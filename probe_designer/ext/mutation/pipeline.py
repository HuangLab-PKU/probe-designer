"""Mutation pipeline orchestrator.

Runs the four phases (thermal-filter -> BLAST -> stitch -> RT primer) inside a
single Python process, sharing the workdir and the chemistry list across
phases. Output paths match the per-experiment script templates this replaces:

  output/step1_thermal_filter/mutation_probe_thermal_<chem>.json
                              best_mutation_probes_thermal_<chem>.xlsx
  output/step2_blast_filter/   mutation_probe_blast_input_<chem>.fasta
                              blast_results_<chem>/batch_<i>.xml
                              mutation_probe_blast_filter_best_<chem>.xlsx
                              summary_<chem>.json
  output/step3_probe_stitch/  <patient>_mut_final_probes.{xlsx,fasta}
  output/<patient>_mut_final_probes.{xlsx,fasta}     (mirror)
  output/<patient>_mut_cDNA_RT_primers.{xlsx,fasta}  (only if cDNA in chemistries)
  output/last_no.txt
  <experiment>/probe_info.xlsx                       (mirror; one row per probe)

Most callers just want :func:`run_mutation_pipeline`.
"""
from __future__ import annotations

import json
import logging
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from tqdm import tqdm

from probe_designer.blast.gene_aware import (
    apply_gene_aware_filter,
    extract_binding_query,
    find_best_per_mutation_record,
    parse_blast_results,
    run_batch_blast_search,
)
from probe_designer.config import DatabaseConfig, FilterConfig, SearchConfig
from probe_designer.database import DatabaseInterface
from probe_designer.ext.mutation.config import MutationConfig
from probe_designer.ext.mutation.probe import (
    MutationProbeDesigner,
    MutationProbeDesigner3End,
    MutationProbeDesignerInvader,
    assemble_ilock,
)
from probe_designer.genome.local_fasta import build_genome_accessor
from probe_designer.io.probe_schema import (
    CHEM_RT,
    apply_final_column_order,
    make_padlock_name,
    make_rt_primer_name,
)
from probe_designer.probe_assembly import assemble_plain_padlock
from probe_designer.rt_primer import design_rt_primer

logger = logging.getLogger(__name__)


MUT_KEY_COLS = ["Gene", "Chr", "Start", "End", "Ref", "Alt"]


def _find_repo_root(*candidates: Path) -> Path | None:
    """Walk up from each candidate path until a directory containing
    ``data/genome/`` is found. Returns the first hit, or ``None`` if none.
    """
    for cand in candidates:
        cur = Path(cand).resolve()
        if cur.is_file():
            cur = cur.parent
        while cur != cur.parent:
            if (cur / "data" / "genome").is_dir():
                return cur
            cur = cur.parent
    return None


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass
class MutationResult:
    """In-memory summary of a finished mutation pipeline run."""
    workdir: Path
    chemistries: List[str]
    mutations_total: int
    mutations_mappable: int
    probes_total: int
    last_no: int
    final_probes_xlsx: Path
    final_probes_fasta: Path
    rt_primers_xlsx: Path | None = None
    unmappable: List[Dict[str, Any]] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Phase 1 — thermal filter
# ---------------------------------------------------------------------------


_DESIGNER_DISPATCH = {
    "invader": MutationProbeDesignerInvader,
    "3end": MutationProbeDesigner3End,
    "default": MutationProbeDesigner,
}


def _build_filter_config(cfg: MutationConfig, chemistry: str) -> FilterConfig:
    chem = cfg.chemistry_params[chemistry]
    # Schema-v2 (Phase 0): GC% is the active gate; the legacy G-only fields
    # are kept for back-compat but no longer gate. The mutation pipeline
    # historically gated only G-runs (via g_consecutive_max) — A/T/C runs
    # were never gated for mutation probes (their length is constrained by
    # mutation context, not free design). Keep A/T/C effectively disabled
    # here so Phase 0's new homopolymer rule doesn't silently shed legacy
    # mutation probes; mRNA pipeline does enable the full ATCG gate.
    return FilterConfig(
        min_gc_content=cfg.gc_content_range[0],
        max_gc_content=cfg.gc_content_range[1],
        min_g_content=cfg.g_content_range[0],
        max_g_content=cfg.g_content_range[1],
        max_consecutive_g=cfg.g_consecutive_max,
        max_consecutive_a=999,   # effectively off — mutation didn't gate A/T/C historically
        max_consecutive_t=999,
        max_consecutive_c=999,
        min_tm=chem.tm_range[0],
        max_tm=chem.tm_range[1],
        max_tm_diff=cfg.max_tm_diff,
        min_free_energy=cfg.mfe_min,
        check_rna_structure=chem.check_rna_structure,
    )


def _build_plp_params(cfg: MutationConfig, chemistry: str) -> Dict[str, Any]:
    chem = cfg.chemistry_params[chemistry]
    return dict(
        G_min=cfg.g_content_range[0], G_max=cfg.g_content_range[1],
        G_consecutive=cfg.g_consecutive_max,
        Tm_low=chem.tm_range[0], Tm_high=chem.tm_range[1],
        max_Tm_dif=cfg.max_tm_diff,
        mfe_thre=cfg.mfe_min,
        ini_arm_len=cfg.ini_arm_len,
        min_arm_len=cfg.arm_len_range[0],
        max_arm_len=cfg.arm_len_range[1],
        max_arm_len_dif=cfg.max_tm_diff,
    )


def _ensure_strand(mut_df: pd.DataFrame, mut_list: List[Dict],
                   db: DatabaseInterface) -> pd.DataFrame:
    """Ensure mut_list rows have a strand value; cache back into mut_df."""
    if "strand" not in mut_df.columns:
        mut_df["strand"] = None
    updated = 0
    for i, mut in enumerate(tqdm(mut_list, desc="Get strand info")):
        gene = mut["Gene.refGene.x"]
        cached = mut_df.iloc[i].get("strand") if i < len(mut_df) else None
        if pd.notna(cached) and cached in (1, -1, 0):
            strand = int(cached)
        else:
            strand = db.get_strand_by_symbol(gene)
            if strand is not None:
                mut_df.iloc[i, mut_df.columns.get_loc("strand")] = strand
                updated += 1
        mut["strand"] = strand if strand is not None else 1
    logger.info("Strand info: %d newly fetched; remainder cached/defaulted", updated)
    return mut_df


def _run_chemistry_thermal(
    chemistry: str, mut_list: List[Dict], genome_accessor,
    cfg: MutationConfig, workdir: Path,
) -> List[Dict]:
    """Run the per-chemistry designer; emit JSON + best-per-mutation xlsx."""
    chem = cfg.chemistry_params[chemistry]
    plp_params = _build_plp_params(cfg, chemistry)
    filter_cfg = _build_filter_config(cfg, chemistry)
    designer_cls = _DESIGNER_DISPATCH[chem.designer]
    designer = designer_cls(SearchConfig(), filter_cfg, target_type=chem.target_type)
    logger.info("=== %s pass (target=%s, tm=%s, designer=%s) ===",
                chemistry, chem.target_type, chem.tm_range, chem.designer)
    results = designer.design_mutation_probes(mut_list, genome_accessor,
                                              plp_params, pad=cfg.pad)
    # Drop probes that did not pass thermal filter — keep the rest of the
    # mutation record so downstream phases can still see "no probe" entries.
    out = []
    for item in results:
        entry = dict(item)
        entry["plp_ref"] = [p for p in (item.get("plp_ref") or [])
                            if p.get("passed_thermal_filter")]
        entry["plp_mut"] = [p for p in (item.get("plp_mut") or [])
                            if p.get("passed_thermal_filter")]
        out.append(entry)
    _save_step1_outputs(out, chemistry, workdir, chem)
    return out


def _save_step1_outputs(results: List[Dict], chemistry: str,
                        workdir: Path, chem) -> None:
    step1 = workdir / "step1_thermal_filter"
    step1.mkdir(parents=True, exist_ok=True)
    json_path = step1 / f"mutation_probe_thermal_{chemistry}.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump({"params1": results, "chemistry": chemistry}, f,
                  indent=2, default=str, ensure_ascii=False)
    logger.info("[%s] step1 JSON -> %s", chemistry, json_path.name)

    rows: List[Dict] = []
    for item in results:
        mut_probes = item.get("plp_mut") or []
        if not mut_probes:
            continue
        best = min(mut_probes, key=lambda x: x.get("score", float("inf")))
        rows.append({
            "Gene": item.get("Gene.refGene.x", "unknown"),
            "Chr": item.get("Chr", ""),
            "Start": item.get("Start", ""),
            "End": item.get("End", ""),
            "Ref": item.get("Ref", ""),
            "Alt": item.get("Alt", ""),
            "strand": item.get("strand", ""),
            "mut_type": item.get("mut_type", ""),
            "chemistry": chemistry,
            "tm_model": chem.tm_model,
            "mutation_position": chem.mutation_position,
            "plp_pos5_local": best.get("plp_pos5_local"),
            "plp_arm5": best.get("probe_arm5", ""),
            "plp_arm3": best.get("probe_arm3", ""),
            "arm5_len": best.get("arm5_len"),
            "arm3_len": best.get("arm3_len"),
            "tm": best.get("tm"),
            "tm_5prime": best.get("tm_5prime"),
            "tm_3prime": best.get("tm_3prime"),
            "tm_diff": best.get("tm_diff"),
            "g_content": best.get("g_content"),
            "gc_content": best.get("gc_content"),  # Schema-v2 (Phase 0)
            "free_energy": best.get("free_energy"),
            "score": best.get("score"),
            "target_sequence": best.get("target_sequence", ""),
        })
    xlsx = step1 / f"best_mutation_probes_thermal_{chemistry}.xlsx"
    if rows:
        pd.DataFrame(rows).sort_values("Gene").to_excel(
            xlsx, index=False, engine="openpyxl",
        )
        logger.info("[%s] best-per-mutation xlsx: %s (%d rows)",
                    chemistry, xlsx.name, len(rows))
    else:
        logger.warning("[%s] no probes passed thermal filter; skipping xlsx",
                       chemistry)


# ---------------------------------------------------------------------------
# Phase 2 — BLAST gene-aware filter
# ---------------------------------------------------------------------------


def _extract_binding_sequences(records: List[Dict], chemistry: str) -> List[Dict]:
    out: List[Dict] = []
    seen = set()
    for item in records:
        gene = item.get("Gene.refGene.x", "unknown")
        for probe in (item.get("plp_mut") or []):
            arm3 = probe.get("probe_arm3", "")
            arm5 = probe.get("probe_arm5", "")
            binding = extract_binding_query(arm5, arm3, chemistry)
            if not binding or binding in seen:
                continue
            seen.add(binding)
            strand = item.get("strand", "unknown")
            mut_type = item.get("mut_type", "unknown")
            mut_start = item.get("Start", "unknown")
            plp_pos5 = probe.get("plp_pos5_local", 0)
            arm5_len = probe.get("arm5_len", len(arm5))
            arm3_len = probe.get("arm3_len", len(arm3))
            probe_id = (f"{chemistry}_{gene}_strand{strand}_{mut_type}_pos{mut_start}_"
                        f"plp{plp_pos5}_arm5{arm5_len}_arm3{arm3_len}")
            out.append({
                "chemistry": chemistry, "gene": gene, "probe_id": probe_id,
                "binding_sequence": binding,
                "arm3": arm3, "arm5": arm5,
                "original_probe": probe,
                "original_item": item,
            })
    logger.info("[%s] %d unique binding sequences for BLAST", chemistry, len(out))
    return out


def _save_step2_outputs(chemistry: str, best_records: List[Dict],
                         workdir: Path) -> Path:
    step2 = workdir / "step2_blast_filter"
    step2.mkdir(parents=True, exist_ok=True)
    rows: List[Dict] = []
    for r in best_records:
        it = r["item"]
        p = r["probe"]
        rows.append({
            "Gene": it.get("Gene.refGene.x", ""),
            "Chr": it.get("Chr", ""), "Start": it.get("Start", ""),
            "End": it.get("End", ""), "Ref": it.get("Ref", ""),
            "Alt": it.get("Alt", ""), "strand": it.get("strand", ""),
            "mut_type": it.get("mut_type", ""),
            "chemistry": chemistry,
            "plp_arm5": p.get("probe_arm5", ""),
            "plp_arm3": p.get("probe_arm3", ""),
            "arm5_len": p.get("arm5_len", ""),
            "arm3_len": p.get("arm3_len", ""),
            "plp_pos5_local": p.get("plp_pos5_local", ""),
            "tm": p.get("tm", ""),
            "tm_5prime": p.get("tm_5prime", ""),
            "tm_3prime": p.get("tm_3prime", ""),
            "g_content": p.get("g_content", ""),
            "free_energy": p.get("free_energy", ""),
            "score": p.get("score", ""),
            "target_sequence": p.get("target_sequence", ""),
            "blast_hits_count": int(r["seq_info"].get("blast_hits", 0)),
            "same_gene_exact": int(r["seq_info"].get("same_gene_exact", 0)),
        })
    xlsx = step2 / f"mutation_probe_blast_filter_best_{chemistry}.xlsx"
    if rows:
        pd.DataFrame(rows).sort_values("Gene").to_excel(
            xlsx, index=False, engine="openpyxl",
        )
        logger.info("[%s] best-per-mutation BLAST xlsx: %s (%d rows)",
                    chemistry, xlsx.name, len(rows))
    return xlsx


def _run_chemistry_blast(chemistry: str, step1_records: List[Dict],
                          workdir: Path, *, skip_blast: bool) -> Path:
    step2 = workdir / "step2_blast_filter"
    step2.mkdir(parents=True, exist_ok=True)
    seqs = _extract_binding_sequences(step1_records, chemistry)
    if not seqs:
        logger.warning("[%s] no binding sequences; nothing to BLAST", chemistry)
        return _save_step2_outputs(chemistry, [], workdir)

    fasta = step2 / f"mutation_probe_blast_input_{chemistry}.fasta"
    with open(fasta, "w") as f:
        for s in seqs:
            f.write(f">{s['probe_id']}\n{s['binding_sequence']}\n")

    if skip_blast:
        # Promote step1 thermal best directly to step2 best (no BLAST round).
        # Each seq passes with placeholder counts so downstream xlsx columns exist.
        passed = [{**s, "blast_hits": 0, "best_match_percentage": 0,
                   "same_gene_exact": 0} for s in seqs]
        logger.info("[%s] --skip-blast: promoting %d seqs unfiltered",
                    chemistry, len(passed))
    else:
        blast_dir = step2 / f"blast_results_{chemistry}"
        blast_dir.mkdir(parents=True, exist_ok=True)
        run_batch_blast_search(seqs, blast_dir, label=chemistry)
        blast = parse_blast_results(seqs, blast_dir)
        if not blast:
            logger.error("[%s] BLAST returned no parseable results", chemistry)
            return _save_step2_outputs(chemistry, [], workdir)
        passed, rejected = apply_gene_aware_filter(seqs, blast)
        logger.info("[%s] gene-aware filter: %d passed, %d rejected",
                    chemistry, len(passed), len(rejected))

    best_records = find_best_per_mutation_record(passed)
    logger.info("[%s] best-per-mutation: %d selected", chemistry, len(best_records))
    return _save_step2_outputs(chemistry, best_records, workdir)


# ---------------------------------------------------------------------------
# Phase 3 — stitch
# ---------------------------------------------------------------------------


def _load_chemistry_bests(workdir: Path, chemistries: List[str]) -> Dict[str, pd.DataFrame]:
    step2 = workdir / "step2_blast_filter"
    frames: Dict[str, pd.DataFrame] = {}
    for chem in chemistries:
        path = step2 / f"mutation_probe_blast_filter_best_{chem}.xlsx"
        if not path.exists():
            raise FileNotFoundError(f"Required step2 best xlsx missing: {path}")
        frames[chem] = pd.read_excel(path)
        logger.info("Loaded %s best: %d rows", chem, len(frames[chem]))
    return frames


def _merge_chemistry_bests(frames: Dict[str, pd.DataFrame],
                            chemistries: List[str]) -> pd.DataFrame:
    """Inner-join all frames on the mutation key; suffix non-key columns by chem."""
    renamed = {}
    for chem, df in frames.items():
        rename_map = {c: f"{c}_{chem}" for c in df.columns if c not in MUT_KEY_COLS}
        renamed[chem] = df.rename(columns=rename_map)
    merged = renamed[chemistries[0]]
    for chem in chemistries[1:]:
        merged = merged.merge(renamed[chem], on=MUT_KEY_COLS, how="inner")
    # Warn about mutations dropped during the inner join
    key_sets = {chem: set(map(tuple, df[MUT_KEY_COLS].to_numpy().tolist()))
                for chem, df in frames.items()}
    all_keys = set().union(*key_sets.values())
    kept = set(map(tuple, merged[MUT_KEY_COLS].to_numpy().tolist()))
    for k in (all_keys - kept):
        missing = [chem for chem, s in key_sets.items() if k not in s]
        logger.warning("Mutation %s dropped during merge — missing from: %s", k, missing)
    logger.info("Merged: %d mutations present in all %d chemistries",
                len(merged), len(frames))
    return merged.sort_values("Gene").reset_index(drop=True)


def _load_backbone(backbone_file: Path) -> pd.DataFrame:
    df = pd.read_excel(backbone_file)
    required = {"No.", "Barcode", "Sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Backbone file {backbone_file.name} missing columns: {missing}")
    return df


def _build_stitch_row(chemistry: str, merged_row, no: int,
                       bc_name: str, bc_seq: str, ilock_flap: str,
                       backbone_file_name: str, codebook: str) -> Dict[str, Any]:
    suffix = f"_{chemistry}"
    arm5 = str(merged_row[f"plp_arm5{suffix}"]).lower()
    arm3 = str(merged_row[f"plp_arm3{suffix}"]).lower()
    if chemistry == "iLock":
        try:
            probe = assemble_ilock(arm5=arm5, arm3=arm3, backbone=bc_seq, flap=ilock_flap)
        except ValueError as e:
            gene = merged_row["Gene"]; start = merged_row["Start"]
            raise ValueError(f"iLock probe for {gene}@{start}: {e}") from e
        flap, linker = ilock_flap, arm3[-1].upper()
    else:
        probe = assemble_plain_padlock(arm5=arm5, arm3=arm3,
                                        backbone=bc_seq.upper())
        flap, linker = "", ""
        gene = merged_row["Gene"]; start = merged_row["Start"]
        assert probe.startswith(arm5) and probe.endswith(arm3), \
            f"{chemistry} probe for {gene}@{start} violates plain-padlock layout"

    position = int(merged_row["Start"])
    probe_name = make_padlock_name(
        name=str(merged_row["Gene"]), position=position,
        chemistry=chemistry, codebook=codebook, no=no,
    )
    return {
        "order": pd.NA,
        "probe_name": probe_name,
        "probe_sequence": probe,
        "gene_name": merged_row["Gene"], "clone_id": "",
        "chemistry": chemistry,
        "position": position, "target_sequence": merged_row.get(f"target_sequence{suffix}", ""),
        "probe_arm5": arm5, "probe_arm3": arm3, "probe_length": len(probe),
        "No.": no, "codebook": codebook,
        "backbone_name": bc_name, "backbone_sequence": bc_seq,
        "backbone_file": backbone_file_name,
        # Schema-v2 (Phase 0): prefer the real gc_content column; fall back to
        # legacy g_content only for old step1 files that pre-date the fix.
        "gc_content": (
            merged_row.get(f"gc_content{suffix}")
            if pd.notna(merged_row.get(f"gc_content{suffix}", pd.NA))
            else merged_row.get(f"g_content{suffix}")
        ),
        "g_content": merged_row.get(f"g_content{suffix}"),
        "tm": merged_row.get(f"tm{suffix}"),
        "tm_5prime": merged_row.get(f"tm_5prime{suffix}"),
        "tm_3prime": merged_row.get(f"tm_3prime{suffix}"),
        "score": merged_row.get(f"score{suffix}"),
        "free_energy": merged_row.get(f"free_energy{suffix}"),
        "blast_hit_count": merged_row.get(f"blast_hits_count{suffix}", 0),
        "blast_top_hits": pd.NA,
        "same_gene_exact": merged_row.get(f"same_gene_exact{suffix}", 0),
        "chr": merged_row["Chr"],
        "genomic_start": merged_row["Start"], "genomic_end": merged_row["End"],
        "ref": merged_row["Ref"], "alt": merged_row["Alt"],
        "strand": merged_row.get(f"strand{suffix}", ""),
        "mutation_type": merged_row.get(f"mut_type{suffix}", ""),
        "tm_model": merged_row.get(f"tm_model{suffix}", ""),
        "mutation_position": merged_row.get(f"mutation_position{suffix}", ""),
        "arm5_len": merged_row.get(f"arm5_len{suffix}"),
        "arm3_len": merged_row.get(f"arm3_len{suffix}"),
        "plp_pos5_local": merged_row.get(f"plp_pos5_local{suffix}"),
        "ilock_flap": flap, "ilock_linker_nt": linker,
    }


def _stitch(merged_df: pd.DataFrame, backbone_df: pd.DataFrame,
             start_no: int, chemistries: List[str], ilock_flap: str,
             backbone_file_name: str, codebook: str) -> pd.DataFrame:
    backbone_by_no = {int(r["No."]): (r["Barcode"], r["Sequence"])
                      for _, r in backbone_df.iterrows()}
    n_chem = len(chemistries)
    rows: List[Dict] = []
    for i, row in merged_df.iterrows():
        nos = {chem: start_no + n_chem * i + offset
               for offset, chem in enumerate(chemistries)}
        for no in nos.values():
            if no not in backbone_by_no:
                raise ValueError(f"Backbone No. {no} not found in {backbone_file_name}")
        for chem in chemistries:
            bc_name, bc_seq = backbone_by_no[nos[chem]]
            rows.append(_build_stitch_row(
                chem, row, nos[chem], bc_name, bc_seq, ilock_flap,
                backbone_file_name, codebook,
            ))
    return pd.DataFrame(rows).sort_values("No.").reset_index(drop=True)


def _save_step3_outputs(combined: pd.DataFrame, workdir: Path,
                         experiment_dir: Path) -> tuple[Path, Path, int]:
    step3 = workdir / "step3_probe_stitch"
    step3.mkdir(parents=True, exist_ok=True)
    combined = apply_final_column_order(combined, kind="padlock")
    patient_prefix = experiment_dir.name.rsplit("_", 1)[0]
    patient_tag = patient_prefix.split("_")[-1] if "_" in patient_prefix else "BZ"
    out_stem = f"{patient_tag}_mut_final_probes"
    xlsx = step3 / f"{out_stem}.xlsx"
    chemistries_in_run = list(dict.fromkeys(combined["chemistry"]))
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        combined.to_excel(w, sheet_name="Final_Probes", index=False)
        summary_rows = []
        for chem in chemistries_in_run:
            sub = combined[combined["chemistry"] == chem]
            summary_rows.append({
                "Chemistry": chem,
                "Probes": len(sub),
                "Genes": sub["gene_name"].nunique(),
                "Tm_model": sub["tm_model"].iloc[0] if len(sub) else "",
                "Mutation_position": sub["mutation_position"].iloc[0] if len(sub) else "",
                "Avg Length": f"{sub['probe_length'].mean():.1f}" if len(sub) else "",
                "No. range": f"{sub['No.'].min()}-{sub['No.'].max()}" if len(sub) else "",
            })
        pd.DataFrame(summary_rows).to_excel(w, sheet_name="Summary", index=False)
    fasta = step3 / f"{out_stem}.fasta"
    with open(fasta, "w") as f:
        for _, r in combined.iterrows():
            f.write(f">{r['probe_name']}\n{r['probe_sequence']}\n")
    # Mirror to workdir root for legacy scripts that look there
    shutil.copy2(xlsx, workdir / xlsx.name)
    shutil.copy2(fasta, workdir / fasta.name)
    # probe_info mirror at experiment root
    probe_info = experiment_dir / "probe_info.xlsx"
    combined.to_excel(probe_info, index=False)
    last_no = int(combined["No."].max())
    (workdir / "last_no.txt").write_text(str(last_no), encoding="utf-8")
    logger.info("Combined probes: %s (%d rows); last_no=%d",
                xlsx.name, len(combined), last_no)
    return xlsx, fasta, last_no


# ---------------------------------------------------------------------------
# Phase 4 — RT primers (cDNA chemistry only)
# ---------------------------------------------------------------------------


def _run_rt_primers(combined: pd.DataFrame, cfg: MutationConfig,
                     workdir: Path, experiment_dir: Path,
                     genome_accessor) -> Path | None:
    cdna = combined[combined["chemistry"] == "cDNA"].reset_index(drop=True)
    if cdna.empty:
        logger.info("No cDNA probes; skipping RT primer phase")
        return None
    patient_prefix = experiment_dir.name.rsplit("_", 1)[0]
    patient_tag = patient_prefix.split("_")[-1] if "_" in patient_prefix else "BZ"
    rows: List[Dict] = []
    for _, row in cdna.iterrows():
        gene = row["gene_name"]
        no = row["No."]
        position = int(row["genomic_start"])
        primer = design_rt_primer(
            chrom=row["chr"],
            mut_start=position,
            mut_end=int(row["genomic_end"]),
            strand=row.get("strand", 1),
            probe_arm5_len=len(str(row.get("probe_arm5", ""))),
            genome_accessor=genome_accessor,
            gap=cfg.rt_primer_gap,
            min_primer_len=cfg.rt_primer_min_len,
            max_primer_len=cfg.rt_primer_max_len,
            target_primer_len=cfg.rt_primer_target_len,
            tm_min=cfg.rt_primer_tm_range[0],
            tm_max=cfg.rt_primer_tm_range[1],
        )
        rows.append({
            "order": pd.NA,
            "probe_name": make_rt_primer_name(name=str(gene), position=position),
            "probe_sequence": primer["primer_seq"],
            "gene_name": gene, "clone_id": "",
            "chemistry": CHEM_RT,
            "position": position, "target_sequence": pd.NA,
            "probe_length": primer["primer_length"],
            "gc_content": primer["gc_percent"],
            "tm": primer["tm"],
            "paired_padlock_probe_name": row["probe_name"],
            "paired_padlock_no": no,
            "notes": primer["notes"],
            "chr": row["chr"],
            "genomic_start": row["genomic_start"], "genomic_end": row["genomic_end"],
            "ref": row["ref"], "alt": row["alt"],
            "strand": primer["strand"],
            "mutation_type": row.get("mutation_type", ""),
            "primer_genomic_start": primer["primer_genomic_start"],
            "primer_genomic_end": primer["primer_genomic_end"],
            "gap_nt": primer["gap_nt"],
        })
    res_df = apply_final_column_order(pd.DataFrame(rows), kind="rt")
    xlsx = workdir / f"{patient_tag}_mut_cDNA_RT_primers.xlsx"
    fasta = workdir / f"{patient_tag}_mut_cDNA_RT_primers.fasta"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        res_df.to_excel(w, sheet_name="RT_Primers", index=False)
        summary = pd.DataFrame({
            "Metric": [
                "Total primers", "Tm range", "Length range",
                f"Primers with Tm < {cfg.rt_primer_tm_range[0]}C",
                f"Primers with Tm > {cfg.rt_primer_tm_range[1]}C",
                "Primers with notes",
            ],
            "Value": [
                len(res_df),
                f"{res_df['tm'].min():.1f} - {res_df['tm'].max():.1f} C",
                f"{res_df['probe_length'].min()} - {res_df['probe_length'].max()} nt",
                int((res_df["tm"] < cfg.rt_primer_tm_range[0]).sum()),
                int((res_df["tm"] > cfg.rt_primer_tm_range[1]).sum()),
                int((res_df["notes"].fillna("") != "").sum()),
            ],
        })
        summary.to_excel(w, sheet_name="Summary", index=False)
    with open(fasta, "w") as f:
        for _, r in res_df.iterrows():
            if r["probe_sequence"]:
                f.write(f">{r['probe_name']}\n{r['probe_sequence']}\n")
    logger.info("RT primers: %s (%d primers)", xlsx.name, len(res_df))
    return xlsx


# ---------------------------------------------------------------------------
# Unmappable list
# ---------------------------------------------------------------------------


def _emit_unmappable_csv(mut_df: pd.DataFrame, mappable_keys: set,
                          workdir: Path) -> List[Dict]:
    unmappable_rows = []
    for _, row in mut_df.iterrows():
        key = (row.get("Gene.refGene.x"), row.get("Chr"), row.get("Start"),
               row.get("End"), row.get("Ref"), row.get("Alt"))
        if key in mappable_keys:
            continue
        unmappable_rows.append({
            "Gene": row.get("Gene.refGene.x", ""),
            "Chr": row.get("Chr", ""),
            "Start": row.get("Start", ""),
            "End": row.get("End", ""),
            "Ref": row.get("Ref", ""),
            "Alt": row.get("Alt", ""),
            "Reason": "no chemistry produced a passing probe",
        })
    if unmappable_rows:
        path = workdir / "unmappable.csv"
        pd.DataFrame(unmappable_rows).to_csv(path, index=False)
        logger.warning("%d unmappable mutation(s) -> %s",
                       len(unmappable_rows), path.name)
    return unmappable_rows


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------


def run_mutation_pipeline(cfg: MutationConfig) -> MutationResult:
    """Run the four-phase mutation pipeline end-to-end."""
    workdir = cfg.output_dir
    workdir.mkdir(parents=True, exist_ok=True)
    experiment_dir = cfg.input_xlsx.parent.parent  # input/ -> experiment/

    # Phase 0 — strand info + genome accessor + load mutations
    if cfg.genome_path is not None:
        genome_path = cfg.genome_path
    else:
        repo_root = _find_repo_root(workdir, cfg.input_xlsx, Path.cwd())
        if repo_root is None:
            raise FileNotFoundError(
                "Could not auto-locate the genome FASTA. Pass --genome explicitly "
                "or run from a directory whose ancestor contains data/genome/."
            )
        genome_path = cfg.resolve_genome_path(repo_root)
    if not genome_path.exists():
        raise FileNotFoundError(f"Genome FASTA not found: {genome_path}")
    genome_accessor = build_genome_accessor(genome_path)

    mut_df = pd.read_excel(cfg.input_xlsx)
    logger.info("Loaded %d mutations from %s", len(mut_df), cfg.input_xlsx.name)
    mut_list = mut_df[
        ["Gene.refGene.x", "Chr", "Start", "End", "Ref", "Alt"]
    ].to_dict(orient="records")
    db = DatabaseInterface(DatabaseConfig(
        database_type="ensembl", organism="human",
        coord_system_version="GRCh38", max_retries=3,
    ))
    mut_df = _ensure_strand(mut_df, mut_list, db)
    mut_df.to_excel(cfg.input_xlsx, index=False)

    # Phase 1 — thermal filter per chemistry
    step1_results: Dict[str, List[Dict]] = {}
    for chem in cfg.chemistries:
        step1_results[chem] = _run_chemistry_thermal(
            chem, mut_list, genome_accessor, cfg, workdir,
        )

    # Phase 2 — BLAST gene-aware filter per chemistry
    for chem in cfg.chemistries:
        _run_chemistry_blast(chem, step1_results[chem], workdir,
                              skip_blast=cfg.skip_blast)

    # Phase 3 — stitch (inner-merges across chemistries)
    frames = _load_chemistry_bests(workdir, cfg.chemistries)
    merged_df = _merge_chemistry_bests(frames, cfg.chemistries)
    if merged_df.empty:
        raise RuntimeError(
            "No mutations are mappable in ALL selected chemistries; "
            "nothing to stitch. Check per-chemistry step1/step2 outputs."
        )
    backbone_df = _load_backbone(cfg.backbone_file)
    start_no = cfg.resolve_start_no()
    codebook = cfg.resolve_codebook()
    combined = _stitch(merged_df, backbone_df, start_no, cfg.chemistries,
                        cfg.ilock_flap, cfg.backbone_file.name, codebook)
    xlsx, fasta, last_no = _save_step3_outputs(combined, workdir, experiment_dir)

    # Phase 4 — RT primers (cDNA only)
    rt_xlsx = None
    if "cDNA" in cfg.chemistries:
        rt_xlsx = _run_rt_primers(combined, cfg, workdir,
                                   experiment_dir, genome_accessor)

    # Unmappable accounting
    mappable_keys = set(map(tuple, merged_df[MUT_KEY_COLS].to_numpy().tolist()))
    unmappable = _emit_unmappable_csv(mut_df, mappable_keys, workdir)

    return MutationResult(
        workdir=workdir,
        chemistries=cfg.chemistries,
        mutations_total=len(mut_df),
        mutations_mappable=len(merged_df),
        probes_total=len(combined),
        last_no=last_no,
        final_probes_xlsx=xlsx,
        final_probes_fasta=fasta,
        rt_primers_xlsx=rt_xlsx,
        unmappable=unmappable,
    )
