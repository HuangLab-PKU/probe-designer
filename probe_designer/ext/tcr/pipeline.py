"""TCR padlock pipeline orchestrator.

Runs the four phases (find BDS → BLAST QC → assemble → RT primer) inside a
single Python process. Output layout matches the templates this replaces:

  output/<YYYYMMDD_HHMMSS>_TCR/
    bds_candidate.xlsx          (all scanned positions)
    bds_selected.xlsx           (post-Tm-filter, non-overlap-selected)
    bds_selected_blast.xlsx     (with BLAST hit annotations) — only if not skip_blast
    tm_landscape.pdf            — only if plot_tm_landscape
    probes_combined.{xlsx,fasta} (final probes, all selected chemistries)
    rt_primers.{xlsx,fasta}     — only if 'cDNA' in chemistries
  output/latest_run.txt         (pointer to the latest <ts>_TCR subdir)
  output/last_no.txt            (last-assigned No.; chain via --last-no-from)
  <experiment>/probe_info.xlsx  (mirror of probes_combined)

Most callers just want :func:`run_tcr_pipeline`.
"""
from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from probe_designer.blast.gene_aware import (
    parse_blast_results,
    run_batch_blast_search,
)
from probe_designer.ext.tcr.config import TcrConfig
from probe_designer.ext.tcr.landscape import plot_tm_landscape
from probe_designer.ext.tcr.probe import (
    TcrProbeDesigner,
    validate_subseq_in_target,
)
from probe_designer.io.probe_schema import (
    CHEM_RT,
    apply_final_column_order,
    make_padlock_name,
    make_rt_primer_name,
)
from probe_designer.probe_assembly import assemble_plain_padlock
from probe_designer.rt_primer import design_rt_primer_from_target

logger = logging.getLogger(__name__)


@dataclass
class TcrResult:
    """In-memory summary of a finished TCR pipeline run."""
    workdir: Path
    run_dir: Path
    chemistries: List[str]
    clones_total: int
    clones_with_sites: int
    sites_selected: int
    probes_total: int
    last_no: int
    final_probes_xlsx: Path
    final_probes_fasta: Path
    rt_primers_xlsx: Optional[Path] = None
    landscape_pdf: Optional[Path] = None
    skipped_clones: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Input loading — auto-detect TRA_nt / TRB_nt / chain_nt
# ---------------------------------------------------------------------------


_CHAIN_COLUMN_CANDIDATES = ("TRB_nt", "TRA_nt", "chain_nt", "TR_nt", "sequence_nt")


def _load_clones(input_xlsx: Path) -> pd.DataFrame:
    """Load TCR clones xlsx; require consensus_id + cdr3_nt + a chain column.

    Auto-detects the chain column from
    ``TRB_nt / TRA_nt / chain_nt / TR_nt / sequence_nt`` and stores it in a
    standardized ``chain_nt`` column on the returned DataFrame, while
    preserving the original column under its original name.
    """
    df = pd.read_excel(input_xlsx)
    if "consensus_id" not in df.columns:
        raise ValueError(
            f"{input_xlsx.name} missing required column 'consensus_id'"
        )
    if "cdr3_nt" not in df.columns:
        raise ValueError(
            f"{input_xlsx.name} missing required column 'cdr3_nt' (the "
            "must-fall-in-here sub-region)"
        )
    chain_col = next((c for c in _CHAIN_COLUMN_CANDIDATES if c in df.columns), None)
    if chain_col is None:
        raise ValueError(
            f"{input_xlsx.name} missing a chain-sequence column. Expected one of: "
            f"{_CHAIN_COLUMN_CANDIDATES}"
        )
    df = df.copy()
    df["chain_nt"] = df[chain_col].astype(str)
    df["_chain_col_used"] = chain_col
    if "Clone_name" not in df.columns:
        df["Clone_name"] = df["consensus_id"].astype(str)
    if "Type" not in df.columns:
        df["Type"] = ""
    return df


# ---------------------------------------------------------------------------
# Phase 1 — find BDS + filter + select
# ---------------------------------------------------------------------------


def _phase1_find_select(
    clones_df: pd.DataFrame, cfg: TcrConfig, run_dir: Path,
) -> tuple[Dict[str, List[Dict]], Dict[str, List[Dict]], List[str]]:
    """Returns (all_sites_per_clone, selected_per_clone, skipped_clones)."""
    designer = TcrProbeDesigner(
        bds_len=cfg.bds_len, arm_len=cfg.arm_len,
        min_arm_tm=cfg.tm_range[0], max_arm_tm=cfg.tm_range[1],
    )
    all_sites: Dict[str, List[Dict]] = {}
    selected_per_clone: Dict[str, List[Dict]] = {}
    skipped: List[str] = []

    for _, row in clones_df.iterrows():
        clone_name = str(row["consensus_id"])
        target_seq = str(row["chain_nt"])
        allowed = str(row["cdr3_nt"])
        clone_id = str(row["Clone_name"])

        try:
            validate_subseq_in_target(target_seq, allowed)
        except ValueError as e:
            logger.warning("[%s] %s — skipping clone", clone_name, e)
            skipped.append(clone_name)
            continue

        sites = designer.find_binding_sites(target_seq, allowed, clone_name)
        if not sites:
            logger.warning("[%s] no scan candidates produced", clone_name)
            skipped.append(clone_name)
            continue
        all_sites[clone_name] = sites

        passed = designer.filter_by_mrna_tm(sites)
        if not passed:
            logger.warning("[%s] %d candidates but none passed Tm filter %s",
                           clone_name, len(sites), cfg.tm_range)
            continue
        # Score = mfe - tm_diff_mRNA (higher is better; less negative MFE OR
        # smaller arm-Tm-imbalance both raise the score).
        for s in passed:
            s["score"] = round(s.get("mfe", 0.0) - s.get("tm_diff_mRNA", 0.0), 3)
            s["tm_cDNA_warning"] = not (50.0 <= s["tm_cDNA"] <= 75.0)
            s["clone_id"] = clone_id

        selected = designer.select_non_overlapping_sites(
            passed, max_sites=cfg.sites_per_clone,
        )
        if selected:
            selected_per_clone[clone_name] = selected
            logger.info("[%s] %d candidates → %d passed Tm → %d selected sites at %s",
                        clone_name, len(sites), len(passed), len(selected),
                        [s["st"] for s in selected])

    # Save bds_candidate.xlsx (all scanned positions)
    flat_all = [s for sites in all_sites.values() for s in sites]
    if flat_all:
        pd.DataFrame(flat_all).to_excel(
            run_dir / "bds_candidate.xlsx", index=False, engine="openpyxl",
        )
    # Save bds_selected.xlsx
    flat_sel = [s for sites in selected_per_clone.values() for s in sites]
    if flat_sel:
        pd.DataFrame(flat_sel).to_excel(
            run_dir / "bds_selected.xlsx", index=False, engine="openpyxl",
        )
        logger.info("Phase 1: %d clones, %d selected sites → bds_selected.xlsx",
                    len(selected_per_clone), len(flat_sel))
    else:
        raise RuntimeError("No sites selected across all clones — check input + Tm window.")

    return all_sites, selected_per_clone, skipped


# ---------------------------------------------------------------------------
# Phase 2 — BLAST QC (informational; never rejects)
# ---------------------------------------------------------------------------


def _phase2_blast(
    selected_per_clone: Dict[str, List[Dict]], cfg: TcrConfig, run_dir: Path,
) -> pd.DataFrame:
    """Run BLAST on mRNA arms, annotate hit count + top hits, write
    bds_selected_blast.xlsx. Does not drop sites — TCR clonotypes are
    patient-specific and BLAST is informational only.
    """
    flat = [s for sites in selected_per_clone.values() for s in sites]
    if cfg.skip_blast:
        logger.info("Phase 2: --skip-blast set, leaving BLAST annotations empty")
        for s in flat:
            s.setdefault("blast_hit_count", 0)
            s.setdefault("blast_top_hits", "")
        df = pd.DataFrame(flat)
        return df

    blast_dir = run_dir / "blast_results"
    blast_dir.mkdir(parents=True, exist_ok=True)

    # Build BLAST queries: one per selected site, using the mRNA arm
    # concatenation as the binding sequence.
    seqs = []
    for s in flat:
        binding = (s["arm_3prime_dRNA"] + s["arm_5prime_dRNA"]).upper()
        seqs.append({
            "probe_id": s["gene_name"] + f"_pos{s['pos']}",
            "binding_sequence": binding,
            "gene": s["gene_name"],
        })
    logger.info("Phase 2: BLAST QC on %d selected sites", len(seqs))
    run_batch_blast_search(seqs, blast_dir, label="TCR")
    blast = parse_blast_results(seqs, blast_dir)

    # Annotate (informational)
    for s, q in zip(flat, seqs):
        hits = blast.get(q["probe_id"], [])
        s["blast_hit_count"] = len(hits)
        # Keep the top 3 subject IDs (truncated) for quick diagnostic
        s["blast_top_hits"] = "; ".join(
            h.get("subject_id", "")[:60] for h in hits[:3]
        )

    df = pd.DataFrame(flat)
    df.to_excel(run_dir / "bds_selected_blast.xlsx",
                index=False, engine="openpyxl")
    logger.info("Phase 2: bds_selected_blast.xlsx written")
    return df


# ---------------------------------------------------------------------------
# Phase 3 — assemble (per chemistry)
# ---------------------------------------------------------------------------


def _load_backbone_map(backbone_file: Path) -> Dict[int, tuple]:
    df = pd.read_excel(backbone_file)
    required = {"No.", "Barcode", "Sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Backbone file {backbone_file.name} missing columns: {missing}")
    return {int(r["No."]): (r["Barcode"], r["Sequence"]) for _, r in df.iterrows()}


def _assemble_one(chemistry: str, bds_row: Dict, no: int,
                   bc_name: str, bc_seq: str, backbone_file_name: str,
                   codebook: str) -> Dict:
    arm5 = str(bds_row[f"arm_5prime_{chemistry}"]).lower()
    arm3 = str(bds_row[f"arm_3prime_{chemistry}"]).lower()
    probe = assemble_plain_padlock(arm5=arm5, arm3=arm3,
                                    backbone=str(bc_seq).upper())
    clone_name = str(bds_row["gene_name"])
    clone_id = bds_row.get("clone_id", "")
    pos = int(bds_row["pos"])
    return {
        "order": pd.NA,
        "probe_name": make_padlock_name(
            name=str(clone_id) or clone_name, position=pos,
            chemistry=chemistry, codebook=codebook, no=no,
        ),
        "probe_sequence": probe,
        "gene_name": bds_row["gene_name"], "clone_id": clone_id,
        "chemistry": chemistry,
        "position": pos, "target_sequence": bds_row["target_sequence"],
        "probe_arm5": arm5, "probe_arm3": arm3, "probe_length": len(probe),
        "No.": no, "codebook": codebook,
        "backbone_name": bc_name, "backbone_sequence": bc_seq,
        "backbone_file": backbone_file_name,
        "gc_content": bds_row.get("g_content"),
        "tm": bds_row.get(f"tm_{chemistry}"),
        "tm_5prime": bds_row.get(f"tm_5prime_{chemistry}"),
        "tm_3prime": bds_row.get(f"tm_3prime_{chemistry}"),
        "score": bds_row.get("score"),
        "free_energy": bds_row.get("mfe"),
        "blast_hit_count": bds_row.get("blast_hit_count", 0),
        "blast_top_hits": bds_row.get("blast_top_hits", ""),
        "tm_cdna_warning": bds_row.get("tm_cDNA_warning", False),
    }


def _phase3_assemble(
    bds_df: pd.DataFrame, cfg: TcrConfig, run_dir: Path, experiment_dir: Path,
) -> tuple[pd.DataFrame, Path, Path, int]:
    backbone_map = _load_backbone_map(cfg.backbone_file)
    start_no = cfg.resolve_start_no()
    codebook = cfg.resolve_codebook()
    n_chem = len(cfg.chemistries)
    rows: List[Dict] = []
    for i, r in bds_df.reset_index(drop=True).iterrows():
        for offset, chem in enumerate(cfg.chemistries):
            no = start_no + n_chem * i + offset
            if no not in backbone_map:
                raise ValueError(
                    f"Backbone No. {no} not found in {cfg.backbone_file.name}"
                )
            bc_name, bc_seq = backbone_map[no]
            rows.append(_assemble_one(chem, r.to_dict(), no, bc_name, bc_seq,
                                       cfg.backbone_file.name, codebook))

    combined = pd.DataFrame(rows).sort_values("No.").reset_index(drop=True)
    combined = apply_final_column_order(combined, kind="padlock")

    xlsx = run_dir / "probes_combined.xlsx"
    combined.to_excel(xlsx, index=False, engine="openpyxl")
    fasta = run_dir / "probes_combined.fasta"
    with open(fasta, "w") as f:
        for _, r in combined.iterrows():
            f.write(f">{r['probe_name']}\n{r['probe_sequence']}\n")
    # probe_info mirrors at run_dir AND experiment_dir
    combined.to_excel(run_dir / "probe_info.xlsx", index=False, engine="openpyxl")
    combined.to_excel(experiment_dir / "probe_info.xlsx",
                       index=False, engine="openpyxl")

    last_no = int(combined["No."].max())
    (cfg.output_dir / "last_no.txt").write_text(str(last_no), encoding="utf-8")
    logger.info("Phase 3: %d probes (Nos %d-%d) → probes_combined.xlsx",
                len(combined), int(combined["No."].min()), last_no)
    return combined, xlsx, fasta, last_no


# ---------------------------------------------------------------------------
# Phase 4 — RT primers (cDNA only)
# ---------------------------------------------------------------------------


def _phase4_rt_primers(
    combined: pd.DataFrame, clones_df: pd.DataFrame, cfg: TcrConfig,
    run_dir: Path,
) -> Path:
    cdna = combined[combined["chemistry"] == "cDNA"].reset_index(drop=True)
    if cdna.empty:
        raise RuntimeError("RT-primer phase invoked with no cDNA probes")
    chain_by_consensus = dict(zip(
        clones_df["consensus_id"].astype(str),
        clones_df["chain_nt"].astype(str),
    ))
    rows: List[Dict] = []
    for _, probe in cdna.iterrows():
        clone_name = str(probe["gene_name"])
        clone_id = str(probe["clone_id"])
        no = int(probe["No."])
        pos = int(probe["position"])
        target_seq = chain_by_consensus.get(clone_name, "")
        if not target_seq:
            logger.warning("[%s] no chain sequence found for RT primer", clone_name)
            continue
        primer = design_rt_primer_from_target(
            target_seq=target_seq,
            target_end=pos + cfg.bds_len,
            gap=cfg.rt_primer_gap,
            min_primer_len=cfg.rt_primer_min_len,
            max_primer_len=cfg.rt_primer_max_len,
            target_primer_len=cfg.rt_primer_target_len,
            tm_min=cfg.rt_primer_tm_range[0],
            tm_max=cfg.rt_primer_tm_range[1],
        )
        rows.append({
            "order": pd.NA,
            "probe_name": make_rt_primer_name(
                name=str(clone_id) or clone_name, position=pos,
            ),
            "probe_sequence": primer["primer_seq"],
            "gene_name": clone_name, "clone_id": clone_id,
            "chemistry": CHEM_RT,
            "position": pos, "target_sequence": probe.get("target_sequence", ""),
            "probe_length": primer["primer_length"],
            "gc_content": primer["gc_percent"],
            "tm": primer["tm"],
            "paired_padlock_probe_name": probe["probe_name"],
            "paired_padlock_no": no,
            "notes": primer["notes"],
            "primer_start_on_target": primer["primer_target_start"],
            "primer_end_on_target": primer["primer_target_end"],
            "gap_nt": primer["gap_nt"],
        })

    res_df = apply_final_column_order(pd.DataFrame(rows), kind="rt")
    xlsx = run_dir / "rt_primers.xlsx"
    fasta = run_dir / "rt_primers.fasta"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        res_df.to_excel(writer, sheet_name="RT_Primers", index=False)
        if len(res_df):
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
            summary.to_excel(writer, sheet_name="Summary", index=False)
    with open(fasta, "w") as f:
        for _, r in res_df.iterrows():
            if r["probe_sequence"]:
                f.write(f">{r['probe_name']}\n{r['probe_sequence']}\n")
    logger.info("Phase 4: %d RT primers -> rt_primers.xlsx", len(res_df))
    return xlsx


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------


def run_tcr_pipeline(cfg: TcrConfig) -> TcrResult:
    """Run the four-phase TCR pipeline end-to-end."""
    workdir = cfg.output_dir
    workdir.mkdir(parents=True, exist_ok=True)
    experiment_dir = cfg.input_xlsx.parent.parent  # input/ -> experiment/

    # Timestamped run subdir + latest_run.txt
    ts = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    run_dir = workdir / f"{ts}_TCR"
    run_dir.mkdir(parents=True, exist_ok=True)
    (workdir / "latest_run.txt").write_text(run_dir.name, encoding="utf-8")
    logger.info("Run directory: %s", run_dir)

    # Phase 0 — load + validate
    clones_df = _load_clones(cfg.input_xlsx)
    logger.info("Loaded %d clones from %s (chain column: %s)",
                len(clones_df), cfg.input_xlsx.name,
                clones_df["_chain_col_used"].iloc[0])

    # Phase 1 — find + select
    all_sites, selected_per_clone, skipped = _phase1_find_select(
        clones_df, cfg, run_dir,
    )

    # Phase 1b — Tm landscape PDF (optional)
    landscape_pdf: Optional[Path] = None
    if cfg.plot_tm_landscape:
        clone_metadata = {
            row["consensus_id"]: {
                "clone_id": row["Clone_name"], "type": row.get("Type", ""),
            }
            for _, row in clones_df.iterrows()
        }
        landscape_pdf = plot_tm_landscape(
            all_sites, selected_per_clone,
            pdf_path=run_dir / "tm_landscape.pdf",
            tm_min=cfg.tm_range[0], tm_max=cfg.tm_range[1],
            clone_metadata=clone_metadata,
        )

    # Phase 2 — BLAST QC (informational)
    bds_df = _phase2_blast(selected_per_clone, cfg, run_dir)

    # Phase 3 — assemble
    combined, xlsx, fasta, last_no = _phase3_assemble(
        bds_df, cfg, run_dir, experiment_dir,
    )

    # Phase 4 — RT primers (only when cDNA in chemistries)
    rt_xlsx: Optional[Path] = None
    if cfg.has_cdna():
        rt_xlsx = _phase4_rt_primers(combined, clones_df, cfg, run_dir)

    return TcrResult(
        workdir=workdir,
        run_dir=run_dir,
        chemistries=cfg.chemistries,
        clones_total=len(clones_df),
        clones_with_sites=len(selected_per_clone),
        sites_selected=sum(len(s) for s in selected_per_clone.values()),
        probes_total=len(combined),
        last_no=last_no,
        final_probes_xlsx=xlsx,
        final_probes_fasta=fasta,
        rt_primers_xlsx=rt_xlsx,
        landscape_pdf=landscape_pdf,
        skipped_clones=skipped,
    )
