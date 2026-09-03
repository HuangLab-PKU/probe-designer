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
from probe_designer.ext.tcr.nilsson import nilsson_terms
from probe_designer.ext.tcr.repertoire import (
    build_repertoire_index,
    is_clone_specific,
)
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
    """In-memory summary of a finished TCR pipeline run.

    Each chemistry is selected independently. Probes share ``No.`` per
    ``(clone, site_index)`` — so ``probes_total`` can exceed the number of
    unique Nos when dual-chemistry pairs collapse. ``first_no`` / ``last_no``
    bound the contiguous No. range actually consumed.
    """
    workdir: Path
    run_dir: Path
    chemistries: List[str]
    clones_total: int
    clones_with_sites: int
    sites_selected: int
    probes_total: int
    first_no: int
    last_no: int
    final_probes_xlsx: Path
    final_probes_fasta: Path
    rt_primers_xlsx: Optional[Path] = None
    landscape_pdf: Optional[Path] = None
    skipped_clones: List[str] = field(default_factory=list)
    sites_selected_per_chem: Dict[str, int] = field(default_factory=dict)


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
) -> tuple[Dict[str, List[Dict]], Dict[str, Dict[str, List[Dict]]], List[str]]:
    """Independent BDS selection per chemistry.

    Returns
    -------
    all_sites
        ``{clone_name: [site_dict, ...]}`` — every scanned position with
        both chemistries' Tm/arm data attached. Used by the Tm landscape PDF.
    selected_per_chem
        ``{chemistry: {clone_name: [selected_site, ...]}}``. Each chemistry's
        site list is filtered + scored + non-overlap-selected ON ITS OWN
        target (mRNA for dRNA, RT-cDNA for cDNA). The two chemistries may
        emit DIFFERENT BDS positions on the same clone — by design.
    skipped
        Clones that produced no candidates (CDR3 not found, empty scan, etc).
    """
    designer = TcrProbeDesigner(
        bds_len=cfg.bds_len, arm_len=cfg.arm_len,
        # min/max_arm_tm on the designer object are unused — pipeline calls
        # filter_by_chemistry_tm with per-chemistry tm_range from cfg.
        min_arm_tm=cfg.tm_range[0], max_arm_tm=cfg.tm_range[1],
    )
    all_sites: Dict[str, List[Dict]] = {}
    selected_per_chem: Dict[str, Dict[str, List[Dict]]] = {
        chem: {} for chem in cfg.chemistries
    }
    skipped: List[str] = []

    repertoire_index = {}
    if cfg.repertoire:
        repertoire_index = build_repertoire_index(cfg.repertoire, kmer=cfg.bds_len)
        logger.info("Repertoire screen: %d clones, %d indexed %d-mers",
                    len(cfg.repertoire), len(repertoire_index), cfg.bds_len)

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
        for s in sites:
            s["clone_id"] = clone_id
        all_sites[clone_name] = sites

        # Clone specificity, measured: drop sites whose footprint also occurs
        # in another clone of this patient's repertoire. Runs before the margin
        # rule because it is the direct evidence the margin only approximates.
        if repertoire_index:
            unique = [s for s in sites
                      if is_clone_specific(s["target_sequence"],
                                           repertoire_index, allowed)]
            if not unique:
                logger.warning(
                    "[%s] every candidate site also occurs in another clone of "
                    "the repertoire — skipping clone", clone_name,
                )
                skipped.append(clone_name)
                continue
            if len(unique) != len(sites):
                logger.info("[%s] repertoire screen: %d/%d candidates unique",
                            clone_name, len(unique), len(sites))
            sites = unique

        # Clone specificity: keep the ligase's ~3 nt clamp inside the CDR3.
        # Applied after all_sites so the Tm landscape still shows the full
        # scan, and reported so a shrunken candidate pool is never silent.
        if cfg.min_ligation_margin:
            m = cfg.min_ligation_margin
            eligible = [
                s for s in sites
                if s["ligation_point"] - s["cdr3_start"] >= m
                and s["cdr3_end"] - s["ligation_point"] >= m
            ]
            if not eligible:
                logger.warning(
                    "[%s] no site keeps %d nt of CDR3 either side of the nick "
                    "(CDR3 is %d nt) — skipping clone",
                    clone_name, m, len(allowed),
                )
                skipped.append(clone_name)
                continue
            logger.info("[%s] ligation margin >= %d nt: %d/%d candidates kept",
                        clone_name, m, len(eligible), len(sites))
            sites = eligible

        # Soft Tm design (2026-07-19): no hard Tm gate. Every scanned site is a
        # candidate; arm-Tm proximity to the reaction-anchored target enters the
        # SCORE (two-sided) and select_non_overlapping_sites picks the score
        # peaks (min_gap = bds_len => non-overlapping). Each chemistry-pass is
        # independent, so we score+select per-chemistry on COPIES of the sites.
        target_tm = designer.reaction.min_arm_tm
        for chem in cfg.chemistries:
            ranked = []
            for src in sites:
                s = dict(src)  # shallow copy so per-chem fields don't leak
                s["chemistry"] = chem
                avg_tm = (s[f"tm_5prime_{chem}"] + s[f"tm_3prime_{chem}"]) / 2.0
                tm_dev = abs(avg_tm - target_tm)  # 0 at the target, grows both ways
                # Higher is better: less structure (mfe), balanced arms
                # (-tm_diff), and arm Tm near the target (-tm_dev).
                # `or 0.0`: mfe is None when ViennaRNA is unavailable (audit R8),
                # in which case the structure term drops out rather than raising.
                s["score"] = round(
                    (s.get("mfe") or 0.0) - s.get(f"tm_diff_{chem}", 0.0) - tm_dev, 3,
                )
                # Magoulopoulou/Nilsson site quality (opt-in): re-ranks only,
                # never rejects, so a clone with no compliant site still ships.
                if cfg.nilsson_weights is not None:
                    terms = nilsson_terms(
                        s["target_sequence"], cfg.arm_len, chem,
                        weights=cfg.nilsson_weights,
                    )
                    s.update(terms)
                    s["score"] = round(s["score"] + terms["nilsson_score"], 3)
                # tm_cDNA_warning kept for downstream readers (informational)
                s["tm_cDNA_warning"] = not (50.0 <= s["tm_cDNA"] <= 75.0)
                ranked.append(s)
            selected = designer.select_non_overlapping_sites(
                ranked, max_sites=cfg.sites_per_clone,
            )
            if selected:
                selected_per_chem[chem][clone_name] = selected
                logger.info("[%s/%s] %d candidates → %d peak sites at %s",
                            clone_name, chem, len(sites),
                            len(selected), [s["st"] for s in selected])

    # Save bds_candidate.xlsx (all scanned positions; chemistry-agnostic)
    flat_all = [s for sites in all_sites.values() for s in sites]
    if flat_all:
        pd.DataFrame(flat_all).to_excel(
            run_dir / "bds_candidate.xlsx", index=False, engine="openpyxl",
        )

    # Save bds_selected.xlsx — flatten per-chemistry selections, tagged.
    flat_sel = [
        s
        for chem in cfg.chemistries
        for sites in selected_per_chem[chem].values()
        for s in sites
    ]
    if not flat_sel:
        raise RuntimeError(
            "No sites selected across all clones/chemistries — check input + Tm windows."
        )
    pd.DataFrame(flat_sel).to_excel(
        run_dir / "bds_selected.xlsx", index=False, engine="openpyxl",
    )
    per_chem_counts = {chem: sum(len(v) for v in selected_per_chem[chem].values())
                       for chem in cfg.chemistries}
    logger.info("Phase 1: %d sites selected → bds_selected.xlsx (%s)",
                len(flat_sel),
                ", ".join(f"{c}: {n}" for c, n in per_chem_counts.items()))

    return all_sites, selected_per_chem, skipped


# ---------------------------------------------------------------------------
# Phase 2 — BLAST QC (informational; never rejects)
# ---------------------------------------------------------------------------


def _phase2_blast(
    selected_per_chem: Dict[str, Dict[str, List[Dict]]],
    cfg: TcrConfig, run_dir: Path,
) -> Dict[str, pd.DataFrame]:
    """Run BLAST per chemistry. Each chemistry BLASTs its OWN selected arms
    against the transcriptome; annotations are written back onto the
    chemistry-specific site dicts. Informational only — no rejections.

    Returns a ``{chemistry: bds_df}`` mapping for phase 3 to assemble from.
    """
    out: Dict[str, pd.DataFrame] = {}

    if cfg.skip_blast:
        logger.info("Phase 2: --skip-blast set, leaving BLAST annotations empty")
        for chem in cfg.chemistries:
            flat = [s for sites in selected_per_chem[chem].values() for s in sites]
            for s in flat:
                s.setdefault("blast_hit_count", 0)
                s.setdefault("blast_top_hits", "")
            out[chem] = pd.DataFrame(flat)
        return out

    blast_dir = run_dir / "blast_results"
    blast_dir.mkdir(parents=True, exist_ok=True)

    for chem in cfg.chemistries:
        flat = [s for sites in selected_per_chem[chem].values() for s in sites]
        if not flat:
            out[chem] = pd.DataFrame()
            continue
        # BLAST the actual arm DNA strands of this chemistry against the
        # human transcriptome reference. For dRNA, the arms are RC(target);
        # for cDNA, they are sense target. Either way, BLAST against the
        # mRNA reference detects off-target hybridization risk.
        seqs = []
        for s in flat:
            binding = (s[f"arm_3prime_{chem}"] + s[f"arm_5prime_{chem}"]).upper()
            seqs.append({
                "probe_id": f"{s['gene_name']}_pos{s['pos']}_{chem}",
                "binding_sequence": binding,
                "gene": s["gene_name"],
            })
        logger.info("Phase 2 [%s]: BLAST QC on %d selected sites", chem, len(seqs))
        run_batch_blast_search(seqs, blast_dir, label=f"TCR_{chem}")
        blast = parse_blast_results(seqs, blast_dir)
        for s, q in zip(flat, seqs):
            hits = blast.get(q["probe_id"], [])
            s["blast_hit_count"] = len(hits)
            s["blast_top_hits"] = "; ".join(
                h.get("subject_id", "")[:60] for h in hits[:3]
            )
        out[chem] = pd.DataFrame(flat)

    # Combined bds_selected_blast.xlsx with chemistry-tagged rows.
    combined = pd.concat([df for df in out.values() if not df.empty], ignore_index=True)
    if not combined.empty:
        combined.to_excel(run_dir / "bds_selected_blast.xlsx",
                          index=False, engine="openpyxl")
        logger.info("Phase 2: bds_selected_blast.xlsx written")
    return out


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
        # Present only when the Nilsson contribution was enabled; NA otherwise
        # so the column never lies about a site that was never assessed.
        "nilsson_pass": bds_row.get("nilsson_pass", pd.NA),
        "homopolymer_run": bds_row.get("homopolymer_run", pd.NA),
        "junction_base": bds_row.get("junction_base", pd.NA),
        "nilsson_gc_content": bds_row.get("gc_content", pd.NA),
    }


def _phase3_assemble(
    bds_per_chem: Dict[str, pd.DataFrame], clones_df: pd.DataFrame,
    cfg: TcrConfig, run_dir: Path, experiment_dir: Path,
) -> tuple[pd.DataFrame, Path, Path, int]:
    """Assemble probes per clone, pairing chemistries by site_index.

    For each clone (in input order), assign one ``No.`` per ``site_index``
    that the clone has across its chemistries. Within each (clone, site_idx),
    every chemistry that has a site at that index emits a probe row, all
    sharing the SAME ``No.`` + backbone barcode — they decode to the same
    codebook ID (different physical probes, redundant readout).

    Worked example — clonotype106 with dRNA=[256, 301], cDNA=[279]:
        No.X:   dRNA pos256 (site#0) + cDNA pos279 (site#0)  — paired
        No.X+1: dRNA pos301 (site#1) alone (no cDNA at site#1)
    """
    backbone_map = _load_backbone_map(cfg.backbone_file)
    next_no = cfg.resolve_start_no()
    codebook = cfg.resolve_codebook()

    # Group each chemistry's bds rows by clone (= gene_name in this codebase).
    # Within a clone, the rows are already sorted by `pos` ascending from
    # select_non_overlapping_sites — preserve that order for stable site_index.
    per_chem_per_clone: Dict[str, Dict[str, List[Dict]]] = {}
    for chem in cfg.chemistries:
        df = bds_per_chem.get(chem, pd.DataFrame())
        per_chem_per_clone[chem] = {}
        if df.empty:
            continue
        for clone_name, sub in df.groupby("gene_name", sort=False):
            per_chem_per_clone[chem][str(clone_name)] = (
                sub.sort_values("pos").to_dict("records")
            )

    clones_ordered = list(clones_df["consensus_id"].astype(str))
    rows: List[Dict] = []
    for clone_name in clones_ordered:
        per_chem_sites = {
            chem: per_chem_per_clone[chem].get(clone_name, [])
            for chem in cfg.chemistries
        }
        max_idx = max((len(v) for v in per_chem_sites.values()), default=0)
        for site_idx in range(max_idx):
            no = next_no
            if no not in backbone_map:
                raise ValueError(
                    f"Backbone No. {no} not found in {cfg.backbone_file.name}"
                )
            bc_name, bc_seq = backbone_map[no]
            for chem in cfg.chemistries:
                sites = per_chem_sites[chem]
                if site_idx < len(sites):
                    rows.append(_assemble_one(
                        chem, sites[site_idx], no, bc_name, bc_seq,
                        cfg.backbone_file.name, codebook,
                    ))
            next_no += 1

    if not rows:
        raise RuntimeError("No probes assembled — check phase 1/2 outputs.")
    combined = pd.DataFrame(rows).sort_values("No.").reset_index(drop=True)
    combined = apply_final_column_order(combined, kind="padlock")

    xlsx = run_dir / "probes_combined.xlsx"
    combined.to_excel(xlsx, index=False, engine="openpyxl")
    fasta = run_dir / "probes_combined.fasta"
    with open(fasta, "w") as f:
        for _, r in combined.iterrows():
            f.write(f">{r['probe_name']}\n{r['probe_sequence']}\n")
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

    # Phase 1 — find + select (independent per chemistry)
    all_sites, selected_per_chem, skipped = _phase1_find_select(
        clones_df, cfg, run_dir,
    )

    # Phase 1b — Tm landscape PDF (optional). Uses dRNA selection for the
    # marker overlay when present; cDNA-only panels overlay cDNA instead.
    landscape_pdf: Optional[Path] = None
    if cfg.plot_tm_landscape:
        clone_metadata = {
            row["consensus_id"]: {
                "clone_id": row["Clone_name"], "type": row.get("Type", ""),
            }
            for _, row in clones_df.iterrows()
        }
        marker_chem = "dRNA" if "dRNA" in cfg.chemistries else cfg.chemistries[0]
        if not selected_per_chem[marker_chem]:
            logger.warning(
                "Phase 1b: marker chemistry %r has zero selected sites; "
                "Tm landscape will be drawn without selected-site overlays.",
                marker_chem,
            )
        landscape_pdf = plot_tm_landscape(
            all_sites, selected_per_chem[marker_chem],
            pdf_path=run_dir / "tm_landscape.pdf",
            tm_min=cfg.tm_gate_for(marker_chem)[0],
            tm_max=cfg.tm_gate_for(marker_chem)[1],
            clone_metadata=clone_metadata,
        )

    # Phase 2 — BLAST QC per chemistry (informational)
    bds_per_chem = _phase2_blast(selected_per_chem, cfg, run_dir)

    # Phase 3 — assemble per clone, pairing chemistries by site_index
    combined, xlsx, fasta, last_no = _phase3_assemble(
        bds_per_chem, clones_df, cfg, run_dir, experiment_dir,
    )

    # Phase 4 — RT primers (only when cDNA in chemistries)
    rt_xlsx: Optional[Path] = None
    if cfg.has_cdna():
        rt_xlsx = _phase4_rt_primers(combined, clones_df, cfg, run_dir)

    clones_with_any_site = {
        clone for chem in cfg.chemistries
        for clone in selected_per_chem[chem].keys()
    }
    return TcrResult(
        workdir=workdir,
        run_dir=run_dir,
        chemistries=cfg.chemistries,
        clones_total=len(clones_df),
        clones_with_sites=len(clones_with_any_site),
        sites_selected=sum(
            len(s) for chem in cfg.chemistries
            for s in selected_per_chem[chem].values()
        ),
        sites_selected_per_chem={
            chem: sum(len(s) for s in selected_per_chem[chem].values())
            for chem in cfg.chemistries
        },
        probes_total=len(combined),
        first_no=int(combined["No."].min()),
        last_no=last_no,
        final_probes_xlsx=xlsx,
        final_probes_fasta=fasta,
        rt_primers_xlsx=rt_xlsx,
        landscape_pdf=landscape_pdf,
        skipped_clones=skipped,
    )
