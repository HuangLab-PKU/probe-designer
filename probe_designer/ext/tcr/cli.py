"""``probe-design tcr`` — Typer subcommand for the TCR pipeline.

Thin CLI wrapper around :func:`run_tcr_pipeline`. The CLI builds a
:class:`TcrConfig` from typed flags (or loads one from YAML via ``--config``)
and dispatches.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import typer

from probe_designer.ext.tcr.config import ALL_CHEMISTRIES, TcrConfig
from probe_designer.ext.tcr.repertoire import load_repertoire_file
from probe_designer.ext.tcr.pipeline import run_tcr_pipeline

logger = logging.getLogger(__name__)


def _parse_triple(value: str, sep: str = ",") -> tuple:
    parts = value.split(sep)
    if len(parts) != 3:
        raise typer.BadParameter(
            f"Expected three values separated by {sep!r}, got {value!r}"
        )
    try:
        return tuple(float(p) for p in parts)
    except ValueError as e:
        raise typer.BadParameter(f"Could not parse {value!r}: {e}") from e


def _parse_pair(value: str, *, kind: type, sep: str = ",") -> tuple:
    parts = value.split(sep)
    if len(parts) != 2:
        raise typer.BadParameter(
            f"Expected two values separated by {sep!r}, got {value!r}"
        )
    try:
        return (kind(parts[0]), kind(parts[1]))
    except ValueError as e:
        raise typer.BadParameter(f"Could not parse {value!r}: {e}") from e


def command(
    input_xlsx: Path = typer.Option(
        ..., "--input", "-i", exists=True, dir_okay=False, readable=True,
        help="TCR clones Excel: columns consensus_id, cdr3_nt, and a chain "
             "sequence column (TRB_nt / TRA_nt / chain_nt — auto-detected). "
             "Optional: Clone_name, Type."
    ),
    output_dir: Path = typer.Option(
        ..., "--output", "-o",
        help="Output directory (created if absent). Each run writes to "
             "output/<YYYYMMDD_HHMMSS>_TCR/. last_no.txt + latest_run.txt at "
             "output root."
    ),
    backbone: Path = typer.Option(
        ..., "--backbone", exists=True, dir_okay=False, readable=True,
        help="Backbone Excel: columns No., Barcode, Sequence."
    ),
    chemistries_csv: str = typer.Option(
        "dRNA", "--chemistries",
        help="Comma-separated subset of {dRNA, cDNA}. Default: dRNA only. "
             "Add cDNA for both chemistries; RT primer phase auto-runs "
             "whenever cDNA is in the selection. Legacy 'mRNA' is accepted "
             "and normalized to 'dRNA'."
    ),
    start_no: Optional[int] = typer.Option(
        None, "--start-no",
        help="No. of the first probe; mutually exclusive with --last-no-from."
    ),
    last_no_from: Optional[Path] = typer.Option(
        None, "--last-no-from", exists=True, dir_okay=False, readable=True,
        help="Read last_no.txt from a previous panel; first probe = read+1."
    ),
    codebook: Optional[str] = typer.Option(
        None, "--codebook",
        help="Codebook name (e.g. SP369). Goes into every padlock probe_name "
             "and the `codebook` column. If omitted, inferred from the backbone "
             "filename (e.g. backbone_SP369.xlsx -> SP369)."
    ),
    bds_len: int = typer.Option(
        40, "--bds-len",
        help="Total binding site length on target (default: 40). arm_len = bds_len/2."
    ),
    arm_len: int = typer.Option(
        20, "--arm-len",
        help="Arm length (default: 20). Must equal bds_len/2."
    ),
    sites_per_clone: int = typer.Option(
        3, "--sites-per-clone",
        help="Maximum non-overlapping sites selected per clone (default: 3)."
    ),
    tm_range: str = typer.Option(
        "50,60", "--tm-range",
        help="dRNA arm Tm window 'low,high' (default: 50,60)."
    ),
    cdna_tm_range: str = typer.Option(
        "50,60", "--cdna-tm-range",
        help="cDNA arm Tm window 'low,high' (default: 50,60). cDNA chemistry "
             "has its OWN BDS selection — gated by this window independently "
             "of --tm-range."
    ),
    rt_primer_gap: int = typer.Option(
        15, "--rt-primer-gap",
        help="Gap in nt between BDS 3' end and RT primer 5' end on target "
             "sequence (default: 15)."
    ),
    rt_primer_tm_range: str = typer.Option(
        "50,65", "--rt-primer-tm-range",
        help="RT primer Tm window 'low,high' (default: 50,65)."
    ),
    repertoire: Optional[Path] = typer.Option(
        None, "--repertoire",
        help="Repertoire table (.csv/.xlsx) for the SAME patient, carrying a "
             "chain column and a CDR3 column. Any candidate site whose 40-mer "
             "also occurs in another clone is dropped: both arms would "
             "hybridise there and the nick would be fully paired, so the "
             "ligase closes it. Neither BLAST nor a Tm filter sees that, "
             "because the off-target is another TCR the panel is deliberately "
             "imaging.",
    ),
    min_ligation_margin: int = typer.Option(
        0, "--min-cdr3-margin",
        help="Require this many nt of CDR3 on EACH side of the ligation nick. "
             "SplintR clamps ~3 nt either side, so a nick nearer the CDR3 edge "
             "seals on germline V/J and no longer tells clones apart. Default 0 "
             "reproduces earlier runs; 3 is the SplintR-grounded value.",
    ),
    nilsson_bonus: bool = typer.Option(
        False, "--nilsson-bonus",
        help="Let the Magoulopoulou 2026 (Nilsson lab) padlock rules re-rank "
             "candidate sites: no homopolymer >= 4 nt, G/C at the ligation "
             "junction, GC 40-60%. Scoring only — no site is ever rejected, so "
             "clones whose CDR3 admits no compliant site still get a probe. "
             "Off by default so earlier runs reproduce exactly.",
    ),
    nilsson_weights: str = typer.Option(
        "2,2,2", "--nilsson-weights",
        help="Weights 'homopolymer,junction_gc,gc' for --nilsson-bonus. "
             "Ignored unless --nilsson-bonus is passed.",
    ),
    skip_blast: bool = typer.Option(
        False, "--skip-blast",
        help="Skip the NCBI BLAST QC round. BLAST is informational for TCR "
             "(annotates hit counts, never rejects), so skipping it is safe "
             "when iterating on Tm parameters."
    ),
    plot_tm_landscape: bool = typer.Option(
        True, "--plot-tm-landscape/--no-plot-tm-landscape",
        help="Render the per-clone Tm landscape PDF (default: on)."
    ),
    config_path: Optional[Path] = typer.Option(
        None, "--config", exists=True, dir_okay=False, readable=True,
        help="YAML config; takes precedence over individual flags except "
             "input/output/backbone/skip_blast."
    ),
    check_cross_lig: bool = typer.Option(
        False, "--check-cross-lig",
        help="Run primer3 DNA cross-ligation screen on the final probes. "
             "Pair with --target-pool to also screen against an existing pool. "
             "Writes cross_lig_report.tsv next to the final probes xlsx.",
    ),
    target_pool: Optional[str] = typer.Option(
        None, "--target-pool",
        help="Pool ID (e.g. `TNBCmarker_main_v1`) to screen new candidates "
             "against. Requires --check-cross-lig.",
    ),
    repo_root: Optional[Path] = typer.Option(
        None, "--repo-root",
        help="Project root containing pools/ and probes/. Defaults to cwd.",
    ),
    reject_cross_lig: bool = typer.Option(
        False, "--reject-cross-lig",
        help="Drop flagged probes from the final xlsx. Requires --check-cross-lig.",
    ),
    xlig_tm_threshold: float = typer.Option(
        27.0, "--xlig-tm-threshold",
        help="Minimum overall dimer Tm (°C) to flag a pair.",
    ),
    xlig_end_dg_threshold: Optional[float] = typer.Option(
        None, "--xlig-end-dg-threshold",
        help="DEPRECATED — v2 uses junction-geometry instead of end-stability. "
             "Accepted as a no-op; will be removed in a future release.",
    ),
) -> None:
    """Run the TCR padlock probe design pipeline (4 phases)."""
    if config_path is not None:
        cfg = TcrConfig.from_yaml(config_path)
        # CLI overrides for deployment-specific paths
        cfg.input_xlsx = input_xlsx
        cfg.output_dir = output_dir
        cfg.backbone_file = backbone
        if skip_blast:
            cfg.skip_blast = True
        if codebook is not None:
            cfg.codebook = codebook
    else:
        chemistries: List[str] = [
            c.strip() for c in chemistries_csv.split(",") if c.strip()
        ]
        cfg = TcrConfig(
            input_xlsx=input_xlsx,
            output_dir=output_dir,
            backbone_file=backbone,
            chemistries=chemistries,
            bds_len=bds_len,
            arm_len=arm_len,
            sites_per_clone=sites_per_clone,
            tm_range=_parse_pair(tm_range, kind=float),
            cdna_tm_range=_parse_pair(cdna_tm_range, kind=float),
            start_no=start_no,
            last_no_from=last_no_from,
            rt_primer_gap=rt_primer_gap,
            rt_primer_tm_range=_parse_pair(rt_primer_tm_range, kind=float),
            skip_blast=skip_blast,
            plot_tm_landscape=plot_tm_landscape,
            codebook=codebook,
            min_ligation_margin=min_ligation_margin,
            repertoire=(load_repertoire_file(repertoire)
                        if repertoire is not None else None),
            nilsson_weights=(
                _parse_triple(nilsson_weights) if nilsson_bonus else None
            ),
        )

    if target_pool and not check_cross_lig:
        raise typer.BadParameter("--target-pool requires --check-cross-lig")
    if reject_cross_lig and not check_cross_lig:
        raise typer.BadParameter("--reject-cross-lig requires --check-cross-lig")

    logger.info("Running TCR pipeline with chemistries=%s", cfg.chemistries)
    result = run_tcr_pipeline(cfg)

    msg_parts = [
        f"Done. {result.clones_with_sites}/{result.clones_total} clones produced sites; "
        f"{result.sites_selected} sites; {result.probes_total} probes (Nos "
        f"{result.first_no}-{result.last_no})"
    ]
    if result.rt_primers_xlsx:
        msg_parts.append(f"RT primers → {result.rt_primers_xlsx.name}")
    if result.skipped_clones:
        msg_parts.append(f"Skipped {len(result.skipped_clones)} clone(s): "
                          + ",".join(result.skipped_clones))
    typer.echo(" | ".join(msg_parts))
    typer.echo(f"Final probes: {result.final_probes_xlsx}")

    if xlig_end_dg_threshold is not None:
        import warnings
        warnings.warn(
            "--xlig-end-dg-threshold is deprecated in v2 (junction-geometry "
            "criterion replaces end-stability). Accepted as no-op.",
            DeprecationWarning, stacklevel=2,
        )

    if check_cross_lig:
        from probe_designer.qc.assemble_hook import apply_cross_lig_check_to_xlsx

        splint_probes = None
        if target_pool:
            from probe_book.root import resolve_root
            from probe_designer.ext.pool.loader import load_pool_as_probes_for_screen
            root = resolve_root(repo_root)
            try:
                splint_probes = load_pool_as_probes_for_screen(target_pool, root)
            except (ImportError, FileNotFoundError) as exc:
                raise typer.BadParameter(str(exc)) from exc
            logger.info("Loaded %d splint probes from pool %s.",
                         len(splint_probes), target_pool)

        n_hits, n_dropped = apply_cross_lig_check_to_xlsx(
            result.final_probes_xlsx,
            splint_probes=splint_probes,
            reject=reject_cross_lig,
            tm_threshold_c=xlig_tm_threshold,
        )
        typer.echo(
            f"Cross-lig v2: {n_hits} pair(s) flagged; {n_dropped} probe(s) dropped"
            f" (pool={target_pool or '(within-batch only)'})"
        )
