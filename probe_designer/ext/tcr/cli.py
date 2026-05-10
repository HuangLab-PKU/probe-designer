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
from probe_designer.ext.tcr.pipeline import run_tcr_pipeline

logger = logging.getLogger(__name__)


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
    rt_primer_gap: int = typer.Option(
        15, "--rt-primer-gap",
        help="Gap in nt between BDS 3' end and RT primer 5' end on target "
             "sequence (default: 15)."
    ),
    rt_primer_tm_range: str = typer.Option(
        "50,65", "--rt-primer-tm-range",
        help="RT primer Tm window 'low,high' (default: 50,65)."
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
            start_no=start_no,
            last_no_from=last_no_from,
            rt_primer_gap=rt_primer_gap,
            rt_primer_tm_range=_parse_pair(rt_primer_tm_range, kind=float),
            skip_blast=skip_blast,
            plot_tm_landscape=plot_tm_landscape,
            codebook=codebook,
        )

    logger.info("Running TCR pipeline with chemistries=%s", cfg.chemistries)
    result = run_tcr_pipeline(cfg)

    msg_parts = [
        f"Done. {result.clones_with_sites}/{result.clones_total} clones produced sites; "
        f"{result.sites_selected} sites; {result.probes_total} probes (Nos "
        f"{result.last_no - result.probes_total + 1}-{result.last_no})"
    ]
    if result.rt_primers_xlsx:
        msg_parts.append(f"RT primers → {result.rt_primers_xlsx.name}")
    if result.skipped_clones:
        msg_parts.append(f"Skipped {len(result.skipped_clones)} clone(s): "
                          + ",".join(result.skipped_clones))
    typer.echo(" | ".join(msg_parts))
    typer.echo(f"Final probes: {result.final_probes_xlsx}")
