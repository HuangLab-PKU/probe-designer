"""``probe-design mutation`` — Typer subcommand for the mutation pipeline.

Thin CLI wrapper around :func:`run_mutation_pipeline`. The CLI builds a
:class:`MutationConfig` from typed flags (or loads one from YAML via
``--config``) and dispatches.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import typer

from probe_designer.ext.mutation.config import (
    ALL_CHEMISTRIES,
    MutationConfig,
    default_chemistry_params,
)
from probe_designer.ext.mutation.pipeline import run_mutation_pipeline

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
        help="Mutation Excel: columns Gene.refGene.x, Chr, Start, End, Ref, Alt; "
             "optional strand, source_Probe_ID."
    ),
    output_dir: Path = typer.Option(
        ..., "--output", "-o",
        help="Output directory (created if absent). Step1/2/3 subfolders + "
             "final probe xlsx/fasta land here."
    ),
    backbone: Path = typer.Option(
        ..., "--backbone", exists=True, dir_okay=False, readable=True,
        help="Backbone Excel: columns No., Barcode, Sequence."
    ),
    chemistries_csv: str = typer.Option(
        ",".join(ALL_CHEMISTRIES), "--chemistries",
        help="Comma-separated subset of {iLock, mRNA_noiLock, cDNA}."
    ),
    start_no: Optional[int] = typer.Option(
        None, "--start-no",
        help="No. of the first probe; mutually exclusive with --last-no-from."
    ),
    last_no_from: Optional[Path] = typer.Option(
        None, "--last-no-from", exists=True, dir_okay=False, readable=True,
        help="Read last_no.txt from a previous panel; first probe = read+1."
    ),
    ilock_tm_range: Optional[str] = typer.Option(
        None, "--ilock-tm-range",
        help="Override iLock chemistry tm_range as 'low,high' (default: 50,70)."
    ),
    mrna_noilock_tm_range: Optional[str] = typer.Option(
        None, "--mrna-noilock-tm-range",
        help="Override mRNA_noiLock tm_range (default: 55,75)."
    ),
    cdna_tm_range: Optional[str] = typer.Option(
        None, "--cdna-tm-range",
        help="Override cDNA tm_range (default: 45,75)."
    ),
    arm_len_range: str = typer.Option(
        "10,30", "--arm-len-range",
        help="Min,max arm length (nt). Default: 10,30. "
             "Loose-thermal rescue typically uses 8,35."
    ),
    mfe_min: float = typer.Option(
        -15.0, "--mfe-min",
        help="Minimum acceptable RNA free energy (kcal/mol). "
             "Default: -15. Loose-thermal rescue typically uses -25."
    ),
    ilock_flap: str = typer.Option(
        "CGTTGCTGTGGCG", "--ilock-flap",
        help="iLock 5' flap. Default: HuangLab 2026+ (CGTTGCTGTGGCG). "
             "Legacy BZ09 reference: TATATCCCTATAT."
    ),
    rt_primer_gap: int = typer.Option(
        15, "--rt-primer-gap",
        help="Gap in nt between padlock target 3' end and RT primer 5' end on mRNA."
    ),
    rt_primer_tm_range: str = typer.Option(
        "50,65", "--rt-primer-tm-range",
        help="RT primer Tm window 'low,high' (default: 50,65)."
    ),
    skip_blast: bool = typer.Option(
        False, "--skip-blast",
        help="Skip the NCBI BLAST round; promote step1 thermal best to step2 "
             "best directly. Use only when the same target region was BLAST-"
             "validated in a prior panel."
    ),
    config_path: Optional[Path] = typer.Option(
        None, "--config", exists=True, dir_okay=False, readable=True,
        help="YAML config; takes precedence over individual flags."
    ),
    genome_path: Optional[Path] = typer.Option(
        None, "--genome", exists=True, dir_okay=False, readable=True,
        help="Local genome FASTA (default: <repo>/data/genome/GRCh38.fa)."
    ),
) -> None:
    """Run the mutation padlock probe design pipeline (4 phases)."""
    if config_path is not None:
        cfg = MutationConfig.from_yaml(config_path)
        # Allow CLI to override input/output/backbone even when YAML is used,
        # so the YAML can be a per-experiment template.
        cfg.input_xlsx = input_xlsx
        cfg.output_dir = output_dir
        cfg.backbone_file = backbone
        if genome_path is not None:
            cfg.genome_path = genome_path
        if skip_blast:
            cfg.skip_blast = True
    else:
        chemistries: List[str] = [c.strip() for c in chemistries_csv.split(",") if c.strip()]
        chem_params = default_chemistry_params()
        if ilock_tm_range:
            chem_params["iLock"].tm_range = _parse_pair(ilock_tm_range, kind=float)
        if mrna_noilock_tm_range:
            chem_params["mRNA_noiLock"].tm_range = _parse_pair(
                mrna_noilock_tm_range, kind=float,
            )
        if cdna_tm_range:
            chem_params["cDNA"].tm_range = _parse_pair(cdna_tm_range, kind=float)

        cfg = MutationConfig(
            input_xlsx=input_xlsx,
            output_dir=output_dir,
            backbone_file=backbone,
            chemistries=chemistries,
            chemistry_params=chem_params,
            arm_len_range=_parse_pair(arm_len_range, kind=int),
            mfe_min=mfe_min,
            ilock_flap=ilock_flap,
            rt_primer_gap=rt_primer_gap,
            rt_primer_tm_range=_parse_pair(rt_primer_tm_range, kind=float),
            skip_blast=skip_blast,
            start_no=start_no,
            last_no_from=last_no_from,
            genome_path=genome_path,
        )

    logger.info("Running mutation pipeline with chemistries=%s", cfg.chemistries)
    result = run_mutation_pipeline(cfg)
    typer.echo(
        f"Done. {result.mutations_mappable}/{result.mutations_total} mutations mappable; "
        f"{result.probes_total} probes (Nos {result.last_no - result.probes_total + 1}-"
        f"{result.last_no}); -> {result.final_probes_xlsx}"
    )
