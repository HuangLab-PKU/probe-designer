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
        help="Comma-separated subset of {iLock, dRNA, cDNA}. "
             "Legacy 'mRNA_noiLock' is accepted and normalized to 'dRNA'."
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
             "and the `codebook` column. If omitted, inferred from the "
             "backbone filename (e.g. backbone_SP369.xlsx -> SP369)."
    ),
    ilock_tm_range: Optional[str] = typer.Option(
        None, "--ilock-tm-range",
        help="Override iLock chemistry tm_range as 'low,high' (default: 50,70)."
    ),
    drna_tm_range: Optional[str] = typer.Option(
        None, "--drna-tm-range", "--mrna-noilock-tm-range",
        help="Override dRNA (direct-RNA padlock) tm_range (default: 55,75). "
             "The legacy --mrna-noilock-tm-range alias is still accepted."
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
        help="Drop flagged probes from the final xlsx (saved to dropped_cross_lig.tsv). "
             "Requires --check-cross-lig.",
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
        if codebook is not None:
            cfg.codebook = codebook
    else:
        chemistries: List[str] = [c.strip() for c in chemistries_csv.split(",") if c.strip()]
        chem_params = default_chemistry_params()
        if ilock_tm_range:
            chem_params["iLock"].tm_range = _parse_pair(ilock_tm_range, kind=float)
        if drna_tm_range:
            chem_params["dRNA"].tm_range = _parse_pair(drna_tm_range, kind=float)
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
            codebook=codebook,
        )

    if target_pool and not check_cross_lig:
        raise typer.BadParameter("--target-pool requires --check-cross-lig")
    if reject_cross_lig and not check_cross_lig:
        raise typer.BadParameter("--reject-cross-lig requires --check-cross-lig")

    logger.info("Running mutation pipeline with chemistries=%s", cfg.chemistries)
    result = run_mutation_pipeline(cfg)
    typer.echo(
        f"Done. {result.mutations_mappable}/{result.mutations_total} mutations mappable; "
        f"{result.probes_total} probes (Nos {result.last_no - result.probes_total + 1}-"
        f"{result.last_no}); -> {result.final_probes_xlsx}"
    )

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
            from probe_book.root import find_root
            from probe_designer.ext.pool.loader import load_pool_as_probes_for_screen
            root = find_root(repo_root)
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
