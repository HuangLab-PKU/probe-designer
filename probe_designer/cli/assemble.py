"""``probe-design assemble`` — attach backbones to binding sites.

Thin Typer wrapper around :class:`probe_designer.probe_assembly.ProbeAssembler`.
All loader / format-dispatch logic lives in
``probe_designer.probe_assembly.load_binding_sites``; the CLI only wires
flags + calls the library.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import typer

from probe_designer.config import ProbeConfig
from probe_designer.probe_assembly import ProbeAssembler, load_binding_sites


logger = logging.getLogger(__name__)


def assemble(
    binding_sites: Path = typer.Option(
        ..., "--binding-sites", exists=True, dir_okay=False, readable=True,
        help="Binding sites JSON or XLSX (from probe-design design / score)."
    ),
    gene_info: Path = typer.Option(
        ..., "--gene-info", exists=True, dir_okay=False, readable=True,
        help="Gene info XLSX (must have 'gene_name' and 'No.' columns)."
    ),
    backbone: Path = typer.Option(
        ..., "--backbone", exists=True, dir_okay=False, readable=True,
        help="Backbone XLSX (must have 'No.' and 'Sequence' columns)."
    ),
    output: Path = typer.Option(
        ..., "--output", "-o",
        help="Output directory for assembled probe files."
    ),
    ilock: Optional[str] = typer.Option(
        None, "--ilock", case_sensitive=False,
        help="Add iLock modification on 3prime / 5prime / both arms."
    ),
    codebook: Optional[str] = typer.Option(
        None, "--codebook",
        help="Codebook name (e.g. SP369). Goes into every padlock probe_name "
             "and the `codebook` column. If omitted, inferred from the backbone "
             "filename (e.g. backbone_SP369.xlsx -> SP369)."
    ),
    check_cross_lig: bool = typer.Option(
        False, "--check-cross-lig",
        help="Run primer3 DNA cross-ligation screen on the assembled probes. "
             "Always screens within-batch; pair with --target-pool to also "
             "screen against an existing pool. Writes `cross_lig_report.tsv` "
             "and adds a `cross_lig_partners` column to the output xlsx.",
    ),
    target_pool: Optional[str] = typer.Option(
        None, "--target-pool",
        help="Pool ID (e.g. `TNBCmarker_main_v1`) whose existing probes the "
             "new candidates will be screened against. Requires --check-cross-lig. "
             "Reads from `pools/<pool_id>/manifest.yaml` and `probes/registry.tsv` "
             "under the current working directory (or --repo-root).",
    ),
    repo_root: Optional[Path] = typer.Option(
        None, "--repo-root",
        help="Project root containing pools/ and probes/. Defaults to cwd.",
    ),
    reject_cross_lig: bool = typer.Option(
        False, "--reject-cross-lig",
        help="Drop flagged probes from the assembled output. Dropped rows are "
             "saved to `dropped_cross_lig.tsv`. Requires --check-cross-lig.",
    ),
    xlig_tm_threshold: float = typer.Option(
        27.0, "--xlig-tm-threshold",
        help="Minimum overall dimer Tm (°C) to flag a pair (primer3 default).",
    ),
    xlig_end_dg_threshold: Optional[float] = typer.Option(
        None, "--xlig-end-dg-threshold",
        help="DEPRECATED — v2 uses junction-geometry instead of end-stability. "
             "Accepted as a no-op; will be removed in a future release.",
    ),
) -> None:
    """Attach backbones to binding sites and save the assembled probes."""
    if ilock and ilock not in ("3prime", "5prime", "both"):
        raise typer.BadParameter("--ilock must be one of: 3prime, 5prime, both")

    sites = load_binding_sites(binding_sites)
    if not sites:
        typer.echo("No binding sites to assemble.", err=True)
        raise typer.Exit(code=1)

    output.mkdir(parents=True, exist_ok=True)

    assembler = ProbeAssembler(ProbeConfig(
        backbone_file=str(backbone), codebook=codebook,
    ))
    probes_df = assembler.assemble_probes(sites, str(gene_info))
    probes_df["backbone_file"] = backbone.name

    if probes_df.empty:
        typer.echo("No probes were assembled (no valid gene/backbone matches).", err=True)
        raise typer.Exit(code=1)

    if ilock:
        probes_df = assembler.add_ilock_modification(probes_df, position=ilock)

    # Cross-ligation screen (opt-in via --check-cross-lig).
    if target_pool and not check_cross_lig:
        raise typer.BadParameter("--target-pool requires --check-cross-lig")
    if reject_cross_lig and not check_cross_lig:
        raise typer.BadParameter("--reject-cross-lig requires --check-cross-lig")

    if xlig_end_dg_threshold is not None:
        import warnings
        warnings.warn(
            "--xlig-end-dg-threshold is deprecated in v2 (junction-geometry "
            "criterion replaces end-stability). Flag accepted as no-op.",
            DeprecationWarning, stacklevel=2,
        )

    if check_cross_lig:
        from probe_designer.qc.assemble_hook import apply_cross_lig_check

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

        annotated, dropped, report = apply_cross_lig_check(
            probes_df,
            splint_probes=splint_probes,
            reject=reject_cross_lig,
            tm_threshold_c=xlig_tm_threshold,
        )
        output.mkdir(parents=True, exist_ok=True)
        report.to_csv(output / "cross_lig_report.tsv", sep="\t", index=False)
        n_hits = len(report)
        n_dropped = len(dropped)
        if n_dropped:
            dropped.to_csv(output / "dropped_cross_lig.tsv", sep="\t", index=False)
        logger.info(
            "Cross-lig screen v2: %d confirmed pair(s); %d probe(s) dropped; pool=%s",
            n_hits, n_dropped, target_pool or "(within-batch only)",
        )
        typer.echo(
            f"Cross-lig: {n_hits} pair(s) flagged; "
            f"{n_dropped} probe(s) dropped"
            f" (report: {output / 'cross_lig_report.tsv'})"
        )
        probes_df = annotated

    assembler.save_probes(probes_df, str(output))
    logger.info("Assembled %d probes across %d genes -> %s",
                len(probes_df), probes_df["gene_name"].nunique(), output)
    typer.echo(f"Wrote {len(probes_df)} probes to {output}")
