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
) -> None:
    """Attach backbones to binding sites and save the assembled probes."""
    if ilock and ilock not in ("3prime", "5prime", "both"):
        raise typer.BadParameter("--ilock must be one of: 3prime, 5prime, both")

    sites = load_binding_sites(binding_sites)
    if not sites:
        typer.echo("No binding sites to assemble.", err=True)
        raise typer.Exit(code=1)

    output.mkdir(parents=True, exist_ok=True)

    assembler = ProbeAssembler(ProbeConfig(backbone_file=str(backbone)))
    probes_df = assembler.assemble_probes(sites, str(gene_info))
    probes_df["backbone_file"] = backbone.name

    if probes_df.empty:
        typer.echo("No probes were assembled (no valid gene/backbone matches).", err=True)
        raise typer.Exit(code=1)

    if ilock:
        probes_df = assembler.add_ilock_modification(probes_df, position=ilock)

    assembler.save_probes(probes_df, str(output))
    logger.info("Assembled %d probes across %d genes -> %s",
                len(probes_df), probes_df["gene_name"].nunique(), output)
    typer.echo(f"Wrote {len(probes_df)} probes to {output}")
