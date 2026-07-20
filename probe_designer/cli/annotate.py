"""``probe-design annotate`` — write reusable thermodynamic annotation tracks.

For every record in a reference FASTA, emit per-position bedGraph tracks (arm Tm
+ RNAplfold accessibility) computed at the given hybridization conditions, to a
directory (e.g. a NAS path beside the genome). IGV/UCSC read bedGraph directly.
"""
from __future__ import annotations

import logging
from pathlib import Path

import typer
from Bio import SeqIO

from probe_designer.annotate import build_reference_annotations
from probe_designer.chemistry import ReactionConditions

logger = logging.getLogger(__name__)


def annotate(
    fasta: Path = typer.Option(
        ..., "--fasta", exists=True, dir_okay=False, readable=True,
        help="Reference FASTA (transcript/mRNA sequences, 5'->3'); one track set per record.",
    ),
    out_dir: Path = typer.Option(
        ..., "--out-dir", "-o",
        help="Directory for the .bedgraph tracks (e.g. a NAS path beside the genome).",
    ),
    arm_length: int = typer.Option(20, "--arm-length", min=1, help="Arm/window length for the Tm track."),
    chemistry: str = typer.Option("dRNA", "--chemistry", help="dRNA/iLock (DNA:RNA) or cDNA (DNA:DNA)."),
    no_accessibility: bool = typer.Option(
        False, "--no-accessibility", help="Skip the RNAplfold accessibility track (Tm only)."
    ),
    # --- current hybridization conditions (defaults = protocol rca.md v5.3) ---
    monovalent_mm: float = typer.Option(75.0, "--monovalent-mm", help="Monovalent K+ (mM)."),
    mg_mm: float = typer.Option(10.0, "--mg-mm", help="Mg2+ (mM)."),
    dntp_mm: float = typer.Option(0.0, "--dntp-mm", help="dNTP (mM)."),
    strand_nm: float = typer.Option(100.0, "--strand-nm", help="Probe/oligo conc (nM)."),
    formamide_pct: float = typer.Option(20.0, "--formamide-pct", help="Formamide (percent v/v)."),
    lab_temp_c: float = typer.Option(45.0, "--lab-temp-c", help="Hyb anneal temperature (C)."),
    saltcorr: int = typer.Option(5, "--saltcorr", min=1, max=7, help="Biopython salt method."),
) -> None:
    """Write per-position thermodynamic annotation tracks for each FASTA record."""
    reaction = ReactionConditions(
        monovalent_mM=monovalent_mm, mg_mM=mg_mm, dntp_mM=dntp_mm,
        strand_nM=strand_nm, formamide_pct=formamide_pct, lab_temp_c=lab_temp_c,
        saltcorr=saltcorr,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for rec in SeqIO.parse(str(fasta), "fasta"):
        paths = build_reference_annotations(
            str(rec.seq), rec.id, reaction,
            out_dir=out_dir, arm_length=arm_length, chemistry=chemistry,
            accessibility=not no_accessibility,
        )
        typer.echo(f"{rec.id}: {len(paths)} track(s)")
        n += 1
    typer.echo(f"Done: {n} reference(s) at buffer {reaction.signature()} -> {out_dir}")
