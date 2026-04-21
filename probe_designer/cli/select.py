"""``probe-design select`` — apply top-N + min-gap selection to scored sites JSON.

Mirrors the webapp's post-design auto-selection step (task_runner._auto_select_top_n)
so CLI users get the same spatial-diversity-aware picks.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import typer

from probe_designer.scoring import select_top_n_with_gap


logger = logging.getLogger(__name__)


def select(
    input_path: Path = typer.Option(
        ..., "--input", "-i",
        exists=True, dir_okay=False, readable=True,
        help="Scored binding_sites JSON (per-gene lists of site dicts)."
    ),
    output_path: Path = typer.Option(
        None, "--output", "-o",
        help="Output JSON path. Defaults to <input>.selected.json"
    ),
    top_n: int = typer.Option(3, "--top-n", "-n", min=1, help="Max sites per gene to select."),
    min_gap: int = typer.Option(40, "--min-gap", help="Minimum bp gap between selected site midpoints."),
) -> None:
    """Select the top-N highest-scoring sites per gene with a spatial spacing constraint."""
    data = json.loads(input_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise typer.BadParameter("Input JSON must be a {gene: [sites]} mapping.")

    selected_by_gene = {
        gene: select_top_n_with_gap(sites, top_n=top_n, min_gap=min_gap)
        for gene, sites in data.items()
    }

    total_in = sum(len(v) for v in data.values())
    total_out = sum(len(v) for v in selected_by_gene.values())
    logger.info(
        "Selected %d of %d sites across %d genes (top_n=%d, min_gap=%d).",
        total_out, total_in, len(selected_by_gene), top_n, min_gap,
    )

    out = output_path or input_path.with_name(input_path.stem + ".selected.json")
    out.write_text(
        json.dumps(selected_by_gene, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    typer.echo(f"Wrote {out}")
