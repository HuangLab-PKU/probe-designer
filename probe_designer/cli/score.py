"""``probe-design score`` — re-score and peak-rank an existing binding_sites JSON.

Useful when the user has a cached design run and wants to re-apply scoring
parameters without re-running BLAST.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import typer

from probe_designer.scoring import compute_target_score, peak_rank


logger = logging.getLogger(__name__)


def score(
    input_path: Path = typer.Option(
        ..., "--input", "-i",
        exists=True, dir_okay=False, readable=True,
        help="binding_sites JSON ({gene: [sites]} mapping)."
    ),
    output_path: Path = typer.Option(
        None, "--output", "-o",
        help="Output JSON path. Defaults to <input>.scored.json"
    ),
    min_arm_tm: float = typer.Option(50.0, "--min-arm-tm", help="Tm proximity target for scoring."),
    max_tm_diff: float = typer.Option(10.0, "--max-tm-diff", help="Max Tm diff for balance scoring."),
    region_size: int = typer.Option(80, "--region-size", help="bp window for peak_rank grouping."),
    min_gap: int = typer.Option(40, "--min-gap", help="Minimum gap within a peak_rank region."),
) -> None:
    """Compute per-site target scores and peak-rank them within each gene."""
    data = json.loads(input_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise typer.BadParameter("Input JSON must be a {gene: [sites]} mapping.")

    scored: dict[str, list] = {}
    total_sites = 0
    for gene, sites in data.items():
        for site in sites:
            total_isoforms = site.get("_total_isoforms", 1)
            site["score"] = compute_target_score(
                site,
                min_arm_tm=min_arm_tm,
                max_tm_diff=max_tm_diff,
                total_isoforms=total_isoforms,
            )
        ranked = peak_rank(sites, region_size=region_size, min_gap=min_gap)
        # Attach peak_rank index for downstream consumers
        for idx, site in enumerate(ranked):
            site["peak_rank"] = idx
        scored[gene] = ranked
        total_sites += len(ranked)

    logger.info("Scored and peak-ranked %d sites across %d genes.", total_sites, len(scored))

    out = output_path or input_path.with_name(input_path.stem + ".scored.json")
    out.write_text(json.dumps(scored, indent=2, ensure_ascii=False), encoding="utf-8")
    typer.echo(f"Wrote {out}")
