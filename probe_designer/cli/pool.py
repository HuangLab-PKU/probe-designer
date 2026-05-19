"""``probe-design pool`` — pool-level operations (currently: cross-lig audit).

This is the periodic / standalone entry point for cross-ligation
screening of an existing pool. Use the design CLIs' ``--check-cross-lig
--target-pool <pool_id>`` flow for new probes; use this command for
auditing a pool that has already shipped.

Lazy-imports ``probe_book`` so the designer remains installable without
the bank package.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import typer


logger = logging.getLogger(__name__)


pool_app = typer.Typer(
    name="pool",
    help="Pool-level operations (cross-ligation audit).",
    no_args_is_help=True,
)


@pool_app.command("check")
def check(
    pool_id: str = typer.Argument(
        ..., help="Pool ID, e.g. `TNBCmarker_main_v1`."
    ),
    repo_root: Optional[Path] = typer.Option(
        None, "--repo-root",
        help="Project root containing `pools/` and `probes/`. Defaults to cwd."
    ),
    xlig_tm_threshold: float = typer.Option(
        27.0, "--xlig-tm-threshold",
        help="Minimum overall dimer Tm (°C) to flag a pair.",
    ),
    xlig_end_dg_threshold: float = typer.Option(
        -5.0, "--xlig-end-dg-threshold",
        help="Maximum end-stability ΔG (kcal/mol) for 3'-end-annealed flag.",
    ),
) -> None:
    """Run the full Tier-1 + Tier-2 cross-lig screen on an existing pool.

    Writes ``tier1_hits.tsv`` and ``confirmed_hits.tsv`` under
    ``<repo_root>/pools/<pool_id>/pool_check/``. Use for periodic
    re-audits or to re-screen legacy pools with the primer3 DNA backend.
    """
    try:
        from probe_book.pool.manifest import PoolManifest
        from probe_book.probe.registry import ProbeRegistry
        from probe_book.pool.loader import compute_effective_probes
        from probe_book.pool.cross_ligation import screen_pool_cross_ligation
        from probe_book.pool.report import write_tier1_hits, write_confirmed_hits
    except ImportError as exc:
        raise typer.BadParameter(
            f"pool check requires the `probe-book` package: {exc}"
        ) from exc

    root = (repo_root or Path.cwd()).resolve()
    manifest_path = root / "pools" / pool_id / "manifest.yaml"
    if not manifest_path.exists():
        raise typer.BadParameter(f"Pool manifest not found at {manifest_path}")

    pool = PoolManifest.from_yaml(manifest_path)
    registry = ProbeRegistry(root)
    probes = compute_effective_probes(pool, probe_registry=registry)
    binding = [p for p in probes if p.probe_arm5 and p.probe_arm3]
    logger.info("Pool %s — %d effective probes (%d binding-arm).",
                pool.pool_id, len(probes), len(binding))

    tier1, confirmed = screen_pool_cross_ligation(
        binding,
        tm_threshold_c=xlig_tm_threshold,
        end_dg_threshold_kcal=xlig_end_dg_threshold,
    )

    out_dir = root / "pools" / pool_id / "pool_check"
    out_dir.mkdir(parents=True, exist_ok=True)
    write_tier1_hits(tier1, out_dir / "tier1_hits.tsv")
    write_confirmed_hits(confirmed, out_dir / "confirmed_hits.tsv")

    typer.echo(
        f"pool={pool_id}: probes={len(binding)} pairs={len(binding)*(len(binding)-1)//2} "
        f"tier1={len(tier1)} confirmed_xlig={len(confirmed)} -> {out_dir}"
    )
