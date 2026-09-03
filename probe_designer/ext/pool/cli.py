"""``probe-design pool`` — Typer subcommand for pool-level operations.

Lives under ``probe_designer.ext.pool`` because it depends on probe-book
for pool data loading. Registered into the main ``probe-design`` app via
``app.add_typer(pool_app, name="pool")`` at ``cli/app.py``.
"""
from __future__ import annotations

import csv
import logging
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Set

import typer

from probe_designer.qc.cross_ligation import (
    DEFAULT_TM_THRESHOLD_C,
    screen_cross_ligation_v2,
    write_dimer_report,
)


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
        DEFAULT_TM_THRESHOLD_C, "--xlig-tm-threshold",
        help="Minimum overall dimer Tm (°C) to flag a pair.",
    ),
    xlig_end_dg_threshold: Optional[float] = typer.Option(
        None, "--xlig-end-dg-threshold",
        help="DEPRECATED — v2 uses junction-geometry instead of end-stability. "
             "Accepted as a no-op; will be removed in a future release.",
    ),
    nupack: bool = typer.Option(
        False, "--nupack",
        help="Run NUPACK 4 ternary-complex Tier 3 confirmation (v2.3) on top "
             "of v2.2 (slow: 3-h overnight on TNBC-size pool).",
    ),
    nupack_mode: str = typer.Option(
        "tier1", "--nupack-mode",
        help="Prefilter mode for NUPACK Tier 3 — tier1 (default, max sensitivity), "
             "v22-nick (only v2.2 nick-adjacent hits), or all-pairs.",
    ),
    nupack_ensemble: str = typer.Option(
        "stacking", "--nupack-ensemble",
        help="NUPACK Model ensemble — stacking (default, coaxial stacking on) "
             "or nostacking (ablation test).",
    ),
) -> None:
    """Run the v2 cross-lig screen on an existing pool.

    Loads the pool via probe-book, runs the asymmetric v2 algorithm
    (correct ligation geometry + ASCII-parse junction check), and writes
    ``confirmed_dimers_v2.tsv`` under ``<repo_root>/pools/<pool_id>/pool_check/``.
    """
    if xlig_end_dg_threshold is not None:
        warnings.warn(
            "--xlig-end-dg-threshold is deprecated and ignored in v2 "
            "(end-stability has been replaced by junction-geometry).",
            DeprecationWarning,
            stacklevel=2,
        )

    try:
        from probe_book.root import resolve_root
        from probe_designer.ext.pool.loader import load_pool_as_probes_for_screen
    except ImportError as exc:
        raise typer.BadParameter(str(exc)) from exc

    root = resolve_root(repo_root)
    try:
        probes = load_pool_as_probes_for_screen(pool_id, root)
    except (ImportError, FileNotFoundError) as exc:
        raise typer.BadParameter(str(exc)) from exc

    logger.info("Pool %s — %d binding-arm probes loaded.", pool_id, len(probes))

    tier1, dimers = screen_cross_ligation_v2(
        probes,
        tm_threshold_c=xlig_tm_threshold,
    )
    confirmed = [d for d in dimers if d.flagged_overall and d.a_can_ligate_on_b]

    out_dir = root / "pools" / pool_id / "pool_check"
    out_dir.mkdir(parents=True, exist_ok=True)
    write_dimer_report(confirmed, out_dir / "confirmed_dimers_v2.tsv")
    write_dimer_report(dimers, out_dir / "all_dimers_v2.tsv")

    n_pairs = len(probes) * (len(probes) - 1) // 2
    typer.echo(
        f"pool={pool_id}: probes={len(probes)} pairs={n_pairs} "
        f"tier1_hits={len(tier1)} primer3_dimers={len(dimers)} "
        f"confirmed_v2={len(confirmed)} -> {out_dir}"
    )

    if nupack:
        # Lazy import — fails cleanly if NUPACK 4 isn't installed
        try:
            from probe_designer.ext.nupack.cli import _run_ternary
        except ImportError as exc:
            raise typer.BadParameter(
                f"--nupack requires NUPACK 4: {exc}"
            ) from exc
        typer.echo("\nRunning NUPACK Tier 3 …")
        _run_ternary(
            pool_id=pool_id, repo_root=root, mode=nupack_mode,
            ensemble=nupack_ensemble,
            out_dir=out_dir,
        )


@pool_app.command("create-design")
def create_design(
    pool_id: str = typer.Argument(
        ..., help="New pool id, e.g. `TNBCmarker_cand_v1` (<aim>_<purpose>_v<rev>)."
    ),
    probe: List[str] = typer.Option(
        [], "--probe", help="A probe_id to include (repeatable)."
    ),
    probes_file: Optional[Path] = typer.Option(
        None, "--probes-file",
        help="File with one probe_id per line (a leading 'probe_id' header is skipped).",
    ),
    pool_name: Optional[str] = typer.Option(None, "--pool-name"),
    codebook: str = typer.Option("SP369", "--codebook"),
    well_volume_uL: float = typer.Option(10.0, "--well-volume"),
    repo_root: Optional[Path] = typer.Option(None, "--repo-root"),
    force: bool = typer.Option(False, "--force", help="Overwrite an existing manifest."),
) -> None:
    """Create a DESIGN pool (candidate composition) from registry probe_ids.

    Design pools carry no order binding — they exist to validate a composition
    (`probe-design pool check`) BEFORE the probes are ordered. Once the probes
    are ordered, `pool promote` binds them to orders and flips the pool to
    physical (which the pool-mixing-guide can then pipette).
    """
    try:
        from probe_book.pool import PoolKind, PoolManifest
        from probe_book.probe import ProbeRegistry
        from probe_book.root import resolve_root
    except ImportError as exc:
        raise typer.BadParameter(str(exc)) from exc

    root = resolve_root(repo_root)

    probe_ids: List[str] = list(probe)
    if probes_file:
        for line in Path(probes_file).read_text(encoding="utf-8").splitlines():
            s = line.strip().split("\t")[0].strip()
            if s and s.lower() != "probe_id":
                probe_ids.append(s)
    seen: Set[str] = set()
    probe_ids = [p for p in probe_ids if not (p in seen or seen.add(p))]
    if not probe_ids:
        raise typer.BadParameter("no probe_ids given (use --probe and/or --probes-file)")

    known = {r.get("probe_id") for r in ProbeRegistry(root).read_all()}
    unknown = [p for p in probe_ids if p not in known]
    if unknown:
        raise typer.BadParameter(
            f"{len(unknown)} probe_id(s) not in probes/registry.tsv: {unknown[:8]}"
        )

    dest = root / "pools" / pool_id / "manifest.yaml"
    if dest.exists() and not force:
        raise typer.BadParameter(f"{dest} already exists (pass --force to overwrite)")

    pm = PoolManifest(
        pool_id=pool_id, pool_name=pool_name or f"{pool_id} (design)",
        codebook=codebook, target_well_volume_uL=well_volume_uL,
        kind=PoolKind.DESIGN.value,
    )
    for pid in probe_ids:
        pm.add_probe(probe_id=pid)
    pm.to_yaml(dest)
    typer.echo(
        f"created design pool {pool_id}: {len(probe_ids)} probes -> {dest}\n"
        f"validate with:  probe-design pool check {pool_id}"
    )


@pool_app.command("promote")
def promote(
    pool_id: str = typer.Argument(..., help="Design pool id to promote to physical."),
    order: Optional[str] = typer.Option(
        None, "--order", help="Bind ALL probes to this single order_id."
    ),
    bindings_file: Optional[Path] = typer.Option(
        None, "--bindings-file", help="TSV with columns probe_id, order_id."
    ),
    auto: bool = typer.Option(
        False, "--auto",
        help="Resolve each probe's order from orders/*/inventory.tsv (default when "
             "neither --order nor --bindings-file is given; errors if a probe is in "
             "0 or >1 orders).",
    ),
    concentration_M: float = typer.Option(
        1e-8, "--concentration", help="Concentration (M) applied to every probe."
    ),
    repo_root: Optional[Path] = typer.Option(None, "--repo-root"),
) -> None:
    """Promote a DESIGN pool to PHYSICAL by binding each probe to its source order.

    Binding source (pick one): ``--bindings-file`` (per-probe), ``--order`` (all
    from one order), or auto-resolution from the order inventories (the default).
    A physical pool feeds the pool-mixing-guide.
    """
    try:
        from probe_book.pool import PoolKind, PoolManifest
        from probe_book.root import resolve_root
    except ImportError as exc:
        raise typer.BadParameter(str(exc)) from exc

    root = resolve_root(repo_root)
    mpath = root / "pools" / pool_id / "manifest.yaml"
    if not mpath.exists():
        raise typer.BadParameter(f"pool manifest not found: {mpath}")
    pm = PoolManifest.from_yaml(mpath)
    if pm.kind != PoolKind.DESIGN.value:
        raise typer.BadParameter(
            f"pool {pool_id} is not a design pool (kind={pm.kind}); nothing to promote"
        )

    probe_ids = pm.probe_ids()
    bindings: Dict[str, str] = {}
    if bindings_file:
        with open(bindings_file, encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                bindings[r["probe_id"].strip()] = r["order_id"].strip()
    elif order:
        bindings = {p: order for p in probe_ids}
    else:  # auto-resolve from order inventories
        prov: Dict[str, Set[str]] = {}
        for inv in (root / "orders").glob("*/inventory.tsv"):
            oid = inv.parent.name
            with open(inv, encoding="utf-8") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    prov.setdefault(r["probe_id"].strip(), set()).add(oid)
        missing, ambiguous = [], []
        for p in probe_ids:
            orders = prov.get(p, set())
            if len(orders) == 1:
                bindings[p] = next(iter(orders))
            elif not orders:
                missing.append(p)
            else:
                ambiguous.append((p, sorted(orders)))
        if missing or ambiguous:
            msg = []
            if missing:
                msg.append(f"not ordered (in no inventory): {missing}")
            if ambiguous:
                msg.append(f"in multiple orders — use --order or --bindings-file: {ambiguous}")
            raise typer.BadParameter("auto-resolve failed:\n  " + "\n  ".join(msg))

    try:
        pm.promote(bindings=bindings, default_concentration_M=concentration_M)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    pm.to_yaml(mpath)
    n_orders = len({e.order_id for e in pm.recipe})
    typer.echo(
        f"promoted {pool_id} -> physical: {len(pm.recipe)} probes across "
        f"{n_orders} order(s)\n"
        f"mixing guide:  python .claude/skills/pool-mixing-guide/scripts/mixing_guide.py "
        f"--pool {pool_id}"
    )
