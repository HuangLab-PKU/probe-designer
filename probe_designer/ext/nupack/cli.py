"""Optional `probe-design nupack` Typer subcommand for standalone ternary
analysis runs (without going through `probe-design pool check`)."""
from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Optional

try:
    import typer
except ImportError:    # pragma: no cover - typer is a hard dep of the CLI
    typer = None


from probe_designer.ext.nupack.check import screen_ternary, write_ternary_report
from probe_designer.ext.nupack.config import (
    DEFAULT_CELSIUS,
    DEFAULT_ENSEMBLE,
    DEFAULT_MAGNESIUM_M,
    DEFAULT_SODIUM_M,
    DEFAULT_STRAND_CONC_M,
    DEFAULT_VICINITY_N,
)


nupack_app = typer.Typer(
    help="NUPACK 4 ternary-complex cross-lig confirmation (v2.3).",
    no_args_is_help=True,
) if typer is not None else None


def _run_ternary(
    pool_id: str,
    repo_root: Optional[Path] = None,
    mode: str = "tier1",
    sodium_m: float = DEFAULT_SODIUM_M,
    magnesium_m: float = DEFAULT_MAGNESIUM_M,
    celsius: float = DEFAULT_CELSIUS,
    ensemble: str = DEFAULT_ENSEMBLE,
    strand_conc_m: float = DEFAULT_STRAND_CONC_M,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
    out_dir: Optional[Path] = None,
) -> None:
    """Standalone ternary run on a bank-tracked pool."""
    from probe_book.root import resolve_root
    from probe_designer.ext.pool.loader import load_pool_as_probes_for_screen
    from probe_designer.qc.cross_ligation import screen_cross_ligation

    repo_root = resolve_root(repo_root)
    out_dir = out_dir or (repo_root / "pools" / pool_id / "pool_check")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading pool {pool_id} from {repo_root}")
    probes = load_pool_as_probes_for_screen(pool_id, repo_root)
    print(f"  {len(probes)} probes loaded")

    # Build the prefilter pair list. The register scan is the gate; NUPACK only
    # re-scores what it already flagged, so the modes differ by how much of the
    # scan's output to carry forward.
    #   all-registers — every ligation-competent register, Tm-clearing or not
    #   confirmed     — only those over the Tm threshold (the default)
    #   all-pairs     — skip the gate entirely; O(N^2) NUPACK calls
    print(f"Running prefilter (mode={mode}) ...")
    aliases = {"tier1": "all-registers", "v22-nick": "confirmed"}
    resolved = aliases.get(mode, mode)
    if resolved == "all-pairs":
        prefilter_pairs = None
    elif resolved in ("all-registers", "confirmed"):
        registers = screen_cross_ligation(probes)
        if resolved == "confirmed":
            registers = [d for d in registers if d.flagged_overall]
        prefilter_pairs = [(d.seq_a_id, d.seq_b_id) for d in registers]
    else:
        raise ValueError(
            f"unknown --nupack-mode {mode!r}; expected "
            f"all-registers | confirmed | all-pairs"
        )
    print(f"  prefilter pairs: {'(all-pairs)' if prefilter_pairs is None else len(prefilter_pairs)}")

    t0 = time.time()
    last_print = [t0]

    def _progress(idx: int, total: int) -> None:
        now = time.time()
        if now - last_print[0] > 30 or idx == total:
            elapsed = now - t0
            rate = idx / elapsed if elapsed > 0 else 0
            eta = (total - idx) / rate if rate > 0 else 0
            print(f"  [{idx}/{total}] {rate:.1f} pair/s   ETA {eta / 60:.1f} min")
            last_print[0] = now

    hits = screen_ternary(
        probes, prefilter_pairs=prefilter_pairs,
        sodium_m=sodium_m, magnesium_m=magnesium_m, celsius=celsius,
        ensemble=ensemble, strand_conc_m=strand_conc_m,
        vicinity_n_each_side=vicinity_n_each_side,
        on_progress=_progress,
    )
    print(f"\nNUPACK ternary run complete in {(time.time() - t0)/60:.1f} min")
    print(f"  total ternary hits: {len(hits)}")

    confirmed = [
        h for h in hits
        if h.mfe_nick_adjacent and h.mfe_vicinity_contiguous
        and h.ternary_fraction_of_b > 1e-4
    ]
    print(f"  NUPACK-confirmed (MFE nick-adjacent + vicinity + fraction > 1e-4): {len(confirmed)}")

    write_ternary_report(hits, out_dir / "nupack_ternary_dimers.tsv")
    write_ternary_report(confirmed, out_dir / "nupack_confirmed.tsv")
    print(f"  wrote {out_dir}/nupack_ternary_dimers.tsv ({len(hits)} rows)")
    print(f"  wrote {out_dir}/nupack_confirmed.tsv ({len(confirmed)} rows)")


if typer is not None:
    @nupack_app.command("ternary")
    def ternary(
        pool_id: str = typer.Argument(..., help="Pool ID to screen"),
        repo_root: Optional[Path] = typer.Option(None, "--repo-root", help="Project root (auto-detect if omitted)"),
        mode: str = typer.Option("tier1", "--mode", help="tier1 | v22-nick | all-pairs"),
        sodium_m: float = typer.Option(DEFAULT_SODIUM_M, "--sodium-m"),
        magnesium_m: float = typer.Option(DEFAULT_MAGNESIUM_M, "--magnesium-m"),
        celsius: float = typer.Option(DEFAULT_CELSIUS, "--celsius"),
        ensemble: str = typer.Option(DEFAULT_ENSEMBLE, "--ensemble", help="stacking | nostacking"),
        strand_conc_m: float = typer.Option(DEFAULT_STRAND_CONC_M, "--strand-conc"),
        vicinity_n: int = typer.Option(DEFAULT_VICINITY_N, "--vicinity-n"),
        out_dir: Optional[Path] = typer.Option(None, "--out-dir"),
    ) -> None:
        """Run NUPACK 4 ternary-complex screen on a pool (Tier 3 v2.3)."""
        _run_ternary(
            pool_id=pool_id, repo_root=repo_root, mode=mode,
            sodium_m=sodium_m, magnesium_m=magnesium_m, celsius=celsius,
            ensemble=ensemble, strand_conc_m=strand_conc_m,
            vicinity_n_each_side=vicinity_n, out_dir=out_dir,
        )


__all__ = ["nupack_app", "_run_ternary"]
