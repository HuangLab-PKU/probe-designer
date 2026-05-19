"""Cross-ligation screening shim — bridges the designer to bank's screen.

The screen primitive lives in :mod:`bank.probe_book.pool.cross_ligation`
(k-mer junction seed + primer3 DNA NN model + 3'-end stability check).
We lazy-import it so the designer remains installable without bank.

Two screening modes are exposed:

* **Within-batch** (default): the candidate probes are screened against
  themselves. This is always on — every design run should be safe.
* **Pool-aware**: when ``pool_id`` is given, candidates are additionally
  screened against the existing probes in
  ``pools/<pool_id>/manifest.yaml`` (resolved via the global probe
  registry at ``probes/registry.tsv``).

The hits returned only involve at least one candidate — pool-vs-pool
pairs are filtered out (they were the prior pool author's concern, not
this design run's).

Caller responsibility: convert results into per-probe annotations,
decide whether to drop probes, write reports.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# Thresholds — match bank defaults. Exposed for plumbing through CLI.
DEFAULT_TM_THRESHOLD_C: float = 27.0
DEFAULT_END_DG_THRESHOLD_KCAL: float = -5.0
DEFAULT_MV_CONC_MM: float = 50.0
DEFAULT_DV_CONC_MM: float = 10.0
DEFAULT_DNA_CONC_M: float = 1.0e-8
DEFAULT_TEMP_C: float = 37.0


@dataclass
class CandidateProbe:
    """A new probe being designed — possibly without an assigned probe_id yet."""
    probe_id: str
    chemistry: str           # "iLock" | "dRNA" | "cDNA"
    probe_arm5: str
    probe_arm3: str
    sequence: str = ""       # full assembled probe (with flap if iLock); optional
    target: str = ""         # gene_name; informational


@dataclass
class CrossLigHit:
    probe_a_id: str
    probe_b_id: str
    kind: str                 # "Y_junc_in_X_full" | "X_junc_in_Y_full" | "both_directions"
    overall_tm_c: float
    overall_dg_kcal: float
    a_end_dg_kcal: float
    b_end_dg_kcal: float
    a_target: str = ""
    b_target: str = ""
    a_is_existing_pool: bool = False
    b_is_existing_pool: bool = False
    flagged_end_a: bool = False
    flagged_end_b: bool = False

    def as_short_str(self) -> str:
        return (f"{self.probe_b_id}(Tm={self.overall_tm_c:.1f}°C,"
                f"dG={self.overall_dg_kcal:.1f})")


def screen_candidates_against_pool(
    candidates: list[CandidateProbe],
    *,
    pool_id: Optional[str] = None,
    repo_root: Optional[Path] = None,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    end_dg_threshold_kcal: float = DEFAULT_END_DG_THRESHOLD_KCAL,
    mv_conc_mm: float = DEFAULT_MV_CONC_MM,
    dv_conc_mm: float = DEFAULT_DV_CONC_MM,
    dna_conc_m: float = DEFAULT_DNA_CONC_M,
    temp_c: float = DEFAULT_TEMP_C,
) -> tuple[list[CrossLigHit], dict[str, list[CrossLigHit]]]:
    """Run two-tier cross-lig screen on candidates ± pool.

    Args:
        candidates: new probes being designed.
        pool_id: if set, also screen against this pool's existing probes
            (loaded from ``<repo_root>/pools/<pool_id>/manifest.yaml``
            and resolved via ``<repo_root>/probes/registry.tsv``).
        repo_root: project root containing ``pools/`` and ``probes/``.
            Defaults to current working directory.
        tm_threshold_c / end_dg_threshold_kcal: Tier-2 thresholds. Match
            bank defaults; ``flagged_overall AND (flagged_end_a OR
            flagged_end_b)`` is the criterion.

    Returns:
        Tuple ``(all_hits, by_candidate)``:

        * ``all_hits`` — every confirmed cross-lig pair where at least
          one endpoint is a candidate. Pool-vs-pool hits are dropped.
        * ``by_candidate`` — ``{candidate_probe_id: [CrossLigHit, ...]}``
          for quick annotation lookup.
    """
    if not candidates:
        return [], {}

    # Lazy import bank — defer dependency until the function is called.
    try:
        from probe_book.pool.cross_ligation import (
            build_kmer_seed_index,
            _calc_dimer_thermo,
            _tier1_pair_hits,
            is_cross_ligation_risk,
        )
        from probe_book.pool.loader import EffectiveProbe
    except ImportError as exc:
        raise RuntimeError(
            "Cross-lig screen requires the `probe-book` package "
            f"(install via `pip install -e bank/`). Original error: {exc}"
        ) from exc

    # Convert candidates to EffectiveProbe.
    cand_eps: list[EffectiveProbe] = []
    cand_ids: set[str] = set()
    cand_targets: dict[str, str] = {}
    for c in candidates:
        ep = EffectiveProbe(
            probe_id=c.probe_id,
            order_id="__candidate__",
            concentration_M=0.0,
            chemistry=c.chemistry,
            probe_arm5=c.probe_arm5,
            probe_arm3=c.probe_arm3,
            sequence=c.sequence,
            target=c.target,
        )
        cand_eps.append(ep)
        cand_ids.add(c.probe_id)
        cand_targets[c.probe_id] = c.target

    # Pool probes (optional)
    pool_eps: list[EffectiveProbe] = []
    pool_targets: dict[str, str] = {}
    if pool_id is not None:
        pool_eps, pool_targets = _load_pool_effective_probes(
            pool_id, repo_root or Path.cwd(),
        )

    all_eps = cand_eps + pool_eps
    seed_index = build_kmer_seed_index(all_eps)
    targets_index: dict[str, str] = {**cand_targets, **pool_targets}

    # Tier 1 + Tier 2 pairwise — but skip pool-vs-pool.
    from itertools import combinations
    hits: list[CrossLigHit] = []
    by_cand: dict[str, list[CrossLigHit]] = {}
    for a, b in combinations(all_eps, 2):
        a_is_pool = a.probe_id not in cand_ids
        b_is_pool = b.probe_id not in cand_ids
        if a_is_pool and b_is_pool:
            continue   # caller doesn't care about pool-vs-pool here
        sa = seed_index[a.probe_id]
        sb = seed_index[b.probe_id]
        seed_hit = _tier1_pair_hits(sa, sb)
        if seed_hit is None:
            continue
        thermo = _calc_dimer_thermo(
            sa.primer3_seq, sb.primer3_seq,
            a_id=a.probe_id, b_id=b.probe_id,
            tm_threshold_c=tm_threshold_c,
            end_dg_threshold_kcal=end_dg_threshold_kcal,
            mv_conc_mm=mv_conc_mm, dv_conc_mm=dv_conc_mm,
            dna_conc_m=dna_conc_m, temp_c=temp_c,
        )
        if thermo is None or not is_cross_ligation_risk(thermo):
            continue
        hit = CrossLigHit(
            probe_a_id=a.probe_id, probe_b_id=b.probe_id,
            kind=seed_hit.kind,
            overall_tm_c=thermo.overall_tm_c,
            overall_dg_kcal=thermo.overall_dg_kcal,
            a_end_dg_kcal=thermo.a_end_dg_kcal,
            b_end_dg_kcal=thermo.b_end_dg_kcal,
            flagged_end_a=thermo.flagged_end_a,
            flagged_end_b=thermo.flagged_end_b,
            a_target=targets_index.get(a.probe_id, ""),
            b_target=targets_index.get(b.probe_id, ""),
            a_is_existing_pool=a_is_pool,
            b_is_existing_pool=b_is_pool,
        )
        hits.append(hit)
        # Always index by candidate endpoint(s)
        if not a_is_pool:
            by_cand.setdefault(a.probe_id, []).append(hit)
        if not b_is_pool:
            by_cand.setdefault(b.probe_id, []).append(hit)

    return hits, by_cand


def _load_pool_effective_probes(pool_id: str, repo_root: Path
                                  ) -> tuple[list, dict[str, str]]:
    """Resolve ``pools/<pool_id>/manifest.yaml`` to an EffectiveProbe list."""
    from probe_book.pool.manifest import PoolManifest
    from probe_book.probe.registry import ProbeRegistry
    from probe_book.pool.loader import compute_effective_probes

    manifest_path = repo_root / "pools" / pool_id / "manifest.yaml"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Pool manifest not found at {manifest_path}. "
            f"Pass --repo-root if `pools/{pool_id}/` lives elsewhere."
        )
    pool = PoolManifest.from_yaml(manifest_path)
    registry = ProbeRegistry(repo_root)
    eps = compute_effective_probes(pool, probe_registry=registry)
    # Keep only those with binding arms (skip RT primers etc.)
    eps = [e for e in eps if e.probe_arm5 and e.probe_arm3]
    targets = {e.probe_id: e.target for e in eps}
    return eps, targets


__all__ = [
    "CandidateProbe", "CrossLigHit",
    "screen_candidates_against_pool",
    "DEFAULT_TM_THRESHOLD_C", "DEFAULT_END_DG_THRESHOLD_KCAL",
    "DEFAULT_MV_CONC_MM", "DEFAULT_DV_CONC_MM",
    "DEFAULT_DNA_CONC_M", "DEFAULT_TEMP_C",
]
