"""Target scoring, peak-finding ranking, and top-N selection for padlock probe design.

Public API (preserved from previous flat scoring.py):
- compute_target_score(): composite quality score for a binding site
- peak_rank(): position-aware ranking that spreads top targets across the gene

Added in Phase 1:
- select_top_n_with_gap(): DB-agnostic top-N selection with minimum spatial gap
  between selections (migrated from webapp task_runner._auto_select_top_n).
"""

from typing import Any, Dict, List, Optional

from ..policy import DEFAULT_POLICY, ThermoPolicy
from .scorer import compute_off_target_score
from .selection import select_score_peaks, select_top_n_with_gap


__all__ = [
    "compute_target_score",
    "compute_off_target_score",
    "peak_rank",
    "select_score_peaks",
    "select_top_n_with_gap",
]


def compute_target_score(
    site: Dict[str, Any],
    policy: Optional[ThermoPolicy] = None,
    total_isoforms: int = 1,
) -> float:
    """Compute composite quality score for a target region (0-10 scale).

    Components (higher = better) and their weights:
    1. isoform_coverage (0-3): fraction of isoforms covered
    2. blast_support (0-2): same-gene mRNA hit count (after specificity filter)
    3. tm_proximity (0-2): arm Tm near the policy target
    4. tm_balance (0-1): smaller |tm_5prime - tm_3prime|
    5. terminal_gc (0-1): G/C at padlock ligation point (arm5 5'-end, arm3 3'-end)
    6. delta_g (0-1): lower MFE = more stable = better

    The weights live here; the Tm *thresholds and shapes* live on
    :class:`~probe_designer.policy.ThermoPolicy` so the ranker and the filter
    cannot disagree about them (audit R7). Pre-2026-07-21 this took loose
    ``min_arm_tm`` / ``max_tm_diff`` floats that defaulted independently of
    FilterConfig — pass ``ThermoPolicy.resolve(config)`` instead.
    """
    policy = policy or DEFAULT_POLICY
    score = 0.0

    overlap = site.get("isoform_overlap_num", 1)
    if total_isoforms > 0:
        score += min(3.0, (overlap / total_isoforms) * 3.0)
    else:
        score += 1.5  # neutral if no isoform data

    alignments = site.get("blast_alignments", [])
    mrna_hits = sum(
        1 for a in alignments
        if "mrna" in a.get("hit_def", "").lower()
        or "nm_" in a.get("hit_id", "").lower()
        or "xm_" in a.get("hit_id", "").lower()
    )
    score += min(2.0, mrna_hits / 5.0 * 2.0)

    tm5 = site.get("tm_5prime", 0)
    tm3 = site.get("tm_3prime", 0)
    if tm5 > 0 and tm3 > 0:
        score += 2.0 * policy.tm_range_score((tm5 + tm3) / 2)

    score += 1.0 * policy.tm_balance_score(abs(tm5 - tm3))

    arm5 = site.get("arm_5prime", "")
    arm3 = site.get("arm_3prime", "")
    if arm5 and arm5[0] in "GCgc":
        score += 0.5
    if arm3 and arm3[-1] in "GCgc":
        score += 0.5

    # NOTE: docstring claims 'dG=0 -> 1.0 (best)' but the guard is `if mfe < 0`
    # so mfe=0 actually contributes 0. Phase 1 characterization tests lock
    # this behavior; behavior change is a separate decision.
    mfe = site.get("free_energy", 0)
    if mfe < 0:
        score += max(0.0, 1.0 - abs(mfe) / 10.0)

    return round(score, 3)


def peak_rank(
    sites: List[Dict[str, Any]],
    region_size: int = 80,
    min_gap: int = 40,
) -> List[Dict[str, Any]]:
    """Peak-finding ranking: spread top targets across the gene.

    Algorithm:
    1. Group sites into regions (region_size bp each)
    2. Round 1: pick the highest-score "peak" from each region
    3. Round 2: pick next-best from each region (>= min_gap from peak)
    4. Continue until all sites are ranked
    5. Within each round, sort by score descending

    Returns: re-ordered sites list (spatially distributed + high score first)
    """
    if not sites:
        return []

    for site in sites:
        site["_region"] = site.get("st", 0) // region_size

    regions: Dict[int, List[Dict]] = {}
    for site in sites:
        regions.setdefault(site["_region"], []).append(site)
    for region_sites in regions.values():
        region_sites.sort(key=lambda s: s.get("score", 0), reverse=True)

    result = []
    while any(len(v) > 0 for v in regions.values()):
        round_sites = []
        for region_id in sorted(regions.keys()):
            candidates = regions[region_id]
            if not candidates:
                continue

            selected_in_region = [
                s.get("st", 0) for s in result if s.get("_region") == region_id
            ]
            picked = False
            for i, c in enumerate(candidates):
                pos = c.get("st", 0)
                if not any(abs(pos - sp) < min_gap for sp in selected_in_region):
                    round_sites.append(candidates.pop(i))
                    picked = True
                    break

            if not picked:
                continue

        if not round_sites:
            remaining = True
            while remaining:
                remaining = False
                tail_round = []
                for rid in sorted(regions.keys()):
                    if regions[rid]:
                        tail_round.append(regions[rid].pop(0))
                        remaining = True
                tail_round.sort(key=lambda s: s.get("score", 0), reverse=True)
                result.extend(tail_round)
            break

        round_sites.sort(key=lambda s: s.get("score", 0), reverse=True)
        result.extend(round_sites)

    for site in result:
        site.pop("_region", None)

    return result
