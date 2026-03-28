"""Target scoring and peak-finding ranking for padlock probe design.

Provides:
- compute_target_score(): composite quality score for a binding site
- peak_rank(): position-aware ranking that spreads top targets across the gene
"""

from typing import Dict, List, Any, Optional


def compute_target_score(
    site: Dict[str, Any],
    min_arm_tm: float = 50.0,
    max_tm_diff: float = 10.0,
    total_isoforms: int = 1,
) -> float:
    """Compute composite quality score for a target region (0-10 scale).

    Components (higher = better):
    1. isoform_coverage (0-3): fraction of isoforms covered
    2. blast_support (0-2): same-gene mRNA hit count (after specificity filter)
    3. tm_proximity (0-2): arms closer to min_arm_tm preferred
    4. tm_balance (0-1): smaller |tm_5prime - tm_3prime|
    5. terminal_gc (0-1): G/C at padlock ligation point (arm5 5'-end, arm3 3'-end)
    6. delta_g (0-1): lower MFE = more stable = better
    """
    score = 0.0

    # 1. Isoform coverage (0-3)
    overlap = site.get("isoform_overlap_num", 1)
    if total_isoforms > 0:
        score += min(3.0, (overlap / total_isoforms) * 3.0)
    else:
        score += 1.5  # neutral if no isoform data

    # 2. BLAST support (0-2)
    alignments = site.get("blast_alignments", [])
    mrna_hits = sum(
        1 for a in alignments
        if "mrna" in a.get("hit_def", "").lower()
        or "nm_" in a.get("hit_id", "").lower()
        or "xm_" in a.get("hit_id", "").lower()
    )
    score += min(2.0, mrna_hits / 5.0 * 2.0)

    # 3. Tm proximity to lower bound (0-2)
    tm5 = site.get("tm_5prime", 0)
    tm3 = site.get("tm_3prime", 0)
    if tm5 > 0 and tm3 > 0:
        avg_tm = (tm5 + tm3) / 2
        # Prefer Tm close to min_arm_tm (lower = better for specificity)
        # Score decreases as Tm moves away from min_arm_tm
        tm_excess = max(0, avg_tm - min_arm_tm)
        tm_range = 20.0  # normalize over 20°C range
        score += max(0, 2.0 * (1.0 - tm_excess / tm_range))

    # 4. Tm balance (0-1)
    tm_diff = abs(tm5 - tm3)
    if max_tm_diff > 0:
        score += max(0, 1.0 * (1.0 - tm_diff / max_tm_diff))

    # 5. Terminal GC clamp (0-1)
    arm5 = site.get("arm_5prime", "")
    arm3 = site.get("arm_3prime", "")
    if arm5 and arm5[0] in "GCgc":
        score += 0.5
    if arm3 and arm3[-1] in "GCgc":
        score += 0.5

    # 6. Delta G / free energy (0-1)
    mfe = site.get("free_energy", 0)
    if mfe < 0:
        # More negative = better (more stable secondary structure avoidance)
        # Normalize: -10 kcal/mol = 1.0, 0 = 0.0
        score += min(1.0, abs(mfe) / 10.0)

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

    # Assign region IDs
    for site in sites:
        site["_region"] = site.get("st", 0) // region_size

    # Group by region, sort each group by score descending
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

            # Find best candidate not overlapping with already-selected in this region
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

            # If all remaining overlap, skip this region for now
            if not picked:
                continue

        if not round_sites:
            # All remaining sites overlap with selections — append them at the end
            for v in regions.values():
                result.extend(v)
                v.clear()
            break

        # Within this round, sort by score descending
        round_sites.sort(key=lambda s: s.get("score", 0), reverse=True)
        result.extend(round_sites)

    # Clean up temp field
    for site in result:
        site.pop("_region", None)

    return result
