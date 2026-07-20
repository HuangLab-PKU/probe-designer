"""Top-N selection with round-robin-across-clusters spacing.

Binding sites tend to cluster — there's usually a "peak" of nearby acceptable
sites with very similar scores, separated from the next peak by hundreds of
bases. A pure score-greedy pick would stack the first few selections in the
densest peak and leave other peaks uncovered. Users designing padlock panels
want the opposite: spread the selections across as many peaks as possible for
gene-wide coverage, then only double-dip into already-visited peaks if more
sites are needed.

Algorithm:
  1. Bucket sites into regions of ``region_size`` bp using ``position_fn``.
  2. Within each region, order by score DESC.
  3. Round-robin across regions: each pass picks one site from each region
     whose position is >= ``min_gap`` from every already-selected position.
     Round picks are committed in score-DESC order so the highest-scoring
     picks always come first in the output.
  4. Stop when ``top_n`` picks made or no region can contribute without
     violating ``min_gap``.
  5. Sites without a derivable position bypass the spacing check and fill
     any remaining slots at the end (matches the webapp's historical
     "allow selection when position info is missing" branch).
"""
from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional


def _default_position(site: Dict[str, Any]) -> Optional[float]:
    """Midpoint of (st, en) if both present; else st; else None."""
    st = site.get("st")
    en = site.get("en")
    if st is None:
        return None
    if en is None:
        return float(st)
    return (float(st) + float(en)) / 2.0


def select_score_peaks(
    sites: List[Dict[str, Any]],
    min_distance: int,
    *,
    max_n: Optional[int] = None,
    score_key: str = "score",
    position_fn: Optional[Callable[[Dict[str, Any]], Optional[float]]] = None,
) -> List[Dict[str, Any]]:
    """Greedy local-maximum (non-maximum-suppression) peak selection.

    Repeatedly take the highest-scoring remaining site and suppress every other
    site within ``min_distance`` bp of it, so the selections are the score
    "peaks" of the landscape spaced far enough apart that their binding sites
    do not overlap. Equivalent to ``scipy.signal.find_peaks(distance=...)`` on a
    score-vs-position signal, but operates directly on sparse candidate sites.

    Args:
        sites: candidate sites; each opaque except ``st``/``en`` and ``score_key``.
        min_distance: minimum bp between selected positions — set to the binding-
            site length so selected probes cannot overlap.
        max_n: optional cap on the number of peaks returned (None = all peaks).
        score_key: dict key used for the DESC ordering.
        position_fn: extract a numeric position; defaults to the ``(st, en)``
            midpoint. Unpositioned sites are dropped (they have no overlap to
            suppress and no defined peak location).

    Returns:
        Selected peaks in score-DESC order (highest first). Length <= ``max_n``.
    """
    if not sites or min_distance < 0:
        return []
    pos_fn = position_fn or _default_position

    positioned = [(s, pos_fn(s)) for s in sites]
    positioned = [(s, p) for s, p in positioned if p is not None]
    # Highest score first; the greedy pick then defines each peak and its
    # suppression zone. Ties keep input order (Python sort is stable).
    positioned.sort(key=lambda sp: sp[0].get(score_key, 0), reverse=True)

    selected: List[Dict[str, Any]] = []
    selected_pos: List[float] = []
    for site, pos in positioned:
        if any(abs(pos - q) < min_distance for q in selected_pos):
            continue  # inside an already-claimed peak's exclusion zone
        selected.append(site)
        selected_pos.append(pos)
        if max_n is not None and len(selected) >= max_n:
            break
    return selected


def select_top_n_with_gap(
    sites: List[Dict[str, Any]],
    top_n: int,
    min_gap: int = 40,
    *,
    region_size: int = 80,
    position_fn: Optional[Callable[[Dict[str, Any]], Optional[float]]] = None,
    score_key: str = "score",
) -> List[Dict[str, Any]]:
    """Round-robin-across-clusters top-N selection (see module docstring).

    Args:
        sites: candidate sites; each dict is opaque except for the optional
            ``st`` / ``en`` position keys and ``score_key``.
        top_n: maximum number of sites to return.
        min_gap: minimum allowed bp distance between selected site positions.
        region_size: bucket width for clustering (kept as a keyword so callers
            wanting the historical "no clustering" semantic can pass a large
            value).
        position_fn: extract a numeric position from a site. Defaults to
            midpoint of ``(st, en)``, falling back to ``st``; returns None
            if neither is present.
        score_key: dict key used for the DESC ordering.

    Returns:
        Selected sites in selection order (first picks = highest-scoring of
        the first round). Length <= ``top_n``.
    """
    if top_n <= 0 or not sites:
        return []

    pos_fn = position_fn or _default_position

    # --- Step 1: bucket by region ---
    buckets: Dict[int, List[Dict[str, Any]]] = {}
    unpositioned: List[Dict[str, Any]] = []
    for site in sites:
        pos = pos_fn(site)
        if pos is None:
            unpositioned.append(site)
        else:
            region_id = int(pos) // region_size
            buckets.setdefault(region_id, []).append(site)

    # --- Step 2: sort each bucket by score DESC ---
    for bucket in buckets.values():
        bucket.sort(key=lambda s: s.get(score_key, 0), reverse=True)

    # --- Step 3: round-robin across buckets ---
    selected: List[Dict[str, Any]] = []
    selected_positions: List[float] = []

    while len(selected) < top_n and buckets:
        round_picks: List[tuple[Dict[str, Any], float]] = []
        empty_regions: List[int] = []

        for region_id in sorted(buckets.keys()):
            bucket = buckets[region_id]
            picked_index = None
            picked_pos: Optional[float] = None

            for i, candidate in enumerate(bucket):
                pos = pos_fn(candidate)
                # position_fn can't return None here (we filtered those out
                # into `unpositioned`) but defend anyway.
                if pos is None:
                    continue
                if any(abs(pos - sp) < min_gap for sp in selected_positions):
                    continue
                picked_index = i
                picked_pos = pos
                break

            if picked_index is not None:
                round_picks.append((bucket.pop(picked_index), picked_pos))  # type: ignore[arg-type]
            if not bucket:
                empty_regions.append(region_id)

        if not round_picks:
            break  # every remaining bucket head violates min_gap — give up

        # Commit this round's picks ordered by score DESC so the output
        # starts with the best-scoring site of the run.
        round_picks.sort(key=lambda sp: sp[0].get(score_key, 0), reverse=True)
        for site, pos in round_picks:
            if len(selected) >= top_n:
                break
            selected.append(site)
            selected_positions.append(pos)

        for rid in empty_regions:
            buckets.pop(rid, None)

    # --- Step 4: fill any remaining slots with unpositioned sites ---
    if len(selected) < top_n and unpositioned:
        unpositioned.sort(key=lambda s: s.get(score_key, 0), reverse=True)
        for site in unpositioned:
            if len(selected) >= top_n:
                break
            selected.append(site)

    return selected
