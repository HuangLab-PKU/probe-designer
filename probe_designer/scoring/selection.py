"""Top-N selection with spatial spacing.

DB-agnostic extraction of webapp's `_auto_select_top_n` (task_runner.py:177-241).
Used by both the CLI `probe-design select` subcommand and (in Phase 2) by the
Pipeline class. The webapp's task_runner will call this function and project
the returned list onto `JobTargetLink.selected` rows.
"""
from __future__ import annotations

from typing import Callable, List, Dict, Any, Optional


def _default_position(site: Dict[str, Any]) -> Optional[float]:
    """Midpoint of (st, en) if both present; else st; else None."""
    st = site.get("st")
    en = site.get("en")
    if st is None:
        return None
    if en is None:
        return float(st)
    return (float(st) + float(en)) / 2.0


def select_top_n_with_gap(
    sites: List[Dict[str, Any]],
    top_n: int,
    min_gap: int = 40,
    position_fn: Optional[Callable[[Dict[str, Any]], Optional[float]]] = None,
    score_key: str = "score",
) -> List[Dict[str, Any]]:
    """Greedy top-N selection ordered by score with a minimum spatial gap.

    Matches webapp semantics:
      - Sort by score DESC (missing score treated as 0)
      - Walk the sorted list; include a site only if its position is at least
        `min_gap` from all already-selected positions
      - Stop at `top_n` selections
      - Sites with no derivable position bypass the spacing check (matches
        webapp's `except (json.JSONDecodeError, IndexError, TypeError): pass`)

    Args:
        sites: candidate sites; each dict is opaque except for optional keys
            `st`, `en`, and `score`.
        top_n: maximum number of sites to return.
        min_gap: minimum allowed distance between selected site positions.
        position_fn: extract a numeric position from a site. Defaults to
            midpoint of (st, en); falls back to `st`; returns None if neither.
        score_key: sort key for descending order.

    Returns:
        Selected sites in selection order (highest-scoring first).
    """
    if top_n <= 0 or not sites:
        return []

    pos_fn = position_fn or _default_position

    sorted_sites = sorted(
        sites, key=lambda s: s.get(score_key, 0), reverse=True
    )

    selected: List[Dict[str, Any]] = []
    selected_positions: List[float] = []

    for site in sorted_sites:
        if len(selected) >= top_n:
            break

        pos = pos_fn(site)
        if pos is not None:
            too_close = any(abs(pos - sp) < min_gap for sp in selected_positions)
            if too_close:
                continue
            selected_positions.append(pos)

        selected.append(site)

    return selected
