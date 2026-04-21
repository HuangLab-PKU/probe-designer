"""Off-target scoring formulas.

Phase 2 migration of webapp/services/ingest.py:252-256 — the
``10.0 - n_alns * 0.5`` formula now lives in the designer library so
both CLI and webapp share the same off-target score definition.

compute_target_score + peak_rank remain in scoring/__init__.py to preserve
the existing ``from probe_designer.scoring import ...`` surface.
"""
from __future__ import annotations

from typing import Any, Dict, Iterable, Optional


def compute_off_target_score(
    alignments: Iterable[Dict[str, Any]],
    *,
    target_organism: Optional[str] = None,
    slope: float = 0.5,
    max_score: float = 10.0,
) -> float:
    """Off-target penalty based on number of BLAST hits to the target organism.

    Formula (matches webapp/ingest.py:252-256):
        max(0, max_score - n_hits * slope)

    Args:
        alignments: iterable of BLAST alignment dicts with ``hit_def`` keys.
        target_organism: if given, only alignments whose ``hit_def`` contains
            the organism string (case-insensitive) are counted. If None,
            all alignments are counted (caller has pre-filtered).
        slope: points deducted per alignment (default 0.5).
        max_score: starting score with zero alignments (default 10.0).

    Returns:
        Score in [0, max_score], rounded to 3 decimals.
    """
    if target_organism:
        needle = target_organism.lower()
        filtered = [
            a for a in alignments
            if needle in (a.get("hit_def", "") or "").lower()
        ]
    else:
        filtered = list(alignments)

    n = len(filtered)
    return round(max(0.0, max_score - n * slope), 3)
