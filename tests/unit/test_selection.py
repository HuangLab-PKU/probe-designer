"""TDD for probe_designer.scoring.selection.select_top_n_with_gap.

Phase 1 will add this function as a DB-agnostic extraction of the webapp's
_auto_select_top_n (task_runner.py:177). These tests lock the semantics
that both webapp-side and CLI-side callers rely on.

Expected signature:
    select_top_n_with_gap(
        sites: list[dict],
        top_n: int,
        min_gap: int = 40,
        position_fn: Callable[[dict], float] | None = None,
        score_key: str = "score",
    ) -> list[dict]

Semantics:
  - Sorts by score DESC (ties allowed; implementation-defined for ties)
  - Greedy selection: pick top-score, then next-top whose position is
    >= min_gap away from any already-picked position, skip otherwise
  - Returns up to top_n items, in selection order (highest-scoring first)
  - position_fn default: midpoint of site["st"] and site["en"] if both
    present, else site["st"], else no position constraint applied
  - Sites without derivable positions are eligible (no spacing check)
"""
from __future__ import annotations

import pytest

from probe_designer.scoring.selection import select_top_n_with_gap


def _site(score: float, st: int, en: int | None = None) -> dict:
    return {"score": score, "st": st, "en": en if en is not None else st + 40}


class TestSelectTopN:
    def test_empty_list_returns_empty(self):
        assert select_top_n_with_gap([], top_n=3) == []

    def test_top_n_limit_respected(self):
        sites = [_site(score=float(i), st=i * 100) for i in range(10)]
        result = select_top_n_with_gap(sites, top_n=3)
        assert len(result) == 3

    def test_returns_sorted_by_score_descending(self):
        sites = [
            _site(score=5.0, st=0),
            _site(score=9.0, st=200),
            _site(score=7.0, st=400),
        ]
        result = select_top_n_with_gap(sites, top_n=3)
        assert [s["score"] for s in result] == [9.0, 7.0, 5.0]


class TestMinGapConstraint:
    def test_min_gap_40_skips_nearby(self):
        sites = [
            _site(score=10.0, st=100),
            _site(score=9.0, st=120),  # 20nt away - too close
            _site(score=8.0, st=200),  # 100nt away - ok
        ]
        result = select_top_n_with_gap(sites, top_n=3, min_gap=40)
        # Should pick score 10 (st=100) then score 8 (st=200), skip score 9 (st=120)
        assert len(result) == 2
        picked_sts = [s["st"] for s in result]
        assert 100 in picked_sts
        assert 200 in picked_sts
        assert 120 not in picked_sts

    def test_min_gap_zero_allows_all(self):
        sites = [_site(score=float(i), st=100 + i) for i in range(5)]
        result = select_top_n_with_gap(sites, top_n=5, min_gap=0)
        assert len(result) == 5

    def test_min_gap_uses_midpoint(self):
        # Site with st=100, en=200 has midpoint 150
        # Site with st=120, en=180 has midpoint 150 (same midpoint!) - too close
        # Site with st=300, en=400 has midpoint 350 - far enough
        sites = [
            {"score": 10.0, "st": 100, "en": 200},
            {"score": 9.0, "st": 120, "en": 180},
            {"score": 8.0, "st": 300, "en": 400},
        ]
        result = select_top_n_with_gap(sites, top_n=3, min_gap=40)
        assert len(result) == 2


class TestWebappParity:
    """Matches behavior of webapp/task_runner.py:_auto_select_top_n lines 177-241."""

    def test_greedy_not_optimal(self):
        # Greedy picks 10 first, blocking 9; 8 at far position selected.
        sites = [
            _site(score=10.0, st=100),
            _site(score=9.0, st=130),   # too close to 100
            _site(score=8.0, st=200),
        ]
        result = select_top_n_with_gap(sites, top_n=2, min_gap=40)
        assert [s["score"] for s in result] == [10.0, 8.0]

    def test_sites_without_position_allowed_unconditionally(self):
        # Webapp code: except (json.JSONDecodeError, IndexError, TypeError): pass
        # i.e. position-less sites bypass the spacing check
        sites = [
            {"score": 10.0},  # no st/en
            {"score": 9.0},
            {"score": 8.0},
        ]
        result = select_top_n_with_gap(sites, top_n=3, min_gap=40)
        assert len(result) == 3


class TestReturnOrder:
    def test_result_in_selection_order(self):
        # Highest-scoring first; spacing-skipped items don't appear
        sites = [
            _site(score=5.0, st=0),
            _site(score=9.0, st=1000),
            _site(score=7.0, st=2000),
        ]
        result = select_top_n_with_gap(sites, top_n=5, min_gap=40)
        assert [s["score"] for s in result] == [9.0, 7.0, 5.0]


class TestTopNEdgeCases:
    def test_top_n_zero_returns_empty(self):
        sites = [_site(score=5.0, st=0)]
        assert select_top_n_with_gap(sites, top_n=0) == []

    def test_top_n_larger_than_available(self):
        sites = [_site(score=float(i), st=i * 200) for i in range(3)]
        result = select_top_n_with_gap(sites, top_n=100, min_gap=40)
        assert len(result) == 3
