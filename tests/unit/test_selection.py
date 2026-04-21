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


class TestRoundRobinClustering:
    """Round-robin-across-clusters is the core selection semantic.

    Binding sites come in clusters of nearby high-scoring positions; the
    selector covers one-per-cluster first, then double-dips only after
    every cluster has contributed. Greedy-by-score alone would pile picks
    into the densest cluster.
    """

    def test_round_robin_prefers_cluster_spread_over_raw_score(self):
        # Cluster A (mids 120 & 130, both in region 1 at region_size=80)
        # packs scores 10 and 9.9; cluster B is a far lone site at mid=1020
        # with score 5. Pure greedy would pick [10, 9.9] and miss the far
        # cluster entirely (because 9.9's midpoint is 10nt from 10's midpoint
        # — still >= min_gap=10 in the greedy check). Round-robin picks one
        # from each region first.
        sites = [
            _site(score=10.0, st=100),    # cluster A, mid=120
            _site(score=9.9, st=110),     # cluster A, mid=130 (same region as above)
            _site(score=5.0, st=1000),    # cluster B, mid=1020 (region 12)
        ]
        # Use min_gap=5 so a pure greedy would happily take 120 then 130
        # (|130-120|=10 >= 5); the round-robin version must still prefer spread.
        result = select_top_n_with_gap(
            sites, top_n=2, min_gap=5, region_size=80,
        )
        scores = sorted(s["score"] for s in result)
        assert 5.0 in scores, f"round-robin should reach cluster B; got {scores}"
        assert 10.0 in scores, f"cluster A's best must be picked; got {scores}"
        assert 9.9 not in scores, f"greedy-not-round-robin would wrongly pick 9.9; got {scores}"

    def test_second_round_dips_back_into_cluster_when_more_needed(self):
        # top_n=3 with 2 clusters (3 sites vs 1): round 1 picks one each,
        # round 2 goes back into the bigger cluster for the next-best.
        sites = [
            _site(score=10.0, st=100),    # cluster A
            _site(score=9.0, st=200),     # cluster A (different region: 220//80=2 vs 120//80=1)
            _site(score=8.5, st=300),     # cluster A (mid=320, region 4)
            _site(score=6.0, st=2000),    # cluster B, far away
        ]
        # region_size=500 coalesces 100-300 into cluster "A" for this test
        result = select_top_n_with_gap(
            sites, top_n=3, min_gap=40, region_size=500,
        )
        scores = sorted(s["score"] for s in result)
        # Round 1: best of A (10) + best of B (6).
        # Round 2: 2nd-best of A whose midpoint is >=40 away — st=200 mid=220, |220-120|=100 OK -> 9.
        assert scores == [6.0, 9.0, 10.0]

    def test_min_gap_blocks_dip_back_into_tight_cluster(self):
        # Tight cluster: all 3 sites within min_gap of each other.
        # After round 1 picks one, the rest can't be added without violating.
        sites = [
            _site(score=10.0, st=100),   # mid=120
            _site(score=9.0, st=110),    # mid=130; 10 from 120 < 40
            _site(score=8.0, st=120),    # mid=140; 20 from 120 < 40
            _site(score=5.0, st=2000),   # separate cluster
        ]
        result = select_top_n_with_gap(
            sites, top_n=4, min_gap=40, region_size=500,
        )
        scores = sorted(s["score"] for s in result)
        # Expect exactly 2 picks: best of each cluster.
        assert scores == [5.0, 10.0]

    def test_cluster_best_always_wins_over_distant_lower(self):
        # Cluster A has 10.0 AND 9.9 (adjacent); one far lone site at 9.0.
        # Round-robin round 1: A best (10) + B best (9). Committed in score
        # order = [10, 9]. 9.9 stays unpicked.
        sites = [
            _site(score=10.0, st=100),    # cluster A
            _site(score=9.9, st=150),     # cluster A (same region with st=100 when region_size=80: 120//80=1, 170//80=2)
            _site(score=9.0, st=5000),    # cluster B far away
        ]
        # Use a larger region_size to be sure the 2 cluster-A sites share a region
        result = select_top_n_with_gap(
            sites, top_n=2, min_gap=40, region_size=200,
        )
        scores = [s["score"] for s in result]
        assert scores == [10.0, 9.0]


class TestPositionlessSites:
    def test_sites_without_position_fill_after_clustered(self):
        # Positioned sites participate in round-robin; positionless fill at end.
        sites = [
            _site(score=10.0, st=100),    # positioned
            _site(score=9.0, st=1000),    # positioned, different region
            {"score": 7.0},               # no st/en
            {"score": 6.0},
        ]
        result = select_top_n_with_gap(sites, top_n=4, min_gap=40)
        assert len(result) == 4
        # First two must be the positioned pair (both regions covered)
        first_two_scores = {result[0]["score"], result[1]["score"]}
        assert first_two_scores == {9.0, 10.0}

    def test_only_positionless_sites_returned_by_score_desc(self):
        sites = [
            {"score": 10.0},
            {"score": 9.0},
            {"score": 8.0},
        ]
        result = select_top_n_with_gap(sites, top_n=3, min_gap=40)
        assert [s["score"] for s in result] == [10.0, 9.0, 8.0]


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
