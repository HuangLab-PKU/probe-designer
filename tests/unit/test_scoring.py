"""Characterization tests for probe_designer.scoring.

Locks current behavior of compute_target_score() and peak_rank() before they
move into the new scoring/ subpackage in Phase 1. Any numerical drift
(> 1e-6) here must be explicitly approved.

Known discrepancy between code and docstring (NOT fixed here — characterization
only): scoring.py:71-75 says "dG=0 -> 1.0 (best)" but the guard is
`if mfe < 0`, so mfe=0 actually contributes 0. Flagged for Phase 1 review.
"""
from __future__ import annotations

import pytest

from probe_designer.scoring import compute_target_score, peak_rank


TOL = 1e-6


def _neutral_site(**overrides):
    """Site with all components at neutral values so overrides are additive.

    Base components (with total_isoforms=1 default):
      - isoform_coverage: 3.0 (1/1 * 3 capped at 3)
      - blast_support:    0   (no alignments)
      - tm_proximity:     2.0 (avg=50 == min_arm_tm, no excess)
      - tm_balance:       1.0 (diff=0)
      - terminal_gc:      0   (A...T arms)
      - delta_g:          0   (mfe=0 fails `if mfe<0`)
      Total: 6.0
    """
    base = {
        "arm_5prime": "ATATATATAT",  # neither end is G/C
        "arm_3prime": "ATATATATAT",
        "tm_5prime": 50.0,
        "tm_3prime": 50.0,
        "isoform_overlap_num": 1,
        "blast_alignments": [],
        "free_energy": 0.0,
    }
    base.update(overrides)
    return base


BASE_SCORE = 6.0


class TestComputeTargetScoreBaseline:
    def test_neutral_site_base_score(self):
        assert compute_target_score(_neutral_site()) == pytest.approx(BASE_SCORE, abs=TOL)


class TestComputeTargetScoreIsoform:
    def test_full_coverage_gives_three(self):
        s = compute_target_score(_neutral_site(isoform_overlap_num=3), total_isoforms=3)
        assert s == pytest.approx(BASE_SCORE, abs=TOL)  # 3.0 iso same as base 1/1

    def test_partial_coverage_half(self):
        # 2/4 * 3 = 1.5 (vs base 3.0) => -1.5
        s = compute_target_score(_neutral_site(isoform_overlap_num=2), total_isoforms=4)
        assert s == pytest.approx(BASE_SCORE - 1.5, abs=TOL)

    def test_no_isoform_data_gives_neutral_one_point_five(self):
        # total_isoforms=0 => isoform component=1.5 (vs default 3.0) => -1.5
        s = compute_target_score(_neutral_site(), total_isoforms=0)
        assert s == pytest.approx(BASE_SCORE - 1.5, abs=TOL)

    def test_coverage_capped_at_three(self):
        # overlap > total still caps at 3.0
        s = compute_target_score(_neutral_site(isoform_overlap_num=10), total_isoforms=1)
        assert s == pytest.approx(BASE_SCORE, abs=TOL)


class TestComputeTargetScoreBlast:
    def test_five_mrna_hits_saturate_at_two(self):
        site = _neutral_site(
            blast_alignments=[
                {"hit_id": "NM_001", "hit_def": "Homo sapiens ACTB, mRNA"}
                for _ in range(5)
            ]
        )
        # +2.0 blast
        assert compute_target_score(site) == pytest.approx(BASE_SCORE + 2.0, abs=TOL)

    def test_over_five_still_caps_at_two(self):
        site = _neutral_site(
            blast_alignments=[{"hit_id": f"NM_{i}", "hit_def": "mRNA"} for i in range(20)]
        )
        assert compute_target_score(site) == pytest.approx(BASE_SCORE + 2.0, abs=TOL)

    def test_non_mrna_hit_ids_give_zero(self):
        site = _neutral_site(
            blast_alignments=[{"hit_id": "XR_1", "hit_def": "ncRNA or lnc"}]
        )
        assert compute_target_score(site) == pytest.approx(BASE_SCORE, abs=TOL)

    def test_xm_hit_id_counts_as_mrna(self):
        site = _neutral_site(
            blast_alignments=[{"hit_id": f"XM_{i}", "hit_def": "predicted"} for i in range(5)]
        )
        assert compute_target_score(site) == pytest.approx(BASE_SCORE + 2.0, abs=TOL)

    def test_mrna_in_hit_def_counts(self):
        site = _neutral_site(
            blast_alignments=[
                {"hit_id": "gi|123|", "hit_def": "Homo sapiens actin mRNA"}
                for _ in range(5)
            ]
        )
        assert compute_target_score(site) == pytest.approx(BASE_SCORE + 2.0, abs=TOL)


class TestComputeTargetScoreTm:
    def test_tm_proximity_penalty_at_60(self):
        # avg=60, excess=10, 2*(1-10/20)=1.0 (vs base 2.0) => -1.0
        s = compute_target_score(_neutral_site(tm_5prime=60.0, tm_3prime=60.0))
        assert s == pytest.approx(BASE_SCORE - 1.0, abs=TOL)

    def test_tm_proximity_zero_at_far_excess(self):
        # avg=80, dev=30, clamped to 0 (vs base 2.0) => -2.0
        s = compute_target_score(_neutral_site(tm_5prime=80.0, tm_3prime=80.0))
        assert s == pytest.approx(BASE_SCORE - 2.0, abs=TOL)

    def test_tm_proximity_penalized_below_target(self):
        # 2026-07-19: two-sided. avg=40 is 10 BELOW target(50) => dev=10 =>
        # 2*(1-10/20)=1.0 (vs base 2.0) => -1.0. (Old one-sided gave full 2.0.)
        s = compute_target_score(_neutral_site(tm_5prime=40.0, tm_3prime=40.0))
        assert s == pytest.approx(BASE_SCORE - 1.0, abs=TOL)

    def test_tm_balance_linear(self):
        # diff=5 over max_tm_diff=10 => 0.5 (vs base 1.0) => -0.5
        # avg=52.5, excess=2.5 => tm_prox = 2*(1-0.125) = 1.75 (vs 2.0) => -0.25
        # Net = -0.75
        s = compute_target_score(_neutral_site(tm_5prime=50.0, tm_3prime=55.0))
        assert s == pytest.approx(BASE_SCORE - 0.75, abs=TOL)


class TestComputeTargetScoreTerminalGc:
    def test_arm5_start_g_adds_half(self):
        s = compute_target_score(_neutral_site(arm_5prime="GTATATAT"))
        assert s == pytest.approx(BASE_SCORE + 0.5, abs=TOL)

    def test_arm3_end_c_adds_half(self):
        s = compute_target_score(_neutral_site(arm_3prime="TATATATATC"))
        assert s == pytest.approx(BASE_SCORE + 0.5, abs=TOL)

    def test_both_gc_adds_one(self):
        s = compute_target_score(
            _neutral_site(arm_5prime="CATATA", arm_3prime="TATATG")
        )
        assert s == pytest.approx(BASE_SCORE + 1.0, abs=TOL)

    def test_gc_case_insensitive(self):
        s_upper = compute_target_score(_neutral_site(arm_5prime="GATATA"))
        s_lower = compute_target_score(_neutral_site(arm_5prime="gatata"))
        assert s_upper == pytest.approx(s_lower, abs=TOL)


class TestComputeTargetScoreDeltaG:
    """Locks current (buggy-per-docstring) delta_g behavior: only mfe<0 gets a bonus.

    The docstring at scoring.py:71-75 claims 'dG=0 → 1.0 (best)' but the
    implementation at line 74 uses `if mfe < 0:` so mfe=0 contributes 0.
    Phase 1 TODO: decide whether to fix docstring or implementation.
    """

    def test_mfe_zero_gives_no_bonus(self):
        # mfe=0 -> not <0 -> 0 bonus (current behavior)
        assert compute_target_score(_neutral_site(free_energy=0.0)) == pytest.approx(BASE_SCORE, abs=TOL)

    def test_mfe_minus_five_gives_half_bonus(self):
        # 1 - 5/10 = 0.5
        s = compute_target_score(_neutral_site(free_energy=-5.0))
        assert s == pytest.approx(BASE_SCORE + 0.5, abs=TOL)

    def test_mfe_minus_ten_gives_zero_bonus(self):
        # 1 - 10/10 = 0 (boundary)
        s = compute_target_score(_neutral_site(free_energy=-10.0))
        assert s == pytest.approx(BASE_SCORE, abs=TOL)

    def test_mfe_minus_one_near_max_bonus(self):
        # 1 - 1/10 = 0.9
        s = compute_target_score(_neutral_site(free_energy=-1.0))
        assert s == pytest.approx(BASE_SCORE + 0.9, abs=TOL)

    def test_mfe_very_negative_clamped_to_zero_bonus(self):
        # 1 - 20/10 = -1, clamped to 0
        s = compute_target_score(_neutral_site(free_energy=-20.0))
        assert s == pytest.approx(BASE_SCORE, abs=TOL)

    def test_mfe_positive_treated_as_no_data(self):
        # mfe>0 -> not <0 -> 0 bonus (same as mfe=0)
        s = compute_target_score(_neutral_site(free_energy=1.5))
        assert s == pytest.approx(BASE_SCORE, abs=TOL)


class TestComputeTargetScoreMax:
    def test_max_possible_score_is_ten(self):
        # All components maxed: iso 3 + blast 2 + tm_prox 2 + tm_bal 1 + gc 1 + dG 0.9 = 9.9
        # Note: dG max is 0.9 at mfe=-1 (mfe=0 gets 0 per current bug)
        # So the 'documented' 10.0 ceiling is actually 9.9
        s = compute_target_score(
            _neutral_site(
                isoform_overlap_num=5,
                blast_alignments=[{"hit_id": f"NM_{i}", "hit_def": "mRNA"} for i in range(5)],
                arm_5prime="GATAT", arm_3prime="TATAC",
                free_energy=-1.0,
            ),
            total_isoforms=5,
        )
        assert s == pytest.approx(9.9, abs=TOL)


class TestComputeTargetScoreRounding:
    def test_output_rounded_to_three_decimals(self):
        s = compute_target_score(
            _neutral_site(tm_5prime=50.123456, tm_3prime=50.987654)
        )
        assert s == round(s, 3)


class TestPeakRank:
    def test_empty_returns_empty(self):
        assert peak_rank([]) == []

    def test_single_site_returned(self):
        sites = [{"st": 100, "score": 5.0}]
        ranked = peak_rank(sites)
        assert len(ranked) == 1
        assert ranked[0]["st"] == 100
        assert "_region" not in ranked[0]  # temp field cleaned

    def test_region_grouping_and_round_one_sorts_by_score(self):
        # Region 0 (0-79) has best score 8 at st=50; region 1 (80-159) has 5 at st=100
        sites = [
            {"st": 10, "score": 3.0},
            {"st": 50, "score": 8.0},
            {"st": 100, "score": 5.0},
        ]
        ranked = peak_rank(sites, region_size=80, min_gap=40)
        # First two picks = one from each region, sorted by score desc
        assert [s["score"] for s in ranked[:2]] == [8.0, 5.0]
        # All three returned eventually
        assert len(ranked) == 3

    def test_all_sites_returned_exactly_once(self):
        sites = [{"st": i * 10, "score": float(i)} for i in range(5)]
        ranked = peak_rank(sites)
        assert len(ranked) == 5
        assert sorted(s["st"] for s in ranked) == [0, 10, 20, 30, 40]

    def test_temp_region_field_cleaned_up(self):
        sites = [{"st": i * 30, "score": float(i)} for i in range(6)]
        ranked = peak_rank(sites)
        for s in ranked:
            assert "_region" not in s

    def test_score_missing_treated_as_zero(self):
        sites = [{"st": 0}, {"st": 100, "score": 5.0}]
        ranked = peak_rank(sites, region_size=50, min_gap=10)
        # Score 5.0 should come before score-missing in round 1
        assert ranked[0]["score"] == 5.0

    def test_four_regions_each_represented_first_round(self):
        sites = []
        for region_idx in range(4):
            for i in range(3):
                sites.append({
                    "st": region_idx * 100 + i * 10,
                    "score": float(region_idx * 3 + i),
                })
        ranked = peak_rank(sites, region_size=100, min_gap=40)
        first_four_regions = {s["st"] // 100 for s in ranked[:4]}
        assert first_four_regions == {0, 1, 2, 3}

    def test_scores_all_equal(self):
        sites = [{"st": i * 100, "score": 5.0} for i in range(3)]
        ranked = peak_rank(sites)
        assert len(ranked) == 3
        assert all(s["score"] == 5.0 for s in ranked)
