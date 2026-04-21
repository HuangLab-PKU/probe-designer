"""Tests for probe_designer.scoring.compute_off_target_score.

Phase 2 migration target. Locks the webapp formula ``10 - n * 0.5``.
"""
from __future__ import annotations

import pytest

from probe_designer.scoring import compute_off_target_score


class TestBasicFormula:
    def test_zero_alignments_gives_max_score(self):
        assert compute_off_target_score([]) == 10.0

    def test_one_alignment_subtracts_half(self):
        assert compute_off_target_score([{"hit_def": "x"}]) == 9.5

    def test_ten_alignments_gives_five(self):
        alns = [{"hit_def": "x"} for _ in range(10)]
        assert compute_off_target_score(alns) == 5.0

    def test_many_alignments_clamped_to_zero(self):
        alns = [{"hit_def": "x"} for _ in range(100)]
        assert compute_off_target_score(alns) == 0.0


class TestTargetOrganismFiltering:
    def test_filters_by_organism_substring(self):
        alns = [
            {"hit_def": "Homo sapiens ACTB mRNA"},
            {"hit_def": "Mus musculus Actb mRNA"},
            {"hit_def": "Danio rerio actb mRNA"},
        ]
        # Counting only Homo sapiens hits: 1 -> 10 - 0.5 = 9.5
        score = compute_off_target_score(alns, target_organism="Homo sapiens")
        assert score == 9.5

    def test_organism_match_case_insensitive(self):
        alns = [
            {"hit_def": "HOMO SAPIENS ACTB mRNA"},
            {"hit_def": "homo sapiens different gene"},
        ]
        score = compute_off_target_score(alns, target_organism="Homo sapiens")
        assert score == 9.0  # 2 matches

    def test_no_organism_counts_all(self):
        alns = [{"hit_def": "x"} for _ in range(4)]
        assert compute_off_target_score(alns) == 8.0

    def test_missing_hit_def_treated_as_empty(self):
        alns = [{}, {"hit_def": None}, {"hit_def": "Homo sapiens hit"}]
        score = compute_off_target_score(alns, target_organism="Homo sapiens")
        assert score == 9.5  # only the third matches


class TestParameterOverrides:
    def test_custom_slope(self):
        alns = [{"hit_def": "x"} for _ in range(4)]
        # slope=1.0 -> 10 - 4 = 6.0
        assert compute_off_target_score(alns, slope=1.0) == 6.0

    def test_custom_max_score(self):
        alns = [{"hit_def": "x"} for _ in range(2)]
        # max=5, slope=0.5 -> 5 - 1 = 4.0
        assert compute_off_target_score(alns, max_score=5.0) == 4.0


class TestRounding:
    def test_output_rounded_to_three_decimals(self):
        alns = [{"hit_def": "x"} for _ in range(3)]
        # 10 - 3 * 0.3 = 9.1
        assert compute_off_target_score(alns, slope=0.3) == 9.1
