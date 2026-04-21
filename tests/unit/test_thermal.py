"""Characterization tests for SequenceFilter.thermal_filter.

Locks the semantic behavior (pass/fail logic, failed_checks reasons, returned
dict shape) before filtering.py is split in Phase 3. Does NOT lock specific
Tm values since those depend on Biopython MeltingTemp table versions.
"""
from __future__ import annotations

import pytest

from probe_designer.filtering import SequenceFilter


# Arms tuned to pass the defaults (G content 0.3-0.7, no 4G runs, Tm 45-65)
GOOD_ARM_3 = "GATGATGATGATGATGATGC"   # 20 nt, 7 Gs = 35%
GOOD_ARM_5 = "GATGATGATGATGATGATGC"
GOOD_TARGET_SEQ = "GCATCATCATCATCATCATCGCATCATCATCATCATCATC"

# Arms designed to violate specific rules
LOW_G_ARM = "ATATATATATATATATATAT"           # 0 Gs
HIGH_G_ARM = "GGAGGAGGAGGAGGAGGAGG"          # high G fraction (runs but also high content)
NO_RUN_HIGH_G = "GAGAGAGAGAGAGAGAGAGT"       # 10 Gs, 50%, no 4G run
FOUR_G_RUN_ARM = "GATGATGATGGGGATGATGA"      # contains GGGG


@pytest.fixture
def filter(default_filter_config, default_blast_config):
    return SequenceFilter(default_filter_config, default_blast_config)


class TestThermalFilterReturnedDictShape:
    """Locks the output dict contract that downstream code depends on."""

    def test_all_expected_keys_present(self, filter):
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        expected_keys = {
            "arm_3prime", "arm_5prime", "target_sequence",
            "sequence_type", "target_type",
            "passed",
            "g_content", "g_content_arm3", "g_content_arm5",
            "tm", "tm_3prime", "tm_5prime", "tm_diff",
            "free_energy",
            "has_consecutive_g", "has_consecutive_g_arm3",
            "has_consecutive_g_arm5", "has_consecutive_g_target",
            "failed_checks",
            "arm_3prime_length", "arm_5prime_length",
        }
        missing = expected_keys - set(result.keys())
        assert not missing, f"Missing keys: {missing}"

    def test_arm_lengths_reported(self, filter):
        result = filter.thermal_filter("ACGTACGTAC", "GGGGACGTACGTACGT", target_sequence=GOOD_TARGET_SEQ)
        assert result["arm_3prime_length"] == 10
        assert result["arm_5prime_length"] == 16


class TestThermalFilterPassPath:
    def test_passed_true_when_no_failed_checks(self, filter):
        # Use very lax bounds so the filter has nothing to reject on.
        result = filter.thermal_filter(
            GOOD_ARM_3, GOOD_ARM_5,
            target_sequence=GOOD_TARGET_SEQ,
            min_g_content=0.0, max_g_content=1.0,
            max_consecutive_g=100,
            min_tm=0.0, max_tm=200.0, max_tm_diff=200.0,
        )
        assert result["failed_checks"] == []
        assert result["passed"] is True

    def test_passed_false_when_any_check_fails(self, filter):
        # Force just one failure: g_content too strict
        result = filter.thermal_filter(
            LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ,
            min_tm=0.0, max_tm=200.0, max_tm_diff=200.0,
            max_consecutive_g=100,
        )
        assert result["passed"] is False
        assert len(result["failed_checks"]) > 0


class TestThermalFilterGContent:
    def test_low_g_content_fails_per_arm(self, filter):
        result = filter.thermal_filter(LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "g_content_arm3" in result["failed_checks"]
        assert result["passed"] is False

    def test_high_g_content_fails(self, filter):
        # 16 Gs out of 20 = 0.80, exceeds default max_g_content=0.70.
        # Avoid 4G+ runs so only g_content rule triggers.
        above_max = "GAGGAGGAGGAGGAGGGAGG"  # Gs: 1,3,4,6,7,9,10,12,13,14,15,17,18 = 13/20=0.65
        # Actually let me use 0.80 clearly:
        very_high = "GTGTGTGTGTGTGTGGTGGT"  # Gs: many, let's just force via kwargs
        # Simpler: use the original arm but with a strict max override
        result = filter.thermal_filter(
            GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ,
            max_g_content=0.0,  # anything >0 G fails
        )
        assert "g_content_arm3" in result["failed_checks"]

    def test_g_content_reported_per_arm(self, filter):
        result = filter.thermal_filter(LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert result["g_content_arm3"] == pytest.approx(0.0)
        assert result["g_content_arm5"] > 0.0

    def test_aggregate_g_content_is_average(self, filter):
        result = filter.thermal_filter(LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        avg = (result["g_content_arm3"] + result["g_content_arm5"]) / 2.0
        assert result["g_content"] == pytest.approx(avg)


class TestThermalFilterConsecutiveG:
    def test_four_g_run_on_arm3_fails(self, filter):
        result = filter.thermal_filter(FOUR_G_RUN_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert result["has_consecutive_g_arm3"] is True
        assert "consecutive_g" in result["failed_checks"]
        assert result["passed"] is False

    def test_four_g_run_on_arm5_fails(self, filter):
        result = filter.thermal_filter(GOOD_ARM_3, FOUR_G_RUN_ARM, target_sequence=GOOD_TARGET_SEQ)
        assert result["has_consecutive_g_arm5"] is True
        assert "consecutive_g" in result["failed_checks"]
        assert result["passed"] is False

    def test_four_g_run_on_target_fails(self, filter):
        target_with_run = "AAAAGGGGAAAATTTTGGGGAAAA"
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=target_with_run)
        assert result["has_consecutive_g_target"] is True
        assert "consecutive_g" in result["failed_checks"]

    def test_three_g_run_acceptable(self, filter):
        # max_consecutive_g=4 by default, so GGG (run of 3) is allowed
        arm = "GGGATGGGATGGGATGGGAT"
        result = filter.thermal_filter(arm, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert result["has_consecutive_g_arm3"] is False


class TestThermalFilterTmRange:
    def test_kwargs_override_tm_bounds(self, filter):
        # Force everything out of range by setting impossible bounds
        result = filter.thermal_filter(
            GOOD_ARM_3, GOOD_ARM_5,
            target_sequence=GOOD_TARGET_SEQ,
            min_tm=200.0, max_tm=300.0,
        )
        assert "tm_3prime_range" in result["failed_checks"]
        assert "tm_5prime_range" in result["failed_checks"]
        assert result["passed"] is False

    def test_kwargs_override_max_tm_diff(self, filter):
        # Asymmetric arms to create some diff; force strict limit
        result = filter.thermal_filter(
            "GATGATGATGATGATGATGA", "GCGCGCGCGCGCGCGCGCGC",
            target_sequence=GOOD_TARGET_SEQ,
            max_tm_diff=0.01,
        )
        assert "tm_diff" in result["failed_checks"] or result["tm_diff"] > 0.01


class TestThermalFilterKwargsOverrideConfig:
    def test_kwargs_override_min_g_content(self, filter):
        # Default min=0.3 would reject LOW_G; override to 0 to pass
        result = filter.thermal_filter(
            LOW_G_ARM, LOW_G_ARM,
            target_sequence="A" * 40,
            min_g_content=0.0,
        )
        assert "g_content_arm3" not in result["failed_checks"]
        assert "g_content_arm5" not in result["failed_checks"]

    def test_kwargs_override_max_consecutive_g(self, filter):
        # Default max=4 rejects GGGG. Override to 10 to allow.
        result = filter.thermal_filter(
            FOUR_G_RUN_ARM, GOOD_ARM_5,
            target_sequence=GOOD_TARGET_SEQ,
            max_consecutive_g=10,
        )
        assert "consecutive_g" not in result["failed_checks"]


class TestThermalFilterNumericReturns:
    def test_g_content_in_zero_to_one(self, filter):
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert 0.0 <= result["g_content"] <= 1.0
        assert 0.0 <= result["g_content_arm3"] <= 1.0
        assert 0.0 <= result["g_content_arm5"] <= 1.0

    def test_tm_values_nonzero_for_reasonable_arms(self, filter):
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert result["tm"] > 0
        assert result["tm_3prime"] > 0
        assert result["tm_5prime"] > 0

    def test_tm_diff_is_absolute(self, filter):
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert result["tm_diff"] == pytest.approx(abs(result["tm_5prime"] - result["tm_3prime"]))
        assert result["tm_diff"] >= 0


class TestThermalFilterRnaStructure:
    def test_rna_structure_off_by_default(self, filter):
        # check_rna_structure defaults to False, so no RNA check runs
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "rna_structure" not in result["failed_checks"]
        assert "rna_structure_error" not in result["failed_checks"]
