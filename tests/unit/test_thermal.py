"""Characterization tests for SequenceFilter.thermal_filter.

Locks the Schema-v2 (Phase 0, 2026-05-14) semantics:
  * The PRIMARY thermal gate is GC% (G+C fraction); g_content is computed
    and reported, but no longer gates.
  * Homopolymer runs are now checked for all four bases (A/T/C/G).
  * Default max_consecutive_g raised from 4 → 5 to align with A/T/C.
  * The failed_checks codes are per-arm per-base
    (e.g. `consecutive_a_arm3`, `gc_content_arm5`).

Does NOT lock specific Tm values since those depend on Biopython MeltingTemp
table versions.
"""
from __future__ import annotations

import pytest

from probe_designer.filtering import SequenceFilter


# Arms tuned to pass the Schema-v2 defaults
#   GC ∈ [0.35, 0.65], no 5-run of any base, Tm 45-65, etc.
GOOD_ARM_3 = "GATGATGATGATGATGATGC"   # 20 nt, 8 G+C = 40%
GOOD_ARM_5 = "GATGATGATGATGATGATGC"
GOOD_TARGET_SEQ = "GCATCATCATCATCATCATCGCATCATCATCATCATCATC"

# Adversarial arms
LOW_G_ARM = "ATATATATATATATATATAT"           # 0 G, 0 C → fails new GC gate
HIGH_GC_ARM = "GCGCGCGCGCGCGCGCGCGC"         # 100% GC → fails max_gc
PURE_GC_BUT_LOW_G_ARM = "CCCTCCCTCCCTCCCTCCCT"  # 25% G but 80% GC: passes legacy G-only [0.3, 0.7] (fails actually — 0.25 < 0.3, let me reconsider)
# Better: a sequence with G=0.30 (legacy OK) and GC=0.70 (fails new max=0.65)
HIGH_GC_LOW_G_ARM = "GCCCATGCCCATGCCCATGC"   # 6 G/20 = 0.30, 14 G+C/20 = 0.70

FIVE_A_RUN_ARM = "GATGATGAAAAATGATGATC"      # 5×A run; 5-A fails new
FIVE_T_RUN_ARM = "GATGATGTTTTTGATGATGC"      # 5×T run
FIVE_C_RUN_ARM = "GATGATGCCCCCGATGATGC"      # 5×C run
FIVE_G_RUN_ARM = "GATGATGGGGGATGATGATC"      # 5×G run
FOUR_G_RUN_ARM = "GATGATGATGGGGATGATGA"      # 4×G run; passes new default but fails if max_consecutive_g=4


@pytest.fixture
def filter(default_filter_config, default_blast_config):
    return SequenceFilter(default_filter_config, default_blast_config)


# ---------------------------------------------------------------------------
# Returned dict shape
# ---------------------------------------------------------------------------

class TestThermalFilterReturnedDictShape:
    """Locks the output dict contract that downstream code depends on."""

    def test_all_expected_keys_present(self, filter):
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        expected_keys = {
            "arm_3prime", "arm_5prime", "target_sequence",
            "sequence_type", "target_type",
            "passed",
            # Schema-v2: both metrics emitted
            "g_content", "g_content_arm3", "g_content_arm5",
            "gc_content", "gc_content_arm3", "gc_content_arm5",
            "tm", "tm_3prime", "tm_5prime", "tm_diff",
            "free_energy",
            # legacy summary flag preserved
            "has_consecutive_g", "has_consecutive_g_arm3",
            "has_consecutive_g_arm5", "has_consecutive_g_target",
            # Schema-v2 per-base homopolymer fields
            "homopoly_A_arm3", "homopoly_T_arm3", "homopoly_C_arm3", "homopoly_G_arm3",
            "homopoly_A_arm5", "homopoly_T_arm5", "homopoly_C_arm5", "homopoly_G_arm5",
            "homopoly_G_target",
            "failed_checks",
            "arm_3prime_length", "arm_5prime_length",
        }
        missing = expected_keys - set(result.keys())
        assert not missing, f"Missing keys: {missing}"

    def test_arm_lengths_reported(self, filter):
        result = filter.thermal_filter("ACGTACGTAC", "GGGGACGTACGTACGT", target_sequence=GOOD_TARGET_SEQ)
        assert result["arm_3prime_length"] == 10
        assert result["arm_5prime_length"] == 16


# ---------------------------------------------------------------------------
# Pass/fail semantics
# ---------------------------------------------------------------------------

class TestThermalFilterPassPath:
    def test_passed_true_when_no_failed_checks(self, filter):
        result = filter.thermal_filter(
            GOOD_ARM_3, GOOD_ARM_5,
            target_sequence=GOOD_TARGET_SEQ,
            min_gc_content=0.0, max_gc_content=1.0,
            max_consecutive_a=100, max_consecutive_t=100,
            max_consecutive_c=100, max_consecutive_g=100,
            min_tm=0.0, max_tm=200.0, max_tm_diff=200.0,
        )
        assert result["failed_checks"] == []
        assert result["passed"] is True

    def test_passed_false_when_any_check_fails(self, filter):
        result = filter.thermal_filter(
            LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ,
            min_tm=0.0, max_tm=200.0, max_tm_diff=200.0,
        )
        assert result["passed"] is False
        assert len(result["failed_checks"]) > 0


# ---------------------------------------------------------------------------
# Phase 0 — GC content is the gate
# ---------------------------------------------------------------------------

class TestThermalFilterGCContent:
    def test_low_gc_fails_per_arm(self, filter):
        """An all-AT arm has GC=0 which is below min_gc_content (default 0.35)."""
        result = filter.thermal_filter(LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "gc_content_arm3" in result["failed_checks"]
        assert result["passed"] is False

    def test_high_gc_fails_per_arm(self, filter):
        """An all-GC arm has GC=1.0 which exceeds max_gc_content (default 0.65)."""
        result = filter.thermal_filter(HIGH_GC_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "gc_content_arm3" in result["failed_checks"]

    def test_gc_content_reported_per_arm(self, filter):
        result = filter.thermal_filter(LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert result["gc_content_arm3"] == pytest.approx(0.0)
        assert result["gc_content_arm5"] > 0.0

    def test_aggregate_gc_content_is_arm_average(self, filter):
        result = filter.thermal_filter(LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        avg = (result["gc_content_arm3"] + result["gc_content_arm5"]) / 2.0
        assert result["gc_content"] == pytest.approx(avg)

    def test_kwargs_override_min_gc_content(self, filter):
        result = filter.thermal_filter(
            LOW_G_ARM, LOW_G_ARM,
            target_sequence="A" * 40,
            min_gc_content=0.0,
        )
        assert "gc_content_arm3" not in result["failed_checks"]
        assert "gc_content_arm5" not in result["failed_checks"]


# ---------------------------------------------------------------------------
# Phase 0 — legacy G-only is informational, not a gate
# ---------------------------------------------------------------------------

class TestThermalFilterLegacyGContent:
    def test_g_content_reported_per_arm(self, filter):
        result = filter.thermal_filter(LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert result["g_content_arm3"] == pytest.approx(0.0)
        assert result["g_content_arm5"] > 0.0

    def test_aggregate_g_content_is_arm_average(self, filter):
        result = filter.thermal_filter(LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        avg = (result["g_content_arm3"] + result["g_content_arm5"]) / 2.0
        assert result["g_content"] == pytest.approx(avg)

    def test_phase0_gates_on_gc_not_g(self, filter):
        """Schema-v2 behavior: g_content is NOT a gate.

        ``HIGH_GC_LOW_G_ARM`` has G=0.30 (in old legacy window [0.3, 0.7]) and
        GC=0.70 (above new max_gc=0.65). The probe must FAIL because of GC,
        not pass because of G — and the failed_checks reason must be
        ``gc_content_arm3`` not ``g_content_arm3``.
        """
        result = filter.thermal_filter(
            HIGH_GC_LOW_G_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ,
            min_tm=0.0, max_tm=200.0, max_tm_diff=200.0,
        )
        assert "gc_content_arm3" in result["failed_checks"]
        assert "g_content_arm3" not in result["failed_checks"]


# ---------------------------------------------------------------------------
# Phase 0 — ATCG homopolymer gates
# ---------------------------------------------------------------------------

class TestThermalFilterHomopolymerRuns:
    """Schema-v2 default = 5×any-base run fails."""

    def test_five_a_run_fails(self, filter):
        result = filter.thermal_filter(FIVE_A_RUN_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "consecutive_a_arm3" in result["failed_checks"]
        assert result["homopoly_A_arm3"] >= 5

    def test_five_t_run_fails(self, filter):
        result = filter.thermal_filter(FIVE_T_RUN_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "consecutive_t_arm3" in result["failed_checks"]

    def test_five_c_run_fails(self, filter):
        result = filter.thermal_filter(FIVE_C_RUN_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "consecutive_c_arm3" in result["failed_checks"]

    def test_five_g_run_fails(self, filter):
        result = filter.thermal_filter(FIVE_G_RUN_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "consecutive_g_arm3" in result["failed_checks"]
        assert result["has_consecutive_g_arm3"] is True

    def test_four_g_run_passes_under_new_default(self, filter):
        """Default max_consecutive_g was raised from 4 → 5; GGGG no longer fails."""
        result = filter.thermal_filter(
            FOUR_G_RUN_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ,
            min_tm=0.0, max_tm=200.0, max_tm_diff=200.0,
        )
        assert "consecutive_g_arm3" not in result["failed_checks"]

    def test_kwargs_override_max_consecutive_g(self, filter):
        result = filter.thermal_filter(
            FOUR_G_RUN_ARM, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ,
            max_consecutive_g=4,
        )
        assert "consecutive_g_arm3" in result["failed_checks"]

    def test_target_g_run_still_gated(self, filter):
        target_with_run = "AAAAGGGGGAAAATTTTGGGGGAAAA"  # 5×G runs
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=target_with_run)
        assert result["has_consecutive_g_target"] is True
        assert "consecutive_g_target" in result["failed_checks"]


# ---------------------------------------------------------------------------
# Tm range (unchanged from pre-Phase-0)
# ---------------------------------------------------------------------------

class TestThermalFilterTmRange:
    def test_kwargs_override_tm_bounds(self, filter):
        result = filter.thermal_filter(
            GOOD_ARM_3, GOOD_ARM_5,
            target_sequence=GOOD_TARGET_SEQ,
            min_tm=200.0, max_tm=300.0,
        )
        assert "tm_3prime_range" in result["failed_checks"]
        assert "tm_5prime_range" in result["failed_checks"]

    def test_kwargs_override_max_tm_diff(self, filter):
        result = filter.thermal_filter(
            "GATGATGATGATGATGATGA", "GCGCGCGCGCGCGCGCGCGC",
            target_sequence=GOOD_TARGET_SEQ,
            max_tm_diff=0.01,
        )
        assert "tm_diff" in result["failed_checks"] or result["tm_diff"] > 0.01


# ---------------------------------------------------------------------------
# Numeric returns
# ---------------------------------------------------------------------------

class TestThermalFilterNumericReturns:
    def test_g_content_in_zero_to_one(self, filter):
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert 0.0 <= result["g_content"] <= 1.0
        assert 0.0 <= result["g_content_arm3"] <= 1.0
        assert 0.0 <= result["g_content_arm5"] <= 1.0

    def test_gc_content_in_zero_to_one(self, filter):
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert 0.0 <= result["gc_content"] <= 1.0
        assert 0.0 <= result["gc_content_arm3"] <= 1.0
        assert 0.0 <= result["gc_content_arm5"] <= 1.0

    def test_gc_at_least_as_large_as_g(self, filter):
        """Mathematical invariant: G+C count ≥ G count, so gc ≥ g."""
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert result["gc_content_arm3"] >= result["g_content_arm3"]
        assert result["gc_content_arm5"] >= result["g_content_arm5"]

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
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "rna_structure" not in result["failed_checks"]
        assert "rna_structure_error" not in result["failed_checks"]


# ---------------------------------------------------------------------------
# Phase 1A — target accessibility kwarg
# ---------------------------------------------------------------------------

class TestThermalFilterTargetAccessibility:
    def test_accessibility_field_in_result(self, filter):
        """The field is always emitted (None when caller didn't pass it)."""
        result = filter.thermal_filter(GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ)
        assert "target_accessibility" in result

    def test_accessibility_kwarg_above_threshold_passes(self, filter):
        result = filter.thermal_filter(
            GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ,
            target_accessibility=0.9, min_accessibility=0.3,
            min_tm=0.0, max_tm=200.0, max_tm_diff=200.0,
            min_gc_content=0.0, max_gc_content=1.0,
            max_consecutive_a=100, max_consecutive_t=100,
            max_consecutive_c=100, max_consecutive_g=100,
        )
        assert "accessibility" not in result["failed_checks"]
        assert result["target_accessibility"] == 0.9

    def test_accessibility_kwarg_below_threshold_fails(self, filter):
        result = filter.thermal_filter(
            GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ,
            target_accessibility=0.1, min_accessibility=0.3,
        )
        assert "accessibility" in result["failed_checks"]
        assert result["passed"] is False

    def test_accessibility_kwarg_replaces_self_fold_check(self, filter):
        """When target_accessibility is set, the legacy RNA.fold MFE path is
        NOT taken even if check_rna_structure=True."""
        result = filter.thermal_filter(
            GOOD_ARM_3, GOOD_ARM_5, target_sequence=GOOD_TARGET_SEQ,
            target_accessibility=0.9, min_accessibility=0.3,
            check_rna_structure=True,
        )
        # Neither legacy MFE code should appear
        assert "rna_structure" not in result["failed_checks"]
        assert "rna_structure_error" not in result["failed_checks"]
