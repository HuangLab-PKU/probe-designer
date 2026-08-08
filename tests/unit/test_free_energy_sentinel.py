"""`free_energy` distinguishes "not computed" from "computed as zero" (audit R8).

The field was initialised to ``0.0`` and only overwritten inside the legacy
``check_rna_structure`` branch, which runs solely when no pre-computed
accessibility was supplied *and* the target is RNA. So on the cDNA path, and on
any run using the RNAplfold accessibility gate, ``free_energy`` stayed ``0.0``
— indistinguishable from a genuinely unstructured target, and enough to make
the mutation score's ``mfe_penalty = -free_energy`` silently inert.

``None`` is the honest sentinel: consumers can tell the two apart, and the
numeric behaviour of every existing consumer is unchanged (both ``0.0`` and
``None`` contribute nothing to a score).
"""
from __future__ import annotations

import pytest

from probe_designer.filtering import SequenceFilter
from probe_designer.scoring import compute_target_score

ARM_3 = "GATGATGATGATGATGATGC"
ARM_5 = "GATGATGATGATGATGATGC"
RNA_TARGET = "GCATCATCATCATCATCATCGCATCATCATCATCATCATC"


@pytest.fixture
def filt(default_filter_config, default_blast_config):
    return SequenceFilter(default_filter_config, default_blast_config)


class TestSentinel:
    def test_none_when_structure_check_is_off(self, filt):
        result = filt.thermal_filter(
            ARM_3, ARM_5, target_sequence=RNA_TARGET, check_rna_structure=False,
        )
        assert result["free_energy"] is None

    def test_none_when_accessibility_supersedes_the_mfe_path(self, filt):
        """A pre-computed accessibility replaces the self-fold check entirely."""
        result = filt.thermal_filter(
            ARM_3, ARM_5, target_sequence=RNA_TARGET,
            check_rna_structure=True, target_accessibility=0.8,
        )
        assert result["free_energy"] is None

    def test_none_on_a_dna_target_even_with_the_check_on(self, filt):
        """The MFE branch is RNA-only; cDNA never populated it."""
        result = filt.thermal_filter(
            ARM_3, ARM_5, target_sequence=RNA_TARGET,
            target_type="DNA", check_rna_structure=True,
        )
        assert result["free_energy"] is None

    def test_populated_on_the_rna_self_fold_path(self, filt):
        pytest.importorskip("RNA")
        result = filt.thermal_filter(
            ARM_3, ARM_5, target_sequence=RNA_TARGET,
            target_type="RNA", check_rna_structure=True,
        )
        assert isinstance(result["free_energy"], float)
        assert result["free_energy"] <= 0.0


class TestConsumersHandleNone:
    """None must be inert, not a crash — and score exactly as 0.0 used to."""

    def _site(self, **overrides):
        base = {
            "arm_5prime": "ATATATATAT", "arm_3prime": "ATATATATAT",
            "tm_5prime": 50.0, "tm_3prime": 50.0,
            "isoform_overlap_num": 1, "blast_alignments": [],
        }
        base.update(overrides)
        return base

    def test_scorer_treats_none_like_the_old_zero(self):
        assert compute_target_score(self._site(free_energy=None)) == pytest.approx(
            compute_target_score(self._site(free_energy=0.0))
        )

    def test_scorer_still_rewards_a_real_negative_mfe(self):
        neutral = compute_target_score(self._site(free_energy=None))
        structured = compute_target_score(self._site(free_energy=-5.0))
        assert structured > neutral

    def test_mutation_candidate_score_survives_none(self):
        """ext/mutation computed `-free_energy`, which raises on None."""
        from probe_designer.ext.mutation.probe import mfe_penalty_from

        assert mfe_penalty_from(None) == 0.0
        assert mfe_penalty_from(-4.0) == 4.0
