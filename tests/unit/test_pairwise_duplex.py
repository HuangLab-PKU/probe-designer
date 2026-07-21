"""Tests for filtering/pairwise_duplex.py — primer3 DNA backend (audit R5/P3).

Rewritten 2026-07-20: the backend moved from ViennaRNA ``duplexfold`` (RNA
Turner parameters — the wrong physics for a DNA:DNA padlock dimer, and the
source of the dot-bracket/T→U behaviour these tests used to pin) to primer3's
SantaLucia DNA nearest-neighbor model evaluated at the real buffer. ΔG
magnitudes are not comparable to the old RNA-parameter values, so assertions are
relational where possible.
"""
from __future__ import annotations

import pytest

pytest.importorskip("primer3")

from probe_designer.chemistry import ReactionConditions           # noqa: E402
from probe_designer.filtering.pairwise_duplex import (            # noqa: E402
    DEFAULT_DG_THRESHOLD,
    predict_pairwise_duplex,
    screen_all_pairs,
)

REV_COMP_PAIR = (
    "GATGATGATGATGATGATGC",        # 20-mer
    "GCATCATCATCATCATCATC",        # exact RC of above
)
RANDOM_PAIR = (
    "ACGTACGTACGTACGTACGT",
    "AAAGAAAGAAAGAAAGAAAG",
)


class TestPredictPairwiseDuplex:
    def test_reverse_complement_pair_returns_strong_duplex(self):
        hit = predict_pairwise_duplex(*REV_COMP_PAIR)
        assert hit is not None
        assert hit.delta_g < -10.0, f"expected strong duplex, got {hit.delta_g}"

    def test_random_pair_much_weaker_than_rc(self):
        rc = predict_pairwise_duplex(*REV_COMP_PAIR)
        rnd = predict_pairwise_duplex(*RANDOM_PAIR)
        rnd_dg = rnd.delta_g if rnd is not None else 0.0
        assert rc.delta_g < rnd_dg - 5.0

    def test_returns_none_on_empty_strand(self):
        assert predict_pairwise_duplex("", "ACGT") is None
        assert predict_pairwise_duplex("ACGT", "") is None

    def test_ids_propagate_to_hit(self):
        hit = predict_pairwise_duplex(
            *REV_COMP_PAIR, probe_a_id="ProbeA", probe_b_id="ProbeB",
        )
        assert hit.probe_a_id == "ProbeA"
        assert hit.probe_b_id == "ProbeB"

    def test_structure_is_primer3_ascii(self):
        hit = predict_pairwise_duplex(*REV_COMP_PAIR)
        assert "SEQ" in hit.structure or "STR" in hit.structure

    def test_spans_within_bounds_and_nonempty(self):
        a, b = REV_COMP_PAIR
        hit = predict_pairwise_duplex(a, b)
        (a0, a1), (b0, b1) = hit.span_a, hit.span_b
        assert 0 <= a0 < a1 <= len(a)
        assert 0 <= b0 < b1 <= len(b)


class TestBufferAwareness:
    """The screen now sees the real buffer (audit P5/R5), not a bare temperature."""

    def test_higher_temperature_weakens_duplex(self):
        cold = predict_pairwise_duplex(
            *REV_COMP_PAIR,
            reaction=ReactionConditions(lab_temp_c=25.0, formamide_pct=0.0),
        )
        hot = predict_pairwise_duplex(
            *REV_COMP_PAIR,
            reaction=ReactionConditions(lab_temp_c=70.0, formamide_pct=0.0),
        )
        assert hot.delta_g > cold.delta_g

    def test_formamide_raises_effective_temperature(self):
        no_fa = predict_pairwise_duplex(
            *REV_COMP_PAIR, reaction=ReactionConditions(formamide_pct=0.0))
        with_fa = predict_pairwise_duplex(
            *REV_COMP_PAIR, reaction=ReactionConditions(formamide_pct=20.0))
        # formamide -> hotter effective temperature -> weaker duplex
        assert with_fa.delta_g > no_fa.delta_g


class TestScreenAllPairs:
    def test_empty_input_returns_empty_list(self):
        assert screen_all_pairs([]) == []

    def test_single_probe_returns_empty(self):
        assert screen_all_pairs([("only", "ACGTACGTACGTACGTACGT")]) == []

    def test_flags_complementary_pair(self):
        probes = [
            ("A", REV_COMP_PAIR[0]),
            ("B", REV_COMP_PAIR[1]),
            ("C", "AAAGAAAGAAAGAAAGAAAG"),
        ]
        hits = screen_all_pairs(probes, dg_threshold=-10.0)
        ab_hits = [h for h in hits
                   if {h.probe_a_id, h.probe_b_id} == {"A", "B"}]
        assert len(ab_hits) == 1
        assert ab_hits[0].delta_g <= -10.0

    def test_threshold_kwarg_respected(self):
        probes = [("A", REV_COMP_PAIR[0]), ("B", REV_COMP_PAIR[1])]
        # Unreachable threshold -> nothing passes.
        assert screen_all_pairs(probes, dg_threshold=-9999.0) == []
        # Loose threshold -> the pair is returned.
        assert len(screen_all_pairs(probes, dg_threshold=0.0)) == 1

    def test_default_threshold_constant(self):
        assert DEFAULT_DG_THRESHOLD == -12.0
