"""Tests for filtering/pairwise_duplex.py (Phase 1B, 2026-05-14)."""
from __future__ import annotations

import pytest

from probe_designer.filtering.pairwise_duplex import (
    DEFAULT_DG_THRESHOLD,
    DuplexHit,
    _parse_duplex_spans,
    predict_pairwise_duplex,
    screen_all_pairs,
)


REV_COMP_PAIR = (
    "GATGATGATGATGATGATGC",        # 20-mer
    "GCATCATCATCATCATCATC",        # exact RC of above
)
RANDOM_PAIR = (
    "ACGTACGTACGTACGTACGT",
    "TTTAAACCCGGGAAATTTCC",
)


class TestPredictPairwiseDuplex:
    def test_reverse_complement_pair_returns_strong_duplex(self):
        hit = predict_pairwise_duplex(*REV_COMP_PAIR)
        assert hit is not None
        assert hit.delta_g < -20.0, f"expected strong duplex, got {hit.delta_g}"
        assert "(" in hit.structure and ")" in hit.structure
        assert "&" in hit.structure

    def test_random_pair_returns_weak_or_none(self):
        hit = predict_pairwise_duplex(*RANDOM_PAIR)
        # Either no duplex or a weak one — let the orthogonality threshold
        # be the actual gate; here we just assert it's nothing dramatic.
        if hit is not None:
            assert hit.delta_g > -10.0, (
                f"random pair shouldn't bind strongly; got {hit.delta_g}"
            )

    def test_returns_none_on_empty_strand(self):
        assert predict_pairwise_duplex("", "ACGT") is None
        assert predict_pairwise_duplex("ACGT", "") is None

    def test_ids_propagate_to_hit(self):
        hit = predict_pairwise_duplex(
            *REV_COMP_PAIR, probe_a_id="ProbeA", probe_b_id="ProbeB",
        )
        assert hit.probe_a_id == "ProbeA"
        assert hit.probe_b_id == "ProbeB"

    def test_dna_letters_normalized_to_rna(self):
        """T→U substitution is transparent at the API level: same energy
        either way."""
        hit_dna = predict_pairwise_duplex(
            "GATGATGATGATGATGATGC", "GCATCATCATCATCATCATC",
        )
        hit_rna = predict_pairwise_duplex(
            "GAUGAUGAUGAUGAUGAUGC", "GCAUCAUCAUCAUCAUCAUC",
        )
        assert hit_dna is not None and hit_rna is not None
        assert hit_dna.delta_g == pytest.approx(hit_rna.delta_g)

    def test_temperature_change_weakens_duplex(self):
        """At higher temperature, duplex ΔG should be less negative
        (weaker pairing). This is the 'formamide approximation' knob."""
        hit_cold = predict_pairwise_duplex(*REV_COMP_PAIR, temperature=20.0)
        hit_hot = predict_pairwise_duplex(*REV_COMP_PAIR, temperature=60.0)
        assert hit_cold is not None and hit_hot is not None
        assert hit_hot.delta_g > hit_cold.delta_g, (
            f"hot ({hit_hot.delta_g}) should be > cold ({hit_cold.delta_g})"
        )


class TestParseDuplexSpans:
    def test_full_match_spans_both_strands(self):
        # structure with 20 brackets each side, i=20, j=1
        struct = "(" * 20 + "&" + ")" * 20
        span_a, span_b = _parse_duplex_spans(struct, i=20, j=1)
        assert span_a == (0, 20)
        assert span_b == (0, 20)

    def test_partial_match_in_middle(self):
        # 5-bp duplex; left brackets are 5, structure positions 6..10 on A
        # (1-indexed). i = 10. Right is 5 brackets, j = 1 on B.
        struct = "(((((&)))))"
        span_a, span_b = _parse_duplex_spans(struct, i=10, j=1)
        assert span_a == (5, 10)
        assert span_b == (0, 5)


class TestScreenAllPairs:
    def test_empty_input_returns_empty_list(self):
        assert screen_all_pairs([]) == []

    def test_single_probe_returns_empty(self):
        assert screen_all_pairs([("only", "ACGT")]) == []

    def test_flags_only_below_threshold(self):
        # 4 probes:
        # A and B are RC of each other → strong duplex
        # C and D are random → weak/no duplex
        probes = [
            ("A", "GATGATGATGATGATGATGC"),
            ("B", "GCATCATCATCATCATCATC"),
            ("C", "ACGTACGTACGTACGTACGT"),
            ("D", "AAAAATTTTGGGGCCCCTTT"),
        ]
        # default threshold -12 kcal/mol
        hits = screen_all_pairs(probes)
        # AB pair must be flagged
        ab_hits = [h for h in hits
                   if {h.probe_a_id, h.probe_b_id} == {"A", "B"}]
        assert len(ab_hits) == 1
        assert ab_hits[0].delta_g <= -12.0

    def test_threshold_kwarg_respected(self):
        probes = [
            ("A", "GATGATGATGATGATGATGC"),
            ("B", "GCATCATCATCATCATCATC"),
        ]
        # A very tight threshold (impossible to reach) returns no hits.
        assert screen_all_pairs(probes, dg_threshold=-9999.0) == []
        # A very loose threshold (everything passes) returns the pair.
        loose = screen_all_pairs(probes, dg_threshold=0.0)
        assert len(loose) == 1
