"""Ligation geometry read off a BLAST alignment.

The gate these helpers replace was a coverage ratio, which cannot tell a
trimmed terminal mismatch (harmless) from a mismatch sitting on the ligation
junction (a probe that anneals and never circularises). Both score 39/40.
"""
from __future__ import annotations

import pytest

from probe_designer.blast.geometry import (
    DEFAULT_VICINITY_N,
    binding_hit_frame,
    frame_is_bindable,
    hit_is_bindable,
    hsp_junction_verdict,
    hsp_spans_junction,
    junction_offset,
    midpoint_junction_offset,
)

QUERY = "A" * 40
NICK = midpoint_junction_offset(40)          # 20


def _midline(*, mismatch_at=None, length=40, begin=0):
    line = ["|"] * length
    if mismatch_at is not None:
        line[mismatch_at - begin] = " "
    return "".join(line)


class TestJunctionVerdict:
    def test_clamp_default_is_the_splintr_width(self):
        assert DEFAULT_VICINITY_N == 3

    def test_fully_paired_alignment_spans_the_junction(self):
        assert hsp_spans_junction(_midline(), QUERY, 1, NICK)

    @pytest.mark.parametrize("pos", [17, 18, 19, 20, 21, 22])
    def test_any_mismatch_inside_the_clamp_fails(self, pos):
        assert not hsp_spans_junction(_midline(mismatch_at=pos), QUERY, 1, NICK)

    @pytest.mark.parametrize("pos", [0, 16, 23, 39])
    def test_mismatch_outside_the_clamp_passes(self, pos):
        assert hsp_spans_junction(_midline(mismatch_at=pos), QUERY, 1, NICK)

    def test_alignment_ending_before_the_junction_fails(self):
        # A one-arm hit: real evidence that this HSP cannot be ligated on.
        assert hsp_junction_verdict(_midline(length=15), QUERY[:15], 1, NICK) is False

    def test_five_prime_trim_shifts_the_coordinates(self):
        # BLAST trimmed 2 nt off the query start; the nick is still paired.
        assert hsp_spans_junction(_midline(length=38, begin=2), QUERY[2:], 3, NICK)

    def test_subject_insertion_inside_the_clamp_fails(self):
        # A gap column in the query means the subject carries an extra base
        # between two otherwise-paired probe positions: not contiguous.
        query = QUERY[:20] + "-" + QUERY[20:]
        midline = "|" * 20 + " " + "|" * 20
        assert not hsp_spans_junction(midline, query, 1, NICK)

    def test_no_midline_is_no_evidence(self):
        assert hsp_junction_verdict("", QUERY, 1, NICK) is None

    def test_ragged_midline_is_no_evidence(self):
        assert hsp_junction_verdict("|" * 10, QUERY, 1, NICK) is None

    def test_clamp_running_off_the_query_start_is_no_evidence(self):
        assert hsp_junction_verdict(_midline(), QUERY, 1, 2) is None

    def test_zero_vicinity_clamps_only_the_flanking_bases(self):
        assert hsp_spans_junction(
            _midline(mismatch_at=18), QUERY, 1, NICK, vicinity_n=0
        )
        assert not hsp_spans_junction(
            _midline(mismatch_at=19), QUERY, 1, NICK, vicinity_n=0
        )

    def test_strict_wrapper_treats_no_evidence_as_a_failure(self):
        assert hsp_junction_verdict("", QUERY, 1, NICK) is None
        assert not hsp_spans_junction("", QUERY, 1, NICK)


class TestBindingStrand:
    def test_arm_antisense_chemistries_bind_the_minus_strand(self):
        assert binding_hit_frame("dRNA") == -1
        assert binding_hit_frame("iLock") == -1

    def test_sense_chemistry_binds_the_plus_strand(self):
        # cDNA arms are mRNA-sense and template off the cDNA, so the transcript
        # that matters is the one carrying the query verbatim.
        assert binding_hit_frame("cDNA") == 1

    def test_antisense_transcript_is_not_bindable(self):
        assert not hit_is_bindable(1, "dRNA")
        assert not hit_is_bindable(-1, "cDNA")

    def test_real_hits_are_bindable(self):
        assert hit_is_bindable(-1, "dRNA")
        assert hit_is_bindable(1, "cDNA")

    @pytest.mark.parametrize("frame", [None, 0, 2, ""])
    def test_unreported_frame_stays_bindable(self, frame):
        # A missing measurement must never be the thing that clears a probe.
        assert frame_is_bindable(frame, -1)


class TestJunctionOffset:
    def test_plain_chemistry_nicks_after_arm3(self):
        assert junction_offset("A" * 20, "dRNA") == 20

    def test_ilock_nick_follows_fen1_cleavage(self):
        # FEN1 takes arm5[0] with the flap, so the whole arm3 is ligated.
        assert junction_offset("A" * 20, "iLock") == 20

    @pytest.mark.parametrize("chemistry", ["dRNA", "cDNA", "iLock"])
    def test_offset_lands_on_the_first_ligated_base_of_arm5(self, chemistry):
        from probe_designer.blast.gene_aware import extract_binding_query

        arm3, arm5 = "T" * 19 + "C", "G" + "A" * 19
        query = extract_binding_query(arm5, arm3, chemistry)
        first_ligated = arm5[1] if chemistry == "iLock" else arm5[0]
        assert query[junction_offset(arm3, chemistry)] == first_ligated

    def test_midpoint_matches_the_isoform_walk_convention(self):
        # arm_3prime is the first length // 2 nt of the probe.
        assert midpoint_junction_offset(40) == 20
        assert midpoint_junction_offset(41) == 20
