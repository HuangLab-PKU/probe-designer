"""Unit tests for ``probe_designer.qc.dimer_ascii.parse_primer3_dimer_pairing``.

The parser ingests primer3-py's 4-line ASCII dimer (SEQ_top/SEQ_mid/STR_mid/STR_bot
prefixed by ``SEQ\\t`` or ``STR\\t``) and returns base-pair information in
**canonical 5'→3' coordinates** for both strands. v2 cross-lig uses it to
verify the ligation-junction geometry (positions of A's real 3'-OH and
5'-P being paired on adjacent B positions in B's 5'→3').

Test approach: hand-construct ASCII fixtures matching primer3's exact format
(no live primer3 dependency in this layer) so we test the parser in
isolation.
"""
from __future__ import annotations

import pytest

from probe_designer.qc.dimer_ascii import parse_primer3_dimer_pairing


# ----------------------------------------------------------------------
# Fixtures — primer3 ASCII strings constructed by hand to match the
# exact format calc_heterodimer emits (verified empirically).
# ----------------------------------------------------------------------


# Perfect 8-bp blunt duplex:
#   seq_a = "GGGAAACC" (5'→3', 8 nt)
#   seq_b = "GGTTTCCC" (5'→3', 8 nt) — RC of seq_a, so they pair perfectly
# primer3 displays B in 3'→5' direction left-to-right, which is reverse(seq_b) = "CCCTTTGG"
PERFECT_8BP = "SEQ\t\nSEQ\tGGGAAACC\nSTR\tCCCTTTGG\nSTR\t\n"


# Internal bulge fixture: A has a 2-nt unpaired bulge in the middle; B is straight.
#   col:  0 1 2 3 4 5 6 7 8 9
#   topA: . . . a a . . . . .   ← A's unpaired bulge at cols 3-4
#   midA: G G G . . C C C C C   ← A's 8 paired bases
#   midB: C C C . . G G G G G   ← B's 8 paired bases
#   botB: . . . . . . . . . .   ← (empty — B has no unpaired bases in this fixture)
# A total length = 10 (8 paired + 2 unpaired bulge); B total length = 8 (all paired).
BULGED_DUPLEX = (
    "SEQ\t   AA     \n"          # top: 3 spaces, "AA", 5 spaces
    "SEQ\tGGG  CCCCC\n"          # mid: GGG + 2 spaces + CCCCC
    "STR\tCCC  GGGGG\n"          # mid: CCC + 2 spaces + GGGGG
    "STR\t          \n"          # bot: all spaces
)


# Empty / no-dimer fixture — primer3 sometimes emits this shape when no stable dimer found
EMPTY_DIMER = "SEQ\t\nSEQ\t\nSTR\t\nSTR\t\n"


# ----------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------


def test_perfect_8bp_lengths():
    r = parse_primer3_dimer_pairing(PERFECT_8BP)
    assert r["a_len"] == 8
    assert r["b_len"] == 8


def test_perfect_8bp_all_paired():
    r = parse_primer3_dimer_pairing(PERFECT_8BP)
    # All 8 columns should be paired (mid lines both have letters in cols 0-7)
    assert r["paired_columns"] == set(range(8))


def test_perfect_8bp_pair_positions_5to3_canonical():
    """A's pos i should pair with B's pos (7-i) in 5'→3' canonical coords.

    Because B is displayed 3'→5' (= reverse of B's 5'→3'), B's display index 0
    corresponds to B's 5'→3' position (b_len - 1) = 7. So A[0] pairs with B[7].
    A[7] pairs with B[0].
    """
    r = parse_primer3_dimer_pairing(PERFECT_8BP)
    pairs = r["pair_positions_a_to_b"]
    expected = {i: 7 - i for i in range(8)}
    assert pairs == expected


def test_perfect_8bp_ligation_geometry():
    """If we treat seq_a as a rotated ``arm3 + arm5`` with arm3 = first 4, arm5 = last 4,
    then A's pos 3 (= real 3'-OH) and pos 4 (= real 5'-P) should pair with adjacent
    B positions, with b_5to3 decreasing by 1.
    """
    r = parse_primer3_dimer_pairing(PERFECT_8BP)
    pairs = r["pair_positions_a_to_b"]
    assert pairs[3] == 4  # A position 3 → B position 4
    assert pairs[4] == 3  # A position 4 → B position 3
    # Adjacency check: pairs[3] - pairs[4] should == 1 (B's 5'→3' decreases as A advances)
    assert pairs[3] - pairs[4] == 1


def test_bulged_duplex_lengths():
    r = parse_primer3_dimer_pairing(BULGED_DUPLEX)
    assert r["a_len"] == 10   # 8 paired + 2 unpaired bulge
    assert r["b_len"] == 8    # 8 paired


def test_bulged_duplex_paired_columns():
    r = parse_primer3_dimer_pairing(BULGED_DUPLEX)
    # Paired cols are where mid-A AND mid-B both have a letter:
    # cols 0,1,2 (GGG/CCC) and cols 5,6,7,8,9 (CCCCC/GGGGG)
    assert r["paired_columns"] == {0, 1, 2, 5, 6, 7, 8, 9}


def test_bulged_duplex_pair_positions():
    r = parse_primer3_dimer_pairing(BULGED_DUPLEX)
    pairs = r["pair_positions_a_to_b"]
    # A's indices in 5'→3' (walking cols 0-9, picking mid first then top):
    #   col 0,1,2 → A[0,1,2] paired (mid)
    #   col 3,4   → A[3,4] from top (unpaired bulge)
    #   col 5-9   → A[5,6,7,8,9] paired (mid)
    # B's display indices walking cols 0-9 (mid only, no bot):
    #   col 0,1,2 → B_disp 0,1,2; col 5-9 → B_disp 3,4,5,6,7
    # Convert B_disp to canonical: b_5to3 = (b_len - 1 - b_disp) = 7 - b_disp
    # So B_disp 0,1,2 → b_5to3 7,6,5; B_disp 3-7 → b_5to3 4,3,2,1,0
    assert pairs[0] == 7   # A[0]=G ↔ B at col 0, b_5to3 = 7
    assert pairs[1] == 6
    assert pairs[2] == 5
    # A[3], A[4] are unpaired (top line) — not in pairs
    assert 3 not in pairs
    assert 4 not in pairs
    assert pairs[5] == 4
    assert pairs[6] == 3
    assert pairs[7] == 2
    assert pairs[8] == 1
    assert pairs[9] == 0


def test_empty_dimer_returns_empty():
    r = parse_primer3_dimer_pairing(EMPTY_DIMER)
    assert r["a_len"] == 0
    assert r["b_len"] == 0
    assert r["paired_columns"] == set()
    assert r["pair_positions_a_to_b"] == {}


def test_malformed_structure_returns_empty():
    """Garbage input should not crash; returns empty structure."""
    r = parse_primer3_dimer_pairing("not a primer3 ascii structure")
    assert r["a_len"] == 0
    assert r["b_len"] == 0


def test_lookup_maps_consistent():
    """``a_idx_to_col`` and ``col_to_a_idx`` should be inverses of each other."""
    r = parse_primer3_dimer_pairing(PERFECT_8BP)
    for a_idx, col in r["a_idx_to_col"].items():
        assert r["col_to_a_idx"][col] == a_idx
    for b_idx, col in r["b_idx_to_col"].items():
        assert r["col_to_b_idx"][col] == b_idx


def test_dashes_treated_as_gaps_not_bases():
    """primer3 sometimes emits ``-`` in lines to indicate alignment gaps. Those
    are NOT base positions and must not advance the running index."""
    # A occupies cols 0-2 ("AAA"); col 3 is a dash (gap, no A base); cols 4-6 ("CCC")
    ascii_with_dashes = (
        "SEQ\t          \n"
        "SEQ\tAAA-CCC   \n"
        "STR\tTTT GGG   \n"
        "STR\t          \n"
    )
    r = parse_primer3_dimer_pairing(ascii_with_dashes)
    assert r["a_len"] == 6   # 3 + 3 letters, dash skipped
    assert r["b_len"] == 6   # 3 + 3 letters
