"""Tests for the TCR clone-specificity screen: exact hosts, and near ones.

A CDR3 padlock's job is to tell one clone from hundreds of near-identical ones, and the
repertoire screen is the only direct test of that. Screening on the EXACT 40-mer answers a
narrower question than the ligase asks: two clones can differ by a single base far from the
nick, and there both arms still hybridise and the nick is still perfectly paired, so the
probe closes on the other clone too. Judged here the way the BLAST path judges an
off-target hit -- the clamp around the nick must pair, the rest may carry a few mismatches.
"""
from __future__ import annotations

from probe_designer.ext.tcr.repertoire import (
    RepertoireScreen,
    build_repertoire_index,
    is_clone_specific,
)

V = "CACAGTGGCTCAAGCTTCTGATGGCTCAAG"
J = "AACTGTGAGTCTGGTGCCTTGTCCAAAGAA"
OWN = "TGTGCCAGCAGCTTAGGGACAGGGACTTAT"
NICK_IN_CDR3 = 15
SIZE = 40


def chain_of(cdr3: str) -> str:
    return V + cdr3 + J


def footprint_of(cdr3: str) -> str:
    start = len(V) + NICK_IN_CDR3 - SIZE // 2
    return chain_of(cdr3)[start:start + SIZE]


def mutate(cdr3: str, *positions: int) -> str:
    out = list(cdr3)
    for p in positions:
        out[p] = "A" if out[p] != "A" else "C"
    return "".join(out)


def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


ONE_OFF = mutate(OWN, 5)          # 10 nt from the nick: outside the clamp
AT_NICK = mutate(OWN, NICK_IN_CDR3)   # the nick base itself
NEAR_NICK = mutate(OWN, NICK_IN_CDR3 + 2)
FIVE_OFF = mutate(OWN, 2, 4, 6, 8, 10)
OTHER = "TGTGCCAGCAGTGGGACTAGCGGGAGCTAC"


def screen(*cdr3s, **kw):
    return RepertoireScreen([(c, chain_of(c)) for c in cdr3s], **kw)


def test_exact_screen_cannot_see_a_one_base_difference():
    index = build_repertoire_index([(ONE_OFF, chain_of(ONE_OFF))], kmer=SIZE)
    assert is_clone_specific(footprint_of(OWN), index, own_cdr3=OWN)


def test_near_screen_rejects_a_site_one_base_off_another_clone():
    s = screen(OWN, ONE_OFF)
    hosts = s.hosts(footprint_of(OWN), nick=SIZE // 2)
    assert hosts == {OWN: 0, ONE_OFF: 1}   # hosts() lists the designed clone too
    assert not s.is_clone_specific(footprint_of(OWN), own_cdr3=OWN, nick=SIZE // 2)


def test_a_difference_inside_the_clamp_keeps_the_probe_specific():
    for cdr3 in (AT_NICK, NEAR_NICK):
        s = screen(OWN, cdr3)
        assert s.hosts(footprint_of(OWN), nick=SIZE // 2) == {OWN: 0}
        assert s.is_clone_specific(footprint_of(OWN), own_cdr3=OWN, nick=SIZE // 2)


def test_mismatches_beyond_the_budget_are_not_a_host():
    s = screen(OWN, FIVE_OFF, max_mismatches=4)
    assert s.hosts(footprint_of(OWN), nick=SIZE // 2) == {OWN: 0}


def test_zero_budget_reproduces_the_exact_screen():
    s = screen(OWN, ONE_OFF, max_mismatches=0)
    assert s.is_clone_specific(footprint_of(OWN), own_cdr3=OWN, nick=SIZE // 2)
    same = RepertoireScreen([(OWN, chain_of(OWN))], max_mismatches=0)
    assert same.hosts(footprint_of(OWN), nick=SIZE // 2) == {OWN: 0}


def test_one_clone_filed_twice_is_not_cross_reactivity():
    s = RepertoireScreen([(OWN, chain_of(OWN)), (OWN, chain_of(OWN)[2:])])
    assert s.is_clone_specific(footprint_of(OWN), own_cdr3=OWN, nick=SIZE // 2)


def test_an_unrelated_clone_is_not_a_host():
    s = screen(OWN, OTHER)
    assert s.hosts(footprint_of(OWN), nick=SIZE // 2) == {OWN: 0}


def test_a_chain_stored_on_the_other_strand_is_still_found():
    s = RepertoireScreen([(ONE_OFF, rc(chain_of(ONE_OFF)))])
    assert s.hosts(footprint_of(OWN), nick=SIZE // 2) == {ONE_OFF: 1}


def test_the_default_budget_and_clamp_are_the_ligase_grounded_ones():
    s = RepertoireScreen([])
    assert (s.clamp, s.max_mismatches, s.kmer) == (3, 4, 40)
    assert s.is_clone_specific(footprint_of(OWN), own_cdr3=OWN, nick=SIZE // 2)
