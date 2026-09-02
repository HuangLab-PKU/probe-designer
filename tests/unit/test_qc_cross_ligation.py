"""Unit tests for ``probe_designer.qc.cross_ligation`` — the v3 register scan.

The tests that matter here are the **adversarial** ones. v2.2 passed a suite of
perfect-match fixtures while being blind to any productive register that was
not also the globally most stable alignment, which is the normal case in a real
pool. ``test_productive_register_hidden_behind_*`` reproduces exactly that and
is the regression guard for the 2026-09-02 audit; a perfect-match fixture alone
proves almost nothing.

This module has NO bank dependency (the static check in
``test_no_bank_import_in_qc.py`` enforces that for the package).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from probe_designer.chemistry import ReactionConditions, reverse_complement as rc
from probe_designer.qc.cross_ligation import (
    DEFAULT_TM_THRESHOLD_C,
    DEFAULT_VICINITY_N,
    LigationDimer,
    ProbeForScreen,
    build_v2_geom,
    find_ligation_registers,
    junction_block,
    screen_cross_ligation_v2,
    screen_self_circularisation,
    write_dimer_report,
    write_self_circ_report,
    _build_geometry,
    _fray_trim,
    _register_window,
)


BB = "TCCCTACACGACGCTCTTCCG"
# Fillers with no complementarity to the arms below, and none to each other.
F1 = "GAGGAAGAGGAAGAGGAAGAGGA"
F2 = "CTCCTTCTCCTTCTCCTTCTCCT"
F3 = "AGAGGAAGAGGAAGAGGAAGAG"

ARM5_A = "ACGTGCAGTCATGGCTAACG"
ARM3_A = "TTGCCAGATCGTACGATCCA"
#: The rotated ligator: arm3 then arm5, with the ligation nick in the middle.
ROTATED_A = ARM3_A + ARM5_A
#: A splint perfectly templating A's whole 40 nt ligation footprint.
SPLINT40 = rc(ROTATED_A)
#: 12 nt each side of the nick — comfortably ligation-competent on its own.
CORE24 = SPLINT40[8:32]

PROBE_A = ProbeForScreen("A", "dRNA", ARM5_A, ARM3_A, ARM5_A + BB + ARM3_A, "GENE_A")


def _probe(pid: str, seq: str, target: str = "GENE_B") -> ProbeForScreen:
    """A splint-side probe carrying ``seq``; its own arms are incidental."""
    return ProbeForScreen(pid, "dRNA", seq[:20], seq[-20:], seq, target)


def _ligator_hits(probes, ligator_id="A"):
    _, dimers = screen_cross_ligation_v2(probes, include_self_pairs=False)
    return [
        d for d in dimers
        if d.seq_a_id == ligator_id and d.seq_b_id != ligator_id
    ]


def _is_flagged(probes, ligator_id="A") -> bool:
    return any(
        d.flagged_overall and d.a_can_ligate_on_b
        for d in _ligator_hits(probes, ligator_id)
    )


# ----------------------------------------------------------------------
# Geometry primitives
# ----------------------------------------------------------------------


def test_junction_block_is_last_of_arm3_then_first_of_arm5():
    block = junction_block(ARM3_A, ARM5_A, 3)
    assert len(block) == 2 * (3 + 1)
    assert block == ARM3_A[-4:] + ARM5_A[:4]


def test_junction_block_empty_when_an_arm_is_too_short():
    assert junction_block("AC", ARM5_A, 3) == ""
    assert junction_block(ARM3_A, "AC", 3) == ""


def test_find_registers_is_exact_on_a_constructed_splint():
    """The needle sits at a known offset, so ``j`` is known exactly."""
    geom = _build_geometry(PROBE_A, DEFAULT_VICINITY_N)
    splint = F1 + rc(geom.junction_block) + F2
    registers = find_ligation_registers(geom, splint, DEFAULT_VICINITY_N)
    assert registers == [len(F1) + DEFAULT_VICINITY_N]


def test_register_window_puts_the_nick_between_j_and_j_plus_one():
    geom = _build_geometry(PROBE_A, DEFAULT_VICINITY_N)
    splint = F1 + SPLINT40 + F2
    j = find_ligation_registers(geom, splint, DEFAULT_VICINITY_N)[0]
    lig_sub, c_sub, nick_idx = _register_window(
        geom.ligator, len(geom.arm3_effective), len(geom.arm5), splint, j,
    )
    # arm3's 3'-OH is the base just before the nick; arm5's 5'-P just after.
    assert lig_sub[nick_idx - 1] == geom.arm3_effective[-1]
    assert lig_sub[nick_idx] == geom.arm5[0]
    # And on a perfect splint the whole 40 nt pairs.
    assert all(a != b for a, b in zip(lig_sub, c_sub))  # complementary, never equal
    assert len(lig_sub) == len(ROTATED_A)


def test_fray_trim_stops_at_two_consecutive_mismatches():
    lig = "AAAA" + "GG" + "AAAA"
    #     paired    mism   paired-again-but-unreachable
    c = "TTTT" + "AA" + "TTTT"
    left, right = _fray_trim(lig, c, nick_idx=2)
    assert (left, right) == (0, 4), "helix must stop at the 2-mismatch run"


def test_fray_trim_walks_through_a_single_mismatch():
    lig = "AAAA" + "G" + "AAAA"
    c = "TTTT" + "A" + "TTTT"
    left, right = _fray_trim(lig, c, nick_idx=2)
    assert (left, right) == (0, 9)


# ----------------------------------------------------------------------
# Core screen behaviour
# ----------------------------------------------------------------------


def test_empty_probe_list_returns_empty():
    assert screen_cross_ligation_v2([]) == ([], [])


def test_single_probe_has_no_cross_pair():
    tier1, dimers = screen_cross_ligation_v2([PROBE_A], include_self_pairs=False)
    assert tier1 == [] and dimers == []


def test_safe_pair_produces_no_register_at_all():
    safe_a = ProbeForScreen(
        "safe_A", "dRNA", "GCATAGCAGCAGCAGCATAG", "TGTGTGTGCACGCACGCATG",
        "GCATAGCAGCAGCAGCATAG" + BB + "TGTGTGTGCACGCACGCATG", "G1",
    )
    safe_b = ProbeForScreen(
        "safe_B", "dRNA", "AAGAAGAAGAAGAAGAAGAA", "TTCTTCTTCTTCTTCTTCTT",
        "AAGAAGAAGAAGAAGAAGAA" + BB + "TTCTTCTTCTTCTTCTTCTT", "G2",
    )
    tier1, dimers = screen_cross_ligation_v2(
        [safe_a, safe_b], include_self_pairs=False,
    )
    assert tier1 == [] and dimers == []


def test_perfect_splint_is_flagged():
    assert _is_flagged([PROBE_A, _probe("B", F1 + SPLINT40 + F2)])


def test_flagged_dimer_has_nick_adjacent_geometry():
    """``b_3oh - b_5p == +1``: both arms anneal antiparallel to the splint."""
    for d in _ligator_hits([PROBE_A, _probe("B", F1 + SPLINT40 + F2)]):
        assert d.b_3oh_pos is not None and d.b_5p_pos is not None
        assert d.b_3oh_pos - d.b_5p_pos == 1
        assert d.b_5p_pos == d.nick_pos_on_b


def test_partial_core_alone_is_flagged():
    """24 bp across the nick is ligation-competent without the outer 16."""
    assert _is_flagged([PROBE_A, _probe("B", F1 + CORE24 + F2)])


# ----------------------------------------------------------------------
# The regression that v2.2 failed — a productive register that is not argmax
# ----------------------------------------------------------------------


@pytest.mark.parametrize("decoy_len", [12, 16, 20])
def test_productive_register_hidden_behind_a_more_stable_decoy(decoy_len):
    """The audit's killer case: same productive site, plus a better non-productive one.

    ``CORE24`` is untouched and still ligation-competent. Elsewhere on the SAME
    splint sits a longer perfect match for arm3 alone, which any argmax-based
    method reports instead — v2.2 called all three of these safe. Since every
    probe in a pool shares a backbone, a competing >=12 nt match is the norm,
    not a contrivance.
    """
    decoy = rc(ARM3_A)[:decoy_len]
    splint = F1 + decoy + F2 + CORE24 + F3
    assert _is_flagged([PROBE_A, _probe("B", splint)])


def test_decoy_without_a_productive_register_stays_clean():
    """The mirror image: a strong arm3 match but no nick — must not fire."""
    splint = F1 + rc(ARM3_A) + F2 + F3
    assert not _is_flagged([PROBE_A, _probe("B", splint)])
    assert _ligator_hits([PROBE_A, _probe("B", splint)]) == []


def test_register_tm_beats_either_arm_alone():
    """The productive complex is one duplex across the nick, not two halves.

    Scoring the arms separately is what made a real site look weak enough to
    lose its argmax; the joint duplex is the quantity that decides occupancy.
    """
    d = max(
        _ligator_hits([PROBE_A, _probe("B", F1 + CORE24 + F2)]),
        key=lambda x: x.overall_tm_c,
    )
    assert d.overall_tm_c > d.arm3_tm_c
    assert d.overall_tm_c > d.arm5_tm_c


def test_limiting_arm_tm_is_the_min_not_the_max():
    """Ligation needs BOTH arms annealed; v2.2 aggregated with ``max``."""
    for d in _ligator_hits([PROBE_A, _probe("B", F1 + CORE24 + F2)]):
        assert d.limiting_arm_tm_c == min(d.arm3_tm_c, d.arm5_tm_c)


# ----------------------------------------------------------------------
# Self-pairs and self-circularisation — unscreened before v3
# ----------------------------------------------------------------------


def test_self_pair_is_screened_when_requested():
    """A probe templating another copy of itself: ``combinations`` excluded it."""
    geom = _build_geometry(PROBE_A, DEFAULT_VICINITY_N)
    # Give A a backbone carrying the RC of its own junction block, so a second
    # copy of A can act as its splint.
    seq = ARM5_A + "TCCCTA" + rc(geom.junction_block) + "CGCTCTTCCG" + ARM3_A
    self_lig = ProbeForScreen("S", "dRNA", ARM5_A, ARM3_A, seq, "GENE_S")

    _, without = screen_cross_ligation_v2([self_lig], include_self_pairs=False)
    _, with_self = screen_cross_ligation_v2([self_lig], include_self_pairs=True)
    assert without == []
    assert [d for d in with_self if d.is_self_pair], "self-pair must be reachable"


def test_self_circularisation_is_detected():
    """One molecule presenting its own 3'-OH and 5'-P — no splint at all."""
    geom = _build_geometry(PROBE_A, DEFAULT_VICINITY_N)
    spacer = "TCCCTACACGACGCTCTTCCGATCTACGT"   # long enough to clear both arms
    seq = ARM5_A + spacer + rc(geom.junction_block) + spacer + ARM3_A
    folder = ProbeForScreen("SC", "dRNA", ARM5_A, ARM3_A, seq, "GENE_SC")

    hits = screen_self_circularisation([folder, PROBE_A])
    assert [h.probe_id for h in hits] == ["SC"]
    assert hits[0].junction_run_nt >= 2 * (DEFAULT_VICINITY_N + 1)
    assert "^" in hits[0].alignment


def test_self_circularisation_rejects_registers_overlapping_the_arms():
    """A base cannot pair with itself, so an overlapping register is not a fold."""
    geom = _build_geometry(PROBE_A, DEFAULT_VICINITY_N)
    # Same block, but with no room between the arms — the 40 nt register would
    # have to reuse arm bases as its own template.
    seq = ARM5_A + "TCCCTA" + rc(geom.junction_block) + "CGCTCTTCCG" + ARM3_A
    cramped = ProbeForScreen("CR", "dRNA", ARM5_A, ARM3_A, seq, "GENE_CR")
    assert screen_self_circularisation([cramped]) == []


# ----------------------------------------------------------------------
# iLock chemistry
# ----------------------------------------------------------------------


#: Satisfies the invader invariant ``arm5[0] == arm3[-1]`` that
#: ``ext.mutation.probe.verify_iLock_probe`` enforces — both 'G' here.
ILOCK = ProbeForScreen(
    probe_id="ilock_A", chemistry="iLock",
    probe_arm5="GCAGCAGCAGCAGCAGCAGC",
    probe_arm3="ATGCTGCTGCTGCTGCTGCTG",     # 21 nt, ends in the shared 'G'
    sequence=("CGTTGCTGTGGCG"               # flap
              + "gcagcagcagcagcagcagc"      # arm5
              + BB
              + "atgctgctgctgctgctgctg").upper(),
    target="GIL",
)


def test_ilock_invader_invariant_holds_for_the_fixture():
    assert ILOCK.probe_arm5[0].upper() == ILOCK.probe_arm3[-1].upper()


def test_ilock_removes_arm5_first_base_not_arm3_first_base():
    """FEN1 cleaves the flap + ``arm5[0]`` block; ``arm3[-1]`` is the 3'-OH.

    v2.2 returned ``arm3[1:]`` with arm5 whole — right 3'-OH base by accident,
    5'-P one base early, so every iLock register was out of phase.
    """
    arm3_eff, arm5_eff, _ = build_v2_geom(ILOCK)
    assert arm3_eff == ILOCK.probe_arm3.upper(), "arm3 must stay whole"
    assert arm5_eff == ILOCK.probe_arm5[1:].upper(), "arm5 loses its first base"
    assert arm3_eff[-1] == ILOCK.probe_arm3[-1].upper()   # the 3'-OH


def test_ilock_junction_block_does_not_duplicate_the_snp_base():
    """The shared base appears once at the nick, on the arm3 side."""
    geom = _build_geometry(ILOCK, DEFAULT_VICINITY_N)
    arm3 = ILOCK.probe_arm3.upper()
    arm5 = ILOCK.probe_arm5.upper()
    assert geom.junction_block == arm3[-4:] + arm5[1:5]
    # The wrong v2.2 form would have carried arm3[-1] and arm5[0] — the same
    # nucleotide — on both sides of the nick.
    assert geom.junction_block != arm3[-4:] + arm5[:4]


def test_non_ilock_arms_are_untouched():
    arm3_eff, arm5_eff, _ = build_v2_geom(PROBE_A)
    assert arm3_eff == ARM3_A and arm5_eff == ARM5_A


# ----------------------------------------------------------------------
# Guards and reports
# ----------------------------------------------------------------------


def test_disabling_the_clamp_is_rejected():
    """n=0 would admit every position on every splint — fail loudly instead."""
    with pytest.raises(ValueError, match="vicinity_n_each_side"):
        screen_cross_ligation_v2([PROBE_A], vicinity_n_each_side=0)


def test_arm_missing_from_sequence_raises():
    bad = ProbeForScreen("bad", "dRNA", ARM5_A, ARM3_A, "ACGT" * 10, "G")
    with pytest.raises(ValueError, match="registry inconsistent"):
        screen_cross_ligation_v2([bad, PROBE_A])


def test_threshold_is_honoured():
    probes = [PROBE_A, _probe("B", F1 + SPLINT40 + F2)]
    assert _is_flagged(probes)
    _, dimers = screen_cross_ligation_v2(
        probes, tm_threshold_c=200.0, include_self_pairs=False,
    )
    assert dimers, "the register still exists"
    assert not any(d.flagged_overall for d in dimers), "but nothing clears 200 C"


def test_write_dimer_report_smoke(tmp_path: Path):
    _, dimers = screen_cross_ligation_v2(
        [PROBE_A, _probe("B", F1 + SPLINT40 + F2)], include_self_pairs=False,
    )
    out = tmp_path / "dimers.tsv"
    write_dimer_report(dimers, out)
    text = out.read_text(encoding="utf-8")
    header, *rows = text.splitlines()
    for column in ("seq_a_id", "overall_tm_c", "limiting_arm_tm_c",
                   "junction_run_nt", "nick_pos_on_b", "alignment"):
        assert column in header
    assert len(rows) == len(dimers)
    assert "\\n" in rows[0], "the alignment must be escaped onto one line"


def test_write_dimer_report_empty_writes_header_only(tmp_path: Path):
    out = tmp_path / "empty.tsv"
    write_dimer_report([], out)
    assert out.read_text(encoding="utf-8").strip().startswith("seq_a_id")


def test_write_self_circ_report_smoke(tmp_path: Path):
    out = tmp_path / "selfcirc.tsv"
    write_self_circ_report([], out)
    assert "probe_id" in out.read_text(encoding="utf-8")


def test_reaction_conditions_change_the_score():
    """The buffer is a parameter, not a constant baked into the screen."""
    probes = [PROBE_A, _probe("B", F1 + SPLINT40 + F2)]
    hot = ReactionConditions(monovalent_mM=500.0)
    cold = ReactionConditions(monovalent_mM=10.0)
    _, hot_dimers = screen_cross_ligation_v2(
        probes, reaction=hot, include_self_pairs=False,
    )
    _, cold_dimers = screen_cross_ligation_v2(
        probes, reaction=cold, include_self_pairs=False,
    )
    assert hot_dimers[0].overall_tm_c > cold_dimers[0].overall_tm_c
