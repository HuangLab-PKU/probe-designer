"""Tests for ``probe_designer.ext.nupack.check``.

Two things are pinned here that a future reader will otherwise get wrong:

1. ``test_padlock_wrap_is_a_pseudoknot`` — the arms CANNOT be modelled as one
   tethered strand, because the wrapped complex is a pseudoknot and NUPACK's
   ensemble is nested. This is why the arms are split, and the test exists so
   that "obvious improvement" is not attempted a third time.
2. The tether is paid back as effective concentration, not ignored.
"""
from __future__ import annotations

from pathlib import Path

import pytest

# Skip the whole module if NUPACK 4 is not installed
pytest.importorskip("nupack")


from probe_designer.qc.cross_ligation import ProbeForScreen  # noqa: E402
from probe_designer.ext.nupack.check import (  # noqa: E402
    TernaryHit,
    screen_ternary,
    tether_effective_concentration,
    write_ternary_report,
    _check_mfe_nick_geometry,
    _check_mfe_vicinity_contiguous,
    _loop_length_nt,
    _parse_dotbracket_pairs,
)


#: Non-self-complementary backbone. It must NOT be an ``ACGT`` repeat: ACGT is
#: its own reverse complement, so two probes carrying such a backbone form a
#: full-length backbone duplex that dominates the MFE and hides the arms.
BB = "TCCCTACACGACGCTCTTCCGATCTACACTCTTTCCCTACACGACGCTCTTCCGA"

#: B's arm5 is RC(arm3_A + arm5_A), so B templates A's whole 40 nt footprint.
XLIG_A = ProbeForScreen(
    probe_id="xlig_A", chemistry="dRNA",
    probe_arm5="CCCAATTGCGCAATATCATG",
    probe_arm3="GGGGGAAAATTTCCCAAGGG",
    sequence=("CCCAATTGCGCAATATCATG" + BB + "GGGGGAAAATTTCCCAAGGG").upper(),
    target="GX",
)
XLIG_B = ProbeForScreen(
    probe_id="xlig_B", chemistry="dRNA",
    probe_arm5="CATGATATTGCGCAATTGGGCCCTTGGGAAATTTTCCCCC",
    probe_arm3="AAGCTTAACTGGCCATAAGT",
    sequence=("CATGATATTGCGCAATTGGGCCCTTGGGAAATTTTCCCCC"
              + BB + "AAGCTTAACTGGCCATAAGT").upper(),
    target="GY",
)

PAIR = [("xlig_A", "xlig_B")]


# ----------------------------------------------------------------------
# The pseudoknot constraint — the reason the arms are split
# ----------------------------------------------------------------------


def test_padlock_wrap_is_a_pseudoknot_so_the_arms_must_stay_split():
    """One tethered ligator strand + splint cannot express the wrap.

    Minimal case: two 10 nt arms perfectly complementary to ADJACENT blocks of
    a 20 nt splint. If the nested ensemble could hold the padlock geometry, the
    MFE would be the full 20 bp wrap with both ligator termini paired. It is
    not — one arm binds and the other is left dangling, because the two duplex
    segments cross.
    """
    import nupack

    def rc(s: str) -> str:
        return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    arm5, arm3 = "ACGTGCAGTC", "TTGCCAGATC"
    ligator = arm5 + "TTTTTTTTTTTT" + arm3
    splint = rc(arm5) + rc(arm3)

    model = nupack.Model(material="dna", celsius=37, sodium=0.075,
                         magnesium=0.01, ensemble="stacking")
    cplx = nupack.Complex([
        nupack.Strand(ligator, name="lig"), nupack.Strand(splint, name="spl"),
    ])
    result = nupack.complex_analysis(
        complexes=[cplx], model=model, compute=["mfe"],
    )
    pairs = _parse_dotbracket_pairs(str(result[cplx].mfe[0].structure))

    both_termini_paired = (0 in pairs) and (len(ligator) - 1 in pairs)
    assert not both_termini_paired, (
        "NUPACK expressed a padlock wrap — if this ever passes, the nested-model "
        "limitation is gone and check.py can go back to a 2-strand tethered tube"
    )


# ----------------------------------------------------------------------
# The tether, as effective concentration
# ----------------------------------------------------------------------


def test_tether_concentration_is_far_above_bulk():
    """A 55 nt backbone holds the second arm at ~1e-4 M, not ~1e-7 M."""
    c_eff = tether_effective_concentration(55)
    assert 1e-4 < c_eff < 1e-2
    assert c_eff / 1e-7 > 1000, "the tether must be worth orders of magnitude"


def test_tether_concentration_falls_with_loop_length():
    assert (
        tether_effective_concentration(20)
        > tether_effective_concentration(55)
        > tether_effective_concentration(200)
    )


def test_tether_concentration_rejects_a_zero_loop():
    with pytest.raises(ValueError, match="loop"):
        tether_effective_concentration(0)


def test_loop_length_is_the_backbone_between_the_arms():
    assert _loop_length_nt(XLIG_A) == len(BB)


# ----------------------------------------------------------------------
# Dot-bracket parser
# ----------------------------------------------------------------------


def test_parse_dotbracket_simple_pair():
    pairs = _parse_dotbracket_pairs("((..))")
    assert pairs[0] == 5 and pairs[5] == 0
    assert pairs[1] == 4 and pairs[4] == 1
    assert 2 not in pairs and 3 not in pairs


def test_parse_dotbracket_multi_strand_separator():
    """'+' separators must NOT advance the position counter."""
    assert _parse_dotbracket_pairs("..+..+..") == {}
    pairs = _parse_dotbracket_pairs("((+))")
    assert pairs[0] == 3 and pairs[1] == 2


# ----------------------------------------------------------------------
# MFE geometry helpers — concat order (arm3 | splint | arm5)
# ----------------------------------------------------------------------


ARM3_LEN, SPLINT_LEN, ARM5_LEN = 20, 100, 20
ARM5_BASE = ARM3_LEN + SPLINT_LEN


def _sp(pos: int) -> int:
    return ARM3_LEN + pos


def test_nick_geometry_true_when_termini_are_adjacent():
    pairs = {ARM3_LEN - 1: _sp(51), ARM5_BASE: _sp(50)}
    ok, b_3oh, b_5p = _check_mfe_nick_geometry(pairs, ARM3_LEN, SPLINT_LEN)
    assert ok and (b_3oh, b_5p) == (51, 50)


def test_nick_geometry_false_when_not_adjacent():
    pairs = {ARM3_LEN - 1: _sp(70), ARM5_BASE: _sp(50)}
    ok, _, _ = _check_mfe_nick_geometry(pairs, ARM3_LEN, SPLINT_LEN)
    assert not ok


def test_nick_geometry_false_when_a_terminus_is_unpaired():
    ok, b_3oh, b_5p = _check_mfe_nick_geometry(
        {ARM3_LEN - 1: _sp(51)}, ARM3_LEN, SPLINT_LEN,
    )
    assert not ok and b_3oh is None and b_5p is None


def test_nick_geometry_false_when_a_terminus_pairs_outside_the_splint():
    """Arm-to-arm pairing is not a splinted nick."""
    pairs = {ARM3_LEN - 1: ARM5_BASE + 5, ARM5_BASE: _sp(50)}
    ok, _, _ = _check_mfe_nick_geometry(pairs, ARM3_LEN, SPLINT_LEN)
    assert not ok


def _clamped(n: int, b_5p: int = 50) -> dict:
    pairs = {}
    for i in range(n + 1):
        pairs[ARM3_LEN - 1 - i] = _sp(b_5p + 1 + i)      # arm3 side
        pairs[ARM5_BASE + i] = _sp(b_5p - i)             # arm5 side
    return pairs


def test_vicinity_accepts_a_clean_clamp():
    assert _check_mfe_vicinity_contiguous(
        _clamped(3), ARM3_LEN, SPLINT_LEN, ARM5_LEN, 3,
    )


def test_vicinity_rejects_a_bulge():
    pairs = _clamped(3)
    pairs[ARM5_BASE + 2] = _sp(10)
    assert not _check_mfe_vicinity_contiguous(
        pairs, ARM3_LEN, SPLINT_LEN, ARM5_LEN, 3,
    )


def test_vicinity_rejects_an_unpaired_clamp_base():
    pairs = _clamped(3)
    del pairs[ARM3_LEN - 2]
    assert not _check_mfe_vicinity_contiguous(
        pairs, ARM3_LEN, SPLINT_LEN, ARM5_LEN, 3,
    )


def test_vicinity_disabled_when_n_is_zero():
    assert not _check_mfe_vicinity_contiguous(
        _clamped(3), ARM3_LEN, SPLINT_LEN, ARM5_LEN, 0,
    )


# ----------------------------------------------------------------------
# End-to-end against NUPACK
# ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def hit() -> TernaryHit:
    hits = screen_ternary([XLIG_A, XLIG_B], prefilter_pairs=PAIR)
    assert len(hits) == 1
    return hits[0]


def test_ternary_is_enumerated(hit):
    assert hit.seq_a_id == "xlig_A" and hit.seq_b_id == "xlig_B"
    assert hit.arm3_len == 20 and hit.arm5_len == 20
    assert hit.splint_len == len(XLIG_B.sequence)
    assert hit.mfe_dotbracket.count("+") == 2, "three strands"


def test_ternary_is_bound_and_favourable(hit):
    assert hit.ternary_dg_kcal < 0.0
    assert hit.ternary_concentration_m > 0.0
    assert 0.0 <= hit.ternary_fraction_of_b <= 1.0


def test_perfect_splint_gives_nick_adjacent_geometry(hit):
    assert hit.mfe_nick_adjacent, hit.mfe_dotbracket
    assert hit.b_3oh_pos - hit.b_5p_pos == 1
    assert hit.mfe_vicinity_contiguous


def test_tether_is_recorded_in_the_hit(hit):
    assert hit.tether_loop_nt == len(BB)
    assert hit.arm_conc_m == tether_effective_concentration(len(BB))


def test_tether_raises_the_productive_fraction():
    """The v2.3 behaviour (arms at bulk) is available, and gives less complex."""
    tethered = screen_ternary([XLIG_A, XLIG_B], prefilter_pairs=PAIR)[0]
    untethered = screen_ternary(
        [XLIG_A, XLIG_B], prefilter_pairs=PAIR, arm_conc_m=1.0e-7,
    )[0]
    assert tethered.ternary_fraction_of_b > untethered.ternary_fraction_of_b


def test_missing_probe_id_raises_rather_than_skipping():
    """A silently dropped pair is indistinguishable from a clean one."""
    with pytest.raises(KeyError, match="absent"):
        screen_ternary([XLIG_A], prefilter_pairs=[("xlig_A", "nope")])


def test_write_report_smoke(tmp_path: Path, hit):
    out = tmp_path / "ternary.tsv"
    write_ternary_report([hit], out)
    text = out.read_text(encoding="utf-8")
    for column in (
        "seq_a_id", "ternary_dg_kcal", "ternary_concentration_m",
        "ternary_fraction_of_b", "mfe_nick_adjacent", "mfe_vicinity_contiguous",
        "flagged", "tether_loop_nt", "arm_conc_m", "ensemble_used",
        "mfe_dotbracket",
    ):
        assert column in text, f"missing column {column}"
    assert len(text.splitlines()) == 2


def test_write_report_empty_writes_header_only(tmp_path: Path):
    out = tmp_path / "empty.tsv"
    write_ternary_report([], out)
    assert out.read_text(encoding="utf-8").strip().startswith("seq_a_id")


def test_coaxial_stacking_makes_the_nicked_ternary_more_stable():
    """``ensemble='stacking'`` is the reason to reach for NUPACK at all."""
    stacked = screen_ternary(
        [XLIG_A, XLIG_B], prefilter_pairs=PAIR, ensemble="stacking",
    )[0]
    plain = screen_ternary(
        [XLIG_A, XLIG_B], prefilter_pairs=PAIR, ensemble="nostacking",
    )[0]
    assert stacked.ternary_dg_kcal < plain.ternary_dg_kcal
