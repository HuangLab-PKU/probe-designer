"""Tests for probe_designer.ext.tcr.probe."""
from __future__ import annotations

import pytest

from probe_designer.ext.tcr.probe import (
    TcrProbeDesigner,
    _reverse_complement,
    validate_cdr3_in_trb,
)


# ---------------------------------------------------------------------------
# validate_cdr3_in_trb
# ---------------------------------------------------------------------------


def test_validate_cdr3_returns_position_when_substring():
    trb = "ATGCATGCAAACDR3SEQATGCATGC"
    cdr3 = "AAACDR3SEQ"
    # _reverse_complement isn't used here, just need a real-looking pair
    pos = validate_cdr3_in_trb(trb, cdr3)
    assert pos == 8


def test_validate_cdr3_raises_when_not_in_trb():
    with pytest.raises(ValueError, match="CDR3.*not found"):
        validate_cdr3_in_trb("ATGCATGCAT", "GGGGGGGG")


def test_validate_cdr3_strips_n_and_uppercases():
    trb = "atgcatNNNcdr3seqNATGC"
    cdr3 = "cdr3seq"
    pos = validate_cdr3_in_trb(trb, cdr3)
    assert pos == 6  # after uppercasing + N-strip


# ---------------------------------------------------------------------------
# Backwards-compat shim still works
# ---------------------------------------------------------------------------


def test_legacy_import_path_still_works():
    from probe_designer.ext.tcr_probe import TcrProbeDesigner as Legacy
    assert Legacy is TcrProbeDesigner


# ---------------------------------------------------------------------------
# Designer construction
# ---------------------------------------------------------------------------


def test_designer_default_arm_len_is_half_bds():
    d = TcrProbeDesigner(bds_len=40)
    assert d.arm_len == 20


def test_designer_rejects_arm_len_not_half_bds():
    with pytest.raises(ValueError, match="arm_len.*bds_len"):
        TcrProbeDesigner(bds_len=40, arm_len=15)


# ---------------------------------------------------------------------------
# find_cdr3_binding_sites — dual chemistry output
# ---------------------------------------------------------------------------


@pytest.fixture
def trb_cdr3_pair():
    """A small synthetic TRB with a clearly-marked CDR3 region."""
    # 100 nt TRB, CDR3 starts at position 30, length 30
    upstream = "AAAAATTTTAAAACGCGCATGCATGCATGCA"  # 30+1 nt
    cdr3 = "GGGCATCGATCGATCGATCGATCGATCGATCG"  # 32 nt
    downstream = "TGCATGCAATGCATGAACGTACGTACGTACGTACGTACGT"
    return {
        "trb": upstream + cdr3 + downstream,
        "cdr3": cdr3,
        "cdr3_start_expected": len(upstream),
    }


def test_dual_chemistry_arms_present_for_each_site(trb_cdr3_pair):
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=0.0, max_arm_tm=200.0)
    sites = d.find_cdr3_binding_sites(
        trb_cdr3_pair["trb"], trb_cdr3_pair["cdr3"], "TestClone1",
    )
    assert sites, "expected at least one BDS candidate"
    s = sites[0]
    # Both chemistry arm pairs present
    for key in ("arm_5prime_dRNA", "arm_3prime_dRNA",
                "arm_5prime_cDNA", "arm_3prime_cDNA"):
        assert key in s, f"missing {key}"
        assert len(s[key]) == 20
    # dRNA arm = RC of target slice; cDNA arm = literal target slice
    target = s["target_sequence"]
    assert s["arm_3prime_cDNA"].upper() == target[:20]
    assert s["arm_5prime_cDNA"].upper() == target[20:]
    assert s["arm_3prime_dRNA"].upper() == _reverse_complement(target)[:20]
    assert s["arm_5prime_dRNA"].upper() == _reverse_complement(target)[20:]


def test_dual_chemistry_tms_present(trb_cdr3_pair):
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=0.0, max_arm_tm=200.0)
    sites = d.find_cdr3_binding_sites(
        trb_cdr3_pair["trb"], trb_cdr3_pair["cdr3"], "TestClone1",
    )
    s = sites[0]
    for key in ("tm_dRNA", "tm_5prime_dRNA", "tm_3prime_dRNA",
                "tm_cDNA", "tm_5prime_cDNA", "tm_3prime_cDNA"):
        assert key in s
        assert isinstance(s[key], float)
    # cDNA Tm should be slightly higher than dRNA Tm for a balanced sequence
    assert s["tm_cDNA"] > 0
    assert s["tm_dRNA"] > 0


def test_ligation_point_is_inside_cdr3(trb_cdr3_pair):
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=0.0, max_arm_tm=200.0)
    sites = d.find_cdr3_binding_sites(
        trb_cdr3_pair["trb"], trb_cdr3_pair["cdr3"], "TestClone1",
    )
    cdr3_start = trb_cdr3_pair["cdr3_start_expected"]
    cdr3_end = cdr3_start + len(trb_cdr3_pair["cdr3"])
    for s in sites:
        lp = s["ligation_point"]
        assert cdr3_start <= lp <= cdr3_end, \
            f"ligation_point {lp} outside CDR3 [{cdr3_start},{cdr3_end}]"


def test_returns_empty_when_cdr3_not_in_trb():
    d = TcrProbeDesigner(bds_len=40, arm_len=20)
    sites = d.find_cdr3_binding_sites("ATGCATGCATGC", "GGGGGGGG", "X")
    assert sites == []


# ---------------------------------------------------------------------------
# Filter + select
# ---------------------------------------------------------------------------


def test_filter_by_chemistry_tm_drops_outside_window(trb_cdr3_pair):
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=200.0, max_arm_tm=300.0)
    sites = d.find_cdr3_binding_sites(
        trb_cdr3_pair["trb"], trb_cdr3_pair["cdr3"], "X",
    )
    # impossible window ⇒ filter drops everything
    assert d.filter_by_chemistry_tm(sites, chemistry="dRNA") == []


def test_select_non_overlapping_respects_min_gap():
    d = TcrProbeDesigner(bds_len=40, arm_len=20)
    candidates = [
        {"st": 0, "score": 1.0},
        {"st": 10, "score": 5.0},   # within bds_len of #1 ⇒ should be skipped
        {"st": 50, "score": 3.0},   # 50 >= 40 from #1 ⇒ keepable
        {"st": 100, "score": 2.0},  # 100 - 50 = 50 ⇒ keepable
    ]
    selected = d.select_non_overlapping_sites(candidates, max_sites=3)
    starts = sorted(s["st"] for s in selected)
    # Greedy by score: picks 10 first (best score). 0 violates gap. 50 OK. 100 OK.
    # But 10 and 50 differ by 40 = bds_len ≥ min_gap ⇒ both kept.
    # Then 100 - 50 = 50 ≥ 40 ⇒ kept.
    assert starts == [10, 50, 100]


# ---------------------------------------------------------------------------
# Tm input convention — independent per chemistry
# ---------------------------------------------------------------------------

def test_drna_tm_computed_on_rna_target_strand(trb_cdr3_pair):
    """dRNA Tm must feed R_DNA_NN1 the RNA (target) strand — the reverse
    complement of each DNA probe arm — per Biopython's convention ("seq must be
    the RNA sequence"), evaluated at the reaction buffer + formamide.

    Before 2026-07-17 the code passed the DNA probe strand (the target_rc
    halves) directly, which is the wrong convention for the asymmetric
    R_DNA_NN1 table (see experiments/20260715_tm_deltag_methods_audit/).
    """
    from Bio.SeqUtils import MeltingTemp as mt
    from probe_designer.chemistry import ReactionConditions, dna_revcomp_to_rna

    rc = ReactionConditions()
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=0.0, max_arm_tm=200.0)
    sites = d.find_cdr3_binding_sites(
        trb_cdr3_pair["trb"], trb_cdr3_pair["cdr3"], "TestClone1",
    )
    s = sites[0]
    target = s["target_sequence"]
    target_rc = _reverse_complement(target)

    def expected(dna_arm):
        tm = mt.Tm_NN(dna_revcomp_to_rna(dna_arm), nn_table=mt.R_DNA_NN1,
                      **rc.tm_nn_kwargs())
        return round(rc.apply_formamide(tm), 2)

    # Probe arms are the RC of the target halves (target_rc halves).
    assert s["tm_5prime_dRNA"] == expected(target_rc[20:]), (
        f"tm_5prime_dRNA should be Tm of the RNA-sense strand; got "
        f"{s['tm_5prime_dRNA']}, expected {expected(target_rc[20:])}"
    )
    assert s["tm_3prime_dRNA"] == expected(target_rc[:20])
    assert s["tm_dRNA"] == expected(target_rc)


def test_cdna_tm_uses_target_half_directly(trb_cdr3_pair):
    """cDNA arm = target half (sense DNA, since cDNA = RC(mRNA) and the padlock
    arms read sense again). DNA_NN4 is symmetric; computed at the reaction
    buffer + formamide."""
    from Bio.SeqUtils import MeltingTemp as mt
    from probe_designer.chemistry import ReactionConditions

    rc = ReactionConditions()
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=0.0, max_arm_tm=200.0)
    sites = d.find_cdr3_binding_sites(
        trb_cdr3_pair["trb"], trb_cdr3_pair["cdr3"], "TestClone1",
    )
    s = sites[0]
    target = s["target_sequence"]

    def expected(seq):
        tm = mt.Tm_NN(seq, nn_table=mt.DNA_NN4, **rc.tm_nn_kwargs())
        return round(rc.apply_formamide(tm), 2)

    assert s["tm_5prime_cDNA"] == expected(target[20:])
    assert s["tm_3prime_cDNA"] == expected(target[:20])
    assert s["tm_cDNA"] == expected(target)


def test_per_chemistry_tm_diff_and_score_present(trb_cdr3_pair):
    """Each site carries chemistry-specific tm_diff and score columns so the
    pipeline can filter/select per chemistry."""
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=0.0, max_arm_tm=200.0)
    sites = d.find_cdr3_binding_sites(
        trb_cdr3_pair["trb"], trb_cdr3_pair["cdr3"], "TestClone1",
    )
    s = sites[0]
    for key in ("tm_diff_dRNA", "tm_diff_cDNA"):
        assert key in s and isinstance(s[key], float)


# ---------------------------------------------------------------------------
# Per-chemistry filter helper
# ---------------------------------------------------------------------------

def test_filter_by_chemistry_tm_drna_gate(trb_cdr3_pair):
    """filter_by_chemistry_tm with chemistry='dRNA' uses the dRNA arm Tms."""
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=200.0, max_arm_tm=300.0)
    sites = d.find_cdr3_binding_sites(
        trb_cdr3_pair["trb"], trb_cdr3_pair["cdr3"], "X",
    )
    # impossible window ⇒ drops everything
    assert d.filter_by_chemistry_tm(sites, chemistry="dRNA") == []


def test_filter_by_chemistry_tm_cdna_gate_independent():
    """Synthetic sites with very different dRNA vs cDNA Tms — filter passes
    one chemistry but not the other."""
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=50.0, max_arm_tm=60.0)
    sites = [
        {"tm_5prime_dRNA": 55.0, "tm_3prime_dRNA": 55.0,
         "tm_5prime_cDNA": 80.0, "tm_3prime_cDNA": 80.0},  # passes dRNA, fails cDNA
        {"tm_5prime_dRNA": 80.0, "tm_3prime_dRNA": 80.0,
         "tm_5prime_cDNA": 55.0, "tm_3prime_cDNA": 55.0},  # passes cDNA, fails dRNA
    ]
    passed_drna = d.filter_by_chemistry_tm(sites, chemistry="dRNA")
    passed_cdna = d.filter_by_chemistry_tm(sites, chemistry="cDNA")
    assert len(passed_drna) == 1 and passed_drna[0]["tm_5prime_dRNA"] == 55.0
    assert len(passed_cdna) == 1 and passed_cdna[0]["tm_5prime_cDNA"] == 55.0
