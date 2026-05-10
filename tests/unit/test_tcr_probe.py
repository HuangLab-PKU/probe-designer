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


def test_filter_by_mrna_tm_drops_outside_window(trb_cdr3_pair):
    d = TcrProbeDesigner(bds_len=40, arm_len=20, min_arm_tm=200.0, max_arm_tm=300.0)
    sites = d.find_cdr3_binding_sites(
        trb_cdr3_pair["trb"], trb_cdr3_pair["cdr3"], "X",
    )
    # impossible window ⇒ filter drops everything
    assert d.filter_by_mrna_tm(sites) == []


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
