"""Tests for the Magoulopoulou/Nilsson padlock site-quality terms.

Rules come from Magoulopoulou et al. 2026, eBioMedicine 127:106264 (Fig. 1B +
Methods): a candidate target k-mer should carry *less than four repetitive
bases*, *G or C at the ligation position*, and *40-60% GC*.

The paper states them on a 30-mer with the ligation point at the 16th base.
Here they are expressed on our 40-mer BDS (20+20 arms), anchored on probe
geometry rather than on a fixed index so the same definition serves both
chemistries.
"""
from __future__ import annotations

import pytest

from probe_designer.ext.tcr.nilsson import (
    NilssonWeights,
    gc_fraction,
    junction_base,
    max_homopolymer_run,
    nilsson_terms,
)


# ---------------------------------------------------------------------------
# max_homopolymer_run
# ---------------------------------------------------------------------------


def test_homopolymer_run_counts_longest_stretch():
    assert max_homopolymer_run("ACGTACGT") == 1
    assert max_homopolymer_run("AACCGGTT") == 2
    assert max_homopolymer_run("ACGGGGTA") == 4
    assert max_homopolymer_run("GGGGGACGT") == 5


def test_homopolymer_run_is_case_insensitive():
    assert max_homopolymer_run("acgggGta") == 4


def test_homopolymer_run_of_empty_sequence_is_zero():
    assert max_homopolymer_run("") == 0


# ---------------------------------------------------------------------------
# gc_fraction
# ---------------------------------------------------------------------------


def test_gc_fraction_is_a_fraction_not_a_percent():
    assert gc_fraction("GCGCATAT") == pytest.approx(0.5)
    assert gc_fraction("AAAA") == pytest.approx(0.0)
    assert gc_fraction("gcgc") == pytest.approx(1.0)


def test_gc_fraction_of_empty_sequence_is_zero():
    assert gc_fraction("") == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# junction_base — the target base paired with the probe's 3'-OH terminus
# ---------------------------------------------------------------------------


def test_junction_base_drna_is_first_base_of_three_prime_half():
    # dRNA arms are RC of the target halves; the probe's 3' end pairs with
    # target[arm_len] — the paper's "16th position" of a 15+15 30-mer.
    target = "A" * 20 + "G" + "C" * 19
    assert junction_base(target, arm_len=20, chemistry="dRNA") == "G"


def test_junction_base_cdna_is_last_base_of_five_prime_half():
    # cDNA arms are the literal target halves, so the probe's 3' end pairs
    # with the complement of target[arm_len - 1]; G/C-ness is preserved under
    # complementation, so the reported base is target[arm_len - 1] itself.
    target = "A" * 19 + "C" + "T" * 20
    assert junction_base(target, arm_len=20, chemistry="cDNA") == "C"


def test_junction_base_rejects_unknown_chemistry():
    with pytest.raises(KeyError):
        junction_base("ACGT" * 10, arm_len=20, chemistry="iLock")


def test_junction_base_uppercases():
    target = "a" * 20 + "g" + "c" * 19
    assert junction_base(target, arm_len=20, chemistry="dRNA") == "G"


# ---------------------------------------------------------------------------
# nilsson_terms — the scoring contribution
# ---------------------------------------------------------------------------


def _perfect_target() -> str:
    """40-mer: no run >= 4, G at the dRNA junction, GC = 50%."""
    target = "ACGTTGCAACGTTGCAACGT" + "GATGCCATGATGCCATGATG"
    assert len(target) == 40
    return target


def test_perfect_site_passes_all_three_rules():
    t = _perfect_target()
    terms = nilsson_terms(t, arm_len=20, chemistry="dRNA")
    assert terms["homopolymer_pass"] is True
    assert terms["junction_gc_pass"] is True
    assert terms["gc_pass"] is True
    assert terms["nilsson_pass"] is True


def test_perfect_site_scores_higher_than_homopolymer_site():
    good = nilsson_terms(_perfect_target(), arm_len=20, chemistry="dRNA")
    bad_seq = "ACGTTGCAACGTTGCAAAAA" + "GATGCCATGATGCCATGCAT"
    bad = nilsson_terms(bad_seq, arm_len=20, chemistry="dRNA")
    assert bad["homopolymer_pass"] is False
    assert bad["nilsson_score"] < good["nilsson_score"]


def test_at_junction_loses_the_junction_bonus():
    t = _perfect_target()
    at_junction = t[:20] + "T" + t[21:]
    terms = nilsson_terms(at_junction, arm_len=20, chemistry="dRNA")
    assert terms["junction_base"] == "T"
    assert terms["junction_gc_pass"] is False
    assert terms["nilsson_pass"] is False


def test_gc_outside_window_is_penalised_proportionally():
    near = nilsson_terms("GC" * 4 + "AT" * 16, arm_len=20, chemistry="dRNA")
    far = nilsson_terms("GC" * 2 + "AT" * 18, arm_len=20, chemistry="dRNA")
    assert near["gc_pass"] is False and far["gc_pass"] is False
    assert far["nilsson_score"] < near["nilsson_score"]


def test_gc_window_boundaries_are_inclusive():
    lo = nilsson_terms("GC" * 8 + "AT" * 12, arm_len=20, chemistry="dRNA")
    hi = nilsson_terms("GC" * 12 + "AT" * 8, arm_len=20, chemistry="dRNA")
    assert lo["gc_content"] == pytest.approx(0.40)
    assert hi["gc_content"] == pytest.approx(0.60)
    assert lo["gc_pass"] is True and hi["gc_pass"] is True


def test_zero_weights_make_the_contribution_vanish():
    t = "AAAA" + "T" * 16 + "T" * 20  # violates all three rules
    terms = nilsson_terms(t, arm_len=20, chemistry="dRNA",
                          weights=NilssonWeights(homopolymer=0.0,
                                                 junction_gc=0.0, gc=0.0))
    assert terms["nilsson_score"] == pytest.approx(0.0)
    # Flags still report the truth even when the weights silence the score.
    assert terms["nilsson_pass"] is False


def test_score_is_never_positive_beyond_the_junction_bonus():
    """The only positive term is the junction bonus, so a site can never buy
    its way past a real Tm/structure problem."""
    t = _perfect_target()
    w = NilssonWeights()
    terms = nilsson_terms(t, arm_len=20, chemistry="dRNA", weights=w)
    assert terms["nilsson_score"] == pytest.approx(w.junction_gc)


# ---------------------------------------------------------------------------
# Config + pipeline wiring
# ---------------------------------------------------------------------------


def test_config_defaults_to_no_nilsson_contribution(tmp_path):
    """Off by default: every pre-2026-08 TCR run must stay reproducible."""
    from probe_designer.ext.tcr.config import TcrConfig

    cfg = TcrConfig(input_xlsx=tmp_path / "i.xlsx", output_dir=tmp_path,
                    backbone_file=tmp_path / "b.xlsx", start_no=1)
    assert cfg.nilsson_weights is None


def test_config_accepts_explicit_weights(tmp_path):
    from probe_designer.ext.tcr.config import TcrConfig

    cfg = TcrConfig(input_xlsx=tmp_path / "i.xlsx", output_dir=tmp_path,
                    backbone_file=tmp_path / "b.xlsx", start_no=1,
                    nilsson_weights=(1.0, 2.0, 3.0))
    assert cfg.nilsson_weights == NilssonWeights(homopolymer=1.0,
                                                 junction_gc=2.0, gc=3.0)


def test_config_rejects_negative_weights(tmp_path):
    from probe_designer.ext.tcr.config import TcrConfig

    with pytest.raises(ValueError, match="nilsson_weights"):
        TcrConfig(input_xlsx=tmp_path / "i.xlsx", output_dir=tmp_path,
                  backbone_file=tmp_path / "b.xlsx", start_no=1,
                  nilsson_weights=(-1.0, 2.0, 3.0))


def _site_scores(tmp_path, weights):
    """Run phase 1 on a two-clone frame and return {clone: [(pos, score)]}."""
    import pandas as pd

    from probe_designer.ext.tcr.config import TcrConfig
    from probe_designer.ext.tcr.pipeline import _phase1_find_select

    # CDR3 long enough to admit several ligation points.
    cdr3 = "TGTGCCAGCAGTTTATGGGACACAGATACGCAGTATTTT"
    trb = ("ACCTGGAGCCCCCAGAACTGGCAGACACCTGCCTGATGCTGCCATGGGCCCCCAGCTCC"
           + cdr3
           + "GGCCCAGGCACCCGGCTGACAGTGCTCGAGGACCTGAAAAACGTGTTCCCACCC")
    clones = pd.DataFrame([{"consensus_id": "ct1_consensus1", "chain_nt": trb,
                            "cdr3_nt": cdr3, "Clone_name": "TEST_ct1"}])
    cfg = TcrConfig(input_xlsx=tmp_path / "i.xlsx", output_dir=tmp_path,
                    backbone_file=tmp_path / "b.xlsx", start_no=1,
                    chemistries=["cDNA"], nilsson_weights=weights,
                    plot_tm_landscape=False, skip_blast=True)
    run_dir = tmp_path / "run"
    run_dir.mkdir(exist_ok=True)
    _all, selected, _skipped = _phase1_find_select(clones, cfg, run_dir)
    return selected["cDNA"]["ct1_consensus1"]


def test_pipeline_attaches_nilsson_flags_when_enabled(tmp_path):
    sites = _site_scores(tmp_path, (2.0, 2.0, 2.0))
    assert sites, "expected at least one selected site"
    for s in sites:
        assert "nilsson_pass" in s
        assert "junction_base" in s
        assert "homopolymer_run" in s


def test_pipeline_omits_nilsson_flags_when_disabled(tmp_path):
    sites = _site_scores(tmp_path, None)
    assert sites
    for s in sites:
        assert "nilsson_pass" not in s


def test_enabling_nilsson_changes_the_selection(tmp_path):
    """The whole point of the contribution: it re-ranks. Either a different
    position wins, or the same one wins with a different score."""
    off = {s["pos"]: s["score"] for s in _site_scores(tmp_path, None)}
    on = {s["pos"]: s["score"] for s in _site_scores(tmp_path, (2.0, 2.0, 2.0))}
    assert off != on


def test_selected_site_score_includes_its_nilsson_term(tmp_path):
    for s in _site_scores(tmp_path, (2.0, 2.0, 2.0)):
        base = nilsson_terms(s["target_sequence"], 20, "cDNA",
                             weights=NilssonWeights(2.0, 2.0, 2.0))
        assert s["nilsson_score"] == pytest.approx(base["nilsson_score"])
        assert s["junction_base"] == base["junction_base"]


def test_zero_weights_leave_the_score_untouched(tmp_path):
    off = {s["pos"]: s["score"] for s in _site_scores(tmp_path, None)}
    zeroed = {s["pos"]: s["score"] for s in _site_scores(tmp_path, (0.0, 0.0, 0.0))}
    assert off == zeroed


# ---------------------------------------------------------------------------
# Ligation-margin constraint — clone specificity, not a Nilsson rule
# ---------------------------------------------------------------------------


def _margin_phase1(tmp_path, margin):
    """Phase 1 over a short-CDR3 clone plus a long-CDR3 companion.

    The companion keeps the run non-empty so the short clone's fate is
    observable as a per-clone skip rather than a whole-run error.
    """
    import pandas as pd

    from probe_designer.ext.tcr.config import TcrConfig
    from probe_designer.ext.tcr.pipeline import _phase1_find_select

    cdr3 = "TGCAGTGCTCGTTCGGGGGGTGCCATGACTGAAGCTTTCTTT"
    trb = ("GACAGCAGCTTCTACATCTGCAGTGCTCGCACGTAGCATCGAAGCTTAGCTAAGCTT"
           + cdr3
           + "GGCAGTGGAACCCAGCTCTCTGTCTTGGAGGACCTGAACAAGGTGTTCCCACCCGAG")
    long_cdr3 = cdr3 + "ACGTTGCAACGTTGCAACGTGATGCCATGATGCCATGATGCA"
    long_trb = ("GACAGCAGCTTCTACATCTGCAGTGCTCGCACGTAGCATCGAAGCTTAGCTAAGCTT"
                + long_cdr3
                + "GGCAGTGGAACCCAGCTCTCTGTCTTGGAGGACCTGAACAAGGTGTTCCCACCCGAG")
    clones = pd.DataFrame([
        {"consensus_id": "ct74_consensus1", "chain_nt": trb,
         "cdr3_nt": cdr3, "Clone_name": "TEST_ct74"},
        {"consensus_id": "ctlong_consensus1", "chain_nt": long_trb,
         "cdr3_nt": long_cdr3, "Clone_name": "TEST_ctlong"},
    ])
    cfg = TcrConfig(input_xlsx=tmp_path / "i.xlsx", output_dir=tmp_path,
                    backbone_file=tmp_path / "b.xlsx", start_no=1,
                    chemistries=["cDNA"], sites_per_clone=3,
                    min_ligation_margin=margin,
                    plot_tm_landscape=False, skip_blast=True)
    run_dir = tmp_path / f"run{margin}"
    run_dir.mkdir(exist_ok=True)
    return _phase1_find_select(clones, cfg, run_dir)


def _margin_sites(tmp_path, margin):
    _all, selected, _skipped = _margin_phase1(tmp_path, margin)
    return selected["cDNA"].get("ct74_consensus1", [])


def test_min_ligation_margin_defaults_to_zero(tmp_path):
    from probe_designer.ext.tcr.config import TcrConfig

    cfg = TcrConfig(input_xlsx=tmp_path / "i.xlsx", output_dir=tmp_path,
                    backbone_file=tmp_path / "b.xlsx", start_no=1)
    assert cfg.min_ligation_margin == 0


def test_zero_margin_admits_nicks_at_the_cdr3_edge(tmp_path):
    """Unconstrained, the scan offers nicks flush with the CDR3 boundary —
    the defect the margin exists to remove."""
    all_sites, _selected, _skipped = _margin_phase1(tmp_path, 0)
    margins = [min(s["ligation_point"] - s["cdr3_start"],
                   s["cdr3_end"] - s["ligation_point"])
               for s in all_sites["ct74_consensus1"]]
    assert min(margins) == 0


def test_margin_three_keeps_the_splintr_clamp_inside_cdr3(tmp_path):
    sites = _margin_sites(tmp_path, 3)
    assert sites, "constraining the nick must not wipe out the clone"
    for s in sites:
        assert s["ligation_point"] - s["cdr3_start"] >= 3
        assert s["cdr3_end"] - s["ligation_point"] >= 3


def test_margin_too_wide_for_the_cdr3_skips_the_clone(tmp_path):
    _all, selected, skipped = _margin_phase1(tmp_path, 40)
    assert "ct74_consensus1" in skipped
    assert "ct74_consensus1" not in selected["cDNA"]


def test_config_rejects_negative_margin(tmp_path):
    from probe_designer.ext.tcr.config import TcrConfig

    with pytest.raises(ValueError, match="min_ligation_margin"):
        TcrConfig(input_xlsx=tmp_path / "i.xlsx", output_dir=tmp_path,
                  backbone_file=tmp_path / "b.xlsx", start_no=1,
                  min_ligation_margin=-1)


# ---------------------------------------------------------------------------
# Repertoire uniqueness — a site must not occur in a second clone
# ---------------------------------------------------------------------------


def test_repertoire_index_maps_kmers_to_distinct_cdr3s():
    from probe_designer.ext.tcr.repertoire import build_repertoire_index

    chain_a = "AAAACCCCGGGGTTTT" * 4          # 64 nt
    chain_b = "TTTTGGGGCCCCAAAA" + chain_a[16:]
    idx = build_repertoire_index(
        [("CDR3A", chain_a), ("CDR3B", chain_b)], kmer=40)
    # A 40-mer shared by both chains maps to two distinct CDR3s.
    shared = chain_a[24:64]
    assert shared in chain_b
    assert idx[shared] == {"CDR3A", "CDR3B"}


def test_repertoire_index_covers_both_strands():
    from probe_designer.ext.tcr.repertoire import build_repertoire_index

    chain = "ACGT" * 20
    idx = build_repertoire_index([("C1", chain)], kmer=40)
    fwd = chain[:40]
    rev = fwd.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    assert idx[fwd] == {"C1"}
    assert idx[rev] == {"C1"}


def test_repertoire_index_is_case_insensitive():
    from probe_designer.ext.tcr.repertoire import build_repertoire_index

    idx = build_repertoire_index([("C1", ("acgt" * 20))], kmer=40)
    assert ("ACGT" * 10) in idx


def _repertoire_phase1(tmp_path, repertoire):
    """Phase 1 for one clone whose CDR3 also occurs in a decoy chain."""
    import pandas as pd

    from probe_designer.ext.tcr.config import TcrConfig
    from probe_designer.ext.tcr.pipeline import _phase1_find_select

    cdr3 = "TGTGCCAGCAGTTTATGGGACACAGATACGCAGTATTTT"
    trb = ("ACCTGGAGCCCCCAGAACTGGCAGACACCTGCCTGATGCTGCCATGGGCCCCCAGCTCC"
           + cdr3
           + "GGCCCAGGCACCCGGCTGACAGTGCTCGAGGACCTGAAAAACGTGTTCCCACCC")
    # A companion clone the screen never touches, so the run stays non-empty
    # and ct1's fate shows up as a per-clone skip rather than a global error.
    other_cdr3 = "TGCAGTGCTCGTTCGGGGGGTGCCATGACTGAAGCTTTCTTT"
    other_trb = ("GACAGCAGCTTCTACATCTGCAGTGCTCGCACGTAGCATCGAAGCTTAGCTAAGCTT"
                 + other_cdr3
                 + "GGCAGTGGAACCCAGCTCTCTGTCTTGGAGGACCTGAACAAGGTGTTCCCACCCGAG")
    clones = pd.DataFrame([
        {"consensus_id": "ct1_consensus1", "chain_nt": trb,
         "cdr3_nt": cdr3, "Clone_name": "TEST_ct1"},
        {"consensus_id": "ctother_consensus1", "chain_nt": other_trb,
         "cdr3_nt": other_cdr3, "Clone_name": "TEST_ctother"},
    ])
    if repertoire is not None:
        repertoire = list(repertoire) + [(other_cdr3, other_trb)]
    cfg = TcrConfig(input_xlsx=tmp_path / "i.xlsx", output_dir=tmp_path,
                    backbone_file=tmp_path / "b.xlsx", start_no=1,
                    chemistries=["cDNA"], sites_per_clone=3,
                    min_ligation_margin=3, repertoire=repertoire,
                    plot_tm_landscape=False, skip_blast=True)
    run_dir = tmp_path / "rep_run"
    run_dir.mkdir(exist_ok=True)
    return trb, cdr3, _phase1_find_select(clones, cfg, run_dir)


def test_repertoire_none_keeps_every_site(tmp_path):
    _trb, _cdr3, (_all, selected, _skipped) = _repertoire_phase1(tmp_path, None)
    assert selected["cDNA"]["ct1_consensus1"]


def test_site_shared_with_another_clone_is_dropped(tmp_path):
    """A decoy chain carrying the whole CDR3 region under a DIFFERENT CDR3
    makes every site in that stretch non-specific, so the clone is skipped.

    The repertoire also lists the designed clone itself, as a real one would;
    its own entry must not be what disqualifies it.
    """
    trb, cdr3, _ = _repertoire_phase1(tmp_path, None)
    decoy = [(cdr3, trb), ("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", trb)]
    _t, _c, (_all, selected, skipped) = _repertoire_phase1(tmp_path, decoy)
    assert "ct1_consensus1" in skipped
    assert "ct1_consensus1" not in selected["cDNA"]


def test_own_repertoire_entry_alone_never_disqualifies(tmp_path):
    """The designed clone is normally IN the repertoire; that must be a no-op."""
    trb, cdr3, _ = _repertoire_phase1(tmp_path, None)
    _t, _c, (_all, selected, skipped) = _repertoire_phase1(
        tmp_path, [(cdr3, trb)])
    assert "ct1_consensus1" not in skipped
    assert selected["cDNA"]["ct1_consensus1"]


def test_repertoire_entry_with_the_same_cdr3_is_not_cross_reactivity(tmp_path):
    """The same clone filed twice (identical CDR3, differently trimmed chain)
    must not count against its own probes."""
    trb, cdr3, _ = _repertoire_phase1(tmp_path, None)
    same_clone_again = [(cdr3, trb), (cdr3, trb[5:])]
    _t, _c, (_all, selected, skipped) = _repertoire_phase1(
        tmp_path, same_clone_again)
    assert "ct1_consensus1" not in skipped
    assert selected["cDNA"]["ct1_consensus1"]


def test_load_repertoire_accepts_a_bare_sequence_column(tmp_path):
    """The 2026-09-02 exports name the chain column `sequence`."""
    import pandas as pd

    from probe_designer.ext.tcr.repertoire import load_repertoire_file

    path = tmp_path / "BZ10.csv"
    pd.DataFrame({"clonotype_id": ["BZ10_clonotype1"],
                  "consensus_id": ["clonotype1_consensus_1"],
                  "cdr3_nt": ["TGTGCC"], "sequence": ["ACGT" * 20],
                  "cdr3_in_sequence": ["yes"]}).to_csv(path, index=False)
    assert load_repertoire_file(path) == [("TGTGCC", "ACGT" * 20)]


def test_load_repertoire_prefers_a_named_chain_column_over_sequence(tmp_path):
    """`sequence` is generic, so a specific name must win when both exist."""
    import pandas as pd

    from probe_designer.ext.tcr.repertoire import load_repertoire_file

    path = tmp_path / "both.csv"
    pd.DataFrame({"cdr3_nt": ["TGTGCC"], "TRB_nt": ["AAAA"],
                  "sequence": ["CCCC"]}).to_csv(path, index=False)
    assert load_repertoire_file(path) == [("TGTGCC", "AAAA")]
