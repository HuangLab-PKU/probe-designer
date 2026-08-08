"""ThermoPolicy — arm-Tm design rules shared by the gate and the ranker.

Audit R4 + R7. Two problems with one root cause: ``min_tm`` and ``max_tm_diff``
each existed twice — once in ``filtering`` (which rejects) and once in
``scoring`` (which ranks) — and the copies had drifted (scoring defaulted
``max_tm_diff=10`` and ``min_arm_tm=50`` independently of FilterConfig; the
config default ``min_tm`` was 45 against the YAMLs' 50). A rule that lives in
two places is a rule that eventually disagrees with itself.

ThermoPolicy owns the thresholds *and* the 0-1 scoring shapes, so the gate and
the ranker cannot diverge, and the composite weights stay in ``scoring`` where
they belong.

R4 proper: the balance tolerance tightens 10 -> 5 C (Larsson 2010, Krzywkowski
2017) but becomes a *soft* rule by default. Measured on 232 shipped probes, a
hard 5 C gate would reject 17.7% (19.7% of cDNA); the same 5 C used as a
scoring reference costs nothing and sharpens ranking exactly where field
practice says the discrimination matters.
"""
from __future__ import annotations

import pytest

from probe_designer.chemistry import ReactionConditions
from probe_designer.config import FilterConfig
from probe_designer.policy import ThermoPolicy


class TestDefaults:
    def test_balance_tolerance_follows_field_practice(self):
        assert ThermoPolicy().max_tm_diff == 5.0

    def test_both_tm_rules_are_soft_by_default(self):
        policy = ThermoPolicy()
        assert policy.enforce_tm_gate is False
        assert policy.enforce_tm_diff_gate is False


class TestResolveFromConfig:
    def test_maps_filter_config_field_names(self):
        cfg = FilterConfig(min_tm=52.0, max_tm=68.0, max_tm_diff=4.0)
        policy = ThermoPolicy.resolve(cfg)
        assert (policy.min_arm_tm, policy.max_arm_tm) == (52.0, 68.0)
        assert policy.max_tm_diff == 4.0

    def test_per_call_overrides_win(self):
        cfg = FilterConfig(max_tm_diff=5.0)
        policy = ThermoPolicy.resolve(cfg, max_tm_diff=0.5, enforce_tm_diff_gate=True)
        assert policy.max_tm_diff == 0.5
        assert policy.enforce_tm_diff_gate is True

    def test_none_override_does_not_clobber_the_config(self):
        """kwargs.get(...) returning None must not wipe a real config value."""
        cfg = FilterConfig(max_tm_diff=7.0)
        assert ThermoPolicy.resolve(cfg, max_tm_diff=None).max_tm_diff == 7.0


class TestFromReaction:
    def test_target_tm_tracks_the_reaction(self):
        """Change the hyb temperature and the Tm target moves with it."""
        cold = ThermoPolicy.from_reaction(ReactionConditions(lab_temp_c=45.0))
        hot = ThermoPolicy.from_reaction(ReactionConditions(lab_temp_c=60.0))
        assert cold.min_arm_tm == 50.0          # 45 + 5 margin
        assert hot.min_arm_tm == 65.0
        assert hot.max_arm_tm > cold.max_arm_tm


class TestBalanceScoring:
    def test_perfect_balance_scores_one(self):
        assert ThermoPolicy().tm_balance_score(0.0) == pytest.approx(1.0)

    def test_score_reaches_zero_at_the_tolerance(self):
        assert ThermoPolicy(max_tm_diff=5.0).tm_balance_score(5.0) == pytest.approx(0.0)

    def test_score_is_clamped_not_negative_beyond_tolerance(self):
        assert ThermoPolicy(max_tm_diff=5.0).tm_balance_score(50.0) == 0.0

    def test_tightening_the_tolerance_sharpens_discrimination(self):
        """The R4 point: 5 C separates probes that 10 C treats as similar."""
        loose, tight = ThermoPolicy(max_tm_diff=10.0), ThermoPolicy(max_tm_diff=5.0)
        assert tight.tm_balance_score(2.0) < loose.tm_balance_score(2.0)

    def test_zero_tolerance_disables_the_term(self):
        assert ThermoPolicy(max_tm_diff=0.0).tm_balance_score(3.0) == 0.0


class TestRangeScoring:
    def test_peaks_at_the_target(self):
        policy = ThermoPolicy(min_arm_tm=50.0)
        assert policy.tm_range_score(50.0) == pytest.approx(1.0)

    def test_two_sided_falloff(self):
        """Too cold is as bad as too hot — the Phase 1 two-sided shape."""
        policy = ThermoPolicy(min_arm_tm=50.0, tm_falloff_c=20.0)
        assert policy.tm_range_score(40.0) == pytest.approx(policy.tm_range_score(60.0))
        assert policy.tm_range_score(60.0) == pytest.approx(0.5)

    def test_clamped_at_zero(self):
        assert ThermoPolicy(min_arm_tm=50.0).tm_range_score(200.0) == 0.0


class TestGateSemantics:
    def test_soft_by_default_nothing_is_rejected(self):
        policy = ThermoPolicy(max_tm_diff=5.0)
        assert policy.tm_balance_ok(50.0) is True
        assert policy.tm_range_ok(5.0) is True

    def test_hard_gate_rejects_when_enabled(self):
        policy = ThermoPolicy(
            min_arm_tm=50.0, max_arm_tm=70.0, max_tm_diff=5.0,
            enforce_tm_gate=True, enforce_tm_diff_gate=True,
        )
        assert policy.tm_balance_ok(6.0) is False
        assert policy.tm_balance_ok(4.0) is True
        assert policy.tm_range_ok(45.0) is False
        assert policy.tm_range_ok(75.0) is False
        assert policy.tm_range_ok(60.0) is True

    def test_gate_boundary_is_inclusive(self):
        policy = ThermoPolicy(
            min_arm_tm=50.0, max_arm_tm=70.0, max_tm_diff=5.0,
            enforce_tm_gate=True, enforce_tm_diff_gate=True,
        )
        assert policy.tm_balance_ok(5.0) is True
        assert policy.tm_range_ok(50.0) is True
        assert policy.tm_range_ok(70.0) is True


class TestSingleSourceOfTruth:
    """The regression that motivated the object: gate and ranker agreeing."""

    def test_scorer_and_filter_read_the_same_tolerance(self):
        from probe_designer.scoring import compute_target_score

        cfg = FilterConfig(max_tm_diff=5.0, enforce_tm_diff_gate=True)
        policy = ThermoPolicy.resolve(cfg)

        # A probe exactly at the tolerance: scores 0 on balance, still passes.
        site = {"arm_5prime": "ATATATATAT", "arm_3prime": "ATATATATAT",
                "tm_5prime": 47.5, "tm_3prime": 52.5,
                "isoform_overlap_num": 1, "blast_alignments": [], "free_energy": 0.0}
        assert policy.tm_balance_score(5.0) == pytest.approx(0.0)
        assert policy.tm_balance_ok(5.0) is True
        # and the composite reflects exactly that zero, not some other tolerance
        with_balance = compute_target_score(site, policy=ThermoPolicy(max_tm_diff=100.0))
        without = compute_target_score(site, policy=policy)
        assert with_balance > without

    def test_filter_gate_and_scorer_agree_on_the_same_config(self):
        """One config in, one tolerance out — for the gate and the ranker alike."""
        from probe_designer.config import BlastConfig
        from probe_designer.filtering import SequenceFilter

        cfg = FilterConfig(max_tm_diff=5.0, enforce_tm_diff_gate=True)
        policy = ThermoPolicy.resolve(cfg)
        filt = SequenceFilter(cfg, BlastConfig())

        # 5' arm ~ AT-rich (low Tm), 3' arm GC-rich (high Tm) -> large tm_diff.
        result = filt.thermal_filter(
            arm_3prime="GCGCGCGCGCGCGCGCGCGC",
            arm_5prime="ATATATATATATATATATAT",
            sequence_type="DNA", target_type="DNA", target_sequence=None,
        )
        assert result["tm_diff"] > policy.max_tm_diff
        # The gate rejected it, and the scorer independently scores it 0.
        assert "tm_diff" in result["failed_checks"]
        assert policy.tm_balance_score(result["tm_diff"]) == 0.0
