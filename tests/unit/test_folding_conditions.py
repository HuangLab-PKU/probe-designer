"""FoldingConditions — RNAplfold geometry as a first-class, reaction-linked object.

Audit R9 (Lange 2012, NAR 40:5215): the ViennaRNA default W = L = 70 is
artifact-prone for accessibility — too short a span truncates long-range pairing
and systematically over-estimates how open a site is. The recommendation is
L ~ 100 with W = L + 50.

The window/span/temperature triple was previously three loose numbers repeated
across FilterConfig, search_strategies, annotate, accessibility and the YAMLs —
the same duplication the buffer had before ReactionConditions. These tests pin
the object, its validation, and the rule that folding temperature *tracks the
reaction* unless explicitly overridden (so the accessibility used to filter is
the accessibility written into the reference annotation).
"""
from __future__ import annotations

import pytest

from probe_designer.chemistry import FoldingConditions, ReactionConditions
from probe_designer.config import FilterConfig


class TestDefaults:
    def test_defaults_follow_lange_2012(self):
        fc = FoldingConditions()
        assert fc.span == 100
        assert fc.window == fc.span + 50 == 150

    def test_temperature_defaults_to_tracking_the_reaction(self):
        assert FoldingConditions().temperature_c is None


class TestValidation:
    def test_span_above_window_is_rejected(self):
        with pytest.raises(ValueError, match="span"):
            FoldingConditions(window=50, span=100)

    @pytest.mark.parametrize("kwargs", [{"window": 0}, {"span": 0}, {"window": -10}])
    def test_non_positive_geometry_is_rejected(self, kwargs):
        with pytest.raises(ValueError):
            FoldingConditions(**kwargs)


class TestPlfoldKwargs:
    def test_temperature_tracks_the_reaction_when_unset(self):
        rc = ReactionConditions(lab_temp_c=45.0, formamide_pct=20.0)
        kw = FoldingConditions().plfold_kwargs(rc)
        # 45 lab + 20% * 0.5 C/% formamide = 55 C effective
        assert kw["temperature"] == pytest.approx(rc.effective_celsius)
        assert kw["temperature"] == pytest.approx(55.0)

    def test_explicit_temperature_overrides_the_reaction(self):
        rc = ReactionConditions(lab_temp_c=45.0, formamide_pct=20.0)
        kw = FoldingConditions(temperature_c=37.0).plfold_kwargs(rc)
        assert kw["temperature"] == pytest.approx(37.0)

    def test_without_a_reaction_falls_back_to_the_viennarna_convention(self):
        kw = FoldingConditions().plfold_kwargs()
        assert kw["temperature"] == pytest.approx(37.0)

    def test_kwargs_match_compute_plfold_profile_signature(self):
        """The mapping must be splat-able straight into the primitive."""
        import inspect

        from probe_designer.filtering.accessibility import compute_plfold_profile

        params = inspect.signature(compute_plfold_profile).parameters
        for key in FoldingConditions().plfold_kwargs():
            assert key in params, f"{key} is not a compute_plfold_profile parameter"


class TestFilterConfigIntegration:
    def test_filter_config_defaults_are_lange_compliant(self):
        cfg = FilterConfig()
        assert (cfg.plfold_window, cfg.plfold_span) == (150, 100)

    def test_filter_config_builds_folding_conditions(self):
        cfg = FilterConfig(plfold_window=120, plfold_span=80, plfold_temperature=42.0)
        fc = cfg.folding_conditions()
        assert (fc.window, fc.span, fc.temperature_c) == (120, 80, 42.0)

    def test_accessibility_module_constants_track_the_object(self):
        """One source of truth: the module defaults are not a second copy."""
        from probe_designer.filtering import accessibility

        fc = FoldingConditions()
        assert accessibility.DEFAULT_WINDOW == fc.window
        assert accessibility.DEFAULT_SPAN == fc.span
