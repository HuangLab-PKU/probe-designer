"""The design CLI --*-buffer flags override FilterConfig and re-anchor Tm."""
from __future__ import annotations

import pytest

from probe_designer.cli.design import _apply_buffer_overrides
from probe_designer.config import FilterConfig

_NONE = dict(monovalent_mm=None, mg_mm=None, dntp_mm=None, strand_nm=None,
             formamide_pct=None, lab_temp_c=None, tm_margin_c=None, saltcorr=None)


def test_non_none_values_override():
    fc = FilterConfig()
    _apply_buffer_overrides(fc, **{**_NONE, "mg_mm": 4.0, "formamide_pct": 0.0})
    assert fc.mg_mM == 4.0
    assert fc.formamide_pct == 0.0
    assert fc.monovalent_mM == 75.0  # untouched


def test_lab_temp_reanchors_window():
    fc = FilterConfig()
    _apply_buffer_overrides(fc, **{**_NONE, "lab_temp_c": 37.0, "tm_margin_c": 8.0})
    assert fc.min_tm == 45.0   # 37 + 8
    assert fc.max_tm == 62.0   # 37 + 25


def test_window_untouched_without_temp_flags():
    fc = FilterConfig()
    _apply_buffer_overrides(fc, **{**_NONE, "mg_mm": 4.0})
    assert fc.min_tm == 50.0   # default anchor unchanged
    assert fc.max_tm == 70.0


def test_bad_buffer_fails_fast():
    fc = FilterConfig()
    with pytest.raises(ValueError):
        _apply_buffer_overrides(fc, **{**_NONE, "mg_mm": -1.0})
