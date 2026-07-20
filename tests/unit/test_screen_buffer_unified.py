"""The cross-lig / NUPACK screen buffers derive from ReactionConditions (P5).

Locks the single-source-of-truth: the screen constants must equal the shared
ReactionConditions defaults, and keep their historical numeric values so this
unification is behavior-preserving.
"""
from __future__ import annotations

import pytest

from probe_designer.chemistry import ReactionConditions

_LAB = ReactionConditions()


class TestCrossLigationBufferDerived:
    def test_matches_reaction_conditions(self):
        from probe_designer.qc import cross_ligation as xl
        assert xl.DEFAULT_MV_CONC_MM == _LAB.monovalent_mM
        assert xl.DEFAULT_DV_CONC_MM == _LAB.mg_mM
        assert xl.DEFAULT_DNA_CONC_M == _LAB.strand_nM * 1e-9
        assert xl.DEFAULT_TEMP_C == _LAB.effective_celsius

    def test_historical_values_unchanged(self):
        from probe_designer.qc import cross_ligation as xl
        assert xl.DEFAULT_MV_CONC_MM == 75.0
        assert xl.DEFAULT_DV_CONC_MM == 10.0
        assert xl.DEFAULT_DNA_CONC_M == pytest.approx(1.0e-7)
        assert xl.DEFAULT_TEMP_C == 55.0


class TestNupackBufferDerived:
    def test_matches_reaction_conditions(self):
        from probe_designer.ext.nupack import config as nc
        nk = _LAB.nupack_kwargs()
        assert nc.DEFAULT_SODIUM_M == nk["sodium"]
        assert nc.DEFAULT_MAGNESIUM_M == nk["magnesium"]
        assert nc.DEFAULT_CELSIUS == nk["celsius"]
        assert nc.DEFAULT_STRAND_CONC_M == _LAB.strand_nM * 1e-9

    def test_historical_values_unchanged(self):
        from probe_designer.ext.nupack import config as nc
        assert nc.DEFAULT_SODIUM_M == pytest.approx(0.075)
        assert nc.DEFAULT_MAGNESIUM_M == pytest.approx(0.010)
        assert nc.DEFAULT_CELSIUS == 55.0
        assert nc.DEFAULT_STRAND_CONC_M == pytest.approx(1.0e-7)
