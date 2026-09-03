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
    """The screen holds the buffer as one object, not four loose scalars.

    The four ``DEFAULT_*_CONC`` constants were restatements of fields on the same
    ``ReactionConditions`` instance, which is what this test existed to police.
    Now the screen takes a ``reaction`` and defaults it to ``LAB_HYBRIDISATION``,
    so drift is structurally impossible — what is left worth locking is that the
    default really is the lab protocol, at its historical numbers.
    """

    def test_screen_default_is_the_lab_reaction(self):
        from probe_designer.qc import cross_ligation as xl
        assert xl.LAB_HYBRIDISATION == _LAB

    def test_historical_values_unchanged(self):
        from probe_designer.qc import cross_ligation as xl
        lab = xl.LAB_HYBRIDISATION
        assert lab.monovalent_mM == 75.0
        assert lab.mg_mM == 10.0
        assert lab.strand_nM * 1e-9 == pytest.approx(1.0e-7)
        assert lab.effective_celsius == 55.0

    def test_clamp_is_the_shared_definition(self):
        """One ligase clamp across the register scan, BLAST and NUPACK."""
        from probe_designer.blast import gene_aware
        from probe_designer.chemistry import LIGASE_CLAMP_NT
        from probe_designer.ext.nupack import config as nc
        from probe_designer.qc import cross_ligation as xl
        assert xl.DEFAULT_VICINITY_N == LIGASE_CLAMP_NT
        assert nc.DEFAULT_VICINITY_N == LIGASE_CLAMP_NT
        assert gene_aware.DEFAULT_CLAMP_NT == LIGASE_CLAMP_NT


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
