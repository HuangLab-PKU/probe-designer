"""Per-position Tm profile / reference-annotation cache."""
from __future__ import annotations

import numpy as np
from Bio.SeqUtils import MeltingTemp as mt

from probe_designer.chemistry import ReactionConditions
from probe_designer.filtering.thermo_profile import (
    compute_tm_profile,
    cached_tm_profile,
)

SEQ = "GCATTCAGGTCACCTTGATGCATTCAGGTCACCTTGATG"  # 39 nt, RNA-sense reference


class TestSignature:
    def test_buffer_signature_changes_with_mg(self):
        assert ReactionConditions().signature() != ReactionConditions(mg_mM=4.0).signature()

    def test_signature_ignores_margin_and_temp(self):
        # tm_margin_c / lab_temp_c do not change the Tm value -> not in the key.
        a = ReactionConditions().signature()
        b = ReactionConditions(lab_temp_c=37.0, tm_margin_c=8.0).signature()
        assert a == b


class TestComputeTmProfile:
    def test_shape_and_trailing_nan(self):
        prof = compute_tm_profile(SEQ, ReactionConditions(), arm_length=20)
        assert prof.shape == (len(SEQ),)
        assert np.isnan(prof[-1])                 # no full 20-mer at the end
        assert not np.isnan(prof[0])

    def test_drna_value_matches_direct_hybrid_tm(self):
        rc = ReactionConditions()
        prof = compute_tm_profile(SEQ, rc, arm_length=20, chemistry="dRNA")
        # For dRNA the reference window IS the RNA-sense strand for R_DNA_NN1.
        expected = rc.apply_formamide(
            mt.Tm_NN(SEQ[0:20], nn_table=mt.R_DNA_NN1, **rc.tm_nn_kwargs())
        )
        assert prof[0] == expected

    def test_cdna_uses_dna_table(self):
        rc = ReactionConditions()
        prof = compute_tm_profile(SEQ, rc, arm_length=20, chemistry="cDNA")
        expected = rc.apply_formamide(
            mt.Tm_NN(SEQ[0:20], nn_table=mt.DNA_NN4, **rc.tm_nn_kwargs())
        )
        assert prof[0] == expected

    def test_mg_raises_profile(self):
        no_mg = compute_tm_profile(SEQ, ReactionConditions(mg_mM=0.0, formamide_pct=0.0))
        with_mg = compute_tm_profile(SEQ, ReactionConditions(mg_mM=10.0, formamide_pct=0.0))
        assert with_mg[0] > no_mg[0]


class TestCachedTmProfile:
    def test_writes_and_reloads(self, tmp_path):
        rc = ReactionConditions()
        a = cached_tm_profile(SEQ, "TX1", rc, cache_dir=tmp_path, arm_length=20)
        files = list(tmp_path.glob("*.npy"))
        assert len(files) == 1
        b = cached_tm_profile(SEQ, "TX1", rc, cache_dir=tmp_path, arm_length=20)
        assert np.allclose(a, b, equal_nan=True)
        assert len(list(tmp_path.glob("*.npy"))) == 1   # reused, not regenerated

    def test_different_buffer_is_different_annotation(self, tmp_path):
        cached_tm_profile(SEQ, "TX1", ReactionConditions(), cache_dir=tmp_path)
        cached_tm_profile(SEQ, "TX1", ReactionConditions(mg_mM=4.0), cache_dir=tmp_path)
        assert len(list(tmp_path.glob("*.npy"))) == 2   # keyed by buffer signature
