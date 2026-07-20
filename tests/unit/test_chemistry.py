"""Tests for probe_designer.chemistry — the ReactionConditions buffer object.

Pins the buffer defaults (current canonical protocol rca.md v5.3: 75 mM K+,
10 mM Mg2+, 0.1 uM strand, 20% formamide, 45 C anneal) and the derived
quantities that the Tm/ΔG call sites depend on. See
experiments/20260715_tm_deltag_methods_audit/ for the motivating audit.
"""
from __future__ import annotations

import math

import pytest
from Bio.SeqUtils import MeltingTemp as mt

from probe_designer.chemistry import (
    ReactionConditions,
    dna_revcomp_to_rna,
    register_solvent,
)


# ---------------------------------------------------------------------------
# dna_revcomp_to_rna — the shared RNA-sense helper (also the strand fix)
# ---------------------------------------------------------------------------

class TestDnaRevcompToRna:
    def test_reverse_complement_keeps_thymine(self):
        # RC of the DNA probe arm = the RNA-sense (mRNA) strand; Biopython's
        # R_DNA_NN1 is keyed on T, so T is intentionally NOT transcribed to U.
        assert dna_revcomp_to_rna("AAAC") == "GTTT"

    def test_is_reverse_complement_not_just_complement(self):
        assert dna_revcomp_to_rna("ACGG") == "CCGT"

    def test_uppercases_and_handles_n(self):
        assert dna_revcomp_to_rna("acgn") == "NCGT"

    def test_roundtrip_double_rc_is_identity(self):
        arm = "GATCGGATCCATTGCA"
        assert dna_revcomp_to_rna(dna_revcomp_to_rna(arm)) == arm


# ---------------------------------------------------------------------------
# ReactionConditions defaults + validation
# ---------------------------------------------------------------------------

class TestReactionConditionsDefaults:
    def test_protocol_defaults(self):
        rc = ReactionConditions()
        assert rc.monovalent_mM == 75.0
        assert rc.mg_mM == 10.0
        assert rc.dntp_mM == 0.0
        assert rc.strand_nM == 100.0
        assert rc.formamide_pct == 20.0
        assert rc.formamide_deg_per_pct == 0.5
        assert rc.lab_temp_c == 45.0
        assert rc.tm_margin_c == 5.0
        assert rc.saltcorr == 5

    def test_min_arm_tm_is_lab_temp_plus_margin(self):
        assert ReactionConditions().min_arm_tm == 50.0
        assert ReactionConditions(lab_temp_c=37.0, tm_margin_c=8.0).min_arm_tm == 45.0

    def test_effective_celsius_adds_formamide_offset(self):
        # 45 + 20 * 0.5 = 55 (matches the cross-lig / NUPACK screen model)
        assert ReactionConditions().effective_celsius == 55.0


class TestReactionConditionsValidation:
    @pytest.mark.parametrize("field,value", [
        ("monovalent_mM", -1.0),
        ("mg_mM", -1.0),
        ("dntp_mM", -1.0),
        ("formamide_pct", -1.0),
    ])
    def test_negative_concentrations_rejected(self, field, value):
        with pytest.raises(ValueError):
            ReactionConditions(**{field: value})

    def test_nonpositive_strand_rejected(self):
        with pytest.raises(ValueError):
            ReactionConditions(strand_nM=0.0)

    def test_zero_monovalent_rejected(self):
        # salt_correction methods 1-6 divide by ln(monovalent); require > 0.
        with pytest.raises(ValueError):
            ReactionConditions(monovalent_mM=0.0)

    def test_bad_saltcorr_rejected(self):
        with pytest.raises(ValueError):
            ReactionConditions(saltcorr=9)


# ---------------------------------------------------------------------------
# Helper payloads for the thermodynamic engines
# ---------------------------------------------------------------------------

class TestTmNnKwargs:
    def test_monovalent_mapped_to_potassium_not_sodium(self):
        # Protocol buffer is KCl-based; Na stays 0 so von Ahsen sums K correctly.
        kw = ReactionConditions().tm_nn_kwargs()
        assert kw["Na"] == 0.0
        assert kw["K"] == 75.0
        assert kw["Mg"] == 10.0
        assert kw["dNTPs"] == 0.0
        assert kw["dnac1"] == 100.0 and kw["dnac2"] == 100.0
        assert kw["saltcorr"] == 5

    def test_kwargs_make_tm_nn_run_and_raise_tm_vs_biopython_default(self):
        # Passing 10 mM Mg (via von Ahsen in saltcorr=5) must raise Tm vs the
        # Na=50/Mg=0 Biopython default — the whole point of the fix.
        seq = "GCATTCAGGTCACCTTGA"  # RNA-sense strand, R_DNA_NN1
        default_tm = mt.Tm_NN(seq, nn_table=mt.R_DNA_NN1)
        buffered_tm = mt.Tm_NN(seq, nn_table=mt.R_DNA_NN1,
                               **ReactionConditions().tm_nn_kwargs())
        assert buffered_tm > default_tm + 5.0


class TestFormamideCorrection:
    def test_formamide_lowers_tm_linearly(self):
        rc = ReactionConditions()  # 20% * 0.5 = 10 C depression
        assert rc.apply_formamide(60.0) == pytest.approx(50.0)

    def test_zero_formamide_is_noop(self):
        rc = ReactionConditions(formamide_pct=0.0)
        assert rc.apply_formamide(60.0) == pytest.approx(60.0)


class TestFromBuffer:
    def test_granular_monovalents_fold_to_total(self):
        # Na/K/Tris are equivalent for Tm; total = Na + K + Tris/2.
        rc = ReactionConditions.from_buffer(na_mM=25.0, k_mM=50.0)
        assert rc.monovalent_mM == 75.0

    def test_tris_counts_half(self):
        rc = ReactionConditions.from_buffer(k_mM=50.0, tris_mM=50.0)
        assert rc.monovalent_mM == 75.0  # 50 + 50/2

    def test_other_kwargs_pass_through(self):
        rc = ReactionConditions.from_buffer(k_mM=75.0, mg_mM=4.0, formamide_pct=0.0)
        assert rc.mg_mM == 4.0 and rc.formamide_pct == 0.0


class TestSolvents:
    def test_default_solvents_empty_matches_formamide_only(self):
        rc = ReactionConditions()
        # 20% formamide * 0.5 = 10 C, no extra solvents.
        assert rc.solvent_tm_depression == pytest.approx(10.0)
        assert rc.effective_celsius == pytest.approx(55.0)

    def test_extra_solvent_adds_depression(self):
        rc = ReactionConditions(solvents={"dmso": 10.0})  # +0.6*10 = 6
        assert rc.solvent_tm_depression == pytest.approx(16.0)
        assert rc.apply_solvents(60.0) == pytest.approx(44.0)
        assert rc.effective_celsius == pytest.approx(61.0)

    def test_apply_formamide_alias_still_works(self):
        assert ReactionConditions().apply_formamide(60.0) == pytest.approx(50.0)

    def test_unknown_solvent_raises_until_registered(self):
        rc = ReactionConditions(solvents={"betaine": 20.0})
        with pytest.raises(ValueError, match="register_solvent"):
            _ = rc.solvent_tm_depression
        register_solvent("betaine", 0.0)  # isostabilizer ~ no net Tm shift here
        assert rc.solvent_tm_depression == pytest.approx(10.0)  # formamide only

    def test_negative_solvent_percent_rejected(self):
        with pytest.raises(ValueError):
            ReactionConditions(solvents={"dmso": -1.0})

    def test_signature_distinguishes_solvents(self):
        assert (ReactionConditions().signature()
                != ReactionConditions(solvents={"dmso": 5.0}).signature())


class TestEngineKwargs:
    def test_primer3_kwargs(self):
        kw = ReactionConditions().primer3_kwargs()
        assert kw["mv_conc"] == 75.0
        assert kw["dv_conc"] == 10.0
        assert kw["dntp_conc"] == 0.0
        assert kw["dna_conc"] == 100.0   # nM
        assert kw["temp_c"] == 55.0      # effective (formamide-shifted)

    def test_nupack_kwargs_molar(self):
        kw = ReactionConditions().nupack_kwargs()
        assert kw["sodium"] == pytest.approx(0.075)
        assert kw["magnesium"] == pytest.approx(0.010)
        assert kw["celsius"] == 55.0
