"""Banerjee 2020 hybrid NN parameters + the one place that picks an NN table.

Audit R6. Two classes of thing are tested here, and the first matters more than
the second: a wrong thermodynamic parameter silently produces wrong Tm for every
probe the lab orders, so the table is validated against the paper's own internal
relations rather than merely transcribed.

1. **Table integrity** — ΔG°37 = ΔH° − 310.15·ΔS° for all 16 dimers and both
   initiation terms (catches a transcription slip in any single number), and the
   key set maps exactly onto Biopython's, which pins the notation conversion
   (the paper writes both strands 5'->3'; Biopython wants the DNA strand 3'->5').

2. **Salt reference** — Banerjee's parameters were measured *at* 100 mM NaCl, so
   Biopython's 1 M-referenced salt correction must not be applied on top of them.
"""
from __future__ import annotations

import pytest
from Bio.SeqUtils import MeltingTemp as mt

from probe_designer.chemistry import ReactionConditions
from probe_designer.nn_tables import (
    DEFAULT_HYBRID_MODEL,
    HYBRID_BANERJEE2020,
    HYBRID_SUGIMOTO1995,
    R_DNA_BANERJEE2020,
    melting_temperature,
    nn_model_for,
    nn_model_for_chemistry,
)

T37 = 310.15

# (Biopython key, published ΔG°37) from Banerjee 2020 Table 2 as corrected in
# NAR 2021 49:10796. ΔH/ΔS live in the module; ΔG is the independent check.
PUBLISHED_DG37 = {
    "AA/TT": -0.7, "AC/TG": -1.5, "AG/TC": -1.3, "AT/TA": -0.4,
    "CA/GT": -1.2, "CC/GG": -1.7, "CG/GC": -1.4, "CT/GA": -0.4,
    "GA/CT": -1.5, "GC/CG": -2.0, "GG/CC": -2.3, "GT/CA": -1.4,
    "TA/AT": -0.5, "TC/AG": -1.4, "TG/AC": -1.6, "TT/AA": 0.2,
    "init_G/C": 2.0, "init_A/T": 2.6,
}

ARMS = [
    "GATGATGATGATGATGATGC", "GCTAGCTAGCTAGCTAGCTA", "ATGCATGCATGCATGCATGC",
    "GGCCGGCCGGAAGGCCGGCC", "ATATATATATCGCGATATAT", "TTGACCAGTGCATTCGAAGT",
]


class TestTableIntegrity:
    @pytest.mark.parametrize("key,dg_published", sorted(PUBLISHED_DG37.items()))
    def test_gibbs_relation_holds(self, key, dg_published):
        """ΔG°37 = ΔH° − T·ΔS° — catches a slip in any single transcribed value."""
        dh, ds = R_DNA_BANERJEE2020[key]
        assert dh - T37 * ds / 1000.0 == pytest.approx(dg_published, abs=0.02)

    def test_key_set_matches_biopython_hybrid_table(self):
        """Pins the notation conversion, not just the numbers."""
        ours = {k for k in R_DNA_BANERJEE2020 if "/" in k and not k.startswith("init")}
        theirs = {k for k in mt.R_DNA_NN1 if "/" in k and not k.startswith("init")}
        assert ours == theirs

    def test_uses_the_corrected_initiation_entropies(self):
        """The 2021 correction changed these; the originals fail the Gibbs check."""
        assert R_DNA_BANERJEE2020["init_G/C"] == (0.0, -6.4)   # not -4.9
        assert R_DNA_BANERJEE2020["init_A/T"] == (0.0, -8.4)   # not -7.0

    def test_every_dimer_is_stabilizing(self):
        for key, (dh, ds) in R_DNA_BANERJEE2020.items():
            if "/" in key and not key.startswith("init"):
                assert dh < 0 and ds < 0, f"{key} is not a stabilizing NN step"


class TestSaltReference:
    def test_banerjee_declares_its_measurement_salt(self):
        assert HYBRID_BANERJEE2020.reference_monovalent_mM == 100.0
        assert HYBRID_SUGIMOTO1995.reference_monovalent_mM is None  # 1 M convention

    def test_no_double_correction_at_the_reference_condition(self):
        """At 100 mM monovalent / no Mg the parameters need no salt adjustment."""
        reaction = ReactionConditions(
            monovalent_mM=100.0, mg_mM=0.0, formamide_pct=0.0, saltcorr=5,
        )
        for arm in ARMS:
            raw = mt.Tm_NN(arm, nn_table=R_DNA_BANERJEE2020,
                           **{**reaction.tm_nn_kwargs(), "saltcorr": 0})
            assert melting_temperature(arm, HYBRID_BANERJEE2020, reaction) == pytest.approx(raw, abs=0.01)

    def test_disagreement_with_sugimoto_matches_the_published_magnitude(self):
        """Banerjee reports Sugimoto errs ~4.9 C at 100 mM. Applying Biopython's
        1 M-referenced correction on top of Banerjee's 100 mM parameters instead
        puts them ~14.5 C apart — the bug this reference handling exists to
        prevent, and the reason a raw table swap would have been wrong."""
        reaction = ReactionConditions(
            monovalent_mM=100.0, mg_mM=0.0, formamide_pct=0.0, saltcorr=5,
        )
        diffs = [
            melting_temperature(a, HYBRID_BANERJEE2020, reaction)
            - melting_temperature(a, HYBRID_SUGIMOTO1995, reaction)
            for a in ARMS
        ]
        mean_diff = sum(diffs) / len(diffs)
        assert -8.0 < mean_diff < 0.0, f"mean {mean_diff:.2f} C is not the published magnitude"

    def test_added_magnesium_raises_tm_relative_to_the_reference(self):
        """The salt model is still used — for the shift from 100 mM, not from 1 M."""
        bare = ReactionConditions(monovalent_mM=100.0, mg_mM=0.0, formamide_pct=0.0)
        with_mg = ReactionConditions(monovalent_mM=100.0, mg_mM=10.0, formamide_pct=0.0)
        assert (melting_temperature(ARMS[0], HYBRID_BANERJEE2020, with_mg)
                > melting_temperature(ARMS[0], HYBRID_BANERJEE2020, bare))


class TestModelSelection:
    def test_sugimoto_remains_the_default(self):
        assert DEFAULT_HYBRID_MODEL == "sugimoto1995"
        assert ReactionConditions().hybrid_nn_model == "sugimoto1995"
        assert nn_model_for("DNA", "RNA").table is mt.R_DNA_NN1

    @pytest.mark.parametrize("seq_t,tgt_t,expected", [
        ("DNA", "RNA", mt.R_DNA_NN1),
        ("DNA", "DNA", mt.DNA_NN4),
        ("RNA", "RNA", mt.RNA_NN1),
    ])
    def test_strand_pair_selects_the_right_table(self, seq_t, tgt_t, expected):
        assert nn_model_for(seq_t, tgt_t).table is expected

    @pytest.mark.parametrize("chemistry,expected", [
        ("dRNA", mt.R_DNA_NN1), ("iLock", mt.R_DNA_NN1),
        ("cDNA", mt.DNA_NN4), ("CDNA", mt.DNA_NN4),
    ])
    def test_chemistry_label_selects_the_right_table(self, chemistry, expected):
        assert nn_model_for_chemistry(chemistry).table is expected

    def test_hybrid_model_is_selectable(self):
        assert nn_model_for("DNA", "RNA", "banerjee2020").table is R_DNA_BANERJEE2020

    def test_unknown_model_fails_fast(self):
        """Never silently fall back to different physics."""
        with pytest.raises(ValueError, match="unknown hybrid_nn_model"):
            nn_model_for("DNA", "RNA", "sugimoto1996")
        with pytest.raises(ValueError, match="hybrid_nn_model"):
            ReactionConditions(hybrid_nn_model="nope")

    def test_model_choice_does_not_leak_into_non_hybrid_chemistries(self):
        """cDNA is DNA:DNA — the hybrid selection must not touch it."""
        assert nn_model_for("DNA", "DNA", "banerjee2020").table is mt.DNA_NN4
