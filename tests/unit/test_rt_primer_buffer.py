"""RT primer Tm: strand-fixed (RNA-sense) and computed at the Maxima H Minus buffer."""
from __future__ import annotations

from Bio.SeqUtils import MeltingTemp as mt

from probe_designer.chemistry import ReactionConditions, dna_revcomp_to_rna
from probe_designer.rt_primer import MAXIMA_H_MINUS_RT, _compute_tm


def test_maxima_buffer_composition():
    rc = MAXIMA_H_MINUS_RT
    assert rc.monovalent_mM == 75.0   # 1X KCl
    assert rc.mg_mM == 3.0            # 1X MgCl2
    assert rc.formamide_pct == 0.0   # no formamide in RT
    assert rc.lab_temp_c == 50.0     # Maxima H Minus reaction temp


def test_compute_tm_uses_rna_sense_strand_and_buffer():
    # The primer is DNA; R_DNA_NN1 must see the RNA template = RC(primer).
    primer = "GCATTCAGGTCACCTTGATG"
    rc = MAXIMA_H_MINUS_RT
    expected = round(
        mt.Tm_NN(dna_revcomp_to_rna(primer), nn_table=mt.R_DNA_NN1, **rc.tm_nn_kwargs()),
        1,
    )
    assert _compute_tm(primer, rc) == expected


def test_compute_tm_differs_from_dna_strand_convention():
    # Regression guard for the strand bug: DNA-strand Tm != RNA-sense Tm.
    primer = "GCATTCAGGTCACCTTGATG"
    rc = MAXIMA_H_MINUS_RT
    dna_strand_tm = round(
        mt.Tm_NN(primer, nn_table=mt.R_DNA_NN1, **rc.tm_nn_kwargs()), 1
    )
    assert _compute_tm(primer, rc) != dna_strand_tm
