"""The thermal filter computes Tm at the configured reaction buffer.

Before 2026-07-17 every Tm used Biopython defaults (Na=50, Mg=0, 50 nM) and no
formamide. These tests pin that (a) Mg2+ raises the gated Tm, (b) formamide
lowers it, and (c) both reach Tm through the public thermal_filter, via the
ReactionConditions on FilterConfig.
"""
from __future__ import annotations

from probe_designer.config import BlastConfig, FilterConfig
from probe_designer.filtering import SequenceFilter

ARM = "GCATTCAGGTCACCTTGA"  # 18 nt, ~50% GC


def _tm(filter_config):
    f = SequenceFilter(filter_config, BlastConfig())
    # wide gate so we read the value, not the pass/fail
    res = f.thermal_filter(ARM, ARM, min_tm=0.0, max_tm=200.0)
    return res["tm_3prime"]


class TestBufferReachesTm:
    def test_magnesium_raises_tm(self):
        base = FilterConfig(monovalent_mM=50.0, mg_mM=0.0, strand_nM=25.0,
                            formamide_pct=0.0)
        with_mg = FilterConfig(monovalent_mM=50.0, mg_mM=10.0, strand_nM=25.0,
                               formamide_pct=0.0)
        assert _tm(with_mg) > _tm(base) + 5.0

    def test_formamide_lowers_tm(self):
        no_fmd = FilterConfig(monovalent_mM=50.0, mg_mM=0.0, strand_nM=25.0,
                              formamide_pct=0.0)
        with_fmd = FilterConfig(monovalent_mM=50.0, mg_mM=0.0, strand_nM=25.0,
                                formamide_pct=20.0, formamide_deg_per_pct=0.5)
        # 20% * 0.5 = 10 C depression
        assert _tm(with_fmd) == _tm(no_fmd) - 10.0

    def test_default_buffer_differs_from_biopython_default(self):
        # The whole point: default config no longer equals the naked Tm_NN call.
        from Bio.SeqUtils import MeltingTemp as mt
        from probe_designer.chemistry import dna_revcomp_to_rna

        naked = mt.Tm_NN(dna_revcomp_to_rna(ARM), nn_table=mt.R_DNA_NN1)
        buffered = _tm(FilterConfig())
        assert abs(buffered - naked) > 1.0
