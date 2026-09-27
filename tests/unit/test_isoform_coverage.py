"""Junction-core window coverage geometry.

Regression home for ``DESIGNER_ISSUE_isoform_coverage.md`` (2026-08-02).

``isoform_consensus`` used to credit an isoform for merely OVERLAPPING a
candidate window (``IsoformAwareness.isoforms_union``). Because that count is
also the sort key, the defect did not just fail to reject bad sites — it
PREFERRED them: a window straddling a region boundary touches more regions,
unions more isoforms, and outranks a window sitting cleanly inside the
canonical transcript.

The lab rule an isoform must satisfy to count as targeted:

    the junction-centred core of the binding site (roughly +/-16 nt around the
    ligation junction) must be present AND contiguous in that isoform's mRNA

because the ligase reads the junction, and distal-arm shortfall is tolerable
only down to a bound (an effective arm below ~16 nt will not hold at the lab's
37 C / 10 nM hybridization).

Coordinates in ``TestRealFabp1Straddler`` are the real CRC-panel FABP1 site 0
(``experiments/20260727_CRC_BRAF_round1/output_v2``).
"""
from __future__ import annotations

import pytest

from probe_designer.isoform_coverage import (
    PRODUCTIVE_BIOTYPES,
    WindowCoverage,
    credit_isoforms,
    is_productive,
    window_coverage,
)


def _iso(name, exons, *, biotype="protein_coding", strand=1):
    """Minimal isoform dict in the legacy Ensembl-Transcript shape."""
    return {
        "id": f"ENST_{name}",
        "display_name": name,
        "biotype": biotype,
        "strand": strand,
        "seq_region_name": "chrT",
        "start": min(s for s, _ in exons),
        "end": max(e for _, e in exons),
        "Exon": [{"start": s, "end": e} for s, e in exons],
    }


# ---------------------------------------------------------------------------
# The simple case: window entirely inside one exon
# ---------------------------------------------------------------------------

class TestWindowInsideOneExon:
    def test_full_coverage_both_sides(self):
        cov = window_coverage([(1000, 1039)], [(900, 1200)], strand=1)
        assert cov.lower_nt == 20
        assert cov.upper_nt == 20

    def test_full_coverage_credits_at_default_arm(self):
        cov = window_coverage([(1000, 1039)], [(900, 1200)], strand=1)
        assert cov.credited(16)

    def test_arms_are_full_length(self):
        cov = window_coverage([(1000, 1039)], [(900, 1200)], strand=1)
        assert (cov.arm_3prime_nt, cov.arm_5prime_nt) == (20, 20)

    def test_window_flush_with_exon_edges_is_still_full(self):
        cov = window_coverage([(1000, 1039)], [(1000, 1039)], strand=1)
        assert (cov.lower_nt, cov.upper_nt) == (20, 20)


# ---------------------------------------------------------------------------
# The defect: an exon boundary INSIDE the window
# ---------------------------------------------------------------------------

class TestRealFabp1Straddler:
    """FABP1 site 0 of the CRC panel, verbatim.

    Window chr2:88,123,073-88,123,112 (minus strand), built on FABP1-202 where
    it sits inside one exon. FABP1-201's exon ends at 88,123,104, so only 32 of
    the 40 nt are contiguous in that mRNA — BLAST measured 32/40, query
    positions 9-40. The union rule credited FABP1-201 anyway.
    """

    WINDOW = [(88123073, 88123112)]
    FABP1_202 = [(88122982, 88123138)]                    # one exon, spans it all
    FABP1_201 = [(88122982, 88123104), (88124000, 88124100)]  # exon ends mid-window

    def test_host_isoform_is_fully_covered(self):
        cov = window_coverage(self.WINDOW, self.FABP1_202, strand=-1)
        assert (cov.lower_nt, cov.upper_nt) == (20, 20)

    def test_straddled_isoform_loses_eight_nt(self):
        cov = window_coverage(self.WINDOW, self.FABP1_201, strand=-1)
        # 88123105..88123112 (8 nt) are spliced out of FABP1-201.
        assert cov.lower_nt + cov.upper_nt == 32

    def test_shortfall_lands_on_the_distal_end_of_arm5(self):
        cov = window_coverage(self.WINDOW, self.FABP1_201, strand=-1)
        # Minus strand: the genomic-upper half is the target's 5' half, which
        # the probe reads with arm5. 8 nt lost there leaves arm5 at 12.
        assert cov.arm_3prime_nt == 20
        assert cov.arm_5prime_nt == 12

    def test_straddled_isoform_is_not_credited(self):
        cov = window_coverage(self.WINDOW, self.FABP1_201, strand=-1)
        assert not cov.credited(16)

    def test_straddled_isoform_would_pass_a_lenient_bound(self):
        # Locks the knob's meaning: 12 nt survives a 12-nt bound, not a 16-nt one.
        cov = window_coverage(self.WINDOW, self.FABP1_201, strand=-1)
        assert cov.credited(12)


class TestPlusStrandArmMapping:
    """Same geometry, plus strand: the two arms swap sides."""

    WINDOW = [(1000, 1039)]
    CLIPPED = [(900, 1031)]  # exon ends 8 nt before the window's top

    def test_minus_strand_shortfall_is_arm5(self):
        cov = window_coverage(self.WINDOW, self.CLIPPED, strand=-1)
        assert (cov.arm_3prime_nt, cov.arm_5prime_nt) == (20, 12)

    def test_plus_strand_shortfall_is_arm3(self):
        cov = window_coverage(self.WINDOW, self.CLIPPED, strand=1)
        assert (cov.arm_3prime_nt, cov.arm_5prime_nt) == (12, 20)

    def test_side_lengths_are_strand_independent(self):
        plus = window_coverage(self.WINDOW, self.CLIPPED, strand=1)
        minus = window_coverage(self.WINDOW, self.CLIPPED, strand=-1)
        assert (plus.lower_nt, plus.upper_nt) == (minus.lower_nt, minus.upper_nt)


class TestJunctionItselfMissing:
    """If the ligation junction is not paired, the isoform is not bound at all."""

    def test_intron_over_the_junction_gives_zero(self):
        # Window 1000..1039; junction between 1019 and 1020. The isoform has an
        # intron at 1015..1024, so neither side reaches the junction.
        cov = window_coverage(
            [(1000, 1039)], [(900, 1014), (1025, 1100)], strand=1
        )
        assert not cov.ligatable
        assert (cov.lower_nt, cov.upper_nt) == (0, 0)

    def test_isoform_absent_from_the_locus_gives_zero(self):
        cov = window_coverage([(1000, 1039)], [(5000, 6000)], strand=1)
        assert not cov.ligatable
        assert (cov.lower_nt, cov.upper_nt) == (0, 0)

    def test_run_not_containing_the_junction_is_not_counted(self):
        # Present at 1000..1007 and 1017..1039; the junction (1019/1020) sits in
        # the SECOND run, so only that run counts: 3 nt below, 20 above.
        cov = window_coverage(
            [(1000, 1039)], [(900, 1007), (1017, 1100)], strand=1
        )
        assert (cov.lower_nt, cov.upper_nt) == (3, 20)
        assert not cov.credited(16)


# ---------------------------------------------------------------------------
# Splice junctions inside the window (contiguity, not mere presence)
# ---------------------------------------------------------------------------

class TestSplicedWindow:
    """A window that spans a splice junction of its host.

    ``find_target_region`` only ever returns a multi-segment window when the
    host's own exon ends inside it, so the segments are contiguous in the HOST
    mRNA by construction. Another isoform binds the same probe only when it
    shares that exact junction.
    """

    WINDOW = [(1000, 1019), (2000, 2019)]  # 20 + 20 nt across one host intron

    def test_isoform_sharing_the_junction_is_fully_covered(self):
        cov = window_coverage(self.WINDOW, [(900, 1019), (2000, 2100)], strand=1)
        assert (cov.lower_nt, cov.upper_nt) == (20, 20)

    def test_isoform_with_a_different_acceptor_is_not_credited(self):
        # Same donor (1019) but the acceptor is 2050, so the probe's upper half
        # is absent from this mRNA entirely.
        cov = window_coverage(self.WINDOW, [(900, 1019), (2050, 2100)], strand=1)
        assert not cov.ligatable
        assert (cov.lower_nt, cov.upper_nt) == (0, 0)

    def test_retained_intron_isoform_is_not_credited(self):
        # One long exon covering both halves AND the intron between them. Every
        # base is PRESENT, but they are not ADJACENT in the mRNA — the intron
        # sits in between — so the probe cannot ligate on it. This is exactly
        # the case `isoforms_union` got wrong.
        cov = window_coverage(self.WINDOW, [(900, 2100)], strand=1)
        assert not cov.ligatable
        assert (cov.lower_nt, cov.upper_nt) == (0, 0)
        assert not cov.credited(16)

    def test_unligatable_is_rejected_even_at_a_zero_bound(self):
        # `credited` must not degenerate into "0 >= 0" when the junction is
        # unpaired; ligatability is a separate, prior condition.
        cov = window_coverage(self.WINDOW, [(900, 2100)], strand=1)
        assert not cov.credited(0)

    def test_genomically_adjacent_exons_stay_contiguous(self):
        # Degenerate zero-length intron: exon ends at 1019, next starts at 1020.
        # Contiguous in the mRNA either way.
        cov = window_coverage([(1000, 1039)], [(900, 1019), (1020, 1100)], strand=1)
        assert (cov.lower_nt, cov.upper_nt) == (20, 20)


# ---------------------------------------------------------------------------
# The bound itself
# ---------------------------------------------------------------------------

class TestMinArmBound:
    @pytest.mark.parametrize(
        "exon_end, expected_upper",
        [(1035, 16), (1034, 15)],
    )
    def test_upper_arm_length_tracks_the_exon_end(self, exon_end, expected_upper):
        cov = window_coverage([(1000, 1039)], [(900, exon_end)], strand=1)
        assert cov.upper_nt == expected_upper

    def test_sixteen_nt_arm_is_credited_at_the_default_bound(self):
        cov = window_coverage([(1000, 1039)], [(900, 1035)], strand=1)
        assert cov.credited(16)

    def test_fifteen_nt_arm_is_rejected_at_the_default_bound(self):
        cov = window_coverage([(1000, 1039)], [(900, 1034)], strand=1)
        assert not cov.credited(16)

    def test_zero_bound_accepts_anything_that_pairs_the_junction(self):
        cov = window_coverage([(1000, 1039)], [(1019, 1020)], strand=1)
        assert (cov.lower_nt, cov.upper_nt) == (1, 1)
        assert cov.credited(0)


# ---------------------------------------------------------------------------
# credit_isoforms — the set the strategy consumes
# ---------------------------------------------------------------------------

class TestCreditIsoforms:
    ISOFORMS = [
        _iso("G-201", [(1000, 1199), (1300, 1399)]),                    # host
        _iso("G-204", [(1000, 1100), (1150, 1199), (1300, 1399)]),      # boundary at 1100
        _iso("G-202", [(1000, 1399)], biotype="retained_intron"),
        _iso("G-203", [(1200, 1299)], biotype="processed_transcript"),
    ]

    def test_clean_window_credits_every_isoform_that_contains_it(self):
        credited = credit_isoforms([(1000, 1039)], self.ISOFORMS, strand=1, min_arm=16)
        assert set(credited) == {"G-201", "G-204", "G-202"}

    def test_straddling_window_drops_the_boundary_isoform(self):
        # 1080..1119 straddles G-204's exon end at 1100: 21 nt present, only
        # 1 nt of it above the junction.
        credited = credit_isoforms([(1080, 1119)], self.ISOFORMS, strand=1, min_arm=16)
        assert "G-204" not in credited
        assert set(credited) == {"G-201", "G-202"}

    def test_credited_count_never_exceeds_the_isoform_count(self):
        for start in range(1000, 1360, 40):
            credited = credit_isoforms(
                [(start, start + 39)], self.ISOFORMS, strand=1, min_arm=16
            )
            assert len(credited) <= len(self.ISOFORMS)

    def test_result_carries_per_isoform_arm_lengths(self):
        credited = credit_isoforms([(1000, 1039)], self.ISOFORMS, strand=1, min_arm=16)
        assert isinstance(credited["G-201"], WindowCoverage)
        assert credited["G-201"].arm_3prime_nt == 20

    def test_window_in_the_hosts_intron_credits_no_coding_isoform(self):
        # The GNLY / MZB1 shape: a site built on a minor transcript whose exon
        # is intronic in every protein-coding transcript.
        credited = credit_isoforms([(1230, 1269)], self.ISOFORMS, strand=1, min_arm=16)
        assert set(credited) == {"G-202", "G-203"}
        assert not any(
            iso["biotype"] == "protein_coding"
            for iso in self.ISOFORMS
            if iso["display_name"] in credited
        )


# ---------------------------------------------------------------------------
# Which biotypes count as "codes for a protein"
# ---------------------------------------------------------------------------

class TestProductiveBiotypes:
    """`protein_coding` alone is too narrow.

    Immune-receptor genes carry their own Ensembl biotypes — TRAC and TRBC1 are
    ``TR_C_gene``, IGHG1 is ``IG_C_gene`` — and never appear as
    ``protein_coding``. The 2026-08-02 bank audit surfaced this by flagging six
    perfectly good TRAC/TRBC1 sites as binding no coding transcript.
    """

    def test_protein_coding_is_productive(self):
        assert is_productive("protein_coding")

    @pytest.mark.parametrize(
        "biotype",
        ["TR_C_gene", "TR_V_gene", "TR_J_gene", "TR_D_gene",
         "IG_C_gene", "IG_V_gene", "IG_J_gene", "IG_D_gene"],
    )
    def test_immune_receptor_segments_are_productive(self, biotype):
        assert is_productive(biotype)

    @pytest.mark.parametrize(
        "biotype",
        ["retained_intron", "processed_transcript", "nonsense_mediated_decay",
         "lncRNA", "protein_coding_CDS_not_defined", "IG_C_pseudogene",
         "TR_V_pseudogene", "processed_pseudogene", ""],
    )
    def test_unproductive_biotypes_are_rejected(self, biotype):
        assert not is_productive(biotype)

    def test_pseudogenes_are_not_admitted_by_a_prefix_match(self):
        # The obvious implementation — "starts with IG_ or TR_" — would let
        # every immune pseudogene in. GENCODE has hundreds of them.
        assert not (PRODUCTIVE_BIOTYPES & {"IG_C_pseudogene", "TR_V_pseudogene"})
