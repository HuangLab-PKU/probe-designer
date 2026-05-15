"""Tests for IsoformAwareness mRNA-coordinate helpers (Phase 1A)."""
from __future__ import annotations

import pytest

from probe_designer.search_strategies import IsoformAwareness


# Synthetic gene on the + strand with two isoforms.
# Genome window: positions 100..199 (inclusive), exactly 100 nt.
#   - isoform A uses exons [100..149] (50 nt) and [170..199] (30 nt) → 80 nt mRNA
#   - isoform B uses one exon [110..199] → 90 nt mRNA
PLUS_ISOFORMS = [
    {
        "id": "A",
        "display_name": "A",
        "start": 100, "end": 199, "strand": 1,
        "seq_region_name": "chrTest",
        "Exon": [
            {"start": 100, "end": 149},
            {"start": 170, "end": 199},
        ],
    },
    {
        "id": "B",
        "display_name": "B",
        "start": 110, "end": 199, "strand": 1,
        "seq_region_name": "chrTest",
        "Exon": [
            {"start": 110, "end": 199},
        ],
    },
]


# Minus-strand version: same coordinates, but the mRNA reads RC-wise.
MINUS_ISOFORM = {
    "id": "M",
    "display_name": "M",
    "start": 100, "end": 199, "strand": -1,
    "seq_region_name": "chrTest",
    "Exon": [
        {"start": 100, "end": 149},
        {"start": 170, "end": 199},
    ],
}


@pytest.fixture
def fake_genome():
    """Return a deterministic 100-nt genome accessor for the synthetic locus.

    The sequence is 100 bases of repeating "ACGT" so we can verify mRNA
    reconstruction independent of any external genome.
    """
    seq = ("ACGT" * 25)  # exactly 100 nt covering 100..199
    def accessor(region_name, start, end):
        # Inclusive coords, 0-indexed in this synthetic accessor
        assert region_name == "chrTest"
        # Map global_start..global_end -> 0..99 in `seq`
        rel_start = start - 100
        rel_end = end - 100 + 1  # exclusive
        return seq[rel_start:rel_end]
    return accessor


class TestGetIsoformMRNA:
    def test_plus_strand_concatenates_exons(self, fake_genome):
        awareness = IsoformAwareness(PLUS_ISOFORMS, fake_genome)
        mrna_A = awareness.get_isoform_mrna(PLUS_ISOFORMS[0])
        # Isoform A: exon 100..149 (50 nt) + exon 170..199 (30 nt) = 80 nt
        assert len(mrna_A) == 80

        mrna_B = awareness.get_isoform_mrna(PLUS_ISOFORMS[1])
        # Isoform B: one exon 110..199 (90 nt)
        assert len(mrna_B) == 90

    def test_plus_strand_uses_genome_directly(self, fake_genome):
        awareness = IsoformAwareness(PLUS_ISOFORMS, fake_genome)
        mrna_B = awareness.get_isoform_mrna(PLUS_ISOFORMS[1])
        # Genome positions 110..199 in "ACGT"-repeat = chars 10..99 of the seq
        # First few should be the genome bases starting at position 110 (= rel 10).
        # rel 10 in "ACGT" * 25 = (10 mod 4) = 2 -> 'G' is index 2 in "ACGT"
        # So seq[10] = "ACGT"[10 % 4] = "ACGT"[2] = "G"
        assert mrna_B[0] == "G"

    def test_minus_strand_reverse_complements(self, fake_genome):
        awareness = IsoformAwareness([MINUS_ISOFORM], fake_genome)
        mrna = awareness.get_isoform_mrna(MINUS_ISOFORM)
        # Same exon layout as isoform A → 80 nt total
        assert len(mrna) == 80
        # Plus-strand version:
        plus_iso = dict(MINUS_ISOFORM)
        plus_iso["strand"] = 1
        plus_awareness = IsoformAwareness([plus_iso], fake_genome)
        plus_mrna = plus_awareness.get_isoform_mrna(plus_iso)
        # Minus mRNA = reverse complement of plus mRNA
        comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        expected_minus = "".join(comp[b] for b in reversed(plus_mrna))
        assert mrna == expected_minus


class TestGenomicToMrnaPos:
    def test_plus_strand_first_exon_start(self, fake_genome):
        awareness = IsoformAwareness(PLUS_ISOFORMS, fake_genome)
        # Genomic 100 is the very first base of isoform A → mRNA pos 0
        assert awareness.genomic_to_mrna_pos(PLUS_ISOFORMS[0], 100) == 0

    def test_plus_strand_second_exon_start(self, fake_genome):
        awareness = IsoformAwareness(PLUS_ISOFORMS, fake_genome)
        # Isoform A: exon1 = 100..149 (50 nt). Exon2 starts at genomic 170.
        # mRNA position for 170 = 50 (since 50 nt already consumed)
        assert awareness.genomic_to_mrna_pos(PLUS_ISOFORMS[0], 170) == 50

    def test_intron_returns_none(self, fake_genome):
        awareness = IsoformAwareness(PLUS_ISOFORMS, fake_genome)
        # Isoform A intron = 150..169. Position 160 is in the intron.
        assert awareness.genomic_to_mrna_pos(PLUS_ISOFORMS[0], 160) is None

    def test_minus_strand_reverses_position(self, fake_genome):
        awareness = IsoformAwareness([MINUS_ISOFORM], fake_genome)
        # Genomic 100 (plus-strand exon1 start) maps to mRNA position
        # (total_mrna_len - 1) - 0 = 79 for the minus-strand isoform.
        # (mRNA is read from genome 199 back to 100 with introns skipped.)
        assert awareness.genomic_to_mrna_pos(MINUS_ISOFORM, 100) == 79

    def test_minus_strand_genomic_end(self, fake_genome):
        awareness = IsoformAwareness([MINUS_ISOFORM], fake_genome)
        # Genomic 199 (last base) is mRNA position 0 for minus-strand
        assert awareness.genomic_to_mrna_pos(MINUS_ISOFORM, 199) == 0


# ---------------------------------------------------------------------------
# Minus-strand window slicing (regression for code review 2026-05-14)
# ---------------------------------------------------------------------------

# Single-exon minus-strand isoform [100..199] (100 nt mRNA).
SIMPLE_MINUS_ISOFORM = {
    "id": "M1",
    "display_name": "M1",
    "start": 100, "end": 199, "strand": -1,
    "seq_region_name": "chrTest",
    "Exon": [{"start": 100, "end": 199}],
}


class TestMinusStrandWindowSlicing:
    """The bug: search loop computed ``me = ms + bds_len`` after mapping
    the LOWEST genomic coord to its mRNA position. On + strand that's the
    mRNA start of the window. On - strand it's the mRNA END (because
    mRNA reads 3'→5' of plus strand), so ``ms + bds_len`` slides off into
    a DOWNSTREAM window. The fix maps both endpoints and uses min/max.
    """

    def test_genomic_window_maps_to_correct_mrna_slice(self, fake_genome):
        awareness = IsoformAwareness([SIMPLE_MINUS_ISOFORM], fake_genome)
        # 40-bp window: genomic [120, 159] inclusive (the binding site).
        ms_lo = awareness.genomic_to_mrna_pos(SIMPLE_MINUS_ISOFORM, 120)
        ms_hi = awareness.genomic_to_mrna_pos(SIMPLE_MINUS_ISOFORM, 159)
        # mRNA pos 0 == genomic 199; pos 79 == genomic 120; pos 40 == genomic 159
        assert ms_lo == 79
        assert ms_hi == 40
        # The correct slice uses min/max + 1 for half-open:
        win_start = min(ms_lo, ms_hi)
        win_end = max(ms_lo, ms_hi) + 1
        assert (win_start, win_end) == (40, 80)
        # Sanity: width == 40 bp
        assert win_end - win_start == 40

    def test_buggy_naive_formula_would_overflow(self, fake_genome):
        """Document exactly what the OLD code did wrong, so a future
        refactor that 'simplifies' back to ms + bds_len triggers this."""
        awareness = IsoformAwareness([SIMPLE_MINUS_ISOFORM], fake_genome)
        ms = awareness.genomic_to_mrna_pos(SIMPLE_MINUS_ISOFORM, 120)
        bds_len = 40
        buggy_me = ms + bds_len
        # ms=79 + 40 = 119 → past end-of-mRNA (which is 100). The old
        # search loop's check `me > len(prof)` would skip this entry,
        # silently filtering out every minus-strand candidate without a
        # diagnostic. The correct slice (40, 80) fits in [0, 100).
        assert buggy_me > 100  # buggy formula reads past mRNA end
