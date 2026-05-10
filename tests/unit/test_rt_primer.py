"""Tests for probe_designer.rt_primer.design_rt_primer."""
from __future__ import annotations

import pytest

from probe_designer.rt_primer import design_rt_primer


def make_constant_accessor(motif: str):
    """Genome accessor stub: returns a repeated motif of any length asked.

    The point is that the function under test only cares about the
    nucleotide composition of the slice, not its absolute coordinates.
    """
    def acc(chrom, start, end):
        n = max(0, end - (start - 1))
        if n == 0 or not motif:
            return ""
        # Repeat the motif enough times then trim.
        repeats = (n // len(motif)) + 1
        return (motif * repeats)[:n]
    return acc


def test_designs_primer_within_tm_range_on_plus_strand():
    """A balanced 50/50 GC seed should yield a primer in 50-65 °C window."""
    # Constant 50% GC sequence — Tm scales with length.
    seed = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"  # 40 nt, 50% GC
    acc = make_constant_accessor(seed)
    out = design_rt_primer(
        chrom="chr1",
        mut_start=100, mut_end=100,
        strand=1, probe_arm5_len=20,
        genome_accessor=acc,
    )
    assert out["primer_length"] >= 15
    assert 50.0 <= out["tm"] <= 65.0
    assert out["primer_seq"]
    # +strand: primer is RC of fetched genomic.
    assert set(out["primer_seq"]).issubset(set("ACGT"))


def test_minus_strand_does_not_reverse_complement():
    """For strand=-1 the primer is the literal genomic slice (sense strand)."""
    seed = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    acc = make_constant_accessor(seed)
    out = design_rt_primer(
        chrom="chr1",
        mut_start=100, mut_end=100,
        strand=-1, probe_arm5_len=20,
        genome_accessor=acc,
    )
    assert out["primer_seq"] != ""
    # No RC step ⇒ primer is upper-cased substring of genomic
    assert out["primer_seq"] in seed.upper()


def test_failed_genome_lookup_returns_failure_record():
    """If genome accessor returns empty, primer fails gracefully."""
    out = design_rt_primer(
        chrom="chr99",
        mut_start=1, mut_end=1,
        strand=1, probe_arm5_len=20,
        genome_accessor=lambda c, s, e: "",
    )
    assert out["primer_seq"] == ""
    assert "FAILED" in out["notes"]


def test_strand_default_on_nan():
    """NaN strand value should default to +1 (no error)."""
    seed = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    acc = make_constant_accessor(seed)
    out = design_rt_primer(
        chrom="chr1",
        mut_start=100, mut_end=100,
        strand=float("nan"), probe_arm5_len=20,
        genome_accessor=acc,
    )
    assert out["strand"] == 1


def test_out_of_range_tm_records_warning_in_notes():
    """When no length yields Tm in range, fallback path picks closest and notes it."""
    # Force fallback by demanding a Tm window no real primer can hit
    seed = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"  # 50% GC
    acc = make_constant_accessor(seed)
    out = design_rt_primer(
        chrom="chr1", mut_start=100, mut_end=100,
        strand=1, probe_arm5_len=20,
        genome_accessor=acc,
        tm_min=90.0, tm_max=100.0,  # impossible window
    )
    assert out["primer_seq"] != ""
    assert out["tm"] < 90.0
    assert "below" in out["notes"]


# ---------------------------------------------------------------------------
# design_rt_primer_from_target — TCR entry point
# ---------------------------------------------------------------------------


from probe_designer.rt_primer import design_rt_primer_from_target  # noqa: E402


def test_from_target_basic_primer_in_tm_range():
    """A balanced TRB-like sequence yields a primer in [50, 65]."""
    # 50% GC repeating motif, plenty of length
    target = ("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"  # 40 nt BDS region
              "AAAATTTTAAAA"                              # 12 nt gap region
              "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC")  # 40 nt primer region
    out = design_rt_primer_from_target(
        target_seq=target, target_end=40, gap=12,
    )
    assert out["primer_seq"]
    assert 50.0 <= out["tm"] <= 65.0
    assert out["primer_length"] >= 15
    # primer is RC of target slice → should NOT equal the target slice
    slice_taken = target[40 + 12: 40 + 12 + out["primer_length"]]
    assert out["primer_seq"] != slice_taken


def test_from_target_fails_when_target_too_short():
    """If target_end + gap + min_primer_len > len(target), fail gracefully."""
    target = "ATGCATGCATGCATGCATGCATGCATGCATGC"  # 32 nt
    out = design_rt_primer_from_target(
        target_seq=target, target_end=20, gap=15,
    )
    # 20 + 15 + 15 = 50 > 32 ⇒ no primer fits
    assert out["primer_seq"] == ""
    assert "FAILED" in out["notes"]


def test_from_target_returns_target_offsets_in_record():
    """The result records the (start, end) the primer occupies on target_seq."""
    target = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC" * 4
    out = design_rt_primer_from_target(
        target_seq=target, target_end=40, gap=15,
    )
    assert out["primer_target_start"] == 55     # 40 + 15
    assert out["primer_target_end"] == 55 + out["primer_length"]
    assert out["strand"] == 1                   # always +1 for TCR target
