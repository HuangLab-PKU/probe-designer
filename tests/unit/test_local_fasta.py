"""Tests for probe_designer.genome.local_fasta.

Skips when GRCh38.fa is not present locally so CI machines without the
genome reference still pass.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from probe_designer.genome.local_fasta import build_genome_accessor


REPO_ROOT = Path(__file__).resolve().parents[3]
GRCH38_FASTA = REPO_ROOT / "data" / "genome" / "GRCh38.fa"


def _require_genome():
    if not GRCH38_FASTA.exists():
        pytest.skip(f"GRCh38 FASTA not present at {GRCH38_FASTA}")


def test_round_trip_known_locus():
    """Fetch a known TP53 promoter slice and check length + alphabet."""
    _require_genome()
    acc = build_genome_accessor(GRCH38_FASTA)
    seq = acc("chr17", 7676000, 7676100)
    assert len(seq) == 101
    assert set(seq.upper()).issubset(set("ACGTN"))


def test_chr_prefix_normalization_on_unprefixed_chrom():
    """Calling with bare '17' must hit the same slice as 'chr17'."""
    _require_genome()
    acc = build_genome_accessor(GRCH38_FASTA)
    a = acc("17", 7676000, 7676100)
    b = acc("chr17", 7676000, 7676100)
    assert a == b
    assert len(a) == 101


def test_unknown_chromosome_returns_empty_string():
    """Lookup of a non-existent chromosome must not raise; returns ''."""
    _require_genome()
    acc = build_genome_accessor(GRCH38_FASTA)
    assert acc("chr99", 1, 100) == ""
