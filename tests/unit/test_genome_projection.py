"""Transcript-profile -> genome-coordinate projection."""
from __future__ import annotations

import numpy as np

from probe_designer.genome_projection import (
    pick_canonical_isoform,
    project_profile_to_genome,
    transcript_to_genome_positions,
    write_genome_bedgraph,
)

# Two exons, Ensembl 1-based inclusive.
EXONS = [{"start": 100, "end": 102}, {"start": 200, "end": 201}]  # 3 + 2 = 5 nt


class TestTranscriptToGenome:
    def test_plus_strand_ascending(self):
        # transcript 5'->3' = ascending genomic
        assert transcript_to_genome_positions(EXONS, 1) == [100, 101, 102, 200, 201]

    def test_minus_strand_descending(self):
        # mRNA 5' end = highest genomic coordinate
        assert transcript_to_genome_positions(EXONS, -1) == [201, 200, 102, 101, 100]


class TestProjectProfile:
    def test_plus_strand_mapping_0based_sorted(self):
        prof = np.array([50.0, 51.0, 52.0, 53.0, 54.0])
        pairs = project_profile_to_genome(prof, EXONS, 1)
        assert pairs == [(99, 50.0), (100, 51.0), (101, 52.0), (199, 53.0), (200, 54.0)]

    def test_minus_strand_sorted_ascending(self):
        prof = np.array([50.0, 51.0, 52.0, 53.0, 54.0])  # t-pos0 -> genomic 201
        pairs = project_profile_to_genome(prof, EXONS, -1)
        # sorted by genomic pos: genomic 100(0based 99)=t4=54 ... 201(200)=t0=50
        assert pairs == [(99, 54.0), (100, 53.0), (101, 52.0), (199, 51.0), (200, 50.0)]

    def test_nan_skipped(self):
        prof = np.array([50.0, np.nan, 52.0, np.nan, 54.0])
        pairs = project_profile_to_genome(prof, EXONS, 1)
        assert [p[0] for p in pairs] == [99, 101, 200]


class TestWriteGenomeBedgraph:
    def test_coalesces_contiguous_equal(self, tmp_path):
        pairs = [(99, 50.0), (100, 50.0), (101, 52.0)]
        path = write_genome_bedgraph(pairs, "chr7", tmp_path / "g.bedgraph")
        lines = path.read_text(encoding="utf-8").splitlines()
        assert lines[0].startswith("track type=bedGraph")
        assert lines[1] == "chr7\t99\t101\t50.00"   # contiguous equal coalesced
        assert lines[2] == "chr7\t101\t102\t52.00"

    def test_gap_not_coalesced(self, tmp_path):
        pairs = [(99, 50.0), (200, 50.0)]  # same value but non-contiguous (intron)
        path = write_genome_bedgraph(pairs, "chr7", tmp_path / "g.bedgraph")
        lines = path.read_text(encoding="utf-8").splitlines()
        assert len(lines) == 3  # track + two separate intervals


class TestPickCanonical:
    def test_prefers_canonical_flag(self):
        isos = [
            {"display_name": "t1", "Exon": [{"start": 1, "end": 10}]},
            {"display_name": "t2", "is_canonical": 1, "Exon": [{"start": 1, "end": 5}]},
        ]
        assert pick_canonical_isoform(isos)["display_name"] == "t2"

    def test_falls_back_to_longest(self):
        isos = [
            {"display_name": "short", "Exon": [{"start": 1, "end": 10}]},
            {"display_name": "long", "Exon": [{"start": 1, "end": 50}]},
        ]
        assert pick_canonical_isoform(isos)["display_name"] == "long"
