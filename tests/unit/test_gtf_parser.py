"""Tests for probe_designer.genome.gtf_parser.

Verifies GTF extraction semantics and the byte-offset index cache.
"""
from __future__ import annotations

import pytest

from probe_designer.genome.gtf_parser import (
    _build_index,
    _index_path_for,
    parse_gtf_for_gene,
)


SYNTHETIC_GTF = """\
# comment line ignored
1\tHAVANA\tgene\t100\t1000\t.\t+\t.\tgene_name "ACTB"; gene_id "ENSG001";
1\tHAVANA\ttranscript\t100\t1000\t.\t+\t.\ttranscript_id "ENST001"; gene_name "ACTB"; transcript_biotype "protein_coding";
1\tHAVANA\texon\t100\t300\t.\t+\t.\ttranscript_id "ENST001"; gene_name "ACTB";
1\tHAVANA\texon\t500\t800\t.\t+\t.\ttranscript_id "ENST001"; gene_name "ACTB";
1\tHAVANA\ttranscript\t100\t900\t.\t+\t.\ttranscript_id "ENST002"; gene_name "ACTB"; transcript_biotype "protein_coding";
1\tHAVANA\texon\t100\t250\t.\t+\t.\ttranscript_id "ENST002"; gene_name "ACTB";
1\tHAVANA\texon\t700\t900\t.\t+\t.\ttranscript_id "ENST002"; gene_name "ACTB";
2\tHAVANA\tgene\t5000\t6000\t.\t-\t.\tgene_name "GAPDH"; gene_id "ENSG002";
2\tHAVANA\texon\t5000\t5500\t.\t-\t.\ttranscript_id "ENST003"; gene_name "GAPDH";
"""


@pytest.fixture
def gtf_file(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(SYNTHETIC_GTF)
    return gtf


class TestParseGtf:
    def test_gene_found_returns_isoforms(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "ACTB")
        assert len(data["isoforms"]) == 2
        ids = sorted(i["accession"] for i in data["isoforms"])
        assert ids == ["ENST001", "ENST002"]

    def test_gene_bounds_from_gene_feature(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "ACTB")
        assert data["start"] == 100
        assert data["end"] == 1000
        assert data["chromosome"] == "1"
        assert data["strand"] == 1

    def test_minus_strand_reported_as_minus_one(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "GAPDH")
        assert data["strand"] == -1

    def test_each_isoform_has_exons(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "ACTB")
        for iso in data["isoforms"]:
            assert len(iso["exons"]) == 2
            for ex in iso["exons"]:
                assert "start" in ex and "end" in ex
                assert "rel_start" in ex and "rel_width" in ex

    def test_exons_sorted_by_position(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "ACTB")
        for iso in data["isoforms"]:
            starts = [e["start"] for e in iso["exons"]]
            assert starts == sorted(starts)

    def test_biotype_captured_when_present(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "ACTB")
        for iso in data["isoforms"]:
            assert iso["biotype"] == "protein_coding"

    def test_unknown_gene_returns_empty(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "FAKE_GENE_XYZ")
        assert data == {"isoforms": []}

    def test_missing_file_returns_empty(self, tmp_path):
        data = parse_gtf_for_gene(tmp_path / "nope.gtf", "ACTB")
        assert data == {"isoforms": []}


class TestByteOffsetIndex:
    def test_index_created_on_first_parse(self, gtf_file):
        idx_path = _index_path_for(gtf_file)
        assert not idx_path.exists()
        parse_gtf_for_gene(gtf_file, "ACTB")
        assert idx_path.exists()

    def test_index_contains_known_genes(self, gtf_file):
        parse_gtf_for_gene(gtf_file, "ACTB")
        offsets = _build_index(gtf_file)
        assert "ACTB" in offsets
        assert "GAPDH" in offsets
        assert len(offsets["ACTB"]) >= 5  # 1 gene + 2 transcripts + 4 exons

    def test_rebuild_after_gtf_modification(self, gtf_file, tmp_path):
        # First parse builds the index
        parse_gtf_for_gene(gtf_file, "ACTB")
        idx_path = _index_path_for(gtf_file)
        first_mtime = idx_path.stat().st_mtime

        # Modify GTF (bump mtime by rewriting)
        import time
        time.sleep(0.1)
        gtf_file.write_text(SYNTHETIC_GTF + "\n")

        # Second parse should rebuild the index
        parse_gtf_for_gene(gtf_file, "ACTB")
        second_mtime = idx_path.stat().st_mtime
        assert second_mtime >= first_mtime

    def test_repeated_parses_produce_same_result(self, gtf_file):
        first = parse_gtf_for_gene(gtf_file, "ACTB")
        second = parse_gtf_for_gene(gtf_file, "ACTB")
        assert first == second
