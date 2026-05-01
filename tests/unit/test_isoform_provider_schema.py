"""Regression tests: isoform providers must return a schema that
search_strategies.IsoformAwareness can consume.

Background: in v0.2.1, both fetch_ensembl_isoforms and parse_gtf_for_gene
returned isoform dicts shaped as
``{"accession", "biotype", "start", "end", "exons": [...]}``,
but search_strategies.IsoformAwareness reads
``iso["seq_region_name"]``, ``iso["display_name"]``, ``iso["id"]``,
``iso["strand"]`` and ``iso["Exon"]`` (capital E). Calling
``probe-design design --strategy isoform_consensus`` therefore raised
``KeyError: 'seq_region_name'`` for every gene. These tests pin the
schema so the bug can't come back silently.
"""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from probe_designer.genome.ensembl_client import fetch_ensembl_isoforms
from probe_designer.genome.gtf_parser import parse_gtf_for_gene
from probe_designer.search_strategies import IsoformAwareness


# Fields the search strategy reads off each isoform dict (legacy shape from
# raw Ensembl Transcript records, kept by the unfiltered _db.get_isoform_info
# path that webapp + early CLI relied on).
LEGACY_ISO_FIELDS = {
    "seq_region_name", "display_name", "id", "strand", "Exon",
    "start", "end",
}

LEGACY_EXON_FIELDS = {"start", "end"}


# --- ensembl_client.fetch_ensembl_isoforms ----------------------------------


@pytest.fixture
def fake_ensembl_responses():
    """Mock the two Ensembl REST calls fetch_ensembl_isoforms makes:
    /lookup/symbol/<species>/<gene>  and  /lookup/id/<gene_id>?expand=1.
    """
    gene_lookup = {
        "id": "ENSG_TEST",
        "seq_region_name": "11",
        "start": 1000,
        "end": 5000,
        "strand": 1,
    }
    expand = {
        "Transcript": [
            {
                "id": "ENST_A",
                "display_name": "TEST-201",
                "biotype": "protein_coding",
                "seq_region_name": "11",
                "start": 1000,
                "end": 5000,
                "strand": 1,
                "Exon": [
                    {"start": 1000, "end": 1200, "strand": 1},
                    {"start": 4000, "end": 5000, "strand": 1},
                ],
            },
            {
                "id": "ENST_B",
                "display_name": "TEST-202",
                "biotype": "retained_intron",
                "seq_region_name": "11",
                "start": 1000,
                "end": 4500,
                "strand": 1,
                "Exon": [
                    {"start": 1000, "end": 4500, "strand": 1},
                ],
            },
        ]
    }
    sequence = [
        MagicMock(status_code=200, json=lambda: gene_lookup),
        MagicMock(status_code=200, json=lambda: expand),
    ]
    with patch("probe_designer.genome.ensembl_client.requests.get") as mocked:
        mocked.side_effect = sequence
        yield mocked


class TestEnsemblClientSchema:
    def test_each_isoform_has_legacy_fields(self, fake_ensembl_responses):
        data = fetch_ensembl_isoforms("TEST", "human")
        assert data["isoforms"], "expected non-empty isoforms"
        for iso in data["isoforms"]:
            missing = LEGACY_ISO_FIELDS - set(iso)
            assert not missing, f"missing legacy fields on isoform: {missing}"

    def test_each_exon_has_start_end(self, fake_ensembl_responses):
        data = fetch_ensembl_isoforms("TEST", "human")
        for iso in data["isoforms"]:
            for exon in iso["Exon"]:
                missing = LEGACY_EXON_FIELDS - set(exon)
                assert not missing, f"missing exon fields: {missing}"

    def test_isoform_awareness_can_consume(self, fake_ensembl_responses):
        """The integration contract: feeding fetch_ensembl_isoforms output
        directly to IsoformAwareness must not raise (this is exactly what
        DefaultIsoformProvider does inside the CLI Pipeline).
        """
        data = fetch_ensembl_isoforms("TEST", "human")
        # Passing dummy genome accessor — IsoformAwareness only calls it on
        # demand from get_genome_seq, which __init__ doesn't touch.
        IsoformAwareness(
            isoforms=data["isoforms"],
            genome_accessor=lambda *a, **kw: "",
        )

    def test_new_shape_fields_still_present(self, fake_ensembl_responses):
        """Backward compat: callers using the v0.2.1 shape (accession,
        lowercase exons) keep working.
        """
        data = fetch_ensembl_isoforms("TEST", "human")
        for iso in data["isoforms"]:
            assert "accession" in iso
            assert "exons" in iso
            assert "biotype" in iso


# --- gtf_parser.parse_gtf_for_gene -----------------------------------------


SYNTHETIC_GTF = """\
1\tHAVANA\tgene\t100\t1000\t.\t+\t.\tgene_name "ACTB"; gene_id "ENSG001";
1\tHAVANA\ttranscript\t100\t1000\t.\t+\t.\ttranscript_id "ENST001"; gene_name "ACTB"; transcript_biotype "protein_coding";
1\tHAVANA\texon\t100\t300\t.\t+\t.\ttranscript_id "ENST001"; gene_name "ACTB";
1\tHAVANA\texon\t500\t800\t.\t+\t.\ttranscript_id "ENST001"; gene_name "ACTB";
2\tHAVANA\ttranscript\t5000\t6000\t.\t-\t.\ttranscript_id "ENST002"; gene_name "GAPDH"; transcript_biotype "protein_coding";
2\tHAVANA\texon\t5000\t6000\t.\t-\t.\ttranscript_id "ENST002"; gene_name "GAPDH";
"""


@pytest.fixture
def gtf_file(tmp_path):
    p = tmp_path / "mini.gtf"
    p.write_text(SYNTHETIC_GTF)
    return p


class TestGtfParserSchema:
    def test_each_isoform_has_legacy_fields(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "ACTB")
        assert data["isoforms"]
        for iso in data["isoforms"]:
            missing = LEGACY_ISO_FIELDS - set(iso)
            assert not missing, f"missing legacy fields on isoform: {missing}"

    def test_each_exon_has_start_end(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "ACTB")
        for iso in data["isoforms"]:
            for exon in iso["Exon"]:
                missing = LEGACY_EXON_FIELDS - set(exon)
                assert not missing

    def test_isoform_awareness_can_consume(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "ACTB")
        IsoformAwareness(
            isoforms=data["isoforms"],
            genome_accessor=lambda *a, **kw: "",
        )

    def test_minus_strand_isoforms_keep_strand(self, gtf_file):
        data = parse_gtf_for_gene(gtf_file, "GAPDH")
        for iso in data["isoforms"]:
            assert iso["strand"] == -1

    def test_new_shape_fields_still_present(self, gtf_file):
        """Same backward compat as Ensembl test."""
        data = parse_gtf_for_gene(gtf_file, "ACTB")
        for iso in data["isoforms"]:
            assert "accession" in iso
            assert "exons" in iso
            assert "biotype" in iso
