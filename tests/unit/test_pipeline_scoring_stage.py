"""`Pipeline._run_single` must actually reach its scoring stage in a test.

Every test in ``test_pipeline.py`` mocks ``get_gene_sequences`` to return
``{"sequences": {}}``, so ``_run_single`` returns at the "No reference sequence"
guard long before search / filter / score ever run. A full-green suite therefore
said nothing about the back half of the per-gene body — and that is exactly where
a `ThermoPolicy` call shipped without its import, a NameError on every real
`probe-design design` run that 634 passing tests did not notice.

These tests drive a gene all the way through to selection with the DB and BLAST
stubbed, so the scoring stage is executed rather than skipped.
"""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from probe_designer.config import ConfigManager
from probe_designer.pipeline import Pipeline

# 240 nt of mixed-composition sequence: long enough for the single_sequence
# strategy to slide out several candidate windows.
TARGET = ("GATGATGATGATGATGATGCTAGCTAGCTAGCTAGCTAGCATGCATGCATGCATGCATGC"
          "TTGACCAGTGCATTCGAAGTGGCCGGCCGGAAGGCCGGCCATATATCGCGATATATGCTA"
          "GCATCATCATCATCATCATCGGATCCGGATCCGGATCCAAGCTTAAGCTTAAGCTTGGCC"
          "ACGTACGTACGTACGTACGTCCGGAATTCCGGAATTCCGGTTAACCGGTTAACCGGTTAA")


@pytest.fixture
def pipeline_with_sequence():
    cfg = ConfigManager()
    cfg.database.organism = "human"
    cfg.search.search_strategy = "single_sequence"
    cfg.search.binding_site_length = 40
    with patch("probe_designer.pipeline.pipeline.DatabaseInterface") as db_cls:
        db = db_cls.return_value
        # Shape matters: _get_sequences_for_gene keys by GENE and reads
        # seq_data["sequence"]. Getting this wrong is why the existing tests
        # all bail at the "No reference sequence" guard.
        db.get_gene_sequences.return_value = {
            "sequences": {
                "TESTGENE": {
                    "sequence": TARGET,
                    # single_sequence builds a genomic_context from these four
                    # and raises KeyError without them.
                    "seq_region_name": "1",
                    "start": 1_000_000,
                    "end": 1_000_000 + len(TARGET) - 1,
                    "strand": 1,
                },
            },
        }
        yield Pipeline(cfg)


class TestScoringStageIsReached:
    def test_run_single_scores_and_selects_sites(self, pipeline_with_sequence):
        """The regression: this path NameError'd on a missing ThermoPolicy import."""
        result = pipeline_with_sequence.run(
            ["TESTGENE"], strategy="single_sequence", skip_blast=True,
        )
        gene_result = result.per_gene["TESTGENE"]
        assert not gene_result.errors, f"pipeline reported errors: {gene_result.errors}"
        assert gene_result.sites, "no sites survived to selection"
        for site in gene_result.sites:
            assert isinstance(site["score"], float)
            assert "peak_rank" in site

    def test_scores_respond_to_the_configured_policy(self, pipeline_with_sequence):
        """Proves the config actually reaches the scorer, not just that it ran."""
        pipeline = pipeline_with_sequence
        baseline = pipeline.run(
            ["TESTGENE"], strategy="single_sequence", skip_blast=True,
        ).per_gene["TESTGENE"].sites

        # Move the Tm target far from where these arms sit; proximity must drop.
        pipeline.config.filter.min_tm = 95.0
        pipeline.config.filter.max_tm = 115.0
        shifted = pipeline.run(
            ["TESTGENE"], strategy="single_sequence", skip_blast=True,
        ).per_gene["TESTGENE"].sites

        assert baseline and shifted
        mean_base = sum(s["score"] for s in baseline) / len(baseline)
        mean_shifted = sum(s["score"] for s in shifted) / len(shifted)
        assert mean_shifted < mean_base, (
            "moving the Tm target 45 C away left scores unchanged — the "
            f"config is not reaching the scorer ({mean_base} vs {mean_shifted})"
        )

