"""Tests for the biotype filter on DefaultIsoformProvider + CLI flag.

The filter lets callers say "I only want protein_coding transcripts to
participate in isoform_consensus", which avoids the well-known dilution
problem where Ensembl's retained_intron / NMD / processed_transcript
isoforms drag down isoform_overlap_num for genuine consensus regions.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from probe_designer.cli.app import app
from probe_designer.genome.isoform_provider import DefaultIsoformProvider


# --- Provider-level tests --------------------------------------------------


@pytest.fixture
def fake_ensembl_two_biotypes():
    """Two transcripts, one protein_coding + one retained_intron."""
    gene_lookup = {
        "id": "ENSG_X", "seq_region_name": "1",
        "start": 100, "end": 5000, "strand": 1,
    }
    expand = {
        "Transcript": [
            {"id": "ENST_PC", "display_name": "X-201",
             "biotype": "protein_coding",
             "seq_region_name": "1", "start": 100, "end": 5000, "strand": 1,
             "Exon": [{"start": 100, "end": 5000, "strand": 1}]},
            {"id": "ENST_RI", "display_name": "X-202",
             "biotype": "retained_intron",
             "seq_region_name": "1", "start": 100, "end": 4500, "strand": 1,
             "Exon": [{"start": 100, "end": 4500, "strand": 1}]},
        ]
    }
    sequence = [
        MagicMock(status_code=200, json=lambda: gene_lookup),
        MagicMock(status_code=200, json=lambda: expand),
    ]
    with patch("probe_designer.genome.ensembl_client.requests.get") as mocked:
        mocked.side_effect = sequence
        yield mocked


class TestDefaultProviderBiotypeFilter:
    def test_no_filter_returns_all_isoforms(self, fake_ensembl_two_biotypes):
        provider = DefaultIsoformProvider()
        isoforms = provider.get_isoforms("X", "human")
        assert {iso["biotype"] for iso in isoforms} == {"protein_coding", "retained_intron"}

    def test_filter_keeps_only_matching(self, fake_ensembl_two_biotypes):
        provider = DefaultIsoformProvider(keep_biotypes={"protein_coding"})
        isoforms = provider.get_isoforms("X", "human")
        assert len(isoforms) == 1
        assert isoforms[0]["biotype"] == "protein_coding"

    def test_filter_with_multiple_biotypes(self, fake_ensembl_two_biotypes):
        # Both biotypes kept — same as no filter, but explicit.
        provider = DefaultIsoformProvider(
            keep_biotypes={"protein_coding", "retained_intron"}
        )
        isoforms = provider.get_isoforms("X", "human")
        assert len(isoforms) == 2

    def test_filter_dropping_all_falls_back_to_unfiltered(self, fake_ensembl_two_biotypes):
        # No isoform matches "lncRNA", so we'd get []. Fallback returns the
        # full set so the pipeline can still attempt the run; the user sees
        # the warning in the gene_report errors.
        provider = DefaultIsoformProvider(keep_biotypes={"lncRNA"})
        isoforms = provider.get_isoforms("X", "human")
        assert len(isoforms) == 2

    def test_filter_with_empty_set_acts_like_no_filter(self, fake_ensembl_two_biotypes):
        # Empty set should not silently drop everything; treat as "no filter".
        provider = DefaultIsoformProvider(keep_biotypes=set())
        isoforms = provider.get_isoforms("X", "human")
        assert len(isoforms) == 2


# --- CLI-level test --------------------------------------------------------


class TestCliIsoformBiotypeFlag:
    """Smoke test: `probe-design design --isoform-biotype protein_coding` is
    accepted and forwarded to DefaultIsoformProvider. We mock Pipeline so the
    test stays a unit test (no network, no BLAST).
    """

    def test_flag_constructs_provider_with_keep_biotypes(
        self, tmp_path: Path
    ):
        genes_file = tmp_path / "g.txt"
        genes_file.write_text("FOO\n")
        out_dir = tmp_path / "out"

        captured: dict = {}

        class _StubPipeline:
            def __init__(self, *args, **kw):
                captured["isoform_provider"] = kw.get("isoform_provider")

            def run(self, genes, **kw):
                # Build the smallest possible result with the right shape.
                from probe_designer.pipeline.pipeline import (
                    GeneResult, PipelineResult,
                )
                gr = GeneResult(gene=genes[0])
                return PipelineResult(
                    per_gene={genes[0]: gr},
                    config_snapshot={},
                    run_id="test",
                    output_dir=out_dir,
                )

        runner = CliRunner()
        with patch("probe_designer.cli.design.Pipeline", _StubPipeline):
            result = runner.invoke(
                app,
                [
                    "design",
                    "--genes-file", str(genes_file),
                    "--strategy", "isoform_consensus",
                    "--species", "human",
                    "--output", str(out_dir),
                    "--isoform-biotype", "protein_coding",
                    "--isoform-biotype", "protein_coding_LoF",
                ],
            )

        assert result.exit_code == 0, result.output
        provider = captured["isoform_provider"]
        assert provider is not None
        assert isinstance(provider, DefaultIsoformProvider)
        assert provider.keep_biotypes == {"protein_coding", "protein_coding_LoF"}

    def test_no_flag_provider_has_no_filter(self, tmp_path: Path):
        genes_file = tmp_path / "g.txt"
        genes_file.write_text("FOO\n")
        out_dir = tmp_path / "out"

        captured: dict = {}

        class _StubPipeline:
            def __init__(self, *args, **kw):
                captured["isoform_provider"] = kw.get("isoform_provider")

            def run(self, genes, **kw):
                from probe_designer.pipeline.pipeline import (
                    GeneResult, PipelineResult,
                )
                gr = GeneResult(gene=genes[0])
                return PipelineResult(
                    per_gene={genes[0]: gr},
                    config_snapshot={},
                    run_id="test",
                    output_dir=out_dir,
                )

        runner = CliRunner()
        with patch("probe_designer.cli.design.Pipeline", _StubPipeline):
            result = runner.invoke(
                app,
                [
                    "design",
                    "--genes-file", str(genes_file),
                    "--strategy", "isoform_consensus",
                    "--species", "human",
                    "--output", str(out_dir),
                ],
            )

        assert result.exit_code == 0, result.output
        provider = captured["isoform_provider"]
        # Either keep_biotypes is None / empty — both mean "no filter".
        assert not getattr(provider, "keep_biotypes", None)
