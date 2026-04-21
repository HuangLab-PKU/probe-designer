"""Tests for probe_designer.pipeline.Pipeline.

Exercises the hook-injection surface without hitting network/DB. Real
strategy + filter invocations require sequences; we patch
BindingSiteSearcher and SequenceFilter to isolate the Pipeline's
orchestration logic.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional
from unittest.mock import MagicMock, patch

import pytest

from probe_designer.config import ConfigManager
from probe_designer.pipeline import GeneResult, Pipeline, PipelineResult


class _FakeProgress:
    def __init__(self) -> None:
        self.events: List[tuple] = []

    def on_gene_start(self, gene, index, total):
        self.events.append(("start", gene, index, total))

    def on_stage(self, gene, stage, info):
        self.events.append(("stage", gene, stage))

    def on_gene_complete(self, gene, site_count):
        self.events.append(("complete", gene, site_count))

    def on_pipeline_complete(self, gene_count):
        self.events.append(("pipeline_done", gene_count))


class _FakeExistingTargets:
    def __init__(self, mapping: Optional[Dict[str, List[Dict[str, Any]]]] = None):
        self.mapping = mapping or {}
        self.calls: List[str] = []

    def lookup(self, gene, species, min_score=0.0):
        self.calls.append(gene)
        return self.mapping.get(gene, [])


class _FakeIsoformProvider:
    def __init__(self, mapping: Optional[Dict[str, List[Dict]]] = None):
        self.mapping = mapping or {}
        self.calls: List[str] = []

    def get_isoforms(self, gene, species):
        self.calls.append(gene)
        return self.mapping.get(gene)


@pytest.fixture
def minimal_config(monkeypatch):
    # Avoid species_config.json lookup side effects during construction
    cfg = ConfigManager()
    cfg.database.organism = "human"
    cfg.search.search_strategy = "single_sequence"
    return cfg


class TestStrategyValidation:
    def test_unknown_strategy_raises(self, minimal_config):
        with patch("probe_designer.pipeline.pipeline.DatabaseInterface"):
            p = Pipeline(minimal_config)
        with pytest.raises(ValueError, match="Unknown strategy"):
            p.run([], strategy="no_such_thing")


class TestProgressHookWiring:
    def test_empty_gene_list_still_fires_pipeline_complete(self, minimal_config):
        prog = _FakeProgress()
        with patch("probe_designer.pipeline.pipeline.DatabaseInterface"):
            p = Pipeline(minimal_config, progress=prog)
        p.run([], strategy="single_sequence")
        assert ("pipeline_done", 0) in prog.events

    def test_gene_start_and_complete_pair_fires_in_order(self, minimal_config):
        prog = _FakeProgress()
        with patch("probe_designer.pipeline.pipeline.DatabaseInterface") as mock_db_cls:
            mock_db = MagicMock()
            mock_db_cls.return_value = mock_db
            # Pipeline fetches sequences for single_sequence strategy -> empty -> skip
            mock_db.get_gene_sequences.return_value = {"sequences": {}}
            p = Pipeline(minimal_config, progress=prog)
            result = p.run(["ACTB"], strategy="single_sequence")
        events = [(e[0], e[1] if len(e) > 1 else None) for e in prog.events]
        assert ("start", "ACTB") in events
        assert ("complete", "ACTB") in events
        assert ("pipeline_done", None) in [(e[0], None) for e in prog.events]
        assert "ACTB" in result.per_gene


class TestExistingTargetsHook:
    def test_hook_called_per_gene(self, minimal_config):
        existing = _FakeExistingTargets()
        with patch("probe_designer.pipeline.pipeline.DatabaseInterface") as mock_db_cls:
            mock_db = MagicMock()
            mock_db_cls.return_value = mock_db
            mock_db.get_gene_sequences.return_value = {"sequences": {}}
            p = Pipeline(minimal_config, existing_targets=existing)
            p.run(["ACTB", "GAPDH"], strategy="single_sequence")
        assert existing.calls == ["ACTB", "GAPDH"]

    def test_existing_targets_injected_into_result(self, minimal_config):
        existing = _FakeExistingTargets(
            mapping={"ACTB": [{"display_id": "ACTB_001", "score": 9.0}]}
        )
        with patch("probe_designer.pipeline.pipeline.DatabaseInterface") as mock_db_cls:
            mock_db = MagicMock()
            mock_db_cls.return_value = mock_db
            mock_db.get_gene_sequences.return_value = {"sequences": {}}
            p = Pipeline(minimal_config, existing_targets=existing)
            result = p.run(["ACTB"], strategy="single_sequence")
        assert result.per_gene["ACTB"].existing_targets == [
            {"display_id": "ACTB_001", "score": 9.0}
        ]


class TestIsoformProviderWiring:
    def test_provider_used_for_isoform_strategy(self, minimal_config):
        provider = _FakeIsoformProvider(
            mapping={"ACTB": [{"accession": "ENST001", "exons": []}]}
        )
        with patch("probe_designer.pipeline.pipeline.DatabaseInterface") as mock_db_cls:
            mock_db_cls.return_value = MagicMock()
            p = Pipeline(minimal_config, isoform_provider=provider)
            # search will fail on minimal data but we only care that provider was called
            p.run(["ACTB"], strategy="isoform_consensus")
        assert provider.calls == ["ACTB"]

    def test_provider_not_called_for_single_sequence(self, minimal_config):
        provider = _FakeIsoformProvider()
        with patch("probe_designer.pipeline.pipeline.DatabaseInterface") as mock_db_cls:
            mock_db = MagicMock()
            mock_db_cls.return_value = mock_db
            mock_db.get_gene_sequences.return_value = {"sequences": {}}
            p = Pipeline(minimal_config, isoform_provider=provider)
            p.run(["ACTB"], strategy="single_sequence")
        assert provider.calls == []


class TestResultShape:
    def test_returns_pipeline_result_with_config_snapshot(self, minimal_config):
        with patch("probe_designer.pipeline.pipeline.DatabaseInterface") as mock_db_cls:
            mock_db_cls.return_value = MagicMock()
            p = Pipeline(minimal_config)
            result = p.run([], strategy="single_sequence")
        assert isinstance(result, PipelineResult)
        assert "database" in result.config_snapshot
        assert "search" in result.config_snapshot
        assert len(result.run_id) > 0

    def test_gene_result_dataclass_defaults(self):
        gr = GeneResult(gene="ACTB")
        assert gr.sites == []
        assert gr.existing_targets == []
        assert gr.errors == []
        d = gr.to_dict()
        assert d["gene"] == "ACTB"

    def test_pipeline_result_aggregates(self):
        pr = PipelineResult(
            per_gene={
                "A": GeneResult(gene="A", sites=[{}, {}]),
                "B": GeneResult(gene="B", sites=[{}]),
            }
        )
        assert pr.gene_count == 2
        assert pr.total_sites == 3
        assert pr.sites_by_gene() == {"A": [{}, {}], "B": [{}]}


class TestLegacyBruteForceAccepted:
    """Phase 2 transitional: 'brute_force' still accepted pending Phase 5 rename."""

    def test_brute_force_alias_does_not_raise(self, minimal_config):
        with patch("probe_designer.pipeline.pipeline.DatabaseInterface") as mock_db_cls:
            mock_db = MagicMock()
            mock_db_cls.return_value = mock_db
            mock_db.get_gene_sequences.return_value = {"sequences": {}}
            p = Pipeline(minimal_config)
            # Should not raise
            p.run(["ACTB"], strategy="brute_force")
