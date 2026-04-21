"""Unified padlock probe design pipeline.

This package provides the central :class:`Pipeline` orchestration layer shared
by both the ``probe-design`` CLI and the webapp. Callers inject hooks to
customize progress reporting, DB lookups, and isoform data sources without
coupling the designer library to any webapp-specific state.
"""
from probe_designer.pipeline.hooks import (
    ExistingTargetsHook,
    IsoformProvider,
    ProgressHook,
)
from probe_designer.pipeline.pipeline import Pipeline
from probe_designer.pipeline.result import GeneResult, PipelineResult


__all__ = [
    "Pipeline",
    "PipelineResult",
    "GeneResult",
    "ProgressHook",
    "ExistingTargetsHook",
    "IsoformProvider",
]
