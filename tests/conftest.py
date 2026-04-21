"""Shared fixtures for designer tests."""
from __future__ import annotations

import pytest

from probe_designer.config import FilterConfig, BlastConfig


@pytest.fixture
def default_filter_config() -> FilterConfig:
    """FilterConfig with pipeline defaults; tests can override individual fields."""
    return FilterConfig()


@pytest.fixture
def default_blast_config() -> BlastConfig:
    return BlastConfig()


@pytest.fixture
def sample_binding_site():
    """Representative binding site dict produced by the search+filter pipeline."""
    return {
        "gene_name": "ACTB",
        "sequence": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "target_sequence": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "arm_3prime": "ACGTACGTACGTACGTACGT",
        "arm_5prime": "ACGTACGTACGTACGTACGT",
        "st": 100,
        "en": 140,
        "strand": "+",
        "length": 40,
        "g_content": 0.25,
        "tm": 58.0,
        "tm_3prime": 56.0,
        "tm_5prime": 57.0,
        "tm_diff": 1.0,
        "free_energy": -2.5,
        "isoform_overlap_num": 1,
        "blast_alignments": [],
    }
