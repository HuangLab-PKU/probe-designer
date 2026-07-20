"""Smoke test for probe_designer.ext.tcr.landscape.plot_tm_landscape."""
from __future__ import annotations

from pathlib import Path

import pytest

from probe_designer.ext.tcr.landscape import plot_tm_landscape


@pytest.fixture
def fake_clone_sites():
    """Two clones, ~20 sites each, with all required Tm + CDR3 fields."""
    out = {}
    for clone in ("CloneA", "CloneB"):
        sites = []
        for pos in range(0, 20):
            sites.append({
                "ligation_point": pos + 20,
                "tm_5prime_mRNA": 55.0 + (pos % 5),
                "tm_3prime_mRNA": 56.0 + (pos % 4),
                "tm_5prime_cDNA": 60.0 + (pos % 5),
                "tm_3prime_cDNA": 61.0 + (pos % 4),
                "cdr3_start": 30,
                "cdr3_end": 60,
            })
        out[clone] = sites
    return out


def test_plot_tm_landscape_writes_pdf(tmp_path: Path, fake_clone_sites):
    out = tmp_path / "landscape.pdf"
    selected = {"CloneA": [fake_clone_sites["CloneA"][5]]}
    result = plot_tm_landscape(fake_clone_sites, selected, pdf_path=out)
    if result is None:
        pytest.skip("matplotlib unavailable")
    assert out.exists()
    assert out.stat().st_size > 1024


def test_plot_tm_landscape_handles_empty_clone(tmp_path: Path):
    """An empty input dict still creates the PDF (or skips cleanly)."""
    out = tmp_path / "empty.pdf"
    result = plot_tm_landscape({}, {}, pdf_path=out)
    # Either matplotlib unavailable (None) or PDF was created with no pages.
    # Both are acceptable — no exception.
    assert result is None or out.exists()
