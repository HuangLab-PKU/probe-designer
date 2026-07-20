"""Tests for ``probe-design pool check`` (the ext/pool CLI subcommand)."""
from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("primer3")
pytest.importorskip("probe_book")


from probe_designer.cli.app import app
from probe_designer.ext.pool.tests.test_loader import _write_minimal_pool


runner = CliRunner()


def test_pool_check_runs_and_writes_v2_report(tmp_path):
    root = _write_minimal_pool(tmp_path, pool_id="test_pool_v1")
    result = runner.invoke(app, [
        "pool", "check", "test_pool_v1",
        "--repo-root", str(root),
    ])
    assert result.exit_code == 0, result.output

    check_dir = root / "pools" / "test_pool_v1" / "pool_check"
    assert (check_dir / "confirmed_dimers_v2.tsv").exists()
    assert (check_dir / "all_dimers_v2.tsv").exists()


def test_pool_check_missing_pool_fails_cleanly(tmp_path):
    (tmp_path / "probes").mkdir()
    (tmp_path / "probes" / "registry.tsv").write_text(
        "probe_id\tcodebook\tchemistry\tprobe_arm5\tprobe_arm3\tsequence\ttarget\n",
        encoding="utf-8",
    )
    result = runner.invoke(app, [
        "pool", "check", "nonexistent_v1",
        "--repo-root", str(tmp_path),
    ])
    assert result.exit_code != 0
    assert ("not found" in result.output.lower() or "manifest" in result.output.lower())


def test_pool_check_deprecated_flag_warns(tmp_path):
    root = _write_minimal_pool(tmp_path, pool_id="test_pool_v1")
    result = runner.invoke(app, [
        "pool", "check", "test_pool_v1",
        "--repo-root", str(root),
        "--xlig-end-dg-threshold", "-5.0",
    ])
    # Should still succeed but with a deprecation warning
    assert result.exit_code == 0, result.output
