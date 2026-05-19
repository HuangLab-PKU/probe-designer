"""Tests for ``probe-design pool check <pool_id>``."""
from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("primer3")
pytest.importorskip("probe_book")


from probe_designer.cli.app import app


runner = CliRunner()


def _write_minimal_pool(tmp_path: Path) -> Path:
    """Build a minimal repo-root with one pool + matching probes registry.

    The two probes form a known cross-lig pair (GC repeat).
    """
    (tmp_path / "probes").mkdir()
    (tmp_path / "probes" / "registry.tsv").write_text(
        "probe_id\tcodebook\tchemistry\tprobe_arm5\tprobe_arm3\tsequence\ttarget\n"
        "SP369_AAAA\tSP369\tdRNA\tGCAGCAGCAGCAGCAGCAGC\tTGCTGCTGCTGCTGCTGCTG\t\tGENEX\n"
        "SP369_BBBB\tSP369\tdRNA\tCAGCAGCAGCAGCAGCAGCA\tGCTGCTGCTGCTGCTGCTGC\t\tGENEY\n",
        encoding="utf-8",
    )
    pool_dir = tmp_path / "pools" / "test_pool_v1"
    pool_dir.mkdir(parents=True)
    (pool_dir / "manifest.yaml").write_text(
        "pool_id: test_pool_v1\n"
        "pool_name: test\n"
        "codebook: SP369\n"
        "target_well_volume_uL: 10.0\n"
        "recipe:\n"
        "- order_id: ORD1\n"
        "  probe_id: SP369_AAAA\n"
        "  concentration_M: 1.0e-08\n"
        "- order_id: ORD1\n"
        "  probe_id: SP369_BBBB\n"
        "  concentration_M: 1.0e-08\n",
        encoding="utf-8",
    )
    return tmp_path


def test_pool_check_runs_and_writes_report(tmp_path):
    root = _write_minimal_pool(tmp_path)
    result = runner.invoke(app, [
        "pool", "check", "test_pool_v1",
        "--repo-root", str(root),
    ])
    assert result.exit_code == 0, result.output

    check_dir = root / "pools" / "test_pool_v1" / "pool_check"
    assert (check_dir / "tier1_hits.tsv").exists()
    assert (check_dir / "confirmed_hits.tsv").exists()
    # The two dimer probes should be in confirmed_hits.
    confirmed = (check_dir / "confirmed_hits.tsv").read_text(encoding="utf-8")
    assert "SP369_AAAA" in confirmed and "SP369_BBBB" in confirmed


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
    assert "not found" in result.output.lower() or "manifest" in result.output.lower()
