"""Tests for ``probe_designer.ext.pool.loader``."""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("probe_book")


from probe_designer.ext.pool.loader import load_pool_as_probes_for_screen
from probe_designer.qc.cross_ligation import ProbeForScreen


def _write_minimal_pool(tmp_path: Path, *, pool_id: str = "test_pool_v1") -> Path:
    """Build a tiny pool + matching probes registry. Returns the repo_root."""
    (tmp_path / "probes").mkdir()
    (tmp_path / "probes" / "registry.tsv").write_text(
        "probe_id\tcodebook\tchemistry\tprobe_arm5\tprobe_arm3\tsequence\ttarget\n"
        "SP369_AAAA\tSP369\tdRNA\tGCATAGCAGCAGCAGCATAG\tTGTGTGTGCACGCACGCATG\t"
        "GCATAGCAGCAGCAGCATAGTCCCTACACGACGCTCTTCCGTGTGTGTGCACGCACGCATG\tGENEA\n"
        "SP369_BBBB\tSP369\tdRNA\tAAGAAGAAGAAGAAGAAGAA\tTTCTTCTTCTTCTTCTTCTT\t"
        "AAGAAGAAGAAGAAGAAGAATCCCTACACGACGCTCTTCCGTTCTTCTTCTTCTTCTTCTT\tGENEB\n",
        encoding="utf-8",
    )
    pool_dir = tmp_path / "pools" / pool_id
    pool_dir.mkdir(parents=True)
    (pool_dir / "manifest.yaml").write_text(
        f"pool_id: {pool_id}\n"
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


def test_loads_pool_and_returns_probes_for_screen(tmp_path):
    root = _write_minimal_pool(tmp_path)
    probes = load_pool_as_probes_for_screen("test_pool_v1", root)
    assert len(probes) == 2
    assert all(isinstance(p, ProbeForScreen) for p in probes)
    ids = {p.probe_id for p in probes}
    assert ids == {"SP369_AAAA", "SP369_BBBB"}
    a = next(p for p in probes if p.probe_id == "SP369_AAAA")
    assert a.chemistry == "dRNA"
    assert a.probe_arm5 == "GCATAGCAGCAGCAGCATAG"
    assert "GCATAGCAGCAGCAGCATAG" in a.sequence.upper()
    assert a.target == "GENEA"


def test_loads_design_pool_without_order_ids(tmp_path):
    """A DESIGN pool (kind: design, recipe entries with no order_id) must still
    resolve into ProbeForScreen so cross-ligation / 'can I add this probe'
    validation works BEFORE the probes are ordered."""
    (tmp_path / "probes").mkdir()
    (tmp_path / "probes" / "registry.tsv").write_text(
        "probe_id\tcodebook\tchemistry\tprobe_arm5\tprobe_arm3\tsequence\ttarget\n"
        "SP369_AAAA\tSP369\tdRNA\tGCATAGCAGCAGCAGCATAG\tTGTGTGTGCACGCACGCATG\t"
        "GCATAGCAGCAGCAGCATAGTCCCTACACGACGCTCTTCCGTGTGTGTGCACGCACGCATG\tGENEA\n",
        encoding="utf-8",
    )
    pool_dir = tmp_path / "pools" / "cand_design_v1"
    pool_dir.mkdir(parents=True)
    (pool_dir / "manifest.yaml").write_text(
        "pool_id: cand_design_v1\n"
        "pool_name: candidate composition\n"
        "codebook: SP369\n"
        "target_well_volume_uL: 10.0\n"
        "kind: design\n"
        "recipe:\n"
        "- probe_id: SP369_AAAA\n"
        "  order_id: null\n"
        "  concentration_M: 0.0\n",
        encoding="utf-8",
    )
    probes = load_pool_as_probes_for_screen("cand_design_v1", tmp_path)
    assert len(probes) == 1
    assert probes[0].probe_id == "SP369_AAAA"
    assert probes[0].probe_arm5 == "GCATAGCAGCAGCAGCATAG"


def test_missing_pool_manifest_raises_clear_error(tmp_path):
    (tmp_path / "probes").mkdir()
    (tmp_path / "probes" / "registry.tsv").write_text(
        "probe_id\tcodebook\tchemistry\tprobe_arm5\tprobe_arm3\tsequence\ttarget\n",
        encoding="utf-8",
    )
    with pytest.raises(FileNotFoundError) as exc:
        load_pool_as_probes_for_screen("does_not_exist_v1", tmp_path)
    assert "manifest" in str(exc.value).lower() or "not found" in str(exc.value).lower()
