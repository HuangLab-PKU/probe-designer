"""Tests for ``probe-design pool create-design`` and ``pool promote``."""
from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("probe_book")

from probe_book.pool import PoolKind, PoolManifest
from probe_designer.cli.app import app

runner = CliRunner()

REG_HEADER = "probe_id\tcodebook\tchemistry\tprobe_arm5\tprobe_arm3\tsequence\ttarget\n"


def _write_registry(root: Path, probe_ids):
    (root / "probes").mkdir(parents=True, exist_ok=True)
    rows = REG_HEADER
    for i, pid in enumerate(probe_ids):
        rows += f"{pid}\tSP369\tdRNA\tARM5{i:04d}AAAAAAAAAAAA\tARM3{i:04d}TTTTTTTTTTTT\tSEQ{i}\tGENE{i}\n"
    (root / "probes" / "registry.tsv").write_text(rows, encoding="utf-8")


def _write_order_inventory(root: Path, order_id: str, probe_ids):
    d = root / "orders" / order_id
    d.mkdir(parents=True, exist_ok=True)
    body = "probe_id\tquantity_ordered_M\tquantity_remaining_M\tstatus\n"
    for pid in probe_ids:
        body += f"{pid}\t1e-06\t1e-06\tin_stock\n"
    (d / "inventory.tsv").write_text(body, encoding="utf-8")


def _write_design_pool(root: Path, pool_id: str, probe_ids):
    pm = PoolManifest(pool_id=pool_id, pool_name="cand", codebook="SP369",
                      target_well_volume_uL=10.0, kind=PoolKind.DESIGN.value)
    for pid in probe_ids:
        pm.add_probe(probe_id=pid)
    pm.to_yaml(root / "pools" / pool_id / "manifest.yaml")


# ---- create-design ----

def test_create_design_writes_design_pool(tmp_path):
    _write_registry(tmp_path, ["SP369_p1", "SP369_p2"])
    result = runner.invoke(app, [
        "pool", "create-design", "TNBCmarker_cand_v1",
        "--probe", "SP369_p1", "--probe", "SP369_p2",
        "--repo-root", str(tmp_path),
    ])
    assert result.exit_code == 0, result.output
    pm = PoolManifest.from_yaml(tmp_path / "pools" / "TNBCmarker_cand_v1" / "manifest.yaml")
    assert pm.kind == PoolKind.DESIGN.value
    assert pm.probe_ids() == ["SP369_p1", "SP369_p2"]
    assert all(e.order_id is None for e in pm.recipe)


def test_create_design_reads_probes_file(tmp_path):
    _write_registry(tmp_path, ["SP369_p1", "SP369_p2"])
    pf = tmp_path / "probes.txt"
    pf.write_text("probe_id\nSP369_p1\nSP369_p2\n", encoding="utf-8")
    result = runner.invoke(app, [
        "pool", "create-design", "TNBCmarker_cand_v1",
        "--probes-file", str(pf), "--repo-root", str(tmp_path),
    ])
    assert result.exit_code == 0, result.output
    pm = PoolManifest.from_yaml(tmp_path / "pools" / "TNBCmarker_cand_v1" / "manifest.yaml")
    assert len(pm.recipe) == 2


def test_create_design_rejects_unknown_probe(tmp_path):
    _write_registry(tmp_path, ["SP369_p1"])
    result = runner.invoke(app, [
        "pool", "create-design", "TNBCmarker_cand_v1",
        "--probe", "SP369_p1", "--probe", "SP369_UNKNOWN",
        "--repo-root", str(tmp_path),
    ])
    assert result.exit_code != 0
    assert "unknown" in result.output.lower() or "registry" in result.output.lower()


def test_create_design_refuses_overwrite(tmp_path):
    _write_registry(tmp_path, ["SP369_p1"])
    _write_design_pool(tmp_path, "TNBCmarker_cand_v1", ["SP369_p1"])
    result = runner.invoke(app, [
        "pool", "create-design", "TNBCmarker_cand_v1",
        "--probe", "SP369_p1", "--repo-root", str(tmp_path),
    ])
    assert result.exit_code != 0
    assert "exist" in result.output.lower()


# ---- promote ----

def test_promote_with_order_binds_all(tmp_path):
    _write_registry(tmp_path, ["SP369_p1", "SP369_p2"])
    _write_design_pool(tmp_path, "TNBCmarker_cand_v1", ["SP369_p1", "SP369_p2"])
    result = runner.invoke(app, [
        "pool", "promote", "TNBCmarker_cand_v1",
        "--order", "O0012_SP369_TNBCmarker_20260613",
        "--concentration", "1e-8", "--repo-root", str(tmp_path),
    ])
    assert result.exit_code == 0, result.output
    pm = PoolManifest.from_yaml(tmp_path / "pools" / "TNBCmarker_cand_v1" / "manifest.yaml")
    assert pm.kind == PoolKind.PHYSICAL.value
    assert all(e.order_id == "O0012_SP369_TNBCmarker_20260613" for e in pm.recipe)
    assert all(e.concentration_M == 1e-8 for e in pm.recipe)


def test_promote_auto_resolves_from_inventories(tmp_path):
    _write_registry(tmp_path, ["SP369_p1", "SP369_p2"])
    _write_design_pool(tmp_path, "TNBCmarker_cand_v1", ["SP369_p1", "SP369_p2"])
    _write_order_inventory(tmp_path, "O0011_SP369_TNBCmarker_20260612", ["SP369_p1"])
    _write_order_inventory(tmp_path, "O0012_SP369_TNBCmarker_20260613", ["SP369_p2"])
    result = runner.invoke(app, [
        "pool", "promote", "TNBCmarker_cand_v1",
        "--auto", "--concentration", "1e-8", "--repo-root", str(tmp_path),
    ])
    assert result.exit_code == 0, result.output
    pm = PoolManifest.from_yaml(tmp_path / "pools" / "TNBCmarker_cand_v1" / "manifest.yaml")
    binding = {e.probe_id: e.order_id for e in pm.recipe}
    assert binding["SP369_p1"] == "O0011_SP369_TNBCmarker_20260612"
    assert binding["SP369_p2"] == "O0012_SP369_TNBCmarker_20260613"


def test_promote_auto_errors_on_unordered_probe(tmp_path):
    _write_registry(tmp_path, ["SP369_p1"])
    _write_design_pool(tmp_path, "TNBCmarker_cand_v1", ["SP369_p1"])
    # no inventory contains SP369_p1
    result = runner.invoke(app, [
        "pool", "promote", "TNBCmarker_cand_v1",
        "--auto", "--repo-root", str(tmp_path),
    ])
    assert result.exit_code != 0
    assert "no order" in result.output.lower() or "unresolved" in result.output.lower() or "not" in result.output.lower()


def test_promote_rejects_physical_pool(tmp_path):
    _write_registry(tmp_path, ["SP369_p1"])
    pm = PoolManifest(pool_id="TNBCmarker_phys_v1", pool_name="p", codebook="SP369",
                      target_well_volume_uL=10.0)  # physical default
    pm.add_probe(order_id="O1", probe_id="SP369_p1", concentration_M=1e-8)
    pm.to_yaml(tmp_path / "pools" / "TNBCmarker_phys_v1" / "manifest.yaml")
    result = runner.invoke(app, [
        "pool", "promote", "TNBCmarker_phys_v1",
        "--order", "O1", "--repo-root", str(tmp_path),
    ])
    assert result.exit_code != 0
    assert "design" in result.output.lower()
