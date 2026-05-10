"""Tests for probe_designer.ext.tcr.config."""
from __future__ import annotations

from pathlib import Path

import pytest

from probe_designer.ext.tcr.config import ALL_CHEMISTRIES, TcrConfig


def make_kwargs(tmp_path: Path, **overrides):
    base = dict(
        input_xlsx=tmp_path / "in.xlsx",
        output_dir=tmp_path / "out",
        backbone_file=tmp_path / "backbone.xlsx",
        start_no=100,
    )
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# Defaults reflect TCR convention: mRNA-only by default
# ---------------------------------------------------------------------------


def test_default_chemistry_is_mrna_only(tmp_path):
    cfg = TcrConfig(**make_kwargs(tmp_path))
    assert cfg.chemistries == ["mRNA"]
    assert not cfg.has_cdna()


def test_chemistries_with_cdna_triggers_has_cdna(tmp_path):
    cfg = TcrConfig(**make_kwargs(tmp_path, chemistries=["mRNA", "cDNA"]))
    assert cfg.has_cdna()


def test_default_geometry_40_20_3(tmp_path):
    cfg = TcrConfig(**make_kwargs(tmp_path))
    assert cfg.bds_len == 40
    assert cfg.arm_len == 20
    assert cfg.sites_per_clone == 3


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_unknown_chemistry_raises(tmp_path):
    with pytest.raises(ValueError, match="Unknown chemistry"):
        TcrConfig(**make_kwargs(tmp_path, chemistries=["bogus"]))


def test_empty_chemistries_raises(tmp_path):
    with pytest.raises(ValueError, match="must not be empty"):
        TcrConfig(**make_kwargs(tmp_path, chemistries=[]))


def test_arm_len_must_be_half_bds(tmp_path):
    with pytest.raises(ValueError, match="arm_len.*bds_len"):
        TcrConfig(**make_kwargs(tmp_path, bds_len=40, arm_len=15))


def test_inverted_tm_range_raises(tmp_path):
    with pytest.raises(ValueError, match="tm_range"):
        TcrConfig(**make_kwargs(tmp_path, tm_range=(60.0, 50.0)))


def test_must_provide_exactly_one_of_start_no_or_last_no_from(tmp_path):
    with pytest.raises(ValueError, match="Exactly one"):
        TcrConfig(**make_kwargs(tmp_path, last_no_from=tmp_path / "x.txt"))
    with pytest.raises(ValueError, match="Exactly one"):
        kw = make_kwargs(tmp_path)
        kw.pop("start_no")
        TcrConfig(**kw)


def test_resolve_start_no_from_file(tmp_path):
    last = tmp_path / "last_no.txt"
    last.write_text("369", encoding="utf-8")
    kw = make_kwargs(tmp_path)
    kw.pop("start_no")
    kw["last_no_from"] = last
    cfg = TcrConfig(**kw)
    assert cfg.resolve_start_no() == 370


# ---------------------------------------------------------------------------
# YAML round-trip
# ---------------------------------------------------------------------------


def test_yaml_loader_with_dual_chemistry(tmp_path):
    yaml_text = """
input_xlsx: in.xlsx
output_dir: out/
backbone_file: backbone.xlsx
chemistries: [mRNA, cDNA]
bds_len: 40
arm_len: 20
tm_range: [55, 65]
sites_per_clone: 5
start_no: 200
rt_primer_gap: 12
"""
    p = tmp_path / "config.yaml"
    p.write_text(yaml_text, encoding="utf-8")
    cfg = TcrConfig.from_yaml(p)
    assert cfg.chemistries == ["mRNA", "cDNA"]
    assert cfg.has_cdna()
    assert cfg.tm_range == (55.0, 65.0)
    assert cfg.sites_per_clone == 5
    assert cfg.rt_primer_gap == 12
