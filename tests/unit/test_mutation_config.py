"""Tests for probe_designer.ext.mutation.config."""
from __future__ import annotations

from pathlib import Path

import pytest

from probe_designer.ext.mutation.config import (
    ALL_CHEMISTRIES,
    ChemistryParams,
    MutationConfig,
    default_chemistry_params,
)


def make_kwargs(tmp_path: Path, **overrides):
    """Minimal valid kwargs for MutationConfig."""
    base = dict(
        input_xlsx=tmp_path / "in.xlsx",
        output_dir=tmp_path / "out",
        backbone_file=tmp_path / "backbone.xlsx",
        start_no=100,
    )
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------


def test_default_chemistries_are_all_three(tmp_path):
    cfg = MutationConfig(**make_kwargs(tmp_path))
    assert cfg.chemistries == list(ALL_CHEMISTRIES)


def test_default_ilock_tm_range_is_50_70(tmp_path):
    """Per 2026-05-10 audit, default iLock Tm window is 50-70 (was 55-75)."""
    cfg = MutationConfig(**make_kwargs(tmp_path))
    assert cfg.chemistry_params["iLock"].tm_range == (50.0, 70.0)


def test_default_chemistry_params_independent_per_call():
    a = default_chemistry_params()
    b = default_chemistry_params()
    a["iLock"].tm_range = (1.0, 2.0)
    assert b["iLock"].tm_range != (1.0, 2.0)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_unknown_chemistry_raises(tmp_path):
    with pytest.raises(ValueError, match="Unknown chemistry"):
        MutationConfig(**make_kwargs(tmp_path, chemistries=["bogus"]))


def test_empty_chemistries_raises(tmp_path):
    with pytest.raises(ValueError, match="must not be empty"):
        MutationConfig(**make_kwargs(tmp_path, chemistries=[]))


def test_inverted_tm_range_raises(tmp_path):
    bad = default_chemistry_params()
    bad["iLock"].tm_range = (80.0, 50.0)
    with pytest.raises(ValueError, match="tm_range"):
        MutationConfig(**make_kwargs(tmp_path, chemistry_params=bad))


def test_must_provide_exactly_one_of_start_no_or_last_no_from(tmp_path):
    # Both → error
    with pytest.raises(ValueError, match="Exactly one"):
        MutationConfig(**make_kwargs(
            tmp_path, last_no_from=tmp_path / "x.txt",
        ))
    # Neither → error
    with pytest.raises(ValueError, match="Exactly one"):
        kw = make_kwargs(tmp_path)
        kw.pop("start_no")
        MutationConfig(**kw)


def test_resolve_start_no_from_file(tmp_path):
    last = tmp_path / "last_no.txt"
    last.write_text("369", encoding="utf-8")
    kw = make_kwargs(tmp_path)
    kw.pop("start_no")
    kw["last_no_from"] = last
    cfg = MutationConfig(**kw)
    assert cfg.resolve_start_no() == 370


def test_resolve_start_no_from_explicit(tmp_path):
    cfg = MutationConfig(**make_kwargs(tmp_path, start_no=500))
    assert cfg.resolve_start_no() == 500


# ---------------------------------------------------------------------------
# YAML round-trip
# ---------------------------------------------------------------------------


def test_yaml_loader_with_per_chemistry_override(tmp_path):
    yaml_text = """
input_xlsx: in.xlsx
output_dir: out/
backbone_file: backbone.xlsx
chemistries: [iLock]
chemistry_params:
  iLock:
    tm_range: [45, 80]
arm_len_range: [8, 35]
mfe_min: -25.0
start_no: 100
"""
    p = tmp_path / "config.yaml"
    p.write_text(yaml_text, encoding="utf-8")

    cfg = MutationConfig.from_yaml(p)
    assert cfg.chemistries == ["iLock"]
    assert cfg.chemistry_params["iLock"].tm_range == (45.0, 80.0)
    # Other chemistries keep their defaults even if not in this run's chemistries
    assert cfg.chemistry_params["dRNA"].tm_range == (55.0, 75.0)
    assert cfg.arm_len_range == (8, 35)
    assert cfg.mfe_min == -25.0


def test_legacy_chemistry_token_normalized_to_drna(tmp_path):
    """Legacy ``mRNA_noiLock`` is accepted but normalized to ``dRNA``."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        cfg = MutationConfig(
            **make_kwargs(tmp_path, chemistries=["mRNA_noiLock"]),
        )
    assert cfg.chemistries == ["dRNA"]


def test_codebook_field_default_none(tmp_path):
    cfg = MutationConfig(**make_kwargs(tmp_path))
    assert cfg.codebook is None
