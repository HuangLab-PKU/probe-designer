"""Tests for ``probe-design assemble`` subcommand."""
from __future__ import annotations

import json

import pandas as pd
import pytest
from typer.testing import CliRunner

from probe_designer.cli.app import app


runner = CliRunner()


@pytest.fixture
def backbone_xlsx(tmp_path):
    """Backbone file with one sequence row."""
    p = tmp_path / "backbone.xlsx"
    pd.DataFrame({
        "No.": ["1"],
        "Sequence": ["GTTTTTGGTAAGCTTCGGATCCTCAGACGGAAGACTC"],
    }).to_excel(p, index=False)
    return p


@pytest.fixture
def gene_info_xlsx(tmp_path):
    """gene_info with a single assignment: gene ACTB -> backbone No. '1'."""
    p = tmp_path / "gene_info.xlsx"
    pd.DataFrame({
        "gene_name": ["ACTB"],
        "No.": ["1"],
    }).to_excel(p, index=False)
    return p


@pytest.fixture
def binding_sites_json(tmp_path):
    p = tmp_path / "sites.json"
    payload = {
        "ACTB": [
            {
                "gene_name": "ACTB",
                "arm_3prime": "ACGTACGTACGTACGTACGT",
                "arm_5prime": "TGCATGCATGCATGCATGCA",
                "st": 100, "en": 140,
                "g_content": 0.5, "tm": 60.0,
                "tm_3prime": 59.0, "tm_5prime": 61.0,
            },
        ],
    }
    p.write_text(json.dumps(payload), encoding="utf-8")
    return p


class TestAssembleSubcommand:
    def test_happy_path_json_input(
        self, tmp_path, binding_sites_json, gene_info_xlsx, backbone_xlsx
    ):
        out = tmp_path / "out"
        result = runner.invoke(app, [
            "assemble",
            "--binding-sites", str(binding_sites_json),
            "--gene-info", str(gene_info_xlsx),
            "--backbone", str(backbone_xlsx),
            "--output", str(out),
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()
        xlsx_files = list(out.glob("*.xlsx"))
        assert xlsx_files, f"no xlsx in {out}: {list(out.iterdir())}"

    def test_empty_binding_sites_exits_nonzero(
        self, tmp_path, gene_info_xlsx, backbone_xlsx
    ):
        empty = tmp_path / "empty.json"
        empty.write_text("{}", encoding="utf-8")
        out = tmp_path / "out"
        result = runner.invoke(app, [
            "assemble",
            "--binding-sites", str(empty),
            "--gene-info", str(gene_info_xlsx),
            "--backbone", str(backbone_xlsx),
            "--output", str(out),
        ])
        assert result.exit_code != 0

    def test_bad_ilock_value_rejected(
        self, tmp_path, binding_sites_json, gene_info_xlsx, backbone_xlsx
    ):
        result = runner.invoke(app, [
            "assemble",
            "--binding-sites", str(binding_sites_json),
            "--gene-info", str(gene_info_xlsx),
            "--backbone", str(backbone_xlsx),
            "--output", str(tmp_path / "out"),
            "--ilock", "nonsense",
        ])
        assert result.exit_code != 0
        assert "ilock" in result.output.lower() or "nonsense" in result.output.lower()
