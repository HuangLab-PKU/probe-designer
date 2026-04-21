"""Integration-level tests for ``probe-design merge``.

Builds synthetic XLSX fixtures in tmp_path and drives the command through
Typer's CliRunner — no network, no real pandas IO against large files.
"""
from __future__ import annotations

import pandas as pd
import pytest
from typer.testing import CliRunner

from probe_designer.cli.app import app


runner = CliRunner()


def _write_xlsx(path, rows):
    pd.DataFrame(rows).to_excel(path, index=False)


@pytest.fixture
def results_tree(tmp_path):
    """Two subdirs each with a filtered_binding_sites.xlsx.

    - run_a covers ACTB (pos 100) + GAPDH (pos 50)
    - run_b covers ACTB (pos 100 dup) + ACTB (pos 500) + UNKNOWN_GENE (filtered out)
    """
    (tmp_path / "run_a").mkdir()
    (tmp_path / "run_b").mkdir()

    _write_xlsx(tmp_path / "run_a" / "filtered_binding_sites.xlsx", [
        {"gene_name": "ACTB", "st": 100, "en": 140, "tm": 60.0, "g_content": 0.4},
        {"gene_name": "GAPDH", "st": 50, "en": 90, "tm": 62.0, "g_content": 0.45},
    ])
    _write_xlsx(tmp_path / "run_b" / "filtered_binding_sites.xlsx", [
        {"gene_name": "actb", "st": 100, "en": 140, "tm": 60.0, "g_content": 0.4},  # dup (case-insensitive gene)
        {"gene_name": "ACTB", "st": 500, "en": 540, "tm": 61.5, "g_content": 0.5},
        {"gene_name": "UNKNOWN_GENE", "st": 1, "en": 40, "tm": 60.0, "g_content": 0.5},
    ])
    return tmp_path


@pytest.fixture
def gene_info_file(tmp_path):
    p = tmp_path / "gene_info.xlsx"
    pd.DataFrame({
        "No.": [1, 2],
        "gene_name": ["ACTB", "GAPDH"],  # order matters for sort
    }).to_excel(p, index=False)
    return p


class TestMergeSubcommand:
    def test_happy_path(self, results_tree, gene_info_file, tmp_path):
        out = tmp_path / "merged.xlsx"
        result = runner.invoke(app, [
            "merge",
            "--results-dir", str(results_tree),
            "--gene-info", str(gene_info_file),
            "--output", str(out),
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()

        df = pd.read_excel(out)
        # ACTB has 2 unique positions after dedup, GAPDH has 1, UNKNOWN filtered out -> 3 rows
        assert len(df) == 3
        assert set(df["gene_name"].str.upper()) == {"ACTB", "GAPDH"}

    def test_missing_genes_file_written(self, tmp_path, results_tree):
        # gene_info lists a gene that isn't in any results file
        info = tmp_path / "gi.xlsx"
        pd.DataFrame({
            "No.": [1, 2, 3],
            "gene_name": ["ACTB", "GAPDH", "MISSING_GENE"],
        }).to_excel(info, index=False)

        out = tmp_path / "merged.xlsx"
        result = runner.invoke(app, [
            "merge", "-r", str(results_tree), "-g", str(info), "-o", str(out),
        ])
        assert result.exit_code == 0, result.output
        missing_txt = out.parent / "missing_genes.txt"
        assert missing_txt.exists()
        assert "MISSING_GENE" in missing_txt.read_text(encoding="utf-8")
        assert "ACTB" not in missing_txt.read_text(encoding="utf-8")

    def test_sorted_by_gene_info_order(self, tmp_path):
        # Put GAPDH first in gene_info; merged output should list GAPDH rows first.
        (tmp_path / "r1").mkdir()
        _write_xlsx(tmp_path / "r1" / "filtered_binding_sites.xlsx", [
            {"gene_name": "ACTB", "st": 10, "en": 50, "tm": 60.0, "g_content": 0.4},
            {"gene_name": "GAPDH", "st": 20, "en": 60, "tm": 60.0, "g_content": 0.4},
        ])
        info = tmp_path / "gi.xlsx"
        pd.DataFrame({
            "No.": [1, 2],
            "gene_name": ["GAPDH", "ACTB"],
        }).to_excel(info, index=False)
        out = tmp_path / "merged.xlsx"
        result = runner.invoke(app, [
            "merge", "-r", str(tmp_path), "-g", str(info), "-o", str(out),
        ])
        assert result.exit_code == 0, result.output
        df = pd.read_excel(out)
        assert df["gene_name"].str.upper().tolist() == ["GAPDH", "ACTB"]

    def test_empty_results_dir_exits_nonzero(self, tmp_path, gene_info_file):
        empty = tmp_path / "empty"
        empty.mkdir()
        out = tmp_path / "merged.xlsx"
        result = runner.invoke(app, [
            "merge", "-r", str(empty), "-g", str(gene_info_file), "-o", str(out),
        ])
        assert result.exit_code != 0

    def test_custom_missing_output_path(self, results_tree, tmp_path):
        info = tmp_path / "gi.xlsx"
        pd.DataFrame({
            "No.": [1, 2, 3],
            "gene_name": ["ACTB", "GAPDH", "MISSING_GENE"],
        }).to_excel(info, index=False)
        out = tmp_path / "merged.xlsx"
        missing_custom = tmp_path / "subdir" / "my_missing.txt"
        missing_custom.parent.mkdir(parents=True)
        result = runner.invoke(app, [
            "merge", "-r", str(results_tree), "-g", str(info), "-o", str(out),
            "-m", str(missing_custom),
        ])
        assert result.exit_code == 0, result.output
        assert missing_custom.exists()
        assert "MISSING_GENE" in missing_custom.read_text(encoding="utf-8")
