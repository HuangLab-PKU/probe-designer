"""Tests for probe_designer.probe_assembly.load_binding_sites."""
from __future__ import annotations

import json

import pandas as pd
import pytest

from probe_designer.probe_assembly import load_binding_sites


class TestJsonPath:
    def test_reads_json_dict_verbatim(self, tmp_path):
        payload = {
            "ACTB": [{"arm_3prime": "ACGT", "arm_5prime": "TGCA", "st": 100, "en": 140}]
        }
        p = tmp_path / "sites.json"
        p.write_text(json.dumps(payload), encoding="utf-8")

        result = load_binding_sites(p)

        assert result == payload

    def test_json_with_empty_dict_ok(self, tmp_path):
        p = tmp_path / "empty.json"
        p.write_text("{}", encoding="utf-8")
        assert load_binding_sites(p) == {}


class TestExcelPath:
    def _write_xlsx(self, tmp_path, rows):
        df = pd.DataFrame(rows)
        p = tmp_path / "sites.xlsx"
        df.to_excel(p, index=False)
        return p

    def test_reads_xlsx_and_groups_by_gene(self, tmp_path):
        p = self._write_xlsx(tmp_path, [
            {"gene_name": "ACTB", "arm3": "ACGT", "arm5": "TGCA",
             "st": 100, "en": 140, "g_content": 0.5, "tm": 60.0,
             "tm3": 58.0, "tm5": 61.0},
            {"gene_name": "ACTB", "arm3": "AAAA", "arm5": "TTTT",
             "st": 200, "en": 240, "g_content": 0.0, "tm": 55.0,
             "tm3": 54.0, "tm5": 55.5},
            {"gene_name": "GAPDH", "arm3": "CCCC", "arm5": "GGGG",
             "st": 10, "en": 50, "g_content": 1.0, "tm": 72.0,
             "tm3": 71.0, "tm5": 72.5},
        ])

        result = load_binding_sites(p)

        assert set(result.keys()) == {"ACTB", "GAPDH"}
        assert len(result["ACTB"]) == 2
        assert len(result["GAPDH"]) == 1
        # Required fields renamed to canonical keys
        site = result["ACTB"][0]
        assert site["arm_3prime"] == "ACGT"
        assert site["arm_5prime"] == "TGCA"
        assert site["tm_3prime"] == 58.0
        assert site["tm_5prime"] == 61.0

    def test_xlsx_missing_required_columns_raises(self, tmp_path):
        p = self._write_xlsx(tmp_path, [{"gene_name": "ACTB", "arm3": "A"}])
        with pytest.raises(ValueError, match="missing required columns"):
            load_binding_sites(p)

    def test_optional_isoform_overlap_propagates(self, tmp_path):
        p = self._write_xlsx(tmp_path, [
            {"gene_name": "ACTB", "arm3": "A", "arm5": "T",
             "st": 1, "en": 40, "g_content": 0.5, "tm": 60.0,
             "tm3": 58.0, "tm5": 61.0, "isoform_overlap_num": 3},
        ])
        result = load_binding_sites(p)
        assert result["ACTB"][0]["isoform_overlap_num"] == 3


class TestDispatch:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_binding_sites(tmp_path / "nope.json")

    def test_unknown_extension_raises(self, tmp_path):
        p = tmp_path / "sites.txt"
        p.write_text("anything")
        with pytest.raises(ValueError, match="Unsupported"):
            load_binding_sites(p)
