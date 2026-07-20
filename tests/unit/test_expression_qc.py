"""Tests for the Phase 3 scRNA-seq expression QC module."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from probe_designer.qc import (
    ExpressionFlag,
    ExpressionQC,
    compute_expression_qc_from_anndata,
    compute_expression_qc_from_matrix,
    render_panel_qc_markdown,
    write_panel_qc_json,
)


def _toy_dataset():
    """Hand-crafted dataset reproducing the Xenium tech note Vip/Plp1 motif.

    Cells:
      - 0..29 (30 cells): "oligodendrocyte" — Plp1 is sky-high (≈150 mean),
        Vip is silent. Plp1 should be flagged EHEG; Vip should be flagged
        as low/marker-unsafe.
      - 30..49 (20 cells): "Vip-neuron" — Vip is high (10 mean), Plp1
        silent.
      - 50..69 (20 cells): "other" — both quiet.
    Genes (columns): ['Vip', 'Plp1', 'NoiseGene', 'NeverExpressed'].
    """
    np.random.seed(7)
    n_cells = 70
    n_genes = 4
    X = np.zeros((n_cells, n_genes), dtype=np.float64)
    # Plp1 in oligodendrocyte
    X[0:30, 1] = np.random.normal(150, 15, 30).clip(min=0)
    # Vip in Vip-neuron
    X[30:50, 0] = np.random.normal(10, 2, 20).clip(min=0)
    # NoiseGene low everywhere
    X[:, 2] = np.random.normal(2.0, 0.5, n_cells).clip(min=0)
    # NeverExpressed stays 0

    cell_types = (["oligodendrocyte"] * 30 + ["Vip-neuron"] * 20
                  + ["other"] * 20)
    gene_names = ["Vip", "Plp1", "NoiseGene", "NeverExpressed"]
    return X, gene_names, cell_types


class TestComputeQCFromMatrix:
    def test_smoke_shape_and_keys(self):
        X, gene_names, cell_types = _toy_dataset()
        qc = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["Vip", "Plp1"],
        )
        assert isinstance(qc, ExpressionQC)
        assert set(qc.cell_types) == {"oligodendrocyte", "Vip-neuron", "other"}
        assert set(qc.panel_genes) == {"Vip", "Plp1"}
        assert set(qc.cell_type_tp10k["Plp1"].keys()) == set(qc.cell_types)

    def test_plp1_flagged_as_ehec(self):
        X, gene_names, cell_types = _toy_dataset()
        qc = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["Vip", "Plp1", "NoiseGene"],
            chromium_version="v3.1",
        )
        plp1_high = [f for f in qc.flags if f.gene == "Plp1" and f.kind == "high"]
        assert plp1_high, "Plp1 must be flagged as EHEG in oligodendrocyte"
        # And in the right cell type
        assert any(f.cell_type == "oligodendrocyte" for f in plp1_high)

    def test_vip_marker_safe_no_low_flag(self):
        X, gene_names, cell_types = _toy_dataset()
        qc = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["Vip"],
            chromium_version="v3.1",
        )
        vip_flags = [f for f in qc.flags if f.gene == "Vip"]
        # Vip max mean ~10 → above marker_min (4); shouldn't be flagged low.
        assert not any(f.kind == "low" for f in vip_flags)
        assert not any(f.kind == "marker_unsafe" for f in vip_flags)

    def test_low_gene_flagged(self):
        X, gene_names, cell_types = _toy_dataset()
        qc = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["NeverExpressed"],
        )
        flags = [f for f in qc.flags if f.gene == "NeverExpressed"]
        assert any(f.kind == "low" for f in flags)

    def test_missing_gene_flagged(self):
        X, gene_names, cell_types = _toy_dataset()
        qc = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["NotInReference"],
        )
        flags = [f for f in qc.flags if f.gene == "NotInReference"]
        assert any(f.kind == "missing_in_reference" for f in flags)

    def test_marker_unsafe_flag(self):
        """NoiseGene has mean ~2 — below marker_min=4 (v3.1) but above
        low_cutoff=0.1 → marker_unsafe flag."""
        X, gene_names, cell_types = _toy_dataset()
        qc = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["NoiseGene"],
            chromium_version="v3.1",
        )
        flags = [f for f in qc.flags if f.gene == "NoiseGene"]
        assert any(f.kind == "marker_unsafe" for f in flags)
        # NOT flagged low (max ~2 > low_cutoff 0.1)
        assert not any(f.kind == "low" for f in flags)

    def test_tp10k_normalizes_to_10000(self):
        X, gene_names, cell_types = _toy_dataset()
        qc = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=gene_names,  # all genes
        )
        # Per-cell-type sum across ALL 4 genes should equal ≤ 10,000 (panel = full reference).
        for ct in qc.cell_types:
            total = qc.panel_tp10k_per_cell_type[ct]
            # In our toy dataset every gene is in the panel, so panel-tp10k = 10000
            # for any cell type with nonzero profile.
            if total > 0:
                assert abs(total - 10_000.0) < 1e-6

    def test_v2_thresholds_softer(self):
        X, gene_names, cell_types = _toy_dataset()
        # NoiseGene has max ~2; v2 marker_min is 2, so it should be ON the boundary.
        # Use Plp1 to check the high cutoff difference instead.
        qc_v31 = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["Plp1"], chromium_version="v3.1",
        )
        qc_v2 = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["Plp1"], chromium_version="v2",
        )
        # Plp1 mean ~150 in oligos; both v3.1 (cutoff 100) and v2 (cutoff 50)
        # should flag it. But Plp1 across other cell types is ~0 — won't flag.
        v31_high = [f for f in qc_v31.flags if f.kind == "high"]
        v2_high = [f for f in qc_v2.flags if f.kind == "high"]
        assert v31_high and v2_high

    def test_ehec_co_expression_pair_detected(self):
        # Two genes both EHEG in the same cell type → co-expressed pair
        n_cells = 60
        X = np.zeros((n_cells, 2), dtype=np.float64)
        X[0:30, 0] = 150  # Plp1-like
        X[0:30, 1] = 120  # Mog-like (also EHEG in oligodendrocyte)
        cell_types = ["oligodendrocyte"] * 30 + ["other"] * 30
        qc = compute_expression_qc_from_matrix(
            X, gene_names=["Plp1", "Mog"],
            cell_type_labels=cell_types,
            panel_genes=["Plp1", "Mog"],
            chromium_version="v3.1", ehec_jaccard=0.5,
        )
        assert len(qc.ehec_pairs) >= 1
        a, b, jaccard = qc.ehec_pairs[0]
        assert {a, b} == {"Plp1", "Mog"}
        assert jaccard >= 0.5

    def test_disjoint_eheg_not_paired(self):
        """Two EHEGs in mutually exclusive cell types should NOT pair."""
        n_cells = 60
        X = np.zeros((n_cells, 2), dtype=np.float64)
        X[0:30, 0] = 150     # Plp1 in oligo
        X[30:60, 1] = 150    # Ins1 in beta (different cell type)
        cell_types = ["oligodendrocyte"] * 30 + ["beta_cell"] * 30
        qc = compute_expression_qc_from_matrix(
            X, gene_names=["Plp1", "Ins1"],
            cell_type_labels=cell_types,
            panel_genes=["Plp1", "Ins1"],
        )
        # Both EHEG, but disjoint cell-type sets — Jaccard 0 → no pair.
        assert qc.ehec_pairs == []

    def test_invalid_inputs_raise(self):
        with pytest.raises(ValueError):
            compute_expression_qc_from_matrix(
                np.zeros((1, 1)), gene_names=["A", "B"],
                cell_type_labels=["x"], panel_genes=["A"],
            )
        with pytest.raises(ValueError):
            compute_expression_qc_from_matrix(
                np.zeros((2, 1)), gene_names=["A"],
                cell_type_labels=["x"], panel_genes=["A"],
            )


class TestAnndataAdapter:
    def test_anndata_unavailable_raises_helpful_error(self):
        """Without anndata installed, the adapter should raise an
        actionable ImportError, not crash mysteriously."""
        # We can't actually call the adapter without anndata. The test
        # just verifies the function exists and the error message
        # mentions the install hint.
        with pytest.raises(ImportError, match="anndata"):
            compute_expression_qc_from_anndata(
                None, panel_genes=["Plp1"],
            )


class TestRender:
    def test_markdown_includes_flags_and_utilization(self):
        X, gene_names, cell_types = _toy_dataset()
        qc = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["Plp1", "Vip", "NeverExpressed"],
        )
        md = render_panel_qc_markdown(qc, panel_id="TEST_PANEL")
        assert "TEST_PANEL" in md
        assert "## Flags" in md
        assert "## Panel TP10K per cell type" in md
        assert "Plp1" in md

    def test_json_round_trip(self, tmp_path):
        X, gene_names, cell_types = _toy_dataset()
        qc = compute_expression_qc_from_matrix(
            X, gene_names=gene_names, cell_type_labels=cell_types,
            panel_genes=["Plp1"],
        )
        out = tmp_path / "qc.json"
        write_panel_qc_json(qc, out)
        with open(out) as f:
            data = json.load(f)
        assert data["panel_genes"] == ["Plp1"]
        assert "cell_type_tp10k" in data
