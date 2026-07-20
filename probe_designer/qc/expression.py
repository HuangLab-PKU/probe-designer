"""scRNA-seq expression QC for designed panels.

Phase 3 (2026-05-14). Numbers below are anchored to the Xenium tech note
``CG000643_DesigningCustomXeniumPanels_TechnicalNote_RevB`` Table 1:

* Low cutoff: <0.1 mean transcripts per cell → drop (won't be detected).
* High cutoff: >50 (Chromium v2) or >100 (v3.1) mean transcripts → EHEG
  (extremely highly expressed; risks optical crowding).
* Cell-type-marker viability: ≥2 (v2) or ≥4 (v3.1) → safe to use as a
  standalone marker.

Co-expressed EHEGs: per the tech note, two EHEGs co-expressed in the
same cell type knock out detection of up to 20% of co-expressed genes.
We flag pairs whose per-cell-type expression overlap (Jaccard on the
"cell types where expressed" sets) exceeds a threshold.
"""
from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Literal, Optional, Sequence, Tuple

import numpy as np


logger = logging.getLogger(__name__)


# Defaults anchored to Xenium tech note Table 1.
LOW_CUTOFF_DEFAULT: float = 0.1
HIGH_CUTOFF_DEFAULT_V31: float = 100.0
HIGH_CUTOFF_DEFAULT_V2: float = 50.0
MARKER_MIN_DEFAULT_V31: float = 4.0
MARKER_MIN_DEFAULT_V2: float = 2.0


FlagKind = Literal["low", "high", "marker_unsafe", "missing_in_reference"]


@dataclass
class ExpressionFlag:
    """One panel-level concern about a gene."""
    gene: str
    cell_type: Optional[str]   # None for whole-panel-level flags
    kind: FlagKind
    value: float                # the numeric driver (e.g. mean transcripts/cell)
    note: str = ""


@dataclass
class ExpressionQC:
    """Result of an expression QC run.

    Attributes:
        panel_genes: the panel genes we evaluated (preserves caller order).
        cell_types: the cell types found in the reference.
        cell_type_tp10k: dict of dicts ``{gene -> {cell_type -> tp10k}}``.
        mean_transcripts_per_cell: ``{gene -> {cell_type -> mean count}}``.
        panel_tp10k_per_cell_type: ``{cell_type -> sum of panel-gene TP10K}``.
        flags: list of ``ExpressionFlag`` ordered by (cell_type, gene).
        ehec_pairs: list of (gene_a, gene_b, jaccard) for co-expressed EHEGs.
    """
    panel_genes: List[str]
    cell_types: List[str]
    cell_type_tp10k: Dict[str, Dict[str, float]]
    mean_transcripts_per_cell: Dict[str, Dict[str, float]]
    panel_tp10k_per_cell_type: Dict[str, float]
    flags: List[ExpressionFlag]
    ehec_pairs: List[Tuple[str, str, float]]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "panel_genes": self.panel_genes,
            "cell_types": self.cell_types,
            "cell_type_tp10k": self.cell_type_tp10k,
            "mean_transcripts_per_cell": self.mean_transcripts_per_cell,
            "panel_tp10k_per_cell_type": self.panel_tp10k_per_cell_type,
            "flags": [asdict(f) for f in self.flags],
            "ehec_pairs": [{"gene_a": a, "gene_b": b, "jaccard": j}
                           for a, b, j in self.ehec_pairs],
        }


def _validate_inputs(
    matrix: np.ndarray, gene_names: Sequence[str],
    cell_type_labels: Sequence[str],
) -> None:
    if matrix.ndim != 2:
        raise ValueError(f"matrix must be 2D, got shape {matrix.shape}")
    n_cells, n_genes = matrix.shape
    if len(gene_names) != n_genes:
        raise ValueError(
            f"gene_names length {len(gene_names)} != matrix columns {n_genes}"
        )
    if len(cell_type_labels) != n_cells:
        raise ValueError(
            f"cell_type_labels length {len(cell_type_labels)} != matrix rows {n_cells}"
        )


def compute_expression_qc_from_matrix(
    matrix: np.ndarray,
    *,
    gene_names: Sequence[str],
    cell_type_labels: Sequence[str],
    panel_genes: Sequence[str],
    low_cutoff: float = LOW_CUTOFF_DEFAULT,
    high_cutoff: Optional[float] = None,
    marker_min: Optional[float] = None,
    chromium_version: Literal["v2", "v3.1"] = "v3.1",
    ehec_jaccard: float = 0.5,
) -> ExpressionQC:
    """Run panel QC on a raw expression matrix.

    Args:
        matrix: cells × genes counts matrix (NumPy ndarray). Sparse inputs
            should be ``toarray()``-ed by the caller.
        gene_names: the column labels of ``matrix``. Used to look up
            ``panel_genes`` rows.
        cell_type_labels: per-cell category. Used for grouping.
        panel_genes: the panel genes to evaluate.
        low_cutoff: genes with ``mean transcripts/cell < low_cutoff`` in
            every cell type get a ``"low"`` flag.
        high_cutoff: genes above this threshold in ANY cell type get
            an EHEG ``"high"`` flag. Defaults to Xenium's v3.1 (100) /
            v2 (50) based on ``chromium_version``.
        marker_min: standalone-marker viability floor. Defaults to
            Xenium's v3.1 (4) / v2 (2). If a panel gene's maximum
            cross-cell-type mean is below this, we flag ``"marker_unsafe"``.
        chromium_version: drives the high_cutoff / marker_min defaults.
        ehec_jaccard: minimum Jaccard overlap on "cell types where
            expressed at all" for two EHEG genes to be flagged as
            co-expressed (default 0.5).
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    _validate_inputs(matrix, gene_names, cell_type_labels)

    if high_cutoff is None:
        high_cutoff = (HIGH_CUTOFF_DEFAULT_V31
                       if chromium_version == "v3.1"
                       else HIGH_CUTOFF_DEFAULT_V2)
    if marker_min is None:
        marker_min = (MARKER_MIN_DEFAULT_V31
                      if chromium_version == "v3.1"
                      else MARKER_MIN_DEFAULT_V2)

    gene_to_col = {g: i for i, g in enumerate(gene_names)}
    cell_types = sorted(set(cell_type_labels))

    # Cell-type indices into matrix rows.
    cell_type_to_rows: Dict[str, np.ndarray] = {
        ct: np.array([i for i, lab in enumerate(cell_type_labels) if lab == ct],
                     dtype=int)
        for ct in cell_types
    }

    # mean transcripts/cell per (gene, cell_type)
    mean_counts: Dict[str, Dict[str, float]] = {}
    tp10k: Dict[str, Dict[str, float]] = {}

    # Precompute, per cell_type, the mean profile across ALL reference genes
    # (used to normalize to TP10K).
    profile_sums: Dict[str, float] = {}
    full_profile_mean: Dict[str, np.ndarray] = {}
    for ct, rows in cell_type_to_rows.items():
        if rows.size == 0:
            full_profile_mean[ct] = np.zeros(len(gene_names), dtype=np.float64)
            profile_sums[ct] = 0.0
            continue
        mean_vec = matrix[rows, :].mean(axis=0)
        full_profile_mean[ct] = mean_vec
        profile_sums[ct] = float(mean_vec.sum())

    missing_panel: List[str] = []
    for gene in panel_genes:
        col = gene_to_col.get(gene)
        if col is None:
            missing_panel.append(gene)
            mean_counts[gene] = {ct: 0.0 for ct in cell_types}
            tp10k[gene] = {ct: 0.0 for ct in cell_types}
            continue
        mean_counts[gene] = {}
        tp10k[gene] = {}
        for ct in cell_types:
            cmean = float(full_profile_mean[ct][col])
            mean_counts[gene][ct] = cmean
            total = profile_sums[ct]
            tp10k[gene][ct] = (cmean / total * 10_000.0) if total > 0 else 0.0

    # Panel-level utilization per cell type.
    panel_tp10k_per_cell_type: Dict[str, float] = {
        ct: sum(tp10k[g][ct] for g in panel_genes) for ct in cell_types
    }

    # Flags
    flags: List[ExpressionFlag] = []
    for gene in missing_panel:
        flags.append(ExpressionFlag(
            gene=gene, cell_type=None, kind="missing_in_reference",
            value=0.0, note="gene name not present in scRNA-seq reference",
        ))
    for gene in panel_genes:
        if gene in missing_panel:
            continue
        max_mean = max(mean_counts[gene].values()) if mean_counts[gene] else 0.0
        # Low: below cutoff in every cell type
        if max_mean < low_cutoff:
            flags.append(ExpressionFlag(
                gene=gene, cell_type=None, kind="low",
                value=max_mean,
                note=f"max mean transcripts/cell {max_mean:.3f} < low_cutoff {low_cutoff}",
            ))
        # Marker-unsafe: max < marker_min (only meaningful if not also "low")
        elif max_mean < marker_min:
            flags.append(ExpressionFlag(
                gene=gene, cell_type=None, kind="marker_unsafe",
                value=max_mean,
                note=f"max mean transcripts/cell {max_mean:.3f} < marker_min {marker_min}",
            ))
        # High: any single cell type above high_cutoff
        for ct, cmean in mean_counts[gene].items():
            if cmean > high_cutoff:
                flags.append(ExpressionFlag(
                    gene=gene, cell_type=ct, kind="high",
                    value=cmean,
                    note=f"mean transcripts/cell {cmean:.1f} > high_cutoff {high_cutoff}",
                ))

    # Co-expressed EHEG pairs (Jaccard on expressed-cell-type sets).
    ehec_genes = sorted({f.gene for f in flags if f.kind == "high"})
    ehec_pairs: List[Tuple[str, str, float]] = []
    def expressed_set(gene: str) -> set[str]:
        return {ct for ct, c in mean_counts[gene].items() if c > low_cutoff}
    for i, g1 in enumerate(ehec_genes):
        s1 = expressed_set(g1)
        for g2 in ehec_genes[i + 1:]:
            s2 = expressed_set(g2)
            if not s1 and not s2:
                continue
            inter = len(s1 & s2)
            union = len(s1 | s2)
            if union == 0:
                continue
            j = inter / union
            if j >= ehec_jaccard:
                ehec_pairs.append((g1, g2, j))

    return ExpressionQC(
        panel_genes=list(panel_genes),
        cell_types=cell_types,
        cell_type_tp10k=tp10k,
        mean_transcripts_per_cell=mean_counts,
        panel_tp10k_per_cell_type=panel_tp10k_per_cell_type,
        flags=flags,
        ehec_pairs=ehec_pairs,
    )


def compute_expression_qc_from_anndata(
    adata: Any,
    *,
    panel_genes: Sequence[str],
    cell_type_key: str = "celltype",
    **kwargs: Any,
) -> ExpressionQC:
    """Adapter for AnnData / 10x HDF5 inputs.

    Requires anndata to be installed; raises ImportError with an install
    hint otherwise.
    """
    try:
        import anndata  # noqa: F401  (proves the package is importable)
    except ImportError as exc:
        raise ImportError(
            "anndata is required for AnnData input. Install via "
            "`mamba install -n probe-design anndata scanpy`."
        ) from exc

    # Pull the matrix; densify sparse if needed.
    X = adata.X
    if hasattr(X, "toarray"):
        matrix = X.toarray()
    else:
        matrix = np.asarray(X)
    gene_names = list(adata.var_names)
    cell_type_labels = list(adata.obs[cell_type_key].astype(str))
    return compute_expression_qc_from_matrix(
        matrix,
        gene_names=gene_names,
        cell_type_labels=cell_type_labels,
        panel_genes=panel_genes,
        **kwargs,
    )
