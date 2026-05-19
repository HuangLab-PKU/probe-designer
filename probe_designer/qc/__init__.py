"""Panel-level QC against scRNA-seq references.

Phase 3 (2026-05-14). Given a scRNA-seq matrix (cells × genes with a
cell_type annotation) and the target list of a designed panel, the QC
module:

* computes per-cell-type, per-panel-gene TP10K (transcripts per 10,000),
* flags genes that are too low to detect or so high they crowd out
  neighbors (EHEGs),
* surfaces co-expressed EHEG pairs (Xenium tech-note rule #D),
* renders the result as a markdown + json report.

The core (``compute_expression_qc_from_matrix``) operates on a raw numpy
matrix so unit tests do not need anndata. The convenience entry point
``compute_expression_qc_from_anndata`` lazy-imports anndata and adapts.

See ``~/.claude/plans/composed-gliding-sparkle.md`` §Phase 3.
"""
from probe_designer.qc.expression import (
    ExpressionFlag,
    ExpressionQC,
    LOW_CUTOFF_DEFAULT,
    HIGH_CUTOFF_DEFAULT_V31,
    HIGH_CUTOFF_DEFAULT_V2,
    MARKER_MIN_DEFAULT_V31,
    MARKER_MIN_DEFAULT_V2,
    compute_expression_qc_from_matrix,
    compute_expression_qc_from_anndata,
)
from probe_designer.qc.report import (
    render_panel_qc_markdown,
    write_panel_qc_json,
)
from probe_designer.qc.viz import (
    IsoformMeta,
    ProbeAnnotation,
    plot_accessibility_tracks,
    plot_accessibility_heatmap,
    plot_rnafold_structure,
)


__all__ = [
    "ExpressionFlag", "ExpressionQC",
    "LOW_CUTOFF_DEFAULT",
    "HIGH_CUTOFF_DEFAULT_V31", "HIGH_CUTOFF_DEFAULT_V2",
    "MARKER_MIN_DEFAULT_V31", "MARKER_MIN_DEFAULT_V2",
    "compute_expression_qc_from_matrix",
    "compute_expression_qc_from_anndata",
    "render_panel_qc_markdown",
    "write_panel_qc_json",
    "IsoformMeta", "ProbeAnnotation",
    "plot_accessibility_tracks", "plot_accessibility_heatmap",
    "plot_rnafold_structure",
]
