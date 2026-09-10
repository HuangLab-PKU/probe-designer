"""Chemistry-aware BLAST helpers shared between mutation / TCR pipelines.

The mRNA pipeline keeps using ``probe_designer.filtering.SequenceFilter`` for
its species-aware specificity check; this subpackage implements the
gene-aware variant required by mutation panels (and, in the next round,
TCR panels).
"""
from probe_designer.blast.geometry import (
    DEFAULT_VICINITY_N,
    binding_hit_frame,
    frame_is_bindable,
    hit_is_bindable,
    hsp_junction_verdict,
    hsp_spans_junction,
    junction_offset,
    midpoint_junction_offset,
)
from probe_designer.blast.gene_aware import (
    DEFAULT_CLAMP_NT,
    DEFAULT_IDENTITY_THRESHOLD,
    DEFAULT_MAX_OFF_TARGET_BINDING,
    DEFAULT_MIN_ALIGN_LENGTH,
    GENE_SYNONYMS,
    apply_gene_aware_filter,
    binding_query_junction_index,
    extract_binding_query,
    find_best_per_mutation_record,
    hit_can_ligate,
    hsp_covers_junction,
    is_same_gene_hit,
    parse_blast_results,
    run_batch_blast_search,
)

__all__ = [
    "DEFAULT_CLAMP_NT",
    "DEFAULT_IDENTITY_THRESHOLD",
    "DEFAULT_MAX_OFF_TARGET_BINDING",
    "DEFAULT_MIN_ALIGN_LENGTH",
    "DEFAULT_VICINITY_N",
    "GENE_SYNONYMS",
    "apply_gene_aware_filter",
    "binding_hit_frame",
    "binding_query_junction_index",
    "extract_binding_query",
    "find_best_per_mutation_record",
    "frame_is_bindable",
    "hit_can_ligate",
    "hit_is_bindable",
    "hsp_covers_junction",
    "hsp_junction_verdict",
    "hsp_spans_junction",
    "is_same_gene_hit",
    "junction_offset",
    "midpoint_junction_offset",
    "parse_blast_results",
    "run_batch_blast_search",
]
