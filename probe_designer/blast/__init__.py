"""Chemistry-aware BLAST helpers shared between mutation / TCR pipelines.

The mRNA pipeline keeps using ``probe_designer.filtering.SequenceFilter`` for
its species-aware specificity check; this subpackage implements the
gene-aware variant required by mutation panels (and, in the next round,
TCR panels).
"""
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
    hsp_spans_junction,
    is_same_gene_hit,
    parse_blast_results,
    run_batch_blast_search,
)

__all__ = [
    "DEFAULT_CLAMP_NT",
    "DEFAULT_IDENTITY_THRESHOLD",
    "DEFAULT_MAX_OFF_TARGET_BINDING",
    "DEFAULT_MIN_ALIGN_LENGTH",
    "GENE_SYNONYMS",
    "apply_gene_aware_filter",
    "binding_query_junction_index",
    "extract_binding_query",
    "find_best_per_mutation_record",
    "hsp_spans_junction",
    "is_same_gene_hit",
    "parse_blast_results",
    "run_batch_blast_search",
]
