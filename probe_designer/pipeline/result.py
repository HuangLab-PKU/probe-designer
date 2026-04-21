"""In-memory pipeline result containers.

Both CLI and webapp receive :class:`PipelineResult` from ``Pipeline.run()``.
The CLI serializes it to disk via :mod:`probe_designer.io.writers`; the
webapp ingests it directly into the database without a JSON round trip.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass
class GeneResult:
    """Per-gene pipeline output."""

    gene: str
    sites: List[Dict[str, Any]] = field(default_factory=list)
    existing_targets: List[Dict[str, Any]] = field(default_factory=list)
    sequence_source: Optional[str] = None  # "ensembl" | "ncbi" | "local" | None
    sequence_length: Optional[int] = None
    isoform_count: Optional[int] = None
    errors: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "gene": self.gene,
            "sites": self.sites,
            "existing_targets": self.existing_targets,
            "sequence_source": self.sequence_source,
            "sequence_length": self.sequence_length,
            "isoform_count": self.isoform_count,
            "errors": self.errors,
        }


@dataclass
class PipelineResult:
    """Whole-run pipeline output."""

    per_gene: Dict[str, GeneResult] = field(default_factory=dict)
    config_snapshot: Dict[str, Any] = field(default_factory=dict)
    run_id: str = ""
    output_dir: Optional[Path] = None

    @property
    def gene_count(self) -> int:
        return len(self.per_gene)

    @property
    def total_sites(self) -> int:
        return sum(len(r.sites) for r in self.per_gene.values())

    def sites_by_gene(self) -> Dict[str, List[Dict[str, Any]]]:
        """Flattened view matching the old dict-of-lists JSON format."""
        return {gene: result.sites for gene, result in self.per_gene.items()}
