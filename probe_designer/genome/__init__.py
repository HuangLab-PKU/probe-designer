"""Genome reference and isoform data sources.

Lives in the designer library so both CLI and webapp share:
  - :func:`parse_gtf_for_gene` — streamed single-gene extraction from GTF
  - :func:`fetch_ensembl_isoforms` — Ensembl REST symbol -> exons
  - :class:`DefaultIsoformProvider` — local GTF first, Ensembl fallback
"""
from probe_designer.genome.ensembl_client import fetch_ensembl_isoforms
from probe_designer.genome.gtf_parser import parse_gtf_for_gene
from probe_designer.genome.isoform_provider import DefaultIsoformProvider


__all__ = [
    "parse_gtf_for_gene",
    "fetch_ensembl_isoforms",
    "DefaultIsoformProvider",
]
