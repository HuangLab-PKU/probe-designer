"""Default isoform provider: local GTF -> Ensembl REST chain.

Webapp wraps this with a DB cache (not included here to keep designer
library DB-agnostic).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set

from probe_designer.genome.ensembl_client import fetch_ensembl_isoforms
from probe_designer.genome.gtf_parser import parse_gtf_for_gene


logger = logging.getLogger(__name__)


class DefaultIsoformProvider:
    """Resolves per-gene isoforms from a local GTF first, Ensembl REST second.

    Args:
        gtf_path: if provided, local GTF is tried first.
        rest_base: override Ensembl REST endpoint (useful for testing).
        keep_biotypes: if non-empty, drop isoforms whose ``biotype`` is not
            in this set before returning. Useful for ``isoform_consensus``,
            where retained_intron / NMD / processed_transcript isoforms
            otherwise dilute consensus across the productive protein-coding
            transcripts. ``None`` or empty == no filter (keep everything).
            If the filter would drop ALL isoforms for a gene, the unfiltered
            list is returned and a warning is logged.
    """

    def __init__(
        self,
        gtf_path: Optional[str | Path] = None,
        *,
        rest_base: str = "https://rest.ensembl.org",
        keep_biotypes: Optional[Iterable[str]] = None,
    ) -> None:
        self.gtf_path = Path(gtf_path) if gtf_path else None
        self.rest_base = rest_base
        self.keep_biotypes: Set[str] = set(keep_biotypes) if keep_biotypes else set()

    def get_isoforms(
        self, gene: str, species: str
    ) -> Optional[List[Dict[str, Any]]]:
        if self.gtf_path and self.gtf_path.exists():
            data = parse_gtf_for_gene(self.gtf_path, gene)
            isoforms = data.get("isoforms") or []
        else:
            data = fetch_ensembl_isoforms(gene, species, rest_base=self.rest_base)
            isoforms = data.get("isoforms") or []

        return self._apply_biotype_filter(gene, isoforms) if isoforms else None

    def _apply_biotype_filter(
        self, gene: str, isoforms: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        if not self.keep_biotypes:
            return isoforms
        kept = [iso for iso in isoforms if iso.get("biotype") in self.keep_biotypes]
        if not kept:
            seen = sorted({str(iso.get("biotype", "")) for iso in isoforms})
            logger.warning(
                "[%s] biotype filter %s dropped all %d isoforms (present: %s); "
                "falling back to unfiltered list",
                gene, sorted(self.keep_biotypes), len(isoforms), seen,
            )
            return isoforms
        if len(kept) < len(isoforms):
            logger.info(
                "[%s] biotype filter: %d/%d isoforms kept",
                gene, len(kept), len(isoforms),
            )
        return kept
