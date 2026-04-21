"""Default isoform provider: local GTF -> Ensembl REST chain.

Webapp wraps this with a DB cache (not included here to keep designer
library DB-agnostic).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from probe_designer.genome.ensembl_client import fetch_ensembl_isoforms
from probe_designer.genome.gtf_parser import parse_gtf_for_gene


logger = logging.getLogger(__name__)


class DefaultIsoformProvider:
    """Resolves per-gene isoforms from a local GTF first, Ensembl REST second.

    Args:
        gtf_path: if provided, local GTF is tried first.
        rest_base: override Ensembl REST endpoint (useful for testing).
    """

    def __init__(
        self,
        gtf_path: Optional[str | Path] = None,
        *,
        rest_base: str = "https://rest.ensembl.org",
    ) -> None:
        self.gtf_path = Path(gtf_path) if gtf_path else None
        self.rest_base = rest_base

    def get_isoforms(
        self, gene: str, species: str
    ) -> Optional[List[Dict[str, Any]]]:
        if self.gtf_path and self.gtf_path.exists():
            data = parse_gtf_for_gene(self.gtf_path, gene)
            isoforms = data.get("isoforms") or []
            if isoforms:
                return isoforms

        data = fetch_ensembl_isoforms(gene, species, rest_base=self.rest_base)
        isoforms = data.get("isoforms") or []
        return isoforms or None
