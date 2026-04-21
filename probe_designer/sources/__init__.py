"""Sequence data sources (Ensembl, NCBI, local FASTA).

Phase 4 introduces :mod:`credentials` here; Phase 5 will move the rest of
``database.py`` into this subpackage (``ensembl.py``, ``ncbi.py``, etc.).
"""
from probe_designer.sources.credentials import get_entrez_credentials


__all__ = ["get_entrez_credentials"]
