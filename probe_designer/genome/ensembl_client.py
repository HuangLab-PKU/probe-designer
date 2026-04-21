"""Ensembl REST client for isoform/exon lookup.

Migrated from webapp/services/genome.py::_fetch_from_ensembl (no DB access,
pure HTTP). Both CLI and webapp adapters import this.
"""
from __future__ import annotations

import logging
from typing import Any, Dict

import requests


logger = logging.getLogger(__name__)


_SPECIES_TO_ENSEMBL = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
    "rat": "rattus_norvegicus",
    "zebrafish": "danio_rerio",
}


def fetch_ensembl_isoforms(
    gene_symbol: str,
    species: str,
    *,
    rest_base: str = "https://rest.ensembl.org",
    timeout: float = 15.0,
) -> Dict[str, Any]:
    """Fetch isoform and exon metadata from Ensembl REST.

    Returns a dict with keys ``gene_id_ensembl``, ``chromosome``, ``start``,
    ``end``, ``strand``, ``isoforms``. On any HTTP or parse failure returns
    ``{"isoforms": [], "message": ...}`` (non-fatal).
    """
    ensembl_species = _SPECIES_TO_ENSEMBL.get(species, species)

    try:
        resp = requests.get(
            f"{rest_base}/lookup/symbol/{ensembl_species}/{gene_symbol}",
            headers={"Content-Type": "application/json"},
            timeout=timeout,
        )
        if resp.status_code != 200:
            return {"isoforms": [], "message": f"Ensembl lookup failed: {resp.status_code}"}

        gene_data = resp.json()
        gene_id_ensembl = gene_data.get("id")
        gene_start = gene_data.get("start", 0)
        gene_end = gene_data.get("end", 0)
        chrom = gene_data.get("seq_region_name", "")
        strand = gene_data.get("strand", 1)

        resp2 = requests.get(
            f"{rest_base}/lookup/id/{gene_id_ensembl}?expand=1",
            headers={"Content-Type": "application/json"},
            timeout=timeout,
        )
        if resp2.status_code != 200:
            return {"isoforms": [], "message": "Could not expand transcripts"}

        gene_full = resp2.json()
        transcripts = gene_full.get("Transcript", [])

        gene_length = gene_end - gene_start
        if gene_length <= 0:
            gene_length = 1

        isoforms = []
        for tx in transcripts:
            exons = []
            for exon in tx.get("Exon", []):
                ex_start = exon.get("start", 0)
                ex_end = exon.get("end", 0)
                rel_start = max(0, (ex_start - gene_start) / gene_length * 100)
                rel_width = max(0.5, (ex_end - ex_start) / gene_length * 100)
                exons.append({
                    "start": ex_start,
                    "end": ex_end,
                    "rel_start": round(rel_start, 2),
                    "rel_width": round(rel_width, 2),
                })
            isoforms.append({
                "accession": tx.get("id", ""),
                "biotype": tx.get("biotype", ""),
                "start": tx.get("start", 0),
                "end": tx.get("end", 0),
                "exons": exons,
            })

        return {
            "gene_id_ensembl": gene_id_ensembl,
            "chromosome": chrom,
            "start": gene_start,
            "end": gene_end,
            "strand": strand,
            "isoforms": isoforms,
        }

    except requests.RequestException as exc:
        return {"isoforms": [], "message": f"Ensembl fetch error: {exc}"}
    except (ValueError, KeyError) as exc:
        return {"isoforms": [], "message": f"Ensembl parse error: {exc}"}
