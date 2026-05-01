"""Streamed GTF extraction for a single gene.

Migrated from webapp/services/genome.py::_parse_gtf_for_gene. Adds an
optional ``<gtf>.pidx.json`` byte-offset index so many-gene runs don't
re-scan the full GTF per gene.
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional


logger = logging.getLogger(__name__)


def parse_gtf_for_gene(gtf_path: str | os.PathLike, gene_symbol: str) -> Dict[str, Any]:
    """Return isoform + exon structure for ``gene_symbol`` from ``gtf_path``.

    The return shape matches webapp's previous contract:
        {
            "chromosome": str,
            "start": int, "end": int,
            "strand": 1 | -1,
            "isoforms": [
                {"accession": str, "biotype": str,
                 "start": int, "end": int,
                 "exons": [{"start", "end", "rel_start", "rel_width"}]},
                ...
            ],
        }

    On any parse error the function returns ``{"isoforms": []}`` (non-fatal).
    """
    gtf_path = Path(gtf_path)
    if not gtf_path.exists():
        return {"isoforms": []}

    index = _load_or_build_index(gtf_path)

    transcripts: Dict[str, Dict[str, Any]] = {}
    gene_start: Optional[int] = None
    gene_end: Optional[int] = None
    chrom = ""
    strand = "+"

    try:
        with open(gtf_path, "r", encoding="utf-8") as fp:
            offsets = index.get(gene_symbol) if index else None
            if offsets is not None:
                # Fast path: seek to known byte offsets
                for offset in offsets:
                    fp.seek(offset)
                    line = fp.readline()
                    if not line or f'gene_name "{gene_symbol}"' not in line:
                        continue
                    gene_start, gene_end, chrom, strand = _absorb_line(
                        line, transcripts, gene_symbol, gene_start, gene_end, chrom, strand
                    )
            else:
                # Full scan fallback
                for line in fp:
                    if line.startswith("#"):
                        continue
                    if f'gene_name "{gene_symbol}"' not in line:
                        continue
                    gene_start, gene_end, chrom, strand = _absorb_line(
                        line, transcripts, gene_symbol, gene_start, gene_end, chrom, strand
                    )
    except (OSError, ValueError) as exc:
        logger.warning("GTF parse failed for %s: %s", gene_symbol, exc)
        return {"isoforms": []}

    if not transcripts:
        return {"isoforms": []}

    if gene_start is None:
        gene_start = min(t["start"] for t in transcripts.values())
        gene_end = max(t["end"] for t in transcripts.values())

    gene_length = max(1, (gene_end or gene_start) - gene_start)

    strand_int = 1 if strand == "+" else -1

    isoforms = []
    for tx_id, tx_data in transcripts.items():
        exons = []
        exon_legacy = []
        for ex in tx_data["exons"]:
            rel_start = max(0.0, (ex["start"] - gene_start) / gene_length * 100)
            rel_width = max(0.5, (ex["end"] - ex["start"]) / gene_length * 100)
            exons.append({
                "start": ex["start"],
                "end": ex["end"],
                "rel_start": round(rel_start, 2),
                "rel_width": round(rel_width, 2),
            })
            exon_legacy.append({
                "start": ex["start"],
                "end": ex["end"],
                "strand": strand_int,
            })
        isoforms.append({
            # v0.2.1 shape
            "accession": tx_id,
            "biotype": tx_data.get("biotype", ""),
            "start": tx_data["start"],
            "end": tx_data["end"],
            "exons": sorted(exons, key=lambda e: e["start"]),
            # Legacy Ensembl-Transcript shape (search_strategies.IsoformAwareness)
            "id": tx_id,
            "display_name": tx_id,
            "seq_region_name": chrom,
            "strand": strand_int,
            "Exon": sorted(exon_legacy, key=lambda e: e["start"]),
        })

    return {
        "chromosome": chrom,
        "start": gene_start,
        "end": gene_end,
        "strand": strand_int,
        "isoforms": isoforms,
    }


def _absorb_line(
    line: str,
    transcripts: Dict[str, Dict[str, Any]],
    gene_symbol: str,
    gene_start: Optional[int],
    gene_end: Optional[int],
    chrom: str,
    strand: str,
):
    parts = line.strip().split("\t")
    if len(parts) < 9:
        return gene_start, gene_end, chrom, strand

    feature = parts[2]
    start = int(parts[3])
    end = int(parts[4])
    chrom = parts[0]
    strand = parts[6]
    attrs = parts[8]

    tx_id = _extract_attr(attrs, "transcript_id")
    biotype = _extract_attr(attrs, "transcript_biotype") or ""

    if feature == "gene":
        gene_start = start
        gene_end = end
    elif feature == "exon" and tx_id:
        tx = transcripts.setdefault(
            tx_id,
            {"accession": tx_id, "exons": [], "start": start, "end": end, "biotype": biotype},
        )
        tx["exons"].append({"start": start, "end": end})
        tx["start"] = min(tx["start"], start)
        tx["end"] = max(tx["end"], end)
        if biotype and not tx.get("biotype"):
            tx["biotype"] = biotype
    elif feature == "transcript" and tx_id:
        tx = transcripts.setdefault(
            tx_id,
            {"accession": tx_id, "exons": [], "start": start, "end": end, "biotype": biotype},
        )
        tx["biotype"] = biotype or tx.get("biotype", "")
        tx["start"] = min(tx["start"], start)
        tx["end"] = max(tx["end"], end)

    return gene_start, gene_end, chrom, strand


def _extract_attr(attrs: str, key: str) -> Optional[str]:
    for part in attrs.split(";"):
        part = part.strip()
        if part.startswith(f'{key} "'):
            return part.split('"')[1]
    return None


# --- Byte-offset index -----------------------------------------------------
#
# Multi-gene CLI runs call parse_gtf_for_gene() many times. A full file scan
# per call costs O(n_genes × file_size). We maintain an auxiliary JSON
# alongside the GTF mapping gene_symbol → list-of-line-start-offsets, so
# subsequent calls can seek directly to matching lines.
#
# The index is rebuilt whenever the GTF mtime changes.

def _index_path_for(gtf_path: Path) -> Path:
    return gtf_path.with_suffix(gtf_path.suffix + ".pidx.json")


def _load_or_build_index(gtf_path: Path) -> Optional[Dict[str, List[int]]]:
    idx_path = _index_path_for(gtf_path)
    try:
        gtf_mtime = gtf_path.stat().st_mtime
    except OSError:
        return None

    if idx_path.exists():
        try:
            idx_mtime = idx_path.stat().st_mtime
            if idx_mtime >= gtf_mtime:
                with open(idx_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                if data.get("gtf_mtime") == gtf_mtime:
                    return data.get("offsets", {})
        except (OSError, ValueError, json.JSONDecodeError):
            pass  # fall through to rebuild

    try:
        offsets = _build_index(gtf_path)
    except OSError as exc:
        logger.warning("GTF index build failed: %s", exc)
        return None

    try:
        with open(idx_path, "w", encoding="utf-8") as f:
            json.dump({"gtf_mtime": gtf_mtime, "offsets": offsets}, f)
    except OSError as exc:
        logger.warning("Could not write GTF index %s: %s", idx_path, exc)

    return offsets


def _build_index(gtf_path: Path) -> Dict[str, List[int]]:
    """Scan the GTF once, recording byte offset of every line per gene_name."""
    offsets: Dict[str, List[int]] = {}
    with open(gtf_path, "rb") as fp:
        while True:
            line_start = fp.tell()
            raw = fp.readline()
            if not raw:
                break
            if raw.startswith(b"#"):
                continue
            # Scan only for gene_name marker; avoid decoding full line
            marker = raw.find(b'gene_name "')
            if marker < 0:
                continue
            end = raw.find(b'"', marker + len(b'gene_name "'))
            if end < 0:
                continue
            gene = raw[marker + len(b'gene_name "'): end].decode("utf-8", errors="replace")
            offsets.setdefault(gene, []).append(line_start)
    return offsets
