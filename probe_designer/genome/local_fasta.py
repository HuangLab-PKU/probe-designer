"""pyfaidx-backed local genome accessor.

Builds the closure used by mutation / TCR designers to fetch genomic slices
without paying an Ensembl REST round-trip per query. Handles the chr-prefix
mismatch that comes up when the input mutation table uses ``chr3`` but the
indexed FASTA uses ``3`` (or vice versa).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Callable, Union

logger = logging.getLogger(__name__)

GenomeAccessor = Callable[[Union[str, int], int, int], str]


def build_genome_accessor(fasta_path: Union[str, Path]) -> GenomeAccessor:
    """Return ``accessor(chrom, start, end)`` that fetches 1-based inclusive slices.

    Parameters
    ----------
    fasta_path : str or Path
        Path to an indexed FASTA (``pyfaidx`` will create the ``.fai`` if missing).

    Returns
    -------
    accessor : callable
        ``accessor("chr3", 100, 200)`` returns the 101-nt substring at that
        locus as a string. Returns ``""`` and logs at DEBUG level on lookup
        failure (e.g. unknown chromosome) to keep batched callers running.
    """
    from pyfaidx import Fasta

    fasta_path = Path(fasta_path)
    genome = Fasta(str(fasta_path))
    fasta_keys = set(genome.keys())
    has_chr_prefix = any(k.startswith("chr") for k in fasta_keys)

    def accessor(chrom: Union[str, int], start: int, end: int) -> str:
        chrom_str = str(chrom)
        if has_chr_prefix:
            key = chrom_str if chrom_str.startswith("chr") else f"chr{chrom_str}"
        else:
            key = chrom_str[3:] if chrom_str.startswith("chr") else chrom_str
        try:
            return str(genome[key][start - 1:end])
        except Exception as exc:
            logger.debug("FASTA lookup failed for %s:%d-%d via key %r: %s",
                         chrom, start, end, key, exc)
            return ""

    return accessor
