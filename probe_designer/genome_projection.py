"""Project a per-transcript-position profile onto genome coordinates.

A per-position annotation (arm Tm, accessibility, ...) lives in transcript
coordinates. To view it in a genome browser against the assembly, map each
transcript position to its genomic coordinate through the transcript's exon
structure and strand, and write a genome-coordinate bedGraph.

Isoform caveat: this is only well-defined for ONE transcript at a time (isoforms
disagree at shared exon positions and accessibility is per-isoform), so callers
should project the CANONICAL/representative transcript per gene, not all of them.
Junction caveat: for a windowed metric (e.g. arm Tm) the value is placed at the
window's ANCHOR position's genomic coordinate; a window that spans a splice
junction is annotated at its anchor, not split across the junction.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np


def _exon_bounds(exon: Any) -> Tuple[int, int]:
    """(start, end) genomic, Ensembl 1-based inclusive, from a dict or 2-tuple."""
    if isinstance(exon, dict):
        return int(exon["start"]), int(exon["end"])
    return int(exon[0]), int(exon[1])


def transcript_to_genome_positions(exons: List[Any], strand: int) -> List[int]:
    """Genomic (1-based) position for each transcript position, 5'->3'.

    Exons are concatenated in transcript order — ascending genomic order for the
    ``+`` strand, descending for the ``-`` strand (whose mRNA 5' end is the
    highest genomic coordinate).
    """
    ex = sorted((_exon_bounds(e) for e in exons), key=lambda t: t[0])
    t2g: List[int] = []
    if int(strand) < 0:
        for s, e in reversed(ex):
            t2g.extend(range(e, s - 1, -1))
    else:
        for s, e in ex:
            t2g.extend(range(s, e + 1))
    return t2g


def project_profile_to_genome(
    profile: np.ndarray, exons: List[Any], strand: int
) -> List[Tuple[int, float]]:
    """Map ``profile[i]`` (transcript position i) to ``(genomic_pos0, value)``.

    Returns 0-based genomic positions sorted ascending; NaN entries are skipped.
    If the profile is longer than the spliced length it is truncated (and vice
    versa), guarding against a transcript/exon length mismatch.
    """
    t2g = transcript_to_genome_positions(exons, strand)
    n = min(len(profile), len(t2g))
    out: List[Tuple[int, float]] = []
    for i in range(n):
        v = profile[i]
        if v == v:  # not NaN
            out.append((t2g[i] - 1, float(v)))  # Ensembl 1-based -> bedGraph 0-based
    out.sort(key=lambda p: p[0])
    return out


def write_genome_bedgraph(
    pairs: List[Tuple[int, float]],
    chrom: str,
    path: Path,
    *,
    track_name: str = "arm_tm",
) -> Path:
    """Write sorted ``(genomic_pos0, value)`` pairs as a genome bedGraph.

    Contiguous positions carrying an equal value are coalesced into one interval.
    """
    lines = [f'track type=bedGraph name="{track_name}"']
    i, n = 0, len(pairs)
    while i < n:
        pos, v = pairs[i]
        end = pos + 1
        j = i + 1
        while j < n and pairs[j][0] == end and pairs[j][1] == v:
            end += 1
            j += 1
        lines.append(f"{chrom}\t{pos}\t{end}\t{v:.2f}")
        i = j
    path = Path(path)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def pick_canonical_isoform(isoforms: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Choose the representative transcript to project to the genome.

    Prefers an explicit canonical flag (``is_canonical`` / ``canonical`` /
    ``mane_select``); otherwise the longest transcript (by total exon length).
    """
    for key in ("is_canonical", "canonical", "mane_select"):
        flagged = [iso for iso in isoforms if iso.get(key)]
        if flagged:
            return flagged[0]

    def total_exon_len(iso: Dict[str, Any]) -> int:
        return sum(e2 - e1 + 1 for e1, e2 in (_exon_bounds(x) for x in iso.get("Exon", [])))

    return max(isoforms, key=total_exon_len)


__all__ = [
    "transcript_to_genome_positions",
    "project_profile_to_genome",
    "write_genome_bedgraph",
    "pick_canonical_isoform",
]
