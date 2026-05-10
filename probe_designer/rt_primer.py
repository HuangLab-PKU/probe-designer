"""Reverse-transcription primer designer.

Used by mutation cDNA chemistry today (and TCR cDNA next round). Designs a
primer that hybridizes to mRNA a configurable number of nt downstream of the
padlock target so RT extension produces the cDNA the padlock will bind.

The function searches primer lengths between ``min_primer_len`` and
``max_primer_len`` and returns the candidate whose ``R_DNA_NN1`` Tm falls
inside ``[tm_min, tm_max]`` and whose length is closest to
``target_primer_len``. Falls back to the closest-Tm candidate if none is in
range; the caller can decide whether to reject it.
"""
from __future__ import annotations

import math
from typing import Callable, Dict, Optional

from Bio.SeqUtils import MeltingTemp as mt


GenomeAccessor = Callable[[str, int, int], str]


def _reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp[b] for b in reversed(seq.upper()))


def _gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    return round(sum(1 for b in seq.upper() if b in "GC") / len(seq) * 100, 1)


def _compute_tm(seq: str) -> float:
    return round(mt.Tm_NN(seq, nn_table=mt.R_DNA_NN1), 1)


def _resolve_strand(strand) -> int:
    """Coerce strand to ``+1`` or ``-1``; default ``+1`` on None / NaN."""
    if strand is None:
        return 1
    try:
        if math.isnan(float(strand)):
            return 1
    except (ValueError, TypeError):
        return 1
    return int(strand)


def _primer_genomic_region(
    strand: int, mut_start: int, mut_end: int,
    probe_arm5_len: int, gap: int, primer_len: int,
) -> tuple[int, int]:
    """Return the (start, end) genomic coordinates the primer would bind.

    For strand=+1 the primer sits downstream of the padlock target on the
    plus strand; for strand=-1 it sits upstream (because target is on the
    minus strand).
    """
    if strand == 1:
        target_3prime_end = mut_end + probe_arm5_len
        primer_start = target_3prime_end + gap
        primer_end = primer_start + primer_len
    else:
        target_3prime_end = mut_start - probe_arm5_len
        primer_end = target_3prime_end - gap
        primer_start = primer_end - primer_len
    return primer_start, primer_end


def design_rt_primer(
    *,
    chrom: str,
    mut_start: int,
    mut_end: int,
    strand,
    probe_arm5_len: int,
    genome_accessor: GenomeAccessor,
    gap: int = 15,
    min_primer_len: int = 15,
    max_primer_len: int = 30,
    target_primer_len: int = 20,
    tm_min: float = 50.0,
    tm_max: float = 65.0,
) -> Dict:
    """Design one RT primer for the given mutation.

    Returns a dict with:

    ``primer_seq, primer_length, tm, gc_percent, primer_genomic_start,
    primer_genomic_end, gap_nt, strand, notes``

    On failure (no genomic sequence, no candidate found) returns an entry
    with ``primer_seq=""`` and a non-empty ``notes`` field; the caller can
    log / drop / retry as appropriate.
    """
    strand = _resolve_strand(strand)
    notes: list[str] = []

    # Fetch the widest possible genomic window (max primer length) once,
    # then slice locally for each candidate length.
    fetch_start_raw, fetch_end_raw = _primer_genomic_region(
        strand, mut_start, mut_end, probe_arm5_len, gap, max_primer_len,
    )
    fetch_start = min(fetch_start_raw, fetch_end_raw)
    fetch_end = max(fetch_start_raw, fetch_end_raw)
    genomic = genome_accessor(chrom, fetch_start, fetch_end)

    if not genomic or len(genomic) < min_primer_len:
        return {
            "primer_seq": "", "primer_length": 0,
            "tm": 0, "gc_percent": 0,
            "primer_genomic_start": fetch_start,
            "primer_genomic_end": fetch_end,
            "gap_nt": gap, "strand": strand,
            "notes": "FAILED: could not fetch genomic sequence",
        }

    best_in_range: Optional[tuple] = None
    best_fallback: Optional[tuple] = None

    for plen in range(min_primer_len, max_primer_len + 1):
        p_start, p_end = _primer_genomic_region(
            strand, mut_start, mut_end, probe_arm5_len, gap, plen,
        )
        f_start = min(p_start, p_end)
        f_end = max(p_start, p_end)
        offset = f_start - fetch_start
        sub = genomic[offset: offset + plen]
        if len(sub) < plen:
            continue
        primer = (_reverse_complement(sub.upper())
                  if strand == 1 else sub.upper())
        try:
            tm = _compute_tm(primer)
        except Exception:
            continue
        if tm_min <= tm <= tm_max:
            cand = (primer, tm, plen, f_start, f_end)
            if (best_in_range is None
                    or abs(plen - target_primer_len)
                    < abs(best_in_range[2] - target_primer_len)):
                best_in_range = cand
        else:
            dist = max(tm_min - tm, 0) + max(tm - tm_max, 0)
            fb_dist = (float("inf") if best_fallback is None
                       else max(tm_min - best_fallback[1], 0)
                       + max(best_fallback[1] - tm_max, 0))
            if (dist < fb_dist or
                (dist == fb_dist and best_fallback is not None and
                 abs(plen - target_primer_len)
                 < abs(best_fallback[2] - target_primer_len))):
                best_fallback = (primer, tm, plen, f_start, f_end)

    if best_in_range:
        primer, tm, plen, p_s, p_e = best_in_range
    elif best_fallback:
        primer, tm, plen, p_s, p_e = best_fallback
        if tm < tm_min:
            notes.append(f"Tm {tm}C below {tm_min}C (AT-rich region)")
        else:
            notes.append(f"Tm {tm}C above {tm_max}C (GC-rich region)")
    else:
        primer, tm, plen = "", 0, 0
        p_s, p_e = fetch_start, fetch_end
        notes.append("FAILED: no valid primer sequence")

    return {
        "primer_seq": primer, "primer_length": plen,
        "tm": tm, "gc_percent": _gc_percent(primer),
        "primer_genomic_start": p_s, "primer_genomic_end": p_e,
        "gap_nt": gap, "strand": strand,
        "notes": "; ".join(notes) if notes else "",
    }
