"""Reverse-transcription primer designer.

Two entry points:

- :func:`design_rt_primer` — mutation cDNA chemistry. Fetches genomic
  sequence via a callable (``GenomeAccessor``) using chromosomal coordinates
  and a strand. Primer = RC of the genomic slice (+strand) or literal
  (-strand).
- :func:`design_rt_primer_from_target` — TCR cDNA chemistry. Takes a target
  sequence directly (e.g. TRB mRNA) and a 0-indexed offset; primer is
  always RC of the target slice (TRB is mRNA-sense).

Both search primer lengths in ``[min_primer_len, max_primer_len]`` and
return the candidate whose ``R_DNA_NN1`` Tm falls inside ``[tm_min, tm_max]``
and whose length is closest to ``target_primer_len``. Falls back to the
closest-Tm candidate if none is in range; the caller decides whether to
reject it.
"""
from __future__ import annotations

import math
from typing import Callable, Dict, Optional

from Bio.SeqUtils import MeltingTemp as mt

from probe_designer.chemistry import ReactionConditions, dna_revcomp_to_rna


GenomeAccessor = Callable[[str, int, int], str]


def _reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp[b] for b in reversed(seq.upper()))


def _gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    return round(sum(1 for b in seq.upper() if b in "GC") / len(seq) * 100, 1)


def _compute_tm(primer: str, reaction: ReactionConditions) -> float:
    """Tm of a DNA RT primer annealed to its mRNA template.

    The primer:template duplex is DNA:RNA, so ``R_DNA_NN1`` must be fed the RNA
    (template) strand — the reverse complement of the primer — per Biopython's
    convention. Salt/Mg conditions come from ``reaction``. Reverse transcription
    carries no formamide, so the caller passes a formamide-free ReactionConditions.
    """
    rna_template = dna_revcomp_to_rna(primer)
    tm = mt.Tm_NN(rna_template, nn_table=mt.R_DNA_NN1, **reaction.tm_nn_kwargs())
    return round(reaction.apply_formamide(tm), 1)


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
    reaction: Optional[ReactionConditions] = None,
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
    # RT anneals in a formamide-free reverse-transcription buffer. The default
    # reproduces the legacy Biopython salt (Na=50, Mg=0, 25 nM) so the strand-
    # convention fix lands WITHOUT shifting the RT-primer Tm scale — the RT
    # window is not yet re-anchored to the real RT buffer (needs the RT reaction
    # temperature; a Phase-1 follow-up). Pass an explicit ``reaction`` to override.
    reaction = reaction or ReactionConditions(
        monovalent_mM=50.0, mg_mM=0.0, strand_nM=25.0, formamide_pct=0.0
    )
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
            tm = _compute_tm(primer, reaction)
        except ValueError:
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


def design_rt_primer_from_target(
    target_seq: str,
    target_end: int,
    *,
    gap: int = 15,
    min_primer_len: int = 15,
    max_primer_len: int = 30,
    target_primer_len: int = 20,
    tm_min: float = 50.0,
    tm_max: float = 65.0,
    reaction: Optional[ReactionConditions] = None,
) -> Dict:
    """Design one RT primer for a probe whose target ends at ``target_end``
    (0-indexed exclusive) on ``target_seq``.

    Used for TCR cDNA chemistry where the target sequence is the patient's
    TRB mRNA and no genomic accessor / strand math applies. The primer
    hybridizes to ``target_seq[target_end + gap : target_end + gap + plen]``
    and is the reverse complement of that slice.

    Returns a dict with the same keys as :func:`design_rt_primer` plus
    ``primer_target_start`` / ``primer_target_end`` (0-indexed positions on
    ``target_seq``). The ``primer_genomic_start`` / ``primer_genomic_end``
    keys carry the same target offsets (genomic semantics don't apply).
    """
    target_seq = target_seq.upper()
    # RT anneals in a formamide-free reverse-transcription buffer. The default
    # reproduces the legacy Biopython salt (Na=50, Mg=0, 25 nM) so the strand-
    # convention fix lands WITHOUT shifting the RT-primer Tm scale — the RT
    # window is not yet re-anchored to the real RT buffer (needs the RT reaction
    # temperature; a Phase-1 follow-up). Pass an explicit ``reaction`` to override.
    reaction = reaction or ReactionConditions(
        monovalent_mM=50.0, mg_mM=0.0, strand_nM=25.0, formamide_pct=0.0
    )
    notes: list[str] = []
    best_in_range: Optional[tuple] = None
    best_fallback: Optional[tuple] = None

    for plen in range(min_primer_len, max_primer_len + 1):
        start = target_end + gap
        end = start + plen
        if end > len(target_seq):
            continue
        region = target_seq[start:end]
        if len(region) < plen or "N" in region:
            continue
        primer = _reverse_complement(region)
        try:
            tm = _compute_tm(primer, reaction)
        except ValueError:
            continue
        if tm_min <= tm <= tm_max:
            cand = (primer, tm, plen, start, end)
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
                best_fallback = (primer, tm, plen, start, end)

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
        p_s, p_e = target_end + gap, target_end + gap + target_primer_len
        notes.append("FAILED: could not fit primer (target sequence too short or N-rich)")

    return {
        "primer_seq": primer, "primer_length": plen,
        "tm": tm, "gc_percent": _gc_percent(primer),
        "primer_target_start": p_s, "primer_target_end": p_e,
        "primer_genomic_start": p_s, "primer_genomic_end": p_e,
        "gap_nt": gap, "strand": 1,
        "notes": "; ".join(notes) if notes else "",
    }
