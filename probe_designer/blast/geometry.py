"""Ligation geometry read off a BLAST alignment.

Two properties decide whether a BLAST hit is a probe's real binding site, and
neither of them is a length.

**The nick.** A padlock only circularises where its two arms meet on a
contiguous stretch of target; SplintR / PBCV-1 ligase clamps ~3 nt on each side
of the nick and will not close over a mismatch there (Lohman 2014, Krzywkowski
2017 -- the same evidence behind ``DEFAULT_VICINITY_N``). A coverage ratio
cannot see this: a 39/40 alignment scores 0.975 whether the one unmatched base
is a trimmed terminal nucleotide (harmless) or the junction itself (dead
probe). ``Hsp_midline`` -- the ``|``/space line sitting in the same ``<Hsp>``
node the parser already walks -- resolves it, and lets the BLAST side apply the
criterion :mod:`probe_designer.isoform_coverage` already applies to exon
structure.

**The strand.** A hit on the strand the probe cannot template off is not an
off-target. Which strand that is follows from the chemistry, not from the
aligner:

* ``dRNA`` / ``iLock`` -- the arms are RC(mRNA) and bind the mRNA itself, so a
  real hit has the subject carrying RC(query): BLAST reports hit frame -1. A
  plus-strand hit means the transcript carries the probe's own sequence -- an
  antisense transcript -- which cannot base-pair with it.
* ``cDNA`` -- the arms are mRNA-sense and template off the cDNA, so the hit
  that matters is the transcript carrying the query verbatim: frame +1.

Measured on the SP369 marker panel: only 30 of 4299 full-length identities sat
on the non-binding strand, but counting them moved the non-specific set from 40
to 61, and four probes (CD27, ACTA2, KRT7, NCAM1) were condemned purely by their
own gene's antisense transcript. It is easy to miss because an antisense locus
really is a different gene, so the same-gene test says "different" and the
rejection looks correct. That audit resolved one clean gene symbol per hit;
``gene_aware.is_same_gene_hit`` instead matches the symbol as a whole word, which
already succeeds inside ``KRT7-AS`` (the dash is a non-word character), so what
the strand test adds here is antisense transcripts of OTHER loci.
"""
from __future__ import annotations

from typing import Optional

from probe_designer.qc.cross_ligation import DEFAULT_VICINITY_N

__all__ = [
    "DEFAULT_VICINITY_N",
    "binding_hit_frame",
    "frame_is_bindable",
    "hit_is_bindable",
    "hsp_junction_verdict",
    "hsp_spans_junction",
    "junction_offset",
    "midpoint_junction_offset",
]


# Chemistries whose arms are mRNA-sense: the probe templates off the cDNA, so
# the transcript that matters is the one carrying the query verbatim.
_SENSE_CHEMISTRIES = frozenset({"cDNA"})


def binding_hit_frame(chemistry: str) -> int:
    """The ``Hsp_hit-frame`` sign a hit must carry to be bindable by ``chemistry``."""
    return 1 if chemistry in _SENSE_CHEMISTRIES else -1


def frame_is_bindable(hit_frame: Optional[int], binding_frame: int) -> bool:
    """``False`` only when BLAST positively places the hit on the strand the
    probe cannot bind.

    An unreported or unexpected frame stays bindable: a missing measurement must
    never be the thing that lets a probe through.
    """
    if hit_frame not in (1, -1):
        return True
    return hit_frame == binding_frame


def hit_is_bindable(hit_frame: Optional[int], chemistry: str) -> bool:
    """:func:`frame_is_bindable` with the expected frame taken from chemistry."""
    return frame_is_bindable(hit_frame, binding_hit_frame(chemistry))


def junction_offset(arm3: str, chemistry: str) -> int:
    """Query index of the first ligated base 3' of the nick, for an ``arm3 + arm5`` query.

    ``len(arm3)`` for every chemistry. The iLock query
    (:func:`~probe_designer.blast.gene_aware.extract_binding_query`) drops
    ``arm3[-1]``, but ``arm5[0] == arm3[-1]``, so it still reads
    ``arm3 + arm5[1:]`` -- which is the ligated product: FEN1 cleaves the flap
    *plus* ``arm5[0]``, leaving ``arm3[-1]`` as the 3'-OH and ``arm5[1]`` as the
    5'-P (:func:`probe_designer.qc.cross_ligation.ligation_geometry`). A nick at
    ``len(arm3) - 1`` is the out-of-phase iLock geometry of cross-lig v2.2.
    ``chemistry`` stays in the signature so callers need not special-case iLock.
    """
    return len(arm3 or "")


def midpoint_junction_offset(probe_length: int) -> int:
    """Nick position for a probe whose arms have not been split out yet.

    ``arm_3prime`` is the first ``length // 2`` nt of the probe -- the same
    convention :func:`probe_designer.isoform_coverage._junction_index` applies
    to the target window.
    """
    return int(probe_length) // 2


def hsp_junction_verdict(
    midline: str,
    query_seq: str,
    query_start: int,
    junction_offset_nt: int,
    *,
    vicinity_n: int = DEFAULT_VICINITY_N,
) -> Optional[bool]:
    """Is the nick vicinity contiguously Watson-Crick paired in this HSP?

    ``midline`` and ``query_seq`` are ``Hsp_midline`` / ``Hsp_qseq``
    (Biopython: ``hsp.match`` / ``hsp.query``); ``query_start`` is the 1-based
    ``Hsp_query-from`` (Biopython ``hsp.query_start``); ``junction_offset_nt``
    is the 0-based query index of the first base 3' of the nick.

    ``vicinity_n`` bases each side of the nick must be present in the alignment
    and paired, with no gap in either sequence between them -- every column in
    the span has to be ``|``, which rejects a mismatch, a deletion and an
    insertion alike. ``0`` requires only the two bases flanking the nick.

    Returns ``None`` when the alignment carries no midline to judge from. The
    caller decides which way "no evidence" falls; both callers make it fall
    against the probe, in the direction that is conservative for their question.
    """
    if not midline or not query_seq or len(midline) != len(query_seq):
        return None
    span = max(int(vicinity_n), 0) or 1
    lo = int(junction_offset_nt) - span
    hi = int(junction_offset_nt) + span - 1
    if lo < 0:
        # The clamp reaches past the start of the query: nothing to measure,
        # which is not the same as a broken nick. Say so, and let each caller
        # fall the conservative way for its own question.
        return None

    # Walk the alignment columns, skipping the ones where the query carries a
    # gap, to find the columns holding query positions lo and hi.
    col_lo: Optional[int] = None
    col_hi: Optional[int] = None
    pos = int(query_start) - 1
    for col, base in enumerate(query_seq):
        if base == "-":
            continue
        if pos == lo:
            col_lo = col
        if pos == hi:
            col_hi = col
            break
        pos += 1
    if col_lo is None or col_hi is None:
        return False
    return all(char == "|" for char in midline[col_lo:col_hi + 1])


def hsp_spans_junction(
    midline: str,
    query_seq: str,
    query_start: int,
    junction_offset_nt: int,
    *,
    vicinity_n: int = DEFAULT_VICINITY_N,
) -> bool:
    """Strict form of :func:`hsp_junction_verdict`: no midline is not a pass."""
    return hsp_junction_verdict(
        midline, query_seq, query_start, junction_offset_nt,
        vicinity_n=vicinity_n,
    ) is True
