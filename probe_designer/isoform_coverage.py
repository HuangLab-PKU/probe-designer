"""Which isoforms a padlock binding window can actually be ligated on.

``isoform_consensus`` needs one number per candidate window: how many of the
gene's transcripts the probe will bind. It used to answer that with
``IsoformAwareness.isoforms_union`` — the set of isoforms owning any exon
fragment the window OVERLAPS. That is not the same question. An isoform can
own part of the window and still splice the rest away, and because the count is
also a ranking key, over-crediting actively PREFERRED such windows. Six sites
in the 2026-07-27 CRC panel bind no protein-coding transcript of their own gene
as a result (``DESIGNER_ISSUE_isoform_coverage.md``).

The condition a transcript has to meet is physical, not positional. The probe
is one contiguous oligo whose two arms meet at the ligation junction; SplintR /
Ampligase read that junction and will not close a nick that is not properly
paired on both sides. So a transcript is bound only when

    a *contiguous run of its own mRNA* covers the junction, with at least
    ``min_arm`` nt of the window on each side of it.

Contiguity, not mere presence, is the operative word: a retained-intron
transcript can contain every base of a junction-spanning window and still not
bind it, because in that mRNA the two halves are separated by the intron.
Shortfall at the far ends of the arms is tolerated down to ``min_arm`` — the
lab's stated rule is a junction-centred core of about +/-16 nt — because those
bases are distal to the nick and cost affinity rather than fidelity.

Coordinates are Ensembl-style: 1-based, both ends inclusive.
"""
from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple


Segment = Tuple[int, int]

PROTEIN_CODING = "protein_coding"

# Transcript biotypes that make a protein. `protein_coding` alone is too narrow:
# immune-receptor genes carry their own Ensembl biotypes and NEVER appear as
# protein_coding — TRAC and TRBC1 are TR_C_gene, IGHG1 is IG_C_gene — so a
# protein_coding-only rule silently rejects every valid TCR/BCR constant-region
# probe. Found by the 2026-08-02 bank audit, which flagged six good TRAC/TRBC1
# sites as binding nothing.
#
# Listed explicitly rather than matched by an `IG_`/`TR_` prefix, which would
# also admit the several hundred IG_*_pseudogene / TR_*_pseudogene entries.
# `protein_coding_CDS_not_defined` is likewise excluded: it is what newer
# GENCODE releases call transcripts that used to be processed_transcript, and it
# has no CDS.
PRODUCTIVE_BIOTYPES = frozenset({
    PROTEIN_CODING,
    "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
    "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
})


def is_productive(biotype: str) -> bool:
    """True when this transcript biotype makes a protein."""
    return biotype in PRODUCTIVE_BIOTYPES


@dataclass(frozen=True)
class WindowCoverage:
    """How much of one binding window survives in one isoform's mRNA.

    ``lower_nt`` / ``upper_nt`` count the window bases contiguously paired on
    the genomic-lower and genomic-upper side of the ligation junction; the arm
    fields are the same two numbers mapped through the strand onto the padlock's
    arms. All four are zero when the junction itself is not paired, since the
    probe cannot be ligated on that transcript at all — ``ligatable`` says which
    of the two zero cases you are looking at.
    """

    ligatable: bool
    lower_nt: int
    upper_nt: int
    arm_3prime_nt: int
    arm_5prime_nt: int

    def credited(self, min_arm: int) -> bool:
        """True when this isoform counts as targeted at the given bound."""
        return (
            self.ligatable
            and self.lower_nt >= min_arm
            and self.upper_nt >= min_arm
        )


_UNCOVERED = WindowCoverage(
    ligatable=False, lower_nt=0, upper_nt=0, arm_3prime_nt=0, arm_5prime_nt=0
)


def _merge_overlaps(exons: List[Segment]) -> List[Segment]:
    """Collapse exons that overlap; leave abutting ones alone.

    A well-formed transcript has neither, but an overlap would make the exon-end
    list non-monotonic and silently break the bisect. Abutting exons (a
    zero-length intron) are kept separate because ``splices()`` still has to see
    that junction.
    """
    merged: List[Segment] = []
    for start, end in exons:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


class _ExonTrack:
    """One isoform's exons, prepared for repeated window queries."""

    __slots__ = ("name", "biotype", "exons", "_ends", "_acceptor_of")

    def __init__(self, isoform: Dict[str, Any]) -> None:
        self.name = isoform.get("display_name") or isoform.get("id") or ""
        self.biotype = isoform.get("biotype") or ""
        self.exons: List[Segment] = _merge_overlaps(
            sorted((int(e["start"]), int(e["end"])) for e in isoform.get("Exon", []))
        )
        # Ascending because _merge_overlaps leaves no exon containing another —
        # the bisect in pieces() depends on it.
        self._ends = [end for _, end in self.exons]
        # donor -> acceptor of the very next exon; the only way two window
        # halves separated by a host intron can be adjacent in this mRNA.
        self._acceptor_of = {
            end: self.exons[i + 1][0] for i, (_, end) in enumerate(self.exons[:-1])
        }

    def pieces(self, start: int, end: int) -> List[Segment]:
        """Parts of ``[start, end]`` present in this isoform, ascending.

        Exons that abut in the genome are merged, since their bases are
        adjacent in the mRNA too.
        """
        out: List[Segment] = []
        for i in range(bisect_left(self._ends, start), len(self.exons)):
            ex_start, ex_end = self.exons[i]
            if ex_start > end:
                break
            lo, hi = max(ex_start, start), min(ex_end, end)
            if lo > hi:
                continue
            if out and lo <= out[-1][1] + 1:
                out[-1] = (out[-1][0], max(hi, out[-1][1]))
            else:
                out.append((lo, hi))
        return out

    def splices(self, donor: int, acceptor: int) -> bool:
        """True when this isoform joins ``donor`` to ``acceptor`` exactly."""
        return self._acceptor_of.get(donor) == acceptor


def _junction_index(total: int, strand: int) -> int:
    """Window index (genomic-ascending, 0-based) the ligation junction sits at.

    ``arm_3prime`` is the first ``total // 2`` nt of the probe and the probe is
    the reverse complement of the window, so arm3 pairs with the LAST
    ``total // 2`` bases of the target. On the plus strand the target reads
    along ascending genomic coordinates and that is the upper side; on the minus
    strand it is the lower one.
    """
    return total - total // 2 if strand >= 0 else total // 2


def window_coverage(
    segments: Sequence[Segment],
    exons: Sequence[Segment] | Sequence[Dict[str, Any]],
    strand: int,
) -> WindowCoverage:
    """Coverage of one window by one isoform's exon structure.

    ``segments`` are the window's genomic pieces in ascending order — the
    output of ``IsoformAwareness.find_target_region``, contiguous in the HOST
    isoform's mRNA by construction. ``exons`` is the exon list of the isoform
    being asked about, as ``(start, end)`` pairs or Ensembl exon dicts.
    """
    normalized: List[Segment] = sorted(
        (int(s), int(e)) for s, e in
        ((x["start"], x["end"]) if isinstance(x, dict) else x for x in exons)
    )
    track = _ExonTrack({"Exon": [{"start": s, "end": e} for s, e in normalized]})
    return _coverage(list(segments), track, strand)


def _coverage(
    segments: Sequence[Segment], track: _ExonTrack, strand: int,
) -> WindowCoverage:
    segs: List[Segment] = sorted((int(s), int(e)) for s, e in segments)
    total = sum(end - start + 1 for start, end in segs)
    if total < 2 or not track.exons:
        return _UNCOVERED

    junction = _junction_index(total, strand)

    # Walk the window in genomic-ascending order, accumulating maximal runs of
    # window positions that are contiguous in THIS isoform's mRNA.
    runs: List[List[int]] = []
    offset = 0
    prev_end: Optional[int] = None
    prev_reached_end = False
    for index, (seg_start, seg_end) in enumerate(segs):
        parts = track.pieces(seg_start, seg_end)
        for part_index, (lo, hi) in enumerate(parts):
            win_start = offset + (lo - seg_start)
            win_end = offset + (hi - seg_start) + 1
            joins_previous = (
                part_index == 0
                and index > 0
                and lo == seg_start
                and prev_reached_end
                and track.splices(prev_end, seg_start)
            )
            if joins_previous and runs and runs[-1][1] == win_start:
                runs[-1][1] = win_end
            else:
                runs.append([win_start, win_end])
        prev_reached_end = bool(parts) and parts[-1][1] == seg_end
        prev_end = seg_end
        offset += seg_end - seg_start + 1

    for run_start, run_end in runs:
        if run_start < junction < run_end:
            lower = junction - run_start
            upper = run_end - junction
            arm_3prime, arm_5prime = (
                (upper, lower) if strand >= 0 else (lower, upper)
            )
            return WindowCoverage(
                ligatable=True,
                lower_nt=lower,
                upper_nt=upper,
                arm_3prime_nt=arm_3prime,
                arm_5prime_nt=arm_5prime,
            )
    return _UNCOVERED


class IsoformCoverageIndex:
    """Prepared exon tracks for one gene, queried once per candidate window."""

    def __init__(self, isoforms: Sequence[Dict[str, Any]]) -> None:
        self._tracks = [_ExonTrack(iso) for iso in isoforms]
        self._biotypes = {t.name: t.biotype for t in self._tracks}

    @property
    def biotypes(self) -> Dict[str, str]:
        return self._biotypes

    def biotype(self, isoform_name: str) -> str:
        return self._biotypes.get(isoform_name, "")

    def has_protein_coding(self) -> bool:
        return any(is_productive(t.biotype) for t in self._tracks)

    def is_productive_isoform(self, isoform_name: str) -> bool:
        return is_productive(self.biotype(isoform_name))

    def credit(
        self, segments: Sequence[Segment], strand: int, min_arm: int,
    ) -> Dict[str, WindowCoverage]:
        """Isoforms this window is ligatable on, with their arm geometry."""
        credited: Dict[str, WindowCoverage] = {}
        for track in self._tracks:
            coverage = _coverage(segments, track, strand)
            if coverage.credited(min_arm):
                credited[track.name] = coverage
        return credited


def credit_isoforms(
    segments: Sequence[Segment],
    isoforms: Sequence[Dict[str, Any]],
    *,
    strand: int,
    min_arm: int,
) -> Dict[str, WindowCoverage]:
    """One-shot :meth:`IsoformCoverageIndex.credit` for callers without a gene."""
    return IsoformCoverageIndex(isoforms).credit(segments, strand, min_arm)


def coding_first(isoforms: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Protein-coding transcripts first, canonical ahead of the rest.

    The search walks isoforms in order and the position-dedup keeps whichever
    one reached a window first, so this decides which transcript a shipped probe
    is recorded as designed against. Sorting is stable, so the provider's own
    order survives inside each group.
    """
    def priority(isoform: Dict[str, Any]) -> Tuple[int, int]:
        canonical = any(
            isoform.get(key)
            for key in ("is_canonical", "canonical", "mane_select")
        )
        coding = is_productive(isoform.get("biotype") or "")
        return (0 if coding else 1, 0 if canonical else 1)

    return sorted(isoforms, key=priority)


__all__ = [
    "PRODUCTIVE_BIOTYPES",
    "PROTEIN_CODING",
    "IsoformCoverageIndex",
    "WindowCoverage",
    "coding_first",
    "credit_isoforms",
    "is_productive",
    "window_coverage",
]
