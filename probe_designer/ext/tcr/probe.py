"""Single-sequence padlock probe design with an allowed-region constraint.

Designed for the TCR use case (scan a TRB or TRA mRNA, constrain ligation to
fall inside CDR3) but the API is generic — any caller with a single target
sequence and a "must-fall-in-this-substring" requirement can use it.

Each candidate site carries arm sequences and Tm values for BOTH chemistries
(direct mRNA padlock and cDNA-targeting padlock); chemistry selection happens
downstream at assembly time.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from Bio.SeqUtils import MeltingTemp as mt

from probe_designer.chemistry import (
    ReactionConditions,
    dna_revcomp_to_rna,
    reverse_complement,
)

try:
    import RNA  # ViennaRNA bindings — provides MFE via RNA.fold
    _HAS_VIENNARNA = True
except ImportError:
    _HAS_VIENNARNA = False


# Gap-preserving: aligned clone sequences carry '-', which the canonical
# implementation passes through unchanged.
_reverse_complement = reverse_complement


def validate_subseq_in_target(target_seq: str, allowed_subseq: str) -> int:
    """Verify ``allowed_subseq`` is a substring of ``target_seq`` and return
    its 0-indexed start position. Raises ``ValueError`` with a helpful
    message if not found (case-insensitive, N-stripped).
    """
    target = target_seq.replace("N", "").upper()
    sub = allowed_subseq.replace("N", "").upper()
    pos = target.find(sub)
    if pos == -1:
        raise ValueError(
            f"Allowed sub-sequence ({sub!r}) not found in target "
            f"(target length={len(target)}). For TCR usage this typically "
            f"means the CDR3 was not found in the TRB/TRA — check input pairing."
        )
    return pos


# Backwards-compat alias — TCR callers used to call this name.
validate_cdr3_in_trb = validate_subseq_in_target


class TcrProbeDesigner:
    """Design padlock probes for a single sequence with an allowed-region
    constraint on the ligation point.

    Despite the name (kept for the TCR use case it serves), the scan logic is
    generic: any caller with a target sequence and a sub-region the ligation
    must fall in can use it. For TCR, target = full TRA or TRB mRNA,
    allowed sub-region = CDR3.

    Each scan position emits arm sequences + Tm for BOTH the mRNA chemistry
    (arms = RC of target slice; Tm via ``R_DNA_NN1``) and the cDNA chemistry
    (arms = literal target slice; Tm via ``DNA_NN4``). Chemistry selection
    happens downstream — the same scan supports either or both.

    Parameters
    ----------
    bds_len : int
        Total binding site length (default 40).
    arm_len : int
        Length of each arm — half of ``bds_len`` (default 20).
    min_arm_tm, max_arm_tm : float
        mRNA-arm Tm window for filtering candidates (default 50.0–60.0 °C).
        Filtering uses the mRNA Tm only; cDNA Tm is informational.
    """

    def __init__(
        self,
        bds_len: int = 40,
        arm_len: Optional[int] = None,
        min_arm_tm: float = 50.0,
        max_arm_tm: float = 60.0,
        reaction: Optional[ReactionConditions] = None,
    ):
        self.bds_len = bds_len
        self.arm_len = arm_len if arm_len is not None else bds_len // 2
        if 2 * self.arm_len != self.bds_len:
            raise ValueError(
                f"arm_len ({self.arm_len}) * 2 must equal bds_len ({self.bds_len})"
            )
        self.min_arm_tm = min_arm_tm
        self.max_arm_tm = max_arm_tm
        self.reaction = reaction if reaction is not None else ReactionConditions()

    def _tm_hybrid(self, dna_arm: str) -> float:
        """Tm of a DNA probe arm hybridizing to mRNA.

        ``R_DNA_NN1`` requires the RNA (target) strand, so the DNA arm is
        reverse-complemented to its RNA-sense before the lookup. Buffer +
        formamide come from ``self.reaction``.
        """
        rna = dna_revcomp_to_rna(dna_arm)
        tm = mt.Tm_NN(rna, nn_table=mt.R_DNA_NN1, **self.reaction.tm_nn_kwargs())
        return self.reaction.apply_formamide(tm)

    def _tm_dna(self, seq: str) -> float:
        """Tm of a DNA:DNA duplex (cDNA chemistry), ``DNA_NN4`` at the buffer."""
        tm = mt.Tm_NN(seq, nn_table=mt.DNA_NN4, **self.reaction.tm_nn_kwargs())
        return self.reaction.apply_formamide(tm)

    # ------------------------------------------------------------------
    # Core scan
    # ------------------------------------------------------------------

    def find_binding_sites(
        self,
        target_seq: str,
        allowed_subseq: str,
        name: str,
    ) -> List[Dict[str, Any]]:
        """Scan ``target_seq`` for binding sites whose ligation point is
        inside ``allowed_subseq``.

        Returns one dict per candidate site, carrying BOTH mRNA and cDNA
        chemistry variants' arms + Tms. The caller filters / scores / selects
        which sites to ship.

        Returns ``[]`` (no error) if ``allowed_subseq`` is not a substring of
        ``target_seq`` — callers wanting a hard error should use
        :func:`validate_subseq_in_target` first.

        For backwards compatibility the dict keeps the field names
        ``cdr3_start`` / ``cdr3_end`` (used by Tm landscape plotting and old
        downstream xlsx readers) — they hold the bounds of ``allowed_subseq``
        regardless of whether the caller is doing TCR-CDR3 or something else.
        """
        target = target_seq.replace("N", "").upper()
        sub = allowed_subseq.replace("N", "").upper()
        sub_start = target.find(sub)
        if sub_start == -1:
            return []
        sub_end = sub_start + len(sub)

        bds_len = self.bds_len
        arm_len = self.arm_len

        search_start = max(0, sub_start - arm_len)
        search_end = min(len(target) - bds_len, sub_end - arm_len)

        sites: List[Dict[str, Any]] = []
        for pos in range(search_start, search_end + 1):
            if pos + bds_len > len(target):
                continue

            target_seq = target[pos: pos + bds_len]
            target_rc = _reverse_complement(target_seq)

            # mRNA chemistry: probe binds mRNA directly; arms = RC of target
            arm_3prime_dRNA = target_rc[:arm_len]
            arm_5prime_dRNA = target_rc[arm_len:]
            # cDNA chemistry: probe binds RT-derived cDNA; arms = literal target
            # The two chemistries have arms that are reverse-complements of each other
            arm_3prime_cDNA = target_seq[:arm_len]
            arm_5prime_cDNA = target_seq[arm_len:]

            # dRNA Tm: the probe binds mRNA, so R_DNA_NN1 must be fed the RNA
            # (target) strand. _tm_hybrid reverse-complements each DNA arm back
            # to its RNA-sense — Biopython's documented convention ("seq must be
            # the RNA sequence"). cDNA Tm uses the target half directly (DNA_NN4,
            # symmetric). Both at the reaction buffer + formamide.
            tm_full_rna = self._tm_hybrid(target_rc)
            tm_5prime_rna = self._tm_hybrid(arm_5prime_dRNA)
            tm_3prime_rna = self._tm_hybrid(arm_3prime_dRNA)
            tm_full_dna = self._tm_dna(target_seq)
            tm_5prime_dna = self._tm_dna(arm_5prime_cDNA)
            tm_3prime_dna = self._tm_dna(arm_3prime_cDNA)

            g_5 = target_seq[:arm_len].count("G") / arm_len
            g_3 = target_seq[arm_len:].count("G") / arm_len

            # MFE via ViennaRNA; None when unavailable (audit R8 — 0.0 would read
            # as "unstructured" rather than "not computed"). No try/except: a
            # ViennaRNA failure is a real problem and must surface, matching the
            # mRNA path; the pipeline isolates it per gene.
            mfe = RNA.fold(target_seq)[1] if _HAS_VIENNARNA else None

            sites.append({
                "gene_name": name,
                "pos": pos,
                "st": pos,
                "en": pos + bds_len - 1,
                "ligation_point": pos + arm_len,
                "target_sequence": target_seq,
                # mRNA chemistry
                "arm_3prime_dRNA": arm_3prime_dRNA.lower(),
                "arm_5prime_dRNA": arm_5prime_dRNA.lower(),
                "tm_dRNA": round(tm_full_rna, 2),
                "tm_5prime_dRNA": round(tm_5prime_rna, 2),
                "tm_3prime_dRNA": round(tm_3prime_rna, 2),
                "tm_diff_dRNA": round(abs(tm_5prime_rna - tm_3prime_rna), 2),
                # cDNA chemistry
                "arm_3prime_cDNA": arm_3prime_cDNA.lower(),
                "arm_5prime_cDNA": arm_5prime_cDNA.lower(),
                "tm_cDNA": round(tm_full_dna, 2),
                "tm_5prime_cDNA": round(tm_5prime_dna, 2),
                "tm_3prime_cDNA": round(tm_3prime_dna, 2),
                "tm_diff_cDNA": round(abs(tm_5prime_dna - tm_3prime_dna), 2),
                # G content (computed on target → same for both chemistries)
                "g_content": round((g_5 + g_3) / 2, 3),
                "g_content_5prime": round(g_5, 3),
                "g_content_3prime": round(g_3, 3),
                "mfe": None if mfe is None else round(mfe, 2),
                # Allowed-subseq bounds (kept under cdr3_* keys for backwards
                # compatibility with downstream Tm-landscape + xlsx readers).
                "cdr3_start": sub_start,
                "cdr3_end": sub_end,
            })
        return sites

    # Backwards-compat alias — preserves the old method name.
    find_cdr3_binding_sites = find_binding_sites

    # ------------------------------------------------------------------
    # Filtering + selection
    # ------------------------------------------------------------------

    def filter_by_chemistry_tm(
        self,
        candidates: List[Dict[str, Any]],
        *,
        chemistry: str,
        tm_range: Optional[tuple] = None,
    ) -> List[Dict[str, Any]]:
        """Keep only sites whose chemistry-specific 5' AND 3' arm Tms fall
        in the gate. ``chemistry`` is ``"dRNA"`` or ``"cDNA"``. If
        ``tm_range`` is None, falls back to ``[min_arm_tm, max_arm_tm]``
        (originally set on the designer)."""
        lo, hi = tm_range if tm_range is not None else (self.min_arm_tm, self.max_arm_tm)
        k5 = f"tm_5prime_{chemistry}"
        k3 = f"tm_3prime_{chemistry}"
        return [s for s in candidates if lo <= s[k5] <= hi and lo <= s[k3] <= hi]

    def select_non_overlapping_sites(
        self,
        candidates: List[Dict[str, Any]],
        *,
        min_gap: Optional[int] = None,
        max_sites: int = 3,
    ) -> List[Dict[str, Any]]:
        """Greedily select non-overlapping sites by descending ``score``.

        Each input site must carry a ``score`` key. If ``min_gap`` is None,
        defaults to ``self.bds_len`` (no two selected sites overlap).
        """
        if min_gap is None:
            min_gap = self.bds_len
        ranked = sorted(candidates, key=lambda x: x["score"], reverse=True)
        selected: List[Dict[str, Any]] = []
        for site in ranked:
            if all(abs(site["st"] - s["st"]) >= min_gap for s in selected):
                selected.append(site)
                if len(selected) >= max_sites:
                    break
        selected.sort(key=lambda x: x["st"])
        return selected
