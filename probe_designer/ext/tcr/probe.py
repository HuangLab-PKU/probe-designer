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

try:
    import RNA  # ViennaRNA bindings — provides MFE via RNA.fold
    _HAS_VIENNARNA = True
except ImportError:
    _HAS_VIENNARNA = False


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-", "N": "N"}
    return "".join(comp[b] for b in reversed(seq.upper()))


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
    ):
        self.bds_len = bds_len
        self.arm_len = arm_len if arm_len is not None else bds_len // 2
        if 2 * self.arm_len != self.bds_len:
            raise ValueError(
                f"arm_len ({self.arm_len}) * 2 must equal bds_len ({self.bds_len})"
            )
        self.min_arm_tm = min_arm_tm
        self.max_arm_tm = max_arm_tm

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
            arm_3prime_mRNA = target_rc[:arm_len]
            arm_5prime_mRNA = target_rc[arm_len:]
            # cDNA chemistry: probe binds RT-derived cDNA; arms = literal target
            # The two chemistries have arms that are reverse-complements of each other
            arm_3prime_cDNA = target_seq[:arm_len]
            arm_5prime_cDNA = target_seq[arm_len:]

            try:
                tm_full_rna = mt.Tm_NN(target_seq, nn_table=mt.R_DNA_NN1)
                tm_5prime_rna = mt.Tm_NN(target_seq[arm_len:], nn_table=mt.R_DNA_NN1)
                tm_3prime_rna = mt.Tm_NN(target_seq[:arm_len], nn_table=mt.R_DNA_NN1)
                tm_full_dna = mt.Tm_NN(target_seq, nn_table=mt.DNA_NN4)
                tm_5prime_dna = mt.Tm_NN(target_seq[arm_len:], nn_table=mt.DNA_NN4)
                tm_3prime_dna = mt.Tm_NN(target_seq[:arm_len], nn_table=mt.DNA_NN4)
            except Exception:
                tm_full_rna = tm_5prime_rna = tm_3prime_rna = 0.0
                tm_full_dna = tm_5prime_dna = tm_3prime_dna = 0.0

            g_5 = target_seq[:arm_len].count("G") / arm_len
            g_3 = target_seq[arm_len:].count("G") / arm_len

            # MFE via ViennaRNA when available; 0.0 otherwise.
            if _HAS_VIENNARNA:
                try:
                    _, mfe = RNA.fold(target_seq)
                except Exception:
                    mfe = 0.0
            else:
                mfe = 0.0

            sites.append({
                "gene_name": name,
                "pos": pos,
                "st": pos,
                "en": pos + bds_len - 1,
                "ligation_point": pos + arm_len,
                "target_sequence": target_seq,
                # mRNA chemistry
                "arm_3prime_mRNA": arm_3prime_mRNA.lower(),
                "arm_5prime_mRNA": arm_5prime_mRNA.lower(),
                "tm_mRNA": round(tm_full_rna, 2),
                "tm_5prime_mRNA": round(tm_5prime_rna, 2),
                "tm_3prime_mRNA": round(tm_3prime_rna, 2),
                "tm_diff_mRNA": round(abs(tm_5prime_rna - tm_3prime_rna), 2),
                # cDNA chemistry
                "arm_3prime_cDNA": arm_3prime_cDNA.lower(),
                "arm_5prime_cDNA": arm_5prime_cDNA.lower(),
                "tm_cDNA": round(tm_full_dna, 2),
                "tm_5prime_cDNA": round(tm_5prime_dna, 2),
                "tm_3prime_cDNA": round(tm_3prime_dna, 2),
                # G content (computed on target → same for both chemistries)
                "g_content": round((g_5 + g_3) / 2, 3),
                "g_content_5prime": round(g_5, 3),
                "g_content_3prime": round(g_3, 3),
                "mfe": round(mfe, 2),
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

    def filter_by_mrna_tm(
        self, candidates: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Keep only sites whose 5' AND 3' mRNA arm Tms fall in the configured
        ``[min_arm_tm, max_arm_tm]`` window."""
        return [
            s for s in candidates
            if self.min_arm_tm <= s["tm_5prime_mRNA"] <= self.max_arm_tm
            and self.min_arm_tm <= s["tm_3prime_mRNA"] <= self.max_arm_tm
        ]

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
