"""TCR clone-specific padlock probe design module.

Designs padlock probes targeting CDR3 regions of T-cell receptor sequences.
The ligation point is constrained to fall within CDR3, while arms may extend
outside it.
"""

from typing import List, Dict, Any, Optional
from Bio.SeqUtils import MeltingTemp as mt

from probe_designer.filtering import SequenceFilter


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}
    return "".join(comp[b] for b in reversed(seq.upper()))


class TcrProbeDesigner:
    """Design padlock probes for TCR clone CDR3 regions.

    Parameters
    ----------
    bds_len : int
        Total binding site length (default 40).
    arm_len : int
        Length of each arm (default bds_len // 2).
    min_arm_tm : float
        Minimum acceptable arm Tm (default 50.0).
    max_arm_tm : float
        Maximum acceptable arm Tm (default 60.0).
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
        self.min_arm_tm = min_arm_tm
        self.max_arm_tm = max_arm_tm

    # ------------------------------------------------------------------
    # Core helpers
    # ------------------------------------------------------------------

    def find_cdr3_binding_sites(
        self,
        trb_seq: str,
        cdr3_seq: str,
        gene_name: str,
    ) -> List[Dict[str, Any]]:
        """Find binding sites where the ligation point falls inside CDR3.

        The ligation point is defined as ``pos + arm_len`` (the junction
        between the two padlock arms on the target strand).

        Parameters
        ----------
        trb_seq : str
            Full TRB nucleotide sequence.
        cdr3_seq : str
            CDR3 nucleotide sub-sequence (must be found within *trb_seq*).
        gene_name : str
            Identifier for this clone / gene.

        Returns
        -------
        list[dict]
            Candidate binding sites with thermodynamic annotations.
        """
        trb_seq = trb_seq.replace("N", "").upper()
        cdr3_seq = cdr3_seq.replace("N", "").upper()

        cdr3_start = trb_seq.find(cdr3_seq)
        if cdr3_start == -1:
            return []
        cdr3_end = cdr3_start + len(cdr3_seq)

        bds_len = self.bds_len
        arm_len = self.arm_len

        search_start = max(0, cdr3_start - arm_len)
        search_end = min(len(trb_seq) - bds_len, cdr3_end - arm_len)

        sites: List[Dict[str, Any]] = []
        for pos in range(search_start, search_end + 1):
            if pos + bds_len > len(trb_seq):
                continue

            target_seq = trb_seq[pos : pos + bds_len]
            target_rc = _reverse_complement(target_seq)
            arm_3prime = target_rc[:arm_len]
            arm_5prime = target_rc[arm_len:]

            try:
                tm_full = mt.Tm_NN(target_seq, nn_table=mt.R_DNA_NN1)
                tm_5prime = mt.Tm_NN(target_seq[:arm_len], nn_table=mt.R_DNA_NN1)
                tm_3prime = mt.Tm_NN(target_seq[arm_len:], nn_table=mt.R_DNA_NN1)
            except Exception:
                tm_full = tm_5prime = tm_3prime = 0.0

            g_content_5 = target_seq[:arm_len].count("G") / arm_len
            g_content_3 = target_seq[arm_len:].count("G") / arm_len
            g_content = (g_content_5 + g_content_3) / 2.0

            sites.append(
                {
                    "gene_name": gene_name,
                    "pos": pos,
                    "st": pos,
                    "en": pos + bds_len - 1,
                    "ligation_point": pos + arm_len,
                    "sequence": arm_3prime + arm_5prime,
                    "arm_5prime": arm_5prime.lower(),
                    "arm_3prime": arm_3prime.lower(),
                    "target_sequence": target_seq,
                    "tm": round(tm_full, 2),
                    "tm_5prime": round(tm_5prime, 2),
                    "tm_3prime": round(tm_3prime, 2),
                    "tm_diff": round(abs(tm_5prime - tm_3prime), 2),
                    "g_content": round(g_content, 3),
                    "g_content_5prime": round(g_content_5, 3),
                    "g_content_3prime": round(g_content_3, 3),
                    "cdr3_start": cdr3_start,
                    "cdr3_end": cdr3_end,
                }
            )

        return sites

    def select_non_overlapping_sites(
        self,
        candidates: List[Dict[str, Any]],
        min_gap: Optional[int] = None,
        max_sites: int = 3,
    ) -> List[Dict[str, Any]]:
        """Greedily select non-overlapping sites by descending score.

        Parameters
        ----------
        candidates : list[dict]
            Candidate sites; each must contain ``"st"`` and ``"score"`` keys.
        min_gap : int, optional
            Minimum start-position gap between selected sites
            (default ``self.bds_len``).
        max_sites : int
            Maximum number of sites to select.

        Returns
        -------
        list[dict]
            Selected sites sorted by position.
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

    # ------------------------------------------------------------------
    # High-level entry point
    # ------------------------------------------------------------------

    def design_tcr_probes(
        self,
        clone_sequences: Dict[str, str],
        cdr3_regions: Dict[str, str],
        max_sites_per_clone: int = 3,
        min_gap: Optional[int] = None,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Design probes for multiple TCR clones.

        Parameters
        ----------
        clone_sequences : dict[str, str]
            Mapping of clone/gene name to full TRB nucleotide sequence.
        cdr3_regions : dict[str, str]
            Mapping of clone/gene name to CDR3 nucleotide sub-sequence.
        max_sites_per_clone : int
            Maximum probes to select per clone.
        min_gap : int, optional
            Minimum gap between selected binding sites.

        Returns
        -------
        dict[str, list[dict]]
            Mapping of clone name to selected probe sites.
        """
        results: Dict[str, List[Dict[str, Any]]] = {}

        for name in clone_sequences:
            trb = clone_sequences[name]
            cdr3 = cdr3_regions.get(name, "")
            if not cdr3:
                continue

            candidates = self.find_cdr3_binding_sites(trb, cdr3, name)
            if not candidates:
                continue

            # Tm filter
            passed = [
                s
                for s in candidates
                if self.min_arm_tm <= s["tm_5prime"] <= self.max_arm_tm
                and self.min_arm_tm <= s["tm_3prime"] <= self.max_arm_tm
            ]
            if not passed:
                continue

            # Score: higher is better (less negative MFE + small Tm diff)
            for s in passed:
                s["score"] = round(-s["tm_diff"], 3)

            selected = self.select_non_overlapping_sites(
                passed, min_gap=min_gap, max_sites=max_sites_per_clone
            )
            if selected:
                results[name] = selected

        return results
