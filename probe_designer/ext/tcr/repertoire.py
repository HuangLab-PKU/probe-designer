"""Clone-specificity screen against a patient's whole TCR repertoire.

A CDR3 padlock is only useful if it ligates on **one** clone. The nick-margin rule is a good
prior but not a guarantee: CDR3beta almost always opens with the germline-encoded CASS motif,
so a site a few nt inside the CDR3 can still lie in sequence every clone using that V gene
shares. When that happens both arms hybridise on the other clone's transcript, the nick is
fully paired, and the ligase closes it -- a false positive that no thermodynamic or BLAST
filter would flag, because the off-target is a real human transcript the panel is
deliberately looking at.

Two screens live here. ``build_repertoire_index`` / ``is_clone_specific`` ask whether the
EXACT footprint occurs in another clone. That is narrower than what the ligase needs: beyond
~3 nt either side of the nick a mismatch only lowers duplex stability, so a site one or two
bases off another clone's chain still anneals and is still sealed there. Found 2026-09-15
when the final check of the BZ09/BZ10/BZ12 mixing guide aligned every ordered TCR padlock
against its sample's clone list: 10 padlocks designed with the exact screen matched another
clone within 4 mismatches, all >= 10 nt from the nick (BZ12 clonotype58 vs clonotype55: one
base, 20 nt away). ``RepertoireScreen`` is the rule the pipeline uses now: the site counts on
another clone when the clamp around the nick pairs perfectly and the whole footprint carries
at most ``max_mismatches`` mismatches -- the criterion the BLAST path applies to off-target
hits. A near clone is therefore told apart only by a difference INSIDE the clamp, which is
where the screen pushes the chosen nick.

Both key on CDR3: two rows with the same CDR3 are one T-cell clone filed twice (the lab's
repertoire exports routinely carry a few), so they count once.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from probe_designer.chemistry import reverse_complement

# SplintR / PBCV-1 ligase clamps ~3 nt each side of the nick (Lohman 2014, Krzywkowski 2017);
# the same value the BLAST path uses for DEFAULT_CLAMP_NT.
LIGASE_CLAMP_NT = 3


def build_repertoire_index(
    clones: Iterable[Tuple[str, str]], kmer: int = 40,
) -> Dict[str, Set[str]]:
    """Map every ``kmer``-mer of every chain to the CDR3s that carry it.

    ``clones`` is ``(cdr3, chain)`` pairs. Both strands are indexed, so a
    footprint can be looked up as designed regardless of chemistry.
    """
    index: Dict[str, Set[str]] = defaultdict(set)
    for cdr3, chain in clones:
        chain = chain.upper()
        for seq in (chain, reverse_complement(chain)):
            for i in range(len(seq) - kmer + 1):
                index[seq[i:i + kmer]].add(cdr3.upper())
    return index


def host_cdr3s(footprint: str, index: Dict[str, Set[str]]) -> Set[str]:
    """Distinct CDR3s whose chain contains ``footprint`` exactly."""
    return index.get(footprint.upper(), set())


def is_clone_specific(
    footprint: str, index: Dict[str, Set[str]], own_cdr3: str = "",
) -> bool:
    """True when the exact footprint reaches no clone OTHER than ``own_cdr3``.

    Excluding the designed clone explicitly, rather than counting hosts and
    allowing one, keeps the answer correct whether or not the repertoire
    happens to contain the clone being designed. An empty index means "no
    repertoire supplied" -- the screen is then a no-op rather than a filter
    that rejects everything. This is the exact screen; see
    :class:`RepertoireScreen` for the one the pipeline applies.
    """
    if not index:
        return True
    return not (host_cdr3s(footprint, index) - {own_cdr3.upper()})


@dataclass
class RepertoireScreen:
    """Which OTHER clones a candidate site could still ligate on.

    A clone hosts the site when, on either strand of its chain, the ``2 * clamp`` nt around
    the nick match exactly and the whole footprint carries at most ``max_mismatches``
    mismatches. ``max_mismatches=0`` reduces to the exact screen. The seed index is built
    once per repertoire, so screening every candidate of every clone stays cheap.
    """

    clones: Sequence[Tuple[str, str]]
    kmer: int = 40
    clamp: int = LIGASE_CLAMP_NT
    max_mismatches: int = 4
    _strands: List[Tuple[str, str]] = field(default_factory=list, init=False, repr=False)
    _seeds: Dict[str, List[Tuple[int, int]]] = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self) -> None:
        if self.clamp < 1 or self.max_mismatches < 0:
            raise ValueError(f"clamp must be >= 1 and max_mismatches >= 0 (got {self.clamp}, {self.max_mismatches})")
        seeds: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        width = 2 * self.clamp
        for cdr3, chain in self.clones:
            chain = chain.upper()
            for seq in (chain, reverse_complement(chain)):
                k = len(self._strands)
                self._strands.append((cdr3.upper(), seq))
                for i in range(len(seq) - width + 1):
                    seeds[seq[i:i + width]].append((k, i))
        self._seeds = dict(seeds)

    def _alignments(self, footprint: str, nick: Optional[int]):
        fp = footprint.upper()
        nick = len(fp) // 2 if nick is None else nick
        if not self.clamp <= nick <= len(fp) - self.clamp:
            raise ValueError(f"nick {nick} leaves less than {self.clamp} nt of footprint on one side")
        lead = nick - self.clamp
        for k, i in self._seeds.get(fp[lead:nick + self.clamp], ()):
            cdr3, seq = self._strands[k]
            start = i - lead
            if start < 0 or start + len(fp) > len(seq):
                continue
            offsets = [j - nick for j, (a, b) in enumerate(zip(fp, seq[start:start + len(fp)])) if a != b]
            yield cdr3, seq, offsets

    def hosts(self, footprint: str, nick: Optional[int] = None) -> Dict[str, int]:
        """CDR3 -> fewest mismatches, for every clone the site could ligate on."""
        best: Dict[str, int] = {}
        for cdr3, _seq, offsets in self._alignments(footprint, nick):
            if len(offsets) <= self.max_mismatches and len(offsets) < best.get(cdr3, len(offsets) + 1):
                best[cdr3] = len(offsets)
        return best

    def mismatch_offsets(self, footprint: str, chain: str, nick: Optional[int] = None) -> List[int]:
        """Mismatch positions, relative to the nick, of the closest ligatable alignment on ``chain``."""
        wanted = {chain.upper(), reverse_complement(chain.upper())}
        found = [offsets for _cdr3, seq, offsets in self._alignments(footprint, nick) if seq in wanted]
        return min(found, key=len) if found else []

    def is_clone_specific(self, footprint: str, own_cdr3: str = "", nick: Optional[int] = None) -> bool:
        """True when the site could ligate on no clone other than ``own_cdr3``."""
        if not self.clones:
            return True
        return not (set(self.hosts(footprint, nick)) - {own_cdr3.upper()})


# Column names seen in the lab's repertoire exports. `sequence` is last on
# purpose: it is generic, so any specifically-named chain column wins over it.
_CHAIN_COLUMNS = ("TRB_nt", "TRA_nt", "chain_nt", "TR_nt", "sequence_nt",
                  "sequence")
_CDR3_COLUMNS = ("cdr3_nt", "cdr3_TRB_nt", "cdr3_TRA_nt", "cdr3")


def load_repertoire_file(path) -> list:
    """Read a repertoire table into ``[(cdr3, chain), ...]``.

    Accepts .csv or .xlsx and auto-detects the chain and CDR3 columns, so the
    per-patient exports (`ZCH_BZ07_spatial_TRB_sequence.csv` with
    `cdr3_TRB_nt`, `ZCH_BZ10_clone_cell_count_with_consensus_TRB.csv` with
    `cdr3_nt`) both work unchanged.
    """
    from pathlib import Path

    import pandas as pd

    path = Path(path)
    df = pd.read_csv(path) if path.suffix.lower() == ".csv" \
        else pd.read_excel(path)
    chain_col = next((c for c in _CHAIN_COLUMNS if c in df.columns), None)
    cdr3_col = next((c for c in _CDR3_COLUMNS if c in df.columns), None)
    if chain_col is None or cdr3_col is None:
        raise ValueError(
            f"{path.name}: need a chain column {_CHAIN_COLUMNS} and a CDR3 "
            f"column {_CDR3_COLUMNS}; found {list(df.columns)}"
        )
    out = []
    for _, r in df.iterrows():
        chain, cdr3 = r[chain_col], r[cdr3_col]
        if isinstance(chain, str) and isinstance(cdr3, str):
            out.append((cdr3.upper(), chain.upper()))
    if not out:
        raise ValueError(f"{path.name}: no usable (cdr3, chain) rows")
    return out
