"""Clone-specificity screen against a patient's whole TCR repertoire.

A CDR3 padlock is only useful if its 40 nt footprint occurs in **one** clone.
The nick-margin rule is a good prior but not a guarantee: CDR3beta almost
always opens with the germline-encoded CASS motif, so a site a few nt inside
the CDR3 can still lie in sequence every clone using that V gene shares. When
that happens both arms hybridise on the other clone's transcript, the nick is
fully paired, and the ligase closes it — a false positive that no thermodynamic
or BLAST filter would flag, because the off-target is a real human transcript
the panel is deliberately looking at.

Measured directly here: index every k-mer of every chain in the repertoire,
then ask how many DISTINCT CDR3s a candidate footprint reaches. Two rows with
the same CDR3 are one T-cell clone filed twice (the lab's repertoire exports
routinely carry a few), so they count once.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, Set, Tuple

from probe_designer.chemistry import reverse_complement


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
    """Distinct CDR3s whose chain contains ``footprint``."""
    return index.get(footprint.upper(), set())


def is_clone_specific(
    footprint: str, index: Dict[str, Set[str]], own_cdr3: str = "",
) -> bool:
    """True when the footprint reaches no clone OTHER than ``own_cdr3``.

    Excluding the designed clone explicitly, rather than counting hosts and
    allowing one, keeps the answer correct whether or not the repertoire
    happens to contain the clone being designed — and a repertoire that lists
    the same clone twice under different ids still counts once, because the
    index is keyed on CDR3.

    An empty index means "no repertoire supplied" — the screen is then a no-op
    rather than a filter that rejects everything.
    """
    if not index:
        return True
    return not (host_cdr3s(footprint, index) - {own_cdr3.upper()})


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
