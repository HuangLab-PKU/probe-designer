"""Pairwise heteroduplex prediction via primer3 (DNA nearest-neighbor).

Used by the **orthogonality screen** — flags any pair of probes in the same
panel whose duplex ΔG sits below a permissive threshold. Action is *log only*:
the panel manifest grows a warnings TSV, but probes are not auto-dropped.

**2026-07-20 (audit R5/P3): switched from ViennaRNA ``duplexfold`` to primer3.**
Padlock–padlock dimers are DNA:DNA, but ViennaRNA has no first-class DNA
parameter set, so the previous implementation scored them with RNA (Turner)
energies — the wrong physics, which also made the ΔG cutoffs rest on RNA
energetics. primer3 uses the SantaLucia unified DNA nearest-neighbor model and
takes the real buffer (monovalent / Mg2+ / dNTP / strand conc / temperature)
from a :class:`~probe_designer.chemistry.ReactionConditions`, so this screen now
shares the same physics and buffer as the cross-ligation screen.

Note on thresholds: ΔG magnitudes are NOT comparable to the old RNA-parameter
values. ``DEFAULT_DG_THRESHOLD`` is kept at the historical −12 kcal/mol as a
permissive log-only default; re-tune against a real panel if it is ever used to
drop probes.

This module does NOT model multi-strand competition. For concentration-aware
ensemble analysis (one probe competing across many partners) use the NUPACK
wrapper in ``ext/nupack``.
"""
from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, List, Optional, Sequence, Tuple

from probe_designer.chemistry import ReactionConditions

logger = logging.getLogger(__name__)


DEFAULT_DG_THRESHOLD = -12.0  # kcal/mol — orthogonality default (log-only)


@dataclass(frozen=True)
class DuplexHit:
    """One heteroduplex prediction.

    Attributes:
        probe_a_id: identifier for strand A.
        probe_b_id: identifier for strand B.
        delta_g: ΔG of the optimal duplex (kcal/mol; negative = stable).
        structure: primer3 4-line ASCII dimer structure.
        span_a: 0-indexed half-open ``(start, end)`` range on strand A that
            participates in the duplex (from the parsed base pairing).
        span_b: 0-indexed half-open ``(start, end)`` range on strand B.
    """
    probe_a_id: str
    probe_b_id: str
    delta_g: float
    structure: str
    span_a: Tuple[int, int]
    span_b: Tuple[int, int]


def _spans_from_ascii(ascii_structure: str) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Half-open paired spans on A and B from a primer3 ASCII dimer structure.

    Returns ``((0, 0), (0, 0))`` when no base pairs are present (or the ASCII is
    unparseable) — the parser degrades gracefully rather than raising.
    """
    from probe_designer.qc.dimer_ascii import parse_primer3_dimer_pairing  # lazy

    parsed: Dict = parse_primer3_dimer_pairing(ascii_structure or "")
    pairs = parsed.get("pair_positions_a_to_b") or {}
    if not pairs:
        return (0, 0), (0, 0)
    a_idx = sorted(pairs.keys())
    b_idx = sorted(pairs.values())
    return (a_idx[0], a_idx[-1] + 1), (b_idx[0], b_idx[-1] + 1)


def predict_pairwise_duplex(
    probe_a: str, probe_b: str,
    *,
    probe_a_id: str = "A",
    probe_b_id: str = "B",
    reaction: Optional[ReactionConditions] = None,
) -> Optional[DuplexHit]:
    """Optimal DNA:DNA heteroduplex between two probe sequences (primer3).

    Args:
        probe_a, probe_b: DNA sequences (5'->3').
        probe_a_id, probe_b_id: identifiers carried into the hit.
        reaction: buffer conditions; defaults to the protocol
            :class:`ReactionConditions` (its ``effective_celsius`` is the
            simulation temperature, so formamide raises stringency).

    Returns:
        ``DuplexHit``, or ``None`` for empty input or when primer3 finds no
        structure.
    """
    if not probe_a or not probe_b:
        return None
    import primer3  # lazy — keeps import cost off the hot path

    reaction = reaction or ReactionConditions()
    thermo = primer3.calc_heterodimer(
        probe_a.upper(), probe_b.upper(),
        output_structure=True,
        **reaction.primer3_kwargs(),
    )
    if not getattr(thermo, "structure_found", False):
        return None
    ascii_structure = str(getattr(thermo, "ascii_structure", "") or "")
    span_a, span_b = _spans_from_ascii(ascii_structure)
    return DuplexHit(
        probe_a_id=probe_a_id,
        probe_b_id=probe_b_id,
        delta_g=float(thermo.dg) / 1000.0,  # primer3 reports cal/mol
        structure=ascii_structure,
        span_a=span_a,
        span_b=span_b,
    )


def _pair_task(
    args: Tuple[str, str, str, str, ReactionConditions]
) -> Optional[DuplexHit]:
    """ProcessPoolExecutor friendly wrapper (must be importable by name)."""
    a, b, a_id, b_id, reaction = args
    return predict_pairwise_duplex(
        a, b, probe_a_id=a_id, probe_b_id=b_id, reaction=reaction,
    )


def screen_all_pairs(
    probes: Sequence[Tuple[str, str]],
    *,
    dg_threshold: float = DEFAULT_DG_THRESHOLD,
    reaction: Optional[ReactionConditions] = None,
    n_workers: Optional[int] = None,
) -> List[DuplexHit]:
    """All-pairs duplex screen across ``probes``.

    Args:
        probes: sequence of ``(probe_id, probe_sequence)`` tuples.
        dg_threshold: hits with ``delta_g <= dg_threshold`` are returned.
        reaction: buffer conditions (defaults to the protocol conditions).
        n_workers: parallelism. ``None``/0 runs serially (avoids ProcessPool
            startup cost on small panels); a positive int enables workers;
            ``-1`` uses ``os.cpu_count()``.

    Returns:
        List of ``DuplexHit`` for every pair (i, j), i<j, whose ΔG passes the
        threshold. The same pair is reported once.
    """
    items = list(probes)
    pairs = list(combinations(items, 2))
    if not pairs:
        return []

    reaction = reaction or ReactionConditions()
    tasks = [
        (a_seq, b_seq, a_id, b_id, reaction)
        for (a_id, a_seq), (b_id, b_seq) in pairs
    ]

    hits: List[DuplexHit] = []
    if n_workers is None or n_workers == 0:
        for t in tasks:
            hit = _pair_task(t)
            if hit is not None and hit.delta_g <= dg_threshold:
                hits.append(hit)
        return hits

    if n_workers < 0:
        import os
        n_workers = max(1, os.cpu_count() or 1)

    with ProcessPoolExecutor(max_workers=n_workers) as exe:
        futures = [exe.submit(_pair_task, t) for t in tasks]
        for fut in as_completed(futures):
            hit = fut.result()
            if hit is not None and hit.delta_g <= dg_threshold:
                hits.append(hit)
    return hits


__all__ = [
    "DEFAULT_DG_THRESHOLD",
    "DuplexHit",
    "predict_pairwise_duplex",
    "screen_all_pairs",
]
