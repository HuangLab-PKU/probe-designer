"""Pairwise heteroduplex prediction via ViennaRNA RNAduplex.

Phase 1B (2026-05-14) ships the shared backend used by:

* the **orthogonality consumer** — flags any pair of probes in the same
  panel whose duplex ΔG sits below a permissive threshold (~−12 kcal/mol).
  Action is *log only*: the panel manifest grows a warnings TSV, but
  probes are not auto-dropped.
* the Phase 2D **cross-ligation screen** — runs only on pairs already
  flagged by a fast k-mer junction-spanning seed; uses a tighter ΔG cut
  (~−10 kcal/mol) and a junction-span policy on top of the same backend.

The wrapper is intentionally thin around ``RNA.duplexfold`` because:

1. ``duplexfold`` returns the *optimal* duplex of two strands (single best
   alignment). For the panel-level screens this is sufficient — we are
   looking for the strongest unintended pairing, not enumerating all
   sub-optimal hits.
2. ViennaRNA's dot-bracket already encodes the alignment span on each
   strand; the parser below converts it to 0-indexed half-open Python
   ranges.

Temperature plumbing follows the Phase 1A convention: pass
``temperature`` (°C) and the global ``RNA.cvar.temperature`` is set
before calling duplexfold. Raise to 47–55 °C to mimic ~20% formamide
destabilization in our hybridization buffer.

This module does NOT model multi-strand competition. For ensemble
analysis (concentration-aware competition, useful when one probe could
pair with multiple others), use the NUPACK wrapper in Phase 2D's
``pool/ensemble.py``.
"""
from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, List, Optional, Sequence, Tuple

try:
    import RNA  # ViennaRNA bindings (>= 2.7.0)
    _HAS_VIENNARNA = True
except ImportError:
    _HAS_VIENNARNA = False

logger = logging.getLogger(__name__)


DEFAULT_DG_THRESHOLD = -12.0  # kcal/mol — orthogonality default
DEFAULT_TEMPERATURE = 37.0


@dataclass(frozen=True)
class DuplexHit:
    """One heteroduplex prediction.

    Attributes:
        probe_a_id: identifier for strand A.
        probe_b_id: identifier for strand B.
        delta_g: ΔG of the optimal duplex (kcal/mol; negative = stable).
        structure: ViennaRNA dot-bracket; ``&`` separates the two strands.
        span_a: 0-indexed half-open ``(start, end)`` range on strand A that
            participates in the duplex.
        span_b: 0-indexed half-open ``(start, end)`` range on strand B.
    """
    probe_a_id: str
    probe_b_id: str
    delta_g: float
    structure: str
    span_a: Tuple[int, int]
    span_b: Tuple[int, int]


def _parse_duplex_spans(
    structure: str, i: int, j: int,
) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Convert ViennaRNA's ``structure`` + ``i`` + ``j`` to 0-indexed spans.

    The structure looks like ``"(((((((&)))))))"`` where ``&`` separates
    strand A (left) and strand B (right). ViennaRNA reports:

    * ``i`` — 1-indexed END position on strand A (last base of left side
      structure).
    * ``j`` — 1-indexed START position on strand B (first base of right
      side structure).

    Returns half-open Python ranges so they slice directly:
    ``probe_a[span_a[0]:span_a[1]]`` is the duplex region on A.
    """
    if "&" not in structure:
        # Defensive: should never happen for duplexfold output.
        return (0, 0), (0, 0)
    left, right = structure.split("&", 1)
    len_left = len(left)
    len_right = len(right)
    # i is 1-indexed END on A; START_A_1idx = i - len_left + 1; convert to 0-idx half-open
    start_a = i - len_left  # 1-idx start - 1 = 0-idx start
    end_a = i               # 1-idx end inclusive == 0-idx exclusive end
    # j is 1-indexed START on B; END_B_1idx = j + len_right - 1
    start_b = j - 1
    end_b = j - 1 + len_right
    return (start_a, end_a), (start_b, end_b)


def predict_pairwise_duplex(
    probe_a: str, probe_b: str,
    *,
    probe_a_id: str = "A",
    probe_b_id: str = "B",
    temperature: float = DEFAULT_TEMPERATURE,
) -> Optional[DuplexHit]:
    """Compute optimal RNAduplex between two probe sequences.

    Returns ``None`` when ViennaRNA reports no duplex (rare; usually only
    for empty/very short inputs).

    The duplex is computed in the **DNA-DNA hybridization regime** for
    our padlock pool dimerization use case. ViennaRNA's RNA energy
    parameters are used (no first-class DNA parameter set in
    ViennaRNA); this is a known approximation. NUPACK ``material='dna'``
    is the rigorous alternative when concentration-aware competition
    matters (Phase 2D).
    """
    if not _HAS_VIENNARNA:
        raise RuntimeError(
            "ViennaRNA Python bindings not available. Install via "
            "`mamba install -n probe-design viennarna`."
        )
    if not probe_a or not probe_b:
        return None
    # ViennaRNA uses U; uppercase + T→U for safety.
    a = probe_a.upper().replace("T", "U")
    b = probe_b.upper().replace("T", "U")
    RNA.cvar.temperature = float(temperature)
    result = RNA.duplexfold(a, b)
    if result is None:
        return None
    span_a, span_b = _parse_duplex_spans(result.structure, result.i, result.j)
    return DuplexHit(
        probe_a_id=probe_a_id,
        probe_b_id=probe_b_id,
        delta_g=float(result.energy),
        structure=str(result.structure),
        span_a=span_a,
        span_b=span_b,
    )


def _pair_task(args: Tuple[str, str, str, str, float]) -> Optional[DuplexHit]:
    """ProcessPoolExecutor friendly wrapper (must be importable by name)."""
    a, b, a_id, b_id, temp = args
    return predict_pairwise_duplex(
        a, b, probe_a_id=a_id, probe_b_id=b_id, temperature=temp,
    )


def screen_all_pairs(
    probes: Sequence[Tuple[str, str]],
    *,
    dg_threshold: float = DEFAULT_DG_THRESHOLD,
    temperature: float = DEFAULT_TEMPERATURE,
    n_workers: Optional[int] = None,
) -> List[DuplexHit]:
    """All-pairs duplex screen across ``probes``.

    Args:
        probes: sequence of ``(probe_id, probe_sequence)`` tuples.
        dg_threshold: hits with ``delta_g <= dg_threshold`` are returned.
            Default −12 kcal/mol is the orthogonality threshold; cross-
            ligation callers should pass ``-10`` (tighter) plus their own
            junction-span filter.
        temperature: folding temperature in °C.
        n_workers: parallelism. ``None`` uses serial execution (avoids
            ProcessPool startup cost on small panels); set to a positive
            int to enable workers. ``-1`` uses ``os.cpu_count()``.

    Returns:
        List of ``DuplexHit`` for every pair (i, j), i<j, whose ΔG passes
        the threshold. The same pair is reported once.
    """
    items = list(probes)
    pairs = list(combinations(items, 2))
    if not pairs:
        return []

    tasks = [
        (a_seq, b_seq, a_id, b_id, temperature)
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
    "DEFAULT_TEMPERATURE",
    "DuplexHit",
    "predict_pairwise_duplex",
    "screen_all_pairs",
]
