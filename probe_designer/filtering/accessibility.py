"""Per-transcript binding-site accessibility via ViennaRNA RNAplfold.

Phase 1A (2026-05-14) replaces the per-window ``RNA.fold`` self-fold gate in
``thermal_filter`` with a per-base unpaired-probability profile computed
over the full transcript context. The biological rationale: a candidate
binding window may look unstructured when folded in isolation but be
locked behind a long-range stem in the real mRNA, or vice versa. Folding
the entire transcript and asking "is this window open in context?" is
the right physical question.

Implementation: ``ViennaRNA.pfl_fold_up(seq, ulength=1, W, L)`` returns a
length-``len(seq)+1`` tuple of ``(_, p_unpaired_ulen1)`` rows. Index 0 is
a placeholder; positions 1..N map to base 1..N. We expose a clean
``compute_plfold_profile`` returning a 0-indexed numpy array of length N.

Window geometry (2026-07-21, audit R9): W and L come from
``chemistry.FoldingConditions``, defaulting to Lange 2012's W = L + 50 with
L ≈ 100. The previous W = 70 / L = 40 could not represent a base pair spanning
more than 40 nt, so any site locked behind a long-range stem read as open —
on real transcripts the shorter geometry over-states mean p_unpaired by ≈ 0.10.

Experimental-conditions note: ViennaRNA defaults assume buffered conditions
with no denaturant. Our hybridization buffer (~20% formamide) destabilizes RNA
secondary structure, so folding happens at the *solvent-effective* temperature
rather than a nominal 37 °C. ``FoldingConditions.temperature_c`` is ``None`` by
default, meaning "track ``ReactionConditions.effective_celsius``" — this is what
keeps the accessibility used to filter identical to the accessibility written
into the reference annotation. Pin it to a number to decouple the two.

Caching: per-(transcript_id, len, W, L, temperature) profiles are stored
as ``.npy`` files under a caller-supplied directory. The cache key omits
the sequence itself but verifies length on load so we never deliver a
mismatched profile.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np

from probe_designer.chemistry import (
    DEFAULT_FOLD_TEMPERATURE_C,
    FoldingConditions,
)

try:
    import RNA  # ViennaRNA Python bindings (>= 2.7.0)
    _HAS_VIENNARNA = True
except ImportError:
    _HAS_VIENNARNA = False

logger = logging.getLogger(__name__)


# Single source of truth for the geometry: chemistry.FoldingConditions (audit
# R9 / Lange 2012). These aliases stay so callers can keep importing constants,
# but they are derived, not a second copy.
_DEFAULTS = FoldingConditions()
DEFAULT_WINDOW = _DEFAULTS.window
DEFAULT_SPAN = _DEFAULTS.span
DEFAULT_TEMPERATURE = DEFAULT_FOLD_TEMPERATURE_C


def compute_plfold_profile(
    seq: str,
    *,
    window: int = DEFAULT_WINDOW,
    span: int = DEFAULT_SPAN,
    temperature: float = DEFAULT_TEMPERATURE,
) -> np.ndarray:
    """Per-base unpaired probability for ``seq`` via ViennaRNA RNAplfold.

    Args:
        seq: transcript sequence (DNA or RNA letters; U/T are equivalent for
            ViennaRNA's purposes — internally normalised).
        window: sliding window size W (default 150 = L + 50, Lange 2012).
            Lower for short transcripts (e.g. 40 for CDR3 contexts); it is
            clipped to ``len(seq)`` automatically.
        span: maximum base-pair span L within the window (default 100).
        temperature: folding temperature in °C (default 37.0, the bare
            ViennaRNA convention). Callers with a buffer should pass
            ``**FoldingConditions().plfold_kwargs(reaction)`` instead of
            hand-picking this.

    Returns:
        ``np.ndarray`` of length ``len(seq)`` with ``p_unpaired`` ∈ [0, 1]
        for each base (0-indexed).

    Raises:
        RuntimeError: if ViennaRNA is not importable.
        ValueError: if ``seq`` is empty or ``window`` exceeds ``len(seq)``.
    """
    if not _HAS_VIENNARNA:
        raise RuntimeError(
            "ViennaRNA Python bindings not available. Install via "
            "`mamba install -n probe-design viennarna`."
        )
    if not seq:
        raise ValueError("seq must be non-empty")
    n = len(seq)
    # ViennaRNA accepts T/U interchangeably; uppercase for safety.
    rna_seq = seq.upper().replace("T", "U")
    # plfold needs window <= len(seq); clip to avoid silent garbage.
    effective_W = min(window, n)
    effective_L = min(span, effective_W)

    # Set fold parameters via model details.
    md = RNA.md()
    md.temperature = float(temperature)
    # Apply globally — RNA.pfl_fold_up reads global params on this version.
    RNA.cvar.temperature = float(temperature)

    raw = RNA.pfl_fold_up(rna_seq, 1, effective_W, effective_L)
    # raw is a tuple of length n+1: entry 0 is a placeholder, entries 1..n
    # are (p_unpaired_ulen0=0, p_unpaired_ulen1) tuples per position.
    profile = np.fromiter(
        (raw[i + 1][1] for i in range(n)), dtype=np.float64, count=n,
    )
    # Numerical hygiene: clip to [0, 1] in case ViennaRNA returns very
    # small negatives from cancellation.
    np.clip(profile, 0.0, 1.0, out=profile)
    return profile


def mean_open_probability(
    profile: np.ndarray, start: int, end: int,
) -> float:
    """Mean p_unpaired over the half-open window ``[start, end)``.

    Args:
        profile: per-base profile returned by :func:`compute_plfold_profile`.
        start: inclusive 0-based start position.
        end: exclusive 0-based end position. Must satisfy
            ``0 <= start < end <= len(profile)``.

    Returns:
        Mean of ``profile[start:end]``. Returns ``0.0`` when the window
        falls outside the profile (caller should treat this as "fail
        gate"; logged at DEBUG to aid diagnostics).
    """
    n = len(profile)
    if not (0 <= start < end <= n):
        logger.debug(
            "mean_open_probability: window [%d, %d) outside profile of length %d",
            start, end, n,
        )
        return 0.0
    return float(profile[start:end].mean())


def aggregate_open_probability(
    per_isoform: dict[str, float],
    *,
    weights: Optional[dict[str, float]] = None,
) -> float:
    """Combine per-isoform mean-open probabilities into one score.

    Args:
        per_isoform: mapping ``{isoform_id: mean_open_probability}``.
        weights: optional ``{isoform_id: weight}`` — typically per-cell-type
            expression weights derived from scRNA-seq. Missing isoforms
            default to weight 0 (treated as not expressed). When None,
            uniform weighting (simple mean across isoforms).

    Returns:
        Weighted mean open probability. Returns ``0.0`` for empty input.

    Notes:
        Why weighting matters: a binding site might be fully open in a
        low-abundance isoform but folded in the dominant one. We want
        probes that work on whatever isoform the cell actually expresses,
        so the dominant isoform should drive the score.
    """
    if not per_isoform:
        return 0.0
    if weights is None:
        return float(sum(per_isoform.values()) / len(per_isoform))
    total_w = 0.0
    weighted = 0.0
    for iso, p_open in per_isoform.items():
        w = float(weights.get(iso, 0.0))
        if w <= 0:
            continue
        weighted += p_open * w
        total_w += w
    if total_w <= 0:
        return 0.0
    return weighted / total_w


def cached_profile(
    seq: str,
    transcript_id: str,
    *,
    cache_dir: Path,
    window: int = DEFAULT_WINDOW,
    span: int = DEFAULT_SPAN,
    temperature: float = DEFAULT_TEMPERATURE,
) -> np.ndarray:
    """Disk-cached plfold profile.

    Cache key: ``<cache_dir>/<transcript_id>_W{W}_L{L}_T{T}_N{len}.npy``.
    The sequence is NOT in the key — instead the cached file's length is
    checked on load; on mismatch the file is silently regenerated. This
    guards against the case where a transcript ID is reused but the
    underlying sequence changed (e.g. genome assembly update).
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = f"{transcript_id}_W{window}_L{span}_T{temperature:g}_N{len(seq)}.npy"
    path = cache_dir / key
    if path.exists():
        try:
            arr = np.load(path)
            if arr.shape == (len(seq),):
                return arr
            logger.warning(
                "plfold cache %s shape mismatch (%s vs len %d); regenerating",
                path, arr.shape, len(seq),
            )
        except (OSError, ValueError) as exc:
            logger.warning("plfold cache %s unreadable (%s); regenerating",
                           path, exc)
    arr = compute_plfold_profile(
        seq, window=window, span=span, temperature=temperature,
    )
    try:
        np.save(path, arr)
    except OSError as exc:
        logger.warning("failed to write plfold cache %s: %s", path, exc)
    return arr


__all__ = [
    "DEFAULT_WINDOW", "DEFAULT_SPAN", "DEFAULT_TEMPERATURE",
    "compute_plfold_profile", "mean_open_probability",
    "aggregate_open_probability", "cached_profile",
]
