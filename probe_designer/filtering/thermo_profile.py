"""Per-position melting-temperature profile — thermo as a reference annotation.

For a reference transcript/mRNA sequence, a chemistry, an arm length, and a
ReactionConditions buffer, compute the arm Tm at every position (sliding
window). This turns Tm into a reusable annotation on the genome/RNA reference:
compute once per (reference, condition, arm length), cache to disk, and read
during design / filtering / visualization instead of recomputing per candidate.

For DNA:RNA (direct-mRNA / iLock) the reference window IS the RNA-sense strand
that Biopython's ``R_DNA_NN1`` requires (no reverse-complement needed); for
DNA:DNA (cDNA) the window is used directly with ``DNA_NN4``. Caching mirrors
``filtering.accessibility.cached_profile``.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from Bio.SeqUtils import MeltingTemp as mt

from probe_designer.chemistry import ReactionConditions

logger = logging.getLogger(__name__)

DEFAULT_ARM_LENGTH = 20


def _nn_table(chemistry: str):
    """dRNA / iLock -> DNA:RNA hybrid table; cDNA -> DNA:DNA table."""
    return mt.DNA_NN4 if str(chemistry).lower() == "cdna" else mt.R_DNA_NN1


def compute_tm_profile(
    seq: str,
    reaction: ReactionConditions,
    *,
    arm_length: int = DEFAULT_ARM_LENGTH,
    chemistry: str = "dRNA",
) -> np.ndarray:
    """Arm Tm at every start position along ``seq`` (an RNA-sense reference).

    Args:
        seq: Reference transcript / mRNA sequence, 5'->3' (RNA-sense).
        reaction: Buffer conditions the Tm is computed at.
        arm_length: Arm length (window size) in nt.
        chemistry: ``"dRNA"``/``"iLock"`` (DNA:RNA) or ``"cDNA"`` (DNA:DNA).

    Returns:
        Array of length ``len(seq)``; entry ``i`` is the Tm (°C, buffer +
        formamide) of the ``arm_length``-mer arm binding at
        ``seq[i:i+arm_length]``. The trailing ``arm_length-1`` positions and any
        ambiguous window are ``NaN``.
    """
    seq = seq.upper()
    n = len(seq)
    table = _nn_table(chemistry)
    buffer = reaction.tm_nn_kwargs()
    out = np.full(n, np.nan, dtype=np.float64)
    for i in range(n - arm_length + 1):
        window = seq[i:i + arm_length]
        try:
            tm = mt.Tm_NN(window, nn_table=table, **buffer)
        except ValueError:
            continue
        out[i] = reaction.apply_formamide(tm)
    return out


def cached_tm_profile(
    seq: str,
    reference_id: str,
    reaction: ReactionConditions,
    *,
    cache_dir: Path,
    arm_length: int = DEFAULT_ARM_LENGTH,
    chemistry: str = "dRNA",
) -> np.ndarray:
    """Disk-cached :func:`compute_tm_profile`.

    Cache key: ``<reference_id>_<chemistry>_L{arm}_<buffer_sig>_N{len}.npy``, so
    a different buffer or arm length is a different annotation. The sequence
    length is verified on load; a mismatch (reused id, changed sequence)
    silently regenerates.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    chem = str(chemistry).lower()
    key = (f"{reference_id}_{chem}_L{arm_length}"
           f"_{reaction.signature()}_N{len(seq)}.npy")
    path = cache_dir / key
    if path.exists():
        try:
            arr = np.load(path)
            if arr.shape == (len(seq),):
                return arr
            logger.warning("tm profile cache %s shape mismatch; regenerating", path)
        except (OSError, ValueError) as exc:
            logger.warning("tm profile cache %s unreadable (%s); regenerating", path, exc)
    arr = compute_tm_profile(seq, reaction, arm_length=arm_length, chemistry=chem)
    try:
        np.save(path, arr)
    except OSError as exc:
        logger.warning("failed to write tm profile cache %s: %s", path, exc)
    return arr


def write_bedgraph(
    profile: np.ndarray,
    chrom: str,
    start: int,
    path: Path,
    *,
    track_name: str = "arm_tm",
) -> Path:
    """Write a per-position Tm profile as a bedGraph — an IGV/UCSC-compatible
    genome annotation to store on NAS alongside the reference.

    Each non-NaN ``profile[i]`` becomes a single-base interval at
    ``start + i`` carrying the arm Tm anchored there; equal adjacent values are
    coalesced. NaN positions (no full window) are skipped.

    Args:
        profile: 1-D array from :func:`compute_tm_profile`.
        chrom: chromosome / contig / transcript id the profile sits on.
        start: 0-based reference start of ``profile[0]``.
        path: output ``.bedgraph`` path.
        track_name: bedGraph track name shown in IGV.

    Returns:
        ``path``. Convert to indexed **BigWig** for genome-scale tracks with
        UCSC ``bedGraphToBigWig`` (needs a chrom.sizes) or ``pyBigWig``.
    """
    profile = np.asarray(profile, dtype=float)
    n = len(profile)
    lines = [f'track type=bedGraph name="{track_name}"']
    i = 0
    while i < n:
        v = profile[i]
        if np.isnan(v):
            i += 1
            continue
        j = i + 1
        while j < n and profile[j] == v:
            j += 1
        lines.append(f"{chrom}\t{start + i}\t{start + j}\t{v:.2f}")
        i = j
    path = Path(path)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


__all__ = [
    "DEFAULT_ARM_LENGTH",
    "compute_tm_profile",
    "cached_tm_profile",
    "write_bedgraph",
]
