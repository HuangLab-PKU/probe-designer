"""Build reusable thermodynamic annotations for a reference sequence.

First batch of per-position tracks, all computed at the CURRENT hybridization
conditions (a ``ReactionConditions``):

  * **arm Tm** — sliding-window arm melting temperature at the buffer
    (``filtering.thermo_profile``).
  * **accessibility** — per-base unpaired probability from folding the FULL
    mRNA (RNAplfold; ``filtering.accessibility``), i.e. the ΔG-derived
    accessibility, evaluated at the formamide-effective temperature.

Each track is written as bedGraph (IGV/UCSC-compatible) so it can live on NAS
beside the reference and be read during design / filtering / visualization
instead of being recomputed. The set is deliberately extensible — add more
track builders here as the annotation library grows.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

from probe_designer.chemistry import ReactionConditions
from probe_designer.filtering.thermo_profile import compute_tm_profile, write_bedgraph

try:
    from probe_designer.filtering.accessibility import (
        DEFAULT_SPAN,
        DEFAULT_WINDOW,
        compute_plfold_profile,
    )
    _HAS_ACCESSIBILITY = True
except ImportError:  # ViennaRNA missing
    _HAS_ACCESSIBILITY = False
    DEFAULT_WINDOW, DEFAULT_SPAN = 70, 40

logger = logging.getLogger(__name__)


def build_reference_annotations(
    seq: str,
    reference_id: str,
    reaction: ReactionConditions,
    *,
    out_dir: Path,
    arm_length: int = 20,
    chemistry: str = "dRNA",
    chrom: Optional[str] = None,
    start: int = 0,
    accessibility: bool = True,
    plfold_window: int = DEFAULT_WINDOW,
    plfold_span: int = DEFAULT_SPAN,
) -> List[Path]:
    """Write the first-batch thermodynamic annotation tracks for one reference.

    Args:
        seq: Reference transcript / mRNA sequence, 5'->3' (RNA-sense).
        reference_id: Identifier used in filenames and as the bedGraph chrom.
        reaction: The hybridization conditions to compute at (buffer + formamide).
        out_dir: Directory to write the ``.bedgraph`` tracks into (e.g. a NAS
            path beside the reference).
        arm_length: Arm/window length for the Tm track.
        chemistry: ``"dRNA"``/``"iLock"`` (DNA:RNA) or ``"cDNA"`` (DNA:DNA).
        chrom: bedGraph chrom column (defaults to ``reference_id``).
        start: 0-based reference offset of ``seq[0]`` (for genomic placement).
        accessibility: also write the RNAplfold accessibility track.
        plfold_window, plfold_span: RNAplfold W / L.

    Returns:
        Paths of the tracks written.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chrom = chrom or reference_id
    sig = reaction.signature()
    written: List[Path] = []

    # 1. Arm Tm at the buffer.
    tm = compute_tm_profile(seq, reaction, arm_length=arm_length, chemistry=chemistry)
    tm_path = out_dir / f"{reference_id}_armTm_L{arm_length}_{chemistry}_{sig}.bedgraph"
    write_bedgraph(tm, chrom, start, tm_path, track_name=f"{reference_id}_armTm")
    written.append(tm_path)
    logger.info("wrote arm-Tm annotation %s", tm_path.name)

    # 2. Accessibility from full-mRNA folding at the formamide-effective temp.
    if accessibility:
        if not _HAS_ACCESSIBILITY:
            logger.warning("ViennaRNA not available; skipping accessibility track")
        else:
            eff_t = reaction.effective_celsius
            acc = compute_plfold_profile(
                seq, window=plfold_window, span=plfold_span, temperature=eff_t,
            )
            acc_path = out_dir / f"{reference_id}_accessibility_T{eff_t:g}.bedgraph"
            write_bedgraph(
                acc, chrom, start, acc_path,
                track_name=f"{reference_id}_accessibility",
            )
            written.append(acc_path)
            logger.info("wrote accessibility annotation %s", acc_path.name)

    return written


def emit_annotations_for_sequences(
    seqs: dict,
    reaction: ReactionConditions,
    out_dir: Path,
    *,
    arm_length: int = 20,
    chemistry: str = "dRNA",
    accessibility: bool = True,
) -> List[Path]:
    """Write annotation tracks for many references (``{reference_id: seq}``).

    A thin batch wrapper over :func:`build_reference_annotations`; used to emit a
    design run's transcripts in one call. Non-string / empty sequences are
    skipped. Returns all track paths written.
    """
    written: List[Path] = []
    for reference_id, seq in seqs.items():
        if not isinstance(seq, str) or not seq:
            continue
        written.extend(build_reference_annotations(
            seq, str(reference_id), reaction,
            out_dir=out_dir, arm_length=arm_length,
            chemistry=chemistry, accessibility=accessibility,
        ))
    return written


__all__ = ["build_reference_annotations", "emit_annotations_for_sequences"]
