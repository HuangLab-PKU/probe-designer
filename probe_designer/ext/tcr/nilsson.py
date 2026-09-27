"""Magoulopoulou/Nilsson padlock site-quality terms for TCR probe selection.

Source: Magoulopoulou et al. 2026, eBioMedicine 127:106264 (Mats Nilsson lab).
Their TCR-ISS pipeline splits germline V/J/C segments into 30-mers and keeps
only those with

* fewer than four repetitive bases,
* G or C at the 16th position (the ligation point of a 15+15 arm split),
* GC content between 40 and 60%.

We target patient-specific CDR3s rather than germline segments, so the CDR3
window admits only a handful of candidate positions and a hard filter would
leave some clones with no probe at all. The three rules therefore enter as a
*score contribution*: sites that satisfy them sort ahead of sites that do not,
without any site being dropped. Flags travel alongside the score so a run can
report which shipped sites break which rule.

Only the junction bonus is positive, so a compliant site can never outrank one
whose arm Tm or secondary structure is genuinely better by more than the bonus.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

# The paper's GC window, as fractions.
GC_WINDOW = (0.40, 0.60)

# "< 4 repetitive bases" — a run of 4 identical bases is already a violation.
MAX_HOMOPOLYMER_RUN = 3

# GC deviation is measured in units of 0.1 so one full 10-percentage-point miss
# costs the same as one excess homopolymer base.
_GC_DEVIATION_UNIT = 10.0


@dataclass(frozen=True)
class NilssonWeights:
    """Per-rule score weights. Set a weight to 0 to silence that rule."""

    homopolymer: float = 2.0
    junction_gc: float = 2.0
    gc: float = 2.0


def max_homopolymer_run(seq: str) -> int:
    """Length of the longest single-base run in ``seq`` (0 for empty)."""
    seq = seq.upper()
    best = run = 0
    prev = ""
    for base in seq:
        run = run + 1 if base == prev else 1
        prev = base
        if run > best:
            best = run
    return best


def gc_fraction(seq: str) -> float:
    """G+C as a fraction of length (0.0 for empty), not a percent."""
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def junction_base(target: str, arm_len: int, chemistry: str) -> str:
    """Target base at the ligation junction, on the strand the probe binds.

    Anchored on probe geometry: the base paired with the probe's 3'-OH
    terminus, which is what the paper's "16th position" of a 15+15 30-mer
    identifies. ``dRNA`` arms are the reverse complement of the target halves,
    so that base is ``target[arm_len]``; ``cDNA`` arms are the literal halves,
    putting the 3' terminus one base the other side of the nick at
    ``target[arm_len - 1]``. G/C-ness survives complementation, so reporting
    the mRNA-sense base answers the rule for either chemistry.
    """
    offset = {"dRNA": arm_len, "cDNA": arm_len - 1}[chemistry]
    return target.upper()[offset]


def nilsson_terms(
    target: str,
    arm_len: int,
    chemistry: str,
    weights: Optional[NilssonWeights] = None,
) -> Dict[str, Any]:
    """Score contribution + per-rule flags for one candidate binding site.

    ``target`` is the mRNA-sense BDS (``2 * arm_len`` nt). Returns the raw
    measurements, a boolean per rule, and ``nilsson_score`` — add it to the
    site's existing score.
    """
    weights = weights if weights is not None else NilssonWeights()
    target = target.upper()

    run = max_homopolymer_run(target)
    gc = gc_fraction(target)
    base = junction_base(target, arm_len, chemistry)

    gc_lo, gc_hi = GC_WINDOW
    gc_deviation = max(0.0, gc_lo - gc, gc - gc_hi)

    homopolymer_pass = run <= MAX_HOMOPOLYMER_RUN
    junction_gc_pass = base in ("G", "C")
    gc_pass = gc_deviation == 0.0

    score = (
        -weights.homopolymer * max(0, run - MAX_HOMOPOLYMER_RUN)
        + (weights.junction_gc if junction_gc_pass else 0.0)
        - weights.gc * _GC_DEVIATION_UNIT * gc_deviation
    )

    return {
        "homopolymer_run": run,
        "homopolymer_pass": homopolymer_pass,
        "junction_base": base,
        "junction_gc_pass": junction_gc_pass,
        "gc_content": round(gc, 4),
        "gc_pass": gc_pass,
        "nilsson_pass": homopolymer_pass and junction_gc_pass and gc_pass,
        "nilsson_score": round(score, 3),
    }
