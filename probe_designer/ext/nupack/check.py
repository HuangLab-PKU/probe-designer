"""NUPACK cross-ligation Tier 3 — split arms, with the tether restored as
effective concentration.

Given ligator A and splint B, build a 3-strand tube
``{A.arm3_effective, A.arm5_effective, B.sequence}``, enumerate complexes up to
size 3, find the productive (arm3 . splint . arm5) ternary and report its free
energy (including coaxial stacking at the nick), its equilibrium concentration,
the fraction of splint tied up in it, and whether the MFE structure puts A's
3'-OH and 5'-P on adjacent bases of B.

**Why the arms are split — and why they must stay split.** The obvious model is
one tethered ligator strand plus the splint, so the backbone loop is handled by
the partition function. It does not work: a padlock wrapped on a splint is a
**pseudoknot**, and NUPACK's ensemble is nested (pseudoknot-free), so that
complex cannot be represented at all. Verified 2026-09-02 on a minimal case — a
32 nt ligator whose two 10 nt arms perfectly complement adjacent blocks of a
20 nt splint: the MFE pairs one arm (10 bp, -15.1 kcal/mol) and leaves the other
entirely unpaired, because pairing both would cross. Concretely, with the strands
concatenated as (ligator | splint), the pair from A's 5'-P and the pair from its
3'-OH always cross, at any register. Splitting the arms into separate strands is
what makes the geometry expressible; it is a workaround for a model limitation,
not a claim that the arms are independent. **Do not "fix" this back to two
strands.**

**What the split does cost, and how it is paid back.** Free arms lose the
tether: in reality, once one arm binds, the other is held next to the site at an
effective local concentration set by the backbone loop, not at the bulk probe
concentration. That is a large factor — roughly 500 uM for a 55 nt loop against
0.1 uM in the tube, i.e. ~5000x — and v2.3 dropped it entirely by giving the
arms the bulk concentration. Here the arm strands enter the tube at
:func:`tether_effective_concentration` instead, which restores the dominant part
of the cooperativity inside a model that can actually express the complex.

The result is an **upper bound** on the productive complex: both arms get the
tethered concentration, though physically only the second one to bind does.
Upper-bounding is the right direction for a safety screen, and this is Tier 3 —
the screening decision belongs to ``qc.cross_ligation``'s register scan, which is
exhaustive and costs microseconds. This module answers the more expensive
question, *how much* splint is tied up, for the handful of pairs that flags.
Feed it ``prefilter_pairs``; a full N x (N-1) sweep is neither affordable nor
useful.

``import nupack`` is LAZY — importing this module succeeds without NUPACK
installed; only :func:`screen_ternary` raises a clear ImportError.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from probe_designer.qc.cross_ligation import (
    ProbeForScreen,
    backbone_loop_nt,
    ligation_geometry,
)

from probe_designer.ext.nupack.config import (
    DEFAULT_CELSIUS,
    DEFAULT_ENSEMBLE,
    DEFAULT_FRACTION_THRESHOLD,
    DEFAULT_MAGNESIUM_M,
    DEFAULT_SODIUM_M,
    DEFAULT_STRAND_CONC_M,
    DEFAULT_VICINITY_N,
)


# ----------------------------------------------------------------------
# The tether
# ----------------------------------------------------------------------

#: ssDNA persistence length (nm) and rise per nucleotide (nm), for the
#: worm-like-chain mean-square end-to-end distance used below. Values are the
#: standard ones for ssDNA at ~100 mM monovalent salt.
SSDNA_PERSISTENCE_NM: float = 1.5
SSDNA_RISE_PER_NT_NM: float = 0.63

_AVOGADRO_PER_NM3_TO_M: float = 1e24 / 6.02214076e23   # molecules/nm^3 -> mol/L


def tether_effective_concentration(loop_nt: int) -> float:
    """Effective local concentration (M) of one chain end at the other.

    Jacobson-Stockmayer: the probability density of finding the far end of a
    Gaussian chain at the origin, ``(3 / (2 pi <r^2>))^(3/2)``, converted to
    molar. ``<r^2> = 2 * l_p * L_c`` for a worm-like chain, with ``L_c`` the
    contour length of the loop.

    This is what the backbone buys a padlock: for a 55 nt loop it comes out
    around 5e-4 M, five thousand times the ~1e-7 M bulk probe concentration.
    Modelling the arms as free oligos at bulk concentration throws that away.

    Raises ``ValueError`` for a non-positive loop, which would be a probe whose
    arms are contiguous — not a padlock.
    """
    if loop_nt <= 0:
        raise ValueError(
            f"tether loop must be > 0 nt, got {loop_nt}; a padlock's arms are "
            f"joined by a backbone, and with no loop there is nothing to tether"
        )
    contour_nm = SSDNA_RISE_PER_NT_NM * loop_nt
    mean_square_nm2 = 2.0 * SSDNA_PERSISTENCE_NM * contour_nm
    density_per_nm3 = (3.0 / (2.0 * math.pi * mean_square_nm2)) ** 1.5
    return density_per_nm3 * _AVOGADRO_PER_NM3_TO_M


#: The tether's contour, from ``qc.cross_ligation``'s single derivation of where
#: the arms sit in the assembled sequence.
_loop_length_nt = backbone_loop_nt


# ----------------------------------------------------------------------
# Public dataclass
# ----------------------------------------------------------------------


@dataclass
class TernaryHit:
    """One NUPACK ternary evaluation for an (A-ligator, B-splint) pair."""
    seq_a_id: str
    seq_b_id: str
    a_target: str = ""
    b_target: str = ""
    # Concatenation order in the canonical dot-bracket: arm3 + splint + arm5
    arm3_len: int = 0
    splint_len: int = 0
    arm5_len: int = 0
    # NUPACK results
    ternary_dg_kcal: float = 0.0          # full ternary dG at celsius (incl. coaxial)
    ternary_concentration_m: float = 0.0  # equilibrium [arm3.splint.arm5]
    ternary_fraction_of_b: float = 0.0    # ternary / total splint
    mfe_dotbracket: str = ""              # dot-bracket with '+' between strands
    mfe_nick_adjacent: bool = False       # 3'-OH and 5'-P on adjacent splint bases
    mfe_vicinity_contiguous: bool = False  # per-arm +/-n contiguous from terminus
    b_3oh_pos: Optional[int] = None       # splint 5'->3' partner of arm3's last base
    b_5p_pos: Optional[int] = None        # splint 5'->3' partner of arm5's first base
    # Audit / reproducibility
    tether_loop_nt: int = 0
    arm_conc_m: float = 0.0               # effective concentration used for the arms
    ensemble_used: str = "stacking"
    celsius_used: float = 55.0
    sodium_m: float = 0.075
    magnesium_m: float = 0.010
    strand_conc_m: float = 1.0e-7         # bulk concentration of the splint
    vicinity_n_each_side: int = 3

    @property
    def flagged(self) -> bool:
        """Enough splint tied up in a nick-adjacent complex to matter."""
        return (
            self.mfe_nick_adjacent
            and self.mfe_vicinity_contiguous
            and self.ternary_fraction_of_b > DEFAULT_FRACTION_THRESHOLD
        )


# ----------------------------------------------------------------------
# Dot-bracket parsing
# ----------------------------------------------------------------------


def _parse_dotbracket_pairs(dotbracket: str) -> Dict[int, int]:
    """Parse a dot-bracket string (with '+' strand separators) into
    ``{pos_i: pos_j}`` in CONCATENATED coordinates (ignoring '+' chars).
    """
    pairs: Dict[int, int] = {}
    stack: list[int] = []
    pos = 0
    for ch in dotbracket:
        if ch == "+":
            continue
        if ch == "(":
            stack.append(pos)
        elif ch == ")":
            if not stack:
                pos += 1
                continue
            opener = stack.pop()
            pairs[opener] = pos
            pairs[pos] = opener
        pos += 1
    return pairs


def _check_mfe_nick_geometry(
    pairs: Dict[int, int], arm3_len: int, splint_len: int,
) -> Tuple[bool, Optional[int], Optional[int]]:
    """In concatenated coords (arm3 | splint | arm5):

        arm3's 3'-OH = concat pos ``arm3_len - 1``
        arm5's 5'-P  = concat pos ``arm3_len + splint_len``

    Both must pair into the splint, and ``b_3oh - b_5p == 1`` — the same
    antiparallel arithmetic the register scan uses. Splint-local positions are
    returned.
    """
    p3 = pairs.get(arm3_len - 1)
    p5 = pairs.get(arm3_len + splint_len)
    if p3 is None or p5 is None:
        return False, None, None
    lo, hi = arm3_len, arm3_len + splint_len
    if not (lo <= p3 < hi) or not (lo <= p5 < hi):
        return False, None, None
    b_3oh, b_5p = p3 - arm3_len, p5 - arm3_len
    return (b_3oh - b_5p == 1), b_3oh, b_5p


def _check_mfe_vicinity_contiguous(
    pairs: Dict[int, int], arm3_len: int, splint_len: int, arm5_len: int,
    n_each_side: int,
) -> bool:
    """Ligase-clamp contiguity on both sides of the nick in the MFE structure.

    arm3's last ``n+1`` bases and arm5's first ``n+1`` must each pair into the
    splint along a strictly descending run — antiparallel, one base at a time,
    no bulges.
    """
    if n_each_side <= 0:
        return False
    if arm3_len <= n_each_side or arm5_len <= n_each_side:
        return False
    lo, hi = arm3_len, arm3_len + splint_len
    arm5_base = arm3_len + splint_len

    for positions in (
        range(arm3_len - 1 - n_each_side, arm3_len),          # arm3 side
        range(arm5_base, arm5_base + n_each_side + 1),        # arm5 side
    ):
        partners = []
        for a in positions:
            p = pairs.get(a)
            if p is None or not (lo <= p < hi):
                return False
            partners.append(p - arm3_len)
        if any(partners[i + 1] - partners[i] != -1 for i in range(len(partners) - 1)):
            return False
    return True


# ----------------------------------------------------------------------
# Core: single-direction NUPACK call
# ----------------------------------------------------------------------


def _evaluate_ternary_directional(
    *,
    ligator: ProbeForScreen, splint: ProbeForScreen,
    sodium_m: float, magnesium_m: float, celsius: float,
    ensemble: str, strand_conc_m: float, vicinity_n_each_side: int,
    arm_conc_m: Optional[float] = None,
) -> Optional[TernaryHit]:
    """One 3-strand NUPACK analysis. None if the ternary is not enumerated."""
    try:
        import nupack
    except ImportError as exc:
        raise ImportError(
            "NUPACK 4 is required for ext.nupack.check.screen_ternary. "
            "Install via academic registration at https://nupack.org/ "
            "and `pip install nupack-x.x.tar.gz` into the probe-design env."
        ) from exc

    arm3_eff, arm5_eff, _ = ligation_geometry(ligator)
    _, _, splint_seq = ligation_geometry(splint)

    loop_nt = _loop_length_nt(ligator)
    if arm_conc_m is None:
        arm_conc_m = tether_effective_concentration(loop_nt)

    s_arm3 = nupack.Strand(arm3_eff, name="arm3")
    s_splint = nupack.Strand(splint_seq, name="splint")
    s_arm5 = nupack.Strand(arm5_eff, name="arm5")

    tube = nupack.Tube(
        strands={
            s_arm3: arm_conc_m,
            s_splint: strand_conc_m,
            s_arm5: arm_conc_m,
        },
        complexes=nupack.SetSpec(max_size=3),
        name="cross_lig_tube",
    )
    model = nupack.Model(
        material="dna", celsius=celsius,
        sodium=sodium_m, magnesium=magnesium_m, ensemble=ensemble,
    )
    result = nupack.tube_analysis(tubes=[tube], model=model, compute=["mfe"])

    target_strands = {s_arm3, s_splint, s_arm5}
    productive = None
    concentration = 0.0
    free_energy = 0.0
    splint_total = 0.0
    for cplx, conc in result.tubes[tube].complex_concentrations.items():
        strands = list(cplx.strands)
        if s_splint in strands:
            splint_total += float(conc) * strands.count(s_splint)
        if len(strands) == 3 and set(strands) == target_strands:
            productive = cplx
            concentration = float(conc)
            free_energy = float(result.complexes[cplx].free_energy)

    if productive is None:
        return None

    mfe_structures = result.complexes[productive].mfe
    if not mfe_structures:
        return None
    dotbracket = str(mfe_structures[0].structure)

    # Remap NUPACK's own strand ordering onto the canonical (arm3 | splint | arm5).
    name_to_len = {
        "arm3": len(arm3_eff),
        "splint": len(splint_seq),
        "arm5": len(arm5_eff),
    }
    offsets: Dict[str, int] = {}
    cursor = 0
    for strand in productive.strands:
        offsets[strand.name] = cursor
        cursor += name_to_len[strand.name]

    remap: Dict[int, int] = {}
    for i in range(len(arm3_eff)):
        remap[offsets["arm3"] + i] = i
    for i in range(len(splint_seq)):
        remap[offsets["splint"] + i] = len(arm3_eff) + i
    for i in range(len(arm5_eff)):
        remap[offsets["arm5"] + i] = len(arm3_eff) + len(splint_seq) + i

    raw_pairs = _parse_dotbracket_pairs(dotbracket)
    pairs = {
        remap[k]: remap[v]
        for k, v in raw_pairs.items()
        if k in remap and v in remap
    }

    can_ligate, b_3oh, b_5p = _check_mfe_nick_geometry(
        pairs, len(arm3_eff), len(splint_seq),
    )
    vicinity_ok = (
        _check_mfe_vicinity_contiguous(
            pairs, len(arm3_eff), len(splint_seq), len(arm5_eff),
            vicinity_n_each_side,
        )
        if can_ligate else False
    )

    return TernaryHit(
        seq_a_id=ligator.probe_id, seq_b_id=splint.probe_id,
        a_target=ligator.target, b_target=splint.target,
        arm3_len=len(arm3_eff), splint_len=len(splint_seq),
        arm5_len=len(arm5_eff),
        ternary_dg_kcal=free_energy,
        ternary_concentration_m=concentration,
        ternary_fraction_of_b=(
            concentration / splint_total if splint_total > 0 else 0.0
        ),
        mfe_dotbracket=dotbracket,
        mfe_nick_adjacent=can_ligate,
        mfe_vicinity_contiguous=vicinity_ok,
        b_3oh_pos=b_3oh, b_5p_pos=b_5p,
        tether_loop_nt=loop_nt, arm_conc_m=arm_conc_m,
        ensemble_used=ensemble, celsius_used=celsius,
        sodium_m=sodium_m, magnesium_m=magnesium_m,
        strand_conc_m=strand_conc_m,
        vicinity_n_each_side=vicinity_n_each_side,
    )


# ----------------------------------------------------------------------
# Public entry point
# ----------------------------------------------------------------------


def screen_ternary(
    probes: List[ProbeForScreen],
    *,
    prefilter_pairs: Optional[Iterable[Tuple[str, str]]] = None,
    sodium_m: float = DEFAULT_SODIUM_M,
    magnesium_m: float = DEFAULT_MAGNESIUM_M,
    celsius: float = DEFAULT_CELSIUS,
    ensemble: str = DEFAULT_ENSEMBLE,
    strand_conc_m: float = DEFAULT_STRAND_CONC_M,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
    arm_conc_m: Optional[float] = None,
    on_progress: Optional[Any] = None,
) -> List[TernaryHit]:
    """NUPACK the productive ternary for directional (ligator, splint) pairs.

    Args:
        probes: full probe list.
        prefilter_pairs: ``(ligator_id, splint_id)`` DIRECTIONAL pairs to
            evaluate — normally the pairs ``qc.cross_ligation`` flagged. If
            None, every directional pair runs, which is O(N^2) NUPACK calls and
            only sane for a handful of probes.
        sodium_m, magnesium_m, celsius, ensemble: NUPACK Model parameters.
        strand_conc_m: bulk concentration of the splint strand.
        arm_conc_m: concentration for the two arm strands. Defaults to the
            ligator's own :func:`tether_effective_concentration`, restoring the
            backbone's cooperativity; pass ``strand_conc_m`` explicitly to get
            the v2.3 behaviour of untethered free arms.
        vicinity_n_each_side: ligase-clamp contiguity on the MFE structure.
        on_progress: optional ``(idx, total)`` callable for progress reporting.

    Returns one :class:`TernaryHit` per pair whose ternary NUPACK enumerated. A
    pair that raises is **not** skipped: the exception propagates with both
    probe ids attached, because a silently dropped pair is indistinguishable
    from a clean one in the report.

    Raises ``ImportError`` if NUPACK 4 is not installed.
    """
    probe_by_id = {p.probe_id: p for p in probes}

    if prefilter_pairs is None:
        pair_list: list[Tuple[str, str]] = []
        for a, b in combinations(probes, 2):
            pair_list.append((a.probe_id, b.probe_id))
            pair_list.append((b.probe_id, a.probe_id))
    else:
        pair_list = list(prefilter_pairs)

    hits: List[TernaryHit] = []
    total = len(pair_list)
    for idx, (lig_id, splint_id) in enumerate(pair_list):
        ligator = probe_by_id.get(lig_id)
        splint = probe_by_id.get(splint_id)
        if ligator is None or splint is None:
            raise KeyError(
                f"pair ({lig_id}, {splint_id}) references a probe absent from "
                f"the probe list"
            )
        try:
            hit = _evaluate_ternary_directional(
                ligator=ligator, splint=splint,
                sodium_m=sodium_m, magnesium_m=magnesium_m, celsius=celsius,
                ensemble=ensemble, strand_conc_m=strand_conc_m,
                vicinity_n_each_side=vicinity_n_each_side,
                arm_conc_m=arm_conc_m,
            )
        except ImportError:
            raise
        except Exception as exc:
            raise RuntimeError(
                f"NUPACK evaluation failed for ligator={lig_id} splint={splint_id}"
            ) from exc
        if hit is not None:
            hits.append(hit)
        if on_progress is not None:
            on_progress(idx + 1, total)
    return hits


# ----------------------------------------------------------------------
# Report writer
# ----------------------------------------------------------------------


REPORT_COLUMNS: Tuple[str, ...] = (
    "seq_a_id", "seq_b_id", "a_target", "b_target",
    "arm3_len", "splint_len", "arm5_len",
    "ternary_dg_kcal", "ternary_concentration_m", "ternary_fraction_of_b",
    "mfe_nick_adjacent", "mfe_vicinity_contiguous", "flagged",
    "b_3oh_pos", "b_5p_pos",
    "tether_loop_nt", "arm_conc_m",
    "ensemble_used", "celsius_used",
    "sodium_m", "magnesium_m", "strand_conc_m",
    "vicinity_n_each_side",
    "mfe_dotbracket",
)


def write_ternary_report(hits: List[TernaryHit], tsv_path: Path | str) -> None:
    """Write hits to TSV, most splint tied up first."""
    tsv_path = Path(tsv_path)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(hits, key=lambda h: h.ternary_fraction_of_b, reverse=True)

    def _esc(text: str) -> str:
        return (text or "").replace("\n", "\\n").replace("\t", "\\t")

    with open(tsv_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(list(REPORT_COLUMNS))
        for h in ordered:
            writer.writerow([
                h.seq_a_id, h.seq_b_id, h.a_target, h.b_target,
                h.arm3_len, h.splint_len, h.arm5_len,
                f"{h.ternary_dg_kcal:.3f}",
                f"{h.ternary_concentration_m:.3e}",
                f"{h.ternary_fraction_of_b:.3e}",
                h.mfe_nick_adjacent, h.mfe_vicinity_contiguous, h.flagged,
                h.b_3oh_pos if h.b_3oh_pos is not None else "",
                h.b_5p_pos if h.b_5p_pos is not None else "",
                h.tether_loop_nt, f"{h.arm_conc_m:.3e}",
                h.ensemble_used, f"{h.celsius_used:.1f}",
                f"{h.sodium_m:.4f}", f"{h.magnesium_m:.4f}",
                f"{h.strand_conc_m:.3e}",
                h.vicinity_n_each_side,
                _esc(h.mfe_dotbracket),
            ])


__all__ = [
    "TernaryHit",
    "screen_ternary",
    "write_ternary_report",
    "tether_effective_concentration",
    "REPORT_COLUMNS",
]
