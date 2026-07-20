"""NUPACK 4 ternary-complex cross-ligation Tier 3 (v2.3).

Given probes A (ligator) and B (splint), build a 3-strand NUPACK tube
``{A.arm3_effective, A.arm5, B.sequence}``, enumerate complexes up to
size 3 with ``SetSpec(max_size=3)``, find the productive
(arm3·splint·arm5) ternary complex, and report:

* its full ternary ΔG (incl. coaxial stacking at the nick, via
  ``Model(ensemble='stacking')`` — the DEFAULT)
* its equilibrium concentration at lab strand concentrations
* fraction of splint participating in productive ternary
* MFE dot-bracket structure of the ternary
* MFE-geometry flags: nick-adjacent + per-arm vicinity ±n contiguous

This module is the OPTIONAL Tier 3 above qc/cross_ligation.py's v2.2
split-fragment screen. ``import nupack`` is LAZY — top-level
``from probe_designer.ext.nupack.check import screen_ternary`` succeeds
without nupack installed; only :func:`screen_ternary` itself raises a
clear ImportError.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from probe_designer.qc.cross_ligation import ProbeForScreen, build_v2_geom

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
# Public dataclass
# ----------------------------------------------------------------------


@dataclass
class TernaryHit:
    """One NUPACK ternary complex evaluation for an (A-ligator, B-splint) pair."""
    seq_a_id: str
    seq_b_id: str
    a_target: str = ""
    b_target: str = ""
    # Concatenation order in dot-bracket: arm3 + splint + arm5
    arm3_len: int = 0
    splint_len: int = 0
    arm5_len: int = 0
    # NUPACK results
    ternary_dg_kcal: float = 0.0          # full ternary ΔG at celsius (incl. coaxial)
    ternary_concentration_m: float = 0.0  # equilibrium [arm3·splint·arm5]
    ternary_fraction_of_b: float = 0.0    # ternary / total splint
    mfe_dotbracket: str = ""              # dot-bracket with '+' between strands
    mfe_nick_adjacent: bool = False       # arm3 3'-OH + arm5 5'-P pair on adjacent splint bases
    mfe_vicinity_contiguous: bool = False # per-arm ±n bases contiguous from terminus
    b_3oh_pos: Optional[int] = None       # splint canonical 5'→3' partner of arm3 last base
    b_5p_pos: Optional[int] = None        # splint canonical 5'→3' partner of arm5 first base
    # Audit / reproducibility
    ensemble_used: str = "stacking"
    celsius_used: float = 55.0
    sodium_m: float = 0.075
    magnesium_m: float = 0.010
    strand_conc_m: float = 1.0e-7
    vicinity_n_each_side: int = 3


# ----------------------------------------------------------------------
# Helpers: dot-bracket parsing for 3-strand concatenated structure
# ----------------------------------------------------------------------


def _parse_dotbracket_pairs(dotbracket: str) -> Dict[int, int]:
    """Parse a dot-bracket string (with '+' strand separators) into a
    base-pair dict ``{pos_i: pos_j}`` where positions are CONCATENATED
    coordinates (ignoring '+' chars).
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
                # Malformed — silently skip
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
        arm3's 3'-OH = concat pos (arm3_len - 1)
        arm5's 5'-P  = concat pos (arm3_len + splint_len)
    Returns ``(can_ligate, b_3oh, b_5p)`` where b_3oh / b_5p are splint
    canonical 5'→3' positions (concat partner minus arm3_len). can_ligate
    requires both ends pair to splint AND b_3oh - b_5p == 1.
    """
    arm3_3oh_concat = arm3_len - 1
    arm5_5p_concat = arm3_len + splint_len

    p3 = pairs.get(arm3_3oh_concat)
    p5 = pairs.get(arm5_5p_concat)
    if p3 is None or p5 is None:
        return False, None, None
    # Both partners must fall inside the splint range
    if not (arm3_len <= p3 < arm3_len + splint_len):
        return False, None, None
    if not (arm3_len <= p5 < arm3_len + splint_len):
        return False, None, None
    b_3oh = p3 - arm3_len
    b_5p = p5 - arm3_len
    return (b_3oh - b_5p == 1), b_3oh, b_5p


def _check_mfe_vicinity_contiguous(
    pairs: Dict[int, int], arm3_len: int, splint_len: int, arm5_len: int,
    n_each_side: int,
) -> bool:
    """Per-arm vicinity contiguity on the MFE ternary structure.

    arm3 (concat positions ``[0 .. arm3_len-1]``): last n+1 positions
    must all pair to splint AND splint partners form strict descending run.
    arm5 (concat positions ``[arm3_len+splint_len .. end]``): first n+1
    positions same requirement.
    """
    if n_each_side <= 0:
        return False
    if arm3_len <= n_each_side or arm5_len <= n_each_side:
        return False

    # arm3: positions [arm3_len - 1 - n .. arm3_len - 1] in concat coords
    arm3_concat_positions = list(range(arm3_len - 1 - n_each_side, arm3_len))
    for a in arm3_concat_positions:
        p = pairs.get(a)
        if p is None or not (arm3_len <= p < arm3_len + splint_len):
            return False
    arm3_b_seq = [pairs[a] - arm3_len for a in arm3_concat_positions]
    for i in range(len(arm3_b_seq) - 1):
        if arm3_b_seq[i + 1] - arm3_b_seq[i] != -1:
            return False

    # arm5: positions [arm3_len+splint_len .. arm3_len+splint_len + n] in concat
    base = arm3_len + splint_len
    arm5_concat_positions = list(range(base, base + n_each_side + 1))
    for a in arm5_concat_positions:
        p = pairs.get(a)
        if p is None or not (arm3_len <= p < arm3_len + splint_len):
            return False
    arm5_b_seq = [pairs[a] - arm3_len for a in arm5_concat_positions]
    for i in range(len(arm5_b_seq) - 1):
        if arm5_b_seq[i + 1] - arm5_b_seq[i] != -1:
            return False

    return True


# ----------------------------------------------------------------------
# Core: single-direction NUPACK ternary call
# ----------------------------------------------------------------------


def _evaluate_ternary_directional(
    *,
    ligator: ProbeForScreen, splint: ProbeForScreen,
    sodium_m: float, magnesium_m: float, celsius: float,
    ensemble: str, strand_conc_m: float, vicinity_n_each_side: int,
) -> Optional[TernaryHit]:
    """Run a single 3-strand NUPACK ternary analysis (ligator's arm3 +
    ligator's arm5 + splint's full sequence). Returns None if NUPACK
    finds no productive ternary complex (very rare).
    """
    try:
        import nupack
    except ImportError as exc:
        raise ImportError(
            "NUPACK 4 is required for ext.nupack.check.screen_ternary. "
            "Install via academic registration at https://nupack.org/ "
            "and `pip install nupack-x.x.tar.gz` into the probe-design env."
        ) from exc

    arm3_eff, arm5, splint_seq = build_v2_geom(ligator)
    _, _, splint_b_seq = build_v2_geom(splint)
    assert splint_seq != splint_b_seq or splint.probe_id == ligator.probe_id, (
        "build_v2_geom returned different sequences for two probes — should be impossible"
    )
    splint_seq = splint_b_seq   # B's sequence is the splint

    # Strands. NUPACK strand names appear in the dot-bracket order — we put
    # arm3 first, splint second, arm5 third.
    s_arm3 = nupack.Strand(arm3_eff, name="arm3")
    s_splint = nupack.Strand(splint_seq, name="splint")
    s_arm5 = nupack.Strand(arm5, name="arm5")

    tube = nupack.Tube(
        strands={s_arm3: strand_conc_m, s_splint: strand_conc_m, s_arm5: strand_conc_m},
        complexes=nupack.SetSpec(max_size=3),
        name="cross_lig_tube",
    )
    model = nupack.Model(
        material="dna",
        celsius=celsius,
        sodium=sodium_m,
        magnesium=magnesium_m,
        ensemble=ensemble,
    )

    result = nupack.tube_analysis(
        tubes=[tube], model=model,
        compute=["mfe"],
    )

    # Find the productive ternary complex (3-strand, contains arm3 + splint + arm5)
    target_strands = (s_arm3, s_splint, s_arm5)
    productive_complex = None
    productive_conc = 0.0
    productive_dg = 0.0
    for cplx, data in result.tubes[tube].complex_concentrations.items():
        strands_in_cplx = tuple(cplx.strands)
        if (set(strands_in_cplx) == set(target_strands)
                and len(strands_in_cplx) == 3):
            productive_complex = cplx
            productive_conc = float(data)
            productive_dg = float(result.complexes[cplx].free_energy)
            break

    if productive_complex is None:
        # No productive ternary enumerated (essentially never happens with max_size=3)
        return None

    # Compute splint total concentration (sum over all complexes containing splint)
    splint_total = 0.0
    for cplx, data in result.tubes[tube].complex_concentrations.items():
        if s_splint in cplx.strands:
            splint_total += float(data) * cplx.strands.count(s_splint)
    fraction_of_b = (productive_conc / splint_total) if splint_total > 0 else 0.0

    # MFE structure of the productive ternary
    mfe_structures = result.complexes[productive_complex].mfe
    if not mfe_structures:
        return None
    mfe_struct = mfe_structures[0]
    # Order the strands as (arm3, splint, arm5) in the dot-bracket
    # NUPACK returns the MFE structure in the strand order it stored
    # internally; we reorder by re-querying with explicit ordering.
    # The simplest portable approach: get structure as str (which respects
    # the strand order in the Complex), then parse pair indices in that order.
    mfe_dotbracket = str(mfe_struct.structure)

    # The complex's strands attribute tells us the actual ordering used by NUPACK.
    # Map strand-name → known sequence length so we can compute concat offsets
    # without relying on str(Strand) returning the sequence.
    name_to_len = {
        "arm3": len(arm3_eff),
        "splint": len(splint_seq),
        "arm5": len(arm5),
    }
    cplx_strand_order = list(productive_complex.strands)
    strand_lens = [name_to_len[s.name] for s in cplx_strand_order]
    offsets = [0]
    for L in strand_lens:
        offsets.append(offsets[-1] + L)
    # Locate each named strand's offset by name (independent of NUPACK ordering)
    name_offsets = {
        s.name: offsets[i] for i, s in enumerate(cplx_strand_order)
    }
    arm3_offset = name_offsets["arm3"]
    splint_offset = name_offsets["splint"]
    arm5_offset = name_offsets["arm5"]

    # Parse dot-bracket pairs in the complex's ordering
    pairs_concat = _parse_dotbracket_pairs(mfe_dotbracket)
    # Translate to (arm3 | splint | arm5) canonical ordering by remapping concat positions
    # We need pairs as if order were (arm3, splint, arm5).
    # Build a remap: original_concat_pos → canonical_concat_pos
    remap: Dict[int, int] = {}
    # arm3 chunk
    for i in range(len(arm3_eff)):
        remap[arm3_offset + i] = i
    # splint chunk
    for i in range(len(splint_seq)):
        remap[splint_offset + i] = len(arm3_eff) + i
    # arm5 chunk
    for i in range(len(arm5)):
        remap[arm5_offset + i] = len(arm3_eff) + len(splint_seq) + i

    pairs: Dict[int, int] = {}
    for k, v in pairs_concat.items():
        rk = remap.get(k); rv = remap.get(v)
        if rk is not None and rv is not None:
            pairs[rk] = rv

    can_ligate, b_3oh, b_5p = _check_mfe_nick_geometry(
        pairs, len(arm3_eff), len(splint_seq),
    )
    vicinity_ok = (
        _check_mfe_vicinity_contiguous(
            pairs, len(arm3_eff), len(splint_seq), len(arm5), vicinity_n_each_side,
        )
        if can_ligate and vicinity_n_each_side > 0 else False
    )

    return TernaryHit(
        seq_a_id=ligator.probe_id,
        seq_b_id=splint.probe_id,
        a_target=ligator.target,
        b_target=splint.target,
        arm3_len=len(arm3_eff),
        splint_len=len(splint_seq),
        arm5_len=len(arm5),
        ternary_dg_kcal=productive_dg,
        ternary_concentration_m=productive_conc,
        ternary_fraction_of_b=fraction_of_b,
        mfe_dotbracket=mfe_dotbracket,
        mfe_nick_adjacent=can_ligate,
        mfe_vicinity_contiguous=vicinity_ok,
        b_3oh_pos=b_3oh,
        b_5p_pos=b_5p,
        ensemble_used=ensemble,
        celsius_used=celsius,
        sodium_m=sodium_m,
        magnesium_m=magnesium_m,
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
    on_progress: Optional[Any] = None,
) -> List[TernaryHit]:
    """Run NUPACK 4 ternary-complex screen on (ligator, splint) directional pairs.

    Args:
        probes: full probe list (ProbeForScreen).
        prefilter_pairs: iterable of ``(ligator_id, splint_id)`` DIRECTIONAL
            pairs to evaluate. If None, all directional combinations are
            evaluated (~N×(N-1) calls — VERY slow for large pools).
        sodium_m, magnesium_m, celsius, ensemble, strand_conc_m: NUPACK Model
            and tube parameters (see :mod:`probe_designer.ext.nupack.config`).
        vicinity_n_each_side: per-arm contiguity ±n requirement on MFE structure.
        on_progress: optional callable ``(idx, total)`` invoked after each
            pair evaluation, for batch-run progress tracking.

    Returns a list of :class:`TernaryHit`, one per directional pair where
    the productive ternary was enumerated.

    Raises ``ImportError`` if NUPACK 4 is not installed.
    """
    probe_by_id = {p.probe_id: p for p in probes}

    if prefilter_pairs is None:
        # All directional pairs
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
            continue
        try:
            hit = _evaluate_ternary_directional(
                ligator=ligator, splint=splint,
                sodium_m=sodium_m, magnesium_m=magnesium_m, celsius=celsius,
                ensemble=ensemble, strand_conc_m=strand_conc_m,
                vicinity_n_each_side=vicinity_n_each_side,
            )
        except ImportError:
            # NUPACK missing — fatal for the whole audit, not a per-pair issue
            raise
        except Exception as exc:
            # Single-pair failures (e.g. NUPACK numerical issues on
            # degenerate sequences) should not kill the whole audit.
            print(f"  [ternary skip] {lig_id} x {splint_id}: {exc}")
            hit = None
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
    "mfe_nick_adjacent", "mfe_vicinity_contiguous",
    "b_3oh_pos", "b_5p_pos",
    "ensemble_used", "celsius_used",
    "sodium_m", "magnesium_m", "strand_conc_m",
    "vicinity_n_each_side",
    "mfe_dotbracket",
)


def write_ternary_report(hits: List[TernaryHit], tsv_path: Path | str) -> None:
    """Write a list of TernaryHit to a TSV with v2.3 schema. Sorted by
    ternary_fraction_of_b descending (most likely cross-lig first).
    """
    tsv_path = Path(tsv_path)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    sorted_hits = sorted(hits, key=lambda h: h.ternary_fraction_of_b, reverse=True)

    def _esc(s: str) -> str:
        return (s or "").replace("\n", "\\n").replace("\t", "\\t")

    with open(tsv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(list(REPORT_COLUMNS))
        for h in sorted_hits:
            b3 = h.b_3oh_pos if h.b_3oh_pos is not None else ""
            b5 = h.b_5p_pos if h.b_5p_pos is not None else ""
            w.writerow([
                h.seq_a_id, h.seq_b_id, h.a_target, h.b_target,
                h.arm3_len, h.splint_len, h.arm5_len,
                f"{h.ternary_dg_kcal:.3f}",
                f"{h.ternary_concentration_m:.3e}",
                f"{h.ternary_fraction_of_b:.3e}",
                h.mfe_nick_adjacent, h.mfe_vicinity_contiguous,
                b3, b5,
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
    "REPORT_COLUMNS",
]
