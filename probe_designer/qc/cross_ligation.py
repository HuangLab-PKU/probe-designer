"""v2.2 cross-ligation screen — split-fragment ligation geometry.

**Physical model (corrected from v2 / v2.1).**

Cross-ligation between two padlock probes occurs when:

1. Probe A's arm3 fragment (carrying the real 3'-OH) and arm5 fragment
   (carrying the real 5'-P) **independently** hybridize to probe B's
   sequence (acting as splint).
2. Their nick-adjacent ends meet on B at consecutive positions:
   ``b_3oh - b_5p == 1`` (B's 5'→3' index decreases by 1 across the nick
   because both arms anneal antiparallel to B).
3. SplintR-style ligase clamps the nick if WC contiguity holds for the
   ~3 bases on each side of the nick (see [[reference-splintr-fidelity]]).

**Why v2.1's `arm3 + arm5` rotation was wrong.** Rotating arm3 + arm5
into one continuous 40-nt oligo models the **product** of ligation (the
closed padlock circle), not the reactant. The backbone tether between
arm3 and arm5 doesn't participate in WC pairing — primer3 falsely
treats them as covalently connected, producing dimer geometries that
can't actually form before ligation. v2.2 runs primer3 twice, once
per fragment.

**Tier 1 (k-mer prefilter).** For each direction (A as ligator on B as
splint), check that A's arm3 has at least one k-mer (window ±n nt
across arm3's 3'-end) whose reverse complement appears in B's full
k-mer set AND A's arm5 has at least one such k-mer (across arm5's
5'-end). Both conditions necessary — the nick joins arm3-end with
arm5-start.

**Tier 2 (primer3 ×2 + nick check).** For each Tier 1 hit, run two
``primer3.calc_heterodimer`` calls:
  * ``calc_heterodimer(arm3, B.sequence, output_structure=True)``
  * ``calc_heterodimer(arm5, B.sequence, output_structure=True)``
Parse each ASCII with :func:`parse_primer3_dimer_pairing`. From arm3's
ASCII, read ``b_3oh = pair_positions[len(arm3_effective) - 1]``. From
arm5's ASCII, read ``b_5p = pair_positions[0]``. Ligation-competent
geometry requires both keys exist AND ``b_3oh - b_5p == 1``.

**Vicinity (per-arm, ±n=3).** SplintR clamps ~3 bases each side of the
nick. Translated to v2.2 split-fragment: arm3's last n+1 positions
``[len(arm3_eff)-1-n .. len(arm3_eff)-1]`` must all be paired in arm3's
ASCII and their B partners must form a strictly descending antiparallel
run (``b[i+1] = b[i] - 1``). Same for arm5's first n+1 positions
``[0 .. n]``. (Note: vicinity is **per-arm**, not across the nick — the
nick itself is the inter-arm gap and is by construction OK once
``b_3oh - b_5p == 1`` is satisfied.)

**iLock chemistry.** iLock probes lose 1 nt at the 3' end of arm3 after
Taq FEN cleavage (the "overlap" base). ``build_v2_geom`` returns
``arm3_effective = arm3[1:]`` for iLock; that's the actual ligation
substrate. Same logic applies — arm3_effective is the fragment we feed
to primer3.

**Architecture**: this module is BANK-FREE. Pool-loading is the
responsibility of ``probe_designer.ext.pool.loader``. The static check
at ``designer/tests/unit/test_no_bank_import_in_qc.py`` enforces this.

See ``[[project-cross-lig-v2-2]]`` for the design rationale and the
2026-05-22 TNBC audit comparison.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any, List, Optional, Set, Tuple

from probe_designer.chemistry import ReactionConditions
from probe_designer.qc.dimer_ascii import parse_primer3_dimer_pairing


# Buffer parameters DERIVED from the single source of truth (ReactionConditions
# defaults = protocol rca.md v5.3), so the cross-lig screen and the Tm path never
# drift apart (audit P5). ReactionConditions is bank-free, keeping qc bank-free.
#   - 75 mM monovalent (K+), 10 mM Mg2+, 0.1 μM probe, 20% formamide → 55 °C
#     effective in the no-formamide NN model (45 °C lab + 0.5 °C/%FA × 20%).
_LAB = ReactionConditions()
DEFAULT_K: int = 4
DEFAULT_JUNCTION_WINDOW: int = 3
DEFAULT_TM_THRESHOLD_C: float = 27.0
DEFAULT_MV_CONC_MM: float = _LAB.monovalent_mM      # 75 mM K+
DEFAULT_DV_CONC_MM: float = _LAB.mg_mM              # 10 mM Mg2+
DEFAULT_DNA_CONC_M: float = _LAB.strand_nM * 1e-9   # 0.1 μM probe
DEFAULT_TEMP_C: float = _LAB.effective_celsius      # 55 °C formamide-effective

# Vicinity contiguity requirement around the ligation nick. SplintR / PBCV-1
# ligase's active site clamps ~3 nt on each side of the nick (contiguous
# Watson-Crick, no bulge). See [[reference-splintr-fidelity]] for the
# Lohman 2014 / Krzywkowski 2017 evidence. Default 3 nt each side.
# Setting to 0 disables the vicinity check.
DEFAULT_VICINITY_N: int = 3


# ----------------------------------------------------------------------
# Public dataclasses
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class ProbeForScreen:
    """Minimal probe protocol the screen consumes — compute input only.

    Converted from bank's ``EffectiveProbe`` at the ext/pool boundary; the
    compute layer never sees bank's loader types.
    """
    probe_id: str
    chemistry: str               # "iLock" | "dRNA" | "cDNA"
    probe_arm5: str
    probe_arm3: str
    sequence: str                # full assembled probe (case-insensitive)
    target: str = ""


@dataclass
class JunctionSeedHitV2:
    """Tier 1 k-mer seed hit — one direction (a as ligator, b as splint).

    A pair passes Tier 1 only when BOTH arms have a junction-spanning
    k-mer match on B (the nick joins arm3-end with arm5-start, so both
    fragments must template).
    """
    probe_a_id: str              # ligator
    probe_b_id: str              # splint
    direction: str = "a_lig_on_b"
    arm3_seed: str = ""          # k-mer near arm3's 3'-end that hit B
    arm5_seed: str = ""          # k-mer near arm5's 5'-end that hit B


@dataclass
class LigationDimer:
    """One Tier 2 dimer call (one direction: A-as-ligator on B-as-splint).

    Two independent primer3 alignments: arm3 vs B, arm5 vs B. Captured
    as two separate ASCII strings + Tm/dG values.
    """
    seq_a_id: str                # ligator (A)
    seq_b_id: str                # splint  (B)
    a_target: str = ""
    b_target: str = ""
    # Per-arm primer3 outputs
    arm3_tm_c: float = 0.0
    arm3_dg_kcal: float = 0.0
    arm3_dimer_structure: str = ""
    arm5_tm_c: float = 0.0
    arm5_dg_kcal: float = 0.0
    arm5_dimer_structure: str = ""
    # Combined / derived
    overall_tm_c: float = 0.0           # max(arm3_tm, arm5_tm) — used for Tm threshold
    flagged_overall: bool = False       # max Tm > threshold
    a_can_ligate_on_b: bool = False     # both ends paired AND nick-adjacent on B
    b_3oh_pos: Optional[int] = None     # B 5'→3' partner of arm3's 3'-OH
    b_5p_pos: Optional[int] = None      # B 5'→3' partner of arm5's 5'-P
    vicinity_contiguous: bool = False   # per-arm vicinity contiguous on B
    vicinity_n_each_side: int = 0       # n used; 0 means check disabled


# Internal seed-set type
@dataclass(frozen=True)
class _SeedSetV22:
    probe_id: str
    chemistry: str
    arm3_effective: str          # arm3[1:] for iLock, arm3 for others
    arm5: str
    splint_seq: str              # probe.sequence (full assembled, w/ bb)
    arm3_end_seeds: Set[str]     # k-mers spanning arm3's 3' end
    arm5_start_seeds: Set[str]   # k-mers spanning arm5's 5' end
    full_kmers: Set[str]         # all k-mers of splint_seq (RC matching)
    target: str = ""


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "U": "A", "N": "N"}


def _rc(seq: str) -> str:
    return "".join(_COMP.get(b, "N") for b in reversed(seq.upper()))


def _kmers(seq: str, k: int) -> Set[str]:
    seq = seq.upper()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def _extract_end_kmers(seq: str, end: str, k: int, window: int) -> Set[str]:
    """Return all k-mers covering the specified end of seq, ±window nt.

    ``end="3"``: k-mers that cover seq's 3'-terminal base (the last char).
        Returns k-mers starting at positions [len-k-window, ..., len-k].
    ``end="5"``: k-mers that cover seq's 5'-terminal base (the first char).
        Returns k-mers starting at positions [0, ..., window].
    """
    seq = seq.upper()
    L = len(seq)
    out: Set[str] = set()
    if end == "3":
        for start in range(max(0, L - k - window), L - k + 1):
            if start + k <= L:
                out.add(seq[start:start + k])
    elif end == "5":
        for start in range(0, min(L - k, window) + 1):
            out.add(seq[start:start + k])
    return out


def build_v2_geom(probe: ProbeForScreen) -> Tuple[str, str, str]:
    """Compute the v2.2 split-fragment geometry for one probe.

    Returns ``(arm3_effective, arm5, splint_seq)``:
      * ``arm3_effective`` — arm3 with iLock's 1-nt overlap dropped.
        This is the fragment carrying the real 3'-OH; its **last** char is
        the 3'-OH-bearing base.
      * ``arm5`` — the fragment carrying the real 5'-P; its **first** char
        is the 5'-P-bearing base.
      * ``splint_seq`` — the probe's full assembled sequence (backbone
        included). When this probe acts as splint for another probe's
        ligation, this is the template.

    Raises ``ValueError`` on registry inconsistency (arm not in sequence).
    """
    is_ilock = probe.chemistry == "iLock"
    arm5 = probe.probe_arm5.upper()
    arm3 = probe.probe_arm3.upper()
    seq_upper = probe.sequence.upper()

    if not seq_upper:
        raise ValueError(
            f"probe {probe.probe_id}: empty assembled sequence — registry inconsistent"
        )
    if arm5 not in seq_upper:
        raise ValueError(
            f"probe {probe.probe_id}: arm5 ({arm5[:10]}...) not in assembled sequence — "
            f"registry inconsistent"
        )
    if arm3 not in seq_upper:
        raise ValueError(
            f"probe {probe.probe_id}: arm3 ({arm3[:10]}...) not in assembled sequence — "
            f"registry inconsistent"
        )

    arm3_effective = arm3[1:] if is_ilock else arm3
    return arm3_effective, arm5, seq_upper


def _build_seed_set(probe: ProbeForScreen, *, k: int, window: int) -> _SeedSetV22:
    """Construct per-arm k-mer seeds + splint k-mer index for one probe."""
    arm3_effective, arm5, splint_seq = build_v2_geom(probe)
    return _SeedSetV22(
        probe_id=probe.probe_id,
        chemistry=probe.chemistry,
        arm3_effective=arm3_effective,
        arm5=arm5,
        splint_seq=splint_seq,
        arm3_end_seeds=_extract_end_kmers(arm3_effective, end="3", k=k, window=window),
        arm5_start_seeds=_extract_end_kmers(arm5, end="5", k=k, window=window),
        full_kmers=_kmers(splint_seq, k),
        target=probe.target,
    )


def _tier1_directional_seed_hit(
    ligator: _SeedSetV22, splint: _SeedSetV22,
) -> Optional[JunctionSeedHitV2]:
    """Check if BOTH ligator's arm3-end AND arm5-start have RC matches in splint.

    The nick joins arm3's 3'-end with arm5's 5'-end; for cross-lig to be
    possible, both fragments must be templated by the splint. If only one
    arm has a match, ligation can't proceed regardless of geometry.
    """
    arm3_hit_seed = ""
    for seed in ligator.arm3_end_seeds:
        if _rc(seed) in splint.full_kmers:
            arm3_hit_seed = seed
            break
    if not arm3_hit_seed:
        return None
    arm5_hit_seed = ""
    for seed in ligator.arm5_start_seeds:
        if _rc(seed) in splint.full_kmers:
            arm5_hit_seed = seed
            break
    if not arm5_hit_seed:
        return None
    return JunctionSeedHitV2(
        probe_a_id=ligator.probe_id,
        probe_b_id=splint.probe_id,
        direction="a_lig_on_b",
        arm3_seed=arm3_hit_seed,
        arm5_seed=arm5_hit_seed,
    )


def _check_nick_geometry(
    arm3_ascii: str, arm5_ascii: str, arm3_eff_len: int,
) -> Tuple[bool, Optional[int], Optional[int]]:
    """Parse the two ASCII alignments and check nick-adjacency on B.

    Returns ``(can_ligate, b_3oh, b_5p)``:
      * ``b_3oh`` — B 5'→3' partner of arm3's last base (arm3 pos len-1),
        or None if arm3's 3'-OH is unpaired.
      * ``b_5p`` — B 5'→3' partner of arm5's first base (arm5 pos 0),
        or None if arm5's 5'-P is unpaired.
      * ``can_ligate`` — both ends paired AND ``b_3oh - b_5p == 1``
        (B's 5'→3' index decreases by 1 across the nick because both
        arms anneal antiparallel to B).
    """
    arm3_pairs = parse_primer3_dimer_pairing(arm3_ascii)["pair_positions_a_to_b"]
    arm5_pairs = parse_primer3_dimer_pairing(arm5_ascii)["pair_positions_a_to_b"]

    arm3_end_pos = arm3_eff_len - 1   # arm3's 3'-OH = last base of arm3_effective
    arm5_start_pos = 0                # arm5's 5'-P = first base of arm5

    b_3oh = arm3_pairs.get(arm3_end_pos)
    b_5p = arm5_pairs.get(arm5_start_pos)
    if b_3oh is None or b_5p is None:
        return False, b_3oh, b_5p
    return (b_3oh - b_5p == 1), b_3oh, b_5p


def _check_per_arm_vicinity(
    arm3_ascii: str, arm5_ascii: str,
    arm3_eff_len: int, arm5_len: int,
    n_each_side: int,
) -> bool:
    """Check vicinity contiguity inside each arm independently.

    arm3's last ``n+1`` positions (``[arm3_eff_len-1-n .. arm3_eff_len-1]``)
    must all be paired and their B partners must form a strict descending
    antiparallel run. Same for arm5's first ``n+1`` positions
    (``[0 .. n]``).

    Returns False when ``n_each_side <= 0`` (documented disable value).
    """
    if n_each_side <= 0:
        return False
    arm3_pairs = parse_primer3_dimer_pairing(arm3_ascii)["pair_positions_a_to_b"]
    arm5_pairs = parse_primer3_dimer_pairing(arm5_ascii)["pair_positions_a_to_b"]

    # arm3: positions [L-1-n, L-1-n+1, ..., L-1] in 5'→3' (n+1 positions)
    arm3_positions = list(range(arm3_eff_len - 1 - n_each_side, arm3_eff_len))
    if any(a not in arm3_pairs for a in arm3_positions):
        return False
    arm3_b_seq = [arm3_pairs[a] for a in arm3_positions]
    for i in range(len(arm3_b_seq) - 1):
        if arm3_b_seq[i + 1] - arm3_b_seq[i] != -1:
            return False

    # arm5: positions [0, 1, ..., n]
    arm5_positions = list(range(0, n_each_side + 1))
    if arm5_len <= n_each_side:
        return False
    if any(a not in arm5_pairs for a in arm5_positions):
        return False
    arm5_b_seq = [arm5_pairs[a] for a in arm5_positions]
    for i in range(len(arm5_b_seq) - 1):
        if arm5_b_seq[i + 1] - arm5_b_seq[i] != -1:
            return False

    return True


def _calc_dimer(
    ligator: _SeedSetV22, splint: _SeedSetV22, *,
    tm_threshold_c: float,
    mv_conc_mm: float, dv_conc_mm: float, dna_conc_m: float, temp_c: float,
    vicinity_n_each_side: int,
) -> Optional[LigationDimer]:
    """Run primer3 twice (arm3 vs B, arm5 vs B), parse, judge."""
    import primer3   # lazy

    buffer = dict(
        mv_conc=mv_conc_mm, dv_conc=dv_conc_mm,
        dna_conc=dna_conc_m * 1e9,   # primer3 expects nM
        temp_c=temp_c,
    )
    arm3_thermo = primer3.calc_heterodimer(
        ligator.arm3_effective, splint.splint_seq,
        output_structure=True, **buffer,
    )
    arm5_thermo = primer3.calc_heterodimer(
        ligator.arm5, splint.splint_seq,
        output_structure=True, **buffer,
    )
    if not (arm3_thermo.structure_found or arm5_thermo.structure_found):
        return None

    arm3_ascii = arm3_thermo.ascii_structure or "" if arm3_thermo.structure_found else ""
    arm5_ascii = arm5_thermo.ascii_structure or "" if arm5_thermo.structure_found else ""

    can_ligate, b_3oh, b_5p = _check_nick_geometry(
        arm3_ascii, arm5_ascii, len(ligator.arm3_effective),
    )
    vicinity_ok = (
        _check_per_arm_vicinity(
            arm3_ascii, arm5_ascii,
            len(ligator.arm3_effective), len(ligator.arm5),
            vicinity_n_each_side,
        )
        if (can_ligate and vicinity_n_each_side > 0) else False
    )

    arm3_tm = float(arm3_thermo.tm) if arm3_thermo.structure_found else 0.0
    arm5_tm = float(arm5_thermo.tm) if arm5_thermo.structure_found else 0.0
    overall_tm = max(arm3_tm, arm5_tm)

    return LigationDimer(
        seq_a_id=ligator.probe_id,
        seq_b_id=splint.probe_id,
        a_target=ligator.target,
        b_target=splint.target,
        arm3_tm_c=arm3_tm,
        arm3_dg_kcal=(float(arm3_thermo.dg) / 1000.0) if arm3_thermo.structure_found else 0.0,
        arm3_dimer_structure=arm3_ascii,
        arm5_tm_c=arm5_tm,
        arm5_dg_kcal=(float(arm5_thermo.dg) / 1000.0) if arm5_thermo.structure_found else 0.0,
        arm5_dimer_structure=arm5_ascii,
        overall_tm_c=overall_tm,
        flagged_overall=overall_tm > tm_threshold_c,
        a_can_ligate_on_b=can_ligate,
        b_3oh_pos=b_3oh,
        b_5p_pos=b_5p,
        vicinity_contiguous=vicinity_ok,
        vicinity_n_each_side=vicinity_n_each_side,
    )


# ----------------------------------------------------------------------
# Public entry point
# ----------------------------------------------------------------------


def screen_cross_ligation_v2(
    probes: List[ProbeForScreen], *,
    k: int = DEFAULT_K,
    window: int = DEFAULT_JUNCTION_WINDOW,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    mv_conc_mm: float = DEFAULT_MV_CONC_MM,
    dv_conc_mm: float = DEFAULT_DV_CONC_MM,
    dna_conc_m: float = DEFAULT_DNA_CONC_M,
    temp_c: float = DEFAULT_TEMP_C,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
) -> Tuple[List[JunctionSeedHitV2], List[LigationDimer]]:
    """v2.2 split-fragment cross-ligation screen.

    For each unordered pair (a, b), screen in BOTH directions
    (a-as-ligator on b-as-splint; b-as-ligator on a-as-splint). Each
    direction emits at most one :class:`LigationDimer`. A pair is
    "confirmed cross-lig risk" iff:

        flagged_overall AND a_can_ligate_on_b AND vicinity_contiguous

    Args:
        vicinity_n_each_side: per-arm contiguity requirement (see
            :data:`DEFAULT_VICINITY_N`). Set to 0 to disable the check.

    Returns ``(tier1_hits, all_dimers)``:
      * ``tier1_hits``: every k-mer-seed match where BOTH arm3-end and
        arm5-start hit the splint.
      * ``all_dimers``: every primer3-computed dimer that passed Tier 1.

    Raises ``ValueError`` (during seed-set construction) on registry
    inconsistency.
    """
    if not probes:
        return [], []

    seed_index: dict[str, _SeedSetV22] = {
        p.probe_id: _build_seed_set(p, k=k, window=window) for p in probes
    }

    tier1: List[JunctionSeedHitV2] = []
    dimers: List[LigationDimer] = []

    for a, b in combinations(probes, 2):
        sa = seed_index[a.probe_id]
        sb = seed_index[b.probe_id]

        # Direction 1: a as ligator on b as splint
        seed_hit_ab = _tier1_directional_seed_hit(sa, sb)
        if seed_hit_ab is not None:
            tier1.append(seed_hit_ab)
            d = _calc_dimer(
                sa, sb,
                tm_threshold_c=tm_threshold_c,
                mv_conc_mm=mv_conc_mm, dv_conc_mm=dv_conc_mm,
                dna_conc_m=dna_conc_m, temp_c=temp_c,
                vicinity_n_each_side=vicinity_n_each_side,
            )
            if d is not None:
                dimers.append(d)

        # Direction 2: b as ligator on a as splint
        seed_hit_ba = _tier1_directional_seed_hit(sb, sa)
        if seed_hit_ba is not None:
            tier1.append(seed_hit_ba)
            d = _calc_dimer(
                sb, sa,
                tm_threshold_c=tm_threshold_c,
                mv_conc_mm=mv_conc_mm, dv_conc_mm=dv_conc_mm,
                dna_conc_m=dna_conc_m, temp_c=temp_c,
                vicinity_n_each_side=vicinity_n_each_side,
            )
            if d is not None:
                dimers.append(d)

    return tier1, dimers


# ----------------------------------------------------------------------
# Report writer
# ----------------------------------------------------------------------


REPORT_COLUMNS: Tuple[str, ...] = (
    "seq_a_id", "seq_b_id", "a_target", "b_target",
    "overall_tm_c",
    "arm3_tm_c", "arm3_dg_kcal",
    "arm5_tm_c", "arm5_dg_kcal",
    "flagged_overall", "a_can_ligate_on_b",
    "vicinity_contiguous", "vicinity_n_each_side",
    "b_3oh_pos", "b_5p_pos",
    "arm3_dimer_structure", "arm5_dimer_structure",
)


def write_dimer_report(dimers: List[LigationDimer], tsv_path: Path | str) -> None:
    """Write a list of :class:`LigationDimer` to a TSV with v2.2 schema.

    Sorted by overall_tm_c descending. Empty list still writes a
    header-only file. Newlines inside the two ASCII fields are escaped
    so each row stays on one TSV line.
    """
    tsv_path = Path(tsv_path)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    sorted_dimers = sorted(dimers, key=lambda d: d.overall_tm_c, reverse=True)

    def _esc(s: str) -> str:
        return (s or "").replace("\n", "\\n").replace("\t", "\\t")

    with open(tsv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(list(REPORT_COLUMNS))
        for d in sorted_dimers:
            b3 = d.b_3oh_pos if d.b_3oh_pos is not None else ""
            b5 = d.b_5p_pos if d.b_5p_pos is not None else ""
            w.writerow([
                d.seq_a_id, d.seq_b_id, d.a_target, d.b_target,
                f"{d.overall_tm_c:.3f}",
                f"{d.arm3_tm_c:.3f}", f"{d.arm3_dg_kcal:.3f}",
                f"{d.arm5_tm_c:.3f}", f"{d.arm5_dg_kcal:.3f}",
                d.flagged_overall, d.a_can_ligate_on_b,
                d.vicinity_contiguous, d.vicinity_n_each_side,
                b3, b5,
                _esc(d.arm3_dimer_structure), _esc(d.arm5_dimer_structure),
            ])


__all__ = [
    "ProbeForScreen", "JunctionSeedHitV2", "LigationDimer",
    "screen_cross_ligation_v2",
    "build_v2_geom",
    "write_dimer_report",
    "DEFAULT_K", "DEFAULT_JUNCTION_WINDOW",
    "DEFAULT_TM_THRESHOLD_C",
    "DEFAULT_MV_CONC_MM", "DEFAULT_DV_CONC_MM",
    "DEFAULT_DNA_CONC_M", "DEFAULT_TEMP_C",
    "DEFAULT_VICINITY_N",
]
