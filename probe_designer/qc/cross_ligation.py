"""v3 cross-ligation screen — enumerate ligation registers, don't take an argmax.

**What the screen is for.** Probe A's two arms can anneal to probe B instead of
to mRNA. If they land nick-adjacent on B, ligase closes **A into a real circle**
carrying A's own barcode, with no target RNA involved. It amplifies and decodes
as gene A everywhere, and being circular it survives the exonuclease step. 10x
describes the same artifact for Xenium panels.

**Why v2.2 was replaced (audit 2026-09-02).** v2.2 called
``primer3.calc_heterodimer`` once per arm against the splint and asked whether
the two independently-chosen best alignments happened to land nick-adjacent.
That is an argmax coincidence test, not a search. Measured on
``TNBCmarker_main_v1``: arm3's 3'-OH was paired in 16.5% of dimers, arm5's 5'-P
in 14.7%, and *both* in 2.4% — the product of the two, i.e. the two alignments
were uncorrelated events rather than one physical complex. Only 3 of 37 187 came
out adjacent, about the rate chance alone predicts. The screen missed the pool's
strongest ligation-competent pair (CD3E on CD3D, 11 bp contiguous across the
nick) while reporting the second-ranked one.

The root cause is that the productive complex's stability is a property of the
two arms *together* — one duplex spanning the nick — but v2.2 scored each arm
alone, then demanded each weak half win its own global argmax. A 12 bp half
loses to any other >=12 bp match on the splint, and since every probe in a pool
carries the same backbone, such matches are the norm.

**The v3 model.** The nick position *determines* the register. With the nick
sitting between splint bases ``j`` and ``j+1``:

* arm3 occupies ``B[j+1 .. j+L3]`` — its 3'-OH pairs ``B[j+1]``
* arm5 occupies ``B[j-L5+1 .. j]`` — its 5'-P pairs ``B[j]``

Both arms run antiparallel to B, so moving 5'->3' along the probe moves 3'->5'
along B; hence ``b_3oh - b_5p == +1``. The two windows are contiguous, so the
rotated ligator ``arm3_effective + arm5`` maps onto one contiguous splint
window. Enumerating ``j`` is O(len B) per direction.

**Tier 1 is now exact, not a heuristic.** The ligase clamp requires ``n+1``
contiguous Watson-Crick pairs on each side of the nick (see
[[reference-splintr-fidelity]]). In the rotated form that is one fixed
``2(n+1)``-nt block straddling the junction, so Tier 1 is a single substring
search for its reverse complement in the splint — no k-mer seeding, no false
negatives, and every hit hands back the exact ``j`` to score.

**Tier 2 scores the fixed alignment.** ``nn_tables.duplex_melting_temperature``
computes the NN Tm of the rotated ligator against the splint window *at that
register*, mismatches included. Because arm3 and arm5 stack coaxially across a
nick, treating the pair as one duplex is the standard approximation; it slightly
overestimates stability, which is the safe direction for a screen.

This also removes v2.1's error without reintroducing it: v2.1 fed the rotated
40-mer to primer3, which was then free to *slide* the alignment so the arm3/arm5
seam sat somewhere other than the nick — geometries that cannot form before
ligation. Fixing the register pins the seam to the nick by construction.

**What v2.2 got right and v3 keeps**: the nick arithmetic, the +/-3 clamp and
its literature basis, the iLock ``arm3[1:]`` substrate, and the bank-free
architecture (pool loading lives in ``probe_designer.ext.pool.loader``; the
static check at ``designer/tests/unit/test_no_bank_import_in_qc.py`` enforces
it).

**Also screened here, and not screened at all before v3**:

* ``A`` templating on *another copy of A* — ``combinations`` had excluded it,
  so the highest-concentration encounter in any pool was never examined.
* **Self-circularisation**: one molecule folding so its own 3'-OH and 5'-P abut,
  needing no splint whatsoever. See :func:`screen_self_circularisation`.

See ``[[project-cross-lig-v3-register-scan]]`` for the audit that produced this.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from probe_designer.chemistry import ReactionConditions, reverse_complement
from probe_designer.nn_tables import DNA_DNA, duplex_melting_temperature
from probe_designer.qc.dimer_ascii import parse_primer3_dimer_pairing


# ----------------------------------------------------------------------
# Reaction conditions
# ----------------------------------------------------------------------

#: The hybridisation the screen models. ``ReactionConditions`` defaults ARE the
#: lab protocol (rca.md v5.3, mirrored in ``configs/config_template.yaml``):
#: 75 mM K+, 10 mM Mg2+, 100 nM strand, 20% formamide, 45 C — 55 C effective in
#: the no-formamide NN model. Pass an explicit ``reaction`` for a different
#: buffer rather than editing this; in an n-plex pool the per-species
#: concentration is the pool stock divided by n, which is the number that
#: belongs in ``strand_nM`` for a bimolecular screen.
LAB_HYBRIDISATION = ReactionConditions()

#: Both strands are DNA here — the splint is another padlock probe, not RNA.
CROSS_LIG_NN_MODEL = DNA_DNA

#: Ligase clamp: n contiguous Watson-Crick pairs required on each side of the
#: nick, so the gate block is 2(n+1) nt wide in the rotated form. SplintR /
#: PBCV-1 clamps ~3 nt each side (Lohman 2014; Krzywkowski 2017) — see
#: [[reference-splintr-fidelity]]. 0 would disable the gate, which makes the
#: screen meaningless; it is rejected rather than silently honoured.
DEFAULT_VICINITY_N: int = 3

#: Flag a register whose duplex Tm clears this. It is an absolute floor, about
#: 28 C below the 55 C effective reaction temperature — deliberately far below
#: it, because ligation is irreversible: a register populated at a small
#: equilibrium fraction is a kinetic trap that accumulates product over an
#: overnight incubation, not a species you can dismiss by its share. (v2.2
#: carried this same 27 C but documented it as "T_rxn - 10" for a 37 C
#: reaction, which stopped being true when the anchor moved to 45 C.)
DEFAULT_TM_THRESHOLD_C: float = 27.0

#: Consecutive mismatches that end the helix. Walking outward from the nick, a
#: run this long stops the duplex: the helix frays there, and it is also where
#: the NN model runs out — Biopython parameterises *single* internal mismatches
#: (DNA_IMM1) but has no term for two in a row, and would drop them from the sum
#: while warning that the result is wrong. Trimming to a core that starts and
#: ends on a Watson-Crick pair keeps every scored neighbour inside the tables.
MAX_CONSECUTIVE_MISMATCH: int = 2

# Buffer scalars, kept for callers that still reference them.
DEFAULT_MV_CONC_MM: float = LAB_HYBRIDISATION.monovalent_mM
DEFAULT_DV_CONC_MM: float = LAB_HYBRIDISATION.mg_mM
DEFAULT_DNA_CONC_M: float = LAB_HYBRIDISATION.strand_nM * 1e-9
DEFAULT_TEMP_C: float = LAB_HYBRIDISATION.effective_celsius


_COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "G": "C", "C": "G"}
_rc = reverse_complement


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


@dataclass(frozen=True)
class JunctionSeedHitV2:
    """Tier 1 hit — one direction (a as ligator, b as splint).

    Unlike the v2.x k-mer seed this is exact: ``junction_block``'s reverse
    complement was found in the splint at ``nick_positions``, so every hit is a
    register that satisfies the ligase clamp and none is missed.
    """
    probe_a_id: str              # ligator
    probe_b_id: str              # splint
    direction: str = "a_lig_on_b"
    junction_block: str = ""     # the 2(n+1)-nt block straddling A's nick
    nick_positions: Tuple[int, ...] = ()   # every j on B satisfying the clamp


@dataclass
class LigationDimer:
    """The best ligation-competent register for one (ligator, splint) direction.

    ``overall_tm_c`` is the Tm of the rotated ligator against the splint window
    *at this register* — one duplex spanning the nick, and the quantity that
    decides whether the complex is populated. It is the only number used as a
    gate.

    The per-arm Tms are **diagnostics**: they split that same alignment at the
    nick, so an asymmetric register (say 8 paired bases on the arm3 side and 4
    on the arm5 side) is visible rather than hidden inside one figure. A half
    of only a few bases melts below 0 C and will report a large negative Tm —
    that is the model saying "this half cannot hold on its own", not an error.
    A half shorter than 2 nt has no nearest neighbour at all and reports 0.0.
    ``limiting_arm_tm_c`` is the ``min`` of the two, since ligation needs BOTH
    arms annealed; v2.2 aggregated with ``max`` and gated on it, which let one
    strong arm carry a pair whose other arm was useless.
    """
    seq_a_id: str                # ligator (A)
    seq_b_id: str                # splint  (B)
    a_target: str = ""
    b_target: str = ""
    # Register
    nick_pos_on_b: int = -1      # nick sits between B[j] and B[j+1]
    b_3oh_pos: Optional[int] = None   # = j + 1
    b_5p_pos: Optional[int] = None    # = j
    junction_run_nt: int = 0     # contiguous WC run straddling the nick
    paired_nt: int = 0           # WC pairs over the whole register
    # Thermodynamics at that register
    overall_tm_c: float = 0.0
    arm3_tm_c: float = 0.0
    arm5_tm_c: float = 0.0
    limiting_arm_tm_c: float = 0.0
    # Verdict
    flagged_overall: bool = False       # overall_tm_c > threshold
    a_can_ligate_on_b: bool = False     # a clamp-satisfying register exists
    vicinity_contiguous: bool = True    # implied by the Tier 1 gate in v3
    vicinity_n_each_side: int = DEFAULT_VICINITY_N
    is_self_pair: bool = False          # ligator and splint are the same probe
    alignment: str = ""                 # rendered register
    n_registers: int = 0                # clamp-satisfying registers considered


@dataclass
class SelfCircularisationHit:
    """One probe folding so its own 3'-OH and 5'-P abut — no splint needed.

    The same register formalism, with the probe's own sequence as splint and
    registers overlapping the arms themselves excluded (a base cannot pair with
    itself). This is target-independent background at its purest: the circle
    forms from a single molecule.
    """
    probe_id: str
    target: str = ""
    nick_pos: int = -1
    junction_run_nt: int = 0
    paired_nt: int = 0
    tm_c: float = 0.0
    flagged: bool = False
    alignment: str = ""


# Internal per-probe geometry
@dataclass(frozen=True)
class _Geometry:
    probe_id: str
    chemistry: str
    arm3_effective: str          # arm3[1:] for iLock, arm3 otherwise
    arm5: str
    ligator: str                 # arm3_effective + arm5, the rotated form
    splint_seq: str              # full assembled probe, when acting as splint
    junction_block: str          # 2(n+1) nt straddling the nick, in ligator
    needle: str                  # reverse complement of junction_block
    arm_spans: Tuple[Tuple[int, int], ...]   # half-open arm spans in splint_seq
    target: str = ""


# ----------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------


def build_v2_geom(probe: ProbeForScreen) -> Tuple[str, str, str]:
    """Compute the ligation-substrate geometry for one probe.

    Returns ``(arm3_effective, arm5_effective, splint_seq)``:
      * ``arm3_effective`` — its **last** base carries the 3'-OH.
      * ``arm5_effective`` — its **first** base carries the 5'-P.
      * ``splint_seq`` — the full assembled sequence, used when this probe acts
        as the splint for another probe's ligation.

    **iLock drops ``arm5[0]``, not ``arm3[0]``.** The invader designer places
    the same SNP-complement base at both arm tips (``arm5[0] == arm3[-1]``, the
    invariant ``ext.mutation.probe.verify_iLock_probe`` enforces) so they
    compete for the target SNP; FEN1 then cleaves the *flap + arm5[0]* block,
    leaving ``arm3[-1]`` as the 3'-OH and ``arm5[1]`` as the new 5'-P. v2.2
    returned ``arm3[1:]`` with arm5 whole, which kept the right 3'-OH base by
    coincidence (``arm3[1:][-1] == arm3[-1]``) but put the 5'-P one base early.
    The junction block then carried the SNP base twice and every iLock register
    was computed out of phase.

    Raises ``ValueError`` on registry inconsistency (arm not in sequence).

    Named for the v2 geometry it introduced; kept because ``ext.nupack.check``
    and ``qc.viz`` import it.
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

    if is_ilock and len(arm5) < 2:
        raise ValueError(
            f"probe {probe.probe_id}: iLock arm5 must be >= 2 nt — FEN1 removes "
            f"arm5[0] with the flap, leaving nothing to ligate"
        )
    arm5_effective = arm5[1:] if is_ilock else arm5
    return arm3, arm5_effective, seq_upper


def ligation_substrate(probe: ProbeForScreen) -> str:
    """The probe as one strand presenting its 5'-P at index 0 and 3'-OH at -1.

    This is the molecule a ligase actually sees: the assembled sequence trimmed
    to start at the effective arm5 and end at arm3's last base. For dRNA/cDNA
    that is the whole probe; for iLock it drops the flap and ``arm5[0]``, the
    block FEN1 cleaves.

    Unlike the split arms, this keeps the **backbone tether** between the two
    ends — the thing that makes a padlock a padlock. Hand this to a folding
    engine (``ext.nupack``) and the loop entropy and the cooperativity between
    the two arms are modelled instead of discarded.
    """
    arm3_eff, arm5_eff, seq = build_v2_geom(probe)
    start = seq.index(arm5_eff)
    end = seq.rindex(arm3_eff) + len(arm3_eff)
    if end <= start:
        raise ValueError(
            f"probe {probe.probe_id}: arm3 ends at {end} but arm5 starts at "
            f"{start} — arms are out of order in the assembled sequence, so the "
            f"probe cannot circularise (see verify_iLock_probe's swapped-stitch "
            f"failure mode)"
        )
    return seq[start:end]


def junction_block(arm3_effective: str, arm5: str, n_each_side: int) -> str:
    """The ``2(n+1)``-nt block straddling the nick in the rotated ligator.

    ``arm3_effective``'s last ``n+1`` bases followed by ``arm5``'s first
    ``n+1`` — exactly the stretch the ligase must see contiguously paired.
    Returns ``""`` when either arm is too short to supply its half.
    """
    if n_each_side < 0:
        raise ValueError(f"n_each_side must be >= 0, got {n_each_side}")
    half = n_each_side + 1
    if len(arm3_effective) < half or len(arm5) < half:
        return ""
    return arm3_effective[-half:] + arm5[:half]


def _build_geometry(probe: ProbeForScreen, n_each_side: int) -> _Geometry:
    arm3_eff, arm5, splint = build_v2_geom(probe)
    block = junction_block(arm3_eff, arm5, n_each_side)
    arm3_literal = probe.probe_arm3.upper()
    arm5_at = splint.index(arm5)
    arm3_at = splint.rindex(arm3_literal)
    return _Geometry(
        probe_id=probe.probe_id,
        chemistry=probe.chemistry,
        arm3_effective=arm3_eff,
        arm5=arm5,
        ligator=arm3_eff + arm5,
        splint_seq=splint,
        junction_block=block,
        needle=_rc(block) if block else "",
        arm_spans=(
            (arm5_at, arm5_at + len(arm5)),
            (arm3_at, arm3_at + len(arm3_literal)),
        ),
        target=probe.target,
    )


def _find_all(haystack: str, needle: str) -> Iterator[int]:
    """Every start index of ``needle`` in ``haystack``, overlaps included."""
    if not needle:
        return
    start = haystack.find(needle)
    while start != -1:
        yield start
        start = haystack.find(needle, start + 1)


def find_ligation_registers(
    geom: _Geometry, splint_seq: str, n_each_side: int,
) -> List[int]:
    """Nick positions ``j`` on ``splint_seq`` where the ligase clamp is satisfied.

    The clamp block spans ``splint[j-n .. j+1+n]``, so an occurrence of
    ``needle`` at index ``p`` means ``j = p + n``. Exact and exhaustive: any
    register a ligase could act on shows up here.
    """
    return [p + n_each_side for p in _find_all(splint_seq, geom.needle)]


# ----------------------------------------------------------------------
# Register scoring
# ----------------------------------------------------------------------


def _is_wc(a: str, b: str) -> bool:
    return _COMPLEMENT.get(a) == b


def _register_window(
    ligator: str, arm3_len: int, arm5_len: int, splint: str, j: int,
) -> Tuple[str, str, int]:
    """Clip the register to the splint; return ``(lig_sub, c_sub, nick_idx)``.

    ``lig[i]`` pairs ``splint[hi - i]`` where ``hi = j + arm3_len``. Bases whose
    partner falls off either end of the splint are dropped — they simply dangle.
    ``c_sub[i]`` is the partner of ``lig_sub[i]`` (Biopython's ``c_seq``
    convention), and ``nick_idx`` is the index in ``lig_sub`` of arm5's 5'-P, so
    the nick sits immediately before it.
    """
    n_splint = len(splint)
    hi = j + arm3_len                 # splint partner of ligator[0]
    lo = j - arm5_len + 1             # splint partner of ligator[-1]
    clip_hi = min(hi, n_splint - 1)
    clip_lo = max(lo, 0)
    i_start = hi - clip_hi
    i_end = hi - clip_lo              # inclusive
    lig_sub = ligator[i_start:i_end + 1]
    window = splint[clip_lo:clip_hi + 1]
    c_sub = window[::-1]              # c_sub[i] pairs lig_sub[i]
    return lig_sub, c_sub, arm3_len - i_start


def _fray_trim(
    lig_sub: str, c_sub: str, nick_idx: int,
    max_consecutive_mismatch: int = MAX_CONSECUTIVE_MISMATCH,
) -> Tuple[int, int]:
    """Half-open bounds of the helix around the nick, after fraying.

    Walks outward from the nick in both directions and stops at
    ``max_consecutive_mismatch`` mismatches in a row; the returned span always
    begins and ends on a Watson-Crick pair. Everything past the fray point is a
    dangling tail that neither holds the complex together nor has NN parameters
    — scoring it would inflate the Tm on exactly the mismatch-rich registers
    the clamp gate lets through by luck.
    """
    left = nick_idx
    run = 0
    for i in range(nick_idx - 1, -1, -1):
        if _is_wc(lig_sub[i], c_sub[i]):
            run = 0
            left = i
        else:
            run += 1
            if run >= max_consecutive_mismatch:
                break

    right = nick_idx
    run = 0
    for i in range(nick_idx, len(lig_sub)):
        if _is_wc(lig_sub[i], c_sub[i]):
            run = 0
            right = i + 1
        else:
            run += 1
            if run >= max_consecutive_mismatch:
                break
    return left, right


def _junction_run(lig_sub: str, c_sub: str, nick_idx: int) -> int:
    """Contiguous Watson-Crick run straddling the nick, in nt."""
    run = 0
    i = nick_idx - 1
    while i >= 0 and _is_wc(lig_sub[i], c_sub[i]):
        run += 1
        i -= 1
    i = nick_idx
    while i < len(lig_sub) and _is_wc(lig_sub[i], c_sub[i]):
        run += 1
        i += 1
    return run


def render_register(lig_sub: str, c_sub: str, nick_idx: int) -> str:
    """Four-line text alignment of one register; ``^`` marks the nick."""
    bars = "".join("|" if _is_wc(a, b) else "." for a, b in zip(lig_sub, c_sub))
    return (
        f"lig 5'-{lig_sub}-3'\n"
        f"       {bars}\n"
        f"spl 3'-{c_sub}-5'\n"
        f"       {' ' * nick_idx}^"
    )


def score_register(
    geom: _Geometry, splint: str, j: int, *,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
) -> Tuple[float, float, float, int, int, str]:
    """Thermodynamics of one fixed register.

    Returns ``(register_tm, arm3_tm, arm5_tm, junction_run, paired_nt,
    alignment)``, all measured over the frayed helix core (see
    :func:`_fray_trim`) rather than the full arm span. The register Tm treats
    the two arms as one duplex across the nick (coaxial stacking); the per-arm
    Tms split that same alignment at the nick, so a register carried mostly by
    one arm stays visible in the report.
    """
    arm3_len = len(geom.arm3_effective)
    lig_full, c_full, nick_full = _register_window(
        geom.ligator, arm3_len, len(geom.arm5), splint, j,
    )
    left, right = _fray_trim(lig_full, c_full, nick_full)
    lig_sub, c_sub = lig_full[left:right], c_full[left:right]
    nick_idx = nick_full - left
    tm = duplex_melting_temperature(lig_sub, c_sub, CROSS_LIG_NN_MODEL, reaction)
    arm3_tm = (
        duplex_melting_temperature(
            lig_sub[:nick_idx], c_sub[:nick_idx], CROSS_LIG_NN_MODEL, reaction,
        )
        if nick_idx >= 2 else 0.0
    )
    arm5_tm = (
        duplex_melting_temperature(
            lig_sub[nick_idx:], c_sub[nick_idx:], CROSS_LIG_NN_MODEL, reaction,
        )
        if len(lig_sub) - nick_idx >= 2 else 0.0
    )
    paired = sum(1 for a, b in zip(lig_sub, c_sub) if _is_wc(a, b))
    return (
        tm, arm3_tm, arm5_tm,
        _junction_run(lig_sub, c_sub, nick_idx), paired,
        render_register(lig_sub, c_sub, nick_idx),
    )


def _overlaps_arms(geom: _Geometry, j: int, splint: str) -> bool:
    """True if this register would ask arm bases to pair with themselves.

    Only meaningful intramolecularly (self-circularisation), where the splint
    IS the ligator's own molecule.
    """
    lo = max(j - len(geom.arm5) + 1, 0)
    hi = min(j + len(geom.arm3_effective), len(splint) - 1)
    for span_lo, span_hi in geom.arm_spans:
        if lo <= span_hi - 1 and span_lo <= hi:
            return True
    return False


def _best_register_dimer(
    lig_geom: _Geometry, splint_geom: _Geometry, *,
    reaction: ReactionConditions,
    tm_threshold_c: float,
    n_each_side: int,
    exclude_arm_overlap: bool = False,
) -> Optional[LigationDimer]:
    """Highest-Tm clamp-satisfying register for one direction, or None."""
    splint = splint_geom.splint_seq
    registers = find_ligation_registers(lig_geom, splint, n_each_side)
    if exclude_arm_overlap:
        registers = [j for j in registers if not _overlaps_arms(lig_geom, j, splint)]
    if not registers:
        return None

    best: Optional[Tuple[float, int, float, float, int, int, str]] = None
    for j in registers:
        tm, arm3_tm, arm5_tm, run, paired, ascii_ = score_register(
            lig_geom, splint, j, reaction=reaction,
        )
        if best is None or tm > best[0]:
            best = (tm, j, arm3_tm, arm5_tm, run, paired, ascii_)

    tm, j, arm3_tm, arm5_tm, run, paired, ascii_ = best
    return LigationDimer(
        seq_a_id=lig_geom.probe_id,
        seq_b_id=splint_geom.probe_id,
        a_target=lig_geom.target,
        b_target=splint_geom.target,
        nick_pos_on_b=j,
        b_3oh_pos=j + 1,
        b_5p_pos=j,
        junction_run_nt=run,
        paired_nt=paired,
        overall_tm_c=tm,
        arm3_tm_c=arm3_tm,
        arm5_tm_c=arm5_tm,
        limiting_arm_tm_c=min(arm3_tm, arm5_tm),
        flagged_overall=tm > tm_threshold_c,
        a_can_ligate_on_b=True,
        vicinity_contiguous=True,
        vicinity_n_each_side=n_each_side,
        is_self_pair=lig_geom.probe_id == splint_geom.probe_id,
        alignment=ascii_,
        n_registers=len(registers),
    )


# ----------------------------------------------------------------------
# Public entry points
# ----------------------------------------------------------------------


def screen_cross_ligation_v2(
    probes: List[ProbeForScreen], *,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
    include_self_pairs: bool = True,
) -> Tuple[List[JunctionSeedHitV2], List[LigationDimer]]:
    """v3 register-enumerating cross-ligation screen.

    Every pair is screened in both directions, plus — new in v3 — each probe
    against another copy of itself, which ``combinations`` had excluded even
    though it is the highest-concentration encounter in any pool.

    A direction is a confirmed risk iff ``flagged_overall and
    a_can_ligate_on_b``. ``vicinity_contiguous`` is True on every returned
    dimer: in v3 the clamp is the Tier 1 gate, so a register that failed it
    never becomes a dimer at all.

    Args:
        tm_threshold_c: flag a register whose duplex Tm clears this.
        reaction: hybridisation conditions; defaults to the lab protocol.
        vicinity_n_each_side: contiguous WC pairs the ligase needs on each side
            of the nick. Must be >= 1 — 0 would admit every position on every
            splint and is rejected rather than silently honoured.
        include_self_pairs: screen A against another copy of A.

    Returns ``(tier1_hits, dimers)``. ``tier1_hits`` carries every
    clamp-satisfying register found, so the funnel stays inspectable.

    Raises ``ValueError`` on registry inconsistency or a disabled clamp.
    """
    if vicinity_n_each_side < 1:
        raise ValueError(
            "vicinity_n_each_side must be >= 1; the contiguous clamp around the "
            "nick is what makes a register ligation-competent, and disabling it "
            "would flag every position on every splint"
        )
    if not probes:
        return [], []

    geoms: Dict[str, _Geometry] = {
        p.probe_id: _build_geometry(p, vicinity_n_each_side) for p in probes
    }

    tier1: List[JunctionSeedHitV2] = []
    dimers: List[LigationDimer] = []

    def _one_direction(lig: _Geometry, spl: _Geometry) -> None:
        registers = find_ligation_registers(lig, spl.splint_seq, vicinity_n_each_side)
        if not registers:
            return
        tier1.append(JunctionSeedHitV2(
            probe_a_id=lig.probe_id, probe_b_id=spl.probe_id,
            junction_block=lig.junction_block,
            nick_positions=tuple(registers),
        ))
        dimer = _best_register_dimer(
            lig, spl, reaction=reaction, tm_threshold_c=tm_threshold_c,
            n_each_side=vicinity_n_each_side,
        )
        if dimer is not None:
            dimers.append(dimer)

    for a, b in combinations(probes, 2):
        _one_direction(geoms[a.probe_id], geoms[b.probe_id])
        _one_direction(geoms[b.probe_id], geoms[a.probe_id])

    if include_self_pairs:
        for probe in probes:
            _one_direction(geoms[probe.probe_id], geoms[probe.probe_id])

    return tier1, dimers


def screen_self_circularisation(
    probes: List[ProbeForScreen], *,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
) -> List[SelfCircularisationHit]:
    """Probes that can circularise on their own fold, with no splint at all.

    Same register formalism with the probe's own sequence as splint, minus the
    registers that would ask an arm base to pair with itself. A hit means one
    molecule can present its 3'-OH and 5'-P to a ligase unaided — the purest
    form of target-independent background, and unscreened before v3.

    The register test is necessary but not sufficient: it does not check that
    the intervening loop can physically close. Treat hits as candidates for a
    fold check (``ext.nupack``), not as proven circularisers.
    """
    if vicinity_n_each_side < 1:
        raise ValueError("vicinity_n_each_side must be >= 1")

    hits: List[SelfCircularisationHit] = []
    for probe in probes:
        geom = _build_geometry(probe, vicinity_n_each_side)
        dimer = _best_register_dimer(
            geom, geom, reaction=reaction, tm_threshold_c=tm_threshold_c,
            n_each_side=vicinity_n_each_side, exclude_arm_overlap=True,
        )
        if dimer is None:
            continue
        hits.append(SelfCircularisationHit(
            probe_id=probe.probe_id, target=probe.target,
            nick_pos=dimer.nick_pos_on_b,
            junction_run_nt=dimer.junction_run_nt,
            paired_nt=dimer.paired_nt,
            tm_c=dimer.overall_tm_c,
            flagged=dimer.flagged_overall,
            alignment=dimer.alignment,
        ))
    return hits


# ----------------------------------------------------------------------
# Report writers
# ----------------------------------------------------------------------


REPORT_COLUMNS: Tuple[str, ...] = (
    "seq_a_id", "seq_b_id", "a_target", "b_target",
    "overall_tm_c", "limiting_arm_tm_c", "arm3_tm_c", "arm5_tm_c",
    "junction_run_nt", "paired_nt",
    "flagged_overall", "a_can_ligate_on_b", "is_self_pair",
    "nick_pos_on_b", "b_3oh_pos", "b_5p_pos",
    "vicinity_n_each_side", "n_registers", "alignment",
)

SELF_CIRC_COLUMNS: Tuple[str, ...] = (
    "probe_id", "target", "tm_c", "junction_run_nt", "paired_nt",
    "nick_pos", "flagged", "alignment",
)


def _esc(text: str) -> str:
    return (text or "").replace("\n", "\\n").replace("\t", "\\t")


def write_dimer_report(dimers: Sequence[LigationDimer], tsv_path: Path | str) -> None:
    """Write :class:`LigationDimer` rows to TSV, highest Tm first.

    An empty list still writes the header. Newlines inside ``alignment`` are
    escaped so each row stays on one line.
    """
    tsv_path = Path(tsv_path)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(list(REPORT_COLUMNS))
        for d in sorted(dimers, key=lambda x: x.overall_tm_c, reverse=True):
            writer.writerow([
                d.seq_a_id, d.seq_b_id, d.a_target, d.b_target,
                f"{d.overall_tm_c:.3f}", f"{d.limiting_arm_tm_c:.3f}",
                f"{d.arm3_tm_c:.3f}", f"{d.arm5_tm_c:.3f}",
                d.junction_run_nt, d.paired_nt,
                d.flagged_overall, d.a_can_ligate_on_b, d.is_self_pair,
                d.nick_pos_on_b,
                d.b_3oh_pos if d.b_3oh_pos is not None else "",
                d.b_5p_pos if d.b_5p_pos is not None else "",
                d.vicinity_n_each_side, d.n_registers, _esc(d.alignment),
            ])


def write_self_circ_report(
    hits: Sequence[SelfCircularisationHit], tsv_path: Path | str,
) -> None:
    """Write :class:`SelfCircularisationHit` rows to TSV, highest Tm first."""
    tsv_path = Path(tsv_path)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(list(SELF_CIRC_COLUMNS))
        for h in sorted(hits, key=lambda x: x.tm_c, reverse=True):
            writer.writerow([
                h.probe_id, h.target, f"{h.tm_c:.3f}",
                h.junction_run_nt, h.paired_nt, h.nick_pos,
                h.flagged, _esc(h.alignment),
            ])


# ----------------------------------------------------------------------
# primer3-ASCII helpers — diagnostics only, NOT the screening criterion
# ----------------------------------------------------------------------
#
# ``qc.viz.plot_cross_lig_dimer`` renders primer3's chosen alignment for a pair.
# That view answers "where does each arm bind best", which is worth looking at
# and is NOT what decides a cross-lig call — see the module docstring. These two
# helpers back that plot and are kept for it alone.


def _check_nick_geometry(
    arm3_ascii: str, arm5_ascii: str, arm3_eff_len: int,
) -> Tuple[bool, Optional[int], Optional[int]]:
    """Whether primer3's two argmax alignments happen to be nick-adjacent.

    Diagnostic only. v2.2 used this as the screening criterion; measured, it
    fires at chance rate (see module docstring).
    """
    arm3_pairs = parse_primer3_dimer_pairing(arm3_ascii)["pair_positions_a_to_b"]
    arm5_pairs = parse_primer3_dimer_pairing(arm5_ascii)["pair_positions_a_to_b"]
    b_3oh = arm3_pairs.get(arm3_eff_len - 1)
    b_5p = arm5_pairs.get(0)
    if b_3oh is None or b_5p is None:
        return False, b_3oh, b_5p
    return (b_3oh - b_5p == 1), b_3oh, b_5p


def _check_per_arm_vicinity(
    arm3_ascii: str, arm5_ascii: str,
    arm3_eff_len: int, arm5_len: int,
    n_each_side: int,
) -> bool:
    """Per-arm contiguity inside primer3's argmax alignments. Diagnostic only."""
    if n_each_side <= 0:
        return False
    arm3_pairs = parse_primer3_dimer_pairing(arm3_ascii)["pair_positions_a_to_b"]
    arm5_pairs = parse_primer3_dimer_pairing(arm5_ascii)["pair_positions_a_to_b"]

    arm3_positions = list(range(arm3_eff_len - 1 - n_each_side, arm3_eff_len))
    if any(a not in arm3_pairs for a in arm3_positions):
        return False
    arm3_b = [arm3_pairs[a] for a in arm3_positions]
    if any(arm3_b[i + 1] - arm3_b[i] != -1 for i in range(len(arm3_b) - 1)):
        return False

    if arm5_len <= n_each_side:
        return False
    arm5_positions = list(range(0, n_each_side + 1))
    if any(a not in arm5_pairs for a in arm5_positions):
        return False
    arm5_b = [arm5_pairs[a] for a in arm5_positions]
    if any(arm5_b[i + 1] - arm5_b[i] != -1 for i in range(len(arm5_b) - 1)):
        return False
    return True


__all__ = [
    "ProbeForScreen", "JunctionSeedHitV2", "LigationDimer",
    "SelfCircularisationHit",
    "screen_cross_ligation_v2", "screen_self_circularisation",
    "build_v2_geom", "ligation_substrate", "junction_block",
    "find_ligation_registers", "score_register", "render_register",
    "write_dimer_report", "write_self_circ_report",
    "REPORT_COLUMNS", "SELF_CIRC_COLUMNS",
    "LAB_HYBRIDISATION", "CROSS_LIG_NN_MODEL",
    "DEFAULT_TM_THRESHOLD_C", "DEFAULT_VICINITY_N",
    "DEFAULT_MV_CONC_MM", "DEFAULT_DV_CONC_MM",
    "DEFAULT_DNA_CONC_M", "DEFAULT_TEMP_C",
]
