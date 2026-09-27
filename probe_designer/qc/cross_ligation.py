"""Cross-ligation screen — enumerate ligation registers, don't take an argmax.

**What the screen is for.** Probe A's two arms can anneal to some other template
instead of to its mRNA target. If they land nick-adjacent on it, ligase closes
**A into a real circle** carrying A's own barcode, with no target RNA involved.
It amplifies and decodes as gene A everywhere, and being circular it survives the
exonuclease step. 10x describes the same artifact for Xenium panels.

**Why the argmax model was replaced (audit 2026-09-02).** v2.2 called
``primer3.calc_heterodimer`` once per arm against the splint and asked whether
the two independently-chosen best alignments happened to land nick-adjacent.
That is a coincidence test, not a search. Measured on ``TNBCmarker_main_v1``:
arm3's 3'-OH was paired in 16.5% of dimers, arm5's 5'-P in 14.7%, and *both* in
2.4% — the product of the two, i.e. uncorrelated events rather than one physical
complex. Only 3 of 37 187 came out adjacent, about the rate chance predicts. It
missed the pool's strongest ligation-competent pair (CD3E on CD3D) while
reporting the second-ranked one.

The root cause is that the productive complex's stability is a property of the
two arms *together* — one duplex spanning the nick — but v2.2 scored each arm
alone, then demanded each weak half win its own global argmax. A 12 bp half loses
to any other >=12 bp match on the splint, and since every probe in a pool carries
the same backbone, such matches are the norm.

**The model.** The nick position determines the register. With the nick sitting
between splint bases ``j`` and ``j+1``:

* arm3 occupies ``B[j+1 .. j+L3]`` — its 3'-OH pairs ``B[j+1]``
* arm5 occupies ``B[j-L5+1 .. j]`` — its 5'-P pairs ``B[j]``

Both arms run antiparallel to B, so moving 5'->3' along the probe moves 3'->5'
along B; hence ``b_3oh - b_5p == +1``. The two windows are contiguous, so the
rotated ligator ``arm3_effective + arm5_effective`` maps onto one contiguous
splint window.

**Tier 1 is exact, and it is an index lookup.** The ligase clamp requires
``n+1`` contiguous Watson-Crick pairs each side of the nick
(``chemistry.LIGASE_CLAMP_NT``), which in the rotated form is one fixed
``2(n+1)``-nt block straddling the junction. So a ligator has exactly one needle,
and every splint is indexed by its ``2(n+1)``-mers once — no k-mer seeding, no
false negatives, and no N^2 scan: each ligator is a single dict lookup that hands
back every ``(splint, j)`` it could ligate on.

**Tier 2 scores the fixed alignment.** ``nn_tables.duplex_melting_temperature``
computes the NN Tm of the rotated ligator against the splint window *at that
register*, mismatches included, trimmed to the helix core that survives fraying.
Because arm3 and arm5 stack coaxially across a nick, treating the pair as one
duplex is the standard approximation; it slightly overestimates stability, which
is the safe direction for a screen.

Fixing the register also removes v2.1's error without reintroducing it: v2.1 fed
the rotated 40-mer to primer3, which was then free to *slide* the alignment so
the arm3/arm5 seam sat somewhere other than the nick — geometries that cannot
form before ligation.

**A splint is not necessarily a probe.** :class:`Splint` is just a labelled
sequence with optional masked spans, so the same machinery screens probe-vs-probe
(:func:`screen_cross_ligation`), a probe against its own fold
(:func:`screen_self_circularisation`), or a panel against transcript sequences
(:func:`screen_against_splints` with ``model=nn_model_for("DNA", "RNA")``).

**Also screened here, and not screened at all before**: ``A`` templating on
*another copy of A* — ``combinations`` had excluded the highest-concentration
encounter in any pool — and self-circularisation, one molecule folding so its own
3'-OH and 5'-P abut.

**Architecture**: this module is BANK-FREE. Pool loading is the job of
``probe_designer.ext.pool.loader``; the static check at
``designer/tests/unit/test_no_bank_import_in_qc.py`` enforces it.

See ``experiments/20260902_crosslig_register_scan`` and
``[[project-cross-lig-v3-register-scan]]``.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from probe_designer.chemistry import (
    LIGASE_CLAMP_NT,
    ReactionConditions,
    is_watson_crick,
    reverse_complement,
)
from probe_designer.nn_tables import DNA_DNA, NNModel, duplex_melting_temperature


# ----------------------------------------------------------------------
# Screen parameters
# ----------------------------------------------------------------------

#: The hybridisation the screen models. ``ReactionConditions`` defaults ARE the
#: lab protocol (rca.md v5.3, mirrored in ``configs/config_template.yaml``):
#: 75 mM K+, 10 mM Mg2+, 100 nM strand, 20% formamide, 45 C — 55 C effective in
#: the no-formamide NN model. Pass an explicit ``reaction`` for a different
#: buffer rather than editing this; in an n-plex pool the per-species
#: concentration is the pool stock divided by n, which is the number that
#: belongs in ``strand_nM`` for a bimolecular screen.
LAB_HYBRIDISATION = ReactionConditions()

#: Probe-vs-probe is DNA:DNA — the splint is another padlock, not RNA. Screening
#: against transcripts wants ``nn_tables.nn_model_for("DNA", "RNA")`` instead;
#: pass it as ``model=``.
CROSS_LIG_NN_MODEL = DNA_DNA

#: Re-exported from ``chemistry`` so callers of this module need not reach past
#: it; ``chemistry.LIGASE_CLAMP_NT`` remains the single definition.
DEFAULT_VICINITY_N: int = LIGASE_CLAMP_NT

#: Flag a register whose duplex Tm clears this. It is an absolute floor, about
#: 28 C below the 55 C effective reaction temperature — deliberately far below
#: it, because ligation is irreversible: a register populated at a small
#: equilibrium fraction is a kinetic trap that accumulates product over an
#: overnight incubation, not a species you can dismiss by its share. (v2.2
#: carried this same 27 C but documented it as "T_rxn - 10" for a 37 C reaction,
#: which stopped being true when the anchor moved to 45 C.)
DEFAULT_TM_THRESHOLD_C: float = 27.0

#: Consecutive mismatches that end the helix. Walking outward from the nick, a
#: run this long stops the duplex: the helix frays there, and it is also where
#: the NN model runs out — Biopython parameterises *single* internal mismatches
#: (DNA_IMM1) but has no term for two in a row. Trimming to a core that starts
#: and ends on a Watson-Crick pair keeps every scored neighbour inside the tables.
MAX_CONSECUTIVE_MISMATCH: int = 2


# ----------------------------------------------------------------------
# Inputs
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
class Splint:
    """A template a ligator might circularise on.

    Deliberately *not* a probe. A padlock, an mRNA and a synthetic control are
    all just a labelled sequence here, which is what lets the transcriptome scan
    reuse this machinery instead of growing a parallel one.

    ``forbidden_spans`` are half-open ``[start, end)`` regions of ``sequence``
    that may not template a register. Two uses so far: the ligator's own arms
    when a probe folds on itself (a base cannot pair with itself), and, for a
    transcriptome scan, the probe's intended target site — without which every
    probe flags on its own transcript.
    """
    splint_id: str
    sequence: str
    label: str = ""                                   # gene / target, for reports
    forbidden_spans: Tuple[Tuple[int, int], ...] = ()


@dataclass(frozen=True)
class LigatorGeometry:
    """One probe resolved into the parts the register scan needs.

    Public because it is the unit of reuse: build it once per probe and screen it
    against any number of splints.
    """
    probe_id: str
    chemistry: str
    arm3_effective: str          # its LAST base carries the 3'-OH
    arm5_effective: str          # its FIRST base carries the 5'-P
    ligator: str                 # arm3_effective + arm5_effective, rotated form
    sequence: str                # full assembled probe
    junction_block: str          # 2(n+1) nt straddling the nick, in the ligator
    needle: str                  # reverse complement of junction_block
    arm_spans: Tuple[Tuple[int, int], ...]   # half-open arm spans in sequence
    target: str = ""

    def as_splint(self) -> Splint:
        """This probe seen as a template for someone else's ligation."""
        return Splint(self.probe_id, self.sequence, self.target)

    def as_self_splint(self) -> Splint:
        """This probe seen as its own template, with its arms masked."""
        return Splint(self.probe_id, self.sequence, self.target, self.arm_spans)


# ----------------------------------------------------------------------
# Outputs
# ----------------------------------------------------------------------


@dataclass
class LigationDimer:
    """The best ligation-competent register for one (ligator, splint) direction.

    ``overall_tm_c`` is the Tm of the rotated ligator against the splint window
    *at this register* — one duplex spanning the nick, and the quantity that
    decides whether the complex is populated. It is the only number used as a
    gate.

    The per-arm Tms are **diagnostics**: they split that same alignment at the
    nick, so an asymmetric register (say 8 paired bases on the arm3 side and 4 on
    the arm5 side) is visible rather than hidden inside one figure. A half of only
    a few bases melts below 0 C and will report a large negative Tm — that is the
    model saying "this half cannot hold on its own", not an error. A half shorter
    than 2 nt has no nearest neighbour at all and reports 0.0.
    """
    seq_a_id: str                # ligator (A)
    seq_b_id: str                # splint  (B)
    a_target: str = ""
    b_target: str = ""
    # Register
    nick_pos_on_b: int = -1      # nick sits between B[j] and B[j+1]
    nick_positions: Tuple[int, ...] = ()   # every clamp-satisfying j considered
    junction_run_nt: int = 0     # contiguous WC run straddling the nick
    paired_nt: int = 0           # WC pairs over the scored helix
    # Thermodynamics at that register
    overall_tm_c: float = 0.0
    arm3_tm_c: float = 0.0
    arm5_tm_c: float = 0.0
    # Verdict
    flagged_overall: bool = False       # overall_tm_c > threshold
    vicinity_n_each_side: int = DEFAULT_VICINITY_N
    alignment: str = ""                 # rendered register

    @property
    def b_3oh_pos(self) -> int:
        """Splint base paired to the ligator's 3'-OH — one past the nick."""
        return self.nick_pos_on_b + 1

    @property
    def b_5p_pos(self) -> int:
        """Splint base paired to the ligator's 5'-P — the nick position."""
        return self.nick_pos_on_b

    @property
    def is_self_pair(self) -> bool:
        """Ligator and splint are the same probe (a second copy of it)."""
        return self.seq_a_id == self.seq_b_id

    @property
    def limiting_arm_tm_c(self) -> float:
        """The weaker half. Ligation needs BOTH arms annealed, so ``min``.

        v2.2 aggregated with ``max`` *and gated on it*, which let one strong arm
        carry a pair whose other arm was useless.
        """
        return min(self.arm3_tm_c, self.arm5_tm_c)

    @property
    def n_registers(self) -> int:
        return len(self.nick_positions)


@dataclass
class SelfCircularisationHit:
    """One probe folding so its own 3'-OH and 5'-P abut — no splint needed.

    The same register formalism with the probe's own sequence as splint and its
    arms masked. This is target-independent background at its purest: the circle
    forms from a single molecule.

    Necessary but not sufficient — it does not check that the intervening loop
    can physically close. Treat hits as candidates for a fold check
    (``ext.nupack``), not as proven circularisers.
    """
    probe_id: str
    target: str = ""
    nick_pos: int = -1
    junction_run_nt: int = 0
    paired_nt: int = 0
    tm_c: float = 0.0
    flagged: bool = False
    alignment: str = ""


# ----------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------


def arm_offsets(probe: ProbeForScreen) -> Tuple[int, int, int, int]:
    """``(arm5_start, arm5_end, arm3_start, arm3_end)`` in the assembled sequence.

    Half-open. One derivation of "where do the arms sit", so that the ligation
    substrate, the masked spans and the backbone loop length cannot disagree.

    Raises ``ValueError`` on registry inconsistency (arm not in sequence, or arms
    out of order — the swapped-stitch failure mode ``verify_iLock_probe``
    guards, which binds but cannot ligate).
    """
    arm3_eff, arm5_eff, seq = ligation_geometry(probe)
    arm5_start = seq.index(arm5_eff)
    arm3_end = seq.rindex(arm3_eff) + len(arm3_eff)
    if arm3_end <= arm5_start:
        raise ValueError(
            f"probe {probe.probe_id}: arm3 ends at {arm3_end} but arm5 starts at "
            f"{arm5_start} — the arms are out of order in the assembled sequence, "
            f"so the probe binds but cannot circularise"
        )
    return arm5_start, arm5_start + len(arm5_eff), arm3_end - len(arm3_eff), arm3_end


def ligation_geometry(probe: ProbeForScreen) -> Tuple[str, str, str]:
    """``(arm3_effective, arm5_effective, sequence)`` — the ligation substrate.

    ``arm3_effective``'s **last** base carries the 3'-OH; ``arm5_effective``'s
    **first** base carries the 5'-P.

    **iLock drops ``arm5[0]``, not ``arm3[0]``.** The invader designer places the
    same SNP-complement base at both arm tips (``arm5[0] == arm3[-1]``, the
    invariant ``ext.mutation.probe.verify_iLock_probe`` enforces) so they compete
    for the target SNP; FEN1 then cleaves the *flap + arm5[0]* block, leaving
    ``arm3[-1]`` as the 3'-OH and ``arm5[1]`` as the new 5'-P. v2.2 returned
    ``arm3[1:]`` with arm5 whole, which kept the right 3'-OH base by coincidence
    (``arm3[1:][-1] == arm3[-1]``) but put the 5'-P one base early — the junction
    block then carried the SNP base twice and every iLock register was computed
    out of phase.

    Raises ``ValueError`` on registry inconsistency (arm not in sequence).
    """
    arm5 = probe.probe_arm5.upper()
    arm3 = probe.probe_arm3.upper()
    seq_upper = probe.sequence.upper()

    if not seq_upper:
        raise ValueError(
            f"probe {probe.probe_id}: empty assembled sequence — registry inconsistent"
        )
    for name, arm in (("arm5", arm5), ("arm3", arm3)):
        if arm not in seq_upper:
            raise ValueError(
                f"probe {probe.probe_id}: {name} ({arm[:10]}...) not in assembled "
                f"sequence — registry inconsistent"
            )

    if probe.chemistry != "iLock":
        return arm3, arm5, seq_upper
    if len(arm5) < 2:
        raise ValueError(
            f"probe {probe.probe_id}: iLock arm5 must be >= 2 nt — FEN1 removes "
            f"arm5[0] with the flap, leaving nothing to ligate"
        )
    return arm3, arm5[1:], seq_upper


#: Historical name; ``ext.nupack`` and ``qc.viz`` imported it as ``build_v2_geom``
#: when the geometry was v2's. Kept so those imports keep working.
build_v2_geom = ligation_geometry


def junction_block(arm3_effective: str, arm5_effective: str, n_each_side: int) -> str:
    """The ``2(n+1)``-nt block straddling the nick in the rotated ligator.

    ``arm3_effective``'s last ``n+1`` bases followed by ``arm5_effective``'s
    first ``n+1`` — exactly the stretch the ligase must see contiguously paired.
    Returns ``""`` when either arm is too short to supply its half.
    """
    if n_each_side < 0:
        raise ValueError(f"n_each_side must be >= 0, got {n_each_side}")
    half = n_each_side + 1
    if len(arm3_effective) < half or len(arm5_effective) < half:
        return ""
    return arm3_effective[-half:] + arm5_effective[:half]


def build_geometry(
    probe: ProbeForScreen, n_each_side: int = DEFAULT_VICINITY_N,
) -> LigatorGeometry:
    """Resolve one probe into a :class:`LigatorGeometry`."""
    arm3_eff, arm5_eff, seq = ligation_geometry(probe)
    arm5_start, arm5_end, arm3_start, arm3_end = arm_offsets(probe)
    block = junction_block(arm3_eff, arm5_eff, n_each_side)
    return LigatorGeometry(
        probe_id=probe.probe_id,
        chemistry=probe.chemistry,
        arm3_effective=arm3_eff,
        arm5_effective=arm5_eff,
        ligator=arm3_eff + arm5_eff,
        sequence=seq,
        junction_block=block,
        needle=reverse_complement(block) if block else "",
        arm_spans=((arm5_start, arm5_end), (arm3_start, arm3_end)),
        target=probe.target,
    )


def ligation_substrate(probe: ProbeForScreen) -> str:
    """The probe as one strand presenting its 5'-P at index 0 and 3'-OH at -1.

    The assembled sequence trimmed to start at the effective arm5 and end at
    arm3's last base — for dRNA/cDNA the whole probe, for iLock without the flap
    and ``arm5[0]``. Unlike the split arms this keeps the **backbone tether**
    between the two ends.
    """
    arm5_start, _, _, arm3_end = arm_offsets(probe)
    return probe.sequence.upper()[arm5_start:arm3_end]


def backbone_loop_nt(probe: ProbeForScreen) -> int:
    """Backbone nucleotides between the two arms — the tether's contour."""
    _, arm5_end, arm3_start, _ = arm_offsets(probe)
    return arm3_start - arm5_end


# ----------------------------------------------------------------------
# Tier 1 — exact register lookup
# ----------------------------------------------------------------------


def build_splint_index(
    splints: Sequence[Splint], n_each_side: int = DEFAULT_VICINITY_N,
) -> Dict[str, List[Tuple[int, int]]]:
    """Index every splint by its ``2(n+1)``-mers: ``{kmer: [(splint_idx, pos)]}``.

    This is what keeps the screen out of an N^2 substring scan. A ligator has
    exactly one needle (its junction block's reverse complement), so finding
    every register it could ligate on, across every splint, is one dict lookup.
    Keys are bounded by ``4 ** (2 * (n + 1))`` = 65,536 for the default clamp,
    and building costs one pass over the splints.
    """
    width = 2 * (n_each_side + 1)
    index: Dict[str, List[Tuple[int, int]]] = {}
    for splint_idx, splint in enumerate(splints):
        seq = splint.sequence
        for pos in range(len(seq) - width + 1):
            index.setdefault(seq[pos:pos + width], []).append((splint_idx, pos))
    return index


def find_ligation_registers(
    needle: str, splint_seq: str, n_each_side: int = DEFAULT_VICINITY_N,
) -> List[int]:
    """Nick positions ``j`` on one splint where the ligase clamp is satisfied.

    The clamp block spans ``splint[j-n .. j+1+n]``, so an occurrence of ``needle``
    at index ``p`` means ``j = p + n``. Exact and exhaustive.

    The screen itself goes through :func:`build_splint_index`; this is the direct
    form, for a single splint or a caller that has no index.
    """
    if not needle:
        return []
    out: List[int] = []
    start = splint_seq.find(needle)
    while start != -1:
        out.append(start + n_each_side)
        start = splint_seq.find(needle, start + 1)
    return out


def _register_is_allowed(
    geom: LigatorGeometry, splint: Splint, nick: int,
) -> bool:
    """False when the register would use a masked span of the splint."""
    if not splint.forbidden_spans:
        return True
    lo = max(nick - len(geom.arm5_effective) + 1, 0)
    hi = min(nick + len(geom.arm3_effective), len(splint.sequence) - 1)
    return not any(
        lo <= span_hi - 1 and span_lo <= hi
        for span_lo, span_hi in splint.forbidden_spans
    )


# ----------------------------------------------------------------------
# Tier 2 — score the fixed alignment
# ----------------------------------------------------------------------


def _register_window(
    geom: LigatorGeometry, splint_seq: str, nick: int,
) -> Tuple[str, str, int]:
    """Clip the register to the splint; return ``(lig_sub, c_sub, nick_idx)``.

    ``ligator[i]`` pairs ``splint[hi - i]`` where ``hi = nick + len(arm3)``. Bases
    whose partner falls off either end of the splint are dropped — they dangle.
    ``c_sub[i]`` is the partner of ``lig_sub[i]`` (Biopython's ``c_seq``
    convention), and ``nick_idx`` is the index in ``lig_sub`` of arm5's 5'-P, so
    the nick sits immediately before it.
    """
    arm3_len = len(geom.arm3_effective)
    hi = nick + arm3_len                                # partner of ligator[0]
    lo = nick - len(geom.arm5_effective) + 1            # partner of ligator[-1]
    clip_hi = min(hi, len(splint_seq) - 1)
    clip_lo = max(lo, 0)
    i_start = hi - clip_hi
    lig_sub = geom.ligator[i_start:hi - clip_lo + 1]
    c_sub = splint_seq[clip_lo:clip_hi + 1][::-1]       # c_sub[i] pairs lig_sub[i]
    return lig_sub, c_sub, arm3_len - i_start


def _fray_trim(
    lig_sub: str, c_sub: str, nick_idx: int,
    max_consecutive_mismatch: int = MAX_CONSECUTIVE_MISMATCH,
) -> Tuple[int, int]:
    """Half-open bounds of the helix around the nick, after fraying.

    Walks outward from the nick in both directions and stops at
    ``max_consecutive_mismatch`` mismatches in a row; the returned span always
    begins and ends on a Watson-Crick pair. Everything past the fray point is a
    dangling tail that neither holds the complex together nor has NN parameters —
    scoring it would inflate the Tm on exactly the mismatch-rich registers the
    clamp gate lets through by luck.
    """
    left = right = nick_idx
    run = 0
    for i in range(nick_idx - 1, -1, -1):
        if is_watson_crick(lig_sub[i], c_sub[i]):
            run, left = 0, i
        else:
            run += 1
            if run >= max_consecutive_mismatch:
                break

    run = 0
    for i in range(nick_idx, len(lig_sub)):
        if is_watson_crick(lig_sub[i], c_sub[i]):
            run, right = 0, i + 1
        else:
            run += 1
            if run >= max_consecutive_mismatch:
                break
    return left, right


def junction_run(lig_sub: str, c_sub: str, nick_idx: int) -> int:
    """Contiguous Watson-Crick run straddling the nick, in nt.

    The same quantity a BLAST midline would give across the junction, which is
    why it is public: an alignment from anywhere can be judged by this rule.
    """
    total = 0
    for i in range(nick_idx - 1, -1, -1):
        if not is_watson_crick(lig_sub[i], c_sub[i]):
            break
        total += 1
    for i in range(nick_idx, len(lig_sub)):
        if not is_watson_crick(lig_sub[i], c_sub[i]):
            break
        total += 1
    return total


def render_register(lig_sub: str, c_sub: str, nick_idx: int) -> str:
    """Four-line text alignment of one register; ``^`` marks the nick."""
    bars = "".join(
        "|" if is_watson_crick(a, b) else "." for a, b in zip(lig_sub, c_sub)
    )
    return (
        f"lig 5'-{lig_sub}-3'\n"
        f"       {bars}\n"
        f"spl 3'-{c_sub}-5'\n"
        f"       {' ' * nick_idx}^"
    )


@dataclass(frozen=True)
class RegisterScore:
    """Thermodynamics and geometry of one fixed register."""
    nick: int
    tm_c: float
    arm3_tm_c: float
    arm5_tm_c: float
    junction_run_nt: int
    paired_nt: int
    alignment: str


def score_register(
    geom: LigatorGeometry, splint_seq: str, nick: int, *,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
    model: NNModel = CROSS_LIG_NN_MODEL,
) -> RegisterScore:
    """Score one fixed register, measured over the frayed helix core.

    The register Tm treats the two arms as one duplex across the nick (coaxial
    stacking); the per-arm Tms split that same alignment at the nick, so a
    register carried mostly by one arm stays visible in the report.
    """
    lig_full, c_full, nick_full = _register_window(geom, splint_seq, nick)
    left, right = _fray_trim(lig_full, c_full, nick_full)
    lig_sub, c_sub = lig_full[left:right], c_full[left:right]
    nick_idx = nick_full - left

    def tm(seq: str, partner: str) -> float:
        return (
            duplex_melting_temperature(seq, partner, model, reaction)
            if len(seq) >= 2 else 0.0
        )

    return RegisterScore(
        nick=nick,
        tm_c=tm(lig_sub, c_sub),
        arm3_tm_c=tm(lig_sub[:nick_idx], c_sub[:nick_idx]),
        arm5_tm_c=tm(lig_sub[nick_idx:], c_sub[nick_idx:]),
        junction_run_nt=junction_run(lig_sub, c_sub, nick_idx),
        paired_nt=sum(1 for a, b in zip(lig_sub, c_sub) if is_watson_crick(a, b)),
        alignment=render_register(lig_sub, c_sub, nick_idx),
    )


def best_register(
    geom: LigatorGeometry, splint: Splint, nicks: Sequence[int], *,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
    model: NNModel = CROSS_LIG_NN_MODEL,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    n_each_side: int = DEFAULT_VICINITY_N,
) -> Optional[LigationDimer]:
    """Highest-Tm allowed register among ``nicks``, or None if there is none."""
    allowed = [j for j in nicks if _register_is_allowed(geom, splint, j)]
    if not allowed:
        return None
    scores = [
        score_register(geom, splint.sequence, j, reaction=reaction, model=model)
        for j in allowed
    ]
    best = max(scores, key=lambda s: s.tm_c)
    return LigationDimer(
        seq_a_id=geom.probe_id,
        seq_b_id=splint.splint_id,
        a_target=geom.target,
        b_target=splint.label,
        nick_pos_on_b=best.nick,
        nick_positions=tuple(allowed),
        junction_run_nt=best.junction_run_nt,
        paired_nt=best.paired_nt,
        overall_tm_c=best.tm_c,
        arm3_tm_c=best.arm3_tm_c,
        arm5_tm_c=best.arm5_tm_c,
        flagged_overall=best.tm_c > tm_threshold_c,
        vicinity_n_each_side=n_each_side,
        alignment=best.alignment,
    )


# ----------------------------------------------------------------------
# Entry points
# ----------------------------------------------------------------------


def _check_clamp(n_each_side: int) -> None:
    if n_each_side < 1:
        raise ValueError(
            "vicinity_n_each_side must be >= 1; the contiguous clamp around the "
            "nick is what makes a register ligation-competent, and disabling it "
            "would flag every position on every splint"
        )


def screen_against_splints(
    ligators: Sequence[ProbeForScreen],
    splints: Sequence[Splint], *,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
    model: NNModel = CROSS_LIG_NN_MODEL,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
) -> List[LigationDimer]:
    """Every ligation-competent (ligator, splint) register, best one per pair.

    The general form: ``splints`` need not be probes. Screening a panel against
    transcript sequences is this function with ``Splint(transcript_id, seq, gene)``
    and ``model=nn_model_for("DNA", "RNA")``.

    Raises ``ValueError`` on registry inconsistency or a disabled clamp.
    """
    _check_clamp(vicinity_n_each_side)
    if not ligators or not splints:
        return []

    index = build_splint_index(splints, vicinity_n_each_side)
    dimers: List[LigationDimer] = []
    for probe in ligators:
        geom = build_geometry(probe, vicinity_n_each_side)
        by_splint: Dict[int, List[int]] = {}
        for splint_idx, pos in index.get(geom.needle, ()):
            by_splint.setdefault(splint_idx, []).append(pos + vicinity_n_each_side)
        for splint_idx, nicks in by_splint.items():
            dimer = best_register(
                geom, splints[splint_idx], nicks,
                reaction=reaction, model=model,
                tm_threshold_c=tm_threshold_c, n_each_side=vicinity_n_each_side,
            )
            if dimer is not None:
                dimers.append(dimer)
    return dimers


def screen_cross_ligation(
    probes: Sequence[ProbeForScreen], *,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
    include_self_pairs: bool = True,
) -> List[LigationDimer]:
    """Screen a pool against itself, both directions, self-pairs included.

    A probe against another copy of itself is the highest-concentration encounter
    in any pool and ``combinations`` had excluded it; ``include_self_pairs=False``
    restores the old behaviour for comparison.

    A direction is a confirmed risk iff ``flagged_overall`` — the clamp is the
    Tier 1 gate, so a register that failed it never becomes a dimer at all.
    """
    _check_clamp(vicinity_n_each_side)
    if not probes:
        return []

    geoms = [build_geometry(p, vicinity_n_each_side) for p in probes]
    splints = [g.as_splint() for g in geoms]
    dimers = screen_against_splints(
        probes, splints,
        tm_threshold_c=tm_threshold_c, reaction=reaction,
        vicinity_n_each_side=vicinity_n_each_side,
    )
    if include_self_pairs:
        return dimers
    return [d for d in dimers if not d.is_self_pair]


def screen_self_circularisation(
    probes: Sequence[ProbeForScreen], *,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
) -> List[SelfCircularisationHit]:
    """Probes that can circularise on their own fold, with no splint at all.

    Same register formalism with the probe's own sequence as splint and its arms
    masked (a base cannot pair with itself). A hit means one molecule can present
    its 3'-OH and 5'-P to a ligase unaided — the purest form of
    target-independent background.
    """
    _check_clamp(vicinity_n_each_side)
    hits: List[SelfCircularisationHit] = []
    for probe in probes:
        geom = build_geometry(probe, vicinity_n_each_side)
        splint = geom.as_self_splint()
        nicks = find_ligation_registers(
            geom.needle, splint.sequence, vicinity_n_each_side,
        )
        dimer = best_register(
            geom, splint, nicks,
            reaction=reaction, tm_threshold_c=tm_threshold_c,
            n_each_side=vicinity_n_each_side,
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


def screen_cross_ligation_v2(
    probes: Sequence[ProbeForScreen], **kwargs,
) -> Tuple[List[LigationDimer], List[LigationDimer]]:
    """Deprecated shim for the old ``(tier1, dimers)`` two-list return.

    The two lists were element-wise identical — ``best_register`` returned a
    record for exactly the directions Tier 1 hit — so the first was a
    same-length restatement of the second. Returns the same list twice; use
    :func:`screen_cross_ligation`.
    """
    dimers = screen_cross_ligation(probes, **kwargs)
    return dimers, dimers


# ----------------------------------------------------------------------
# Report writers
# ----------------------------------------------------------------------

#: ``(column, getter)`` — one description of the screen's output, so the pool CLI
#: and the assemble hook cannot drift apart. ``qc.assemble_hook`` renames the
#: identity columns at the DataFrame boundary for the xlsx and reuses the rest.
DIMER_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("seq_a_id", "seq_a_id"), ("seq_b_id", "seq_b_id"),
    ("a_target", "a_target"), ("b_target", "b_target"),
    ("overall_tm_c", "overall_tm_c"), ("limiting_arm_tm_c", "limiting_arm_tm_c"),
    ("arm3_tm_c", "arm3_tm_c"), ("arm5_tm_c", "arm5_tm_c"),
    ("junction_run_nt", "junction_run_nt"), ("paired_nt", "paired_nt"),
    ("flagged_overall", "flagged_overall"), ("is_self_pair", "is_self_pair"),
    ("nick_pos_on_b", "nick_pos_on_b"),
    ("b_3oh_pos", "b_3oh_pos"), ("b_5p_pos", "b_5p_pos"),
    ("vicinity_n_each_side", "vicinity_n_each_side"),
    ("n_registers", "n_registers"), ("alignment", "alignment"),
)

REPORT_COLUMNS: Tuple[str, ...] = tuple(name for name, _ in DIMER_COLUMNS)

SELF_CIRC_COLUMNS: Tuple[str, ...] = (
    "probe_id", "target", "tm_c", "junction_run_nt", "paired_nt",
    "nick_pos", "flagged", "alignment",
)

_FLOAT_COLUMNS = frozenset({
    "overall_tm_c", "limiting_arm_tm_c", "arm3_tm_c", "arm5_tm_c", "tm_c",
})


def _cell(value: object, column: str) -> object:
    """One formatting rule: floats to 3 dp, newlines escaped onto one line."""
    if column in _FLOAT_COLUMNS:
        return f"{float(value):.3f}"
    if isinstance(value, str):
        return value.replace("\n", "\\n").replace("\t", "\\t")
    return value


def dimer_row(dimer: LigationDimer) -> Dict[str, object]:
    """One dimer as a ``{column: value}`` mapping, formatted for a report."""
    return {
        name: _cell(getattr(dimer, attr), name) for name, attr in DIMER_COLUMNS
    }


def write_dimer_report(
    dimers: Sequence[LigationDimer], tsv_path: Path | str,
) -> None:
    """Write dimers to TSV, highest Tm first. Empty still writes the header."""
    _write_tsv(
        tsv_path, REPORT_COLUMNS,
        (dimer_row(d) for d in sorted(
            dimers, key=lambda x: x.overall_tm_c, reverse=True,
        )),
    )


def write_self_circ_report(
    hits: Sequence[SelfCircularisationHit], tsv_path: Path | str,
) -> None:
    """Write self-circularisation hits to TSV, highest Tm first."""
    _write_tsv(
        tsv_path, SELF_CIRC_COLUMNS,
        (
            {c: _cell(getattr(h, c), c) for c in SELF_CIRC_COLUMNS}
            for h in sorted(hits, key=lambda x: x.tm_c, reverse=True)
        ),
    )


def _write_tsv(
    tsv_path: Path | str,
    columns: Sequence[str],
    rows: Iterable[Dict[str, object]],
) -> None:
    tsv_path = Path(tsv_path)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tsv_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(list(columns))
        for row in rows:
            writer.writerow([row[c] for c in columns])


__all__ = [
    "ProbeForScreen", "Splint", "LigatorGeometry",
    "LigationDimer", "SelfCircularisationHit", "RegisterScore",
    "screen_cross_ligation", "screen_against_splints",
    "screen_self_circularisation", "screen_cross_ligation_v2",
    "ligation_geometry", "build_v2_geom", "build_geometry",
    "arm_offsets", "ligation_substrate", "backbone_loop_nt",
    "junction_block", "build_splint_index", "find_ligation_registers",
    "score_register", "best_register", "junction_run", "render_register",
    "write_dimer_report", "write_self_circ_report", "dimer_row",
    "DIMER_COLUMNS", "REPORT_COLUMNS", "SELF_CIRC_COLUMNS",
    "LAB_HYBRIDISATION", "CROSS_LIG_NN_MODEL",
    "DEFAULT_TM_THRESHOLD_C", "DEFAULT_VICINITY_N", "MAX_CONSECUTIVE_MISMATCH",
]
