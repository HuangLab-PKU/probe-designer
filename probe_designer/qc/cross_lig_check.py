"""High-level cross-ligation entry point used by the design CLIs.

Takes a list of candidate probes (newly designed in this run) and an optional
list of splint probes (from an external pool), runs the v3 register scan
(:func:`screen_cross_ligation_v2`), and returns hits filtered to those
involving at least one candidate.

**This module is bank-free** — no ``probe_book`` imports. The CLI layer
(``probe_designer.cli.assemble`` etc.) loads the external pool's probes via
``probe_designer.ext.pool.loader`` and passes them in as ``splint_probes``.

Caller responsibility: turn :class:`CrossLigHit` records into per-probe
annotations (the ``cross_lig_partners`` xlsx column) and/or drop flagged probes
from the output.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from probe_designer.chemistry import ReactionConditions
from probe_designer.qc.cross_ligation import (
    DEFAULT_TM_THRESHOLD_C,
    DEFAULT_VICINITY_N,
    LAB_HYBRIDISATION,
    ProbeForScreen,
    screen_cross_ligation_v2,
)


@dataclass
class CandidateProbe:
    """A new probe being designed — possibly without an assigned probe_id yet.

    This is the design-CLI-facing input type. Converted internally to
    :class:`ProbeForScreen`. ``sequence`` should be the FULL assembled probe
    (arm5 + bb + arm3 [+ flap for iLock]) so the splint-side representation has
    the backbone available — cross-ligation on a shared backbone is real.
    """
    probe_id: str
    chemistry: str           # "iLock" | "dRNA" | "cDNA"
    probe_arm5: str
    probe_arm3: str
    sequence: str = ""       # full assembled probe — REQUIRED for the splint side
    target: str = ""         # gene_name; informational


@dataclass
class CrossLigHit:
    """One confirmed cross-ligation direction (a-as-ligator on b-as-splint).

    v3 schema: the numbers describe the best ligation-competent **register**
    found on the splint, not two independently-chosen alignments.
    ``overall_tm_c`` — the duplex spanning the nick — is the gate.
    ``limiting_arm_tm_c`` is the weaker of the two halves of that same
    alignment, a diagnostic for how lopsided the register is; a short half
    melts below 0 C and reports a large negative value.
    """
    probe_a_id: str          # ligator
    probe_b_id: str          # splint
    direction: str = "a_lig_on_b"
    overall_tm_c: float = 0.0
    limiting_arm_tm_c: float = 0.0
    arm3_tm_c: float = 0.0
    arm5_tm_c: float = 0.0
    junction_run_nt: int = 0
    paired_nt: int = 0
    a_can_ligate_on_b: bool = False
    vicinity_contiguous: bool = False
    vicinity_n_each_side: int = 0
    is_self_pair: bool = False
    a_target: str = ""
    b_target: str = ""
    a_is_existing_pool: bool = False
    b_is_existing_pool: bool = False
    nick_pos_on_b: int = -1
    b_3oh_pos: Optional[int] = None
    b_5p_pos: Optional[int] = None
    alignment: str = ""

    def partner_id(self, me: str) -> str:
        return self.probe_b_id if self.probe_a_id == me else self.probe_a_id

    def as_short_str(self) -> str:
        return (f"{self.probe_b_id}(Tm={self.overall_tm_c:.1f}C,"
                f"junction={self.junction_run_nt}nt)")


def _candidate_to_pfs(c: CandidateProbe) -> ProbeForScreen:
    return ProbeForScreen(
        probe_id=c.probe_id, chemistry=c.chemistry,
        probe_arm5=c.probe_arm5, probe_arm3=c.probe_arm3,
        sequence=c.sequence or (c.probe_arm5 + c.probe_arm3),  # fall back if missing
        target=c.target,
    )


def screen_candidates(
    candidates: List[CandidateProbe], *,
    splint_probes: Optional[List[ProbeForScreen]] = None,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
) -> Tuple[List[CrossLigHit], Dict[str, List[CrossLigHit]]]:
    """Run the v3 cross-lig screen on candidates, optionally against a pool.

    Args:
        candidates: probes being designed in this run.
        splint_probes: existing probes from an external pool the candidates will
            be mixed with downstream. Each must already be a
            :class:`ProbeForScreen` (use ``ext.pool.loader`` to convert from
            bank's ``EffectiveProbe``). When None or empty, only within-batch
            screening is performed.
        tm_threshold_c: flag a register whose duplex Tm clears this.
        reaction: hybridisation conditions; defaults to the lab protocol.
        vicinity_n_each_side: ligase clamp width, in contiguous WC pairs.

    Returns:
        ``(all_hits, by_candidate)`` —

        * ``all_hits``: every confirmed cross-lig direction where at least one
          endpoint is a candidate. Pool-vs-pool pairs are dropped. A pool
          probe's *self*-pair is likewise dropped: it is a property of the pool
          as shipped, not of this design run.
        * ``by_candidate``: ``{candidate_probe_id: [CrossLigHit, ...]}`` for
          populating per-probe annotation columns.
    """
    if not candidates:
        return [], {}

    cand_pfs = [_candidate_to_pfs(c) for c in candidates]
    pool_pfs = list(splint_probes) if splint_probes else []
    pool_ids = {p.probe_id for p in pool_pfs}

    _, dimers = screen_cross_ligation_v2(
        cand_pfs + pool_pfs,
        tm_threshold_c=tm_threshold_c,
        reaction=reaction,
        vicinity_n_each_side=vicinity_n_each_side,
    )

    hits: List[CrossLigHit] = []
    by_cand: Dict[str, List[CrossLigHit]] = {}
    for d in dimers:
        if not (d.flagged_overall and d.a_can_ligate_on_b):
            continue
        a_is_pool = d.seq_a_id in pool_ids
        b_is_pool = d.seq_b_id in pool_ids
        if a_is_pool and b_is_pool:
            continue  # includes a pool probe's self-pair — not this run's concern
        hit = CrossLigHit(
            probe_a_id=d.seq_a_id, probe_b_id=d.seq_b_id,
            direction="a_lig_on_b",
            overall_tm_c=d.overall_tm_c,
            limiting_arm_tm_c=d.limiting_arm_tm_c,
            arm3_tm_c=d.arm3_tm_c, arm5_tm_c=d.arm5_tm_c,
            junction_run_nt=d.junction_run_nt, paired_nt=d.paired_nt,
            a_can_ligate_on_b=d.a_can_ligate_on_b,
            vicinity_contiguous=d.vicinity_contiguous,
            vicinity_n_each_side=d.vicinity_n_each_side,
            is_self_pair=d.is_self_pair,
            a_target=d.a_target, b_target=d.b_target,
            a_is_existing_pool=a_is_pool, b_is_existing_pool=b_is_pool,
            nick_pos_on_b=d.nick_pos_on_b,
            b_3oh_pos=d.b_3oh_pos, b_5p_pos=d.b_5p_pos,
            alignment=d.alignment,
        )
        hits.append(hit)
        # Index by candidate endpoint(s). A self-pair has one endpoint, so
        # indexing both would list the probe against itself twice.
        if not a_is_pool:
            by_cand.setdefault(d.seq_a_id, []).append(hit)
        if not b_is_pool and not d.is_self_pair:
            by_cand.setdefault(d.seq_b_id, []).append(hit)

    return hits, by_cand


__all__ = [
    "CandidateProbe", "CrossLigHit",
    "screen_candidates",
    "DEFAULT_TM_THRESHOLD_C", "DEFAULT_VICINITY_N", "LAB_HYBRIDISATION",
]
