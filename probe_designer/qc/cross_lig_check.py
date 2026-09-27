"""High-level cross-ligation entry point used by the design CLIs.

Takes a list of candidate probes (newly designed in this run) and an optional
list of splint probes (from an external pool), runs the register scan
(:func:`screen_cross_ligation`), and returns hits involving at least one
candidate.

**This module is bank-free** — no ``probe_book`` imports. The CLI layer
(``probe_designer.cli.assemble`` etc.) loads the external pool's probes via
``probe_designer.ext.pool.loader`` and passes them in as ``splint_probes``.

Caller responsibility: turn :class:`CrossLigHit` records into per-probe
annotations (the ``cross_lig_partners`` xlsx column) and/or drop flagged probes
from the output.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from probe_designer.chemistry import ReactionConditions
from probe_designer.qc.cross_ligation import (
    DEFAULT_TM_THRESHOLD_C,
    DEFAULT_VICINITY_N,
    LAB_HYBRIDISATION,
    LigationDimer,
    ProbeForScreen,
    screen_cross_ligation,
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
    """A confirmed register, plus where its two ends came from.

    Composed rather than copied: the screen's own :class:`LigationDimer` is held
    whole and read through. Mirroring its fields here meant every number added to
    the model had to be written into four places — the dataclass, this one, the
    copy block, and the report schema — and a field added to three of the four is
    a silently empty report column.
    """
    dimer: LigationDimer
    a_is_existing_pool: bool = False
    b_is_existing_pool: bool = False

    # -- read-through, so callers need not reach into .dimer for the basics --
    @property
    def probe_a_id(self) -> str:
        return self.dimer.seq_a_id

    @property
    def probe_b_id(self) -> str:
        return self.dimer.seq_b_id

    @property
    def a_target(self) -> str:
        return self.dimer.a_target

    @property
    def b_target(self) -> str:
        return self.dimer.b_target

    @property
    def overall_tm_c(self) -> float:
        return self.dimer.overall_tm_c

    @property
    def is_self_pair(self) -> bool:
        return self.dimer.is_self_pair

    def partner_id(self, me: str) -> str:
        return self.probe_b_id if self.probe_a_id == me else self.probe_a_id

    def as_short_str(self) -> str:
        return (f"{self.probe_b_id}(Tm={self.overall_tm_c:.1f}C,"
                f"junction={self.dimer.junction_run_nt}nt)")


def _candidate_to_pfs(c: CandidateProbe) -> ProbeForScreen:
    return ProbeForScreen(
        probe_id=c.probe_id, chemistry=c.chemistry,
        probe_arm5=c.probe_arm5, probe_arm3=c.probe_arm3,
        sequence=c.sequence or (c.probe_arm5 + c.probe_arm3),  # fall back if missing
        target=c.target,
    )


def screen_candidates(
    candidates: Sequence[CandidateProbe], *,
    splint_probes: Optional[Sequence[ProbeForScreen]] = None,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    reaction: ReactionConditions = LAB_HYBRIDISATION,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
) -> Tuple[List[CrossLigHit], Dict[str, List[CrossLigHit]]]:
    """Run the cross-lig screen on candidates, optionally against a pool.

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

        * ``all_hits``: every confirmed direction where at least one endpoint is
          a candidate. Pool-vs-pool pairs are dropped, and so is a pool probe's
          *self*-pair: it is a property of the pool as shipped, not of this run.
        * ``by_candidate``: ``{candidate_probe_id: [CrossLigHit, ...]}`` for
          populating per-probe annotation columns.
    """
    if not candidates:
        return [], {}

    cand_pfs = [_candidate_to_pfs(c) for c in candidates]
    pool_pfs = list(splint_probes) if splint_probes else []
    pool_ids = {p.probe_id for p in pool_pfs}

    dimers = screen_cross_ligation(
        cand_pfs + pool_pfs,
        tm_threshold_c=tm_threshold_c,
        reaction=reaction,
        vicinity_n_each_side=vicinity_n_each_side,
    )

    hits: List[CrossLigHit] = []
    by_cand: Dict[str, List[CrossLigHit]] = {}
    for dimer in dimers:
        if not dimer.flagged_overall:
            continue
        a_is_pool = dimer.seq_a_id in pool_ids
        b_is_pool = dimer.seq_b_id in pool_ids
        if a_is_pool and b_is_pool:
            continue  # includes a pool probe's self-pair — not this run's concern
        hit = CrossLigHit(dimer, a_is_pool, b_is_pool)
        hits.append(hit)
        # Index by candidate endpoint(s). A self-pair has one endpoint, so
        # indexing both would list the probe against itself twice.
        if not a_is_pool:
            by_cand.setdefault(dimer.seq_a_id, []).append(hit)
        if not b_is_pool and not dimer.is_self_pair:
            by_cand.setdefault(dimer.seq_b_id, []).append(hit)

    return hits, by_cand


__all__ = [
    "CandidateProbe", "CrossLigHit",
    "screen_candidates",
    "DEFAULT_TM_THRESHOLD_C", "DEFAULT_VICINITY_N", "LAB_HYBRIDISATION",
]
