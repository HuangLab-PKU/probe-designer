"""High-level v2 cross-ligation entry point used by the design CLIs.

Takes a list of candidate probes (newly designed in this run) and an
optional list of splint probes (from an external pool), runs the v2
asymmetric screen (:func:`screen_cross_ligation_v2`), and returns hits
filtered to only those involving at least one candidate.

**This module is bank-free** — no ``probe_book`` imports. The CLI layer
(``probe_designer.cli.assemble`` etc.) is responsible for loading the
external pool's probes via ``probe_designer.ext.pool.loader`` and
passing them in as ``splint_probes``.

Caller responsibility: turn :class:`CrossLigHit` records into per-probe
annotations (the ``cross_lig_partners`` xlsx column) and/or drop flagged
probes from the output.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from probe_designer.qc.cross_ligation import (
    DEFAULT_DNA_CONC_M,
    DEFAULT_DV_CONC_MM,
    DEFAULT_MV_CONC_MM,
    DEFAULT_TEMP_C,
    DEFAULT_TM_THRESHOLD_C,
    DEFAULT_VICINITY_N,
    LigationDimer,
    ProbeForScreen,
    screen_cross_ligation_v2,
)


@dataclass
class CandidateProbe:
    """A new probe being designed — possibly without an assigned probe_id yet.

    This is the design-CLI-facing input type. Converted internally to
    :class:`ProbeForScreen`. ``sequence`` should be the FULL assembled
    probe (arm5 + bb + arm3 [+ flap for iLock]) so the splint-side
    representation has bb available.
    """
    probe_id: str
    chemistry: str           # "iLock" | "dRNA" | "cDNA"
    probe_arm5: str
    probe_arm3: str
    sequence: str = ""       # full assembled probe — REQUIRED for v2 splint side
    target: str = ""         # gene_name; informational


@dataclass
class CrossLigHit:
    """One v2.2 confirmed cross-lig direction (a-as-ligator on b-as-splint).

    v2.2 split-fragment schema: arm3 and arm5 hybridize to splint
    INDEPENDENTLY, so each has its own primer3 alignment, Tm, and dG.
    """
    probe_a_id: str          # ligator
    probe_b_id: str          # splint
    direction: str = "a_lig_on_b"
    overall_tm_c: float = 0.0
    arm3_tm_c: float = 0.0
    arm3_dg_kcal: float = 0.0
    arm5_tm_c: float = 0.0
    arm5_dg_kcal: float = 0.0
    a_can_ligate_on_b: bool = False
    vicinity_contiguous: bool = False
    vicinity_n_each_side: int = 0
    a_target: str = ""
    b_target: str = ""
    a_is_existing_pool: bool = False
    b_is_existing_pool: bool = False
    b_3oh_pos: Optional[int] = None
    b_5p_pos: Optional[int] = None
    arm3_dimer_structure: str = ""
    arm5_dimer_structure: str = ""

    def partner_id(self, me: str) -> str:
        return self.probe_b_id if self.probe_a_id == me else self.probe_a_id

    def as_short_str(self) -> str:
        return (f"{self.probe_b_id}(Tm={self.overall_tm_c:.1f}°C,"
                f"arm3={self.arm3_tm_c:.1f},arm5={self.arm5_tm_c:.1f})")


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
    mv_conc_mm: float = DEFAULT_MV_CONC_MM,
    dv_conc_mm: float = DEFAULT_DV_CONC_MM,
    dna_conc_m: float = DEFAULT_DNA_CONC_M,
    temp_c: float = DEFAULT_TEMP_C,
    vicinity_n_each_side: int = DEFAULT_VICINITY_N,
) -> Tuple[List[CrossLigHit], Dict[str, List[CrossLigHit]]]:
    """Run v2 cross-lig screen on candidates ± an optional external splint pool.

    Args:
        candidates: probes being designed in this run.
        splint_probes: existing probes from an external pool the new candidates
            will be mixed with downstream. Each must already be a
            :class:`ProbeForScreen` (use ``ext.pool.loader`` to convert from
            bank's ``EffectiveProbe``). When None or empty, only within-batch
            screening is performed.
        tm_threshold_c / buffer params: passed through to ``screen_cross_ligation_v2``.

    Returns:
        ``(all_hits, by_candidate)`` —

        * ``all_hits``: every confirmed v2 cross-lig pair where at least one
          endpoint is a candidate. Pool-vs-pool pairs are dropped. Each
          unordered candidate-pair may have up to two hits (one per
          direction).
        * ``by_candidate``: ``{candidate_probe_id: [CrossLigHit, ...]}`` for
          quick lookup when populating per-probe annotation columns.
    """
    if not candidates:
        return [], {}

    cand_ids = {c.probe_id for c in candidates}
    cand_pfs = [_candidate_to_pfs(c) for c in candidates]
    pool_pfs = list(splint_probes) if splint_probes else []
    pool_ids = {p.probe_id for p in pool_pfs}

    all_probes = cand_pfs + pool_pfs
    _, dimers = screen_cross_ligation_v2(
        all_probes,
        tm_threshold_c=tm_threshold_c,
        mv_conc_mm=mv_conc_mm, dv_conc_mm=dv_conc_mm,
        dna_conc_m=dna_conc_m, temp_c=temp_c,
        vicinity_n_each_side=vicinity_n_each_side,
    )

    hits: List[CrossLigHit] = []
    by_cand: Dict[str, List[CrossLigHit]] = {}
    for d in dimers:
        # Confirmed = high-Tm AND ligation geometry AND vicinity-contiguous
        # (the third clause is a no-op when vicinity_n_each_side == 0)
        is_confirmed = d.flagged_overall and d.a_can_ligate_on_b and (
            d.vicinity_contiguous if vicinity_n_each_side > 0 else True
        )
        if not is_confirmed:
            continue
        a_is_pool = d.seq_a_id in pool_ids
        b_is_pool = d.seq_b_id in pool_ids
        if a_is_pool and b_is_pool:
            continue  # pool-vs-pool not this design run's concern
        hit = CrossLigHit(
            probe_a_id=d.seq_a_id, probe_b_id=d.seq_b_id,
            direction="a_lig_on_b",
            overall_tm_c=d.overall_tm_c,
            arm3_tm_c=d.arm3_tm_c, arm3_dg_kcal=d.arm3_dg_kcal,
            arm5_tm_c=d.arm5_tm_c, arm5_dg_kcal=d.arm5_dg_kcal,
            a_can_ligate_on_b=d.a_can_ligate_on_b,
            vicinity_contiguous=d.vicinity_contiguous,
            vicinity_n_each_side=d.vicinity_n_each_side,
            a_target=d.a_target, b_target=d.b_target,
            a_is_existing_pool=a_is_pool, b_is_existing_pool=b_is_pool,
            b_3oh_pos=d.b_3oh_pos, b_5p_pos=d.b_5p_pos,
            arm3_dimer_structure=d.arm3_dimer_structure,
            arm5_dimer_structure=d.arm5_dimer_structure,
        )
        hits.append(hit)
        # Index by candidate endpoint(s) for annotation lookup
        if not a_is_pool:
            by_cand.setdefault(d.seq_a_id, []).append(hit)
        if not b_is_pool:
            by_cand.setdefault(d.seq_b_id, []).append(hit)

    return hits, by_cand


__all__ = [
    "CandidateProbe", "CrossLigHit",
    "screen_candidates",
    "DEFAULT_TM_THRESHOLD_C", "DEFAULT_MV_CONC_MM",
    "DEFAULT_DV_CONC_MM", "DEFAULT_DNA_CONC_M", "DEFAULT_TEMP_C",
]
