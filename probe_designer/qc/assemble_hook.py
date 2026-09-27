"""Cross-ligation hook for the assembled-probes DataFrame.

The design CLIs (``probe-design assemble`` and the mutation/tcr CLIs)
call this hook AFTER the pipeline produces the final probes DataFrame
and BEFORE writing the xlsx. The hook:

1. Converts the DataFrame to a list of :class:`CandidateProbe` records.
2. Runs :func:`screen_candidates` (within-batch + optional external splint pool).
3. Returns three DataFrames — annotated, dropped (if reject=True), and
   the full v2 dimer report.

The hook is **bank-free** — pool-loading happens in the CLI layer (via
``probe_designer.ext.pool.loader``) and the resolved splint probes are
passed in as ``splint_probes``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd

from probe_designer.qc.cross_lig_check import (
    CandidateProbe,
    CrossLigHit,
    DEFAULT_TM_THRESHOLD_C,
    screen_candidates,
)
from probe_designer.qc.cross_ligation import (
    REPORT_COLUMNS as DIMER_REPORT_COLUMNS,
    ProbeForScreen,
    dimer_row,
)


def probes_df_to_candidates(df: pd.DataFrame) -> list[CandidateProbe]:
    """Pull arm5/arm3/chemistry/sequence/target from a Schema-v2 probes DataFrame."""
    out: list[CandidateProbe] = []
    for _, r in df.iterrows():
        arm5 = str(r.get("probe_arm5") or "")
        arm3 = str(r.get("probe_arm3") or "")
        if not arm5 or not arm3:
            continue
        out.append(CandidateProbe(
            probe_id=str(r["probe_name"]),
            chemistry=str(r.get("chemistry") or "dRNA"),
            probe_arm5=arm5,
            probe_arm3=arm3,
            sequence=str(r.get("probe_sequence") or ""),
            target=str(r.get("gene_name") or r.get("target") or ""),
        ))
    return out


def _format_partner_summary(hits: list[CrossLigHit], me: str) -> str:
    """Compact ``<other>[:POOL]:Tm=XX.X:lig=A|B|AB`` strings, ;-separated.

    Each hit is a one-direction (a-as-ligator → b-as-splint) confirmed v2 dimer.
    The ``lig=`` token indicates whether ``me`` is the ligator (A), splint (B),
    or both directions for the same partner — collapsed across hits.
    """
    # Group hits by partner; record whether me is ligator and/or splint
    partner_info: dict[str, dict] = {}
    for h in hits:
        other = h.partner_id(me)
        is_pool_partner = (
            (h.probe_a_id == other and h.a_is_existing_pool)
            or (h.probe_b_id == other and h.b_is_existing_pool)
        )
        info = partner_info.setdefault(other, {
            "is_pool": is_pool_partner,
            "me_is_ligator": False, "me_is_splint": False,
            "max_tm": h.overall_tm_c,
        })
        info["max_tm"] = max(info["max_tm"], h.overall_tm_c)
        if h.probe_a_id == me:
            info["me_is_ligator"] = True
        else:
            info["me_is_splint"] = True

    parts: list[str] = []
    for partner, info in partner_info.items():
        if info["me_is_ligator"] and info["me_is_splint"]:
            lig = "AB"
        elif info["me_is_ligator"]:
            lig = "A"
        else:
            lig = "B"
        tag = ":POOL" if info["is_pool"] else ""
        parts.append(f"{partner}{tag}:Tm={info['max_tm']:.1f}:lig={lig}")
    return ";".join(sorted(parts))


#: The screen's own columns, plus the two this layer adds. Renamed at the
#: DataFrame boundary only because the xlsx has always said ``probe_*`` / gene;
#: the values come from ``cross_ligation.dimer_row`` so the two report writers
#: cannot drift.
_RENAMES = {
    "seq_a_id": "probe_a_id", "seq_b_id": "probe_b_id",
    "a_target": "probe_a_gene", "b_target": "probe_b_gene",
}
REPORT_COLUMNS: tuple[str, ...] = tuple(
    [_RENAMES.get(c, c) for c in DIMER_REPORT_COLUMNS]
    + ["a_is_existing_pool", "b_is_existing_pool"]
)


def _hits_to_report_df(hits: Iterable[CrossLigHit]) -> pd.DataFrame:
    rows = []
    for h in hits:
        row = {_RENAMES.get(k, k): v for k, v in dimer_row(h.dimer).items()}
        row["a_is_existing_pool"] = h.a_is_existing_pool
        row["b_is_existing_pool"] = h.b_is_existing_pool
        rows.append(row)
    if not rows:
        return pd.DataFrame(columns=list(REPORT_COLUMNS))
    return pd.DataFrame(rows).sort_values(
        "overall_tm_c", ascending=False, ignore_index=True,
    )


def apply_cross_lig_check(
    probes_df: pd.DataFrame,
    *,
    splint_probes: Optional[List[ProbeForScreen]] = None,
    reject: bool = False,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run the cross-lig screen on assembled probes; annotate / optionally reject.

    Args:
        probes_df: assembled probes (Schema-v2 with probe_arm5, probe_arm3,
            probe_sequence, chemistry, probe_name, gene_name columns).
        splint_probes: optional list of :class:`ProbeForScreen` from an
            external pool to screen candidates against. Loaded by the CLI
            layer (``ext.pool.loader``); this function does not import bank.
        reject: when True, probes that act as the **ligator** in a confirmed
            hit are removed from ``annotated_df`` and returned in
            ``dropped_df``. Only the ligator is dropped: it is the probe that
            circularises and produces target-independent signal, while the
            splint merely templated it and is otherwise fine. Dropping both
            ends of every pair (what this did before v3) discarded twice the
            probes needed to break the interaction. Every probe involved is
            still annotated, ligator or not.
        tm_threshold_c: flag a register whose duplex Tm clears this.

    Returns ``(annotated_df, dropped_df, report_df)`` — see module docstring.
    """
    cands = probes_df_to_candidates(probes_df)
    hits, by_cand = screen_candidates(
        cands,
        splint_probes=splint_probes,
        tm_threshold_c=tm_threshold_c,
    )

    annotated = probes_df.copy()
    annotated["cross_lig_partners"] = annotated["probe_name"].map(
        lambda pid: _format_partner_summary(by_cand.get(str(pid), []), str(pid))
    )

    if reject and hits:
        ligator_ids = {h.probe_a_id for h in hits}
        flagged_mask = annotated["probe_name"].astype(str).isin(ligator_ids)
        dropped_df = annotated[flagged_mask].copy().reset_index(drop=True)
        annotated = annotated[~flagged_mask].copy().reset_index(drop=True)
    else:
        dropped_df = annotated.iloc[0:0].copy()

    report_df = _hits_to_report_df(hits)
    return annotated, dropped_df, report_df


def apply_cross_lig_check_to_xlsx(
    xlsx_path: Path,
    *,
    splint_probes: Optional[List[ProbeForScreen]] = None,
    reject: bool = False,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
) -> tuple[int, int]:
    """In-place v2 screen on an existing probes xlsx (mutation / TCR post-pipeline).

    Reads ``xlsx_path``, applies :func:`apply_cross_lig_check`, then writes:

    * the annotated DataFrame back to ``xlsx_path`` (overwrites);
    * ``cross_lig_report.tsv`` in the same directory;
    * ``dropped_cross_lig.tsv`` if ``reject=True`` and probes were removed.

    Returns ``(n_hits, n_dropped)``.
    """
    xlsx_path = Path(xlsx_path)
    df = pd.read_excel(xlsx_path)
    annotated, dropped, report = apply_cross_lig_check(
        df,
        splint_probes=splint_probes,
        reject=reject,
        tm_threshold_c=tm_threshold_c,
    )
    out_dir = xlsx_path.parent
    report.to_csv(out_dir / "cross_lig_report.tsv", sep="\t", index=False)
    if len(dropped) > 0:
        dropped.to_csv(out_dir / "dropped_cross_lig.tsv", sep="\t", index=False)
    annotated.to_excel(xlsx_path, index=False)
    return len(report), len(dropped)


__all__ = [
    "probes_df_to_candidates",
    "apply_cross_lig_check",
    "apply_cross_lig_check_to_xlsx",
]
