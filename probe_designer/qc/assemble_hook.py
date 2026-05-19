"""Cross-ligation hook for the assembled-probes DataFrame.

The ``probe-design assemble`` CLI calls this AFTER ``ProbeAssembler``
produces the final probes DataFrame and BEFORE writing the xlsx. The
hook:

1. Converts the DataFrame to a list of :class:`CandidateProbe` records.
2. Runs :func:`screen_candidates_against_pool` (within-batch + optional
   pool screen).
3. Returns three DataFrames — annotated probes, dropped (if reject), and
   the full cross-lig report.

Callers decide what to do with the three frames (write to disk, present
to the user, etc.). This separation lets us reuse the same hook for the
mutation + TCR design CLIs.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from probe_designer.qc.cross_lig_check import (
    CandidateProbe,
    CrossLigHit,
    screen_candidates_against_pool,
    DEFAULT_TM_THRESHOLD_C,
    DEFAULT_END_DG_THRESHOLD_KCAL,
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
    """Compact ``other:Tm=XX:dG=YY`` strings, ;-separated."""
    parts = []
    for h in hits:
        other = h.probe_b_id if h.probe_a_id == me else h.probe_a_id
        tag = ":POOL" if (h.a_is_existing_pool or h.b_is_existing_pool) else ""
        parts.append(f"{other}{tag}:Tm={h.overall_tm_c:.1f}:dG={h.overall_dg_kcal:.1f}")
    return ";".join(sorted(set(parts)))


def _hits_to_report_df(hits: Iterable[CrossLigHit]) -> pd.DataFrame:
    rows = [{
        "probe_a_id": h.probe_a_id, "probe_b_id": h.probe_b_id,
        "probe_a_gene": h.a_target, "probe_b_gene": h.b_target,
        "kind": h.kind,
        "overall_tm_c": h.overall_tm_c,
        "overall_dg_kcal": h.overall_dg_kcal,
        "a_end_dg_kcal": h.a_end_dg_kcal,
        "b_end_dg_kcal": h.b_end_dg_kcal,
        "flagged_end_a": h.flagged_end_a,
        "flagged_end_b": h.flagged_end_b,
        "a_is_existing_pool": h.a_is_existing_pool,
        "b_is_existing_pool": h.b_is_existing_pool,
    } for h in hits]
    return pd.DataFrame(rows).sort_values(
        "overall_tm_c", ascending=False, ignore_index=True,
    ) if rows else pd.DataFrame(columns=[
        "probe_a_id", "probe_b_id", "probe_a_gene", "probe_b_gene", "kind",
        "overall_tm_c", "overall_dg_kcal", "a_end_dg_kcal", "b_end_dg_kcal",
        "flagged_end_a", "flagged_end_b", "a_is_existing_pool", "b_is_existing_pool",
    ])


def apply_cross_lig_check(
    probes_df: pd.DataFrame,
    *,
    pool_id: Optional[str] = None,
    repo_root: Optional[Path] = None,
    reject: bool = False,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    end_dg_threshold_kcal: float = DEFAULT_END_DG_THRESHOLD_KCAL,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run cross-lig screen on the assembled probes; annotate / optionally reject.

    Returns:
        ``(annotated_df, dropped_df, report_df)``.

        * ``annotated_df`` always contains every row from ``probes_df``
          plus a new ``cross_lig_partners`` column with semicolon-separated
          ``<partner_id>[:POOL]:Tm=XX.X:dG=YY.Y`` summaries (empty when
          there's no flagged partner). When ``reject=True``, the
          flagged rows have been removed.
        * ``dropped_df`` contains the removed rows when ``reject=True``
          (empty otherwise). The same cross_lig_partners column is set.
        * ``report_df`` contains every confirmed pair (the deduplicated
          screen output), sorted by Tm descending.
    """
    cands = probes_df_to_candidates(probes_df)
    hits, by_cand = screen_candidates_against_pool(
        cands,
        pool_id=pool_id,
        repo_root=repo_root,
        tm_threshold_c=tm_threshold_c,
        end_dg_threshold_kcal=end_dg_threshold_kcal,
    )

    annotated = probes_df.copy()
    annotated["cross_lig_partners"] = annotated["probe_name"].map(
        lambda pid: _format_partner_summary(by_cand.get(str(pid), []), str(pid))
    )

    report_df = _hits_to_report_df(hits)

    if reject and by_cand:
        flagged_mask = annotated["cross_lig_partners"] != ""
        dropped_df = annotated[flagged_mask].copy().reset_index(drop=True)
        annotated = annotated[~flagged_mask].copy().reset_index(drop=True)
    else:
        dropped_df = annotated.iloc[0:0].copy()

    return annotated, dropped_df, report_df


def apply_cross_lig_check_to_xlsx(
    xlsx_path: Path,
    *,
    pool_id: Optional[str] = None,
    repo_root: Optional[Path] = None,
    reject: bool = False,
    tm_threshold_c: float = DEFAULT_TM_THRESHOLD_C,
    end_dg_threshold_kcal: float = DEFAULT_END_DG_THRESHOLD_KCAL,
) -> tuple[int, int]:
    """In-place screen on an existing probes xlsx (mutation / TCR post-pipeline).

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
        pool_id=pool_id,
        repo_root=repo_root,
        reject=reject,
        tm_threshold_c=tm_threshold_c,
        end_dg_threshold_kcal=end_dg_threshold_kcal,
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
