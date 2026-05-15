"""Markdown + JSON rendering for ExpressionQC."""
from __future__ import annotations

import json
from pathlib import Path
from typing import List

from probe_designer.qc.expression import ExpressionQC


def render_panel_qc_markdown(qc: ExpressionQC, *, panel_id: str = "") -> str:
    """Return a markdown report for the QC result.

    Layout:
      - Header (panel ID, cell-type count, gene count)
      - Flags summary (grouped by kind)
      - Panel utilization per cell type (sorted desc)
      - Co-expressed EHEG pairs
    """
    lines: List[str] = []
    title = f"Panel expression QC — {panel_id}" if panel_id else "Panel expression QC"
    lines.append(f"# {title}")
    lines.append("")
    lines.append(f"- panel genes: **{len(qc.panel_genes)}**")
    lines.append(f"- cell types in reference: **{len(qc.cell_types)}**")
    lines.append("")

    # Flag counts
    counts: dict[str, int] = {}
    for f in qc.flags:
        counts[f.kind] = counts.get(f.kind, 0) + 1
    lines.append("## Flags")
    if not qc.flags:
        lines.append("No flags raised.")
    else:
        for kind, n in sorted(counts.items()):
            lines.append(f"- **{kind}**: {n}")
        lines.append("")
        lines.append("### Detail")
        for f in qc.flags:
            ct = f"@{f.cell_type}" if f.cell_type else ""
            lines.append(
                f"- `{f.gene}` {ct} — {f.kind}: {f.note}"
            )
    lines.append("")

    # Panel utilization per cell type
    lines.append("## Panel TP10K per cell type")
    sorted_cts = sorted(
        qc.panel_tp10k_per_cell_type.items(),
        key=lambda kv: kv[1], reverse=True,
    )
    for ct, val in sorted_cts:
        lines.append(f"- {ct}: {val:.1f}")
    lines.append("")

    # EHEG pairs
    lines.append("## Co-expressed EHEG pairs")
    if not qc.ehec_pairs:
        lines.append("None.")
    else:
        for a, b, j in qc.ehec_pairs:
            lines.append(f"- `{a}` ↔ `{b}` (cell-type overlap Jaccard = {j:.2f})")
    lines.append("")

    return "\n".join(lines) + "\n"


def write_panel_qc_json(qc: ExpressionQC, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(qc.to_dict(), f, indent=2, ensure_ascii=False)
