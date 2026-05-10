"""Tm landscape PDF for TCR runs — one subplot per clone.

Plots both mRNA and cDNA per-arm Tms across the scanned positions, shades the
allowed sub-region (CDR3 for TCR), and marks selected sites with vertical
lines and triangles.

Self-contained — uses matplotlib (already a transitive dep). Skipped silently
if matplotlib is not importable.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


def plot_tm_landscape(
    all_sites: Dict[str, List[Dict[str, Any]]],
    selected_map: Dict[str, List[Dict[str, Any]]],
    *,
    pdf_path: Path,
    tm_min: float = 50.0,
    tm_max: float = 60.0,
    clone_metadata: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Optional[Path]:
    """Render a multi-page PDF: 2 plots per page, one per clone.

    Parameters
    ----------
    all_sites
        ``{clone_name: [site_dict, ...]}`` — output of :func:`TcrProbeDesigner.find_binding_sites`.
    selected_map
        ``{clone_name: [selected_site_dict, ...]}`` — sites that survived the
        Tm filter and non-overlap selection.
    pdf_path
        Where to write the PDF.
    tm_min, tm_max
        mRNA Tm gate (visualized as a horizontal shaded band).
    clone_metadata
        Optional ``{clone_name: {"clone_id": ..., ...}}`` to put in the
        subplot title. Falls back to ``clone_name`` if absent.

    Returns the PDF path on success, or ``None`` if matplotlib is unavailable
    or rendering fails. Logs a warning in either case.
    """
    try:
        import matplotlib  # noqa: F401
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as e:
        logger.warning("matplotlib not available; skipping Tm landscape PDF (%s)", e)
        return None

    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    clone_metadata = clone_metadata or {}

    COLS, ROWS = 1, 2
    PLOTS_PER_PAGE = COLS * ROWS

    try:
        with PdfPages(str(pdf_path)) as pdf:
            fig = None
            plot_idx = 0
            for clone_name in sorted(all_sites):
                sites = all_sites[clone_name]
                if not sites:
                    continue
                if plot_idx % PLOTS_PER_PAGE == 0:
                    if fig is not None:
                        fig.tight_layout(rect=[0, 0, 1, 0.96])
                        pdf.savefig(fig)
                        plt.close(fig)
                    fig, axes = plt.subplots(ROWS, COLS, figsize=(12, 10))
                    fig.suptitle("TCR Arm Tm Landscape (mRNA vs cDNA)", fontsize=14)
                    if not hasattr(axes, "__len__"):
                        axes = [axes]
                    axes = list(axes)
                ax = axes[plot_idx % PLOTS_PER_PAGE]

                cdr3_start = sites[0]["cdr3_start"]
                cdr3_end = sites[0]["cdr3_end"]
                positions = [s["ligation_point"] for s in sites]
                ax.plot(positions, [s["tm_3prime_dRNA"] for s in sites], "b.-",
                        markersize=3, linewidth=1, label="3' Tm dRNA", alpha=0.8)
                ax.plot(positions, [s["tm_5prime_dRNA"] for s in sites], "r.-",
                        markersize=3, linewidth=1, label="5' Tm dRNA", alpha=0.8)
                ax.plot(positions, [s["tm_3prime_cDNA"] for s in sites], "b.--",
                        markersize=2, linewidth=0.7, label="3' Tm cDNA", alpha=0.5)
                ax.plot(positions, [s["tm_5prime_cDNA"] for s in sites], "r.--",
                        markersize=2, linewidth=0.7, label="5' Tm cDNA", alpha=0.5)
                ax.axhspan(tm_min, tm_max, alpha=0.12, color="green",
                           label=f"dRNA Tm gate [{tm_min:.0f},{tm_max:.0f}]")
                ax.axvspan(cdr3_start, cdr3_end, alpha=0.08, color="orange",
                           label="CDR3 / allowed region")
                if clone_name in selected_map:
                    for s in selected_map[clone_name]:
                        lp = s["ligation_point"]
                        ax.axvline(lp, color="black", linewidth=1.2,
                                   linestyle="-", alpha=0.6)
                        ax.plot(lp, s["tm_5prime_dRNA"], "b^", markersize=8, zorder=5)
                        ax.plot(lp, s["tm_3prime_dRNA"], "rv", markersize=8, zorder=5)
                meta = clone_metadata.get(clone_name, {})
                clone_id = meta.get("clone_id", clone_name)
                n_sel = len(selected_map.get(clone_name, []))
                ax.set_title(
                    f"{clone_name}\n({clone_id}, allowed={cdr3_end - cdr3_start}bp, "
                    f"{n_sel} probes)", fontsize=9,
                )
                ax.set_xlabel("Ligation point position", fontsize=8)
                ax.set_ylabel("Tm (°C)", fontsize=8)
                ax.tick_params(labelsize=7)
                ax.set_ylim(30, 85)
                if plot_idx % PLOTS_PER_PAGE == 0:
                    ax.legend(fontsize=6, loc="upper right")
                plot_idx += 1

            if fig is not None:
                fig.tight_layout(rect=[0, 0, 1, 0.96])
                pdf.savefig(fig)
                plt.close(fig)
    except Exception as e:
        logger.warning("Tm landscape PDF rendering failed: %s", e, exc_info=True)
        return None

    logger.info("Tm landscape PDF: %s", pdf_path)
    return pdf_path
