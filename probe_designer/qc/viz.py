"""Panel QC visualizations.

Three plot families:

* :func:`plot_accessibility_tracks` — one row per isoform, all chosen
  temperatures overlaid as line plots; x-axis = chromosomal coordinate
  (multi-isoform alignment via the exon structure). A bottom rug strip
  marks probe-binding sites colored by an outcome label.

* :func:`plot_accessibility_heatmap` — same data, different shape:
  rows = isoforms, x = chromosomal coordinate binned, color = mean
  p_unpaired. Easier to compare many isoforms at once.

* :func:`plot_rnafold_structure` — 2D mRNA secondary structure via
  ViennaRNA ``RNA.fold`` + ``RNA.simple_xy_coordinates``. Probe-binding
  windows can be highlighted on the structure.

All three take ready inputs (isoforms, profiles, probes) — they do NOT
load from any specific cache path. The caller is responsible for
assembling these from whatever data source (Ensembl REST cache,
experiment folders, etc.). This keeps the QC module reusable across
panels.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# ----------------------------------------------------------------------
# Plain dataclasses — the public input types
# ----------------------------------------------------------------------


@dataclass
class IsoformMeta:
    """Enough info to map a mRNA-coordinate profile back to the genome."""
    transcript_id: str
    display_name: str
    chromosome: str
    strand: int   # +1 or -1
    exons: list[tuple[int, int]]   # genomic (start, end), 1-indexed inclusive
    is_canonical: bool = False
    biotype: str = ""

    def mrna_length(self) -> int:
        return sum(e - s + 1 for s, e in self.exons)


@dataclass
class ProbeAnnotation:
    """One probe's genomic footprint + (optional) wet-lab outcome label."""
    probe_id: str
    arm5: str
    arm3: str
    genomic_start: int
    genomic_end: int
    outcome: str = "unknown"   # one of "signal", "dead", "unknown"
    extra: dict = field(default_factory=dict)


# ----------------------------------------------------------------------
# Coordinate mapping (shared by tracks + heatmap)
# ----------------------------------------------------------------------


def _exon_segments_5to3(iso: IsoformMeta) -> list[tuple[int, int]]:
    """Return exons in mRNA 5'→3' order (genomic-low-to-high for +strand;
    reversed for −strand)."""
    segs = list(iso.exons)
    if iso.strand < 0:
        segs = list(reversed(segs))
    return segs


def _mrna_to_genomic_chunks(
    profile: np.ndarray, iso: IsoformMeta,
) -> list[tuple[int, int, np.ndarray]]:
    """Slice an mRNA-coordinate profile into ``(genomic_start,
    genomic_end, chunk_left_to_right)`` per exon."""
    chunks: list[tuple[int, int, np.ndarray]] = []
    segs = _exon_segments_5to3(iso)
    consumed = 0
    for gs, ge in segs:
        n = ge - gs + 1
        chunk = profile[consumed : consumed + n]
        consumed += n
        if iso.strand < 0:
            chunk = chunk[::-1]
        chunks.append((gs, ge, chunk))
    return chunks


# ----------------------------------------------------------------------
# Common helpers
# ----------------------------------------------------------------------


def _outcome_color(outcome: str) -> str:
    return {"signal": "tab:green", "dead": "tab:red"}.get(outcome, "tab:gray")


def _draw_probe_rug(
    ax, probes: Optional[list[ProbeAnnotation]], gstart: int, gend: int, chrom: str,
) -> None:
    """Draw probe spans on a rug-strip axis."""
    ax.set_ylim(0, 1)
    ax.set_xlim(gstart, gend)
    ax.set_yticks([])
    if probes:
        for p in probes:
            if p.genomic_end < gstart or p.genomic_start > gend:
                continue
            ax.axvspan(p.genomic_start, p.genomic_end, color=_outcome_color(p.outcome),
                        alpha=0.55)
        ax.set_ylabel("probes\n(green=signal,\nred=dead)", fontsize=7)
    else:
        ax.set_ylabel("probes", fontsize=8)
    ax.set_xlabel(f"chr{chrom}: genomic position (bp)")


def _ensure_agg_backend():
    if matplotlib.get_backend().lower() not in ("agg", "module://matplotlib_inline.backend_inline"):
        # Don't override user's choice; only switch if we're in a context
        # where rendering would fail (e.g. headless).
        pass


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def plot_accessibility_tracks(
    isoforms: list[IsoformMeta],
    profiles: dict[str, dict[float, np.ndarray]],
    *,
    probes: Optional[list[ProbeAnnotation]] = None,
    temperatures: tuple[float, ...] = (37.0,),
    title: str = "",
    save_path: Optional[Path] = None,
    max_isoforms: int = 6,
    threshold_line: float = 0.30,
) -> plt.Figure:
    """Multi-isoform accessibility tracks aligned on genomic coordinate.

    Args:
        isoforms: ordered list; the first ``max_isoforms`` are plotted.
            Caller decides ordering (typically canonical first, then by
            length).
        profiles: mapping ``{transcript_id: {plfold_T: per-base
            p_unpaired}}``. Missing entries are silently skipped.
        probes: optional probe overlays; outcome attribute drives color.
        temperatures: which temperatures to render (one line per temp).
        title: figure suptitle (defaults to the gene name if you embed it).
        save_path: PNG output; if None, only returns the Figure.
        max_isoforms: cap on isoform tracks rendered.
        threshold_line: horizontal reference line at this p_unpaired
            (typical Phase 1A default 0.30). Set to None to suppress.
    """
    _ensure_agg_backend()
    isos = isoforms[:max_isoforms]
    if not isos:
        raise ValueError("isoforms list is empty")

    gstart = min(min(e[0] for e in i.exons) for i in isos)
    gend = max(max(e[1] for e in i.exons) for i in isos)
    chrom = isos[0].chromosome

    n_iso = len(isos)
    fig, axes = plt.subplots(
        n_iso + 1, 1,
        figsize=(13, max(4.5, 1.1 * n_iso + 1.2)),
        sharex=True,
        gridspec_kw={"height_ratios": [1] * n_iso + [0.55]},
    )
    if n_iso == 1:
        # axes is shape (2,) — already correctly indexable
        pass

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, max(2, len(temperatures))))[:len(temperatures)]

    for i, iso in enumerate(isos):
        ax = axes[i]
        for j, T in enumerate(temperatures):
            profile = profiles.get(iso.transcript_id, {}).get(T)
            if profile is None:
                continue
            if len(profile) != iso.mrna_length():
                continue
            chunks = _mrna_to_genomic_chunks(profile, iso)
            for k, (gs, _ge, chunk) in enumerate(chunks):
                xs = np.arange(gs, gs + len(chunk))
                ax.plot(xs, chunk, color=colors[j], linewidth=0.5, alpha=0.9,
                         label=f"T={T:g}°C" if k == 0 else None)
        ax.set_ylim(0, 1)
        ax.set_xlim(gstart, gend)
        ax.tick_params(axis="y", labelsize=7)
        if threshold_line is not None:
            ax.axhline(threshold_line, color="gray", linestyle=":", linewidth=0.5)
        label = (f"{iso.display_name}{'*' if iso.is_canonical else ''}"
                  f"\n{iso.biotype}" if iso.biotype else
                  f"{iso.display_name}{'*' if iso.is_canonical else ''}")
        ax.set_ylabel(label, fontsize=8)
        if i == 0:
            ax.legend(loc="upper right", fontsize=7, ncol=len(temperatures))

    _draw_probe_rug(axes[-1], probes, gstart, gend, chrom)
    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=130, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_accessibility_heatmap(
    isoforms: list[IsoformMeta],
    profiles: dict[str, np.ndarray],
    *,
    probes: Optional[list[ProbeAnnotation]] = None,
    temperature: float = 37.0,
    title: str = "",
    save_path: Optional[Path] = None,
    max_isoforms: int = 12,
    bin_size: int = 25,
) -> plt.Figure:
    """Multi-isoform accessibility heatmap (rows = isoforms, x = genomic bins).

    Args:
        profiles: mapping ``{transcript_id: profile}`` at the chosen
            temperature. (Single T per heatmap — choose your favorite or
            call this twice.)
        bin_size: bp per heatmap column. Larger = more cells visible per
            isoform, less resolution.
    """
    _ensure_agg_backend()
    isos = isoforms[:max_isoforms]
    if not isos:
        raise ValueError("isoforms list is empty")

    gstart = min(min(e[0] for e in i.exons) for i in isos)
    gend = max(max(e[1] for e in i.exons) for i in isos)
    chrom = isos[0].chromosome
    n_bins = max(1, (gend - gstart + bin_size) // bin_size)

    matrix = np.full((len(isos), n_bins), np.nan)
    labels: list[str] = []

    for i, iso in enumerate(isos):
        profile = profiles.get(iso.transcript_id)
        is_canon = "*" if iso.is_canonical else ""
        if profile is None or len(profile) != iso.mrna_length():
            labels.append(f"{iso.display_name}{is_canon} (no profile)")
            continue
        for (gs, _ge, chunk) in _mrna_to_genomic_chunks(profile, iso):
            for k in range(len(chunk)):
                gpos = gs + k
                b = (gpos - gstart) // bin_size
                if 0 <= b < n_bins:
                    current = matrix[i, b]
                    matrix[i, b] = chunk[k] if np.isnan(current) else max(current, chunk[k])
        labels.append(f"{iso.display_name}{is_canon}")

    fig, (ax_h, ax_rug) = plt.subplots(
        2, 1, figsize=(13, max(3.5, 0.4 * len(isos) + 1.2)),
        sharex=True, gridspec_kw={"height_ratios": [len(isos), 1]},
    )
    im = ax_h.imshow(
        matrix, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0,
        extent=[gstart, gend, len(isos) - 0.5, -0.5],
        interpolation="nearest",
    )
    ax_h.set_yticks(range(len(isos)))
    ax_h.set_yticklabels(labels, fontsize=8)
    if title:
        ax_h.set_title(title)
    else:
        ax_h.set_title(f"accessibility heatmap (T={temperature:g}°C, bin={bin_size}bp)")
    cbar = fig.colorbar(im, ax=ax_h, fraction=0.025, pad=0.01)
    cbar.set_label("p_unpaired", fontsize=8)

    _draw_probe_rug(ax_rug, probes, gstart, gend, chrom)

    fig.tight_layout()
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=130, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_rnafold_structure(
    seq: str,
    *,
    probe_windows: Optional[list[tuple[int, int]]] = None,
    window_size: Optional[int] = None,
    transcript_id: str = "",
    title_suffix: str = "",
    save_path: Optional[Path] = None,
) -> plt.Figure:
    """2D mRNA structure via ViennaRNA Python bindings.

    Args:
        seq: full transcript sequence (ACGU or ACGT).
        probe_windows: list of ``(start, end)`` half-open mRNA windows to
            highlight in red on the structure.
        window_size: if set, fold only a ``window_size``-nt slice centered
            on the first probe window (useful for >2 kb transcripts where
            the full fold is slow and visually crowded).
        transcript_id: appears in the title.
        title_suffix: extra title text.
    """
    _ensure_agg_backend()
    import RNA  # lazy — only this function needs it

    full = seq.upper().replace("T", "U")
    offset = 0
    if window_size and probe_windows:
        s, e = probe_windows[0]
        half = window_size // 2
        ws = max(0, (s + e) // 2 - half)
        we = min(len(full), ws + window_size)
        sub = full[ws:we]
        offset = ws
        adj_windows = [(s - offset, e - offset) for s, e in probe_windows
                        if 0 <= s - offset and e - offset <= len(sub)]
        seq_to_fold = sub
        probe_windows = adj_windows
    else:
        seq_to_fold = full

    structure, mfe = RNA.fold(seq_to_fold)
    coords = list(RNA.simple_xy_coordinates(structure))
    if len(coords) == len(seq_to_fold) + 1:
        coords = coords[1:]
    xs = [c.X for c in coords]
    ys = [c.Y for c in coords]

    fig, ax = plt.subplots(figsize=(11, 11))
    ax.plot(xs, ys, color="#444", linewidth=0.5, alpha=0.6, zorder=1)
    base_colors = {"A": "#1f77b4", "U": "#ff7f0e", "G": "#2ca02c",
                    "C": "#d62728", "T": "#ff7f0e", "N": "#888"}
    cols = [base_colors.get(b, "#888") for b in seq_to_fold]
    ax.scatter(xs, ys, c=cols, s=14, zorder=3, edgecolors="none")

    # Base pairs
    stack: list[int] = []
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            ax.plot([xs[j], xs[i]], [ys[j], ys[i]],
                     color="#888", linewidth=0.4, alpha=0.7, zorder=2)

    if probe_windows:
        for w_idx, (s, e) in enumerate(probe_windows):
            ax.plot(xs[s:e], ys[s:e], color="red", linewidth=2.5, alpha=0.85,
                     zorder=4,
                     label=f"probe {w_idx+1}: mRNA {s+offset}-{e+offset}")
        ax.legend(loc="upper right", fontsize=9)

    ax.set_aspect("equal")
    ax.axis("off")
    title = f"{transcript_id} "
    if offset:
        title += f"(window {offset}-{offset+len(seq_to_fold)}) "
    title += f"— RNA.fold MFE = {mfe:.2f} kcal/mol, length {len(seq_to_fold)} nt"
    if title_suffix:
        title += f"  |  {title_suffix}"
    ax.set_title(title, fontsize=10)

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=130, bbox_inches="tight")
        plt.close(fig)
    return fig


def _render_colored_dimer_ascii(
    ax, ascii_structure: str, *,
    highlight_a_indices: Optional[set] = None,
    highlight_b_indices: Optional[set] = None,
    green_a_indices: Optional[set] = None,
    green_b_indices: Optional[set] = None,
    yellow_a_indices: Optional[set] = None,
    yellow_b_indices: Optional[set] = None,
    x_left: float = 0.06, y_top: float = 0.72,
    char_w_ax: float = 0.0095, line_h_ax: float = 0.045,
) -> dict:
    """Render primer3's dimer ASCII character-by-character with optional
    foreground (red text) and background (green/yellow cell) highlighting.

    A is the top two ASCII rows (SEQ lines); B is the bottom two (STR
    lines). A reads 5'→3' left-to-right; B reads 3'→5' left-to-right
    (antiparallel). For each row we use the 'middle' line's character
    when present (= paired in dimer), otherwise the 'outer' line's (=
    unpaired bulge/tail). Empty columns are skipped on render.

    Args:
        highlight_{a,b}_indices: red+bold characters (junction nick).
        green_{a,b}_indices: green background cell (vicinity, contiguous OK).
        yellow_{a,b}_indices: yellow background cell (vicinity, broken — bulge / mismatch).

    All indices are in canonical 5'→3' coordinates for the respective strand.
    Returns the parsed dict from :func:`parse_primer3_dimer_pairing`.
    """
    from probe_designer.qc.dimer_ascii import parse_primer3_dimer_pairing

    parsed = parse_primer3_dimer_pairing(ascii_structure)
    a_idx_at_col = parsed["col_to_a_idx"]
    b_idx_at_col = parsed["col_to_b_idx"]

    raw_lines = ascii_structure.split("\n")
    contents: list[str] = []
    labels: list[str] = []
    for l in raw_lines:
        if "\t" not in l:
            continue
        prefix, content = l.split("\t", 1)
        labels.append("A:" if prefix == "SEQ" else "B:")
        contents.append(content)
    if len(contents) < 4:
        return parsed
    max_len = max(len(c) for c in contents)
    contents = [c.ljust(max_len) for c in contents]

    a_red = set(highlight_a_indices) if highlight_a_indices else set()
    b_red = set(highlight_b_indices) if highlight_b_indices else set()
    a_green = set(green_a_indices) if green_a_indices else set()
    b_green = set(green_b_indices) if green_b_indices else set()
    a_yellow = set(yellow_a_indices) if yellow_a_indices else set()
    b_yellow = set(yellow_b_indices) if yellow_b_indices else set()

    # Auto-scale char width so the whole ASCII fits within the axes width
    # (cap left at 0.94 of axes width). primer3 dimers for 40-mer ligator
    # vs 95-mer splint can easily be 120+ chars wide.
    available_width = 0.94
    effective_char_w = min(char_w_ax, available_width / max_len)

    # Container box
    box_w = effective_char_w * max_len + 0.04
    box_h = line_h_ax * 4 + 0.04
    box_x = x_left - 0.04
    box_y = y_top - line_h_ax * 3 - 0.025
    rect = plt.Rectangle(
        (box_x, box_y), box_w, box_h,
        transform=ax.transAxes,
        facecolor="#f7f7f0", edgecolor="#888", linewidth=0.6, zorder=1,
        clip_on=False,
    )
    ax.add_patch(rect)

    def _draw_bg(c: int, y: float, color: str) -> None:
        # A cell centered horizontally on the character glyph.
        bg = plt.Rectangle(
            (x_left + c * effective_char_w - effective_char_w * 0.1,
             y - line_h_ax * 0.45),
            effective_char_w * 1.2, line_h_ax * 0.9,
            transform=ax.transAxes,
            facecolor=color, alpha=0.45, edgecolor="none", zorder=2,
            clip_on=False,
        )
        ax.add_patch(bg)

    # We draw both A-line pairs (top + mid) with the SAME A-index background,
    # and both B-line pairs (mid + bot) with the SAME B-index background — so
    # the vicinity strip extends across both rows of each strand visually.
    a_line_ys: list[float] = []
    b_line_ys: list[float] = []
    for line_idx, label in enumerate(labels[:4]):
        y = y_top - line_idx * line_h_ax
        if label.startswith("A"):
            a_line_ys.append(y)
        else:
            b_line_ys.append(y)

    # Draw backgrounds: for each col, look up its A or B index and color it.
    for c in range(max_len):
        a_idx = a_idx_at_col.get(c)
        if a_idx is not None:
            for y in a_line_ys:
                if a_idx in a_green:
                    _draw_bg(c, y, "#7cd97c")
                elif a_idx in a_yellow:
                    _draw_bg(c, y, "#f7e370")
        b_idx = b_idx_at_col.get(c)
        if b_idx is not None:
            for y in b_line_ys:
                if b_idx in b_green:
                    _draw_bg(c, y, "#7cd97c")
                elif b_idx in b_yellow:
                    _draw_bg(c, y, "#f7e370")

    # Font size scales with effective char width to stay readable when squashed
    fontsize = max(6.5, min(10.0, effective_char_w / 0.0095 * 10.0))

    # Draw characters
    for line_idx, (label, content) in enumerate(zip(labels[:4], contents)):
        y = y_top - line_idx * line_h_ax
        ax.text(x_left - 0.022, y, label, family="monospace", fontsize=fontsize,
                 ha="left", va="center", transform=ax.transAxes, zorder=3)
        is_a_line = label.startswith("A")
        for c, ch in enumerate(content):
            if ch == " ":
                continue
            color = "black"
            weight = "normal"
            if is_a_line and c in a_idx_at_col and a_idx_at_col[c] in a_red:
                color = "red"
                weight = "bold"
            elif (not is_a_line) and c in b_idx_at_col and b_idx_at_col[c] in b_red:
                color = "red"
                weight = "bold"
            ax.text(
                x_left + c * effective_char_w, y, ch,
                family="monospace", fontsize=fontsize, color=color, weight=weight,
                ha="left", va="center", transform=ax.transAxes, zorder=4,
            )

    parsed["_render_max_len"] = max_len
    parsed["_render_x_left"] = x_left
    parsed["_render_y_top"] = y_top
    parsed["_render_char_w_ax"] = effective_char_w
    parsed["_render_line_h_ax"] = line_h_ax
    return parsed


def _compute_per_arm_vicinity_sets(
    arm_ascii: str, arm_len: int, end: str, n_each_side: int,
) -> tuple[set, set, set, set]:
    """v2.2 per-arm vicinity highlight: classify arm positions near the
    nick-side terminus into green (contiguous antiparallel run on B from
    the terminus) vs yellow (vicinity but broken by bulge / mismatch /
    non-adjacency). The terminal position itself (the 3'-OH on arm3 or
    5'-P on arm5) is excluded — rendered red.

    Args:
        arm_ascii: primer3 ASCII of (arm fragment vs splint) binary heterodimer.
        arm_len: full length of the arm fragment (for end="3" indexing).
        end: "3" for arm3 (terminal = last base, vicinity = last n+1);
             "5" for arm5 (terminal = first base, vicinity = first n+1).
        n_each_side: vicinity width (how many bases inward from the
            terminus to enforce contiguity, excluding the terminus itself).

    Returns ``(a_green, a_yellow, b_green, b_yellow)`` in canonical 5'→3'
    for the arm and the splint respectively.
    """
    from probe_designer.qc.dimer_ascii import parse_primer3_dimer_pairing
    parsed = parse_primer3_dimer_pairing(arm_ascii)
    pairs = parsed["pair_positions_a_to_b"]

    if end == "3":
        terminal = arm_len - 1
        # Vicinity walks inward (decreasing index) from the terminus.
        vicinity_range = list(range(terminal - n_each_side, terminal + 1))
        # Green starts from the terminus and extends inward as long as B-adjacent.
        a_green: set = set()
        b_green: set = set()
        if terminal in pairs:
            a = terminal
            while (a - 1) >= terminal - n_each_side and (a - 1) in pairs \
                    and pairs[a - 1] - pairs[a] == 1:
                a -= 1
                a_green.add(a)
                b_green.add(pairs[a])
    elif end == "5":
        terminal = 0
        vicinity_range = list(range(0, n_each_side + 1))
        a_green = set()
        b_green = set()
        if terminal in pairs:
            a = terminal
            while (a + 1) <= n_each_side and (a + 1) in pairs \
                    and pairs[a + 1] - pairs[a] == -1:
                a += 1
                a_green.add(a)
                b_green.add(pairs[a])
    else:
        raise ValueError(f"end must be '3' or '5', got {end!r}")

    a_yellow: set = set()
    b_yellow: set = set()
    for a in vicinity_range:
        if a == terminal or a in a_green:
            continue
        a_yellow.add(a)
        if a in pairs:
            b_yellow.add(pairs[a])

    return a_green, a_yellow, b_green, b_yellow


def _annotate_junction_markers(
    ax, parsed: dict, *,
    a_indices: Optional[set] = None,
    b_indices: Optional[set] = None,
) -> None:
    """Draw a small ▼ above A's junction columns and ▲ below B's junction columns."""
    if not parsed or "_render_x_left" not in parsed:
        return
    a_idx_to_col = parsed["a_idx_to_col"]
    b_idx_to_col = parsed["b_idx_to_col"]
    x_left = parsed["_render_x_left"]
    y_top = parsed["_render_y_top"]
    char_w_ax = parsed["_render_char_w_ax"]
    line_h_ax = parsed["_render_line_h_ax"]

    if a_indices:
        for idx in a_indices:
            c = a_idx_to_col.get(idx)
            if c is None:
                continue
            ax.text(x_left + c * char_w_ax, y_top + line_h_ax * 0.6, "▼",
                     family="monospace", fontsize=9, color="red", weight="bold",
                     ha="left", va="center", transform=ax.transAxes, zorder=4)
    if b_indices:
        for idx in b_indices:
            c = b_idx_to_col.get(idx)
            if c is None:
                continue
            ax.text(x_left + c * char_w_ax, y_top - line_h_ax * 3.5, "▲",
                     family="monospace", fontsize=9, color="red", weight="bold",
                     ha="left", va="center", transform=ax.transAxes, zorder=4)


def plot_cross_lig_dimer(
    probe_a: dict,
    probe_b: dict,
    *,
    save_path: Optional[Path] = None,
    title_extra: str = "",
    mv_conc_mm: float = 50.0,
    dv_conc_mm: float = 10.0,
    dna_conc_m: float = 1.0e-8,
    temp_c: float = 37.0,
    vicinity_n_each_side: int = 3,
) -> plt.Figure:
    """Render the v2.2 split-fragment primer3 ASCII for an A-as-ligator /
    B-as-splint cross-lig pair.

    Runs primer3 twice — (arm3 vs B.sequence) and (arm5 vs B.sequence) —
    and draws TWO stacked ASCII panels. arm3's 3'-OH (last base) + arm5's
    5'-P (first base) are red; per-arm vicinity ±n is green (contiguous)
    or yellow (broken).

    Args:
        probe_a, probe_b: dicts with ``probe_id``, ``arm5`` (or
            ``probe_arm5``), ``arm3`` (or ``probe_arm3``), ``chemistry``,
            ``sequence`` (mandatory for B), optional ``target``.
    """
    import primer3  # lazy

    from probe_designer.qc.cross_ligation import (
        ProbeForScreen, build_v2_geom,
        _check_nick_geometry, _check_per_arm_vicinity,
    )

    def _to_pfs(d: dict) -> ProbeForScreen:
        arm5 = d.get("probe_arm5") or d["arm5"]
        arm3 = d.get("probe_arm3") or d["arm3"]
        return ProbeForScreen(
            probe_id=d["probe_id"],
            chemistry=d.get("chemistry", "dRNA"),
            probe_arm5=arm5,
            probe_arm3=arm3,
            sequence=d["sequence"],
            target=d.get("target", ""),
        )

    pfs_a = _to_pfs(probe_a)
    pfs_b = _to_pfs(probe_b)

    arm3_eff, arm5, _ = build_v2_geom(pfs_a)
    _, _, splint_seq = build_v2_geom(pfs_b)

    buffer = dict(
        mv_conc=mv_conc_mm, dv_conc=dv_conc_mm,
        dna_conc=dna_conc_m * 1e9,
        temp_c=temp_c,
    )
    arm3_thermo = primer3.calc_heterodimer(
        arm3_eff, splint_seq, output_structure=True, **buffer,
    )
    arm5_thermo = primer3.calc_heterodimer(
        arm5, splint_seq, output_structure=True, **buffer,
    )
    arm3_ascii = arm3_thermo.ascii_structure or "" if arm3_thermo.structure_found else ""
    arm5_ascii = arm5_thermo.ascii_structure or "" if arm5_thermo.structure_found else ""

    can_ligate, b_3oh, b_5p = _check_nick_geometry(arm3_ascii, arm5_ascii, len(arm3_eff))
    vicinity_ok = (
        _check_per_arm_vicinity(
            arm3_ascii, arm5_ascii, len(arm3_eff), len(arm5), vicinity_n_each_side,
        )
        if (can_ligate and vicinity_n_each_side > 0) else False
    )

    arm3_tm = float(arm3_thermo.tm) if arm3_thermo.structure_found else 0.0
    arm5_tm = float(arm5_thermo.tm) if arm5_thermo.structure_found else 0.0
    arm3_dg = (float(arm3_thermo.dg) / 1000.0) if arm3_thermo.structure_found else 0.0
    arm5_dg = (float(arm5_thermo.dg) / 1000.0) if arm5_thermo.structure_found else 0.0

    return _draw_cross_lig_dimer_figure_v22(
        probe_a=probe_a, probe_b=probe_b,
        arm3_ascii=arm3_ascii, arm5_ascii=arm5_ascii,
        arm3_eff_len=len(arm3_eff), arm5_len=len(arm5),
        arm3_tm_c=arm3_tm, arm3_dg_kcal=arm3_dg,
        arm5_tm_c=arm5_tm, arm5_dg_kcal=arm5_dg,
        b_3oh_pos=b_3oh, b_5p_pos=b_5p,
        can_ligate=can_ligate,
        vicinity_contiguous=vicinity_ok,
        vicinity_n_each_side=vicinity_n_each_side,
        save_path=save_path,
        title_extra=title_extra,
    )


def plot_cross_lig_dimer_from_record(
    probe_a: dict,
    probe_b: dict,
    *,
    arm3_ascii: str,
    arm5_ascii: str,
    arm3_tm_c: float = 0.0,
    arm3_dg_kcal: float = 0.0,
    arm5_tm_c: float = 0.0,
    arm5_dg_kcal: float = 0.0,
    b_3oh_pos: Optional[int],
    b_5p_pos: Optional[int],
    can_ligate: bool,
    vicinity_contiguous: bool = False,
    vicinity_n_each_side: int = 0,
    save_path: Optional[Path] = None,
    title_extra: str = "",
) -> plt.Figure:
    """Render a v2.2 cross-lig figure from a stored audit record (no primer3 re-run).

    Use this when you already have ``arm3_dimer_structure``, ``arm5_dimer_structure``,
    Tm/dG, ``b_3oh_pos``, ``b_5p_pos`` from ``confirmed_dimers_v2.tsv``.
    """
    from probe_designer.qc.cross_ligation import ProbeForScreen, build_v2_geom

    def _to_pfs(d: dict) -> ProbeForScreen:
        arm5 = d.get("probe_arm5") or d["arm5"]
        arm3 = d.get("probe_arm3") or d["arm3"]
        return ProbeForScreen(
            probe_id=d["probe_id"],
            chemistry=d.get("chemistry", "dRNA"),
            probe_arm5=arm5,
            probe_arm3=arm3,
            sequence=d.get("sequence", "") or "PLACEHOLDER",
            target=d.get("target", ""),
        )

    arm3_eff, arm5, _ = build_v2_geom(_to_pfs(probe_a))

    return _draw_cross_lig_dimer_figure_v22(
        probe_a=probe_a, probe_b=probe_b,
        arm3_ascii=arm3_ascii, arm5_ascii=arm5_ascii,
        arm3_eff_len=len(arm3_eff), arm5_len=len(arm5),
        arm3_tm_c=arm3_tm_c, arm3_dg_kcal=arm3_dg_kcal,
        arm5_tm_c=arm5_tm_c, arm5_dg_kcal=arm5_dg_kcal,
        b_3oh_pos=b_3oh_pos, b_5p_pos=b_5p_pos,
        can_ligate=can_ligate,
        vicinity_contiguous=vicinity_contiguous,
        vicinity_n_each_side=vicinity_n_each_side,
        save_path=save_path,
        title_extra=title_extra,
    )


def _draw_cross_lig_dimer_figure_v22(
    *,
    probe_a: dict, probe_b: dict,
    arm3_ascii: str, arm5_ascii: str,
    arm3_eff_len: int, arm5_len: int,
    arm3_tm_c: float, arm3_dg_kcal: float,
    arm5_tm_c: float, arm5_dg_kcal: float,
    b_3oh_pos: Optional[int], b_5p_pos: Optional[int],
    can_ligate: bool,
    vicinity_contiguous: bool,
    vicinity_n_each_side: int,
    save_path: Optional[Path],
    title_extra: str,
) -> plt.Figure:
    """Render two stacked ASCII panels (arm3 + arm5 vs splint) with verdict box."""
    fig = plt.figure(figsize=(15, 9))
    ax = fig.add_subplot(111)
    ax.axis("off")

    target_a = probe_a.get("target", "")
    target_b = probe_b.get("target", "")
    chem_a = probe_a.get("chemistry", "?")
    chem_b = probe_b.get("chemistry", "?")
    title = (
        f"A=ligator: {probe_a['probe_id']} ({target_a}) [{chem_a}]    "
        f"B=splint:  {probe_b['probe_id']} ({target_b}) [{chem_b}]    [v2.2 split-fragment]"
    )
    if title_extra:
        title += f"   — {title_extra}"
    ax.text(0.5, 0.985, title, ha="center", va="top", fontsize=11.5,
             weight="bold", transform=ax.transAxes)

    arm5_a = probe_a.get("probe_arm5") or probe_a.get("arm5", "")
    arm3_a = probe_a.get("probe_arm3") or probe_a.get("arm3", "")
    arm5_b = probe_b.get("probe_arm5") or probe_b.get("arm5", "")
    arm3_b = probe_b.get("probe_arm3") or probe_b.get("arm3", "")
    ax.text(
        0.02, 0.95,
        f"A fragments (independent in solution before ligation):  "
        f"arm3=[{arm3_a}]  arm5=[{arm5_a}]",
        ha="left", va="top", family="monospace", fontsize=8.0, transform=ax.transAxes,
    )
    ax.text(
        0.02, 0.92,
        f"B splint (full assembled probe):  5'-[arm5: {arm5_b}]-[bb..]-[arm3: {arm3_b}]-3'",
        ha="left", va="top", family="monospace", fontsize=8.0, transform=ax.transAxes,
    )

    # ---- ARM3 panel (top) ----
    ax.text(
        0.02, 0.86,
        f"ARM3 vs SPLINT — arm3's 3'-OH = position {arm3_eff_len - 1} (red).  "
        f"Tm = {arm3_tm_c:.1f} °C    ΔG = {arm3_dg_kcal:.2f} kcal/mol",
        ha="left", va="top", family="monospace", fontsize=9.0, weight="bold",
        transform=ax.transAxes,
    )
    arm3_terminal = arm3_eff_len - 1
    arm3_a_hi = {arm3_terminal}
    arm3_b_hi: set = set()
    if b_3oh_pos is not None:
        arm3_b_hi = {b_3oh_pos}
    arm3_a_g: set = set(); arm3_a_y: set = set()
    arm3_b_g: set = set(); arm3_b_y: set = set()
    if arm3_ascii and vicinity_n_each_side > 0:
        arm3_a_g, arm3_a_y, arm3_b_g, arm3_b_y = _compute_per_arm_vicinity_sets(
            arm3_ascii, arm3_eff_len, end="3", n_each_side=vicinity_n_each_side,
        )

    if arm3_ascii:
        _render_colored_dimer_ascii(
            ax, arm3_ascii,
            highlight_a_indices=arm3_a_hi, highlight_b_indices=arm3_b_hi,
            green_a_indices=arm3_a_g, green_b_indices=arm3_b_g,
            yellow_a_indices=arm3_a_y, yellow_b_indices=arm3_b_y,
            x_left=0.06, y_top=0.81,
        )
    else:
        ax.text(0.5, 0.79, "(primer3: no stable arm3 dimer found)",
                 ha="center", va="center", fontsize=10, style="italic",
                 color="gray", transform=ax.transAxes)

    # ---- ARM5 panel (bottom) ----
    ax.text(
        0.02, 0.54,
        f"ARM5 vs SPLINT — arm5's 5'-P = position 0 (red).  "
        f"Tm = {arm5_tm_c:.1f} °C    ΔG = {arm5_dg_kcal:.2f} kcal/mol",
        ha="left", va="top", family="monospace", fontsize=9.0, weight="bold",
        transform=ax.transAxes,
    )
    arm5_a_hi = {0}
    arm5_b_hi: set = set()
    if b_5p_pos is not None:
        arm5_b_hi = {b_5p_pos}
    arm5_a_g: set = set(); arm5_a_y: set = set()
    arm5_b_g: set = set(); arm5_b_y: set = set()
    if arm5_ascii and vicinity_n_each_side > 0:
        arm5_a_g, arm5_a_y, arm5_b_g, arm5_b_y = _compute_per_arm_vicinity_sets(
            arm5_ascii, arm5_len, end="5", n_each_side=vicinity_n_each_side,
        )

    if arm5_ascii:
        _render_colored_dimer_ascii(
            ax, arm5_ascii,
            highlight_a_indices=arm5_a_hi, highlight_b_indices=arm5_b_hi,
            green_a_indices=arm5_a_g, green_b_indices=arm5_b_g,
            yellow_a_indices=arm5_a_y, yellow_b_indices=arm5_b_y,
            x_left=0.06, y_top=0.49,
        )
    else:
        ax.text(0.5, 0.47, "(primer3: no stable arm5 dimer found)",
                 ha="center", va="center", fontsize=10, style="italic",
                 color="gray", transform=ax.transAxes)

    # ---- Verdict box ----
    verdict_parts: list[str] = []
    verdict_parts.append(
        f"NICK ON SPLINT (B):  arm3 3'-OH → B pos {b_3oh_pos}   |   "
        f"arm5 5'-P → B pos {b_5p_pos}"
    )
    if b_3oh_pos is not None and b_5p_pos is not None:
        delta = b_3oh_pos - b_5p_pos
        if delta == 1:
            verdict_parts.append(
                f"NICK ADJACENCY:  Δ(b_3oh − b_5p) = +1  →  nick-adjacent on B  ✓"
            )
        else:
            verdict_parts.append(
                f"NICK ADJACENCY:  Δ(b_3oh − b_5p) = {delta:+d}  →  NOT nick-adjacent on B  ✗"
            )
    else:
        verdict_parts.append("NICK ADJACENCY:  arm3 3'-OH or arm5 5'-P NOT paired in primer3 MFE  ✗")

    if vicinity_n_each_side > 0:
        if vicinity_contiguous:
            verdict_parts.append(
                f"VICINITY ±{vicinity_n_each_side} (green per-arm):  arm3 last {vicinity_n_each_side + 1} "
                f"AND arm5 first {vicinity_n_each_side + 1} bases each contiguously paired on B  ✓"
            )
        else:
            verdict_parts.append(
                f"VICINITY ±{vicinity_n_each_side} (yellow per-arm):  bulge / mismatch / non-adjacent within "
                f"{vicinity_n_each_side + 1} bases of one or both arm termini  ✗"
            )
    is_risk = can_ligate and (vicinity_contiguous if vicinity_n_each_side > 0 else True)
    verdict_parts.append(
        "VERDICT:  LIGATION-COMPETENT for SplintR (cross-lig risk)" if is_risk
        else "VERDICT:  NOT ligation-competent"
    )

    ax.text(0.02, 0.18, "\n".join(verdict_parts), ha="left", va="top",
             family="monospace", fontsize=9.0, transform=ax.transAxes,
             bbox=dict(boxstyle="round,pad=0.4",
                        fc="#ffe5e5" if is_risk else "#e5ffe5", ec="#bb6666"))

    fig.tight_layout()
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=140, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_ternary_complex(
    probe_a: dict,
    probe_b: dict,
    *,
    mfe_dotbracket: str,
    ternary_dg_kcal: float,
    ternary_concentration_m: float,
    ternary_fraction_of_b: float,
    ensemble_used: str,
    celsius_used: float,
    sodium_m: float,
    magnesium_m: float,
    strand_conc_m: float,
    mfe_nick_adjacent: bool,
    mfe_vicinity_contiguous: bool,
    b_3oh_pos: Optional[int],
    b_5p_pos: Optional[int],
    vicinity_n_each_side: int,
    save_path: Optional[Path] = None,
    title_extra: str = "",
) -> plt.Figure:
    """Render the NUPACK 4 ternary-complex visualization (v2.3).

    Layout:

    * Top: probe metadata + NUPACK Model parameters used.
    * Middle: a wrapped pseudo-alignment ASCII showing splint 5'→3' as the
      central strand, with arm3 and arm5 in their MFE-paired positions
      (extracted from the 3-strand dot-bracket). Arm termini in red,
      vicinity in green (contiguous) or yellow (broken).
    * Bottom: verdict box reporting ternary ΔG, equilibrium concentration,
      fraction-of-splint metric, and the MFE-geometry flags.

    The dot-bracket is assumed to follow NUPACK's complex-name strand order
    — for cross-lig audits this is invariably ``splint+arm5+arm3``. The
    function parses on this assumption; pass a different ordering at your
    peril (the verdict text dump always renders correctly even if the
    pseudo-alignment layout doesn't).
    """
    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_subplot(111)
    ax.axis("off")

    target_a = probe_a.get("target", "")
    target_b = probe_b.get("target", "")
    chem_a = probe_a.get("chemistry", "?")
    chem_b = probe_b.get("chemistry", "?")
    title = (
        f"A=ligator: {probe_a['probe_id']} ({target_a}) [{chem_a}]    "
        f"B=splint:  {probe_b['probe_id']} ({target_b}) [{chem_b}]    "
        f"[v2.3 NUPACK ternary, ensemble={ensemble_used}]"
    )
    if title_extra:
        title += f"   — {title_extra}"
    ax.text(0.5, 0.985, title, ha="center", va="top", fontsize=11.5,
             weight="bold", transform=ax.transAxes)

    arm5_a = probe_a.get("probe_arm5") or probe_a.get("arm5", "")
    arm3_a = probe_a.get("probe_arm3") or probe_a.get("arm3", "")
    arm5_b = probe_b.get("probe_arm5") or probe_b.get("arm5", "")
    arm3_b = probe_b.get("probe_arm3") or probe_b.get("arm3", "")
    ax.text(
        0.02, 0.95,
        f"A fragments:   arm3=[{arm3_a}]   arm5=[{arm5_a}]",
        ha="left", va="top", family="monospace", fontsize=8.0, transform=ax.transAxes,
    )
    ax.text(
        0.02, 0.92,
        f"B splint:      5'-[arm5: {arm5_b}]-[bb..]-[arm3: {arm3_b}]-3'",
        ha="left", va="top", family="monospace", fontsize=8.0, transform=ax.transAxes,
    )
    ax.text(
        0.02, 0.89,
        f"NUPACK Model:  celsius={celsius_used:.1f} °C  Na+={sodium_m:.3f} M  Mg2+={magnesium_m:.4f} M  "
        f"[strand]={strand_conc_m:.1e} M",
        ha="left", va="top", family="monospace", fontsize=8.0, transform=ax.transAxes,
    )

    # --- MFE dot-bracket display ---
    # Strands appear in NUPACK complex order: splint | arm5 | arm3
    # We split on '+' and label each chunk.
    chunks = mfe_dotbracket.split("+") if mfe_dotbracket else []
    chunk_labels = ["splint", "arm5", "arm3"]
    if len(chunks) == 3:
        ax.text(
            0.02, 0.83,
            "NUPACK MFE structure (dot-bracket, strand order = splint | arm5 | arm3):",
            ha="left", va="top", family="monospace", fontsize=9.0, weight="bold",
            transform=ax.transAxes,
        )
        y = 0.79
        for lbl, ch in zip(chunk_labels, chunks):
            # Wrap long splint lines (~120 chars)
            ax.text(
                0.04, y, f"{lbl:>8s}:  {ch[:120]}",
                ha="left", va="top", family="monospace", fontsize=7.5,
                transform=ax.transAxes,
            )
            y -= 0.035
            if len(ch) > 120:
                ax.text(
                    0.04, y, f"          {ch[120:240]}",
                    ha="left", va="top", family="monospace", fontsize=7.5,
                    transform=ax.transAxes,
                )
                y -= 0.035
    else:
        ax.text(0.5, 0.78, "(MFE dot-bracket: ternary not enumerated)",
                 ha="center", va="center", fontsize=10, style="italic",
                 color="gray", transform=ax.transAxes)

    # --- Verdict box ---
    verdict_parts: list[str] = []
    verdict_parts.append(
        f"Ternary ΔG:        {ternary_dg_kcal:+.2f} kcal/mol  (includes coaxial stacking when ensemble=stacking)"
    )
    verdict_parts.append(
        f"Ternary [conc]:    {ternary_concentration_m:.3e} M     "
        f"fraction-of-splint: {ternary_fraction_of_b:.3e}"
    )
    if b_3oh_pos is not None and b_5p_pos is not None:
        delta = b_3oh_pos - b_5p_pos
        verdict_parts.append(
            f"MFE nick:          arm3 3'-OH ↔ splint pos {b_3oh_pos}   "
            f"arm5 5'-P ↔ splint pos {b_5p_pos}   Δ={delta:+d}"
        )
    else:
        verdict_parts.append(
            "MFE nick:          arm3 3'-OH or arm5 5'-P NOT paired to splint in MFE structure"
        )

    nick_ok = mfe_nick_adjacent
    vicinity_ok = mfe_vicinity_contiguous
    is_risk = nick_ok and vicinity_ok and (ternary_fraction_of_b > 1e-4)

    verdict_parts.append(
        f"MFE nick-adjacent (Δ=+1):     {'✓' if nick_ok else '✗'}"
    )
    verdict_parts.append(
        f"MFE per-arm vicinity ±{vicinity_n_each_side}:        {'✓' if vicinity_ok else '✗'}"
    )
    verdict_parts.append(
        f"Cross-lig risk (fraction > 1e-4): {'✓ HIGH RISK' if ternary_fraction_of_b > 1e-4 else '✗ low'}"
    )
    verdict_parts.append("")
    if is_risk:
        verdict_parts.append(
            "VERDICT:  NUPACK-CONFIRMED CROSS-LIG RISK (coaxial-stacking-aware ensemble)"
        )
    else:
        reasons = []
        if not nick_ok:
            reasons.append("MFE not nick-adjacent")
        if not vicinity_ok:
            reasons.append("vicinity broken")
        if ternary_fraction_of_b <= 1e-4:
            reasons.append(f"fraction too low ({ternary_fraction_of_b:.1e})")
        verdict_parts.append(
            "VERDICT:  not ligation-competent in ensemble  (" + ", ".join(reasons) + ")"
        )

    ax.text(0.02, 0.32, "\n".join(verdict_parts), ha="left", va="top",
             family="monospace", fontsize=9.0, transform=ax.transAxes,
             bbox=dict(boxstyle="round,pad=0.5",
                        fc="#ffe5e5" if is_risk else "#e5ffe5", ec="#bb6666"))

    fig.tight_layout()
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=140, bbox_inches="tight")
        plt.close(fig)
    return fig


__all__ = [
    "IsoformMeta",
    "ProbeAnnotation",
    "plot_accessibility_tracks",
    "plot_accessibility_heatmap",
    "plot_rnafold_structure",
    "plot_cross_lig_dimer",
    "plot_cross_lig_dimer_from_record",
    "plot_ternary_complex",
]
