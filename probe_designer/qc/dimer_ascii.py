"""Parse primer3's 4-line dimer ASCII into structured base-pair info.

primer3-py's ``calc_heterodimer(..., output_structure=True)`` returns a
``ThermoResult.ascii_structure`` string with this exact format:

    SEQ\\t<seq_a top line — A's UNPAIRED bases at their columns>
    SEQ\\t<seq_a mid line — A's PAIRED bases at their columns>
    STR\\t<seq_b mid line — B's PAIRED bases at their columns>
    STR\\t<seq_b bot line — B's UNPAIRED bases at their columns>

A reads 5'→3' left-to-right (across both SEQ lines combined). B reads
**3'→5' left-to-right** in the display (= antiparallel), which means
B's 5'→3' index = ``b_len - 1 - display_index``. Per column, mid-line takes
precedence over top/bot; a column is "paired" iff BOTH mid-A and mid-B
have a base (non-space, non-dash) at that column.

This module's sole responsibility: convert the ASCII into structured maps
in **canonical 5'→3' coordinates** for both strands. v2 cross-lig consumes
the result to check the ligation-junction geometry (whether A's real
3'-OH and 5'-P pair on adjacent positions of B in B's 5'→3').

See ``[[project-pool-aware-design]]`` and the Phase 2D.2 section of the
plan for the broader v2 architecture.
"""
from __future__ import annotations

from typing import Any, Dict


def parse_primer3_dimer_pairing(ascii_structure: str) -> Dict[str, Any]:
    """Parse primer3 4-line dimer ASCII; return base-pair info.

    Args:
        ascii_structure: the raw string from
            ``primer3.calc_heterodimer(..., output_structure=True).ascii_structure``.
            Format: 4 newline-separated lines, each prefixed with ``SEQ\\t`` or
            ``STR\\t``. Garbage / empty input yields an empty result rather
            than raising.

    Returns dict with keys:
        - ``a_len`` (int): seq_a length recovered from the ASCII.
        - ``b_len`` (int): seq_b length recovered from the ASCII.
        - ``col_to_a_idx`` (dict[int, int]): ASCII column → A's 5'→3' index.
          Only columns where A has a base (in top OR mid line) are present.
        - ``col_to_b_idx`` (dict[int, int]): ASCII column → B's 5'→3' index
          (canonical, NOT the display direction).
        - ``a_idx_to_col``, ``b_idx_to_col``: inverse lookups.
        - ``paired_columns`` (set[int]): columns where mid-A and mid-B both
          have a base — these are the duplex base pairs.
        - ``pair_positions_a_to_b`` (dict[int, int]): for each paired
          column, A's 5'→3' index → B's 5'→3' index. **The ligation criterion
          checks specific keys of this dict** (e.g., the rotated arm3/arm5
          junction positions).
    """
    empty_result: Dict[str, Any] = {
        "a_len": 0, "b_len": 0,
        "col_to_a_idx": {}, "col_to_b_idx": {},
        "a_idx_to_col": {}, "b_idx_to_col": {},
        "paired_columns": set(),
        "pair_positions_a_to_b": {},
    }

    if not isinstance(ascii_structure, str) or not ascii_structure:
        return empty_result

    # Split and keep only the 4 SEQ/STR prefixed lines
    contents: list[str] = []
    labels: list[str] = []
    for line in ascii_structure.split("\n"):
        if "\t" not in line:
            continue
        prefix, content = line.split("\t", 1)
        labels.append(prefix)
        contents.append(content)
    if len(contents) < 4:
        return empty_result

    # Order: SEQ_top, SEQ_mid, STR_mid, STR_bot
    a_top, a_mid, b_mid, b_bot = contents[0], contents[1], contents[2], contents[3]
    max_len = max(len(a_top), len(a_mid), len(b_mid), len(b_bot))
    a_top = a_top.ljust(max_len)
    a_mid = a_mid.ljust(max_len)
    b_mid = b_mid.ljust(max_len)
    b_bot = b_bot.ljust(max_len)

    def _is_base(ch: str) -> bool:
        return ch not in (" ", "-", "")

    # Walk columns for A — mid takes precedence over top
    a_idx = 0
    col_to_a_idx: Dict[int, int] = {}
    a_paired_cols: set[int] = set()
    for c in range(max_len):
        if _is_base(a_mid[c]):
            col_to_a_idx[c] = a_idx
            a_paired_cols.add(c)
            a_idx += 1
        elif _is_base(a_top[c]):
            col_to_a_idx[c] = a_idx
            a_idx += 1
    a_len = a_idx

    # Walk columns for B — track display index, will convert to canonical 5'→3'
    b_display_idx = 0
    col_to_b_display: Dict[int, int] = {}
    b_paired_cols: set[int] = set()
    for c in range(max_len):
        if _is_base(b_mid[c]):
            col_to_b_display[c] = b_display_idx
            b_paired_cols.add(c)
            b_display_idx += 1
        elif _is_base(b_bot[c]):
            col_to_b_display[c] = b_display_idx
            b_display_idx += 1
    b_len = b_display_idx

    # Convert B display index → canonical 5'→3' index.
    # Display direction is 3'→5' left-to-right, so display 0 = canonical (b_len - 1).
    col_to_b_idx: Dict[int, int] = {
        col: (b_len - 1 - di) for col, di in col_to_b_display.items()
    }

    paired_columns = a_paired_cols & b_paired_cols
    pair_positions_a_to_b: Dict[int, int] = {
        col_to_a_idx[c]: col_to_b_idx[c] for c in paired_columns
    }

    a_idx_to_col = {v: k for k, v in col_to_a_idx.items()}
    b_idx_to_col = {v: k for k, v in col_to_b_idx.items()}

    return {
        "a_len": a_len, "b_len": b_len,
        "col_to_a_idx": col_to_a_idx,
        "col_to_b_idx": col_to_b_idx,
        "a_idx_to_col": a_idx_to_col,
        "b_idx_to_col": b_idx_to_col,
        "paired_columns": paired_columns,
        "pair_positions_a_to_b": pair_positions_a_to_b,
    }


__all__ = ["parse_primer3_dimer_pairing"]
