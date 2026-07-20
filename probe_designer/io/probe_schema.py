"""Unified probe-output schema for mRNA / mutation / TCR pipelines.

Every pipeline that produces a final probe Excel file must funnel its
DataFrame through ``apply_final_column_order`` so the user-facing files
share one column contract:

  1. ``order``, ``probe_name``, ``probe_sequence`` (always first three)
  2. identity, chemistry, target, probe parts
  3. encoding (``No.``, ``codebook``, backbone metadata)
  4. thermodynamics, quality, BLAST
  5. pipeline-specific tail (mutation: genomic/mutation context;
     TCR: cDNA-Tm warning)

``probe_name`` follows ``{name}_{position}_{chemistry}_{codebook}_{No.}``
for padlocks and ``{name}_{position}_RT_primer`` for RT primers — built
via ``make_padlock_name`` / ``make_rt_primer_name`` so callers don't
hand-roll f-strings.

The migration CLI (``probe-design migrate-columns``) and the
``probe-order`` skill consume ``LEGACY_RENAMES`` and the parsers below
to read pre-Schema-v2 files.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Literal, Optional

import pandas as pd


# ---------------------------------------------------------------------------
# chemistry tags
# ---------------------------------------------------------------------------

CHEM_DRNA = "dRNA"
CHEM_CDNA = "cDNA"
CHEM_ILOCK = "iLock"
CHEM_RT = "RT_primer"

PADLOCK_CHEMISTRIES: set[str] = {CHEM_DRNA, CHEM_CDNA, CHEM_ILOCK}
ALL_CHEMISTRIES: set[str] = PADLOCK_CHEMISTRIES | {CHEM_RT}


# ---------------------------------------------------------------------------
# canonical column orders
# ---------------------------------------------------------------------------

FINAL_PADLOCK_COLUMNS: list[str] = [
    # 1. user-mandated leading three
    "order", "probe_name", "probe_sequence",
    # 2. identity
    "gene_name", "clone_id",
    # 3. chemistry (iLock is a value, not a separate column)
    "chemistry",
    # 4. target context
    "position", "target_sequence",
    # 5. probe parts
    "probe_arm5", "probe_arm3", "probe_length",
    # 6. encoding
    "No.", "codebook", "backbone_name", "backbone_sequence", "backbone_file",
    # 7. thermodynamics (Schema-v2 keeps both for back-compat:
    #    gc_content = G+C fraction (the primary gate, set by Phase 0 thermal_filter);
    #    g_content  = G-only fraction (legacy; informational only).
    "gc_content", "g_content", "tm", "tm_5prime", "tm_3prime",
    # 8. quality
    "score", "free_energy",
    # 9. BLAST
    "blast_hit_count", "blast_top_hits", "same_gene_exact",
    # 10. mutation-only tail
    "chr", "genomic_start", "genomic_end", "ref", "alt", "strand",
    "mutation_type", "tm_model", "mutation_position",
    "arm5_len", "arm3_len", "plp_pos5_local",
    "ilock_flap", "ilock_linker_nt",
    # 11. TCR-only tail
    "tm_cdna_warning",
]


FINAL_RT_PRIMER_COLUMNS: list[str] = [
    # leading three
    "order", "probe_name", "probe_sequence",
    # identity + chemistry tag (always RT_primer)
    "gene_name", "clone_id", "chemistry",
    # target
    "position", "target_sequence",
    # primer body
    "probe_length", "gc_content", "tm",
    # link back to the padlock this primer pairs with
    "paired_padlock_probe_name", "paired_padlock_no",
    # design notes
    "notes",
    # mutation extras
    "chr", "genomic_start", "genomic_end", "ref", "alt", "strand",
    "mutation_type",
    "primer_genomic_start", "primer_genomic_end", "gap_nt",
    # TCR extras
    "primer_start_on_target", "primer_end_on_target",
]


# ---------------------------------------------------------------------------
# legacy → canonical rename map (used by migration CLI + skill fallback)
# ---------------------------------------------------------------------------

LEGACY_RENAMES: dict[str, str] = {
    # name / sequence / length
    "Probe_Name": "probe_name",
    "RT_Primer_Name": "probe_name",
    "Probe_Seq": "probe_sequence",
    "probe_seq": "probe_sequence",
    "RT_Primer_Sequence": "probe_sequence",
    "Final_Probe_Sequence": "probe_sequence",
    "Probe_Length": "probe_length",
    "RT_Primer_Length": "probe_length",
    # arms
    "arm_5prime": "probe_arm5",
    "arm_3prime": "probe_arm3",
    # GC + Tm
    # Schema-v2 keeps g_content and gc_content as distinct columns (one counts
    # G only, the other G+C). DO NOT rename g_content → gc_content; the values
    # are not interchangeable. Legacy `GC_Percent` columns, in contrast, are
    # unambiguously G+C and can be safely mapped to gc_content.
    "GC_Percent": "gc_content",
    "Tm": "tm",
    # BLAST
    "blast_hits_count": "blast_hit_count",
    # genomic context (mutation pipeline)
    "Chr": "chr",
    "Start": "genomic_start",
    "End": "genomic_end",
    "Ref": "ref",
    "Alt": "alt",
    "Strand": "strand",
    "Mutation_Type": "mutation_type",
    "tm_cDNA_warning": "tm_cdna_warning",
    # mutation RT primer extras
    "Gene": "gene_name",
    "Mut_Start": "genomic_start",
    "Mut_End": "genomic_end",
    "Gap_nt": "gap_nt",
    "Primer_Genomic_Start": "primer_genomic_start",
    "Primer_Genomic_End": "primer_genomic_end",
    "Notes": "notes",
    # TCR RT primer extras
    "BDS_pos": "position",
    "Primer_Start_on_Target": "primer_start_on_target",
    "Primer_End_on_Target": "primer_end_on_target",
    # legacy iLock-flap column casing
    "iLock_flap": "ilock_flap",
    "iLock_linker_nt": "ilock_linker_nt",
    # legacy mRNA position / end naming
    "st": "position",
    # columns that no longer exist in v2 — flagged as "drop"
    "backbone_No.": "drop",
    "Padlock_Probe_No.": "paired_padlock_no",
    "Padlock_Probe_Name": "paired_padlock_probe_name",
    "BDS_end": "drop",
    "target_type": "drop",
    "iLock": "drop",   # boolean column now folded into `chemistry`
    "modification": "drop",  # mRNA pipeline's old marker for iLock
    "en": "drop",      # paired with `st` -> position; redundant
    # energy column (TCR called it mfe; align with mutation's free_energy)
    "mfe": "free_energy",
    # mRNA pipeline backbone-No. column was misnamed
    "probe_id": "No.",
}


# Chemistry value renames (apply after column renames so the column name
# matches before we operate on values).
LEGACY_CHEMISTRY_VALUES: dict[str, str] = {
    "mRNA_noiLock": CHEM_DRNA,   # mutation pipeline's old name for direct-RNA
    "mRNA": CHEM_DRNA,           # TCR pipeline's old name for direct-RNA
}


# ---------------------------------------------------------------------------
# probe_name builders
# ---------------------------------------------------------------------------


def make_padlock_name(
    *,
    name: str,
    position: int,
    chemistry: str,
    codebook: str,
    no: int,
) -> str:
    """Format `{name}_{position}_{chemistry}_{codebook}_{No.}`.

    `chemistry` must be one of ``PADLOCK_CHEMISTRIES`` (RT primers go
    through :func:`make_rt_primer_name`).
    """
    if chemistry not in PADLOCK_CHEMISTRIES:
        raise ValueError(
            f"unknown padlock chemistry {chemistry!r}; "
            f"allowed: {sorted(PADLOCK_CHEMISTRIES)}"
        )
    return f"{name}_{int(position)}_{chemistry}_{codebook}_{int(no)}"


def make_rt_primer_name(*, name: str, position: int) -> str:
    """Format `{name}_{position}_RT_primer`.

    RT primers don't carry a codebook + No. (they're not part of any
    barcode panel) so the suffix is the literal ``RT_primer``.
    """
    return f"{name}_{int(position)}_RT_primer"


# ---------------------------------------------------------------------------
# probe_name parsers (round-trip helpers for migration + probe-order skill)
# ---------------------------------------------------------------------------


_PADLOCK_RE = re.compile(
    r"^(?P<name>.+?)"
    r"_(?P<position>\d+)"
    r"_(?P<chemistry>dRNA|cDNA|iLock)"
    r"_(?P<codebook>[A-Za-z]+\d+)"
    r"_(?P<no>\d+)$"
)


_RT_RE = re.compile(r"^(?P<name>.+?)_(?P<position>\d+)_RT_primer$")


def parse_padlock_name(name: str) -> Optional[dict]:
    """Parse a Schema-v2 padlock probe name; return None on no match.

    The name may itself contain underscores (e.g. ``BZ23_clone23_TRB``);
    the regex anchors at the tail so the leftover prefix is the name.
    """
    m = _PADLOCK_RE.match(name)
    if m is None:
        return None
    return {
        "name": m.group("name"),
        "position": int(m.group("position")),
        "chemistry": m.group("chemistry"),
        "codebook": m.group("codebook"),
        "no": int(m.group("no")),
    }


def parse_rt_primer_name(name: str) -> Optional[dict]:
    """Parse a Schema-v2 RT primer name; return None on no match."""
    m = _RT_RE.match(name)
    if m is None:
        return None
    return {
        "name": m.group("name"),
        "position": int(m.group("position")),
    }


# ---------------------------------------------------------------------------
# codebook resolution
# ---------------------------------------------------------------------------


_CODEBOOK_FROM_FILENAME = re.compile(r"([A-Za-z]+\d+)")


def resolve_codebook(cli_value: Optional[str], backbone_file: Path) -> str:
    """Return the codebook name (e.g. ``SP369``).

    Resolution order:
      1. ``cli_value`` if not None / not empty (CLI flag or YAML config).
      2. First ``[A-Za-z]+\\d+`` token found in ``backbone_file.name``
         (matches ``SP369``, ``MIT30``, ``PRISM2``, etc.).

    Raises:
        ValueError: when no value is provided and the filename has no
            recognisable codebook token.
    """
    if cli_value:
        return cli_value
    m = _CODEBOOK_FROM_FILENAME.search(backbone_file.name)
    if m is not None:
        return m.group(1)
    raise ValueError(
        f"could not infer codebook from backbone filename {backbone_file.name!r}; "
        f"pass --codebook explicitly (e.g. SP369, MIT30)."
    )


# ---------------------------------------------------------------------------
# column-order applier
# ---------------------------------------------------------------------------


_KIND_TO_SCHEMA = {
    "padlock": FINAL_PADLOCK_COLUMNS,
    "rt": FINAL_RT_PRIMER_COLUMNS,
}


def apply_final_column_order(
    df: pd.DataFrame, *, kind: Literal["padlock", "rt"],
) -> pd.DataFrame:
    """Reorder ``df`` columns to match the canonical schema for ``kind``.

    Behaviour:
      * Schema columns come first, in the order defined by
        ``FINAL_PADLOCK_COLUMNS`` / ``FINAL_RT_PRIMER_COLUMNS``.
      * Schema columns missing from ``df`` are added with NaN values so
        that pipelines can ignore irrelevant fields without breaking
        downstream concat.
      * Any caller-provided columns NOT in the schema are appended at
        the tail in their original order — never silently dropped.
    """
    if kind not in _KIND_TO_SCHEMA:
        raise ValueError(
            f"unknown kind {kind!r}; expected one of {list(_KIND_TO_SCHEMA)}"
        )
    schema_cols = _KIND_TO_SCHEMA[kind]
    out = df.copy()
    for col in schema_cols:
        if col not in out.columns:
            out[col] = pd.NA
    extras = [c for c in df.columns if c not in schema_cols]
    return out[schema_cols + extras]


__all__ = [
    "CHEM_DRNA", "CHEM_CDNA", "CHEM_ILOCK", "CHEM_RT",
    "PADLOCK_CHEMISTRIES", "ALL_CHEMISTRIES",
    "FINAL_PADLOCK_COLUMNS", "FINAL_RT_PRIMER_COLUMNS",
    "LEGACY_RENAMES", "LEGACY_CHEMISTRY_VALUES",
    "make_padlock_name", "make_rt_primer_name",
    "parse_padlock_name", "parse_rt_primer_name",
    "resolve_codebook", "apply_final_column_order",
]
