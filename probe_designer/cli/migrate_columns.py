"""``probe-design migrate-columns`` — upgrade legacy probe Excel files to Schema v2.

Walks an experiment directory (or operates on a single file), renames
columns per :data:`LEGACY_RENAMES`, drops the standalone ``iLock`` column
(folding its information into ``chemistry``), normalizes legacy chemistry
tokens (``mRNA_noiLock``/``mRNA`` -> ``dRNA``), regenerates
``probe_name`` against the user-supplied codebook, and re-orders columns
to match :data:`FINAL_PADLOCK_COLUMNS` / :data:`FINAL_RT_PRIMER_COLUMNS`.

A ``.bak`` file is written alongside each migrated input so the operation
is reversible. The CLI never writes to a file twice — re-running on an
already-migrated tree is a no-op for that file (column names already
canonical, probe_name already in v2 format).
"""
from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd
import typer

from probe_designer.io.probe_schema import (
    CHEM_DRNA,
    CHEM_ILOCK,
    CHEM_RT,
    FINAL_PADLOCK_COLUMNS,
    FINAL_RT_PRIMER_COLUMNS,
    LEGACY_CHEMISTRY_VALUES,
    LEGACY_RENAMES,
    apply_final_column_order,
    make_padlock_name,
    make_rt_primer_name,
    parse_padlock_name,
    parse_rt_primer_name,
)

logger = logging.getLogger(__name__)


# Filenames the CLI sweeps for under <experiment>/output/ and its subdirs.
_PADLOCK_GLOBS = (
    "*final_probes.xlsx",
    "*final_iLock_probes.xlsx",
    "probes_combined.xlsx",
    "probe_info.xlsx",
    "probes.xlsx",
)
_RT_GLOBS = (
    "*RT_primers.xlsx",
    "rt_primers.xlsx",
)


def _discover_files(target: Path) -> tuple[List[Path], List[Path]]:
    """Return (padlock_files, rt_files) under ``target``.

    ``target`` may be a single .xlsx or a directory. For directories, look
    under both ``target`` itself and any nested ``output/`` (TCR runs nest
    one level deeper under ``output/<ts>_TCR/``, so we recurse).
    """
    if target.is_file():
        name = target.name.lower()
        if "rt" in name and "primer" in name:
            return [], [target]
        return [target], []

    padlocks: List[Path] = []
    rts: List[Path] = []
    for pattern in _PADLOCK_GLOBS:
        padlocks.extend(sorted(target.rglob(pattern)))
    for pattern in _RT_GLOBS:
        rts.extend(sorted(target.rglob(pattern)))
    # rglob can pick up files matching multiple patterns; dedup while
    # preserving order.
    seen: set = set()
    padlocks = [p for p in padlocks if not (p in seen or seen.add(p))]
    seen.clear()
    rts = [p for p in rts if not (p in seen or seen.add(p))]
    return padlocks, rts


def _apply_renames(df: pd.DataFrame) -> pd.DataFrame:
    """Apply ``LEGACY_RENAMES`` to ``df.columns``, dropping flagged columns."""
    rename_map: dict = {}
    drop_cols: list = []
    for col in df.columns:
        target = LEGACY_RENAMES.get(col)
        if target == "drop":
            drop_cols.append(col)
        elif target is not None and target != col:
            # Avoid colliding with a column already named the v2 target name
            # (would happen if we tried to migrate an already-migrated file).
            if target in df.columns:
                drop_cols.append(col)
            else:
                rename_map[col] = target
    if rename_map:
        df = df.rename(columns=rename_map)
    if drop_cols:
        df = df.drop(columns=drop_cols)
    return df


def _normalize_chemistry_values(df: pd.DataFrame) -> pd.DataFrame:
    """Map legacy chemistry tokens to v2 + fold ``iLock`` boolean back in."""
    if "chemistry" not in df.columns:
        return df
    df = df.copy()
    df["chemistry"] = df["chemistry"].replace(LEGACY_CHEMISTRY_VALUES)
    return df


def _fold_ilock_column(df: pd.DataFrame) -> pd.DataFrame:
    """If the legacy ``iLock`` column is still present, set
    ``chemistry == "iLock"`` for rows where it is yes-ish, then drop it."""
    if "iLock" not in df.columns:
        return df
    df = df.copy()
    yes_mask = df["iLock"].astype(str).str.strip().str.lower().isin(
        {"yes", "true", "1", "y"}
    )
    if "chemistry" not in df.columns:
        df["chemistry"] = pd.NA
    df.loc[yes_mask, "chemistry"] = CHEM_ILOCK
    return df.drop(columns=["iLock"])


def _regenerate_padlock_names(df: pd.DataFrame, codebook: str) -> pd.DataFrame:
    """Rebuild ``probe_name`` for every row using the v2 format."""
    if df.empty:
        return df
    df = df.copy()
    new_names: list = []
    for _, row in df.iterrows():
        # Already-v2 names round-trip cleanly; skip regeneration to preserve
        # the original codebook in case the caller's `--codebook` differs.
        existing = row.get("probe_name")
        if isinstance(existing, str):
            parsed = parse_padlock_name(existing)
            if parsed is not None:
                new_names.append(existing)
                continue
        # TCR rows use clone_id as the leading token; mutation/mRNA use gene_name.
        clone_id = row.get("clone_id")
        leading = (str(clone_id) if clone_id and not pd.isna(clone_id)
                   else str(row.get("gene_name", "")))
        # mutation rows store transcript-relative position under `position`;
        # for genomic-only rows it equals `genomic_start`.
        position = row.get("position")
        if pd.isna(position):
            position = row.get("genomic_start")
        chemistry = row.get("chemistry", CHEM_DRNA)
        no = row.get("No.")
        new_names.append(make_padlock_name(
            name=leading, position=int(position),
            chemistry=str(chemistry), codebook=codebook, no=int(no),
        ))
    df["probe_name"] = new_names
    if "codebook" not in df.columns or df["codebook"].isna().all():
        df["codebook"] = codebook
    return df


def _regenerate_rt_primer_names(
    df: pd.DataFrame, codebook: Optional[str] = None,
) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    new_names: list = []
    new_paired: list = []
    paired_col = "paired_padlock_probe_name" in df.columns
    for _, row in df.iterrows():
        leading_for_self = str(row.get("clone_id") or row.get("gene_name", ""))
        position = row.get("position")
        if pd.isna(position):
            position = row.get("genomic_start")

        existing = row.get("probe_name")
        if isinstance(existing, str) and parse_rt_primer_name(existing) is not None:
            new_names.append(existing)
        else:
            new_names.append(make_rt_primer_name(
                name=leading_for_self, position=int(position),
            ))

        if paired_col:
            paired = row.get("paired_padlock_probe_name")
            paired_no = row.get("paired_padlock_no")
            # Regenerate the link so it matches the migrated padlock file.
            # cDNA is the only chemistry that pairs with an RT primer.
            if (codebook
                    and not pd.isna(paired_no)
                    and not (isinstance(paired, str)
                             and parse_rt_primer_name(paired) is None
                             and "_cDNA_" + codebook + "_" in paired)):
                # `paired` is either NaN or in legacy format -> rebuild it.
                # Use the gene_name (mutation) or clone_id (TCR) and the
                # primer's own position. For mutation files position == the
                # genomic_start of the mutation site, which is what the
                # paired padlock used too.
                new_paired.append(make_padlock_name(
                    name=leading_for_self, position=int(position),
                    chemistry="cDNA", codebook=codebook, no=int(paired_no),
                ))
            else:
                new_paired.append(paired)
    df["probe_name"] = new_names
    if paired_col:
        df["paired_padlock_probe_name"] = new_paired
    if "chemistry" not in df.columns:
        df["chemistry"] = CHEM_RT
    else:
        df["chemistry"] = df["chemistry"].fillna(CHEM_RT).replace(
            {"": CHEM_RT}
        )
    return df


def _migrate_one(
    fp: Path, *, kind: str, codebook: Optional[str], dry_run: bool,
) -> dict:
    """Migrate a single Excel file. Returns a stats dict.

    Schema-v2 (Phase 0) keeps ``g_content`` (G-only) and ``gc_content`` (G+C)
    as distinct columns. Legacy files have only ``g_content`` and store
    G-only values there; migration preserves them and emits an empty
    ``gc_content`` column (the user must re-run design to populate it). We
    log a one-line warning per such file so this is not silently missed.
    """
    df = pd.read_excel(fp)
    original_cols = df.columns.tolist()
    has_legacy_g_only = (
        kind == "padlock"
        and "g_content" in original_cols
        and "gc_content" not in original_cols
    )
    if has_legacy_g_only:
        logger.warning(
            "%s: legacy file has g_content but no gc_content. The migrated file "
            "will contain g_content (G-only fraction, kept for back-compat) plus "
            "an empty gc_content column. Re-run probe-design to populate gc_content.",
            fp,
        )
    # Legacy mutation RT primer files use the standalone `No.` column to
    # store the *paired* padlock's number (RT primers don't get their own
    # codebook slot). Rewrite it before the generic renames run, so we
    # don't lose the value when v2's "RT primers don't have No." rule kicks
    # in below.
    if kind == "rt" and "No." in df.columns and "Padlock_Probe_No." not in df.columns:
        df = df.rename(columns={"No.": "paired_padlock_no"})
    df = _apply_renames(df)
    df = _fold_ilock_column(df)
    df = _normalize_chemistry_values(df)

    if kind == "padlock":
        if codebook is None:
            raise typer.BadParameter(
                f"--codebook is required to migrate padlock file {fp}; pass "
                "the codebook name (e.g. SP369) used when the panel was designed."
            )
        df = _regenerate_padlock_names(df, codebook)
        df = apply_final_column_order(df, kind="padlock")
    elif kind == "rt":
        df = _regenerate_rt_primer_names(df, codebook=codebook)
        df = apply_final_column_order(df, kind="rt")
    else:
        raise ValueError(f"unknown kind {kind!r}")

    stats = {
        "file": str(fp),
        "rows": len(df),
        "columns_before": len(original_cols),
        "columns_after": len(df.columns),
    }
    if dry_run:
        return stats

    bak = fp.with_suffix(fp.suffix + ".bak")
    if not bak.exists():
        shutil.copy2(fp, bak)
    df.to_excel(fp, index=False, engine="openpyxl")
    return stats


def migrate_columns(
    target: Path = typer.Argument(
        ..., exists=True,
        help="Experiment directory or a single probe Excel file."
    ),
    codebook: Optional[str] = typer.Option(
        None, "--codebook",
        help="Codebook name (e.g. SP369). Required for padlock files since "
             "the legacy `_SP_<n>` substring carries the No. but not the "
             "codebook. RT primer files don't need this."
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run",
        help="Discover files and report planned changes; don't write anything."
    ),
) -> None:
    """Upgrade legacy probe Excel files to Schema v2.

    Writes a ``.bak`` file alongside each migrated input. Re-running on a
    tree that's already been migrated is a no-op for those files.
    """
    target = target.resolve()
    padlocks, rts = _discover_files(target)
    if not padlocks and not rts:
        typer.echo(f"No probe Excel files found under {target}.")
        raise typer.Exit(code=1)

    typer.echo(f"Migrating {len(padlocks)} padlock file(s) + "
               f"{len(rts)} RT primer file(s) under {target}")

    failures: List[str] = []
    for fp in padlocks:
        try:
            stats = _migrate_one(fp, kind="padlock", codebook=codebook,
                                  dry_run=dry_run)
            typer.echo(f"  [padlock] {fp.relative_to(target.parent)} : "
                       f"{stats['rows']} rows, "
                       f"{stats['columns_before']} -> {stats['columns_after']} cols")
        except Exception as e:
            msg = f"  [padlock] FAILED {fp}: {e}"
            failures.append(msg)
            typer.echo(msg, err=True)

    for fp in rts:
        try:
            stats = _migrate_one(fp, kind="rt", codebook=codebook,
                                  dry_run=dry_run)
            typer.echo(f"  [rt]      {fp.relative_to(target.parent)} : "
                       f"{stats['rows']} rows, "
                       f"{stats['columns_before']} -> {stats['columns_after']} cols")
        except Exception as e:
            msg = f"  [rt]      FAILED {fp}: {e}"
            failures.append(msg)
            typer.echo(msg, err=True)

    if dry_run:
        typer.echo("[dry-run] no files written.")
    if failures:
        typer.echo(f"\n{len(failures)} file(s) failed to migrate.", err=True)
        raise typer.Exit(code=1)
