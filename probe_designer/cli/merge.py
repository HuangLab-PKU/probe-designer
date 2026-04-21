"""``probe-design merge`` — combine binding-site XLSX outputs across runs.

Walks subdirectories of a results folder, collects every
``filtered_binding_sites.xlsx``, keeps rows whose ``gene_name`` is listed
in the user-supplied ``gene_info.xlsx`` (case-insensitive), drops
duplicates on ``(gene_name, st, en)``, sorts by gene_info row order then
``(st, en)``, and writes the merged XLSX plus a ``missing_genes.txt``
listing genes with no surviving rows.

Replaces the pre-0.2.1 ``scripts/merge_results.py`` utility at one-third
the line count.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import pandas as pd
import typer


logger = logging.getLogger(__name__)


def _find_xlsx_files(results_dir: Path) -> List[Path]:
    """Every ``filtered_binding_sites.xlsx`` directly inside a subdir of ``results_dir``."""
    if not results_dir.exists():
        raise FileNotFoundError(f"Results directory does not exist: {results_dir}")
    return sorted(
        sub / "filtered_binding_sites.xlsx"
        for sub in results_dir.iterdir()
        if sub.is_dir() and (sub / "filtered_binding_sites.xlsx").exists()
    )


def _load_gene_info(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path)
    if "gene_name" not in df.columns:
        raise ValueError(f"gene_info missing 'gene_name' column: {sorted(df.columns)}")
    df["gene_name"] = df["gene_name"].astype(str).str.strip()
    df["_order"] = range(len(df))  # used for deterministic sort
    df["_gene_upper"] = df["gene_name"].str.upper()
    return df


def merge(
    results_dir: Path = typer.Option(
        ..., "--results-dir", "-r", exists=True, file_okay=False,
        help="Root directory containing per-run subdirs each with a "
             "filtered_binding_sites.xlsx."
    ),
    gene_info: Path = typer.Option(
        ..., "--gene-info", "-g", exists=True, dir_okay=False,
        help="Gene info XLSX (must have 'gene_name' column)."
    ),
    output: Path = typer.Option(
        ..., "--output", "-o",
        help="Output XLSX path for the merged binding sites."
    ),
    missing_output: Optional[Path] = typer.Option(
        None, "--missing-output", "-m",
        help="Missing-genes text file (default: <output's dir>/missing_genes.txt)."
    ),
) -> None:
    """Merge filtered binding sites across multiple design runs."""
    files = _find_xlsx_files(results_dir)
    if not files:
        typer.echo(
            f"No filtered_binding_sites.xlsx files found under {results_dir}",
            err=True,
        )
        raise typer.Exit(code=1)
    logger.info("Merging %d files from %s", len(files), results_dir)

    dfs = []
    for p in files:
        try:
            d = pd.read_excel(p)
        except Exception as exc:
            logger.warning("Skipping unreadable %s: %s", p, exc)
            continue
        if "gene_name" not in d.columns:
            logger.warning("Skipping %s: missing 'gene_name' column", p)
            continue
        d["gene_name"] = d["gene_name"].astype(str).str.strip()
        dfs.append(d)

    if not dfs:
        typer.echo("No loadable binding-sites files.", err=True)
        raise typer.Exit(code=1)

    merged = pd.concat(dfs, ignore_index=True)

    # Case-insensitive gene_info filter + sort key via pandas.merge
    info = _load_gene_info(gene_info)
    merged["_gene_upper"] = merged["gene_name"].str.upper()
    kept = merged.merge(
        info[["_gene_upper", "_order"]], on="_gene_upper", how="inner"
    )

    # Dedup uses the case-normalized gene column so "ACTB" and "actb" collapse.
    dedup_cols = ["_gene_upper", "st", "en"] if {"st", "en"}.issubset(kept.columns) else ["_gene_upper", "sequence"]
    kept = kept.drop_duplicates(subset=dedup_cols, keep="first")

    sort_cols = ["_order"] + [c for c in ("st", "en") if c in kept.columns]
    kept = kept.sort_values(by=sort_cols).drop(columns=["_order", "_gene_upper"])

    output = output.with_suffix(".xlsx") if output.suffix.lower() not in (".xlsx", ".xls") else output
    output.parent.mkdir(parents=True, exist_ok=True)
    kept.to_excel(str(output), index=False)
    logger.info("Wrote %d merged rows (%d genes) -> %s",
                len(kept), kept["gene_name"].nunique(), output)

    # Missing genes — preserve gene_info row order
    kept_upper = set(kept["gene_name"].str.upper())
    missing = [g for g in info["gene_name"] if g.upper() not in kept_upper]
    missing_path = missing_output or (output.parent / "missing_genes.txt")
    if missing:
        missing_path.write_text("\n".join(missing) + "\n", encoding="utf-8")
        typer.echo(
            f"Wrote {len(kept)} rows -> {output}; "
            f"{len(missing)} genes missing -> {missing_path}"
        )
    else:
        typer.echo(f"Wrote {len(kept)} rows -> {output}; all genes covered.")
