"""Root Typer application for ``probe-design`` CLI.

Subcommands are registered from sibling modules so each stays focused.
"""
from __future__ import annotations

import logging

import typer

from probe_designer.cli import assemble as assemble_cmd
from probe_designer.cli import design as design_cmd
from probe_designer.cli import merge as merge_cmd
from probe_designer.cli import migrate_columns as migrate_columns_cmd
from probe_designer.cli import score as score_cmd
from probe_designer.cli import select as select_cmd
from probe_designer.cli import validate as validate_cmd
from probe_designer.ext.pool import pool_app
from probe_designer.ext.mutation.cli import command as mutation_cmd
from probe_designer.ext.tcr.cli import command as tcr_cmd


app = typer.Typer(
    name="probe-design",
    help="Padlock probe design pipeline (HuangLab).",
    no_args_is_help=True,
    pretty_exceptions_enable=False,
)

app.command("design", help="Run the full design pipeline (search, filter, BLAST, score, rank, select).")(design_cmd.design)
app.command("score", help="Score and peak-rank a saved binding_sites JSON (no BLAST).")(score_cmd.score)
app.command("select", help="Apply top-N + min-gap selection to a scored binding_sites JSON.")(select_cmd.select)
app.command("merge", help="Merge filtered_binding_sites.xlsx files from multiple runs into one XLSX.")(merge_cmd.merge)
app.command("assemble", help="Attach backbones to binding sites to produce final probe sequences.")(assemble_cmd.assemble)
app.command("validate", help="Load and validate a YAML config (with env-var expansion).")(validate_cmd.validate)
app.command("migrate-columns", help="Upgrade legacy probe Excel files to Schema v2 (.bak written alongside).")(migrate_columns_cmd.migrate_columns)
app.command("mutation", help="Design padlock probes for somatic mutations (iLock / dRNA / cDNA chemistries).")(mutation_cmd)
app.command("tcr", help="Design padlock probes for TCR clone CDR3 regions (dRNA / cDNA chemistries; default dRNA only).")(tcr_cmd)
app.add_typer(pool_app, name="pool", help="Pool-level operations (cross-ligation audit).")


@app.callback()
def _global_options(
    log_level: str = typer.Option(
        "INFO", "--log-level", "-l",
        help="Python logging level (DEBUG, INFO, WARNING, ERROR)."
    ),
) -> None:
    """Global flags applied before any subcommand runs."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


if __name__ == "__main__":
    app()
