"""Root Typer application for ``probe-design`` CLI.

Subcommands are registered from sibling modules so each stays focused.
"""
from __future__ import annotations

import logging

import typer

from probe_designer.cli import design as design_cmd
from probe_designer.cli import score as score_cmd
from probe_designer.cli import select as select_cmd


app = typer.Typer(
    name="probe-design",
    help="Padlock probe design pipeline (HuangLab).",
    no_args_is_help=True,
    pretty_exceptions_enable=False,
)

app.command("design", help="Run the full design pipeline (search, filter, BLAST, score, rank, select).")(design_cmd.design)
app.command("score", help="Score and peak-rank a saved binding_sites JSON (no BLAST).")(score_cmd.score)
app.command("select", help="Apply top-N + min-gap selection to a scored binding_sites JSON.")(select_cmd.select)


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
