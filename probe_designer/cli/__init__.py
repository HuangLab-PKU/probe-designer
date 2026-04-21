"""probe-design Typer CLI.

The ``probe-design`` console script (registered in pyproject.toml)
points at :data:`probe_designer.cli.app.app`.
"""
from probe_designer.cli.app import app

__all__ = ["app"]
