"""``probe-design validate`` — load config with env-expansion and print it.

Non-zero exit code on any validation error so CI/CD can gate on it.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import typer

from probe_designer.config import ConfigManager
from probe_designer.config.loader import load_yaml_with_env
from probe_designer.utils.errors import ConfigError


logger = logging.getLogger(__name__)


def validate(
    config_path: Path = typer.Option(
        ..., "--config", exists=True, dir_okay=False, readable=True,
        help="YAML config to load and validate."
    ),
    env_file: Optional[Path] = typer.Option(
        None, "--env-file", exists=True, dir_okay=False, readable=True,
        help="Optional .env file to load before expansion."
    ),
) -> None:
    """Load YAML (expanding ${VAR} placeholders), hydrate, and print resolved config."""
    # Load .env if requested
    if env_file is not None:
        try:
            from dotenv import load_dotenv
            load_dotenv(env_file, override=False)
        except ImportError:
            typer.echo("python-dotenv not installed; --env-file ignored", err=True)

    # 1. Placeholder expansion (fail-fast on unresolved ${VAR})
    try:
        expanded = load_yaml_with_env(config_path)
    except ConfigError as exc:
        typer.echo(f"Config error: {exc}", err=True)
        raise typer.Exit(code=2)

    typer.echo("=== Resolved YAML ===")
    typer.echo(json.dumps(expanded, indent=2, ensure_ascii=False, default=str))

    # 2. Dataclass hydration + field validation
    try:
        cfg = ConfigManager(str(config_path))
    except Exception as exc:
        typer.echo(f"Failed to hydrate ConfigManager: {exc}", err=True)
        raise typer.Exit(code=2)

    errors = cfg.validate_config()
    if errors:
        typer.echo("\n=== Validation errors ===", err=True)
        for err in errors:
            typer.echo(f"  - {err}", err=True)
        raise typer.Exit(code=1)

    typer.echo("\nConfig OK.")
