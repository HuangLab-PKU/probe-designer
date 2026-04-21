"""YAML config loader with ${ENV_VAR} and ${ENV_VAR:-default} expansion.

Used by the Phase 1+ CLI and Pipeline to let users keep credentials and host
paths out of YAML files. Expansion rules:

    ${FOO}         -> value of env var FOO; raises ConfigError if missing
    ${FOO:-bar}    -> value of FOO, or literal 'bar' if FOO is unset/empty

Expansion happens at the text level BEFORE YAML parsing so the substituted
text is still valid YAML. Non-string YAML values are left untouched.
"""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Mapping, Optional, Union

import yaml

from probe_designer.utils.errors import ConfigError


_ENV_VAR_PATTERN = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)(?::-(.*?))?\}")


def _expand_env(text: str, env: Mapping[str, str]) -> str:
    """Substitute ${VAR} and ${VAR:-default} markers in text using env.

    Raises ConfigError on any unresolved ${VAR} without a default.
    """
    def repl(match: re.Match[str]) -> str:
        var_name = match.group(1)
        default = match.group(2)  # None if no ':-default' present
        value = env.get(var_name)
        if value is not None and value != "":
            return value
        if default is not None:
            return default
        raise ConfigError(
            f"Unresolved environment variable: ${{{var_name}}}. "
            f"Set it in the shell or in a .env file. "
            f"See designer/.env.example for the expected variables."
        )

    return _ENV_VAR_PATTERN.sub(repl, text)


def load_yaml_with_env(
    path: Union[str, Path],
    env: Optional[Mapping[str, str]] = None,
) -> Any:
    """Load a YAML file expanding ${ENV_VAR} placeholders first.

    Args:
        path: path to YAML file
        env: environment mapping to resolve placeholders against. Defaults to
            os.environ when None (typical production use). Tests typically
            pass a fresh dict for isolation.

    Returns:
        Parsed YAML contents (typically a dict).

    Raises:
        ConfigError: if a ${VAR} without default cannot be resolved.
        FileNotFoundError: if the YAML file does not exist.
    """
    if env is None:
        env = dict(os.environ)

    raw = Path(path).read_text(encoding="utf-8")
    expanded = _expand_env(raw, env)
    return yaml.safe_load(expanded)
