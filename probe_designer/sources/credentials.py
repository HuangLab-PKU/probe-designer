"""NCBI Entrez credential resolution.

Replaces the hardcoded ``Entrez.email`` / ``Entrez.api_key`` lines at
database.py:28-29. Credentials are now resolved from environment variables
(loaded from a ``.env`` file if present) so no secrets are committed.

Usage from consumer code::

    from Bio import Entrez
    from probe_designer.sources.credentials import get_entrez_credentials

    email, api_key = get_entrez_credentials(required=False)
    if email:
        Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
"""
from __future__ import annotations

import logging
import os
from typing import Optional, Tuple


logger = logging.getLogger(__name__)


_LOADED_DOTENV = False


def _load_env_file_once() -> None:
    """Load a local .env file into os.environ if python-dotenv is available.

    Only runs on first call and only mutates os.environ for variables that
    are not already set. Safe to call repeatedly.
    """
    global _LOADED_DOTENV
    if _LOADED_DOTENV:
        return
    try:
        from dotenv import load_dotenv  # type: ignore[import-untyped]
    except ImportError:
        _LOADED_DOTENV = True
        return

    for candidate in (".env", "designer/.env"):
        if os.path.exists(candidate):
            load_dotenv(candidate, override=False)
            logger.debug("Loaded env file: %s", candidate)
            break
    _LOADED_DOTENV = True


def get_entrez_credentials(
    *,
    required: bool = False,
) -> Tuple[Optional[str], Optional[str]]:
    """Return ``(email, api_key)`` from env vars ``ENTREZ_EMAIL`` / ``ENTREZ_API_KEY``.

    Args:
        required: if True, raise :class:`probe_designer.utils.errors.ConfigError`
            when either is missing. Default False lets Entrez run without
            credentials (subject to NCBI rate limits).

    Returns:
        Tuple of (email, api_key). Either may be ``None`` when ``required=False``.
    """
    _load_env_file_once()
    email = os.environ.get("ENTREZ_EMAIL") or None
    api_key = os.environ.get("ENTREZ_API_KEY") or None

    if required and (not email or not api_key):
        from probe_designer.utils.errors import ConfigError
        missing = [
            name for name, val in (("ENTREZ_EMAIL", email), ("ENTREZ_API_KEY", api_key))
            if not val
        ]
        raise ConfigError(
            f"NCBI Entrez credentials missing: {', '.join(missing)}. "
            f"Set them in the shell or in a .env file. "
            f"See designer/.env.example for the expected variables."
        )

    return email, api_key
