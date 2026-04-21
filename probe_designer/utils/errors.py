"""Custom exception types for probe_designer.

Phase 1 introduces ConfigError for the new env-aware config loader.
Future phases will add BlastError, GenomeError, etc. as the library grows.
"""
from __future__ import annotations


class ProbeDesignerError(Exception):
    """Base exception for all probe_designer-specific errors."""


class ConfigError(ProbeDesignerError):
    """Raised when configuration is invalid or cannot be resolved.

    Typical triggers:
      - ${ENV_VAR} placeholder in YAML with no matching env var and no default
      - unknown config key in strict mode
      - YAML parse error
    """
