"""Optional ext/pool extension — bridges bank's pool ledger to designer's qc.

This package is the ONLY place in ``probe_designer`` that imports
``probe_book``. Imports inside the module functions are LAZY so the
top-level ``from probe_designer.ext.pool import ...`` works even if
``probe_book`` is not installed; the actual call sites raise a clear
error in that case.

Public API:
    * :func:`load_pool_as_probes_for_screen` — load a pool's existing
      probes and convert them to :class:`ProbeForScreen` for v2 cross-lig.
    * :data:`pool_app` — Typer subcommand group registered as
      ``probe-design pool``.
"""
from probe_designer.ext.pool.loader import load_pool_as_probes_for_screen
from probe_designer.ext.pool.cli import pool_app

__all__ = ["load_pool_as_probes_for_screen", "pool_app"]
