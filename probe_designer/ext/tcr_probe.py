"""Backwards-compatibility shim.

The TCR extension was promoted from a single module to a package at
``probe_designer.ext.tcr``. New code should import from there. This shim
preserves the legacy import path used by older scripts.
"""
from probe_designer.ext.tcr.probe import (  # noqa: F401
    TcrProbeDesigner,
    validate_cdr3_in_trb,
    validate_subseq_in_target,
)
