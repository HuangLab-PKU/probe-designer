"""Backwards-compatibility shim.

The mutation extension was promoted from a single module to a package at
``probe_designer.ext.mutation``. New code should import from there. This shim
preserves the legacy import path used by older tests and experiment scripts.
"""
from probe_designer.ext.mutation import (  # noqa: F401
    MutationProbeDesigner,
    MutationProbeDesigner3End,
    MutationProbeDesignerInvader,
    design_mutation_probes,
    verify_iLock_probe,
)
