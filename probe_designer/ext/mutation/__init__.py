"""Mutation padlock probe extension.

Public API:
    MutationProbeDesigner          — base designer (default geometry; cDNA pass).
    MutationProbeDesignerInvader   — iLock invader chemistry (arm5[0] == arm3[-1] == RC(SNP)).
    MutationProbeDesigner3End      — dRNA chemistry (direct-RNA padlock; mutation at probe 3' tip).
    verify_iLock_probe             — runtime guard for canonical iLock layout + invader invariant.
    assemble_ilock                 — assemble + verify a canonical iLock probe.
    design_mutation_probes         — module-level convenience wrapper around the designer classes.
"""
from probe_designer.ext.mutation.probe import (
    MutationProbeDesigner,
    MutationProbeDesigner3End,
    MutationProbeDesignerInvader,
    assemble_ilock,
    design_mutation_probes,
    verify_iLock_probe,
)

__all__ = [
    "MutationProbeDesigner",
    "MutationProbeDesigner3End",
    "MutationProbeDesignerInvader",
    "assemble_ilock",
    "design_mutation_probes",
    "verify_iLock_probe",
]
