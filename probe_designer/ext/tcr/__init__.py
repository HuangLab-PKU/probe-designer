"""TCR padlock probe extension.

Public API:
    TcrProbeDesigner       — CDR3-constrained scan emitting dual-chemistry arms.
    validate_cdr3_in_trb   — input check; raises ValueError if CDR3 not in TRB.
    TcrConfig              — driver config dataclass (CLI / YAML).
    run_tcr_pipeline       — end-to-end orchestrator (find BDS → BLAST QC → assemble → RT primers).
"""
from probe_designer.ext.tcr.probe import TcrProbeDesigner, validate_cdr3_in_trb

# config + pipeline are imported lazily so probe.py + the shim work even
# during incremental development of this package.
try:
    from probe_designer.ext.tcr.config import TcrConfig  # noqa: F401
    from probe_designer.ext.tcr.pipeline import run_tcr_pipeline  # noqa: F401
    _PIPELINE_AVAILABLE = True
except ImportError:
    _PIPELINE_AVAILABLE = False

__all__ = ["TcrProbeDesigner", "validate_cdr3_in_trb"]
if _PIPELINE_AVAILABLE:
    __all__.extend(["TcrConfig", "run_tcr_pipeline"])
