"""Probe Designer - Padlock probe design pipeline for spatial transcriptomics."""

from probe_designer.config import ConfigManager
from probe_designer.database import DatabaseInterface
from probe_designer.filtering import SequenceFilter
from probe_designer.search_strategies import (
    BruteForceStrategy,
    IsoformConsensusStrategy,
    IsoformSpecificStrategy,
)
from probe_designer.probe_assembly import ProbeAssembler
from probe_designer.scoring import compute_target_score, peak_rank

SequenceSearchStrategy = BruteForceStrategy  # alias for clarity

__version__ = "0.1.0"
