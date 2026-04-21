"""Pluggable hooks for :class:`probe_designer.pipeline.Pipeline`.

These Protocols let callers inject behavior into the pipeline without
designer depending on webapp or DB internals. Every hook is optional —
the Pipeline provides no-op defaults for each.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Protocol, runtime_checkable


@runtime_checkable
class ProgressHook(Protocol):
    """Notified as genes flow through pipeline stages.

    Stages are strings:
      - ``"sequence"`` fetching per-gene reference sequence
      - ``"isoforms"`` fetching isoform metadata
      - ``"search"`` running the search strategy
      - ``"pre_blast"`` thermal filter
      - ``"blast"`` specificity filter
      - ``"score"`` scoring + peak_rank + selection
    """

    def on_gene_start(self, gene: str, index: int, total: int) -> None: ...

    def on_stage(self, gene: str, stage: str, info: Dict[str, Any]) -> None: ...

    def on_gene_complete(self, gene: str, site_count: int) -> None: ...

    def on_pipeline_complete(self, gene_count: int) -> None: ...


@runtime_checkable
class ExistingTargetsHook(Protocol):
    """Look up already-designed targets for a gene (e.g., from a DB).

    The webapp supplies an adapter that queries the ``TargetRegion`` table;
    the CLI leaves this hook unset, so no prior targets are injected.
    """

    def lookup(
        self, gene: str, species: str, min_score: float = 0.0
    ) -> List[Dict[str, Any]]: ...


@runtime_checkable
class IsoformProvider(Protocol):
    """Return per-gene isoform metadata for isoform-aware strategies.

    The default implementation chains local GTF → Ensembl REST. The webapp
    wraps the default with a DB cache.
    """

    def get_isoforms(self, gene: str, species: str) -> Optional[List[Dict[str, Any]]]: ...


# --- Default no-op implementations -----------------------------------------

class _NullProgressHook:
    def on_gene_start(self, gene: str, index: int, total: int) -> None:
        return

    def on_stage(self, gene: str, stage: str, info: Dict[str, Any]) -> None:
        return

    def on_gene_complete(self, gene: str, site_count: int) -> None:
        return

    def on_pipeline_complete(self, gene_count: int) -> None:
        return


class _NullExistingTargetsHook:
    def lookup(
        self, gene: str, species: str, min_score: float = 0.0
    ) -> List[Dict[str, Any]]:
        return []


NULL_PROGRESS: ProgressHook = _NullProgressHook()
NULL_EXISTING: ExistingTargetsHook = _NullExistingTargetsHook()
