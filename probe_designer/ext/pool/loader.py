"""Load an existing pool's effective probes and convert to ProbeForScreen.

Lazy-imports ``probe_book.pool.{manifest,loader}`` and
``probe_book.probe.registry`` inside the function body so that:

* the top-level ``import probe_designer.ext.pool.loader`` succeeds even
  if probe-book is not installed;
* the actual call raises a clear ImportError with install hints.

This is the ONLY place in ``probe_designer`` where ``probe_book`` is
imported. The static check at ``designer/tests/unit/test_no_bank_import_in_qc.py``
enforces that ``probe_designer.qc.*`` stays free of bank imports.
"""
from __future__ import annotations

from pathlib import Path
from typing import List

from probe_designer.qc.cross_ligation import ProbeForScreen


def load_pool_as_probes_for_screen(
    pool_id: str, repo_root: Path,
) -> List[ProbeForScreen]:
    """Load a pool's effective probe list and convert to ProbeForScreen.

    Reads ``repo_root/pools/<pool_id>/manifest.yaml`` and resolves probes
    via ``repo_root/probes/registry.tsv``. Only probes with both
    ``probe_arm5`` AND ``probe_arm3`` are returned (RT primers and other
    non-binding entries are skipped — they don't participate in
    cross-ligation).

    Args:
        pool_id: e.g. ``"TNBCmarker_main_v1"``.
        repo_root: project root containing ``pools/`` and ``probes/``.

    Returns:
        List of :class:`ProbeForScreen` ready for
        :func:`probe_designer.qc.cross_ligation.screen_cross_ligation_v2`.

    Raises:
        ImportError: when ``probe_book`` is not installed — includes an
            actionable hint pointing at ``bank/``.
        FileNotFoundError: when the pool manifest doesn't exist at the
            expected path.
    """
    try:
        from probe_book.pool.manifest import PoolManifest
        from probe_book.pool.loader import compute_effective_probes
        from probe_book.probe.registry import ProbeRegistry
    except ImportError as exc:
        raise ImportError(
            "probe-book is required for --target-pool / probe-design pool check. "
            "Install via `pip install -e bank/` from the project root."
        ) from exc

    repo_root = Path(repo_root).resolve()
    manifest_path = repo_root / "pools" / pool_id / "manifest.yaml"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Pool manifest not found at {manifest_path}. "
            f"Pass --repo-root if `pools/{pool_id}/` lives elsewhere."
        )

    pool = PoolManifest.from_yaml(manifest_path)
    registry = ProbeRegistry(repo_root)
    effective = compute_effective_probes(pool, probe_registry=registry)

    out: List[ProbeForScreen] = []
    for ep in effective:
        if not (ep.probe_arm5 and ep.probe_arm3):
            continue
        out.append(ProbeForScreen(
            probe_id=ep.probe_id,
            chemistry=ep.chemistry or "dRNA",
            probe_arm5=ep.probe_arm5,
            probe_arm3=ep.probe_arm3,
            sequence=(ep.sequence or "").upper() or (ep.probe_arm5 + ep.probe_arm3),
            target=ep.target,
        ))
    return out


__all__ = ["load_pool_as_probes_for_screen"]
