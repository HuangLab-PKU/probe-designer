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


def load_pool_site_quality(pool_id: str, repo_root: Path) -> List[dict]:
    """The derived false-positive status of every sited probe in a pool.

    Reads ``site_quality/registry.tsv`` -- one row per (site, chemistry), rebuilt
    from observations and pool recipes by ``bank/scripts/rebuild_site_quality.py``
    -- and returns one row per pooled probe that has a site, in that probe's own
    chemistry: a site can be clean in cDNA and a false positive in dRNA, and only
    the probe's chemistry says which applies to it.

    A (site, chemistry) the table does not list is ``untested``. A missing table
    raises FileNotFoundError, so "no table" is never read as "all clean".
    """
    import csv

    try:
        from probe_book.pool.manifest import PoolManifest
        from probe_book.probe.registry import ProbeRegistry
    except ImportError as exc:
        raise ImportError(
            "probe-book is required for probe-design pool check. "
            "Install via `pip install -e bank/` from the project root."
        ) from exc

    repo_root = Path(repo_root).resolve()
    table = repo_root / "site_quality" / "registry.tsv"
    if not table.exists():
        raise FileNotFoundError(
            f"{table} has not been built; run bank/scripts/rebuild_site_quality.py")
    manifest_path = repo_root / "pools" / pool_id / "manifest.yaml"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Pool manifest not found at {manifest_path}.")

    with open(table, encoding="utf-8", newline="") as fh:
        status = {(r["site_id"], r["chemistry"]): r
                  for r in csv.DictReader(fh, delimiter="\t")}
    probes = {r["probe_id"]: r for r in ProbeRegistry(repo_root).read_all()}
    out: List[dict] = []
    for entry in PoolManifest.from_yaml(manifest_path).recipe:
        probe = probes.get(entry.probe_id) or {}
        site, chemistry = probe.get("site_id", ""), probe.get("chemistry", "")
        if not site:
            continue
        row = status.get((site, chemistry), {})
        out.append({
            "probe_id": entry.probe_id,
            "site_id": site,
            "chemistry": chemistry,
            "gene_symbol": row.get("gene_symbol") or probe.get("target", ""),
            "fp_status": row.get("fp_status") or "untested",
            "candidate_set_size": row.get("candidate_set_size", ""),
            "differs_by_chemistry": row.get("differs_by_chemistry", "") == "true",
            "fp_conflict": row.get("fp_conflict", ""),
            "fp_proof": row.get("fp_proof", ""),
        })
    return out


__all__ = ["load_pool_as_probes_for_screen", "load_pool_site_quality"]
