#!/usr/bin/env python3
"""DEPRECATED — use ``probe-design merge``.

Legacy shim for users still calling this script path directly.
Will be removed in probe-designer 0.3.0.
"""
from __future__ import annotations

import argparse
import sys
import warnings


def _parse_legacy_args(argv: list[str]) -> list[str]:
    p = argparse.ArgumentParser(description="Legacy merge_results shim")
    p.add_argument("--results-dir", "-r")
    p.add_argument("--gene-info", "-g")
    p.add_argument("--output", "-o")
    p.add_argument("--missing-output", "-m")
    a, _ = p.parse_known_args(argv)

    args = ["merge"]
    if a.results_dir:
        args += ["--results-dir", a.results_dir]
    if a.gene_info:
        args += ["--gene-info", a.gene_info]
    if a.output:
        args += ["--output", a.output]
    if a.missing_output:
        args += ["--missing-output", a.missing_output]
    return args


def main() -> int:
    warnings.warn(
        "scripts/merge_results.py is deprecated; use "
        "'probe-design merge' instead. "
        "This shim will be removed in probe-designer 0.3.0.",
        DeprecationWarning,
        stacklevel=2,
    )
    print(
        "[deprecation] merge_results.py: please migrate to "
        "'probe-design merge'.",
        file=sys.stderr,
    )
    from probe_designer.cli.app import app

    try:
        app(_parse_legacy_args(sys.argv[1:]))
    except SystemExit as exc:
        return int(exc.code or 0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
