#!/usr/bin/env python3
"""DEPRECATED — use ``probe-design assemble``.

Legacy shim for users still calling this script path directly.
Will be removed in probe-designer 0.3.0.
"""
from __future__ import annotations

import argparse
import sys
import warnings


def _parse_legacy_args(argv: list[str]) -> list[str]:
    p = argparse.ArgumentParser(description="Legacy probe_assemble shim")
    p.add_argument("--binding_sites", "--binding-sites", dest="binding_sites")
    p.add_argument("--gene_info", "--gene-info", dest="gene_info")
    p.add_argument("--backbone_file", "--backbone-file", "--backbone", dest="backbone_file")
    p.add_argument("--output_dir", "--output", "-o", dest="output_dir")
    p.add_argument("--ilock")
    a, _ = p.parse_known_args(argv)

    args = ["assemble"]
    if a.binding_sites:
        args += ["--binding-sites", a.binding_sites]
    if a.gene_info:
        args += ["--gene-info", a.gene_info]
    if a.backbone_file:
        args += ["--backbone", a.backbone_file]
    if a.output_dir:
        args += ["--output", a.output_dir]
    if a.ilock:
        args += ["--ilock", a.ilock]
    return args


def main() -> int:
    warnings.warn(
        "scripts/probe_assemble.py is deprecated; use "
        "'probe-design assemble' instead. "
        "This shim will be removed in probe-designer 0.3.0.",
        DeprecationWarning,
        stacklevel=2,
    )
    print(
        "[deprecation] probe_assemble.py: please migrate to "
        "'probe-design assemble'.",
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
