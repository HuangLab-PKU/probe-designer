#!/usr/bin/env python3
"""DEPRECATED — use ``probe-design design --strategy single_sequence``.

Legacy shim for users still calling this script path directly.
Will be removed in probe-designer 0.3.0.
"""
from __future__ import annotations

import argparse
import sys
import warnings


def _parse_legacy_args(argv: list[str]) -> list[str]:
    p = argparse.ArgumentParser(description="Legacy find_probes shim")
    p.add_argument("--config")
    p.add_argument("--genes_file", "--genes-file", dest="genes_file")
    p.add_argument("--gene_info", "--gene-info", dest="gene_info")
    p.add_argument("--species")
    p.add_argument("--database")
    p.add_argument("--blast_species", "--blast-species", dest="blast_species", nargs="+")
    p.add_argument("--genome_fasta", "--genome-fasta", dest="genome_fasta")
    p.add_argument("--output_dir", "--output", dest="output_dir")
    p.add_argument("--progress", action="store_true")
    a, _ = p.parse_known_args(argv)

    args = ["design", "--strategy", "single_sequence"]
    if a.config:
        args += ["--config", a.config]
    if a.genes_file:
        args += ["--genes-file", a.genes_file]
    if a.species:
        args += ["--species", a.species]
    if a.blast_species:
        for s in a.blast_species:
            args += ["--blast-species", s]
    if a.genome_fasta:
        args += ["--genome-fasta", a.genome_fasta]
    if a.output_dir:
        args += ["--output", a.output_dir]
    return args


def main() -> int:
    warnings.warn(
        "scripts/find_probes.py is deprecated; use "
        "'probe-design design --strategy single_sequence' instead. "
        "This shim will be removed in probe-designer 0.3.0.",
        DeprecationWarning,
        stacklevel=2,
    )
    print(
        "[deprecation] find_probes.py: please migrate to "
        "'probe-design design --strategy single_sequence'.",
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
