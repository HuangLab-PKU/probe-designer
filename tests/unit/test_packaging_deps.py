"""Every third-party library the production code imports must be declared.

Audit R10 (2026-07-15): ``primer3-py`` was imported by the cross-ligation screen
— production code — but appeared in neither ``requirements.txt`` nor
``pyproject.toml``. It only worked because the dev env happened to have it, so a
fresh install of the package would ImportError at run time, not install time.

This is a declaration lint, not a behaviour test: it asserts the *class* of bug
(imported-but-undeclared) stays fixed rather than re-checking one library. A
library may be declared as a hard dependency or as a documented optional extra;
what it may not be is invisible.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PYPROJECT = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
REQUIREMENTS = (ROOT / "requirements.txt").read_text(encoding="utf-8")

# import name -> distribution name, for libraries production modules import.
THERMO_LIBS = {
    "primer3": "primer3-py",   # qc/cross_ligation.py, filtering/pairwise_duplex.py
    "Bio": "biopython",        # filtering/, chemistry.py — every Tm_NN call
    "RNA": "ViennaRNA",        # filtering/accessibility.py (optional extra)
}


def _declared_in_pyproject(dist: str) -> bool:
    """True if ``dist`` appears in [project].dependencies or any extra."""
    return re.search(rf'"{re.escape(dist)}\b', PYPROJECT, re.IGNORECASE) is not None


def _declared_in_requirements(dist: str) -> bool:
    """True if ``dist`` appears as an uncommented requirement line."""
    for line in REQUIREMENTS.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if re.match(rf"{re.escape(dist)}\b", line, re.IGNORECASE):
            return True
    return False


@pytest.mark.parametrize("import_name,dist", sorted(THERMO_LIBS.items()))
def test_thermo_library_is_declared(import_name, dist):
    assert _declared_in_pyproject(dist) or _declared_in_requirements(dist), (
        f"{dist} (imported as `{import_name}`) is used by production code but is "
        f"declared in neither pyproject.toml nor requirements.txt"
    )


def test_primer3_is_a_hard_dependency_not_only_an_extra():
    """primer3 backs the cross-ligation + orthogonality screens — not optional."""
    deps_block = PYPROJECT.split("[project.optional-dependencies]")[0]
    assert "primer3-py" in deps_block, (
        "primer3-py must be in [project].dependencies: qc/cross_ligation.py "
        "imports it on the production path, unguarded by a has_* check"
    )


def test_thermo_libraries_are_lower_bounded():
    """Pin lower bounds: Biopython changed its default saltcorr (5 -> 0)."""
    for dist in ("biopython", "primer3-py"):
        line = next(
            (ln for ln in REQUIREMENTS.splitlines()
             if re.match(rf"{re.escape(dist)}\b", ln.strip(), re.IGNORECASE)),
            "",
        )
        assert ">=" in line, f"{dist} needs a lower bound in requirements.txt, got {line!r}"
