"""Static check: ``probe_designer.qc.*`` must NOT import heavy optional deps.

The qc/ layer is the bank-free + nupack-free compute layer. Heavy / optional
dependencies are confined to:

  * bank (``probe_book``) → ``probe_designer.ext.pool.*``
  * NUPACK 4 (``nupack``) → ``probe_designer.ext.nupack.*``

This test parses every .py under qc/ with ``ast`` and fails if any one imports
the forbidden modules (including lazy/function-scope imports — AST sees those
too).
"""
from __future__ import annotations

import ast
from pathlib import Path

import probe_designer.qc


# Modules forbidden inside probe_designer.qc.*
FORBIDDEN_PREFIXES: tuple[str, ...] = ("probe_book", "nupack")

# Where each forbidden module's calling code is allowed to live
FORBIDDEN_REROUTE = {
    "probe_book": "probe_designer.ext.pool.*",
    "nupack": "probe_designer.ext.nupack.*",
}


def _qc_directory() -> Path:
    return Path(probe_designer.qc.__file__).parent


def _iter_python_files(root: Path):
    for path in root.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        yield path


def _find_forbidden_imports(path: Path) -> list[str]:
    """Return a list of offending import statements (formatted) in this file."""
    source = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(source, filename=str(path))
    except SyntaxError as exc:
        raise AssertionError(f"could not parse {path}: {exc}") from exc

    def _is_forbidden(module_name: str) -> bool:
        for prefix in FORBIDDEN_PREFIXES:
            if module_name == prefix or module_name.startswith(prefix + "."):
                return True
        return False

    offenders: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if _is_forbidden(alias.name):
                    offenders.append(f"line {node.lineno}: import {alias.name}")
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if _is_forbidden(mod):
                names = ", ".join(a.name for a in node.names)
                offenders.append(f"line {node.lineno}: from {mod} import {names}")
    return offenders


def test_qc_modules_do_not_import_forbidden():
    """Walk every .py file under probe_designer/qc/ and assert none of them
    statically import the forbidden heavy/optional dependencies (probe_book,
    nupack). AST detects lazy/function-scope imports too.
    """
    qc_dir = _qc_directory()
    violations: dict[Path, list[str]] = {}
    for py in _iter_python_files(qc_dir):
        offenders = _find_forbidden_imports(py)
        if offenders:
            violations[py] = offenders

    if violations:
        msg = [
            f"qc/ must NOT import {', '.join(FORBIDDEN_PREFIXES)} "
            "(bank-free + nupack-free invariant); found:"
        ]
        for path, hits in violations.items():
            msg.append(f"  {path.relative_to(qc_dir)}:")
            for h in hits:
                msg.append(f"    {h}")
        msg.append("")
        msg.append("Move forbidden-dep calls to:")
        for k, v in FORBIDDEN_REROUTE.items():
            msg.append(f"  {k:12s} → {v}")
        raise AssertionError("\n".join(msg))


# Back-compat alias — the test was previously named test_qc_modules_do_not_import_probe_book.
test_qc_modules_do_not_import_probe_book = test_qc_modules_do_not_import_forbidden
