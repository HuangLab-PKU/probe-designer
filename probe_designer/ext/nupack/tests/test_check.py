"""Tests for ``probe_designer.ext.nupack.check`` — v2.3 NUPACK ternary screen."""
from __future__ import annotations

from pathlib import Path

import pytest

# Skip the whole module if NUPACK 4 is not installed
pytest.importorskip("nupack")


from probe_designer.qc.cross_ligation import ProbeForScreen  # noqa: E402
from probe_designer.ext.nupack.check import (  # noqa: E402
    TernaryHit,
    screen_ternary,
    write_ternary_report,
    _parse_dotbracket_pairs,
    _check_mfe_nick_geometry,
    _check_mfe_vicinity_contiguous,
)


# ----------------------------------------------------------------------
# Fixtures — XLIG_A + XLIG_B where B.arm5 = RC(arm3_A + arm5_A)
# (same convention as test_qc_cross_ligation.py)
# ----------------------------------------------------------------------


BB = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"

XLIG_A = ProbeForScreen(
    probe_id="xlig_A", chemistry="dRNA",
    probe_arm5="CCCAATTGCGCAATATCATG",
    probe_arm3="GGGGGAAAATTTCCCAAGGG",
    sequence=("CCCAATTGCGCAATATCATG" + BB + "GGGGGAAAATTTCCCAAGGG").upper(),
    target="GX",
)
XLIG_B = ProbeForScreen(
    probe_id="xlig_B", chemistry="dRNA",
    probe_arm5="CATGATATTGCGCAATTGGGCCCTTGGGAAATTTTCCCCC",
    probe_arm3="AAGCTTAACTGGCCATAAGT",
    sequence=("CATGATATTGCGCAATTGGGCCCTTGGGAAATTTTCCCCC" + BB + "AAGCTTAACTGGCCATAAGT").upper(),
    target="GY",
)


# ----------------------------------------------------------------------
# Dot-bracket parser
# ----------------------------------------------------------------------


def test_parse_dotbracket_simple_pair():
    """Simple 4-mer hairpin: ``((..))`` — pos 0↔5, 1↔4, 2/3 unpaired."""
    pairs = _parse_dotbracket_pairs("((..))")
    assert pairs[0] == 5 and pairs[5] == 0
    assert pairs[1] == 4 and pairs[4] == 1
    assert 2 not in pairs and 3 not in pairs


def test_parse_dotbracket_multi_strand_separator():
    """'+' separators must NOT advance the position counter."""
    # Three strands of 2 chars each, all unpaired: ``..+..+..``
    pairs = _parse_dotbracket_pairs("..+..+..")
    assert pairs == {}
    # Pair across strands: ``((+))`` — pos 0↔3, 1↔2
    pairs = _parse_dotbracket_pairs("((+))")
    assert pairs[0] == 3 and pairs[1] == 2


# ----------------------------------------------------------------------
# MFE geometry helpers
# ----------------------------------------------------------------------


def test_check_mfe_nick_geometry_clean():
    """In a 3-strand concat (arm3 | splint | arm5) with arm3 last base pairing
    to splint pos 0 and arm5 first base pairing to splint pos -1 (i.e. b_3oh=0,
    b_5p=-1, diff=+1), the helper should return can_ligate=True."""
    arm3_len, splint_len = 20, 100
    # arm3 last base = concat pos 19. arm5 first base = concat pos 120.
    # We pair them inside splint such that b_3oh = 1, b_5p = 0 (splint coords).
    pairs = {
        19: arm3_len + 1,           # arm3 last → splint pos 1
        arm3_len + 1: 19,
        120: arm3_len + 0,          # arm5 first → splint pos 0
        arm3_len + 0: 120,
    }
    can_ligate, b3, b5 = _check_mfe_nick_geometry(pairs, arm3_len, splint_len)
    assert can_ligate
    assert b3 == 1 and b5 == 0


def test_check_mfe_nick_geometry_non_adjacent_fails():
    arm3_len, splint_len = 20, 100
    pairs = {
        19: arm3_len + 50,           # arm3 last → splint pos 50
        arm3_len + 50: 19,
        120: arm3_len + 30,          # arm5 first → splint pos 30 (not adjacent)
        arm3_len + 30: 120,
    }
    can_ligate, _, _ = _check_mfe_nick_geometry(pairs, arm3_len, splint_len)
    assert not can_ligate


def test_check_mfe_vicinity_contiguous_clean():
    """arm3 last 4 + arm5 first 4 all paired and antiparallel on splint."""
    arm3_len, splint_len, arm5_len = 20, 100, 20
    pairs: dict = {}
    # arm3 positions 16-19 → splint positions 4, 3, 2, 1 (descending)
    for i, sp in enumerate([4, 3, 2, 1]):
        pairs[arm3_len - 4 + i] = arm3_len + sp
        pairs[arm3_len + sp] = arm3_len - 4 + i
    # arm5 positions 0-3 → splint positions 0, -1, -2, -3? Need positive.
    # Use splint positions 0 down through -3? They have to be ≥ 0.
    # Redo with arm3 last 4 → splint 7,6,5,4; arm5 first 4 → splint 3,2,1,0.
    pairs.clear()
    for i, sp in enumerate([7, 6, 5, 4]):
        pairs[arm3_len - 4 + i] = arm3_len + sp
        pairs[arm3_len + sp] = arm3_len - 4 + i
    arm5_base = arm3_len + splint_len
    for i, sp in enumerate([3, 2, 1, 0]):
        pairs[arm5_base + i] = arm3_len + sp
        pairs[arm3_len + sp] = arm5_base + i
    ok = _check_mfe_vicinity_contiguous(pairs, arm3_len, splint_len, arm5_len, n_each_side=3)
    assert ok


def test_check_mfe_vicinity_disabled_when_n_zero():
    assert not _check_mfe_vicinity_contiguous({}, 20, 100, 20, n_each_side=0)


# ----------------------------------------------------------------------
# Live NUPACK runs
# ----------------------------------------------------------------------


def test_ternary_complex_enumerated():
    """NUPACK must enumerate the productive (arm3, splint, arm5) ternary
    complex for the XLIG fixture and report a non-zero equilibrium
    concentration (however small)."""
    hits = screen_ternary(
        [XLIG_A, XLIG_B],
        prefilter_pairs=[("xlig_A", "xlig_B")],
    )
    assert len(hits) == 1
    h = hits[0]
    assert h.seq_a_id == "xlig_A" and h.seq_b_id == "xlig_B"
    assert h.arm3_len > 0 and h.splint_len > 0 and h.arm5_len > 0
    assert h.ternary_concentration_m > 0


def test_coaxial_stacking_makes_ternary_more_stable():
    """ensemble='stacking' (default) produces a more negative ternary ΔG
    than ensemble='nostacking' — this is the coaxial-stacking contribution
    that v2.2 / primer3 binary cannot capture."""
    hits_stack = screen_ternary(
        [XLIG_A, XLIG_B], prefilter_pairs=[("xlig_A", "xlig_B")],
        ensemble="stacking",
    )
    hits_nostack = screen_ternary(
        [XLIG_A, XLIG_B], prefilter_pairs=[("xlig_A", "xlig_B")],
        ensemble="nostacking",
    )
    assert len(hits_stack) == 1 and len(hits_nostack) == 1
    dg_stack = hits_stack[0].ternary_dg_kcal
    dg_nostack = hits_nostack[0].ternary_dg_kcal
    # Stacking ensemble adds coaxial-stacking subensembles → more stable
    assert dg_stack < dg_nostack
    # Magnitude difference should be at least ~1 kcal/mol (typical coaxial range)
    assert (dg_nostack - dg_stack) >= 1.0


def test_screen_ternary_with_prefilter_pairs_subset():
    """Only specified prefilter_pairs are evaluated; others skipped."""
    # Three probes — without prefilter would evaluate 3 × 2 = 6 directions
    extra = ProbeForScreen(
        probe_id="extra", chemistry="dRNA",
        probe_arm5="AAAAAAAAAAAAAAAAAAAA",
        probe_arm3="TTTTTTTTTTTTTTTTTTTT",
        sequence=("AAAAAAAAAAAAAAAAAAAA" + BB + "TTTTTTTTTTTTTTTTTTTT").upper(),
        target="GZ",
    )
    hits = screen_ternary(
        [XLIG_A, XLIG_B, extra],
        prefilter_pairs=[("xlig_A", "xlig_B")],   # only one direction
    )
    assert len(hits) == 1
    assert hits[0].seq_a_id == "xlig_A" and hits[0].seq_b_id == "xlig_B"


def test_audit_fields_echo_user_params():
    """TernaryHit echoes the buffer parameters used, so audit TSVs are
    reproducible / inspectable."""
    hits = screen_ternary(
        [XLIG_A, XLIG_B], prefilter_pairs=[("xlig_A", "xlig_B")],
        sodium_m=0.1, magnesium_m=0.005, celsius=50.0,
        ensemble="stacking", strand_conc_m=2e-7, vicinity_n_each_side=4,
    )
    h = hits[0]
    assert h.sodium_m == 0.1
    assert h.magnesium_m == 0.005
    assert h.celsius_used == 50.0
    assert h.ensemble_used == "stacking"
    assert h.strand_conc_m == 2e-7
    assert h.vicinity_n_each_side == 4


# ----------------------------------------------------------------------
# Report writer
# ----------------------------------------------------------------------


def test_write_ternary_report_smoke(tmp_path: Path):
    hits = screen_ternary([XLIG_A, XLIG_B], prefilter_pairs=[("xlig_A", "xlig_B")])
    out = tmp_path / "ternary.tsv"
    write_ternary_report(hits, out)
    assert out.exists() and out.stat().st_size > 0
    text = out.read_text(encoding="utf-8")
    # v2.3 schema headers
    for col in (
        "seq_a_id", "ternary_dg_kcal", "ternary_concentration_m",
        "ternary_fraction_of_b", "mfe_nick_adjacent", "mfe_vicinity_contiguous",
        "ensemble_used", "celsius_used", "sodium_m", "magnesium_m",
        "mfe_dotbracket",
    ):
        assert col in text, f"missing column {col} in TSV header"


def test_write_ternary_report_empty(tmp_path: Path):
    out = tmp_path / "empty.tsv"
    write_ternary_report([], out)
    assert out.exists()
    text = out.read_text(encoding="utf-8")
    assert "seq_a_id" in text
    assert text.count("\n") == 1   # header only


# ----------------------------------------------------------------------
# Lazy-import error path
# ----------------------------------------------------------------------


def test_lazy_import_error_message(monkeypatch):
    """Hiding `nupack` from sys.modules + import-finders simulates a missing
    install; the error must mention NUPACK and the install path."""
    import sys
    import importlib

    # Save and remove nupack from sys.modules so re-import triggers the lazy path
    saved = {k: sys.modules.pop(k) for k in list(sys.modules) if k == "nupack" or k.startswith("nupack.")}

    # Block re-import by injecting None into sys.modules['nupack']
    sys.modules["nupack"] = None

    try:
        # Force a fresh import of check
        if "probe_designer.ext.nupack.check" in sys.modules:
            del sys.modules["probe_designer.ext.nupack.check"]
        check_mod = importlib.import_module("probe_designer.ext.nupack.check")
        with pytest.raises(ImportError) as excinfo:
            check_mod.screen_ternary(
                [XLIG_A, XLIG_B], prefilter_pairs=[("xlig_A", "xlig_B")],
            )
        msg = str(excinfo.value)
        assert "NUPACK" in msg or "nupack" in msg
    finally:
        # Restore
        sys.modules.pop("nupack", None)
        for k, v in saved.items():
            sys.modules[k] = v
        # Reload check to restore good state for subsequent tests
        if "probe_designer.ext.nupack.check" in sys.modules:
            importlib.reload(sys.modules["probe_designer.ext.nupack.check"])
