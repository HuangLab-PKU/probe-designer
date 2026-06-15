"""Tests for the v2 assemble-step cross-lig hook + the CLI flag plumbing.

Tests within-batch screening (the path with NO ``--target-pool``). External
pool integration (the path WITH ``--target-pool``) is tested in
``tests/ext/pool/test_loader.py`` and ``test_cli_check.py``.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

pytest.importorskip("primer3")


from probe_designer.cli.app import app
from probe_designer.qc.cross_lig_check import CandidateProbe
from probe_designer.qc.assemble_hook import (
    apply_cross_lig_check,
    probes_df_to_candidates,
)


runner = CliRunner()


_BB = "TCCCTACACGACGCTCTTCCG"   # 21-nt neutral bb


# v2-compatible cross-lig pair: A_rotated (arm3 + arm5, 40 nt) has its RC
# present contiguously inside B's arm5 (40 nt). This guarantees A-as-ligator
# on B-as-splint is ligation-competent under v2's junction-adjacency criterion.
_A_ARM3 = "GGGGGAAAATTTCCCAAGGG"     # last base = G is A's 3'-OH
_A_ARM5 = "CCCAATTGCGCAATATCATG"     # first base = C is A's 5'-P
_B_ARM5 = "CATGATATTGCGCAATTGGGCCCTTGGGAAATTTTCCCCC"   # = RC(arm3_A + arm5_A)

# A's full probe_sequence — arm5 + bb + arm3, all uppercase (assembled form).
A_SEQ = (_A_ARM5 + _BB + _A_ARM3).upper()
B_SEQ = (_B_ARM5 + _BB + "AAGCTTAACTGGCCATAAGT").upper()

V2_DIMER_PROBE_DICT_A = {
    "probe_name": "DIMER_A_100_dRNA_SP369_1",
    "chemistry": "dRNA",
    "probe_arm5": _A_ARM5,
    "probe_arm3": _A_ARM3,
    "probe_sequence": A_SEQ,
    "gene_name": "GENEX",
}
V2_DIMER_PROBE_DICT_B = {
    "probe_name": "DIMER_B_200_dRNA_SP369_1",
    "chemistry": "dRNA",
    "probe_arm5": _B_ARM5,
    "probe_arm3": "AAGCTTAACTGGCCATAAGT",
    "probe_sequence": B_SEQ,
    "gene_name": "GENEY",
}
SAFE_PROBE_DICT = {
    "probe_name": "SAFE_300_dRNA_SP369_1",
    "chemistry": "dRNA",
    "probe_arm5": "GCATAGCAGCAGCAGCATAG",
    "probe_arm3": "TGTGTGTGCACGCACGCATG",
    "probe_sequence": ("GCATAGCAGCAGCAGCATAG" + _BB + "TGTGTGTGCACGCACGCATG").upper(),
    "gene_name": "SAFEG",
}


# ----------------------------------------------------------------------
# Direct hook tests
# ----------------------------------------------------------------------


def test_probes_df_to_candidates_basic():
    df = pd.DataFrame([V2_DIMER_PROBE_DICT_A, SAFE_PROBE_DICT])
    cands = probes_df_to_candidates(df)
    assert len(cands) == 2
    assert all(isinstance(c, CandidateProbe) for c in cands)
    assert cands[0].probe_id == V2_DIMER_PROBE_DICT_A["probe_name"]
    assert cands[0].chemistry == "dRNA"
    assert cands[0].target == "GENEX"


def test_apply_no_hits_returns_input_unchanged_dropped_empty():
    df = pd.DataFrame([SAFE_PROBE_DICT, SAFE_PROBE_DICT.copy()])
    df.loc[1, "probe_name"] = "SAFE_400_dRNA_SP369_1"
    annotated, dropped, report = apply_cross_lig_check(df, reject=False)
    assert len(annotated) == 2
    assert len(dropped) == 0
    assert "cross_lig_partners" in annotated.columns
    assert (annotated["cross_lig_partners"] == "").all()


def test_apply_finds_within_batch_dimer_annotation_only():
    """v2-compatible pair: A_rotated has RC in B's arm5 contiguously → flag."""
    df = pd.DataFrame([V2_DIMER_PROBE_DICT_A, V2_DIMER_PROBE_DICT_B, SAFE_PROBE_DICT])
    annotated, dropped, report = apply_cross_lig_check(df, reject=False)
    assert len(annotated) == 3
    assert len(dropped) == 0
    assert len(report) >= 1
    flagged_rows = annotated[annotated["cross_lig_partners"] != ""]
    assert len(flagged_rows) >= 1
    safe_row = annotated[annotated["probe_name"] == SAFE_PROBE_DICT["probe_name"]]
    assert safe_row["cross_lig_partners"].iloc[0] == ""


def test_apply_reject_drops_dimer_pair():
    df = pd.DataFrame([V2_DIMER_PROBE_DICT_A, V2_DIMER_PROBE_DICT_B, SAFE_PROBE_DICT])
    annotated, dropped, report = apply_cross_lig_check(df, reject=True)
    safe_names = set(annotated["probe_name"])
    assert SAFE_PROBE_DICT["probe_name"] in safe_names
    # At least one of the dimer endpoints must be in dropped
    dropped_names = set(dropped["probe_name"])
    dimer_names = {V2_DIMER_PROBE_DICT_A["probe_name"], V2_DIMER_PROBE_DICT_B["probe_name"]}
    assert dropped_names & dimer_names, f"expected at least one dimer endpoint dropped; dropped={dropped_names}"


def test_apply_with_splint_pool():
    """Splint probes passed in directly (no bank loading); flagged candidate
    against pool member should mark POOL side."""
    from probe_designer.qc.cross_ligation import ProbeForScreen

    splint = ProbeForScreen(
        probe_id="pool_member_X",
        chemistry="dRNA",
        probe_arm5=_B_ARM5,
        probe_arm3="AAGCTTAACTGGCCATAAGT",
        sequence=B_SEQ,
        target="GENEX_POOL",
    )
    df = pd.DataFrame([V2_DIMER_PROBE_DICT_A, SAFE_PROBE_DICT])
    annotated, dropped, report = apply_cross_lig_check(
        df, reject=False, splint_probes=[splint],
    )
    # A should be flagged (against pool_member_X); SAFE stays clean.
    flagged = annotated[annotated["cross_lig_partners"] != ""]
    assert V2_DIMER_PROBE_DICT_A["probe_name"] in set(flagged["probe_name"])
    a_row = flagged[flagged["probe_name"] == V2_DIMER_PROBE_DICT_A["probe_name"]]
    assert "pool_member_X" in a_row["cross_lig_partners"].iloc[0]
    # Pool tag should appear
    assert ":POOL" in a_row["cross_lig_partners"].iloc[0]


# ----------------------------------------------------------------------
# CLI integration (within-batch only — no --target-pool)
# ----------------------------------------------------------------------


@pytest.fixture
def backbone_xlsx(tmp_path):
    p = tmp_path / "backbone_SP369.xlsx"
    pd.DataFrame({
        "No.": ["1"],
        "Sequence": [_BB],
    }).to_excel(p, index=False)
    return p


@pytest.fixture
def gene_info_xlsx(tmp_path):
    p = tmp_path / "gene_info.xlsx"
    pd.DataFrame({
        "gene_name": ["GENEX", "GENEY", "SAFEG"],
        "No.": ["1", "1", "1"],
    }).to_excel(p, index=False)
    return p


@pytest.fixture
def binding_sites_with_dimer(tmp_path):
    p = tmp_path / "sites.json"
    p.write_text(json.dumps({
        "GENEX": [{
            "gene_name": "GENEX",
            "arm_5prime": _A_ARM5, "arm_3prime": _A_ARM3,
            "st": 100, "en": 140, "g_content": 0.5,
            "tm": 70.0, "tm_3prime": 70.0, "tm_5prime": 70.0,
        }],
        "GENEY": [{
            "gene_name": "GENEY",
            "arm_5prime": _B_ARM5[:20], "arm_3prime": _B_ARM5[20:],  # split B's arm5 to fit assemble
            "st": 200, "en": 240, "g_content": 0.5,
            "tm": 70.0, "tm_3prime": 70.0, "tm_5prime": 70.0,
        }],
        "SAFEG": [{
            "gene_name": "SAFEG",
            "arm_5prime": SAFE_PROBE_DICT["probe_arm5"],
            "arm_3prime": SAFE_PROBE_DICT["probe_arm3"],
            "st": 300, "en": 340, "g_content": 0.5,
            "tm": 60.0, "tm_3prime": 60.0, "tm_5prime": 60.0,
        }],
    }), encoding="utf-8")
    return p


def test_cli_assemble_without_xlig_flag_unchanged(
    tmp_path, binding_sites_with_dimer, gene_info_xlsx, backbone_xlsx,
):
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "assemble",
        "--binding-sites", str(binding_sites_with_dimer),
        "--gene-info", str(gene_info_xlsx),
        "--backbone", str(backbone_xlsx),
        "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    assert not (out / "cross_lig_report.tsv").exists()
    assert not (out / "dropped_cross_lig.tsv").exists()


def test_cli_assemble_with_within_batch_screen_writes_report(
    tmp_path, binding_sites_with_dimer, gene_info_xlsx, backbone_xlsx,
):
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "assemble",
        "--binding-sites", str(binding_sites_with_dimer),
        "--gene-info", str(gene_info_xlsx),
        "--backbone", str(backbone_xlsx),
        "--output", str(out),
        "--check-cross-lig",
    ])
    assert result.exit_code == 0, result.output
    report = out / "cross_lig_report.tsv"
    assert report.exists() and report.stat().st_size > 0
