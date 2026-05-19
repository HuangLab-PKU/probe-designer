"""Tests for the assemble-step cross-lig hook + the CLI flag plumbing."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

pytest.importorskip("primer3")
pytest.importorskip("probe_book")


from probe_designer.cli.app import app
from probe_designer.qc.cross_lig_check import CandidateProbe
from probe_designer.qc.assemble_hook import (
    apply_cross_lig_check,
    probes_df_to_candidates,
)


runner = CliRunner()


# Two known-dimer arm pairs (Tm > 27°C, end_dg < -5 kcal/mol).
DIMER_PROBE_DICT_A = {
    "probe_name": "GCXLIG_A_100_dRNA_SP369_1",
    "chemistry": "dRNA",
    "probe_arm5": "GCAGCAGCAGCAGCAGCAGC",
    "probe_arm3": "TGCTGCTGCTGCTGCTGCTG",
    "probe_sequence": "GCAGCAGCAGCAGCAGCAGCTGCTGCTGCTGCTGCTGCTG",
    "gene_name": "GENEX",
}
DIMER_PROBE_DICT_B = {
    "probe_name": "GCXLIG_B_200_dRNA_SP369_1",
    "chemistry": "dRNA",
    "probe_arm5": "CAGCAGCAGCAGCAGCAGCA",
    "probe_arm3": "GCTGCTGCTGCTGCTGCTGC",
    "probe_sequence": "CAGCAGCAGCAGCAGCAGCAGCTGCTGCTGCTGCTGCTGC",
    "gene_name": "GENEY",
}
SAFE_PROBE_DICT = {
    "probe_name": "SAFE_300_dRNA_SP369_1",
    "chemistry": "dRNA",
    "probe_arm5": "GCATAGCAGCAGCAGCATAG",
    "probe_arm3": "TGTGTGTGCACGCACGCATG",
    "probe_sequence": "GCATAGCAGCAGCAGCATAGTGTGTGTGCACGCACGCATG",
    "gene_name": "SAFEG",
}


def test_probes_df_to_candidates_basic():
    df = pd.DataFrame([DIMER_PROBE_DICT_A, SAFE_PROBE_DICT])
    cands = probes_df_to_candidates(df)
    assert len(cands) == 2
    assert all(isinstance(c, CandidateProbe) for c in cands)
    assert cands[0].probe_id == DIMER_PROBE_DICT_A["probe_name"]
    assert cands[0].chemistry == "dRNA"
    assert cands[0].target == "GENEX"


def test_apply_no_hits_returns_input_unchanged_dropped_empty():
    df = pd.DataFrame([SAFE_PROBE_DICT, SAFE_PROBE_DICT.copy()])
    df.loc[1, "probe_name"] = "SAFE_400_dRNA_SP369_1"
    annotated, dropped, report = apply_cross_lig_check(df, reject=False)
    assert len(annotated) == 2
    assert len(dropped) == 0
    assert "cross_lig_partners" in annotated.columns
    # all empty (no partners)
    assert (annotated["cross_lig_partners"] == "").all()


def test_apply_finds_within_batch_dimer_annotation_only():
    df = pd.DataFrame([DIMER_PROBE_DICT_A, DIMER_PROBE_DICT_B, SAFE_PROBE_DICT])
    annotated, dropped, report = apply_cross_lig_check(df, reject=False)
    assert len(annotated) == 3
    assert len(dropped) == 0
    assert len(report) >= 1
    flagged_rows = annotated[annotated["cross_lig_partners"] != ""]
    assert len(flagged_rows) == 2   # both endpoints of the dimer are flagged
    safe_row = annotated[annotated["probe_name"] == SAFE_PROBE_DICT["probe_name"]]
    assert safe_row["cross_lig_partners"].iloc[0] == ""


def test_apply_reject_drops_dimer_pair():
    df = pd.DataFrame([DIMER_PROBE_DICT_A, DIMER_PROBE_DICT_B, SAFE_PROBE_DICT])
    annotated, dropped, report = apply_cross_lig_check(df, reject=True)
    # The safe probe survives; both dimer endpoints go to dropped
    assert len(annotated) == 1
    assert annotated["probe_name"].iloc[0] == SAFE_PROBE_DICT["probe_name"]
    assert len(dropped) == 2
    assert set(dropped["probe_name"]) == {
        DIMER_PROBE_DICT_A["probe_name"], DIMER_PROBE_DICT_B["probe_name"],
    }
    assert "cross_lig_partners" in dropped.columns


def test_apply_with_pool(tmp_path: Path):
    """Pool member cross-ligs with one new candidate; the candidate is flagged."""
    (tmp_path / "probes").mkdir()
    (tmp_path / "probes" / "registry.tsv").write_text(
        "probe_id\tcodebook\tchemistry\tprobe_arm5\tprobe_arm3\tsequence\ttarget\n"
        f"pool_member_X\tSP369\tdRNA\t{DIMER_PROBE_DICT_A['probe_arm5']}\t"
        f"{DIMER_PROBE_DICT_A['probe_arm3']}\t\tGENEX_POOL\n",
        encoding="utf-8",
    )
    pool_dir = tmp_path / "pools" / "test_pool_v1"
    pool_dir.mkdir(parents=True)
    (pool_dir / "manifest.yaml").write_text(
        "pool_id: test_pool_v1\n"
        "pool_name: test\n"
        "codebook: SP369\n"
        "target_well_volume_uL: 10.0\n"
        "recipe:\n"
        "- order_id: ORD1\n"
        "  probe_id: pool_member_X\n"
        "  concentration_M: 1.0e-08\n",
        encoding="utf-8",
    )

    # New candidate dimers with the pool member but not with itself.
    df = pd.DataFrame([DIMER_PROBE_DICT_B, SAFE_PROBE_DICT])
    annotated, dropped, report = apply_cross_lig_check(
        df, reject=False, pool_id="test_pool_v1", repo_root=tmp_path,
    )
    # B should be flagged (against pool_member_X); safe stays clean.
    flagged = annotated[annotated["cross_lig_partners"] != ""]
    assert DIMER_PROBE_DICT_B["probe_name"] in set(flagged["probe_name"])
    assert "pool_member_X" in flagged["cross_lig_partners"].iloc[0]


# -------- CLI integration --------


@pytest.fixture
def backbone_xlsx(tmp_path):
    p = tmp_path / "backbone_SP369.xlsx"
    pd.DataFrame({
        "No.": ["1"],
        "Sequence": ["GTTTTTGGTAAGCTTCGGATCCTCAGACGGAAGACTC"],
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
            "arm_5prime": DIMER_PROBE_DICT_A["probe_arm5"],
            "arm_3prime": DIMER_PROBE_DICT_A["probe_arm3"],
            "st": 100, "en": 140, "g_content": 0.5,
            "tm": 70.0, "tm_3prime": 70.0, "tm_5prime": 70.0,
        }],
        "GENEY": [{
            "gene_name": "GENEY",
            "arm_5prime": DIMER_PROBE_DICT_B["probe_arm5"],
            "arm_3prime": DIMER_PROBE_DICT_B["probe_arm3"],
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
    """No --target-pool/--reject-cross-lig: behavior identical to pre-hook."""
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "assemble",
        "--binding-sites", str(binding_sites_with_dimer),
        "--gene-info", str(gene_info_xlsx),
        "--backbone", str(backbone_xlsx),
        "--output", str(out),
    ])
    assert result.exit_code == 0, result.output
    # No report file written when no flag
    assert not (out / "cross_lig_report.tsv").exists()
    assert not (out / "dropped_cross_lig.tsv").exists()


def test_cli_assemble_with_within_batch_screen_writes_report(
    tmp_path, binding_sites_with_dimer, gene_info_xlsx, backbone_xlsx,
):
    """`--check-cross-lig` (no pool, no reject) emits an annotated xlsx + report."""
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


def test_cli_assemble_reject_drops_dimers(
    tmp_path, binding_sites_with_dimer, gene_info_xlsx, backbone_xlsx,
):
    """The CLI assigns its own probe_name based on (gene, position, …); verify
    the dimer pair makes it into dropped_cross_lig.tsv and the safe probe
    survives into the kept output xlsx."""
    out = tmp_path / "out"
    result = runner.invoke(app, [
        "assemble",
        "--binding-sites", str(binding_sites_with_dimer),
        "--gene-info", str(gene_info_xlsx),
        "--backbone", str(backbone_xlsx),
        "--output", str(out),
        "--check-cross-lig",
        "--reject-cross-lig",
    ])
    assert result.exit_code == 0, result.output
    dropped = out / "dropped_cross_lig.tsv"
    assert dropped.exists()
    drop_df = pd.read_csv(dropped, sep="\t")
    drop_genes = set(drop_df["gene_name"])
    # Both dimer endpoints (GENEX, GENEY) should be dropped; SAFEG should survive.
    assert {"GENEX", "GENEY"} <= drop_genes
    assert "SAFEG" not in drop_genes

    xlsx_files = list(out.glob("*.xlsx"))
    assert xlsx_files
    kept = pd.read_excel(xlsx_files[0])
    assert "SAFEG" in set(kept["gene_name"])
    assert "GENEX" not in set(kept["gene_name"])
