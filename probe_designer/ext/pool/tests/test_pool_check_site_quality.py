"""``probe-design pool check`` names pooled sites that are proven false positives.

The set-based table (ledger ``site_quality/registry.tsv``) is resolved per
(site, chemistry). The check reads it in each probe's OWN chemistry, warns on an
implicated site or a site the chemistries disagree on, and never treats a
candidate -- a site inside an unresolved false-positive union -- as bad.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("primer3")
pytest.importorskip("probe_book")

from probe_designer.cli.app import app

runner = CliRunner()

SITE_A, SITE_B = "SITE_000000000001", "SITE_000000000002"
TABLE_HEADER = ("site_id\tchemistry\tgene_symbol\tfp_status\tfp_proof\tcandidate_set_size\t"
                "n_pools_clean\tn_pools_fp\tfn_status\tdiffers_by_chemistry\tfp_conflict\t"
                "probe_ids\tlast_recomputed\tinputs_hash\n")


def _pool(root: Path, *, table: str | None) -> Path:
    (root / "probes").mkdir()
    (root / "probes" / "registry.tsv").write_text(
        "probe_id\tcodebook\tchemistry\tprobe_arm5\tprobe_arm3\tsequence\ttarget\tsite_id\n"
        "SP369_AAAA\tSP369\tdRNA\tGCATAGCAGCAGCAGCATAG\tTGTGTGTGCACGCACGCATG\t"
        "GCATAGCAGCAGCAGCATAGTCCCTACACGACGCTCTTCCGTGTGTGTGCACGCACGCATG\tGENEA\t"
        f"{SITE_A}\n"
        "SP369_BBBB\tSP369\tcDNA\tAAGAAGAAGAAGAAGAAGAA\tTTCTTCTTCTTCTTCTTCTT\t"
        "AAGAAGAAGAAGAAGAAGAATCCCTACACGACGCTCTTCCGTTCTTCTTCTTCTTCTTCTT\tGENEB\t"
        f"{SITE_B}\n",
        encoding="utf-8",
    )
    pool_dir = root / "pools" / "test_pool_v1"
    pool_dir.mkdir(parents=True)
    (pool_dir / "manifest.yaml").write_text(
        "pool_id: test_pool_v1\npool_name: test\ncodebook: SP369\n"
        "target_well_volume_uL: 10.0\nrecipe:\n"
        "- order_id: ORD1\n  probe_id: SP369_AAAA\n  concentration_M: 1.0e-08\n"
        "- order_id: ORD1\n  probe_id: SP369_BBBB\n  concentration_M: 1.0e-08\n",
        encoding="utf-8",
    )
    if table is not None:
        (root / "site_quality").mkdir()
        (root / "site_quality" / "registry.tsv").write_text(TABLE_HEADER + table,
                                                          encoding="utf-8")
    return root


def _warnings(output: str) -> list[str]:
    return [line for line in output.splitlines() if "WARNING" in line]


def test_an_implicated_site_in_the_probes_chemistry_is_named(tmp_path):
    root = _pool(tmp_path, table=(
        f"{SITE_A}\tdRNA\tGENEA\timplicated\tFP@P_a_v1[1]\t1\t0\t1\tunknown\tfalse\t\t"
        f"SP369_AAAA\t2026-09-10\tabc\n"
        f"{SITE_B}\tcDNA\tGENEB\tcandidate\tFP@P_a_v1[3 left]\t3\t0\t1\tunknown\tfalse\t\t"
        f"SP369_BBBB\t2026-09-10\tabc\n"))
    result = runner.invoke(app, ["pool", "check", "test_pool_v1", "--repo-root", str(root)])
    assert result.exit_code == 0, result.output
    warnings = _warnings(result.output)
    assert any("SP369_AAAA" in w and "implicated" in w for w in warnings)
    # a candidate is inside an unresolved union: reported, never warned about
    assert not any("SP369_BBBB" in w for w in warnings)
    report = root / "pools" / "test_pool_v1" / "pool_check" / "site_quality.tsv"
    assert "candidate" in report.read_text(encoding="utf-8")


def test_only_the_probes_own_chemistry_decides_its_status(tmp_path):
    """SITE_A is implicated in cDNA but clean in dRNA; the pooled probe is dRNA."""
    root = _pool(tmp_path, table=(
        f"{SITE_A}\tcDNA\tGENEA\timplicated\tFP@P_c_v1[1]\t1\t0\t1\tunknown\ttrue\t\t"
        f"SP369_CCCC\t2026-09-10\tabc\n"
        f"{SITE_A}\tdRNA\tGENEA\tclean\tclean@P_d_v1\t\t1\t0\tsensitive\ttrue\t\t"
        f"SP369_AAAA\t2026-09-10\tabc\n"))
    result = runner.invoke(app, ["pool", "check", "test_pool_v1", "--repo-root", str(root)])
    assert result.exit_code == 0, result.output
    (warning,) = _warnings(result.output)
    assert "SP369_AAAA" in warning and "clean here" in warning


def test_a_site_the_table_does_not_list_is_untested(tmp_path):
    root = _pool(tmp_path, table="")
    result = runner.invoke(app, ["pool", "check", "test_pool_v1", "--repo-root", str(root)])
    assert result.exit_code == 0, result.output
    assert "'untested': 2" in result.output
    assert not _warnings(result.output)


def test_a_missing_table_is_never_read_as_clean(tmp_path):
    root = _pool(tmp_path, table=None)
    result = runner.invoke(app, ["pool", "check", "test_pool_v1", "--repo-root", str(root)])
    assert result.exit_code == 0, result.output
    assert "site_quality: not checked" in result.output
