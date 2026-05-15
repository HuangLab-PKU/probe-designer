"""Tests for ``probe-design migrate-columns``.

The migration upgrades a legacy probe Excel file to Schema v2:
  * renames columns per ``LEGACY_RENAMES``
  * folds the standalone ``iLock`` boolean column into ``chemistry``
  * normalizes legacy chemistry tokens (``mRNA_noiLock``/``mRNA`` -> ``dRNA``)
  * regenerates ``probe_name`` against a user-supplied codebook
  * applies the v2 column order
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from probe_designer.cli.app import app
from probe_designer.io.probe_schema import (
    FINAL_PADLOCK_COLUMNS,
    FINAL_RT_PRIMER_COLUMNS,
)


runner = CliRunner()


# ---------------------------------------------------------------------------
# Fixture: a tiny mutation-style legacy final-probes file
# ---------------------------------------------------------------------------


@pytest.fixture
def legacy_mutation_xlsx(tmp_path):
    """One iLock + one mRNA_noiLock + one cDNA row, in the v0/v1 layout."""
    fp = tmp_path / "BZ_mut_final_probes.xlsx"
    pd.DataFrame([
        {
            "No.": 100,
            "Probe_Name": "BRCA1_iLock_12345_SP_100",
            "target_type": "mut", "chemistry": "iLock", "iLock": "yes",
            "gene_name": "BRCA1", "clone_id": "",
            "position": 12345,
            "Chr": "chr17", "Start": 12345, "End": 12345,
            "Ref": "A", "Alt": "G", "Strand": "+",
            "Mutation_Type": "SNV", "tm_model": "R_DNA_NN1",
            "mutation_position": "5end_and_3tip",
            "arm_5prime": "ACGTACGTACGTACGTACGT",
            "arm_3prime": "TGCATGCATGCATGCATGCA",
            "arm5_len": 20, "arm3_len": 20, "plp_pos5_local": 0,
            "iLock_flap": "CGTTGCTGTGGCG", "iLock_linker_nt": "G",
            "backbone_No.": 100, "backbone_name": "BB1",
            "backbone_sequence": "GCATGCATGCATGCATGCATGCATGCATGCATGCAT",
            "Probe_Seq": "ACGT" * 25, "Probe_Length": 100,
            "tm": 60.0, "tm_5prime": 59.0, "tm_3prime": 61.0,
            "g_content": 0.5, "free_energy": -2.0, "score": 5.5,
            "target_sequence": "ACGTACGT" * 5,
            "blast_hits_count": 1, "same_gene_exact": 1,
            "backbone_file": "backbone_SP100.xlsx",
        },
        {
            "No.": 101,
            "Probe_Name": "BRCA1_mRNA_noiLock_12345_SP_101",
            "target_type": "mut", "chemistry": "mRNA_noiLock", "iLock": "no",
            "gene_name": "BRCA1", "clone_id": "",
            "position": 12345,
            "Chr": "chr17", "Start": 12345, "End": 12345,
            "Ref": "A", "Alt": "G", "Strand": "+",
            "Mutation_Type": "SNV", "tm_model": "R_DNA_NN1",
            "mutation_position": "3end",
            "arm_5prime": "ACGTACGTACGTACGTACGT",
            "arm_3prime": "TGCATGCATGCATGCATGCA",
            "arm5_len": 20, "arm3_len": 20, "plp_pos5_local": 0,
            "iLock_flap": "", "iLock_linker_nt": "",
            "backbone_No.": 101, "backbone_name": "BB2",
            "backbone_sequence": "GCATGCATGCATGCATGCATGCATGCATGCATGCAT",
            "Probe_Seq": "ACGT" * 25, "Probe_Length": 100,
            "tm": 60.0, "tm_5prime": 59.0, "tm_3prime": 61.0,
            "g_content": 0.5, "free_energy": -2.0, "score": 5.5,
            "target_sequence": "ACGTACGT" * 5,
            "blast_hits_count": 1, "same_gene_exact": 1,
            "backbone_file": "backbone_SP100.xlsx",
        },
    ]).to_excel(fp, index=False)
    return fp


@pytest.fixture
def legacy_rt_primer_xlsx(tmp_path):
    fp = tmp_path / "BZ_mut_cDNA_RT_primers.xlsx"
    pd.DataFrame([
        {
            "No.": 102, "Gene": "BRCA1", "Chr": "chr17",
            "Mut_Start": 12345, "Mut_End": 12345,
            "Ref": "A", "Alt": "G", "Strand": "+",
            "Mutation_Type": "SNV",
            "Padlock_Probe_Name": "BRCA1_cDNA_12345_SP_102",
            "RT_Primer_Name": "BRCA1_RT_12345_SP_102",
            "RT_Primer_Sequence": "ACGTACGTACGTACGTACGTACGT",
            "RT_Primer_Length": 24, "Tm": 55.0, "GC_Percent": 0.5,
            "Gap_nt": 15,
            "Primer_Genomic_Start": 12360, "Primer_Genomic_End": 12384,
            "Notes": "",
        },
    ]).to_excel(fp, index=False)
    return fp


# ---------------------------------------------------------------------------
# CLI behaviour
# ---------------------------------------------------------------------------


class TestMigrateColumns:
    def test_padlock_file_columns_match_schema_after_migration(
        self, legacy_mutation_xlsx,
    ):
        result = runner.invoke(app, [
            "migrate-columns", str(legacy_mutation_xlsx),
            "--codebook", "SP369",
        ])
        assert result.exit_code == 0, result.output

        df = pd.read_excel(legacy_mutation_xlsx)
        # Schema-v2 leading three columns + every schema column present.
        assert df.columns.tolist()[:3] == ["order", "probe_name", "probe_sequence"]
        assert set(FINAL_PADLOCK_COLUMNS).issubset(df.columns)
        # legacy columns gone. Schema-v2 keeps g_content as a distinct
        # legacy-reporting column (G-only fraction), so it is NOT in this list.
        for legacy in ("Probe_Name", "Probe_Seq", "arm_5prime",
                        "iLock", "Chr", "Start", "backbone_No."):
            assert legacy not in df.columns

    def test_padlock_chemistry_values_normalized(self, legacy_mutation_xlsx):
        runner.invoke(app, [
            "migrate-columns", str(legacy_mutation_xlsx),
            "--codebook", "SP369",
        ])
        df = pd.read_excel(legacy_mutation_xlsx)
        # mRNA_noiLock should be folded into dRNA; iLock survives because
        # the iLock column had "yes" for it.
        assert set(df["chemistry"]) == {"iLock", "dRNA"}

    def test_padlock_probe_names_regenerated(self, legacy_mutation_xlsx):
        runner.invoke(app, [
            "migrate-columns", str(legacy_mutation_xlsx),
            "--codebook", "SP369",
        ])
        df = pd.read_excel(legacy_mutation_xlsx)
        for _, row in df.iterrows():
            # New format: {gene}_{position}_{chem}_{codebook}_{No.}
            parts = row["probe_name"].split("_")
            assert parts[-2] == "SP369"
            assert parts[-1] == str(row["No."])
            assert parts[-3] in {"dRNA", "cDNA", "iLock"}

    def test_bak_file_written(self, legacy_mutation_xlsx):
        runner.invoke(app, [
            "migrate-columns", str(legacy_mutation_xlsx),
            "--codebook", "SP369",
        ])
        bak = legacy_mutation_xlsx.with_suffix(".xlsx.bak")
        assert bak.exists()

    def test_padlock_requires_codebook(self, legacy_mutation_xlsx):
        result = runner.invoke(app, [
            "migrate-columns", str(legacy_mutation_xlsx),
        ])
        assert result.exit_code != 0
        assert "codebook" in result.output.lower()

    def test_dry_run_writes_nothing(self, legacy_mutation_xlsx):
        before = legacy_mutation_xlsx.read_bytes()
        result = runner.invoke(app, [
            "migrate-columns", str(legacy_mutation_xlsx),
            "--codebook", "SP369", "--dry-run",
        ])
        assert result.exit_code == 0
        assert legacy_mutation_xlsx.read_bytes() == before
        # bak shouldn't be created in dry-run either
        assert not legacy_mutation_xlsx.with_suffix(".xlsx.bak").exists()

    def test_already_migrated_file_is_idempotent(self, legacy_mutation_xlsx):
        # First migration
        runner.invoke(app, [
            "migrate-columns", str(legacy_mutation_xlsx),
            "--codebook", "SP369",
        ])
        first = pd.read_excel(legacy_mutation_xlsx)
        # Second migration on the now-canonical file
        result = runner.invoke(app, [
            "migrate-columns", str(legacy_mutation_xlsx),
            "--codebook", "SP369",
        ])
        assert result.exit_code == 0
        second = pd.read_excel(legacy_mutation_xlsx)
        # Probe names should be identical (already-v2 names round-trip).
        assert first["probe_name"].tolist() == second["probe_name"].tolist()
        assert first.columns.tolist() == second.columns.tolist()

    def test_rt_primer_file_columns_match_schema(self, legacy_rt_primer_xlsx):
        result = runner.invoke(app, [
            "migrate-columns", str(legacy_rt_primer_xlsx),
        ])
        assert result.exit_code == 0, result.output
        df = pd.read_excel(legacy_rt_primer_xlsx)
        assert df.columns.tolist()[:3] == ["order", "probe_name", "probe_sequence"]
        assert set(FINAL_RT_PRIMER_COLUMNS).issubset(df.columns)
        assert df["chemistry"].iloc[0] == "RT_primer"
        assert df["probe_name"].iloc[0].endswith("_RT_primer")
        # Tm column renamed
        assert "tm" in df.columns
        assert "Tm" not in df.columns
        # GC_Percent renamed
        assert "gc_content" in df.columns
        assert "GC_Percent" not in df.columns

    def test_directory_walk_picks_up_padlock_and_rt(
        self, tmp_path, legacy_mutation_xlsx, legacy_rt_primer_xlsx,
    ):
        # Both files already live under tmp_path. Run migration against the
        # directory rather than individual files.
        result = runner.invoke(app, [
            "migrate-columns", str(tmp_path),
            "--codebook", "SP369",
        ])
        assert result.exit_code == 0, result.output
        # Both files migrated.
        df_pad = pd.read_excel(legacy_mutation_xlsx)
        df_rt = pd.read_excel(legacy_rt_primer_xlsx)
        assert df_pad.columns.tolist()[:3] == ["order", "probe_name", "probe_sequence"]
        assert df_rt.columns.tolist()[:3] == ["order", "probe_name", "probe_sequence"]
