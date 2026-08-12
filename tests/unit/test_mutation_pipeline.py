"""Smoke test for run_mutation_pipeline.

Runs the full pipeline (with --skip-blast) on BZ07's single TP53 mutation,
which is the smallest known-good real mutation panel in the repo. The test
copies the input xlsx into tmp_path so the original is not mutated.

Skips when GRCh38.fa or the SPRINTseq backbone is not present.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest

from probe_designer.ext.mutation.config import MutationConfig
from probe_designer.ext.mutation.pipeline import run_mutation_pipeline


REPO_ROOT = Path(__file__).resolve().parents[3]
GENOME = REPO_ROOT / "data" / "genome" / "GRCh38.fa"
# The backbone moved out of the experiments root when that tree was emptied
# ("every loose file moves to its own layer"); the bank is now its only home.
BACKBONE = REPO_ROOT / "backbones" / "SP369_add_A_v1" / "sequences.xlsx"
BZ07_INPUT = (REPO_ROOT / "experiments" / "20260425_ZCH_BZ07_mut"
              / "input" / "BZ07_mutations.xlsx")


def _require_assets():
    if not GENOME.exists():
        pytest.skip(f"GRCh38 not present at {GENOME}")
    if not BACKBONE.exists():
        pytest.skip(f"Backbone not present at {BACKBONE}")
    if not BZ07_INPUT.exists():
        pytest.skip(f"BZ07 input not present at {BZ07_INPUT}")


def test_pipeline_smoke_bz07_skip_blast(tmp_path: Path):
    """End-to-end run with the full chemistry triplet, no real BLAST round."""
    _require_assets()

    # Copy input into tmp so we don't mutate the real BZ07 file
    fake_experiment = tmp_path / "20260510_ZCH_BZ07_smoke"
    (fake_experiment / "input").mkdir(parents=True)
    (fake_experiment / "output").mkdir()
    input_xlsx = fake_experiment / "input" / "BZ07_mutations.xlsx"
    shutil.copy2(BZ07_INPUT, input_xlsx)

    cfg = MutationConfig(
        input_xlsx=input_xlsx,
        output_dir=fake_experiment / "output",
        backbone_file=BACKBONE,
        chemistries=["iLock", "dRNA", "cDNA"],
        start_no=200,
        skip_blast=True,
        genome_path=GENOME,
        codebook="SP369",
    )

    result = run_mutation_pipeline(cfg)

    assert result.mutations_total == 1
    assert result.mutations_mappable == 1
    assert result.probes_total == 3  # one per chemistry
    assert result.last_no == 202     # 200 + 3 - 1
    assert result.final_probes_xlsx.exists()
    assert result.final_probes_fasta.exists()
    assert result.rt_primers_xlsx is not None
    assert result.rt_primers_xlsx.exists()

    # Verify the final xlsx has 3 rows in the canonical column order
    df = pd.read_excel(result.final_probes_xlsx, sheet_name="Final_Probes")
    assert len(df) == 3
    assert set(df["chemistry"]) == {"iLock", "dRNA", "cDNA"}
    # Schema-v2 leading three columns
    assert df.columns.tolist()[:3] == ["order", "probe_name", "probe_sequence"]
    # Numbering must be contiguous
    assert sorted(df["No."].tolist()) == [200, 201, 202]
    # iLock probe must start with the default flap
    ilock_row = df[df["chemistry"] == "iLock"].iloc[0]
    assert ilock_row["probe_sequence"].startswith("CGTTGCTGTGGCG")
    # probe_name format: {gene}_{pos}_{chem}_{codebook}_{No.}
    for _, row in df.iterrows():
        parts = row["probe_name"].rsplit("_", 4)
        assert parts[2] in {"iLock", "dRNA", "cDNA"}
        assert parts[3] == "SP369"

    # last_no.txt was written
    assert (cfg.output_dir / "last_no.txt").read_text().strip() == "202"


def test_pipeline_subset_chemistries_just_ilock(tmp_path: Path):
    """When only iLock runs, no RT primers + only one probe per mutation."""
    _require_assets()

    fake_experiment = tmp_path / "20260510_ZCH_BZ07_ilock_only"
    (fake_experiment / "input").mkdir(parents=True)
    (fake_experiment / "output").mkdir()
    input_xlsx = fake_experiment / "input" / "BZ07_mutations.xlsx"
    shutil.copy2(BZ07_INPUT, input_xlsx)

    cfg = MutationConfig(
        input_xlsx=input_xlsx,
        output_dir=fake_experiment / "output",
        backbone_file=BACKBONE,
        chemistries=["iLock"],
        start_no=300,
        skip_blast=True,
        genome_path=GENOME,
        codebook="SP369",
    )
    result = run_mutation_pipeline(cfg)

    assert result.probes_total == 1
    assert result.last_no == 300
    assert result.rt_primers_xlsx is None  # no cDNA -> no RT primers
