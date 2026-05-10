"""Smoke test for run_tcr_pipeline.

Runs the four-phase TCR pipeline end-to-end on the real BZ23_TCR input with
--skip-blast. Skips when the input + backbone files are not present.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest

from probe_designer.ext.tcr.config import TcrConfig
from probe_designer.ext.tcr.pipeline import run_tcr_pipeline


REPO_ROOT = Path(__file__).resolve().parents[3]
BACKBONE = REPO_ROOT / "experiments" / "SPRINTseq_369_add_A_backbone.xlsx"
BZ23_TCR_INPUT = (
    REPO_ROOT / "experiments" / "20260418_ZCH_BZ23_TCR"
    / "input" / "BZ23_tcr_clones.xlsx"
)


def _require_assets():
    if not BACKBONE.exists():
        pytest.skip(f"Backbone not present at {BACKBONE}")
    if not BZ23_TCR_INPUT.exists():
        pytest.skip(f"BZ23_TCR input not present at {BZ23_TCR_INPUT}")


def test_pipeline_smoke_bz23_mrna_only_skip_blast(tmp_path: Path):
    """Default chemistry (mRNA only) — no RT primers, no BLAST round."""
    _require_assets()

    fake_experiment = tmp_path / "20260510_ZCH_BZ23_TCR_smoke"
    (fake_experiment / "input").mkdir(parents=True)
    (fake_experiment / "output").mkdir()
    input_xlsx = fake_experiment / "input" / "BZ23_tcr_clones.xlsx"
    shutil.copy2(BZ23_TCR_INPUT, input_xlsx)

    cfg = TcrConfig(
        input_xlsx=input_xlsx,
        output_dir=fake_experiment / "output",
        backbone_file=BACKBONE,
        chemistries=["mRNA"],
        start_no=200,
        skip_blast=True,
        plot_tm_landscape=False,  # speed up the test; landscape covered separately
    )
    result = run_tcr_pipeline(cfg)

    assert result.clones_total >= 1
    assert result.clones_with_sites >= 1
    assert result.sites_selected >= result.clones_with_sites
    # mRNA only ⇒ probes_total == sites_selected
    assert result.probes_total == result.sites_selected
    assert result.rt_primers_xlsx is None  # no cDNA → no RT primers
    assert result.final_probes_xlsx.exists()
    assert result.final_probes_fasta.exists()

    # Verify final xlsx
    df = pd.read_excel(result.final_probes_xlsx)
    assert set(df["chemistry"]) == {"mRNA"}
    assert df["No."].tolist() == sorted(df["No."].tolist())
    assert df["No."].iloc[0] == 200
    # Probe_Name follows clone-id convention
    for name in df["Probe_Name"]:
        assert "_mRNA_" in name
        assert "_SP_" in name


def test_pipeline_smoke_bz23_dual_chemistry(tmp_path: Path):
    """--chemistries mRNA,cDNA → 2 probes per site + RT primers auto-run."""
    _require_assets()

    fake_experiment = tmp_path / "20260510_ZCH_BZ23_TCR_dual"
    (fake_experiment / "input").mkdir(parents=True)
    (fake_experiment / "output").mkdir()
    input_xlsx = fake_experiment / "input" / "BZ23_tcr_clones.xlsx"
    shutil.copy2(BZ23_TCR_INPUT, input_xlsx)

    cfg = TcrConfig(
        input_xlsx=input_xlsx,
        output_dir=fake_experiment / "output",
        backbone_file=BACKBONE,
        chemistries=["mRNA", "cDNA"],
        start_no=200,
        skip_blast=True,
        plot_tm_landscape=False,
    )
    result = run_tcr_pipeline(cfg)

    assert result.probes_total == result.sites_selected * 2  # 2 chemistries
    assert result.rt_primers_xlsx is not None
    assert result.rt_primers_xlsx.exists()

    df = pd.read_excel(result.final_probes_xlsx)
    assert set(df["chemistry"]) == {"mRNA", "cDNA"}
