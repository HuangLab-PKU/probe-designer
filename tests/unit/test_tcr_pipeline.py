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
# The ledger became its own repo on 2026-08-16, so backbones/ sits one level
# down while experiments/ stayed put. Probed inline rather than imported from
# probe_book, because designer does not depend on the bank; falling back to
# REPO_ROOT keeps this working against a pre-split checkout.
_NESTED = REPO_ROOT / "ledger"
LEDGER = _NESTED if (_NESTED / ".probe-ledger").is_file() else REPO_ROOT
# The backbone moved out of the experiments root when that tree was emptied
# ("every loose file moves to its own layer"); the bank is now its only home.
BACKBONE = LEDGER / "backbones" / "SP369_add_A_v1" / "sequences.xlsx"
BZ23_TCR_INPUT = (
    REPO_ROOT / "experiments" / "20260418_ZCH_BZ23_TCR"
    / "input" / "BZ23_tcr_clones.xlsx"
)


def _require_assets():
    if not BACKBONE.exists():
        pytest.skip(f"Backbone not present at {BACKBONE}")
    if not BZ23_TCR_INPUT.exists():
        pytest.skip(f"BZ23_TCR input not present at {BZ23_TCR_INPUT}")


def test_pipeline_smoke_bz23_drna_only_skip_blast(tmp_path: Path):
    """Default chemistry (dRNA only) — no RT primers, no BLAST round."""
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
        chemistries=["dRNA"],
        start_no=200,
        skip_blast=True,
        plot_tm_landscape=False,  # speed up the test; landscape covered separately
        codebook="SP369",
    )
    result = run_tcr_pipeline(cfg)

    assert result.clones_total >= 1
    assert result.clones_with_sites >= 1
    assert result.sites_selected >= result.clones_with_sites
    # dRNA only -> probes_total == sites_selected
    assert result.probes_total == result.sites_selected
    assert result.rt_primers_xlsx is None  # no cDNA -> no RT primers
    assert result.final_probes_xlsx.exists()
    assert result.final_probes_fasta.exists()

    # Verify final xlsx
    df = pd.read_excel(result.final_probes_xlsx)
    assert set(df["chemistry"]) == {"dRNA"}
    assert df.columns.tolist()[:3] == ["order", "probe_name", "probe_sequence"]
    assert df["No."].tolist() == sorted(df["No."].tolist())
    assert df["No."].iloc[0] == 200
    # probe_name follows new schema: {clone_id}_{pos}_{chem}_{codebook}_{No.}
    for name in df["probe_name"]:
        assert "_dRNA_" in name
        assert "_SP369_" in name


def test_pipeline_smoke_bz23_dual_chemistry(tmp_path: Path):
    """--chemistries dRNA,cDNA — dRNA + cDNA each independently selected;
    Nos are SHARED per (clone, site_index) across chemistries.
    """
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
        chemistries=["dRNA", "cDNA"],
        start_no=200,
        skip_blast=True,
        plot_tm_landscape=False,
        codebook="SP369",
    )
    result = run_tcr_pipeline(cfg)

    assert result.rt_primers_xlsx is not None
    assert result.rt_primers_xlsx.exists()

    df = pd.read_excel(result.final_probes_xlsx)
    assert set(df["chemistry"]) == {"dRNA", "cDNA"}
    # Shared-No. policy: within a clone, dRNA + cDNA rows of site_index 0
    # share the same `No.` AND the same `backbone_sequence`.
    for clone_id, sub in df.groupby("clone_id"):
        drna = sub[sub["chemistry"] == "dRNA"].sort_values("position")
        cdna = sub[sub["chemistry"] == "cDNA"].sort_values("position")
        if len(drna) > 0 and len(cdna) > 0:
            assert drna.iloc[0]["No."] == cdna.iloc[0]["No."], (
                f"{clone_id}: site#0 dRNA No.{drna.iloc[0]['No.']} != "
                f"cDNA No.{cdna.iloc[0]['No.']}"
            )
            assert drna.iloc[0]["backbone_sequence"] == cdna.iloc[0]["backbone_sequence"]


# ---------------------------------------------------------------------------
# Per-chemistry independent selection — the architectural change
# ---------------------------------------------------------------------------


def test_dual_chemistry_soft_selection_no_hard_gate(tmp_path: Path):
    """Soft Tm design (2026-07-19): there is NO hard Tm gate — arm Tm enters the
    per-chemistry score instead. An impossible tm_range no longer drops a
    chemistry; both dRNA and cDNA still produce probes, scored + selected
    INDEPENDENTLY per chemistry (not piggybacking a shared selection).
    """
    _require_assets()

    fake_experiment = tmp_path / "soft_selection"
    (fake_experiment / "input").mkdir(parents=True)
    (fake_experiment / "output").mkdir()
    input_xlsx = fake_experiment / "input" / "BZ23_tcr_clones.xlsx"
    shutil.copy2(BZ23_TCR_INPUT, input_xlsx)

    cfg = TcrConfig(
        input_xlsx=input_xlsx,
        output_dir=fake_experiment / "output",
        backbone_file=BACKBONE,
        chemistries=["dRNA", "cDNA"],
        start_no=200,
        sites_per_clone=2,
        tm_range=(200.0, 300.0),       # 2026-07-19: no longer a hard gate — ignored
        cdna_tm_range=(50.0, 80.0),
        skip_blast=True,
        plot_tm_landscape=False,
        codebook="SP369",
    )
    result = run_tcr_pipeline(cfg)
    df = pd.read_excel(result.final_probes_xlsx)
    # The impossible dRNA "gate" no longer drops dRNA probes under soft design.
    assert set(df["chemistry"]) == {"dRNA", "cDNA"}, (
        f"Expected BOTH chemistries under soft Tm design; got {set(df['chemistry'])}"
    )
    assert result.sites_selected_per_chem["dRNA"] > 0
    assert result.sites_selected_per_chem["cDNA"] > 0


# ---------------------------------------------------------------------------
# Shared-No. policy: dRNA + cDNA share No. per (clone, site_index)
# ---------------------------------------------------------------------------


def test_shared_no_simple_one_site_per_clone(tmp_path: Path):
    """One clone, both chemistries pick one site each — same No. + backbone."""
    _require_assets()

    fake_experiment = tmp_path / "shared_no_simple"
    (fake_experiment / "input").mkdir(parents=True)
    (fake_experiment / "output").mkdir()
    input_xlsx = fake_experiment / "input" / "BZ23_tcr_clones.xlsx"
    shutil.copy2(BZ23_TCR_INPUT, input_xlsx)

    cfg = TcrConfig(
        input_xlsx=input_xlsx,
        output_dir=fake_experiment / "output",
        backbone_file=BACKBONE,
        chemistries=["dRNA", "cDNA"],
        start_no=200,
        sites_per_clone=1,
        skip_blast=True,
        plot_tm_landscape=False,
        codebook="SP369",
    )
    result = run_tcr_pipeline(cfg)
    df = pd.read_excel(result.final_probes_xlsx)

    # For every clone with both chemistries, the two rows share No. and backbone.
    for clone_id, sub in df.groupby("clone_id"):
        chems = set(sub["chemistry"])
        if chems == {"dRNA", "cDNA"}:
            assert len(set(sub["No."])) == 1, (
                f"{clone_id}: dRNA + cDNA rows should share No.; got {sorted(sub['No.'])}"
            )
            assert len(set(sub["backbone_sequence"])) == 1


def test_shared_no_multi_site_pairs_by_index(tmp_path: Path):
    """Multi-site clone: dRNA[0] pairs with cDNA[0] on one No., dRNA[1] (no
    cDNA partner) takes the next No. alone."""
    import pandas as pd

    # Synthetic clone with a long enough allowed region that 2 non-overlapping
    # sites can be selected per chemistry.
    cdr3 = ("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
    flank = "AAAA" * 30
    chain = flank + cdr3 + flank

    input_xlsx = tmp_path / "input" / "syn.xlsx"
    input_xlsx.parent.mkdir(parents=True)
    pd.DataFrame([{
        "Clone_name": "SynClone",
        "consensus_id": "syn1",
        "cdr3_nt": cdr3,
        "TRB_nt": chain,
    }]).to_excel(input_xlsx, index=False)
    out_dir = tmp_path / "output"
    out_dir.mkdir()

    cfg = TcrConfig(
        input_xlsx=input_xlsx,
        output_dir=out_dir,
        backbone_file=BACKBONE,
        chemistries=["dRNA", "cDNA"],
        start_no=200,
        sites_per_clone=2,
        tm_range=(30.0, 90.0),
        cdna_tm_range=(30.0, 90.0),
        skip_blast=True,
        plot_tm_landscape=False,
        codebook="SP369",
    )
    result = run_tcr_pipeline(cfg)
    df = pd.read_excel(result.final_probes_xlsx).sort_values(["chemistry", "position"])

    # Worked-pairing check: per chemistry, sites sorted by position get
    # site_index 0, 1. Site index 0 across both chemistries → same No.
    drna = df[df["chemistry"] == "dRNA"].sort_values("position").reset_index(drop=True)
    cdna = df[df["chemistry"] == "cDNA"].sort_values("position").reset_index(drop=True)
    if len(drna) >= 1 and len(cdna) >= 1:
        assert drna.iloc[0]["No."] == cdna.iloc[0]["No."]
    # Whenever both chemistries have a second site, those share too.
    if len(drna) >= 2 and len(cdna) >= 2:
        assert drna.iloc[1]["No."] == cdna.iloc[1]["No."]
    # Numbering is contiguous starting from start_no.
    nos = sorted(set(df["No."]))
    assert nos == list(range(min(nos), max(nos) + 1))
    assert min(nos) == 200
