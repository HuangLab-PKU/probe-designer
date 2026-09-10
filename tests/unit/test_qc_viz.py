"""Unit tests for ``probe_designer.qc.viz``.

These are smoke tests — matplotlib visual artifacts are hard to assert on
pixel-by-pixel, so we verify the public API contract:

* All three functions accept minimal valid inputs.
* When a save_path is given, a non-empty PNG file appears at that path.
* Each function returns a matplotlib Figure.
* Edge cases (no probes, 1 isoform, missing profile) don't crash.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


from probe_designer.qc.viz import (
    IsoformMeta,
    ProbeAnnotation,
    plot_accessibility_tracks,
    plot_accessibility_heatmap,
    plot_rnafold_structure,
)


def _toy_isoform(tid: str = "ENST_TOY", *, length: int = 200) -> IsoformMeta:
    return IsoformMeta(
        transcript_id=tid,
        display_name=tid.replace("ENST_", "") + "-201",
        chromosome="X",
        strand=1,
        exons=[(1000, 1000 + length - 1)],
        is_canonical=True,
        biotype="protein_coding",
    )


def _toy_profile(n: int = 200, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.uniform(0.3, 0.95, size=n)


def _toy_probe(start: int = 1050, end: int = 1090, outcome: str = "signal") -> ProbeAnnotation:
    return ProbeAnnotation(
        probe_id="SP_TOY_1",
        arm5="A" * 20, arm3="T" * 20,
        genomic_start=start, genomic_end=end,
        outcome=outcome,
    )


# ---------- plot_accessibility_tracks ----------


def test_tracks_minimal_inputs(tmp_path: Path):
    iso = _toy_isoform()
    profiles = {iso.transcript_id: {37.0: _toy_profile()}}
    out = tmp_path / "tracks.png"
    fig = plot_accessibility_tracks(
        [iso], profiles, temperatures=(37.0,), save_path=out,
    )
    assert fig is not None
    assert out.exists() and out.stat().st_size > 0


def test_tracks_multi_isoform_multi_temp(tmp_path: Path):
    isos = [_toy_isoform(f"ENST_T{i}", length=150 + 50 * i) for i in range(3)]
    profiles = {
        iso.transcript_id: {T: _toy_profile(seq_len, seed=i)
                             for T in (37.0, 47.0)
                             for seq_len in [iso.exons[0][1] - iso.exons[0][0] + 1]}
        for i, iso in enumerate(isos)
    }
    out = tmp_path / "tracks_multi.png"
    fig = plot_accessibility_tracks(
        isos, profiles, temperatures=(37.0, 47.0), save_path=out,
    )
    assert out.exists() and out.stat().st_size > 0


def test_tracks_with_probes(tmp_path: Path):
    iso = _toy_isoform()
    profiles = {iso.transcript_id: {37.0: _toy_profile()}}
    probes = [
        _toy_probe(1050, 1090, "signal"),
        ProbeAnnotation("SP_TOY_2", "C" * 20, "G" * 20, 1100, 1140, "dead"),
    ]
    out = tmp_path / "tracks_probes.png"
    plot_accessibility_tracks(
        [iso], profiles, probes=probes,
        temperatures=(37.0,), save_path=out,
    )
    assert out.exists() and out.stat().st_size > 0


def test_tracks_minus_strand(tmp_path: Path):
    """Minus-strand isoform should still render — profile slice is reversed."""
    iso = IsoformMeta(
        transcript_id="ENST_MINUS",
        display_name="MINUS-201",
        chromosome="2",
        strand=-1,
        exons=[(5000, 5199)],
        is_canonical=True,
    )
    profiles = {iso.transcript_id: {37.0: _toy_profile()}}
    out = tmp_path / "tracks_minus.png"
    plot_accessibility_tracks([iso], profiles, save_path=out)
    assert out.exists() and out.stat().st_size > 0


# ---------- plot_accessibility_heatmap ----------


def test_heatmap_minimal(tmp_path: Path):
    iso = _toy_isoform()
    profiles = {iso.transcript_id: _toy_profile()}
    out = tmp_path / "hm.png"
    fig = plot_accessibility_heatmap(
        [iso], profiles, temperature=37.0, save_path=out,
    )
    assert fig is not None
    assert out.exists() and out.stat().st_size > 0


def test_heatmap_multi_exon_isoform(tmp_path: Path):
    iso = IsoformMeta(
        transcript_id="ENST_MULTI",
        display_name="MULTI-201",
        chromosome="3",
        strand=1,
        exons=[(1000, 1099), (2000, 2099)],  # 200 nt mRNA across two exons
        biotype="protein_coding",
    )
    profile = _toy_profile(200)
    out = tmp_path / "hm_multi.png"
    plot_accessibility_heatmap([iso], {iso.transcript_id: profile},
                                temperature=37.0, save_path=out)
    assert out.exists() and out.stat().st_size > 0


# ---------- plot_rnafold_structure ----------


def test_rnafold_short_seq(tmp_path: Path):
    seq = "GGGAAACCCTTTGGGAAACCC"  # forms a simple hairpin
    out = tmp_path / "rna.png"
    fig = plot_rnafold_structure(seq, save_path=out, transcript_id="TOY")
    assert fig is not None
    assert out.exists() and out.stat().st_size > 0


def test_rnafold_with_probe_window(tmp_path: Path):
    seq = "G" * 50 + "AAACCC" * 10 + "G" * 50
    out = tmp_path / "rna_probe.png"
    plot_rnafold_structure(
        seq, probe_windows=[(60, 100)], save_path=out, transcript_id="TOY2",
    )
    assert out.exists() and out.stat().st_size > 0


def test_rnafold_window_subset(tmp_path: Path):
    """window_size centered on the probe should fold only a slice."""
    seq = "A" * 100 + "GCGCGCGCGCGCGCGCGCGC" + "T" * 100  # 220 nt
    out = tmp_path / "rna_window.png"
    plot_rnafold_structure(
        seq, probe_windows=[(100, 120)], window_size=80,
        save_path=out, transcript_id="TOY3",
    )
    assert out.exists() and out.stat().st_size > 0



def test_plot_ternary_complex_renders(tmp_path: Path):
    """v2.3 NUPACK ternary visualization smoke test — feed fabricated
    dot-bracket + metadata, verify figure writes."""
    from probe_designer.qc.viz import plot_ternary_complex

    bb = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    arm5_a = "GCAGCAGCAGCAGCAGCAGC"
    arm3_a = "TGCTGCTGCTGCTGCTGCTG"
    arm5_b = "CAGCAGCAGCAGCAGCAGCA"
    arm3_b = "GCTGCTGCTGCTGCTGCTGC"

    out = tmp_path / "ternary.png"
    plot_ternary_complex(
        probe_a={"probe_id": "SP_A", "probe_arm5": arm5_a, "probe_arm3": arm3_a,
                  "sequence": arm5_a + bb + arm3_a, "target": "GX", "chemistry": "dRNA"},
        probe_b={"probe_id": "SP_B", "probe_arm5": arm5_b, "probe_arm3": arm3_b,
                  "sequence": arm5_b + bb + arm3_b, "target": "GY", "chemistry": "dRNA"},
        # 3-strand dot-bracket: splint | arm5 | arm3, all unpaired (placeholder)
        mfe_dotbracket="." * 95 + "+" + "." * 20 + "+" + "." * 20,
        ternary_dg_kcal=-40.0,
        ternary_concentration_m=1e-10,
        ternary_fraction_of_b=1e-3,
        ensemble_used="stacking",
        celsius_used=55.0,
        sodium_m=0.075, magnesium_m=0.010, strand_conc_m=1e-7,
        mfe_nick_adjacent=False,
        mfe_vicinity_contiguous=False,
        b_3oh_pos=None, b_5p_pos=None,
        vicinity_n_each_side=3,
        save_path=out,
    )
    assert out.exists() and out.stat().st_size > 0
