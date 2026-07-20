"""Tests for filtering/accessibility.py (Phase 1A, 2026-05-14)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from probe_designer.filtering.accessibility import (
    aggregate_open_probability,
    cached_profile,
    compute_plfold_profile,
    mean_open_probability,
)


class TestComputePlfoldProfile:
    def test_returns_numpy_array_of_seq_length(self):
        seq = "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU"  # 40 nt
        prof = compute_plfold_profile(seq, window=70, span=40)
        assert isinstance(prof, np.ndarray)
        assert prof.shape == (40,)

    def test_probabilities_in_unit_interval(self):
        seq = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"  # forms strong hairpins
        prof = compute_plfold_profile(seq, window=70, span=40)
        assert (prof >= 0.0).all()
        assert (prof <= 1.0).all()

    def test_accepts_dna_letters(self):
        """T should be treated like U (ViennaRNA convention)."""
        seq_dna = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        seq_rna = "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU"
        prof_dna = compute_plfold_profile(seq_dna, window=70, span=40)
        prof_rna = compute_plfold_profile(seq_rna, window=70, span=40)
        # Both should produce the same numbers up to float noise.
        assert np.allclose(prof_dna, prof_rna)

    def test_clip_to_window_when_seq_shorter(self):
        """Window > len(seq) must not crash; effective window is clipped."""
        seq = "ACGUACGUACGUACGU"  # 16 nt
        prof = compute_plfold_profile(seq, window=70, span=40)
        assert prof.shape == (16,)
        assert (prof >= 0.0).all()

    def test_empty_seq_raises(self):
        with pytest.raises(ValueError):
            compute_plfold_profile("", window=70, span=40)

    def test_temperature_change_affects_profile(self):
        """Higher temperature destabilizes structure → higher p_unpaired on
        a hairpin-prone sequence. This is the "formamide" recalibration
        path the module docstring mentions."""
        seq = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
        cold = compute_plfold_profile(seq, window=70, span=40, temperature=20.0)
        hot = compute_plfold_profile(seq, window=70, span=40, temperature=60.0)
        # The middle bases should be more open at higher T on average.
        assert hot.mean() > cold.mean()


class TestMeanOpenProbability:
    def test_correct_window_average(self):
        prof = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        # mean over [1, 4) = mean(0.2, 0.3, 0.4) = 0.3
        assert mean_open_probability(prof, 1, 4) == pytest.approx(0.3)

    def test_full_range(self):
        prof = np.array([0.5, 0.5, 0.5, 0.5])
        assert mean_open_probability(prof, 0, 4) == pytest.approx(0.5)

    def test_out_of_range_returns_zero(self):
        prof = np.array([0.5, 0.5, 0.5, 0.5])
        assert mean_open_probability(prof, 0, 5) == 0.0  # end > len
        assert mean_open_probability(prof, -1, 3) == 0.0  # start < 0
        assert mean_open_probability(prof, 3, 3) == 0.0  # empty window


class TestAggregateOpenProbability:
    def test_empty_input_returns_zero(self):
        assert aggregate_open_probability({}) == 0.0

    def test_uniform_weighting_default(self):
        # Mean of {a: 0.4, b: 0.6} = 0.5
        out = aggregate_open_probability({"a": 0.4, "b": 0.6})
        assert out == pytest.approx(0.5)

    def test_weighted_uses_provided_weights(self):
        # Dominant isoform b (weight 0.9) is folded (p=0.1);
        # rare isoform a (weight 0.1) is open (p=1.0).
        # Weighted score = 0.1*1.0 + 0.9*0.1 = 0.19 — close to dominant.
        out = aggregate_open_probability(
            {"a": 1.0, "b": 0.1},
            weights={"a": 0.1, "b": 0.9},
        )
        assert out == pytest.approx(0.19)

    def test_missing_isoform_in_weights_treated_as_zero(self):
        out = aggregate_open_probability(
            {"a": 1.0, "b": 0.5},
            weights={"a": 1.0},  # b missing → weight 0, excluded
        )
        assert out == pytest.approx(1.0)

    def test_zero_total_weight_returns_zero(self):
        out = aggregate_open_probability(
            {"a": 0.5}, weights={"a": 0.0},
        )
        assert out == 0.0


class TestCachedProfile:
    def test_cache_writes_and_reads(self, tmp_path):
        seq = "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU"
        prof1 = cached_profile(seq, "ENST_test", cache_dir=tmp_path)
        # Second call must hit the cached file (same array)
        prof2 = cached_profile(seq, "ENST_test", cache_dir=tmp_path)
        assert np.array_equal(prof1, prof2)
        # File should exist on disk
        cached_files = list(tmp_path.glob("ENST_test_*.npy"))
        assert len(cached_files) == 1

    def test_cache_invalidates_on_length_mismatch(self, tmp_path):
        """If the cached file's length doesn't match seq length, regenerate."""
        seq_a = "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU"  # 40 nt
        seq_b = "ACGUACGUACGUACGU"  # 16 nt
        prof_a = cached_profile(seq_a, "ENST_test", cache_dir=tmp_path)
        # Same transcript_id but different length — cache key changes
        # because length is in the filename, so we should see two files.
        prof_b = cached_profile(seq_b, "ENST_test", cache_dir=tmp_path)
        assert prof_a.shape != prof_b.shape

    def test_cache_key_distinguishes_window_and_temperature(self, tmp_path):
        seq = "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU"
        cached_profile(seq, "ENST_x", cache_dir=tmp_path, window=70, span=40, temperature=37.0)
        cached_profile(seq, "ENST_x", cache_dir=tmp_path, window=70, span=40, temperature=50.0)
        files = list(tmp_path.glob("ENST_x_*.npy"))
        # Two files because temperature differs
        assert len(files) == 2
