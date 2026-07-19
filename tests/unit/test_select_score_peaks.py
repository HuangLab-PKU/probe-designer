"""Greedy local-maximum (NMS) peak selection for probe picking."""
from __future__ import annotations

from probe_designer.scoring import select_score_peaks


def _site(st, score, bds=40):
    return {"st": st, "en": st + bds - 1, "score": score}


class TestSelectScorePeaks:
    def test_empty_returns_empty(self):
        assert select_score_peaks([], 40) == []

    def test_picks_global_max_first(self):
        sites = [_site(0, 3.0), _site(100, 9.0), _site(200, 5.0)]
        peaks = select_score_peaks(sites, min_distance=40)
        assert peaks[0]["st"] == 100  # score 9 first

    def test_suppresses_within_min_distance(self):
        # Two near sites (10 bp apart) + one far. min_distance 40 => the lower
        # of the near pair is suppressed by the higher.
        sites = [_site(100, 9.0), _site(110, 8.0), _site(300, 7.0)]
        peaks = select_score_peaks(sites, min_distance=40)
        starts = sorted(p["st"] for p in peaks)
        assert starts == [100, 300]  # 110 suppressed by 100's peak

    def test_non_overlap_spacing(self):
        # Sites every 20 bp; min_distance 40 => every other one at most.
        sites = [_site(i * 20, 10.0 - i) for i in range(6)]  # scores 10,9,8,...
        peaks = select_score_peaks(sites, min_distance=40)
        positions = sorted(p["st"] for p in peaks)
        for a, b in zip(positions, positions[1:]):
            assert b - a >= 40

    def test_max_n_caps_count(self):
        sites = [_site(i * 100, float(i)) for i in range(10)]
        peaks = select_score_peaks(sites, min_distance=40, max_n=3)
        assert len(peaks) == 3
        # highest scores are the last (i=9,8,7)
        assert {p["st"] for p in peaks} == {900, 800, 700}

    def test_output_in_score_desc_order(self):
        sites = [_site(0, 5.0), _site(100, 9.0), _site(200, 7.0)]
        peaks = select_score_peaks(sites, min_distance=40)
        scores = [p["score"] for p in peaks]
        assert scores == sorted(scores, reverse=True)

    def test_unpositioned_dropped(self):
        sites = [{"score": 9.0}, _site(100, 5.0)]  # first has no st/en
        peaks = select_score_peaks(sites, min_distance=40)
        assert len(peaks) == 1
        assert peaks[0]["st"] == 100
