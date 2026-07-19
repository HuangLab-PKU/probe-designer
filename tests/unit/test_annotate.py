"""First-batch reference thermodynamic annotations (arm Tm + accessibility)."""
from __future__ import annotations

import pytest

from probe_designer.annotate import build_reference_annotations
from probe_designer.chemistry import ReactionConditions

# A short mRNA-sense reference.
SEQ = "GCATTCAGGTCACCTTGATGCATTCAGGTCACCTTGATGCATTCAGGTCACCTTGATG"


def _read_values(path):
    vals = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("track"):
            continue
        vals.append(float(line.split("\t")[3]))
    return vals


class TestBuildReferenceAnnotations:
    def test_writes_arm_tm_track(self, tmp_path):
        paths = build_reference_annotations(
            SEQ, "TX1", ReactionConditions(),
            out_dir=tmp_path, arm_length=20, accessibility=False,
        )
        tm_tracks = [p for p in paths if "armTm" in p.name]
        assert len(tm_tracks) == 1
        assert tm_tracks[0].exists()
        # Tm values are physical (roughly 40-80 C for these arms).
        vals = _read_values(tm_tracks[0])
        assert vals and all(30 < v < 90 for v in vals)

    def test_filename_encodes_buffer_signature(self, tmp_path):
        p_default = build_reference_annotations(
            SEQ, "TX1", ReactionConditions(), out_dir=tmp_path, accessibility=False,
        )[0]
        p_lowmg = build_reference_annotations(
            SEQ, "TX1", ReactionConditions(mg_mM=4.0), out_dir=tmp_path, accessibility=False,
        )[0]
        assert p_default.name != p_lowmg.name  # different condition -> different file

    def test_accessibility_track_when_viennarna_present(self, tmp_path):
        pytest.importorskip("RNA")
        paths = build_reference_annotations(
            SEQ, "TX1", ReactionConditions(), out_dir=tmp_path, accessibility=True,
        )
        acc = [p for p in paths if "accessibility" in p.name]
        assert len(acc) == 1
        vals = _read_values(acc[0])
        assert vals and all(0.0 <= v <= 1.0 for v in vals)  # p_unpaired in [0,1]
