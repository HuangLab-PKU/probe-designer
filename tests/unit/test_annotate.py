"""First-batch reference thermodynamic annotations (arm Tm + accessibility)."""
from __future__ import annotations

import pytest

from probe_designer.annotate import (
    build_canonical_genome_annotations,
    build_reference_annotations,
    emit_annotations_for_sequences,
)
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


class TestEmitForSequences:
    def test_batch_writes_each_reference(self, tmp_path):
        seqs = {"TX1": SEQ, "TX2": SEQ[:40]}
        paths = emit_annotations_for_sequences(
            seqs, ReactionConditions(), tmp_path, accessibility=False,
        )
        assert {p for p in paths if "TX1_armTm" in p.name}
        assert {p for p in paths if "TX2_armTm" in p.name}

    def test_skips_non_string_and_empty(self, tmp_path):
        seqs = {"TX1": SEQ, "BAD": None, "EMPTY": ""}
        paths = emit_annotations_for_sequences(
            seqs, ReactionConditions(), tmp_path, accessibility=False,
        )
        names = " ".join(p.name for p in paths)
        assert "TX1" in names and "BAD" not in names and "EMPTY" not in names


class TestCanonicalGenomeProjection:
    def test_genome_bedgraph_in_genomic_coords(self, tmp_path):
        # 57-nt SEQ split across two exons (29 + 28 = 57) on chr7 + strand.
        iso = {
            "display_name": "ENST_demo", "seq_region_name": "chr7", "strand": 1,
            "Exon": [{"start": 1000, "end": 1028}, {"start": 2000, "end": 2027}],
        }
        paths = build_canonical_genome_annotations(
            iso, SEQ, ReactionConditions(), out_dir=tmp_path,
            arm_length=20, accessibility=False, gene="DEMO",
        )
        assert len(paths) == 1 and paths[0].name.endswith(".genome.bedgraph")
        lines = [l for l in paths[0].read_text(encoding="utf-8").splitlines()
                 if not l.startswith("track")]
        chroms = {l.split("\t")[0] for l in lines}
        starts = [int(l.split("\t")[1]) for l in lines]
        assert chroms == {"chr7"}
        # first exon starts at genomic 1000 (0-based 999); values land there.
        assert min(starts) == 999
        assert any(s >= 1999 for s in starts)  # second exon block too
