"""Unit tests for ``probe_designer.qc.cross_lig_check``.

Tests the thin wrapper that lazy-imports ``bank.probe_book.pool.*`` to
run cross-ligation screening against (a) the candidate batch itself and
optionally (b) an existing pool. The bank package is editable-installed
in our env so we can drive the real screen end-to-end with a tiny mock
pool.
"""
from __future__ import annotations

from pathlib import Path

import pytest


pytest.importorskip("primer3")
pytest.importorskip("probe_book")


from probe_designer.qc.cross_lig_check import (
    CandidateProbe,
    screen_candidates_against_pool,
)


# A safe (no cross-lig) candidate pair — random 20-mers that share no
# meaningful k-mer with their reverse complement.
SAFE_A = CandidateProbe(
    probe_id="cand_safe_A", chemistry="dRNA", target="GENE_A",
    probe_arm5="GCATAGCAGCAGCAGCATAG", probe_arm3="TGTGTGTGCACGCACGCATG",
)
SAFE_B = CandidateProbe(
    probe_id="cand_safe_B", chemistry="dRNA", target="GENE_B",
    probe_arm5="AAGAAGAAGAAGAAGAAGAA", probe_arm3="TTCTTCTTCTTCTTCTTCTT",
)


# A clear cross-lig pair. CROSSLIG_A's wrap-around junction
# arm3[-3:] + arm5[:3] = "TCG|GCA" — and CROSSLIG_B's arm5 contains the RC
# of that ("TGC|CGA"). Their full duplex Tm should exceed 40 °C.
CROSSLIG_A = CandidateProbe(
    probe_id="cand_xlig_A", chemistry="dRNA", target="GENE_X",
    probe_arm5="GCAGCAGCAGCAGCAGCAGC", probe_arm3="TGCTGCTGCTGCTGCTGCTG",
)
CROSSLIG_B = CandidateProbe(
    probe_id="cand_xlig_B", chemistry="dRNA", target="GENE_Y",
    probe_arm5="CAGCAGCAGCAGCAGCAGCA", probe_arm3="GCTGCTGCTGCTGCTGCTGC",
)


def test_within_batch_no_pool_no_hits():
    hits, by_cand = screen_candidates_against_pool([SAFE_A, SAFE_B])
    assert hits == []
    assert by_cand == {}


def test_within_batch_detects_obvious_dimer():
    hits, by_cand = screen_candidates_against_pool([CROSSLIG_A, CROSSLIG_B])
    assert len(hits) >= 1
    pair_ids = {(h.probe_a_id, h.probe_b_id) for h in hits}
    expected = {("cand_xlig_A", "cand_xlig_B"), ("cand_xlig_B", "cand_xlig_A")}
    assert pair_ids & expected
    # by_candidate index should list this hit under both endpoints
    assert "cand_xlig_A" in by_cand and "cand_xlig_B" in by_cand


def test_with_pool_marks_existing_side(tmp_path: Path):
    """Mock a tiny pool registry; verify pool-vs-candidate hits get marked."""
    # Build a fake project root with one pool + matching probe registry.
    (tmp_path / "probes").mkdir()
    registry_tsv = tmp_path / "probes" / "registry.tsv"
    # Make 1 pool probe that will cross-lig with CROSSLIG_B
    registry_tsv.write_text(
        "probe_id\tcodebook\tchemistry\tprobe_arm5\tprobe_arm3\tsequence\ttarget\n"
        f"pool_xlig_A\tSP369\tdRNA\t{CROSSLIG_A.probe_arm5}\t{CROSSLIG_A.probe_arm3}\t\tGENE_X\n",
        encoding="utf-8",
    )
    pool_dir = tmp_path / "pools" / "test_pool_v1"
    pool_dir.mkdir(parents=True)
    (pool_dir / "manifest.yaml").write_text(
        "pool_id: test_pool_v1\n"
        "pool_name: test\n"
        "codebook: SP369\n"
        "target_well_volume_uL: 10.0\n"
        "recipe:\n"
        "- order_id: ORD1\n"
        "  probe_id: pool_xlig_A\n"
        "  concentration_M: 1.0e-08\n",
        encoding="utf-8",
    )

    hits, by_cand = screen_candidates_against_pool(
        [CROSSLIG_B],
        pool_id="test_pool_v1",
        repo_root=tmp_path,
    )
    assert len(hits) >= 1
    # Each hit involves exactly one candidate and one pool probe.
    for h in hits:
        endpoints = {(h.probe_a_id, h.a_is_existing_pool),
                       (h.probe_b_id, h.b_is_existing_pool)}
        assert ("cand_xlig_B", False) in endpoints
        assert ("pool_xlig_A", True) in endpoints
    assert "cand_xlig_B" in by_cand


def test_empty_candidates_returns_empty():
    hits, by_cand = screen_candidates_against_pool([])
    assert hits == [] and by_cand == {}


def test_single_candidate_no_partner():
    hits, by_cand = screen_candidates_against_pool([SAFE_A])
    assert hits == [] and by_cand == {}
