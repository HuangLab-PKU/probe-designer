"""Unit tests for v2 ``probe_designer.qc.cross_lig_check.screen_candidates``.

Tests the within-batch + splint-pool screening path used by the design CLIs.
Bank-free — splint probes are passed in as ``ProbeForScreen`` directly
(the CLI layer's ``ext.pool.loader`` constructs them from bank data).
"""
from __future__ import annotations

import pytest

pytest.importorskip("primer3")


from probe_designer.qc.cross_lig_check import (
    CandidateProbe,
    screen_candidates,
)
from probe_designer.qc.cross_ligation import ProbeForScreen


_BB = "TCCCTACACGACGCTCTTCCG"   # 21-nt neutral backbone for fixtures


def _candidate(probe_id: str, arm5: str, arm3: str, *,
                chemistry: str = "dRNA", target: str = "TEST") -> CandidateProbe:
    return CandidateProbe(
        probe_id=probe_id, chemistry=chemistry,
        probe_arm5=arm5, probe_arm3=arm3,
        sequence=(arm5 + _BB + arm3).upper(),
        target=target,
    )


def _pfs(probe_id: str, arm5: str, arm3: str, *,
         chemistry: str = "dRNA", target: str = "POOL") -> ProbeForScreen:
    return ProbeForScreen(
        probe_id=probe_id, chemistry=chemistry,
        probe_arm5=arm5, probe_arm3=arm3,
        sequence=(arm5 + _BB + arm3).upper(),
        target=target,
    )


# Safe pair — no junction-templating overlap
SAFE_A = _candidate("safe_A", "GCATAGCAGCAGCAGCATAG", "TGTGTGTGCACGCACGCATG", target="G1")
SAFE_B = _candidate("safe_B", "AAGAAGAAGAAGAAGAAGAA", "TTCTTCTTCTTCTTCTTCTT", target="G2")


# True v2 cross-lig: B's arm5 is RC of A's rotated arm3+arm5 (40 nt) so the
# full A ligator pairs contiguously inside B's arm5 (no bb-spanning).
_A_ARM3 = "GGGGGAAAATTTCCCAAGGG"
_A_ARM5 = "CCCAATTGCGCAATATCATG"
_B_ARM5_FROM_A = "CATGATATTGCGCAATTGGGCCCTTGGGAAATTTTCCCCC"  # = RC(arm3_A + arm5_A)

XLIG_A = _candidate("xlig_A", _A_ARM5, _A_ARM3, target="GX")
XLIG_B = CandidateProbe(
    probe_id="xlig_B", chemistry="dRNA",
    probe_arm5=_B_ARM5_FROM_A, probe_arm3="AAGCTTAACTGGCCATAAGT",
    sequence=(_B_ARM5_FROM_A + _BB + "AAGCTTAACTGGCCATAAGT").upper(),
    target="GY",
)


# ----------------------------------------------------------------------
# Within-batch only (no splint pool)
# ----------------------------------------------------------------------


def test_empty_candidates_returns_empty():
    hits, by_cand = screen_candidates([])
    assert hits == [] and by_cand == {}


def test_single_candidate_no_partner():
    hits, by_cand = screen_candidates([SAFE_A])
    assert hits == [] and by_cand == {}


def test_within_batch_safe_pair_unflagged():
    hits, by_cand = screen_candidates([SAFE_A, SAFE_B])
    assert hits == []
    assert by_cand == {}


def test_within_batch_xlig_pair_flagged():
    """A's rotated arms RC-match B's arm5 contiguously → flag."""
    hits, by_cand = screen_candidates([XLIG_A, XLIG_B])
    assert len(hits) >= 1
    # Both endpoints of any confirmed hit appear in by_cand
    for h in hits:
        assert h.dimer.junction_run_nt >= 2 * (h.dimer.vicinity_n_each_side + 1)
        assert h.probe_a_id in by_cand


# ----------------------------------------------------------------------
# Pool-aware (with splint_probes)
# ----------------------------------------------------------------------


def test_pool_partner_marks_existing_side():
    """One candidate; the pool has a probe whose arm5 is RC of the candidate's
    rotated arms. Hit should be flagged AND mark b_is_existing_pool=True.
    """
    splint = _pfs("pool_xlig_partner", _B_ARM5_FROM_A, "AAGCTTAACTGGCCATAAGT", target="POOL_G")
    hits, by_cand = screen_candidates([XLIG_A], splint_probes=[splint])
    confirmed = [h for h in hits if h.dimer.flagged_overall]
    assert confirmed, f"expected ≥1 confirmed hit; got {hits}"
    for h in confirmed:
        # XLIG_A is candidate (a_is_existing_pool=False); pool_xlig_partner is pool (b_is_existing_pool=True)
        assert h.a_is_existing_pool is False
        assert h.b_is_existing_pool is True
        assert h.probe_a_id == "xlig_A"
        assert h.probe_b_id == "pool_xlig_partner"
    # by_cand keyed only by candidate side (= XLIG_A)
    assert "xlig_A" in by_cand
    assert "pool_xlig_partner" not in by_cand   # pool side not indexed


def test_pool_vs_pool_pairs_are_dropped():
    """Pool-vs-pool dimers must not appear in results (they're not this run's concern)."""
    pool1 = _pfs("pool_1", _B_ARM5_FROM_A, "AAGCTTAACTGGCCATAAGT")
    pool2 = _pfs("pool_2", _A_ARM5, _A_ARM3)   # XLIG_A's arms but as pool member
    hits, by_cand = screen_candidates([SAFE_A], splint_probes=[pool1, pool2])
    # SAFE_A doesn't cross-lig with either pool member → no candidate-involving hits
    # Even though pool1 + pool2 form a dimer, it's filtered out (both endpoints in pool)
    for h in hits:
        assert not (h.a_is_existing_pool and h.b_is_existing_pool)
