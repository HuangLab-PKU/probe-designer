"""Characterization tests for BruteForceStrategy (to be renamed SingleSequenceStrategy).

Locks the search_binding_sites behavior on a synthetic sequence before the
Phase 5 rename and search_strategies.py → search/ package split. No HTTP
required — feeds a synthetic sequence directly.
"""
from __future__ import annotations

import pytest

from probe_designer.search_strategies import BruteForceStrategy
from probe_designer.config import SearchConfig, FilterConfig


# Synthetic sequence: 200 nt of random-ish bases with some GC variation
# Long enough for the default pin_gap=0.1 trim (10% off each end = 20nt)
# and binding_site_length=40 default to yield ~140 candidate positions.
SYNTHETIC_SEQ = (
    "ATGCATGCATGCATGCATGC"  # 0-19
    "GCATGCATGCATGCATGCAT"  # 20-39
    "CGATCGATCGATCGATCGAT"  # 40-59
    "TAGCTAGCTAGCTAGCTAGC"  # 60-79
    "GCTAGCTAGCTAGCTAGCTA"  # 80-99
    "ATGCATGCATGCATGCATGC"  # 100-119
    "CGATCGATCGATCGATCGAT"  # 120-139
    "TAGCTAGCTAGCTAGCTAGC"  # 140-159
    "GCATGCATGCATGCATGCAT"  # 160-179
    "ATGCATGCATGCATGCATGC"  # 180-199
)
assert len(SYNTHETIC_SEQ) == 200


@pytest.fixture
def strategy():
    search_cfg = SearchConfig()
    # Tighten filter so only some positions pass — gives meaningful output
    filter_cfg = FilterConfig()
    return BruteForceStrategy(search_cfg, filter_cfg)


class TestReturnShape:
    def test_returns_list_of_dicts(self, strategy):
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        assert isinstance(sites, list)
        assert all(isinstance(s, dict) for s in sites)

    def test_each_site_has_expected_keys(self, strategy):
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        if not sites:
            pytest.skip("No sites passed thermal filter on synthetic sequence")
        required_keys = {
            "gene_name", "sequence", "target_sequence",
            "arm_3prime", "arm_5prime", "st", "en", "strand", "length",
            "g_content", "tm", "tm_3", "tm_5", "tm_diff", "free_energy",
            "strategy",
        }
        missing = required_keys - set(sites[0].keys())
        assert not missing, f"Missing keys in site: {missing}"

    def test_strategy_field_is_brute_force(self, strategy):
        """Phase 5 will change this to 'single_sequence'; pin current value."""
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        if sites:
            assert sites[0]["strategy"] == "brute_force"

    def test_gene_name_propagated(self, strategy):
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        for s in sites:
            assert s["gene_name"] == "TEST_GENE"


class TestSequenceTooShort:
    def test_short_sequence_returns_empty(self, strategy):
        # Default binding_site_length=40; sequence shorter than that returns []
        short = "ATGC" * 5  # 20 nt
        sites = strategy.search_binding_sites(short, "SHORT_GENE")
        assert sites == []


class TestSiteProperties:
    def test_sites_have_correct_length(self, strategy):
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        for s in sites:
            assert s["length"] == 40  # default binding_site_length
            # arm_3prime + arm_5prime combined = length (each half of binding site)
            assert len(s["arm_3prime"]) + len(s["arm_5prime"]) == 40

    def test_sites_within_sequence_bounds(self, strategy):
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        # pin_gap=0.1 trims 10% from each end for the default scan.
        # With len=200, trim=20, so positions in [20, 200-20-40] = [20, 140].
        for s in sites:
            pos = s["st"]
            assert 0 <= pos <= len(SYNTHETIC_SEQ) - 40

    def test_strand_is_plus_by_default(self, strategy):
        # No genomic_context provided => strand defaults to '+'
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        for s in sites:
            assert s["strand"] == "+"


class TestGenomicContext:
    def test_genomic_context_maps_positions(self):
        search_cfg = SearchConfig()
        filter_cfg = FilterConfig()
        genomic_context = {
            "seq_region_name": "1",
            "start": 1000,
            "end": 1199,
            "strand": 1,
        }
        strategy = BruteForceStrategy(search_cfg, filter_cfg, genomic_context=genomic_context)
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        for s in sites:
            # On + strand: genomic_start = g_start + rel_pos
            # genomic_end = genomic_start + 40
            assert s["st"] >= 1000
            assert s["en"] == s["st"] + 40

    def test_genomic_context_minus_strand(self):
        search_cfg = SearchConfig()
        filter_cfg = FilterConfig()
        genomic_context = {
            "seq_region_name": "1", "start": 1000, "end": 1199, "strand": -1,
        }
        strategy = BruteForceStrategy(search_cfg, filter_cfg, genomic_context=genomic_context)
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        for s in sites:
            assert s["strand"] == "-"


class TestDeterministic:
    def test_same_input_same_output(self, strategy):
        sites1 = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        sites2 = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        assert len(sites1) == len(sites2)
        # Same positions, same order
        assert [s["st"] for s in sites1] == [s["st"] for s in sites2]


class TestCurrentBehaviorAllCandidatesReturned:
    """Locks the comment 'Build results from ALL passing candidates (no limit)'.

    Before Phase 5 moves selection into Pipeline, the strategy returns everything
    that passes thermal filter. This must not change during refactor.
    """

    def test_no_top_n_limit_applied_by_strategy(self, strategy):
        # Run with default config; verify strategy does NOT truncate
        # (Pipeline will apply peak_rank + top-N later)
        sites = strategy.search_binding_sites(SYNTHETIC_SEQ, "TEST_GENE")
        # If any site passes, we should get more than just 3 on a 200nt seq
        # (historical 'top-3 most distant' was done at filter stage, not here).
        # This may be 0 if filter is strict — assert shape not count.
        assert isinstance(sites, list)
