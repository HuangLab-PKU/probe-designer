"""Tests for probe_designer.blast.gene_aware (pure helpers + XML parsing)."""
from __future__ import annotations

from pathlib import Path

import pytest

from probe_designer.blast.gene_aware import (
    apply_gene_aware_filter,
    extract_binding_query,
    find_best_per_mutation_record,
    is_same_gene_hit,
    parse_blast_results,
)


# ---------------------------------------------------------------------------
# extract_binding_query — the iLock dedup correctness
# ---------------------------------------------------------------------------

class TestExtractBindingQuery:
    def test_non_ilock_concatenates_arm3_then_arm5(self):
        q = extract_binding_query("AAAA", "TTTT", "mRNA_noiLock")
        assert q == "TTTTAAAA"

    def test_cdna_uses_same_concatenation(self):
        q = extract_binding_query("AAAA", "TTTT", "cDNA")
        assert q == "TTTTAAAA"

    def test_ilock_drops_shared_snp_nt(self):
        # arm3[-1] and arm5[0] both carry RC(SNP); query must NOT double-count.
        q = extract_binding_query("Caaaaaa", "Tttttttt", "iLock")
        # arm3[:-1] + arm5 = "Tttttttt"[:-1] is wrong; correct is arm3[:-1] of arm3.
        # Wait: arm3 here = "Tttttttt"; arm3[:-1] = "Ttttttt"; + arm5("Caaaaaa") = "TtttttttCaaaaaa"
        # Mistake — arm3[:-1] drops the LAST nt of arm3, then concat with arm5.
        # That's the contiguous binding region.
        assert q == "TTTTTTTCAAAAAA"  # arm3[:-1].upper() + arm5.upper()

    def test_empty_arm_returns_empty(self):
        assert extract_binding_query("", "TTT", "iLock") == ""
        assert extract_binding_query("AAA", "", "iLock") == ""


# ---------------------------------------------------------------------------
# is_same_gene_hit — gene-name matching nuances
# ---------------------------------------------------------------------------

class TestIsSameGeneHit:
    def test_exact_uppercase_symbol_in_subject(self):
        assert is_same_gene_hit("PREDICTED: ACTB transcript variant 1", "ACTB")

    def test_paren_wrapped_symbol(self):
        assert is_same_gene_hit("Homo sapiens beta-actin (ACTB), mRNA", "ACTB")

    def test_full_name_synonym_matches(self):
        # GSN <-> 'gelsolin'
        assert is_same_gene_hit("Mus musculus gelsolin (Gsn), mRNA", "GSN")

    def test_unrelated_gene_does_not_match(self):
        assert not is_same_gene_hit("Homo sapiens TLR1 transcript", "TLR6")

    def test_empty_inputs_safe(self):
        assert not is_same_gene_hit("", "ACTB")
        assert not is_same_gene_hit("ACTB", "")


# ---------------------------------------------------------------------------
# apply_gene_aware_filter — the rejection rule
# ---------------------------------------------------------------------------

class TestApplyGeneAwareFilter:
    @pytest.fixture
    def base_seq(self):
        return {
            "probe_id": "p1", "gene": "TLR6",
            "binding_sequence": "GATTACA",
            "original_item": {"Gene.refGene.x": "TLR6"},
            "original_probe": {"score": 0.5},
        }

    def test_no_alignments_passes(self, base_seq):
        passed, rejected = apply_gene_aware_filter([base_seq], {"p1": []})
        assert len(passed) == 1
        assert len(rejected) == 0
        assert passed[0]["blast_hits"] == 0
        assert passed[0]["best_match_percentage"] == 0

    def test_same_gene_exact_match_passes(self, base_seq):
        blast = {"p1": [{"subject_id": "Homo sapiens TLR6 transcript",
                          "match_percentage": 100.0}]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(passed) == 1
        assert passed[0]["same_gene_exact"] == 1

    def test_off_target_exact_match_rejects(self, base_seq):
        blast = {"p1": [{"subject_id": "Homo sapiens TLR1 transcript",
                          "match_percentage": 100.0}]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(rejected) == 1
        assert rejected[0]["reason"] == "off_target_exact_match"

    def test_two_high_similarity_off_target_passes(self, base_seq):
        # Threshold is "3 or more"; 2 should pass.
        blast = {"p1": [
            {"subject_id": "Homo sapiens FOO", "match_percentage": 96.0},
            {"subject_id": "Homo sapiens BAR", "match_percentage": 95.5},
        ]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(passed) == 1

    def test_three_high_similarity_off_target_rejects(self, base_seq):
        blast = {"p1": [
            {"subject_id": "Homo sapiens FOO", "match_percentage": 96.0},
            {"subject_id": "Homo sapiens BAR", "match_percentage": 95.5},
            {"subject_id": "Homo sapiens BAZ", "match_percentage": 95.1},
        ]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(rejected) == 1
        assert rejected[0]["reason"] == "multiple_off_target_high_similarity"


# ---------------------------------------------------------------------------
# find_best_per_mutation_record — group + min-score selection
# ---------------------------------------------------------------------------

class TestFindBestPerMutation:
    def test_groups_by_mutation_key_and_picks_lowest_score(self):
        seqs = [
            {"probe_id": "a", "gene": "TP53", "blast_hits": 0,
             "original_item": {"Gene.refGene.x": "TP53", "Start": 100, "End": 100,
                               "Ref": "C", "Alt": "T"},
             "original_probe": {"score": 5.0}},
            {"probe_id": "b", "gene": "TP53", "blast_hits": 0,
             "original_item": {"Gene.refGene.x": "TP53", "Start": 100, "End": 100,
                               "Ref": "C", "Alt": "T"},
             "original_probe": {"score": 2.0}},
            {"probe_id": "c", "gene": "TP53", "blast_hits": 0,
             "original_item": {"Gene.refGene.x": "TP53", "Start": 200, "End": 200,
                               "Ref": "G", "Alt": "A"},
             "original_probe": {"score": 1.5}},
        ]
        best = find_best_per_mutation_record(seqs)
        assert len(best) == 2
        # The one for Start=100 should have score 2.0 (lower than 5.0)
        scores = sorted(b["score"] for b in best)
        assert scores == [1.5, 2.0]


# ---------------------------------------------------------------------------
# parse_blast_results — XML parsing round-trip
# ---------------------------------------------------------------------------

MINIMAL_BLAST_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_hits>
        <Hit>
          <Hit_def>Homo sapiens TLR6 transcript variant 1</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_identity>40</Hsp_identity>
              <Hsp_align-len>40</Hsp_align-len>
              <Hsp_evalue>1e-10</Hsp_evalue>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def test_parse_blast_results_round_trip(tmp_path: Path):
    blast_dir = tmp_path
    (blast_dir / "batch_1.xml").write_text(MINIMAL_BLAST_XML, encoding="utf-8")
    seqs = [{"probe_id": "p1", "binding_sequence": "GATTACA"}]
    results = parse_blast_results(seqs, blast_dir, batch_size=100)
    assert "p1" in results
    assert len(results["p1"]) == 1
    hit = results["p1"][0]
    assert hit["match_percentage"] == 100.0
    assert hit["align_length"] == 40
    assert "TLR6" in hit["subject_id"]
