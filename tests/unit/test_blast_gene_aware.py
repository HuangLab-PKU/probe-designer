"""Tests for probe_designer.blast.gene_aware (pure helpers + XML parsing)."""
from __future__ import annotations

from pathlib import Path

import pytest

from probe_designer.blast.gene_aware import (
    apply_gene_aware_filter,
    binding_query_junction_index,
    extract_binding_query,
    find_best_per_mutation_record,
    hsp_spans_junction,
    is_same_gene_hit,
    parse_blast_results,
)


# ---------------------------------------------------------------------------
# extract_binding_query — the iLock dedup correctness
# ---------------------------------------------------------------------------

class TestExtractBindingQuery:
    def test_non_ilock_concatenates_arm3_then_arm5(self):
        q = extract_binding_query("AAAA", "TTTT", "dRNA")
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

    def test_mixed_case_ortholog_symbol_matches(self):
        # NCBI names rodent p53 "Tp53" and the human read-through locus
        # "C12orf57"; case-insensitive matching must recognise both as same-gene.
        assert is_same_gene_hit(
            "PREDICTED: Marmota flaviventris tumor protein p53 (Tp53), mRNA", "TP53")
        assert is_same_gene_hit(
            "Homo sapiens chromosome 12 open reading frame 57 (C12orf57), mRNA",
            "C12orf57")

    def test_empty_inputs_safe(self):
        assert not is_same_gene_hit("", "ACTB")
        assert not is_same_gene_hit("ACTB", "")


# ---------------------------------------------------------------------------
# apply_gene_aware_filter — the rejection rule
# ---------------------------------------------------------------------------

ARM3_20 = "TTGCCAGATCGTACGATCCA"      # 20 nt -> junction at query index 20
ARM5_20 = "ACGTGCAGTCATGGCTAACG"      # query is arm3 + arm5, 40 nt total


def _hit(subject, pct, *, q_from, q_to, align=None):
    """One BLAST HSP. Query coords are 1-based inclusive, as BLAST reports."""
    return {
        "subject_id": subject,
        "match_percentage": pct,
        "align_length": align if align is not None else (q_to - q_from + 1),
        "identities": 0,
        "evalue": 1e-9,
        "query_from": q_from,
        "query_to": q_to,
    }


#: Straddles the junction (query index 20) with room on both sides.
def _junction_hit(subject, pct=100.0):
    return _hit(subject, pct, q_from=1, q_to=40)


#: Sits wholly inside the arm3 half — a ligase can never close on it.
def _arm_internal_hit(subject, pct=100.0):
    return _hit(subject, pct, q_from=3, q_to=17)


class TestJunctionHelpers:
    def test_junction_index_is_arm3_length_for_plain_chemistries(self):
        assert binding_query_junction_index(ARM3_20) == 20

    def test_junction_index_matches_the_ilock_query_layout(self):
        """iLock's query drops the duplicated base, so arm3 still fills [0, 20)."""
        arm3, arm5 = "TTGCCAGATCGTACGATCCG", "GCGTGCAGTCATGGCTAACG"  # share 'G'
        assert arm5[0] == arm3[-1]
        query = extract_binding_query(arm5, arm3, "iLock")
        idx = binding_query_junction_index(arm3)
        assert query[:idx] == arm3.upper()
        assert query[idx:] == arm5[1:].upper()

    def test_span_requires_clamp_on_both_sides(self):
        assert hsp_spans_junction(1, 40, 20, 3)
        assert hsp_spans_junction(18, 23, 20, 3)          # exactly 3 each side
        assert not hsp_spans_junction(19, 23, 20, 3)      # 1 short on the arm3 side
        assert not hsp_spans_junction(18, 22, 20, 3)      # 1 short on the arm5 side

    def test_span_false_for_a_hit_entirely_on_one_side(self):
        assert not hsp_spans_junction(1, 15, 20, 3)
        assert not hsp_spans_junction(25, 40, 20, 3)

    def test_span_false_without_coordinates(self):
        assert not hsp_spans_junction(0, 0, 20, 3)


class TestApplyGeneAwareFilter:
    @pytest.fixture
    def base_seq(self):
        return {
            "probe_id": "p1", "gene": "TLR6",
            "binding_sequence": (ARM3_20 + ARM5_20),
            "arm3": ARM3_20, "arm5": ARM5_20,
            "original_item": {"Gene.refGene.x": "TLR6"},
            "original_probe": {"score": 0.5},
        }

    def test_no_alignments_passes(self, base_seq):
        passed, rejected = apply_gene_aware_filter([base_seq], {"p1": []})
        assert len(passed) == 1 and not rejected
        assert passed[0]["blast_hits"] == 0
        assert passed[0]["best_match_percentage"] == 0

    def test_same_gene_junction_hit_passes(self, base_seq):
        blast = {"p1": [_junction_hit("Homo sapiens TLR6 transcript")]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(passed) == 1
        assert passed[0]["same_gene_exact"] == 1

    def test_off_target_junction_spanning_perfect_hit_rejects(self, base_seq):
        blast = {"p1": [_junction_hit("Homo sapiens TLR1 transcript")]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(rejected) == 1
        assert rejected[0]["reason"] == "off_target_ligation_competent"

    def test_single_imperfect_junction_hit_rejects(self, base_seq):
        """39/40 across the junction: the eLife 107070 failure mode.

        Under the old identity-only rule this scored 97.5%, not 100%, so it took
        three such hits to reject — while being exactly the case that produces
        signal indistinguishable from the real target.
        """
        blast = {"p1": [_junction_hit("Homo sapiens APOBEC3D", pct=97.5)]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(rejected) == 1
        assert rejected[0]["reason"] == "off_target_ligation_competent"

    def test_short_perfect_hit_inside_one_arm_passes(self, base_seq):
        """15/15 within arm3: the old rule's false positive.

        100% identity over a fragment the ligase can never close on. It rejected
        a good probe because identity was measured over the HSP, not the
        footprint.
        """
        blast = {"p1": [_arm_internal_hit("Homo sapiens TLR1 transcript")]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(passed) == 1, rejected
        assert passed[0]["off_target_binding_hits"] == 0   # too short to count

    def test_long_off_target_binder_that_misses_the_junction_is_not_fatal(self, base_seq):
        """It sequesters probe but cannot miscall, so one is tolerated."""
        blast = {"p1": [_hit("Homo sapiens FOO", 100.0, q_from=1, q_to=21)]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(passed) == 1
        assert passed[0]["off_target_binding_hits"] == 1

    def test_two_off_target_binders_pass(self, base_seq):
        blast = {"p1": [
            _hit("Homo sapiens FOO", 96.0, q_from=1, q_to=21),
            _hit("Homo sapiens BAR", 95.5, q_from=1, q_to=22),
        ]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(passed) == 1

    def test_three_off_target_binders_reject_as_promiscuous(self, base_seq):
        # All stop at query 22 or earlier: with junction 20 and clamp 3 a hit
        # must reach 23 to span, so none of these is ligation-competent.
        blast = {"p1": [
            _hit("Homo sapiens FOO", 96.0, q_from=1, q_to=20),
            _hit("Homo sapiens BAR", 95.5, q_from=1, q_to=21),
            _hit("Homo sapiens BAZ", 95.1, q_from=1, q_to=22),
        ]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(rejected) == 1
        assert rejected[0]["reason"] == "multiple_off_target_high_similarity"

    def test_junction_hit_below_the_length_floor_does_not_reject(self, base_seq):
        """Spans the junction, but too short to hold a duplex at the anneal temp."""
        blast = {"p1": [_hit("Homo sapiens TLR1", 100.0, q_from=18, q_to=23)]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(passed) == 1

    def test_low_identity_hit_is_ignored(self, base_seq):
        blast = {"p1": [_junction_hit("Homo sapiens TLR1", pct=80.0)]}
        passed, rejected = apply_gene_aware_filter([base_seq], blast)
        assert len(passed) == 1

    def test_record_without_a_locatable_junction_raises(self):
        """Silently degrading to the identity-only verdict is how this bug lived."""
        naked = {"probe_id": "p1", "gene": "TLR6"}
        with pytest.raises(ValueError, match="junction"):
            apply_gene_aware_filter([naked], {"p1": []})

    def test_explicit_junction_index_is_honoured(self, base_seq):
        seq = {**base_seq}
        seq.pop("arm3")
        seq["junction_index"] = 20
        blast = {"p1": [_junction_hit("Homo sapiens TLR1 transcript")]}
        _, rejected = apply_gene_aware_filter([seq], blast)
        assert rejected[0]["reason"] == "off_target_ligation_competent"


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
