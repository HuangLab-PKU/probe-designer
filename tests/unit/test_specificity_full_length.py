"""The specificity filter must check how MUCH of the probe aligned.

From ``DESIGNER_ISSUE_isoform_coverage.md`` section 5: ``_check_specificity_criteria``
accepted a site as soon as it found a same-gene alignment with identity > 0.95,
and never looked at the alignment length. So a probe that matches its own gene
over only 32 of its 40 nt — the signature of a site designed on a minor isoform
and read against the canonical transcript — passed exactly like a 40/40 one.
That is why BLAST did not catch any of the six dead CRC sites.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from probe_designer.config import BlastConfig, FilterConfig
from probe_designer.filtering import SequenceFilter


PROBE = "CTCTTGAAGACAATGTCACCCAATGTCATGGTCTGAAAGC"  # real FABP1 site 0, 40 nt


def _filter(**filter_kwargs) -> SequenceFilter:
    cfg = FilterConfig(target_organisms=["Homo sapiens"], **filter_kwargs)
    return SequenceFilter(cfg, BlastConfig())


def _site(gene="FABP1", sequence=PROBE):
    return {"gene_name": gene, "sequence": sequence}


def _alignment(
    align_length,
    *,
    hit_def="Homo sapiens FABP1 mRNA",
    identity=1.0,
    mismatch_at=None,
    frame=(1, -1),
    query_start=1,
):
    """One parsed HSP, as ``_parse_single_batch_file`` now emits it.

    ``align_length`` columns starting at ``query_start`` (1-based), every one
    paired unless ``mismatch_at`` -- a 0-based position in the probe -- says
    otherwise. The midline is what lets the filter ask WHERE a mismatch fell
    rather than only how many there were.
    """
    begin = query_start - 1
    midline = ["|"] * align_length
    if mismatch_at is not None:
        midline[mismatch_at - begin] = " "
    return {
        "hit_id": "ENST00000393750",
        "hit_def": hit_def,
        "evalue": 1e-13,
        "identity": identity,
        "align_length": align_length,
        "identities": round(identity * align_length) - (mismatch_at is not None),
        "score": 80,
        "frame": frame,
        "midline": "".join(midline),
        "query_seq": PROBE[begin:begin + align_length],
        "query_start": query_start,
        "query_from": query_start,
        "query_to": begin + align_length,
    }


def _blast_info(*alignments):
    return {
        "alignments": list(alignments),
        "evalue": 1e-13,
        "identity": 1.0,
        "alignment_count": len(alignments),
    }


class TestSameGeneCoverage:
    def test_full_length_same_gene_hit_passes(self):
        assert _filter()._check_specificity_criteria(
            _blast_info(_alignment(40)), _site()
        )

    def test_partial_same_gene_hit_is_rejected(self):
        # The FABP1-201 case: 32 of 40 nt align because the rest is spliced out.
        assert not _filter()._check_specificity_criteria(
            _blast_info(_alignment(32)), _site()
        )

    def test_one_full_length_hit_among_partials_is_enough(self):
        # A probe designed on the canonical transcript legitimately gives short
        # hits on minor isoforms; it only has to be full length SOMEWHERE.
        assert _filter()._check_specificity_criteria(
            _blast_info(_alignment(32), _alignment(40)), _site()
        )

    def test_one_nt_short_still_passes_the_default_tolerance(self):
        # 39/40 = 0.975 >= the 0.95 default: BLAST routinely trims a terminal
        # mismatch, and that is not the failure mode this gate is for.
        assert _filter()._check_specificity_criteria(
            _blast_info(_alignment(39)), _site()
        )

    def test_threshold_is_configurable(self):
        assert _filter(min_same_gene_coverage=0.5)._check_specificity_criteria(
            _blast_info(_alignment(32)), _site()
        )

    def test_alignment_without_a_length_is_not_assumed_full(self):
        no_length = _alignment(40)
        del no_length["align_length"]
        assert not _filter()._check_specificity_criteria(
            _blast_info(no_length), _site()
        )


class TestUnchangedRules:
    """The length check is additive — the pre-existing rules still decide."""

    def test_off_target_gene_hit_still_rejects(self):
        off_target = _alignment(40, hit_def="Homo sapiens TLR6 mRNA")
        assert not _filter()._check_specificity_criteria(
            _blast_info(_alignment(40), off_target), _site()
        )

    def test_no_target_organism_hit_still_rejects(self):
        assert not _filter()._check_specificity_criteria(
            _blast_info(_alignment(40, hit_def="Mus musculus Fabp1 mRNA")), _site()
        )

    def test_disabling_specificity_bypasses_everything(self):
        cfg = FilterConfig(target_organisms=["Homo sapiens"], require_specificity=False)
        filt = SequenceFilter(cfg, BlastConfig())
        assert filt._check_specificity_criteria(_blast_info(_alignment(32)), _site())


# ---------------------------------------------------------------------------
# The parser has to record the length in the first place
# ---------------------------------------------------------------------------

_BLAST_XML = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.17.0+</BlastOutput_version>
  <BlastOutput_reference>ref</BlastOutput_reference>
  <BlastOutput_db>human_gencode_db</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>FABP1|st=88123073|en=88123112|g_content=0.20</BlastOutput_query-def>
  <BlastOutput_query-len>40</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>1000</Parameters_expect>
      <Parameters_sc-match>1</Parameters_sc-match>
      <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
      <Parameters_gap-open>0</Parameters_gap-open>
      <Parameters_gap-extend>0</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>FABP1|st=88123073|en=88123112|g_content=0.20</Iteration_query-def>
      <Iteration_query-len>40</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>ENST00000295834</Hit_id>
          <Hit_def>Homo sapiens FABP1-201</Hit_def>
          <Hit_accession>ENST00000295834</Hit_accession>
          <Hit_len>800</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>64.0</Hsp_bit-score>
              <Hsp_score>32</Hsp_score>
              <Hsp_evalue>1e-10</Hsp_evalue>
              <Hsp_query-from>9</Hsp_query-from>
              <Hsp_query-to>40</Hsp_query-to>
              <Hsp_hit-from>200</Hsp_hit-from>
              <Hsp_hit-to>169</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>-1</Hsp_hit-frame>
              <Hsp_identity>32</Hsp_identity>
              <Hsp_positive>32</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>32</Hsp_align-len>
              <Hsp_qseq>CAATGTCACCCAATGTCATGGTCTGAAAGC</Hsp_qseq>
              <Hsp_hseq>CAATGTCACCCAATGTCATGGTCTGAAAGC</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


class TestParserRecordsAlignLength:
    def test_align_length_survives_the_xml_parse(self, tmp_path: Path):
        xml = tmp_path / "blast_results.xml"
        xml.write_text(_BLAST_XML, encoding="utf-8")

        parsed = _filter()._parse_single_batch_file(str(xml), 1)

        (info,) = parsed.values()
        assert info["alignments"][0]["align_length"] == 32

    def test_query_length_is_recorded_too(self, tmp_path: Path):
        # Needed to express coverage as a fraction without trusting the caller
        # to hand back the same sequence that was submitted.
        xml = tmp_path / "blast_results.xml"
        xml.write_text(_BLAST_XML, encoding="utf-8")

        parsed = _filter()._parse_single_batch_file(str(xml), 1)

        (info,) = parsed.values()
        assert info["query_length"] == 40

    def test_a_real_partial_hit_is_then_rejected(self, tmp_path: Path):
        # End to end: parse -> criteria. 32/40 on its own gene must not pass.
        xml = tmp_path / "blast_results.xml"
        xml.write_text(_BLAST_XML, encoding="utf-8")

        filt = _filter()
        (info,) = filt._parse_single_batch_file(str(xml), 1).values()

        assert not filt._check_specificity_criteria(info, _site())


class TestJunctionPairing:
    """Coverage says how much aligned; it cannot say WHERE the gap fell.

    Both alignments below are 39/40 over the probe's full length. One is the
    terminal mismatch BLAST routinely trims and is harmless; the other puts the
    unpaired base on the ligation junction, where SplintR will not close the
    nick. The old gate scored both 0.975 and passed both.
    """

    def test_mismatch_on_the_nick_is_rejected(self):
        assert not _filter()._check_specificity_criteria(
            _blast_info(_alignment(40, mismatch_at=19)), _site()
        )

    def test_mismatch_at_the_far_end_still_passes(self):
        assert _filter()._check_specificity_criteria(
            _blast_info(_alignment(40, mismatch_at=39)), _site()
        )

    def test_mismatch_just_outside_the_clamp_passes(self):
        # +/-3 is the SplintR clamp width; position 16 is one base beyond it.
        assert _filter()._check_specificity_criteria(
            _blast_info(_alignment(40, mismatch_at=16)), _site()
        )

    def test_mismatch_inside_the_clamp_is_rejected(self):
        assert not _filter()._check_specificity_criteria(
            _blast_info(_alignment(40, mismatch_at=17)), _site()
        )

    def test_alignment_stopping_before_the_nick_is_rejected(self):
        assert not _filter()._check_specificity_criteria(
            _blast_info(_alignment(15)), _site()
        )

    def test_missing_midline_is_not_credited(self):
        # Same stance as a missing align_length: an absent measurement is not
        # a pass. It is also the shape a silently-inert check would take.
        alignment = _alignment(40)
        alignment["midline"] = ""
        assert not _filter()._check_specificity_criteria(
            _blast_info(alignment), _site()
        )

    def test_gate_can_be_turned_off(self):
        assert _filter(require_junction_pairing=False)._check_specificity_criteria(
            _blast_info(_alignment(40, mismatch_at=19)), _site()
        )

    def test_vicinity_width_is_configurable(self):
        # n=1 clamps only the two bases flanking the nick, so 17 clears it.
        assert _filter(junction_vicinity_nt=1)._check_specificity_criteria(
            _blast_info(_alignment(40, mismatch_at=17)), _site()
        )


class TestBindingStrand:
    """An antisense transcript carries the probe's own sequence, not its
    complement, so the probe cannot bind it -- for either chemistry. Easy to
    miss, because an antisense locus really is a different gene: the same-gene
    test says "different" and the rejection looks right. The SP369 audit found
    30 such hits among 4299 full-length identities, enough to move the
    non-specific set from 40 to 61.
    """

    def _antisense(self):
        return _alignment(40, hit_def="Homo sapiens KRT7-AS RNA", frame=(1, 1))

    def _sense_off_target(self):
        return _alignment(40, hit_def="Homo sapiens SPINK1 mRNA", frame=(1, -1))

    def test_antisense_hit_to_another_gene_is_not_an_off_target(self):
        assert _filter()._check_specificity_criteria(
            _blast_info(_alignment(40), self._antisense()), _site()
        )

    def test_bindable_off_target_still_rejects(self):
        assert not _filter()._check_specificity_criteria(
            _blast_info(_alignment(40), self._sense_off_target()), _site()
        )

    def test_unreported_frame_is_still_counted_against_the_probe(self):
        assert not _filter()._check_specificity_criteria(
            _blast_info(
                _alignment(40),
                _alignment(40, hit_def="Homo sapiens SPINK1 mRNA", frame=()),
            ),
            _site(),
        )

    def test_sense_design_flips_which_strand_binds(self):
        # A cDNA-style design stores the site mRNA-sense, so its real hits are
        # plus-strand and the minus-strand hit is the one it cannot bind.
        cfg = _filter(binding_hit_frame=1)
        assert cfg._check_specificity_criteria(
            _blast_info(
                _alignment(40, frame=(1, 1)),
                self._sense_off_target(),
            ),
            _site(),
        )


class TestNickFromTheRealArm:
    """The nick is where the arms actually meet, not necessarily the midpoint.

    The mRNA strategies split at ``bds_len // 2`` (``search_strategies``), so
    the two agree today. Reading ``arm_3prime`` keeps the gate aimed at the real
    junction if that ever changes -- the gene-aware path already does.
    """

    def _unbalanced_site(self):
        site = _site()
        site['arm_3prime'] = PROBE[:22]
        site['arm_5prime'] = PROBE[22:]
        return site

    def test_mismatch_at_the_real_nick_is_rejected(self):
        assert not _filter()._check_specificity_criteria(
            _blast_info(_alignment(40, mismatch_at=22)), self._unbalanced_site()
        )

    def test_mismatch_in_the_midpoint_clamp_only_now_passes(self):
        # Position 17 sits inside the midpoint clamp [17,22] but 5 nt away from
        # the real nick at 22, whose clamp is [19,24]. Judged on the arm, this
        # probe is fine; judged on the midpoint it was rejected.
        assert _filter()._check_specificity_criteria(
            _blast_info(_alignment(40, mismatch_at=17)), self._unbalanced_site()
        )

    def test_falls_back_to_the_midpoint_without_arms(self):
        assert not _filter()._check_specificity_criteria(
            _blast_info(_alignment(40, mismatch_at=19)), _site()
        )


class TestLigationAwareOffTargets:
    """Off-target hits are judged the way the gene-aware path judges them.

    The old rule rejected a site on ANY different-gene hit on the binding strand,
    so a 15/15 perfect match inside one arm -- which neither holds at the
    hybridisation temperature nor ligates -- threw away a good site.
    """

    OFF = "Homo sapiens SPINK1 mRNA"

    def _info(self, *off_targets):
        return _blast_info(_alignment(40), *off_targets)

    def test_short_hit_inside_one_arm_no_longer_rejects(self):
        assert _filter()._check_specificity_criteria(
            self._info(_alignment(15, hit_def=self.OFF)), _site()
        )

    def test_full_length_off_target_still_rejects(self):
        assert not _filter()._check_specificity_criteria(
            self._info(_alignment(40, hit_def=self.OFF)), _site()
        )

    def test_long_off_target_broken_at_the_nick_only_binds(self):
        # 39/40 is credible, but the mismatch sits inside the clamp, so no
        # ligase closes it. A single binder is tolerated.
        assert _filter()._check_specificity_criteria(
            self._info(_alignment(40, hit_def=self.OFF, mismatch_at=19)), _site()
        )

    def test_binders_that_miss_the_nick_reject_only_in_numbers(self):
        binder = _alignment(20, hit_def=self.OFF, query_start=21)  # arm5 only
        assert _filter()._check_specificity_criteria(
            self._info(binder, binder), _site()
        )
        assert not _filter()._check_specificity_criteria(
            self._info(binder, binder, binder), _site()
        )

    def test_hit_without_a_midline_is_placed_by_its_coordinates(self):
        hit = _alignment(40, hit_def=self.OFF)
        hit["midline"] = ""
        assert not _filter()._check_specificity_criteria(self._info(hit), _site())

    def test_flag_off_restores_the_old_rule(self):
        assert not _filter(ligation_aware_off_targets=False)._check_specificity_criteria(
            self._info(_alignment(15, hit_def=self.OFF)), _site()
        )
