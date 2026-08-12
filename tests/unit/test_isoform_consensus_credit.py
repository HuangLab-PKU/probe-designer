"""``isoform_consensus`` crediting, end-to-end through IsoformAwareStrategy.

Companion to ``test_isoform_coverage.py``, which locks the pure geometry. This
file locks what the SEARCH STRATEGY does with it — see
``experiments/20260727_CRC_BRAF_round1/DESIGNER_ISSUE_isoform_coverage.md``:

* an isoform is credited only when the ligation junction is paired and both
  effective arms clear ``min_isoform_arm_nt``;
* a selected site must bind at least one protein-coding transcript, so a
  candidate generated on a minor / non-coding isoform cannot ship on its own;
* the isoform a candidate is BUILT on is the canonical one where there is a
  choice, so the position-dedup keeps the protein-coding host.

Synthetic locus (plus strand, 1000..1399), four transcripts:

    G-201  protein_coding      [1000..1199] [1300..1399]
    G-204  protein_coding      [1000..1100] [1150..1199] [1300..1399]
    G-202  retained_intron     [1000..1399]                 (one long exon)
    G-203  processed_transcript             [1200..1299]

G-204's exon boundary at 1100 sits inside a 40-nt window starting at 1080 —
the synthetic form of the real FABP1-201 straddler. 1200..1299 is intronic in
both protein-coding transcripts — the synthetic form of the GNLY / MZB1 sites
that bind no coding transcript at all.
"""
from __future__ import annotations

import pytest

from probe_designer.config import FilterConfig, SearchConfig
from probe_designer.search_strategies import IsoformAwareness, IsoformAwareStrategy


LOCUS_START, LOCUS_END = 1000, 1399
BDS_LEN = 40


def _iso(name, exons, biotype):
    return {
        "id": f"ENST_{name}",
        "display_name": name,
        "biotype": biotype,
        "start": min(s for s, _ in exons),
        "end": max(e for _, e in exons),
        "strand": 1,
        "seq_region_name": "chrT",
        "Exon": [{"start": s, "end": e} for s, e in exons],
    }


# Deliberately minor-first: the pre-fix code took whichever isoform came first
# as the host, so input order silently decided which transcript a shipped probe
# was designed against.
ISOFORMS = [
    _iso("G-202", [(1000, 1399)], "retained_intron"),
    _iso("G-203", [(1200, 1299)], "processed_transcript"),
    _iso("G-201", [(1000, 1199), (1300, 1399)], "protein_coding"),
    _iso("G-204", [(1000, 1100), (1150, 1199), (1300, 1399)], "protein_coding"),
]

NONCODING_ONLY = [ISOFORMS[0], ISOFORMS[1]]

# Immune-receptor constant regions carry their own Ensembl biotypes, never
# `protein_coding`: TRAC/TRBC1 are TR_C_gene, IGHG1 is IG_C_gene. They do code
# for protein. The 2026-08-02 bank audit found six TRAC/TRBC1 sites flagged by
# a rule that only looked for `protein_coding` — the rule was wrong, not the
# sites.
IMMUNE_RECEPTOR = [
    _iso("TRAC-201", [(1000, 1199), (1300, 1399)], "TR_C_gene"),
    _iso("TRAC-202", [(1200, 1299)], "processed_transcript"),
]

BIOTYPE = {iso["display_name"]: iso["biotype"] for iso in ISOFORMS}


@pytest.fixture
def fake_genome():
    """400 nt of repeating ACGT — 50% GC, no homopolymers, so every window
    clears the thermal filter and the test observes crediting alone."""
    seq = "ACGT" * ((LOCUS_END - LOCUS_START + 1) // 4)

    def accessor(region_name, start, end):
        assert region_name == "chrT"
        return seq[start - LOCUS_START: end - LOCUS_START + 1]

    return accessor


def _search(fake_genome, isoforms=ISOFORMS, **search_kwargs):
    search_cfg = SearchConfig(binding_site_length=BDS_LEN, **search_kwargs)
    strategy = IsoformAwareStrategy(
        search_cfg, FilterConfig(), isoforms,
        genome_accessor=fake_genome, mode="consensus",
    )
    return strategy.search_binding_sites("", "G")


def _at(sites, st):
    matches = [s for s in sites if s["st"] == st]
    return matches[0] if matches else None


# ---------------------------------------------------------------------------
# The count is now bounded and means something
# ---------------------------------------------------------------------------

class TestCreditedCount:
    def test_overlap_num_never_exceeds_the_isoform_count(self, fake_genome):
        sites = _search(fake_genome)
        assert sites
        assert max(s["isoform_overlap_num"] for s in sites) <= len(ISOFORMS)

    def test_overlap_num_matches_the_credited_list(self, fake_genome):
        for site in _search(fake_genome):
            assert site["isoform_overlap_num"] == len(site["target_isoforms"])

    def test_host_isoform_is_always_credited_to_itself(self, fake_genome):
        for site in _search(fake_genome):
            assert site["isoform_name"] in site["target_isoforms"]


# ---------------------------------------------------------------------------
# The straddler — the defect that shipped
# ---------------------------------------------------------------------------

class TestBoundaryStraddler:
    """A window at 1080..1119 covers G-204 only up to its exon end at 1100,
    leaving 1 nt above the ligation junction. It must not be credited."""

    def test_old_union_rule_credited_the_straddler(self, fake_genome):
        # Documents the defect this file exists for: the pre-fix code took the
        # union over every atomic splice region the window TOUCHED.
        awareness = IsoformAwareness(ISOFORMS, fake_genome)
        regions = awareness.get_overlapping_regions(1080, 1119)
        assert "G-204" in awareness.isoforms_union(regions)

    def test_straddler_is_no_longer_credited(self, fake_genome):
        site = _at(_search(fake_genome), 1080)
        assert site is not None
        assert "G-204" not in site["target_isoforms"]

    def test_isoforms_that_do_contain_the_window_are_still_credited(self, fake_genome):
        site = _at(_search(fake_genome), 1080)
        assert set(site["target_isoforms"]) == {"G-201", "G-202"}

    def test_a_lenient_bound_lets_the_straddler_back_in(self, fake_genome):
        # 1 nt above the junction, so only a bound of 1 or 0 admits it. Locks
        # that the knob is what decides, not a hard-coded number.
        site = _at(_search(fake_genome, min_isoform_arm_nt=1), 1080)
        assert "G-204" in site["target_isoforms"]


# ---------------------------------------------------------------------------
# The protein-coding requirement — the 6 dead CRC sites
# ---------------------------------------------------------------------------

class TestProteinCodingRequirement:
    """1200..1299 is intronic in both coding transcripts. Candidates there are
    credited only to G-202 (retained intron) and G-203 (processed transcript),
    which is the shape of GNLY|0, B2M|0, VIM|1 and MZB1|0."""

    def test_site_binding_no_coding_transcript_is_dropped(self, fake_genome):
        assert _at(_search(fake_genome), 1220) is None

    def test_every_surviving_site_binds_a_coding_transcript(self, fake_genome):
        for site in _search(fake_genome):
            credited = site["target_isoforms"]
            assert any(BIOTYPE[name] == "protein_coding" for name in credited)

    def test_the_requirement_can_be_switched_off(self, fake_genome):
        site = _at(_search(fake_genome, require_protein_coding=False), 1220)
        assert site is not None
        assert set(site["target_isoforms"]) == {"G-202", "G-203"}

    def test_gene_without_any_coding_transcript_still_produces_sites(self, fake_genome):
        # A lncRNA-only gene must not be wiped out by a requirement it can
        # never satisfy; the run falls back and says so.
        with pytest.warns(UserWarning, match="protein_coding"):
            sites = _search(fake_genome, isoforms=NONCODING_ONLY)
        assert sites

    def test_immune_receptor_constant_region_counts_as_coding(self, fake_genome):
        # TR_C_gene / IG_C_gene transcripts are productive; a TRAC probe must
        # survive the requirement, and without a fallback warning — the gene
        # is not exempt, it genuinely has a coding transcript.
        import warnings as _warnings
        with _warnings.catch_warnings():
            _warnings.simplefilter("error", UserWarning)
            sites = _search(fake_genome, isoforms=IMMUNE_RECEPTOR)
        assert sites
        assert all("TRAC-201" in s["target_isoforms"] for s in sites)

    def test_immune_receptor_still_drops_its_intronic_sites(self, fake_genome):
        # The rule must still bite: 1220..1259 lies in TRAC-201's intron and is
        # carried only by the processed_transcript.
        sites = _search(fake_genome, isoforms=IMMUNE_RECEPTOR)
        assert not [s for s in sites if s["st"] == 1220]

    def test_distal_shortfall_inside_the_bound_is_still_credited(self, fake_genome):
        # 1161..1200 runs 1 nt past G-201's exon end at 1199: 20 nt below the
        # junction, 19 above. Above the bound, so it counts — this is the
        # tolerance the lab rule deliberately allows.
        site = _at(_search(fake_genome), 1161)
        assert site is not None
        assert "G-201" in site["target_isoforms"]


# ---------------------------------------------------------------------------
# What the site record must carry
# ---------------------------------------------------------------------------

class TestSiteRecord:
    def test_records_per_isoform_effective_arm_lengths(self, fake_genome):
        site = _at(_search(fake_genome), 1080)
        assert site["isoform_coverage"]["G-201"] == {
            "arm_3prime_nt": 20, "arm_5prime_nt": 20,
        }

    def test_coverage_keys_match_the_credited_isoforms(self, fake_genome):
        site = _at(_search(fake_genome), 1080)
        assert set(site["isoform_coverage"]) == set(site["target_isoforms"])

    def test_records_how_many_coding_transcripts_are_bound(self, fake_genome):
        site = _at(_search(fake_genome), 1080)
        assert site["protein_coding_overlap_num"] == 1

    def test_record_is_json_serializable(self, fake_genome):
        import json
        json.dumps(_search(fake_genome))


class TestHostIsoformChoice:
    """Position-dedup keeps whichever isoform reached the window first, so the
    iteration order decides which transcript a shipped probe was designed
    against. It must be a coding one wherever there is a choice."""

    def test_host_is_protein_coding_despite_minor_first_input(self, fake_genome):
        site = _at(_search(fake_genome), 1000)
        assert BIOTYPE[site["isoform_name"]] == "protein_coding"

    def test_host_is_coding_wherever_a_coding_transcript_can_host(self, fake_genome):
        # Only an isoform whose own exon fits the whole window can host it, so
        # this is asserted at positions that sit inside an exon of G-201/G-204
        # as well as inside G-202's long one.
        sites = _search(fake_genome)
        for st in (1000, 1060, 1120):
            assert BIOTYPE[_at(sites, st)["isoform_name"]] == "protein_coding"

    def test_a_window_only_a_minor_isoform_can_host_keeps_that_host(self, fake_genome):
        # 1161..1200 runs past every coding exon, so G-202 is the only
        # transcript that can generate it. Re-homing it would be a lie about
        # where the sequence came from.
        assert _at(_search(fake_genome), 1161)["isoform_name"] == "G-202"


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

class TestConfiguration:
    def test_defaults_are_on(self):
        cfg = SearchConfig()
        assert cfg.min_isoform_arm_nt == 16
        assert cfg.require_protein_coding is True

    def test_bound_wider_than_half_the_window_fails_fast(self, fake_genome):
        # 21 nt on each side of the junction is unsatisfiable for a 40-nt
        # window — not even the host qualifies. Refuse rather than return zero
        # sites for every gene.
        with pytest.raises(ValueError, match="min_isoform_arm_nt"):
            _search(fake_genome, min_isoform_arm_nt=21)
