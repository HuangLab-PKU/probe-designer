"""Contract tests for the unified probe-output schema.

These cover the constants, helpers, and validation in
``probe_designer.io.probe_schema``. Pipeline-level integration tests
(mutation/TCR/mRNA) live alongside their respective pipeline tests.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from probe_designer.io import probe_schema as ps


# ---------------------------------------------------------------------------
# chemistry constants
# ---------------------------------------------------------------------------


def test_chemistry_constants_have_expected_values():
    assert ps.CHEM_DRNA == "dRNA"
    assert ps.CHEM_CDNA == "cDNA"
    assert ps.CHEM_ILOCK == "iLock"
    assert ps.CHEM_RT == "RT_primer"


def test_padlock_chemistries_set_excludes_rt():
    assert ps.PADLOCK_CHEMISTRIES == {"dRNA", "cDNA", "iLock"}


# ---------------------------------------------------------------------------
# canonical column lists — first three columns + presence of key fields
# ---------------------------------------------------------------------------


def test_final_padlock_columns_lead_with_canonical_three():
    assert ps.FINAL_PADLOCK_COLUMNS[:3] == ["order", "probe_name", "probe_sequence"]


def test_final_rt_primer_columns_lead_with_canonical_three():
    assert ps.FINAL_RT_PRIMER_COLUMNS[:3] == ["order", "probe_name", "probe_sequence"]


def test_final_padlock_columns_include_renamed_arms():
    cols = ps.FINAL_PADLOCK_COLUMNS
    assert "probe_arm5" in cols and "probe_arm3" in cols
    assert "arm_5prime" not in cols and "arm_3prime" not in cols


def test_final_padlock_columns_keep_no_dot_and_drop_legacy_aliases():
    cols = ps.FINAL_PADLOCK_COLUMNS
    assert "No." in cols
    # iLock is no longer a separate column — it lives inside `chemistry`.
    assert "iLock" not in cols
    # Legacy capitalised aliases are gone.
    assert "Probe_Name" not in cols
    assert "Probe_Seq" not in cols
    assert "GC_Percent" not in cols
    assert "backbone_No." not in cols
    # Schema-v2 (Phase 0, 2026-05-14): gc_content (G+C) is the primary metric
    # but g_content (G-only) is preserved alongside it as a legacy reporting
    # field. The two are distinct measurements; do not collapse them.
    assert "gc_content" in cols
    assert "g_content" in cols


def test_final_padlock_columns_have_no_duplicates():
    assert len(ps.FINAL_PADLOCK_COLUMNS) == len(set(ps.FINAL_PADLOCK_COLUMNS))


def test_final_rt_primer_columns_have_no_duplicates():
    assert len(ps.FINAL_RT_PRIMER_COLUMNS) == len(set(ps.FINAL_RT_PRIMER_COLUMNS))


# ---------------------------------------------------------------------------
# make_padlock_name — format `{name}_{pos}_{chem}_{codebook}_{no}`
# ---------------------------------------------------------------------------


def test_make_padlock_name_basic_drna():
    name = ps.make_padlock_name(
        name="Actb", position=213, chemistry="dRNA", codebook="SP369", no=300,
    )
    assert name == "Actb_213_dRNA_SP369_300"


def test_make_padlock_name_ilock_chemistry():
    name = ps.make_padlock_name(
        name="BRCA1", position=12345, chemistry="iLock", codebook="SP369", no=42,
    )
    assert name == "BRCA1_12345_iLock_SP369_42"


def test_make_padlock_name_clone_id_with_underscores():
    # TCR-style identifier already contains underscores — the format should
    # still parse unambiguously because position is the first numeric chunk
    # after the trailing underscore.
    name = ps.make_padlock_name(
        name="BZ23_clone23_TRB", position=45, chemistry="cDNA",
        codebook="SP500", no=17,
    )
    assert name == "BZ23_clone23_TRB_45_cDNA_SP500_17"


def test_make_padlock_name_rejects_unknown_chemistry():
    with pytest.raises(ValueError, match="chemistry"):
        ps.make_padlock_name(
            name="X", position=1, chemistry="mRNA", codebook="SP1", no=1,
        )


def test_make_padlock_name_rejects_rt_primer_chemistry():
    # RT_primer must go through make_rt_primer_name, never the padlock helper.
    with pytest.raises(ValueError):
        ps.make_padlock_name(
            name="X", position=1, chemistry="RT_primer", codebook="SP1", no=1,
        )


# ---------------------------------------------------------------------------
# make_rt_primer_name — format `{name}_{pos}_RT_primer`
# ---------------------------------------------------------------------------


def test_make_rt_primer_name_basic():
    name = ps.make_rt_primer_name(name="BRCA1", position=12345)
    assert name == "BRCA1_12345_RT_primer"


def test_make_rt_primer_name_clone_id():
    name = ps.make_rt_primer_name(name="BZ23_clone23_TRB", position=45)
    assert name == "BZ23_clone23_TRB_45_RT_primer"


# ---------------------------------------------------------------------------
# resolve_codebook
# ---------------------------------------------------------------------------


def test_resolve_codebook_uses_explicit_value(tmp_path):
    backbone = tmp_path / "backbone_SP500.xlsx"
    backbone.touch()
    assert ps.resolve_codebook("SP369", backbone) == "SP369"


def test_resolve_codebook_extracts_from_filename(tmp_path):
    backbone = tmp_path / "backbone_SP369.xlsx"
    backbone.touch()
    assert ps.resolve_codebook(None, backbone) == "SP369"


def test_resolve_codebook_supports_non_sp_codebooks(tmp_path):
    # Codebooks aren't always SP-prefixed; confirm regex matches MIT30.
    backbone = tmp_path / "MIT30_panel.xlsx"
    backbone.touch()
    assert ps.resolve_codebook(None, backbone) == "MIT30"


def test_resolve_codebook_errors_when_unresolvable(tmp_path):
    backbone = tmp_path / "panel.xlsx"
    backbone.touch()
    with pytest.raises(ValueError, match="codebook"):
        ps.resolve_codebook(None, backbone)


# ---------------------------------------------------------------------------
# apply_final_column_order
# ---------------------------------------------------------------------------


def test_apply_final_column_order_padlock_reorders_known_columns():
    df = pd.DataFrame({
        "probe_arm5": ["aaaa"], "probe_name": ["X"], "order": [None],
        "probe_sequence": ["aaa"], "gene_name": ["X"], "chemistry": ["dRNA"],
        "No.": [1], "codebook": ["SP1"], "tm": [60.0],
    })
    out = ps.apply_final_column_order(df, kind="padlock")
    assert out.columns.tolist()[:3] == ["order", "probe_name", "probe_sequence"]
    # All input columns survive (positions are reordered).
    assert set(out.columns) >= set(df.columns)


def test_apply_final_column_order_padlock_pads_missing_with_nan():
    df = pd.DataFrame({
        "probe_name": ["X"], "probe_sequence": ["aaa"],
        "gene_name": ["X"], "chemistry": ["dRNA"],
        "No.": [1], "codebook": ["SP1"],
    })
    out = ps.apply_final_column_order(df, kind="padlock")
    # `order` is part of the schema even if not provided by the caller.
    assert "order" in out.columns
    assert pd.isna(out["order"].iloc[0])


def test_apply_final_column_order_drops_legacy_columns_silently_when_renamed():
    # Caller is expected to rename before applying order; any leftover legacy
    # columns get appended at the tail rather than silently dropped, so
    # nothing is lost.
    df = pd.DataFrame({
        "probe_name": ["X"], "probe_sequence": ["aaa"],
        "gene_name": ["X"], "chemistry": ["dRNA"],
        "No.": [1], "codebook": ["SP1"],
        "some_unknown_legacy_col": ["junk"],
    })
    out = ps.apply_final_column_order(df, kind="padlock")
    assert "some_unknown_legacy_col" in out.columns
    assert out.columns.tolist().index("some_unknown_legacy_col") >= len(
        ps.FINAL_PADLOCK_COLUMNS
    )


def test_apply_final_column_order_rejects_unknown_kind():
    with pytest.raises(ValueError):
        ps.apply_final_column_order(pd.DataFrame(), kind="binding_sites")


# ---------------------------------------------------------------------------
# legacy renames — exhaustive map covering everything the migration CLI needs
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("old,new", [
    ("Probe_Name", "probe_name"),
    ("Probe_Seq", "probe_sequence"),
    ("probe_seq", "probe_sequence"),
    ("RT_Primer_Name", "probe_name"),
    ("RT_Primer_Sequence", "probe_sequence"),
    ("Final_Probe_Sequence", "probe_sequence"),
    ("Probe_Length", "probe_length"),
    ("RT_Primer_Length", "probe_length"),
    ("arm_5prime", "probe_arm5"),
    ("arm_3prime", "probe_arm3"),
    # Schema-v2 keeps `g_content` and `gc_content` as DISTINCT columns
    # (G-only vs G+C). `g_content → gc_content` would silently change the
    # column's meaning and must NOT be in LEGACY_RENAMES — see Phase 0 notes.
    ("GC_Percent", "gc_content"),
    ("Tm", "tm"),
    ("blast_hits_count", "blast_hit_count"),
    ("Chr", "chr"),
    ("Start", "genomic_start"),
    ("End", "genomic_end"),
    ("Ref", "ref"),
    ("Alt", "alt"),
    ("Strand", "strand"),
    ("Mutation_Type", "mutation_type"),
    ("tm_cDNA_warning", "tm_cdna_warning"),
    ("Padlock_Probe_Name", "paired_padlock_probe_name"),
    ("Padlock_Probe_No.", "paired_padlock_no"),
    ("mfe", "free_energy"),
    ("probe_id", "No."),
    ("probe_seq", "probe_sequence"),
])
def test_legacy_renames_map_old_to_new(old, new):
    assert ps.LEGACY_RENAMES[old] == new


def test_legacy_renames_drop_redundant_backbone_no():
    # backbone_No. duplicates No. and is dropped during migration.
    assert ps.LEGACY_RENAMES.get("backbone_No.") in (None, "drop")


# ---------------------------------------------------------------------------
# parsing — round-trip helper for the migration CLI and probe-order skill
# ---------------------------------------------------------------------------


def test_parse_padlock_name_round_trip():
    name = ps.make_padlock_name(
        name="Actb", position=213, chemistry="dRNA", codebook="SP369", no=300,
    )
    parts = ps.parse_padlock_name(name)
    assert parts["name"] == "Actb"
    assert parts["position"] == 213
    assert parts["chemistry"] == "dRNA"
    assert parts["codebook"] == "SP369"
    assert parts["no"] == 300


def test_parse_padlock_name_underscored_clone_id():
    name = "BZ23_clone23_TRB_45_cDNA_SP500_17"
    parts = ps.parse_padlock_name(name)
    assert parts["name"] == "BZ23_clone23_TRB"
    assert parts["position"] == 45
    assert parts["chemistry"] == "cDNA"
    assert parts["codebook"] == "SP500"
    assert parts["no"] == 17


def test_parse_padlock_name_rejects_rt_primer():
    assert ps.parse_padlock_name("BRCA1_12345_RT_primer") is None


def test_parse_padlock_name_rejects_legacy_format():
    assert ps.parse_padlock_name("BRCA1_iLock_12345_SP_369") is None


def test_parse_rt_primer_name_round_trip():
    name = ps.make_rt_primer_name(name="BRCA1", position=12345)
    parts = ps.parse_rt_primer_name(name)
    assert parts["name"] == "BRCA1"
    assert parts["position"] == 12345


def test_parse_rt_primer_name_rejects_padlock():
    assert ps.parse_rt_primer_name("Actb_213_dRNA_SP369_1") is None
