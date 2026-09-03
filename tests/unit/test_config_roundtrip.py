"""A config must survive save -> load, and a YAML key must actually do something.

Audit R7, generalized. `save_config` mirrored the `FilterConfig` fields by hand,
so every field added since it was written silently vanished on save: the Phase 1
`enforce_tm_gate`, the whole RNAplfold block, `max_alignments`,
`target_organisms`. Load, meanwhile, applies `if hasattr(...)` and skips
everything else — so an unrecognised key (a typo, or a knob that never existed)
is discarded without a word. `final_probes_per_gene` sat in three shipped YAMLs
that way: ignored on load, echoed back as a hardcoded 3 on save, while the real
knob was the `--top-n` CLI option.

The two tests here close the loop from both ends — nothing is lost on the way
out, and nothing is silently ignored on the way in.
"""
from __future__ import annotations

import dataclasses

import pytest
import yaml

from probe_designer.config import ConfigManager, FilterConfig, SearchConfig

SHIPPED_CONFIGS = [
    "config_consensus.yaml",
    "config_specific.yaml",
    "config_single_sequence.yaml",
    "config_template.yaml",
]


def _configs_dir():
    from pathlib import Path
    return Path(__file__).resolve().parents[2] / "configs"


class TestRoundTrip:
    def test_every_filter_field_survives_save_then_load(self, tmp_path):
        cm = ConfigManager()
        # Set every filter field to a non-default so a dropped one is visible.
        cm.filter.enforce_tm_gate = True
        cm.filter.enforce_tm_diff_gate = True
        cm.filter.max_tm_diff = 3.5
        cm.filter.plfold_window = 200
        cm.filter.plfold_span = 120
        cm.filter.plfold_temperature = 42.0
        cm.filter.min_accessibility = 0.3
        cm.filter.mg_mM = 7.5
        cm.filter.solvents = {"dmso": 5.0}
        cm.filter.max_alignments = 9

        out = tmp_path / "roundtrip.yaml"
        cm.save_config(str(out))

        restored = ConfigManager()
        restored.load_config(str(out))
        assert dataclasses.asdict(restored.filter) == dataclasses.asdict(cm.filter)

    def test_saved_yaml_carries_every_filter_field(self, tmp_path):
        """Guards the specific failure: a field the hand-written dict forgot."""
        cm = ConfigManager()
        out = tmp_path / "full.yaml"
        cm.save_config(str(out))
        written = yaml.safe_load(out.read_text(encoding="utf-8"))["filter"]
        expected = {f.name for f in dataclasses.fields(FilterConfig)}
        assert expected - set(written) == set()

    def test_every_search_field_survives_save_then_load(self, tmp_path):
        # Same hand-written-dict trap as `filter`, one section over. The
        # isoform-crediting knobs added 2026-08-02 are safety gates: losing
        # them on a round trip silently reverts a run to the old behaviour.
        cm = ConfigManager()
        cm.search.min_isoform_arm_nt = 18
        cm.search.require_protein_coding = False
        cm.search.binding_site_length = 36
        cm.search.step_size = 2

        out = tmp_path / "roundtrip_search.yaml"
        cm.save_config(str(out))

        restored = ConfigManager()
        restored.load_config(str(out))
        assert dataclasses.asdict(restored.search) == dataclasses.asdict(cm.search)

    def test_saved_yaml_carries_every_search_field(self, tmp_path):
        cm = ConfigManager()
        out = tmp_path / "full_search.yaml"
        cm.save_config(str(out))
        written = yaml.safe_load(out.read_text(encoding="utf-8"))["search"]
        expected = {f.name for f in dataclasses.fields(SearchConfig)}
        assert expected - set(written) == set()


class TestBufferKnobsAreReachable:
    """Every ReactionConditions knob must be settable from a config file.

    ``FilterConfig`` deliberately mirrors ``ReactionConditions`` as flat fields
    because the generic YAML loader cannot round-trip a nested dataclass. That
    mirror is a duplication, so it drifts: ``hybrid_nn_model`` was added to
    ReactionConditions and was unreachable from YAML until this test existed.
    """

    def test_filter_config_exposes_every_reaction_field(self):
        from probe_designer.chemistry import ReactionConditions

        reaction_fields = {f.name for f in dataclasses.fields(ReactionConditions)}
        filter_fields = {f.name for f in dataclasses.fields(FilterConfig)}
        missing = reaction_fields - filter_fields
        assert not missing, (
            f"ReactionConditions knob(s) {sorted(missing)} cannot be set from a "
            "config file — add the flat field to FilterConfig and forward it in "
            "reaction_conditions()"
        )

    def test_values_actually_reach_the_reaction(self):
        """Present-but-not-forwarded is the same bug one step later."""
        cfg = FilterConfig(hybrid_nn_model="banerjee2020", mg_mM=4.0, lab_temp_c=42.0)
        reaction = cfg.reaction_conditions()
        assert reaction.hybrid_nn_model == "banerjee2020"
        assert reaction.mg_mM == 4.0
        assert reaction.lab_temp_c == 42.0


class TestNoDeadKeys:
    @pytest.mark.parametrize("name", SHIPPED_CONFIGS)
    def test_shipped_yaml_filter_keys_all_map_to_a_field(self, name):
        """A key that maps to no field is silently discarded on load."""
        data = yaml.safe_load((_configs_dir() / name).read_text(encoding="utf-8"))
        known = {f.name for f in dataclasses.fields(FilterConfig)}
        unknown = set(data.get("filter") or {}) - known
        assert not unknown, (
            f"{name} sets filter keys that FilterConfig does not define and the "
            f"loader therefore ignores: {sorted(unknown)}"
        )

    @pytest.mark.parametrize("name", SHIPPED_CONFIGS)
    def test_shipped_yaml_search_keys_all_map_to_a_field(self, name):
        data = yaml.safe_load((_configs_dir() / name).read_text(encoding="utf-8"))
        known = {f.name for f in dataclasses.fields(SearchConfig)}
        unknown = set(data.get("search") or {}) - known
        assert not unknown, (
            f"{name} sets search keys that SearchConfig does not define and the "
            f"loader therefore ignores: {sorted(unknown)}"
        )

    def test_loader_warns_instead_of_silently_dropping(self, tmp_path):
        cfg = tmp_path / "typo.yaml"
        cfg.write_text(
            yaml.safe_dump({"filter": {"formamide_pc": 20.0}}), encoding="utf-8",
        )
        cm = ConfigManager()
        with pytest.warns(UserWarning, match="formamide_pc"):
            cm.load_config(str(cfg))

    def test_loader_warns_on_an_unknown_search_key(self, tmp_path):
        cfg = tmp_path / "typo_search.yaml"
        cfg.write_text(
            yaml.safe_dump({"search": {"min_isoform_arm": 16}}), encoding="utf-8",
        )
        cm = ConfigManager()
        with pytest.warns(UserWarning, match="min_isoform_arm"):
            cm.load_config(str(cfg))
