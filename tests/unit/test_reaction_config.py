"""Tests for the ReactionConditions buffer fields on FilterConfig.

The buffer is exposed as flat fields (so the generic YAML loader round-trips
them unchanged) plus a ``reaction_conditions()`` builder that hands a
ReactionConditions to the Tm/ΔG call sites.
"""
from __future__ import annotations

import os
import tempfile

import pytest

from probe_designer.chemistry import ReactionConditions
from probe_designer.config import ConfigManager, FilterConfig


class TestFilterConfigBufferFields:
    def test_defaults_match_protocol(self):
        fc = FilterConfig()
        assert fc.monovalent_mM == 75.0
        assert fc.mg_mM == 10.0
        assert fc.dntp_mM == 0.0
        assert fc.strand_nM == 100.0
        assert fc.formamide_pct == 20.0
        assert fc.formamide_deg_per_pct == 0.5
        assert fc.lab_temp_c == 45.0
        assert fc.tm_margin_c == 5.0
        assert fc.saltcorr == 5

    def test_reaction_conditions_builds_matching_object(self):
        fc = FilterConfig(mg_mM=4.0, formamide_pct=10.0, lab_temp_c=37.0)
        rc = fc.reaction_conditions()
        assert isinstance(rc, ReactionConditions)
        assert rc.mg_mM == 4.0
        assert rc.formamide_pct == 10.0
        assert rc.lab_temp_c == 37.0
        # untouched fields fall back to the protocol defaults
        assert rc.monovalent_mM == 75.0

    def test_reaction_conditions_validation_propagates(self):
        # A nonsensical buffer must fail fast when materialized.
        with pytest.raises(ValueError):
            FilterConfig(mg_mM=-1.0).reaction_conditions()

    def test_solvents_pass_through(self):
        rc = FilterConfig(solvents={"dmso": 5.0}).reaction_conditions()
        assert rc.solvents == {"dmso": 5.0}


class TestBufferYamlRoundTrip:
    def test_load_config_picks_up_buffer_keys(self):
        yaml_text = (
            "filter:\n"
            "  mg_mM: 4.0\n"
            "  formamide_pct: 0.0\n"
            "  saltcorr: 7\n"
        )
        with tempfile.NamedTemporaryFile(
            "w", suffix=".yaml", delete=False, encoding="utf-8"
        ) as fh:
            fh.write(yaml_text)
            path = fh.name
        try:
            cm = ConfigManager(path)
            assert cm.filter.mg_mM == 4.0
            assert cm.filter.formamide_pct == 0.0
            assert cm.filter.saltcorr == 7
            assert cm.filter.reaction_conditions().mg_mM == 4.0
        finally:
            os.unlink(path)

    def test_save_config_emits_buffer_keys(self):
        import yaml

        cm = ConfigManager()
        cm.filter.mg_mM = 4.0
        with tempfile.NamedTemporaryFile(
            "w", suffix=".yaml", delete=False, encoding="utf-8"
        ) as fh:
            path = fh.name
        try:
            cm.save_config(path)
            with open(path, encoding="utf-8") as fh:
                data = yaml.safe_load(fh)
            assert data["filter"]["mg_mM"] == 4.0
            assert "formamide_pct" in data["filter"]
            assert "saltcorr" in data["filter"]
        finally:
            os.unlink(path)
