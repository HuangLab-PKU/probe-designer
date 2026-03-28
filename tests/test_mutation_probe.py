#!/usr/bin/env python3
"""Tests for the mutation probe design module (probe_designer.ext.mutation_probe)."""

import pytest

from probe_designer.ext.mutation_probe import MutationProbeDesigner, design_mutation_probes
from probe_designer.config import SearchConfig, FilterConfig, DatabaseConfig
from probe_designer.database import DatabaseInterface


# ---------------------------------------------------------------------------
# Test: imports
# ---------------------------------------------------------------------------

class TestImports:
    def test_core_imports(self):
        """Core config and database classes are importable."""
        assert SearchConfig is not None
        assert FilterConfig is not None
        assert DatabaseConfig is not None
        assert DatabaseInterface is not None

    def test_mutation_module_imports(self):
        """MutationProbeDesigner and helper function are importable."""
        assert MutationProbeDesigner is not None
        assert design_mutation_probes is not None


# ---------------------------------------------------------------------------
# Test: MutationProbeDesigner basics
# ---------------------------------------------------------------------------

class TestMutationProbeDesigner:

    @pytest.fixture
    def designer(self):
        return MutationProbeDesigner(SearchConfig(), FilterConfig())

    def test_mutation_type_snp(self, designer):
        assert designer._get_mutation_type("A", "T", 100, 101) == "SNP"

    def test_mutation_type_del(self, designer):
        assert designer._get_mutation_type("AT", "-", 100, 102) == "DEL"

    def test_mutation_type_ins(self, designer):
        assert designer._get_mutation_type("-", "AT", 100, 100) == "INS"


# ---------------------------------------------------------------------------
# Test: genomic sequence mutation
# ---------------------------------------------------------------------------

class TestGenomicSequenceMutation:

    @pytest.fixture
    def designer(self):
        return MutationProbeDesigner(SearchConfig(), FilterConfig())

    def test_snp_mutation(self, designer):
        seq = "ATCGATCGATCG"
        # SNP: replace position 2..2 (C -> G)
        mutated, new_start, new_end = designer._mutate_genomic_sequence(
            seq, 2, 2, "C", "G", strand=1
        )
        assert mutated == "ATGGATCGATCG" or mutated[2] == "G"

    def test_deletion(self, designer):
        seq = "ATCGATCGATCG"
        # DEL: remove positions 2..3 (CG -> -)
        mutated, new_start, new_end = designer._mutate_genomic_sequence(
            seq, 2, 3, "CG", "-", strand=1
        )
        assert len(mutated) == len(seq) - 2

    def test_insertion(self, designer):
        seq = "ATCGATCGATCG"
        # INS: insert GG after position 2
        mutated, new_start, new_end = designer._mutate_genomic_sequence(
            seq, 2, 2, "-", "GG", strand=1
        )
        assert len(mutated) == len(seq) + 2


# ---------------------------------------------------------------------------
# Test: padlock arm extraction
# ---------------------------------------------------------------------------

class TestPadlockArmExtraction:

    @pytest.fixture
    def designer(self):
        return MutationProbeDesigner(SearchConfig(), FilterConfig())

    def test_arm_keys_present(self, designer):
        target_seq = "ATCGATCGATCGATCGATCG"
        arm_info = designer._extract_padlock_arms(target_seq, 10, 5, 5, 1)
        for key in ("probe_arm5", "probe_arm3", "target_sequence"):
            assert key in arm_info


# ---------------------------------------------------------------------------
# Test: DatabaseInterface creation
# ---------------------------------------------------------------------------

class TestDatabaseInterface:

    def test_create_ensembl_interface(self):
        config = DatabaseConfig(
            database_type="ensembl",
            organism="human",
            coord_system_version="GRCh38",
            max_retries=3,
        )
        db = DatabaseInterface(config)
        assert db is not None
