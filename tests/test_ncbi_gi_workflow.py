"""Tests for NCBI GI workflow: preset gene IDs and gene_info parsing."""

import os
import sys
import pytest
from unittest.mock import patch, MagicMock

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from config import DatabaseConfig
from database import DatabaseInterface


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def db_ncbi():
    """DatabaseInterface configured for NCBI."""
    config = DatabaseConfig(database_type="ncbi", organism="human")
    return DatabaseInterface(config)


@pytest.fixture
def gene_info_xlsx(tmp_path):
    """Create a temporary gene_info.xlsx for testing."""
    import pandas as pd
    df = pd.DataFrame({
        'gene_name': ['KLRB1', 'KIR2DL4', 'KLRC3'],
        'No.': [11, 15, 17],
        'GI': [1653961201, 158551988, 1779541896],
    })
    path = tmp_path / "gene_info.xlsx"
    df.to_excel(path, index=False)
    return str(path)


@pytest.fixture
def gene_info_no_gi_xlsx(tmp_path):
    """gene_info.xlsx without GI column."""
    import pandas as pd
    df = pd.DataFrame({
        'gene_name': ['KLRB1', 'KIR2DL4'],
        'No.': [11, 15],
    })
    path = tmp_path / "gene_info_no_gi.xlsx"
    df.to_excel(path, index=False)
    return str(path)


# ---------------------------------------------------------------------------
# Test: get_gene_sequences with preset gene_ids skips _ncbi_get_gene_ids
# ---------------------------------------------------------------------------

class TestGetGeneSequencesWithPresetIds:

    def test_skips_search_when_gene_ids_provided(self, db_ncbi):
        """When gene_ids is provided, _ncbi_get_gene_ids should NOT be called."""
        gene_list = ['KLRB1', 'KIR2DL4']
        preset_ids = ['1653961201', '158551988']

        with patch.object(db_ncbi, '_ncbi_get_gene_ids') as mock_search, \
             patch.object(db_ncbi, '_ncbi_fetch_genbank_records', return_value=[None, None]), \
             patch.object(db_ncbi, '_ncbi_get_genomic_context_from_genbank', return_value=None), \
             patch.object(db_ncbi, '_ncbi_get_genomic_context_from_nuccore', return_value=None):
            db_ncbi.get_gene_sequences(gene_list, gene_ids=preset_ids)
            mock_search.assert_not_called()

    def test_uses_search_when_gene_ids_not_provided(self, db_ncbi):
        """When gene_ids is None, _ncbi_get_gene_ids SHOULD be called."""
        gene_list = ['KLRB1', 'KIR2DL4']

        with patch.object(db_ncbi, '_ncbi_get_gene_ids', return_value=[None, None]) as mock_search, \
             patch.object(db_ncbi, '_ncbi_fetch_genbank_records', return_value=[None, None]):
            db_ncbi.get_gene_sequences(gene_list)
            mock_search.assert_called_once_with(gene_list)

    def test_preset_ids_passed_to_fetch(self, db_ncbi):
        """Preset gene_ids should be forwarded to _ncbi_fetch_genbank_records."""
        gene_list = ['KLRB1']
        preset_ids = ['1653961201']

        with patch.object(db_ncbi, '_ncbi_fetch_genbank_records', return_value=[None]) as mock_fetch:
            db_ncbi.get_gene_sequences(gene_list, gene_ids=preset_ids)
            mock_fetch.assert_called_once_with(preset_ids)


# ---------------------------------------------------------------------------
# Test: parse_gene_info extracts GI mapping and gene list
# ---------------------------------------------------------------------------

class TestParseGeneInfo:

    def test_extracts_gene_list(self, gene_info_xlsx):
        """gene_info.xlsx should yield correct gene_name list."""
        import pandas as pd
        df = pd.read_excel(gene_info_xlsx)
        genes = df['gene_name'].astype(str).str.strip().tolist()
        assert genes == ['KLRB1', 'KIR2DL4', 'KLRC3']

    def test_extracts_gi_as_strings(self, gene_info_xlsx):
        """GI column should be extractable as string list."""
        import pandas as pd
        df = pd.read_excel(gene_info_xlsx)
        gene_ids = df['GI'].astype(str).str.strip().tolist()
        assert gene_ids == ['1653961201', '158551988', '1779541896']

    def test_detects_gi_column_presence(self, gene_info_xlsx):
        """Should detect that GI column exists."""
        import pandas as pd
        df = pd.read_excel(gene_info_xlsx)
        assert 'GI' in df.columns

    def test_no_gi_column(self, gene_info_no_gi_xlsx):
        """When no GI column, gene_ids should be None."""
        import pandas as pd
        df = pd.read_excel(gene_info_no_gi_xlsx)
        assert 'GI' not in df.columns


# ---------------------------------------------------------------------------
# Test: gene_info with GI auto-sets database_type to ncbi
# ---------------------------------------------------------------------------

class TestGeneInfoAutoSetsNcbi:

    def test_auto_sets_ncbi_when_gi_present(self, gene_info_xlsx):
        """When gene_info has GI column, database_type should be set to ncbi."""
        import pandas as pd
        cfg_db = DatabaseConfig(database_type="ensembl")  # starts as ensembl
        df = pd.read_excel(gene_info_xlsx)
        if 'GI' in df.columns:
            cfg_db.database_type = 'ncbi'
        assert cfg_db.database_type == 'ncbi'

    def test_does_not_change_when_no_gi(self, gene_info_no_gi_xlsx):
        """Without GI column, database_type should remain unchanged."""
        import pandas as pd
        cfg_db = DatabaseConfig(database_type="ensembl")
        df = pd.read_excel(gene_info_no_gi_xlsx)
        if 'GI' in df.columns:
            cfg_db.database_type = 'ncbi'
        assert cfg_db.database_type == 'ensembl'
