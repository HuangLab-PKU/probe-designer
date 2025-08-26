# Legacy Modules

This directory contains the original modules from the notebook-based implementation of the DNA probe design software. These modules are kept for reference and backward compatibility.

## Modules Overview

### Database Modules
- **`database_ensembl.py`** - Original Ensembl database interface for fetching gene sequences and transcript information
- **`database_ncbi.py`** - Original NCBI database interface using Entrez API for sequence retrieval
- **`database_local.py`** - Local database interface (minimal implementation)

### Search and Binding Modules
- **`search_binding.py`** - Core binding site search algorithms including position search and optimization
- **`seq_filter.py`** - Sequence filtering utilities for pre-BLAST filtering

### BLAST and Analysis Modules
- **`blast.py`** - BLAST analysis utilities for online and local BLAST execution
- **`mutation.py`** - Mutation analysis and detection utilities

### Isoform Processing Modules
- **`isoform_process.py`** - Isoform analysis and processing utilities
- **`isoform_visualization.py`** - Visualization tools for isoform analysis

### Utility Modules
- **`file_operation.py`** - File operation utilities
- **`barcode.py`** - Barcode generation and management utilities

## Migration Status

These modules have been replaced by the new modular architecture:

| Legacy Module | New Module | Status |
|---------------|------------|--------|
| `database_ensembl.py` | `database.py` | ✅ Replaced |
| `database_ncbi.py` | `database.py` | ✅ Integrated |
| `search_binding.py` | `search_strategies.py` | ✅ Replaced |
| `seq_filter.py` | `filtering.py` | ✅ Replaced |
| `blast.py` | `filtering.py` | ✅ Integrated |
| `barcode.py` | `probe_assembly.py` | ✅ Replaced |
| `file_operation.py` | `utils.py` | ✅ Replaced |

## Usage

These modules are primarily for reference. If you need to use the original functionality, you can import them directly:

```python
from src.legacy import search_binding, database_ensembl
```

However, it's recommended to use the new modular architecture in `src/` for better maintainability and features.

## Notes

- These modules contain the original logic from Jupyter notebooks
- They may have hardcoded parameters and limited error handling
- The new architecture provides better configuration management and error handling
- Some functionality has been enhanced in the new modules
