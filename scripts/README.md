# Automated DNA Probe Design Software

This directory contains the main executable scripts for the automated DNA probe design software.

## Files

- `main.py` - Main entry point for probe design workflow
- `merge_results.py` - Utility script for merging results from multiple runs
- `test_probe_designer.py` - Test script for validating components

## Usage

### Main Script (`main.py`)

The main script orchestrates the complete probe design workflow:

```bash
# Basic usage with gene list file
python main.py gene_list.xlsx

# Use custom configuration
python main.py gene_list.xlsx --config my_config.json

# Create configuration template
python main.py --create-config template.json

# Check dependencies
python main.py --check-deps

# Show help
python main.py --help
```

#### Parameters

- `gene_list_file` - Path to gene list file (.xlsx, .csv, or .txt)
- `--config, -c` - Path to configuration JSON file
- `--create-config` - Create a configuration template file
- `--check-deps` - Check if all required dependencies are installed
- `--version` - Show version information

#### Input Format

The gene list file should contain a column named "gene_name" with gene identifiers:

**Excel (.xlsx):**
| gene_name |
|-----------|
| Gene1     |
| Gene2     |
| Gene3     |

**CSV (.csv):**
```csv
gene_name
Gene1
Gene2
Gene3
```

**Text (.txt):**
```
Gene1
Gene2
Gene3
```

#### Output

The script generates the following output files in a timestamped directory:

- `gene_name_list_tosearch.txt` - Processed gene list
- `blast_results.xml` - BLAST analysis results
- `filtered_binding_sites.xlsx` - Filtered binding sites
- `filtered_binding_sites.json` - Filtered binding sites (JSON format)
- `filtering_stats.json` - Filtering statistics
- `{panel_type}_probes.xlsx` - Final probe sequences
- `{panel_type}_probes.fasta` - Probe sequences in FASTA format
- `{panel_type}_probes.json` - Probe data in JSON format
- `{panel_type}_stats.json` - Probe statistics
- `statistics.json` - Overall workflow statistics
- `probe_design.log` - Execution log

### Merge Results Script (`merge_results.py`)

Utility script for merging results from multiple probe design runs:

```bash
python merge_results.py results_directory output_file.xlsx
```

#### Parameters

- `results_directory` - Directory containing subdirectories with results
- `output_file` - Output file path for merged results

#### Usage

This script searches for `probes_wanted.xlsx` files in subdirectories under the specified results directory and merges them into a single file, removing duplicates and sorting by gene name and position.

## Configuration

### Configuration File Format

The software uses JSON configuration files with the following structure:

```json
{
  "database": {
    "organism": "mouse",
    "database_type": "ensembl",
    "coord_system_version": "GRCh38",
    "max_retries": 3
  },
  "search": {
    "binding_site_length": 40,
    "max_binding_sites": 30,
    "search_strategy": "exon_junction",
    "step_size": null
  },
  "filter": {
    "min_g_content": 0.4,
    "max_g_content": 0.7,
    "max_consecutive_g": 4,
    "min_tm": 45.0,
    "max_tm": 65.0,
    "min_free_energy": -10.0,
    "max_alignments": 5,
    "require_specificity": true,
    "target_organisms": ["Mus musculus", "Homo sapiens"]
  },
  "probe": {
    "panel_type": "PRISM",
    "barcode_file": null,
    "primer_left": null,
    "primer_right": null
  },
  "blast": {
    "blast_type": "local",
    "database": "refseq_rna",
    "task": "megablast",
    "evalue": 1e-5
  },
  "output": {
    "output_dir": "results",
    "create_timestamp": true,
    "save_intermediate": true,
    "file_formats": ["xlsx", "fasta", "json"]
  }
}
```

### Parameter Descriptions

#### Database Configuration
- `organism` - Target organism (mouse, human)
- `database_type` - Database type (ensembl, ncbi)
- `coord_system_version` - Genome assembly version
- `max_retries` - Maximum retry attempts for API calls

#### Search Configuration
- `binding_site_length` - Length of binding sites to search for
- `max_binding_sites` - Maximum number of binding sites per gene
- `search_strategy` - Search strategy (exon_junction, brute_force, isoform_specific)
- `step_size` - Step size for brute force search (null for exon junction)

#### Filter Configuration
- `min_g_content` - Minimum G content threshold
- `max_g_content` - Maximum G content threshold
- `max_consecutive_g` - Maximum consecutive G nucleotides
- `min_tm` - Minimum melting temperature
- `max_tm` - Maximum melting temperature
- `min_free_energy` - Minimum RNA secondary structure free energy
- `max_alignments` - Maximum BLAST alignments allowed
- `require_specificity` - Whether to require target organism specificity
- `target_organisms` - List of target organisms for specificity check

#### Probe Configuration
- `panel_type` - Panel type (PRISM, SPRINTseq, custom)
- `barcode_file` - Path to barcode file (null for defaults)
- `primer_left` - Left primer sequence (null for defaults)
- `primer_right` - Right primer sequence (null for defaults)

#### BLAST Configuration
- `blast_type` - BLAST type (local, online)
- `database` - BLAST database name
- `task` - BLAST task type
- `evalue` - E-value threshold

#### Output Configuration
- `output_dir` - Base output directory
- `create_timestamp` - Whether to create timestamped subdirectories
- `save_intermediate` - Whether to save intermediate results
- `file_formats` - Output file formats

## Search Strategies

### Exon Junction Strategy
Searches for binding sites at exon-exon boundaries, which are often more specific and less prone to off-target binding.

### Brute Force Strategy
Scans the entire sequence with a sliding window to find all potential binding sites.

### Isoform Specific Strategy
Identifies unique regions across different isoforms to ensure probe specificity.

## Panel Types

### PRISM
Assembles probes with the format: right arm + BARCODE + left arm

### SPRINTseq
Assembles probes with the format: binding site + composite barcode (primer + barcode + primer + barcode)

### Custom
Uses raw binding sequences without additional components.

## Dependencies

Required Python packages:
- pandas
- numpy
- biopython
- requests
- tqdm

Optional packages:
- RNAfold (for RNA secondary structure analysis)

Install dependencies:
```bash
pip install pandas numpy biopython requests tqdm
```

## Troubleshooting

### Common Issues

1. **Missing dependencies**
   - Run `python main.py --check-deps` to check installed packages
   - Install missing packages with pip

2. **BLAST errors**
   - For local BLAST: Ensure BLAST+ is installed and database is available
   - For online BLAST: Check internet connection and NCBI service status

3. **API rate limits**
   - Increase `max_retries` in database configuration
   - Add delays between requests if needed

4. **Memory issues**
   - Reduce `max_binding_sites` in search configuration
   - Process genes in smaller batches

### Log Files

Check the `probe_design.log` file in the output directory for detailed execution information and error messages.

### Support

For issues and questions, please refer to the main project documentation or contact the development team.
