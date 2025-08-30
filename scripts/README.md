# DNA Probe Design Scripts

This directory contains the main scripts for DNA probe design with support for multiple databases and automated pipeline workflows.

## Available Scripts

### `run_pipeline.py` (NEW - Recommended)
Multi-database probe design pipeline that runs probe design with multiple databases and merges results.

**Usage:**
```bash
# Run consensus strategy with Ensembl only
python run_pipeline.py --genes-file genes.txt --project-id project_001 --strategy consensus

# Run specific strategy with both databases
python run_pipeline.py --genes-file genes.txt --project-id project_002 --strategy specific --databases ensembl ncbi

# Run bruteforce strategy with custom config
python run_pipeline.py --genes-file genes.txt --project-id project_003 --strategy bruteforce --databases ensembl ncbi --config my_config.yaml
```

**Parameters:**
- `--genes-file`: Path to gene list file
- `--project-id`: Project identifier for organizing results
- `--strategy`: Probe design strategy (consensus, specific, bruteforce)
- `--databases`: Databases to use (ensembl, ncbi)
- `--config`: Custom configuration file (YAML or JSON)
- `--blast-species`: Target species for BLAST search
- `--organism`: Target organism (mouse, human, rat, zebrafish)
- `--skip-merge`: Skip result merging

### `find_consensus_probes.py`
Finds consensus binding sites (probes that bind to maximum number of isoforms per gene) using Ensembl and IsoformConsensusStrategy.

**Usage:**
```bash
python find_consensus_probes.py --genes_file genes.txt --config configs/config_consensus.yaml
```

### `find_specific_probes.py`
Finds isoform-specific binding sites per gene using Ensembl and IsoformSpecificStrategy.

**Usage:**
```bash
python find_specific_probes.py --genes_file genes.txt --config configs/config_specific.yaml
```

### `find_probes.py`
Finds binding sites using brute force strategy (works with both Ensembl and NCBI databases).

**Usage:**
```bash
python find_probes.py --genes_file genes.txt --config configs/config_bruteforce.yaml --database ensembl
```

### `probe_stitcher.py`
Assembles probes into panels with barcodes and primers.

**Usage:**
```bash
python probe_stitcher.py --input filtered_binding_sites.xlsx --output probe_panel.xlsx
```

## Pipeline Logic

The new pipeline supports the following workflow:

1. **Single Database Run**: Run probe design with one database (Ensembl or NCBI)
2. **Multi-Database Run**: Run probe design with multiple databases sequentially
3. **Result Merging**: Automatically merge results from multiple databases
4. **Missing Gene Check**: Identify genes that don't have probes designed
5. **Project Organization**: Organize results by project ID with separate directories

### Pipeline Features

- **Automatic Database Selection**: Choose which databases to use for each run
- **Strategy Support**: Support for consensus, specific, and bruteforce strategies
- **Result Merging**: Automatic merging of results from multiple databases
- **Missing Gene Detection**: Identify genes that need additional probe design
- **Project Management**: Organize results by project with clear directory structure
- **Configuration Flexibility**: Support for both YAML and JSON configuration files

### Example Workflow

```bash
# Step 1: Run with Ensembl database
python run_pipeline.py --genes-file genes.txt --project-id my_project --strategy consensus --databases ensembl

# Step 2: Check for missing genes and run with NCBI if needed
python run_pipeline.py --genes-file missing_genes.txt --project-id my_project_ncbi --strategy bruteforce --databases ncbi

# Step 3: Merge all results
python test/merge_results.py --results-dir results/ --gene-list genes.txt --output final_results.xlsx
```

## Configuration

### YAML Configuration (NEW)

The software now supports YAML configuration files with better readability and comments:

```yaml
# Configuration for finding consensus binding sites
database:
  database_type: "ensembl"
  api_key: "your_ncbi_api_key_here"
  email: "your_email@example.com"

search:
  search_strategy: "isoform_consensus"
  binding_site_length: 40
  max_binding_sites: 30

filter:
  min_g_content: 0.3
  max_g_content: 0.7
  min_tm: 60.0
  max_tm: 80.0
  target_organisms: ["Mus musculus", "Homo sapiens"]

blast:
  blast_type: "local"
  database: "nt"
  species: ["Mus musculus"]

output:
  output_dir: "results"
  save_fasta: true
  save_json: true
  save_excel: true

species: "mouse"
```

### Configuration Files

- `configs/config_template.yaml` - Template configuration with detailed comments
- `configs/config_consensus.yaml` - Configuration for consensus probe design
- `configs/config_specific.yaml` - Configuration for specific probe design
- `configs/config_bruteforce.yaml` - Configuration for brute force search

### Parameter Descriptions

#### Database Configuration
- `database_type`: Database type (ensembl, ncbi)
- `api_key`: NCBI API key for database access
- `email`: Email for API requests

#### Search Configuration
- `search_strategy`: Search strategy (brute_force, exon_junction, isoform_specific, isoform_consensus)
- `binding_site_length`: Length of binding sites in base pairs
- `max_binding_sites`: Maximum binding sites per gene

#### Filter Configuration
- `min_g_content` / `max_g_content`: G content limits (0.0-1.0)
- `max_consecutive_g`: Maximum consecutive G nucleotides
- `min_tm` / `max_tm`: Melting temperature limits (°C)
- `max_tm_diff`: Maximum temperature difference between 3' and 5' arms
- `target_organisms`: Target organisms for specificity

#### BLAST Configuration
- `blast_type`: BLAST type (local, online)
- `database`: BLAST database (nt, refseq_rna, nr)
- `species`: Target species for BLAST search

#### Output Configuration
- `output_dir`: Output directory path
- `save_fasta`: Save FASTA format files
- `save_json`: Save JSON format files
- `save_excel`: Save Excel format files

## Search Strategies

### Consensus Strategy (`isoform_consensus`)
Finds probes that bind to the maximum number of isoforms per gene, maximizing coverage.

### Specific Strategy (`isoform_specific`)
Finds probes unique to individual isoforms for isoform-specific targeting.

### Brute Force Strategy (`brute_force`)
Scans the entire sequence with comprehensive filtering and optimization.

### Exon Junction Strategy (`exon_junction`)
Searches for binding sites at exon-exon boundaries for higher specificity.

## Output Structure

### Project Directory Structure
```
results/
├── project_001/
│   ├── ensembl_results/
│   │   ├── binding_sites.json
│   │   ├── filtered_binding_sites.xlsx
│   │   └── filtering_stats.json
│   ├── ncbi_results/
│   │   ├── binding_sites.json
│   │   ├── filtered_binding_sites.xlsx
│   │   └── filtering_stats.json
│   ├── merged_results.xlsx
│   └── missing_genes.txt
```

### Output Files
- `binding_sites.json`: Raw binding sites data
- `filtered_binding_sites.xlsx`: Filtered binding sites in Excel format
- `filtering_stats.json`: Filtering statistics
- `merged_results.xlsx`: Merged results from multiple databases
- `missing_genes.txt`: Genes that need additional probe design

## Dependencies

Required Python packages:
- pandas
- numpy
- biopython
- requests
- tqdm
- pyyaml

Install dependencies:
```bash
pip install pandas numpy biopython requests tqdm pyyaml
```

## Testing

Run the test suite to validate the pipeline:

```bash
python test/test_pipeline.py
```

## Troubleshooting

### Common Issues

1. **Configuration Loading**
   - Ensure YAML files are properly formatted
   - Check file paths and permissions

2. **Database Access**
   - Verify API keys and email addresses
   - Check internet connection for online databases

3. **BLAST Analysis**
   - For local BLAST: Ensure BLAST+ is installed
   - For online BLAST: Check NCBI service status

4. **Memory Issues**
   - Reduce `max_binding_sites` in configuration
   - Process genes in smaller batches

### Support

For issues and questions, please refer to the main project documentation or contact the development team.
