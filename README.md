# Probe Designer

A comprehensive tool for designing padlock probes for **spatial transcriptomics** and **multiplex RNA FISH**. The program searches for optimal binding sites on mRNA sequences using multiple strategies and databases (Ensembl, NCBI).

## Features

- **Multiple Search Strategies**: Consensus, isoform-specific, and brute force approaches
- **Multi-Database Support**: Ensembl and NCBI databases with automatic fallback
- **Advanced Filtering**: Thermodynamic and BLAST-based specificity filtering
- **YAML Configuration**: Human-readable configuration files with detailed documentation
- **Pipeline Automation**: Multi-database runs with automatic result merging
- **Real-time Progress**: Live progress tracking and detailed logging

## Quick Start

### 1. Environment Setup

```bash
# Clone the repository
git clone --depth 1 https://github.com/HuangLab-PKU/probe-designer probe_designer
cd probe_designer

# Create conda environment
conda env create --file environment.yml
conda activate probe_design
```

### 2. Prepare Gene List

Create a text file with one gene name per line:
```bash
# data/example_genes.txt
Actb
Tuba1a
```

### 3. Run Pipeline

#### Single Strategy Run
```bash
# Consensus strategy (recommended for most cases)
python scripts/run_pipeline.py \
    --genes-file data/example_genes.txt \
    --project-id my_project \
    --strategy consensus \
    --databases ensembl

# Isoform-specific strategy
python scripts/run_pipeline.py \
    --genes-file data/example_genes.txt \
    --project-id my_project \
    --strategy specific \
    --databases ensembl

# Brute force strategy
python scripts/run_pipeline.py \
    --genes-file data/example_genes.txt \
    --project-id my_project \
    --strategy bruteforce \
    --databases ensembl
```

#### Multi-Database Run
```bash
# Run with both Ensembl and NCBI, automatically merge results
python scripts/run_pipeline.py \
    --genes-file data/example_genes.txt \
    --project-id my_project \
    --strategy consensus \
    --databases ensembl ncbi
```

### 4. Check Results

Results are organized in `results/{project_id}/`:
- `{database}_results/`: Individual database results
- `merged_results.xlsx`: Combined results from all databases
- `missing_genes.txt`: Genes without successful probe design

## Configuration

The system uses YAML configuration files with detailed parameter documentation:

### Main Configuration Files
- `configs/config_consensus.yaml`: Consensus strategy (recommended)
- `configs/config_specific.yaml`: Isoform-specific strategy
- `configs/config_bruteforce.yaml`: Brute force strategy
- `configs/config_template.yaml`: Template with all options documented

### Key Configuration Sections

#### Database Settings
```yaml
database:
  database_type: "ensembl"  # ensembl, ncbi
  max_retries: 3            # API retry attempts
  organism: "mouse"         # mouse, human, rat, zebrafish
```

#### Search Parameters
```yaml
search:
  search_strategy: "isoform_consensus"  # consensus, specific, bruteforce
  binding_site_length: 40              # Probe length
  max_binding_sites: 30                # Max sites per gene
  window_size: 50                      # Search window
  step_size: 1                         # Search step
```

#### Filtering Criteria
```yaml
filter:
  min_g_content: 0.3                   # Min G content (0-1)
  max_g_content: 0.7                   # Max G content (0-1)
  max_consecutive_g: 3                 # Max consecutive Gs
  min_tm: 60.0                         # Min melting temperature (°C)
  max_tm: 80.0                         # Max melting temperature (°C)
  max_tm_diff: 5.0                     # Max Tm difference between arms
  min_free_energy: -5.0                # Min RNA structure free energy
  check_rna_structure: true            # Enable RNA structure check
  require_specificity: true            # Enable BLAST specificity check
  final_probes_per_gene: 3             # Final probes to select per gene
```

#### BLAST Settings
```yaml
blast:
  blast_type: "local"                  # local, online
  database: "refseq_rna"               # BLAST database
  hitlist_size: 100                    # Number of hits to return
  evalue: 1.0e-5                       # E-value threshold
  species: ["Mus musculus"]            # Target species for specificity
```

## Search Strategies

### 1. Consensus Strategy (`isoform_consensus`)
- **Purpose**: Find binding sites common to multiple isoforms
- **Best for**: Genes with multiple splice variants
- **Output**: Probes that bind to maximum number of isoforms per gene

### 2. Isoform-Specific Strategy (`isoform_specific`)
- **Purpose**: Find binding sites unique to individual isoforms
- **Best for**: Distinguishing between closely related isoforms
- **Output**: Probes specific to individual transcript variants

### 3. Brute Force Strategy (`brute_force`)
- **Purpose**: Exhaustive search across entire transcript
- **Best for**: Genes with complex structure or when other strategies fail
- **Output**: All possible binding sites meeting criteria

## Pipeline Logic

The system implements a sophisticated multi-database pipeline:

1. **Primary Run**: Execute with primary database (usually Ensembl)
2. **Missing Gene Check**: Identify genes without successful probe design
3. **Secondary Run**: Run with alternative database (NCBI) for missing genes
4. **Result Merging**: Combine results from all databases
5. **Final Check**: Report any remaining missing genes

### Example Workflow
```bash
# Step 1: Run with Ensembl
python scripts/run_pipeline.py \
    --genes-file genes.txt \
    --project-id project_001 \
    --strategy consensus \
    --databases ensembl

# Step 2: Check missing genes
cat results/project_001/missing_genes.txt

# Step 3: Run with NCBI for missing genes
python scripts/run_pipeline.py \
    --genes-file results/project_001/missing_genes.txt \
    --project-id project_001_ncbi \
    --strategy consensus \
    --databases ncbi

# Step 4: Merge all results
python test/merge_results.py \
    --results-dir results/ \
    --output final_results.xlsx
```

## Filtering System

### Pre-BLAST Filters (Thermodynamic)
- **G Content**: 30-70% G nucleotides
- **Consecutive Gs**: Maximum 3 consecutive G nucleotides
- **Melting Temperature**: 60-80°C for both arms
- **Tm Difference**: Maximum 5°C difference between 3' and 5' arms
- **RNA Structure**: Minimum free energy > -5 kcal/mol
- **Secondary Structure**: RNA folding analysis (optional)

### Post-BLAST Filters (Specificity)
- **Organism Specificity**: BLAST hits only to target organisms
- **Gene Specificity**: No alignments to non-target genes
- **Complementarity**: Probe must be complementary to target
- **Alignment Count**: Controlled by hitlist_size parameter

## Output Files

### Main Results
- `probes_wanted.xlsx`: Final selected probes with all parameters
- `probes_candidates.xlsx`: All candidate probes before final selection
- `binding_sites_stats.json`: Statistical summary of binding sites
- `blast_results.xml`: BLAST analysis results
- `config_used.yaml`: Configuration file used for this run

### Sequence Files
- `total_binding_sites.fasta`: All candidate binding sites
- `blast_input.fasta`: Sequences submitted to BLAST
- `filtered_binding_sites.fasta`: Sequences passing all filters

### Logging and Debug
- `binding_sites.json`: Detailed binding site information
- `filtering_stats.json`: Filtering statistics
- `missing_genes.txt`: Genes without successful probe design

## Example Results

### Consensus Strategy Example
```bash
python scripts/run_pipeline.py \
    --genes-file data/example_genes.txt \
    --project-id demo_consensus \
    --strategy consensus \
    --databases ensembl
```

**Output Summary**:
- Original sites: 60
- After pre-BLAST filter: 60
- After specificity filter: 11
- Final probes: 6 (3 per gene)

### Specific Strategy Example
```bash
python scripts/run_pipeline.py \
    --genes-file data/example_genes.txt \
    --project-id demo_specific \
    --strategy specific \
    --databases ensembl
```

**Output Summary**:
- Original sites: 150
- After pre-BLAST filter: 150
- After specificity filter: 0 (strict specificity)

### Brute Force Strategy Example
```bash
python scripts/run_pipeline.py \
    --genes-file data/example_genes.txt \
    --project-id demo_bruteforce \
    --strategy bruteforce \
    --databases ensembl
```

**Output Summary**:
- Original sites: 200
- After pre-BLAST filter: 200
- After specificity filter: 52
- Final probes: 6 (3 per gene)

## Advanced Usage

### Custom Configuration
```bash
# Use custom configuration file
python scripts/run_pipeline.py \
    --genes-file genes.txt \
    --project-id custom_run \
    --strategy consensus \
    --config my_config.yaml
```

### BLAST Species Specification
```bash
# Specify target species for BLAST
python scripts/run_pipeline.py \
    --genes-file genes.txt \
    --project-id species_specific \
    --strategy consensus \
    --blast-species "Mus musculus" "Homo sapiens"
```

### Skip Result Merging
```bash
# Run without automatic merging
python scripts/run_pipeline.py \
    --genes-file genes.txt \
    --project-id no_merge \
    --strategy consensus \
    --databases ensembl ncbi \
    --skip-merge
```

## Troubleshooting

### Common Issues

1. **BLAST Command Not Found**
   - Install BLAST+ locally or use online BLAST
   - Set `blast_type: "online"` in configuration

2. **Ensembl API Timeouts**
   - Increase `max_retries` in database configuration
   - Check network connectivity

3. **No Sequences Retrieved**
   - Verify gene names match database conventions
   - Check organism setting matches gene list

4. **Low Specificity Results**
   - Adjust `hitlist_size` and `evalue` parameters
   - Review target species list

### Performance Optimization

1. **Local BLAST Setup**
   ```bash
   # Download BLAST+ and databases
   # Set blast_type: "local" in configuration
   ```

2. **Parallel Processing**
   - Use multiple databases simultaneously
   - Adjust `concurrency` parameter for online BLAST

3. **Memory Management**
   - Reduce `max_binding_sites` for large gene lists
   - Use smaller `window_size` for memory-constrained systems

## Dependencies

- **Python**: 3.8+
- **Biopython**: Sequence analysis and BLAST
- **Pandas**: Data manipulation
- **PyYAML**: Configuration file parsing
- **Requests**: API communication
- **tqdm**: Progress tracking

## Testing

Run the test suite to verify installation:
```bash
python test/test_basic_pipeline.py
```

This will test:
- YAML configuration loading
- Basic pipeline functionality
- Result merging
- Error handling

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this tool in your research, please cite:
```
Probe Designer: A comprehensive tool for padlock probe design in spatial transcriptomics
[Your Name et al.]
[Journal/Conference]
[Year]
```

