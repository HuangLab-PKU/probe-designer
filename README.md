# Probe Designer

A comprehensive tool for designing padlock probes for **spatial transcriptomics** and **multiplex RNA FISH**. The program searches for optimal binding sites on mRNA sequences using multiple strategies and databases (Ensembl, NCBI).

## 📚 Documentation

**For detailed documentation, please visit the [docs/](docs/) directory:**

- **[Installation Guide](docs/installation.md)** - Setup and installation instructions
- **[Development Guide](docs/DEVELOPMENT.md)** - Contributing and development guidelines  
- **[Changelog](docs/CHANGELOG.md)** - Detailed change history and technical updates
- **[Documentation Index](docs/README.md)** - Complete documentation navigation

## 🚀 Quick Start

### 1. Environment Setup
```bash
# Create conda environment
conda env create --file environment.yml
conda activate probe_design
```

### 2. Run Pipeline
```bash
# Consensus strategy (recommended)
python scripts/run_pipeline.py \
    --genes-file data/example_genes.txt \
    --project-id my_project \
    --strategy consensus \
    --databases ensembl
```

### 3. Check Results
Results are organized in `results/{project_id}/`:
- `{database}_results/`: Individual database results
- `merged_results.xlsx`: Combined results from all databases
- `missing_genes.txt`: Genes without successful probe design

## 🔧 Key Features

- **Multiple Search Strategies**: Consensus, isoform-specific, and brute force approaches
- **Multi-Database Support**: Ensembl and NCBI databases with automatic fallback
- **Advanced Filtering**: Thermodynamic and BLAST-based specificity filtering
- **YAML Configuration**: Human-readable configuration files with detailed documentation
- **Pipeline Automation**: Multi-database runs with automatic result merging
- **Real-time Progress**: Live progress tracking and detailed logging

## 📁 Project Structure

```
code/
├── docs/                    # 📚 Complete documentation
├── src/                     # 🔧 Core source code
├── scripts/                 # 🚀 Pipeline scripts
├── configs/                 # ⚙️ Configuration files
├── data/                    # 📊 Example data
├── notebooks/               # 📓 Jupyter notebooks
├── local/                   # 🔬 Local analysis scripts
└── test/                    # 🧪 Test suite
```

## ⚙️ Configuration

Main configuration files in `configs/`:
- `config_consensus.yaml`: Consensus strategy (recommended)
- `config_specific.yaml`: Isoform-specific strategy  
- `config_bruteforce.yaml`: Brute force strategy
- `config_template.yaml`: Template with all options documented

## 🔬 Search Strategies

1. **Consensus**: Find binding sites common to multiple isoforms
2. **Isoform-Specific**: Find binding sites unique to individual isoforms
3. **Brute Force**: Exhaustive search across entire transcript

## 📊 Output Files

Results are saved in `results/{project_id}/`:
- `probes_wanted.xlsx`: Final selected probes
- `probes_candidates.xlsx`: All candidate probes
- `binding_sites_stats.json`: Statistical summary
- `missing_genes.txt`: Genes without successful probe design

## 📄 License

This project is licensed under the MIT License.

---

**For detailed documentation, configuration examples, troubleshooting, and advanced usage, please visit the [docs/](docs/) directory.**
