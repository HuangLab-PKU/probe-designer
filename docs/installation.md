# Installation Guide

This guide will help you set up the probe design environment on your system.

## Prerequisites

### System Requirements
- **Operating System**: Windows 10/11, macOS 10.14+, or Linux
- **Python**: 3.8 or higher
- **Memory**: At least 8GB RAM (16GB recommended for large datasets)
- **Storage**: At least 5GB free space

### Platform-Specific Notes
✅ **Good News**: ViennaRNA is now available on **Windows**, **Linux**, and **macOS** systems!

- **Windows**: ViennaRNA 2.7.0+ is available via pip
- **Linux/macOS**: ViennaRNA is available via conda or pip
- **All platforms**: RNA structure prediction is fully supported

### Required Software
- **Conda** or **Miniconda**: For environment management
- **Git**: For version control
- **Text Editor**: VS Code, PyCharm, or similar

## Installation Steps

### 1. Clone the Repository

```bash
git clone <repository-url> probe_designer
cd probe_designer
```

### 2. Create Environment

#### Option A: Using Conda (Recommended)
```bash
# Create environment from environment.yml
conda env create -f code/environment.yml

# Activate environment
conda activate probe_design
```

#### Option B: Using pip
```bash
# Create virtual environment
python -m venv probe_design_env

# Activate environment
# On Windows:
probe_design_env\Scripts\activate
# On Linux/macOS:
source probe_design_env/bin/activate

# Install dependencies
pip install -r code/requirements.txt
```

### 3. Verify Installation

```bash
# Test basic functionality
python -c "import src.mutation_probe; print('Installation successful!')"

# Test ViennaRNA installation
python -c "import RNA; print('ViennaRNA version:', RNA.__version__)"

# Run a simple test
python local/probe_design_BZ09_TNBC_mut.py --help
```

## Configuration

### 1. Genome Data Setup

You need to have the reference genome FASTA file available:

```bash
# Example path (adjust as needed)
genome_path = "G:/genome/GRCh38.fa"
```

### 2. Configuration Files

Copy and modify configuration files as needed:

```bash
# Copy template configuration
cp configs/config_template.yaml configs/my_config.yaml

# Edit configuration file
nano configs/my_config.yaml
```

### 3. Environment Variables (Optional)

```bash
# Set genome path
export GENOME_PATH="/path/to/genome.fa"

# Set output directory
export OUTPUT_DIR="/path/to/output"
```

## Troubleshooting

### Common Issues

#### 1. Conda Environment Issues
```bash
# If conda command not found
# Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html

# If environment activation fails
conda deactivate
conda activate probe_design
```

#### 2. Python Import Errors
```bash
# Check Python version
python --version

# Reinstall packages
pip install --force-reinstall -r requirements.txt
```

#### 3. Genome File Issues
```bash
# Check if genome file exists and is readable
ls -la /path/to/genome.fa

# Test genome access
python -c "from pyfaidx import Fasta; f = Fasta('/path/to/genome.fa'); print('Genome accessible')"
```

#### 4. RNA Structure Prediction Issues
```bash
# Check if ViennaRNA is available
python -c "import RNA; print('ViennaRNA version:', RNA.__version__)"

# Test RNA structure prediction
python -c "import RNA; structure, mfe = RNA.fold('GGGGAAAACCCC'); print(f'Structure: {structure}, MFE: {mfe}')"

# If ViennaRNA not available, install it
pip install ViennaRNA

# If still having issues, disable RNA structure checks in config
# Set check_rna_structure: false in your configuration file
```

#### 5. Memory Issues
```bash
# For large datasets, consider:
# - Reducing batch size
# - Using smaller test datasets
# - Increasing system memory
```

### Getting Help

If you encounter issues:

1. Check the [Troubleshooting Guide](troubleshooting.md)
2. Review error messages carefully
3. Check system requirements
4. Create an issue on GitHub with:
   - Operating system
   - Python version
   - Error message
   - Steps to reproduce

## Next Steps

After successful installation:

1. Read the [Quick Start Guide](quick-start.md)
2. Try the [Probe Design Workflow](probe-design-workflow.md)
3. Explore [Configuration Options](configuration.md)

## Uninstallation

To remove the installation:

```bash
# Deactivate environment
conda deactivate

# Remove conda environment
conda env remove -n probe_design

# Remove repository (optional)
rm -rf probe_design
```
