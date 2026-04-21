# Probe Design Pipeline Scripts

> **Deprecation notice (0.2.0):** The three `find_*_probes.py` scripts in this
> directory are now thin shims that forward to the new ``probe-design`` CLI.
> They print a DeprecationWarning on start and will be removed in 0.3.0.
> Migrate to ``probe-design`` ASAP — it is the single entry point installed
> alongside ``probe-designer`` via ``pip install -e .``.

## Preferred interface

```bash
probe-design design --genes-file genes.txt --strategy isoform_consensus
probe-design design --genes-file genes.txt --strategy isoform_specific
probe-design design --genes-file genes.txt --strategy single_sequence   # was brute_force
probe-design score  --input out/binding_sites.json
probe-design select --input out/scored.json --top-n 3 --min-gap 40
probe-design validate --config configs/config_consensus.yaml
```

The CLI auto-invokes scoring + peak_rank + top-N spatial selection for every
design run, so CLI output now matches the webapp's quality metric.

## Legacy scripts (deprecated)

### 1. Find Binding Sites

#### `find_consensus_probes.py` *(deprecated: use `probe-design design --strategy isoform_consensus`)*
Finds consensus binding sites (probes that bind to maximum number of isoforms per gene).

#### `find_specific_probes.py` *(deprecated: use `probe-design design --strategy isoform_specific`)*
Finds isoform-specific binding sites per gene.

#### `find_probes.py` *(deprecated: use `probe-design design --strategy single_sequence`)*
Finds binding sites by scanning one input sequence per gene (mRNA from NCBI,
or a full-length Ensembl isoform).

**Common Parameters:**
- `--config`: Configuration file path (YAML format)
- `--genes_file`: Gene list file path (can also be specified in config)
- `--output_dir`: Output directory (can also be specified in config)
- `--progress`: Show progress bars (optional)

**Output:**
- `filtered_binding_sites.xlsx`: Filtered binding sites in Excel format
- `binding_sites.json`: Raw binding sites data
- `filtering_stats.json`: Filtering statistics

### 2. Merge Results

#### `merge_results.py`
Merges binding sites from multiple result directories, filters by gene_info, deduplicates, and sorts.

**Usage:**
```bash
python merge_results.py --results-dir results/ --gene-info gene_info.xlsx --output merged_binding_sites.xlsx
```

**Parameters:**
- `--results-dir` / `-r`: Results directory (will traverse all subdirectories)
- `--gene-info` / `-g`: Gene info Excel file (must contain `gene_name` and `No.` columns)
- `--output` / `-o`: Output file path (merged binding sites Excel file)
- `--missing-output` / `-m`: Missing genes output file path (optional, default: same directory as output)

**Features:**
- Automatically finds all `filtered_binding_sites.xlsx` files in subdirectories
- Filters by genes in `gene_info`
- Deduplicates based on `gene_name`, `st`, `en`
- Sorts by `gene_info` order, then by `st`, `en`
- Reports missing genes (genes in `gene_info` but not in results)

**Output:**
- Merged binding sites Excel file
- `missing_genes.txt`: List of genes not found in results (if any)

### 3. Assemble Probes

#### `probe_assemble.py`
Assembles probes from binding sites and backbone sequences.

**Usage:**
```bash
python probe_assemble.py --binding_sites merged_binding_sites.xlsx --gene_info gene_info.xlsx --backbone_file backbones.xlsx --output_dir output/
```

**Parameters:**
- `--binding_sites`: Binding sites file path (JSON or Excel)
- `--gene_info`: Gene info Excel file (must contain `gene_name` and `No.` columns)
- `--backbone_file`: Backbone Excel file (must contain `No.` and `Sequence` columns)
- `--output_dir`: Output directory
- `--ilock`: Add iLock modification (optional: `3prime`, `5prime`, or `both`)

**Features:**
- Supports per-gene or per-binding-site probe ID assignment
- Validates probe IDs against backbone file
- Optional iLock modification

**Output:**
- `probes.xlsx`: Probe sequences in Excel format
- `probes.fasta`: Probe sequences in FASTA format
- `probes.json`: Probe data in JSON format
- `probe_stats.json`: Probe statistics

## Complete Workflow Example

```bash
# Step 1: Find binding sites
python find_consensus_probes.py --config configs/config_consensus.yaml

# Step 2: Merge results from multiple runs
python merge_results.py --results-dir results/ --gene-info gene_info.xlsx --output merged_binding_sites.xlsx

# Step 3: Assemble probes
python probe_assemble.py --binding_sites merged_binding_sites.xlsx --gene_info gene_info.xlsx --backbone_file backbones.xlsx --output_dir probes/
```

## Configuration Files

Configuration files are in YAML format. See `configs/` directory for examples:

- `config_consensus.yaml`: Configuration for consensus probe design
- `config_specific.yaml`: Configuration for specific probe design
- `config_bruteforce.yaml`: Configuration for brute force search
- `config_template.yaml`: Template with detailed comments

## File Formats

### Gene Info File
Excel file with columns:
- `gene_name`: Gene name
- `No.`: Probe ID number (corresponds to backbone No.)

### Backbone File
Excel file with columns:
- `No.`: Backbone ID number
- `Sequence`: Backbone sequence

### Binding Sites File
Excel file with columns:
- `gene_name`: Gene name
- `arm3`: 3' arm sequence
- `arm5`: 5' arm sequence
- `tm`: Melting temperature (°C)
- `tm3`: 3' arm Melting temperature (°C)
- `tm5`: 5' arm Melting temperature (°C)
- `st`: Genomic start position
- `en`: Genomic end position
- `g_content`: G content (0.0-1.0)
- (and other optional columns)