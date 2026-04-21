# Installation Guide

## Requirements

| Component | Version |
|---|---|
| OS | Windows 10/11, macOS 10.14+, or Linux |
| Python | 3.10+ |
| Conda/Mamba | Miniforge or Miniconda |
| BLAST+ (optional) | 2.13+ for local BLAST; online fallback works without it |
| ViennaRNA (optional) | 2.7+ for RNA secondary-structure check |

Notes:
- ViennaRNA wheels are on PyPI for Windows, Linux, and macOS
  (`pip install ViennaRNA`).
- The pipeline will auto-fall back to online NCBI BLAST if BLAST+ isn't on PATH.

---

## Setup

```bash
# 1. Clone the repo (or initialize submodules from the parent project)
cd designer/

# 2. Create the conda environment (env name: probe-design, with a hyphen)
mamba env create -f environment.yml
mamba activate probe-design

# 3. Install the package in editable mode. This also registers the
#    `probe-design` console script on PATH.
pip install -e .

# 4. (Optional) install optional extras
pip install -e ".[structure]"   # ViennaRNA for RNA MFE
pip install -e ".[dev]"         # pytest + responses + coverage
```

Verify:

```bash
probe-design --help
python -c "import probe_designer; print(probe_designer.__version__)"
```

---

## Credentials

The pipeline uses NCBI Entrez when fetching mRNA sequences. Without
credentials, NCBI rate-limits aggressively.

```bash
# From the designer/ root:
cp .env.example .env

# Edit .env:
#   ENTREZ_EMAIL=your-email@example.com
#   ENTREZ_API_KEY=<your NCBI API key>
```

The key file is `.gitignored`; only `.env.example` is tracked.

To get a key: log into [ncbi.nlm.nih.gov/account](https://ncbi.nlm.nih.gov/account)
and create one in **API Key Management**.

Never commit credentials. YAML configs can reference env vars via
`${ENTREZ_EMAIL}` / `${ENTREZ_API_KEY}` placeholders — these expand at load
time. Run `probe-design validate --config <your.yaml>` to preview.

---

## Reference Genome

Set `genome_fasta_path` under the `genome:` section of your YAML config, or
use `--genome-fasta /path/to/GRCh38.fa` on the CLI. Without a local genome,
the pipeline falls back to Ensembl REST for sequence retrieval (slower, 15s
timeout per query).

Per-species defaults live in `configs/species_config.json`:

```json
{
  "species": {
    "human":  {"genome_fasta_path": "/data/genomes/GRCh38.fa", ...},
    "mouse":  {"genome_fasta_path": "/data/genomes/GRCm39.fa", ...}
  }
}
```

For local GTF (used by isoform strategies, offline), set
`genome.gtf_path` or pass `--gtf /path/to/gencode.gtf`. First parse builds
a `<gtf>.pidx.json` byte-offset index so subsequent gene lookups are O(1).

---

## Troubleshooting

### `probe-design: command not found`

The console script wasn't registered. Re-run `pip install -e .` from the
`designer/` directory. Confirm with:

```bash
pip show probe-designer
# Editable project location: .../designer
```

### `ConfigError: Unresolved environment variable: ${ENTREZ_EMAIL}`

Your YAML has `${ENTREZ_EMAIL}` but no env var + no default. Either:
- Set it in `.env` / your shell, or
- Change YAML to `${ENTREZ_EMAIL:-}` (empty default if truly optional), or
- Remove the placeholder from the YAML.

### ViennaRNA missing

Thermal filter still runs; only the optional RNA secondary-structure check
is skipped. To enable:

```bash
pip install ViennaRNA
```

Set `filter.check_rna_structure: true` in your config.

### BLAST+ not on PATH

The online fallback kicks in automatically. If both fail, the pipeline
returns pre-BLAST thermal results only and logs a warning — you'll see
`"BLAST failed (...), continuing with pre-BLAST results"` in the log.
To skip BLAST entirely (e.g., in dry runs), pass `--skip-blast`.

### Genome FASTA not found

Check `configs/species_config.json` for the path on your species. For
Windows paths with spaces, use forward slashes: `C:/data/genome.fa`.

### NCBI 429 rate limits

Either set `ENTREZ_API_KEY` in `.env` (raises limit from 3 to 10 req/s),
or reduce concurrency via `blast.concurrency` in your config.

---

## Uninstall

```bash
mamba deactivate
mamba env remove -n probe-design
```
