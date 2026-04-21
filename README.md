# Probe Designer

Padlock probe design pipeline for spatial transcriptomics and multiplex RNA FISH.
Searches for optimal binding sites on mRNA sequences, scores them for thermodynamics
and off-target risk, and returns a spatially-diverse top-N selection.

**Current version: 0.2.0** · See [CHANGELOG.md](CHANGELOG.md) for release notes.

---

## Install

```bash
# Create the conda env (name: probe-design, note the hyphen)
mamba env create -f environment.yml
mamba activate probe-design

# Install in editable mode
pip install -e .
```

After install, the `probe-design` CLI is on your PATH.

Copy `.env.example` to `.env` and fill in your NCBI credentials (or leave unset —
Entrez then runs with tighter rate limits):

```bash
cp .env.example .env
# edit .env and set ENTREZ_EMAIL / ENTREZ_API_KEY
```

---

## Quick Start

### One-shot design run

```bash
probe-design design \
    --genes-file data/genes.txt \
    --strategy isoform_consensus \
    --species human \
    --output results/my_project/
```

Output:
- `selected_binding_sites.json` — top-N probes per gene (spatially spaced)
- `scored_binding_sites.json` — all filtered sites with quality scores
- `binding_sites.json` — raw search results (pre-BLAST)
- `config_used.json` — frozen config snapshot for reproducibility
- `gene_report.json` — per-gene diagnostics (source, errors, isoform count)

### Subcommands

| Command | Purpose |
|---|---|
| `probe-design design`   | Full pipeline: search → filter → BLAST → score → peak_rank → top-N |
| `probe-design score`    | Re-score an existing `binding_sites.json` (no BLAST) |
| `probe-design select`   | Apply top-N + min-gap to a scored JSON |
| `probe-design validate` | Dry-run: expand `${VAR}` in YAML, hydrate config, exit non-zero on error |

Run `probe-design <cmd> --help` for per-command flags.

---

## Strategies

| Strategy | Input | Use case |
|---|---|---|
| `single_sequence` | One sequence per gene (NCBI mRNA, or a full-length Ensembl isoform) | Simple / when only one transcript matters |
| `isoform_consensus` | All isoforms of the gene | Design probes that bind across the **most** isoforms |
| `isoform_specific` | All isoforms of the gene | Design probes **unique to** one isoform |

All three feed into the same post-processing: thermal filter → BLAST specificity →
`compute_target_score` (isoform coverage, BLAST support, Tm, GC clamp, ΔG) →
`peak_rank` (region-aware ranking) → `select_top_n_with_gap` (spatial spacing).

> **Rename note:** `brute_force` was renamed to `single_sequence` in 0.2.0
> without an alias. If you have scripts / YAML using `brute_force`, update them.

---

## Configuration

Configs live in `configs/`:

| File | Purpose |
|---|---|
| `config_consensus.yaml` | isoform_consensus defaults (recommended) |
| `config_specific.yaml`  | isoform_specific defaults |
| `config_single_sequence.yaml` | single_sequence defaults |
| `config_template.yaml`  | annotated template with all keys documented |
| `species_config.json`   | per-species genome FASTA paths + taxonomy |

YAML supports `${VAR}` and `${VAR:-default}` expansion, resolved at load time
from your shell env or `designer/.env`. Run `probe-design validate --config <file>`
to preview the fully-resolved config before running a pipeline.

---

## Library API

For Python callers (e.g. the webapp), the central entry point is
`probe_designer.pipeline.Pipeline`:

```python
from probe_designer.config import ConfigManager
from probe_designer.pipeline import Pipeline
from probe_designer.genome import DefaultIsoformProvider

cfg = ConfigManager("configs/config_consensus.yaml", species="human")
pipeline = Pipeline(
    cfg,
    isoform_provider=DefaultIsoformProvider(gtf_path="/path/to/gencode.gtf"),
)
result = pipeline.run(["ACTB", "GAPDH"], strategy="isoform_consensus")

for gene, gene_result in result.per_gene.items():
    print(gene, len(gene_result.sites), "errors:", gene_result.errors)
```

`Pipeline` accepts three optional hooks:
- `progress: ProgressHook` — per-stage callbacks (CLI uses tqdm/logging; webapp writes DB)
- `existing_targets: ExistingTargetsHook` — DB-backed lookup for already-designed sites
- `isoform_provider: IsoformProvider` — GTF/Ensembl chain (webapp wraps this with a DB cache)

See `probe_designer/pipeline/hooks.py` for the Protocol definitions.

---

## Project Layout

```
designer/
  probe_designer/          Python package
    cli/                   Typer CLI (app.py + design/score/select/validate)
    pipeline/              Pipeline class + hooks + result dataclasses
    config/                ConfigManager + env-var YAML loader
    filtering/             Thermal + BLAST + specificity filters
    search_strategies.py   SingleSequence + IsoformConsensus + IsoformSpecific
    scoring/               compute_target_score, peak_rank, select_top_n_with_gap
    genome/                GTF parser, Ensembl client, IsoformProvider
    sources/               Entrez credentials (env-based)
    database.py            Unified Ensembl / NCBI sequence interface
    probe_assembly.py      Backbone attachment
    ext/                   Domain-specific extensions (mutation_probe, tcr_probe)
    legacy/                Pre-0.1 reference code (used by TCR skill; do not delete)
  configs/                 YAML configs + species_config.json
  scripts/                 Legacy deprecation shims (forward to probe-design CLI)
  tests/                   pytest suite (138 passed + 1 skipped)
  docs/                    Installation / development guides
  CHANGELOG.md             Release notes
  .env.example             Credential template (copy to .env)
```

---

## Testing

```bash
python -m pytest tests/ -v
python -m pytest tests/ --cov=probe_designer
```

See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) for contribution guidelines
and [docs/installation.md](docs/installation.md) for detailed setup / troubleshooting.
