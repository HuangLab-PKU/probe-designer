# Migration Guide: 0.1.x → 0.2.0

0.2.0 is a structural refactor with a small number of breaking changes.
This page is the shortest path from a working 0.1.x setup to a working
0.2.0 setup.

## TL;DR Checklist

- [ ] Rotate your NCBI API key (the hardcoded one in 0.1 was public).
- [ ] `cp .env.example .env` and fill in `ENTREZ_EMAIL` / `ENTREZ_API_KEY`.
- [ ] Find every `brute_force` in your YAML / scripts / DB seeds → change to `single_sequence`.
- [ ] Rename `configs/config_bruteforce.yaml` → `configs/config_single_sequence.yaml` (if you kept your own copy).
- [ ] Reinstall: `pip install -e .` to register the new `probe-design` CLI.
- [ ] Migrate any shell scripts from `python scripts/find_*_probes.py` to `probe-design design --strategy <...>`.

The shims in `scripts/` still work for one release — they print a
DeprecationWarning and forward to the new CLI — but they will be removed
in 0.3.0.

---

## Breaking changes, in detail

### 1. `BruteForceStrategy` → `SingleSequenceStrategy`

| Before | After |
|---|---|
| `from probe_designer.search_strategies import BruteForceStrategy` | `from probe_designer.search_strategies import SingleSequenceStrategy` |
| `--strategy brute_force` on CLI | `--strategy single_sequence` |
| `search_strategy: brute_force` in YAML | `search_strategy: single_sequence` |
| `"strategy": "brute_force"` in output JSON | `"strategy": "single_sequence"` |

**No alias** is provided — the old names raise `ValueError`. The webapp
API endpoint `POST /api/design/start` now returns HTTP 422 with a clear
message when given the legacy value.

**Why:** `BruteForce` described the implementation (exhaustive scan), not
the input contract. The strategy takes a single sequence per gene and
scans it for binding sites. The new name pairs symmetrically with
`isoform_consensus` and `isoform_specific`, which both consume the
multi-isoform pre-processing step.

### 2. Hardcoded Entrez credentials removed

**Before:** `probe_designer/database.py` lines 28-29 contained a real
email and API key committed to the public submodule.

**After:** `DatabaseInterface.__init__` calls
`probe_designer.sources.credentials.get_entrez_credentials()` which reads
`ENTREZ_EMAIL` and `ENTREZ_API_KEY` from the environment or from a
`designer/.env` file (via `python-dotenv`).

**Action required:**
1. Rotate the API key that was publicly committed (log into
   `ncbi.nlm.nih.gov/account` → API Key Management → regenerate).
2. `cp designer/.env.example designer/.env`, fill in your new
   email + key. `.env` is gitignored.
3. If you were running without credentials, nothing breaks —
   `get_entrez_credentials(required=False)` returns `(None, None)` and
   Entrez uses anonymous rate limits.

### 3. CLI entry point

0.1.x had five independent scripts under `scripts/`. 0.2.0 adds a single
console script.

| 0.1.x | 0.2.0 |
|---|---|
| `python scripts/find_probes.py --genes_file X` | `probe-design design --strategy single_sequence --genes-file X` |
| `python scripts/find_consensus_probes.py --genes_file X` | `probe-design design --strategy isoform_consensus --genes-file X` |
| `python scripts/find_specific_probes.py --genes_file X` | `probe-design design --strategy isoform_specific --genes-file X` |
| (manual post-processing) | `probe-design score --input <json>`, `probe-design select --input <json>` |

`pip install -e .` puts `probe-design` on PATH. Run `probe-design --help`.

The old scripts still work as thin deprecation shims — they translate
their old flags (`--genes_file` with underscore) to the new CLI's flags
(`--genes-file` with hyphen) and forward. They print a `DeprecationWarning`
and are scheduled for removal in 0.3.0.

### 4. Default ranking + selection now applied automatically

0.1.x CLI saved unsorted, unscored binding sites. 0.2.0 always scores
(`compute_target_score`), ranks spatially (`peak_rank`), and picks
top-N with a minimum gap (`select_top_n_with_gap`) before writing output.

The webapp was already doing this; the CLI now matches.

Control the final selection with:

```bash
probe-design design --top-n 3 --min-gap 40 ...
```

### 5. YAML env-var expansion

YAML values now support `${VAR}` and `${VAR:-default}` expansion,
resolved at load time from your shell or `designer/.env`. Unresolved
placeholders raise `ConfigError` (fail fast).

Example:

```yaml
database:
  email: ${ENTREZ_EMAIL}          # required; error if unset
  api_key: ${ENTREZ_API_KEY:-}    # optional; empty if unset
```

Preview expansion:

```bash
probe-design validate --config configs/config_consensus.yaml
```

### 6. Subpackage moves (import-compatible)

Some internal modules became subpackages. The public imports are
**unchanged**:

```python
from probe_designer.config import ConfigManager        # still works
from probe_designer.scoring import compute_target_score # still works
from probe_designer.filtering import SequenceFilter    # still works
from probe_designer.utils import load_gene_list        # still works
```

Internally, `config.py`, `scoring.py`, `utils.py`, `filtering.py` are now
directories with `__init__.py`. You only notice if you were importing
private names (e.g., `from probe_designer.scoring import _some_private`).

New public surface:

```python
from probe_designer.pipeline import Pipeline, GeneResult, PipelineResult
from probe_designer.genome import parse_gtf_for_gene, DefaultIsoformProvider
from probe_designer.scoring import select_top_n_with_gap, compute_off_target_score
from probe_designer.config.loader import load_yaml_with_env
from probe_designer.utils.errors import ConfigError
```

---

## Webapp-side migration

If you run the webapp:

- `POST /api/design/start` body: change `"strategy": "brute_force"` to
  `"strategy": "single_sequence"`. The endpoint now validates this at the
  API boundary and returns a 422 with a migration hint if you miss one.
- `DATABASE_URL` default previously resolved to `./probedesign.db` (CWD-dependent;
  produced duplicate empty DBs). Now pinned to
  `webapp/backend/probedesign.db` via an absolute path. `DATABASE_URL`
  env var still wins for production (PostgreSQL in docker-compose).
- Clean up any stray `probedesign.db` files outside `webapp/backend/`.
  They're all empty-schema from the CWD bug — safe to delete.

---

## 0.2.1 addenda

- **Selection semantics changed** (webapp-facing). `select_top_n_with_gap`
  (and therefore the webapp's `_auto_select_top_n`) now does
  round-robin-across-clusters instead of greedy-by-score. In practice this
  means when binding sites cluster in dense peaks, the first picks spread
  across peaks for gene-wide coverage rather than piling into the densest
  peak. Same API; no action needed unless you're writing tests that
  assumed the old greedy output. See CHANGELOG 0.2.1 for the full
  rationale and the new `--region-size` flag on `probe-design select`.
- `scripts/find_probes.py` was renamed to `scripts/find_single_sequence_probes.py`
  so the three find_* shims share a consistent naming scheme. The shim still
  accepts all legacy arguments; only the file name changed.
- Two new CLI subcommands finish the surface promised in the 0.2.0 plan:
  `probe-design merge` (replaces `scripts/merge_results.py`) and
  `probe-design assemble` (replaces `scripts/probe_assemble.py`). Both
  scripts now forward via the usual DeprecationWarning shim.
- `scripts/legacy/` (the two pre-refactor files) was deleted.
- `notebooks/` was renamed to `examples/`; two broken notebooks that imported
  the long-renamed `lib.*` modules were deleted.

## What didn't change

- The underlying search, filter, and BLAST algorithms are byte-identical.
- YAML config keys and their meanings (`filter.min_tm`, `blast.evalue`, etc.)
  are unchanged. Only the `search.search_strategy` value changed.
- Output JSON schema (field names, types) is unchanged except the
  `"strategy"` value.
- `probe_designer/legacy/` is intentionally preserved; the
  `tcr-padlock-design` skill imports from there.

---

## Need help?

- Release notes: [../CHANGELOG.md](../CHANGELOG.md)
- Full installation: [installation.md](installation.md)
- Config reference: [../configs/README.md](../configs/README.md) (zh)
