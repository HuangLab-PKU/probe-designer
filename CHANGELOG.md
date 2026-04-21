# Changelog

## 0.2.1 (2026-04-21)

### Behavior change

- `probe_designer.scoring.select_top_n_with_gap` switched from greedy-by-score
  to **round-robin-across-clusters**. Binding sites tend to come in dense
  clusters (a "peak" of adjacent high-scoring sites); the old greedy algorithm
  — a literal port of the webapp's pre-0.2.1 `_auto_select_top_n` — would pile
  its picks into the densest peak and under-cover the rest of the gene. The
  new algorithm buckets sites into regions (default 80 bp, exposed as
  `--region-size` on `probe-design select` and matched to the upstream
  `peak_rank` setting), picks the highest-scoring site per region per round,
  commits each round in score order, and only dips back into already-visited
  regions once every region has contributed. Same `top_n` / `min_gap`
  interface; no API break. Webapp's `_auto_select_top_n` inherits the fix
  because it now delegates to this function.

### Added
- `probe-design merge` subcommand — finishes the CLI surface promised in the
  0.2.0 plan. Rewrite of the 352-line `scripts/merge_results.py` as ~100
  lines using `pandas.merge` for the case-insensitive gene-info filter.
- `probe-design assemble` subcommand — thin Typer wrapper (~40 lines) over
  `ProbeAssembler` that accepts JSON or XLSX via the shared loader.
- `probe_designer.probe_assembly.load_binding_sites(path)` — dispatches on
  file extension (`.json` / `.xlsx`) so the CLI and Python callers share
  one loader. Absorbs the two script-level helpers that used to live in
  `scripts/probe_assemble.py`.
- `examples/README.md` with an index of retained notebooks.

### Changed
- `scripts/find_probes.py` → `scripts/find_single_sequence_probes.py`
  (via `git mv`). The new name pairs symmetrically with
  `find_consensus_probes.py` / `find_specific_probes.py` and hints at the
  strategy.
- Post-review fixes to 0.2.0's Pipeline + webapp:
  - `Pipeline._blast_output_dir` was creating a tempdir per gene on CLI
    dry-runs, never cleaned up — Windows accumulation risk. Now allocates
    one dir per `Pipeline.run()`, cleaned up in a `finally` block.
  - `_DbCachedIsoformProvider` was write-through only (every call
    re-fetched and appended a new row). Now reads the most recent cached
    row first and evicts stale rows on miss.
  - `ingest_binding_sites` was overwriting the Pipeline-computed
    `peak_rank` with the raw insertion index. Now honors
    `site["peak_rank"]` when present.
  - `task_runner._auto_select_top_n` adds an explicit `session.commit()`
    so selection marks persist even if the enclosing pipeline runner
    raises afterwards.
  - `POST /api/design/start` validates strategy as a `Literal`; the
    retired `brute_force` value now returns a clear 422 instead of a 500.
  - Frontend `DesignWizard.vue` was still shipping `brute_force` as the
    default and `consensus` (never a valid backend value) as an option;
    both replaced with the three valid strategy names.
  - `webapp/backend/app/config.py` default SQLite URL was CWD-dependent
    (`sqlite:///./probedesign.db`) and created duplicate empty DBs in
    every directory where `uvicorn` was launched. Now resolved via
    `Path(__file__)` to `webapp/backend/probedesign.db`.

### Removed
- `scripts/legacy/` — `run_pipeline.py` (replaced by `probe-design design`)
  and an older copy of `probe_assemble.py`.
- Broken notebooks `1. binding_searcher_NCBI.ipynb` and
  `1. binding_searcher_ensembl.ipynb` — both imported the renamed `lib.*`
  modules; fully superseded by `probe-design design` with
  `--database {ncbi,ensembl}`.
- Stray `notebooks/random_seq_filtered.txt` data file.

### Moved
- `notebooks/` → `examples/`. Output cells stripped, saving ~56KB.

## 0.2.0 (2026-04-21)

### Breaking changes
- `BruteForceStrategy` renamed to `SingleSequenceStrategy`.
  CLI `--strategy brute_force` and `search_strategy: brute_force` in YAML
  are no longer accepted; use `single_sequence`. No backward-compat alias.
- `configs/config_bruteforce.yaml` renamed to `configs/config_single_sequence.yaml`.
- Hardcoded Entrez email + API key removed from `database.py`. Set
  `ENTREZ_EMAIL` / `ENTREZ_API_KEY` in the environment or in `designer/.env`
  (copy `.env.example`). If you use the previously-committed API key, rotate
  it on ncbi.nlm.nih.gov/account.

### Added
- `probe-design` console script (`[project.scripts]`). Subcommands:
  `design`, `score`, `select`, `validate`. Auto-invokes scoring + peak_rank +
  top-N-with-gap selection so CLI output finally matches webapp quality.
- `probe_designer.pipeline.Pipeline`: central orchestration class shared by
  CLI and webapp, with pluggable `ProgressHook` / `ExistingTargetsHook` /
  `IsoformProvider` Protocols.
- `probe_designer.genome` package: `parse_gtf_for_gene` (with `.pidx.json`
  byte-offset index cache), `fetch_ensembl_isoforms`, `DefaultIsoformProvider`.
- `probe_designer.scoring.scorer.compute_off_target_score`: webapp's
  `10 - n*0.5` formula now in the shared library.
- `probe_designer.scoring.selection.select_top_n_with_gap`: DB-agnostic
  extraction of webapp's `_auto_select_top_n`.
- `probe_designer.config.loader.load_yaml_with_env`: `${VAR}` /
  `${VAR:-default}` expansion with fail-fast on unresolved placeholders.
- `probe_designer.sources.credentials.get_entrez_credentials`: env-var /
  `.env` lookup for NCBI Entrez credentials.
- `.env.example` scaffolded; `.env` variants in `.gitignore`.
- 172+ tests (138 unit + existing mutation/ncbi tests) with Phase 0
  characterization tests locking scoring / thermal_filter / single_sequence
  search behavior.

### Changed
- `filtering.py` promoted to `filtering/` subpackage (no internal class
  decomposition yet; that is deferred to a focused follow-up). Dead
  `_optimize_subsequence` binary-search selection deleted.
- `config.py`, `utils.py`, `scoring.py` promoted to subpackages. All
  previous `from probe_designer.<name> import ...` imports still resolve.
- Webapp `services/pipeline.py` is now a thin adapter: injects
  DB-backed `_DbExistingTargetsHook` + `_DbCachedIsoformProvider` and
  delegates to `probe_designer.pipeline.Pipeline`.
- Webapp `services/genome.py` delegates heavy parsing to designer's
  `probe_designer.genome`; only the DB cache wrapper remains here.
- Webapp `services/ingest.py:252-256` off-target formula replaced with
  a call to `compute_off_target_score`.
- Webapp `services/task_runner.py::_auto_select_top_n` now delegates
  selection semantics to `select_top_n_with_gap`.
- Webapp `api/design.py` strategy default: `brute_force` → `single_sequence`.

### Deprecated
- `scripts/find_probes.py`, `scripts/find_consensus_probes.py`,
  `scripts/find_specific_probes.py` replaced with thin shims that print
  DeprecationWarning and forward to `probe-design design`.
  Scheduled for removal in 0.3.0.

### Removed
- Base `SearchStrategy._optimize_subsequence` (dead code; only called by
  `probe_designer/legacy` and notebooks).
- Duplicate `save_config` method in `ConfigManager` (was silently
  shadowed by the second definition).
- Hardcoded Entrez email + API key in `database.py`.

## 0.1.0
- Initial release.
