# Changelog

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
