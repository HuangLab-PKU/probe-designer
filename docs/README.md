# Documentation Index

## For users

- **[../README.md](../README.md)** — Top-level overview + CLI quick start
- **[installation.md](installation.md)** — Environment setup, credentials, troubleshooting
- **[migration_0.2.md](migration_0.2.md)** — Upgrading from 0.1.x to 0.2.0

## For contributors

- **[DEVELOPMENT.md](DEVELOPMENT.md)** — Commit conventions, branching, PR flow
- **[../CHANGELOG.md](../CHANGELOG.md)** — Release notes (authoritative; older
  copies that once lived here have been removed)

## For hands-on reference

- **[../configs/README.md](../configs/README.md)** — YAML config schema (zh)
- **[../scripts/README.md](../scripts/README.md)** — Legacy script shims and their
  new `probe-design` equivalents

## Inline API docs

Every public module has a module-level docstring. Start from:
- `probe_designer.pipeline.Pipeline` — the central orchestration API
- `probe_designer.scoring` — `compute_target_score`, `peak_rank`, `select_top_n_with_gap`
- `probe_designer.genome` — GTF + Ensembl isoform sources
- `probe_designer.config.loader` — `load_yaml_with_env` (`${VAR}` / `${VAR:-default}`)

```python
python -c "help(<module>)"
```
