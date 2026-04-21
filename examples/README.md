# Examples

Reference notebooks demonstrating ad-hoc analyses and auxiliary tools around
the `probe_designer` package. These are **not** the primary interface —
production probe design runs should use the `probe-design` CLI (see
[../README.md](../README.md)).

All notebooks here have their output cells stripped; re-run them locally to
regenerate figures / tables.

| Notebook | Uses `probe_designer` | Purpose |
|---|---|---|
| `overlap_detection.ipynb` | Yes | Detect genomic overlap among candidate binding sites within a gene |
| `binding_merge.ipynb` | No | Ad-hoc XLSX merge demo (the CLI's `probe-design merge` now supersedes this) |
| `probe_stitch.ipynb` | No | Manual probe stitching sandbox |
| `barcode_generator.ipynb` | No | Generate random DNA barcodes for PRISM/SPRINTseq rounds |
| `random_seq_generator.ipynb` | No | Generate random DNA sequences for controls; inputs to `random_seq_filtered.txt` downstream |

## Running

```bash
mamba activate probe-design
jupyter lab examples/
```

## Related — removed in v0.2.1

- `notebooks/1. binding_searcher_NCBI.ipynb`, `binding_searcher_ensembl.ipynb` —
  deleted because they imported `lib.*` (renamed to `probe_designer.*` years
  ago) and were fully superseded by
  `probe-design design --strategy single_sequence --database {ncbi,ensembl}`.
- `notebooks/random_seq_filtered.txt` — deleted; was stray data output left
  behind by `random_seq_generator.ipynb`. Regenerate locally as needed.
