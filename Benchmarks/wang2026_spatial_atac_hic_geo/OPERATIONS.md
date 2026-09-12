# wang2026_spatial_atac_hic_geo — Operations

Shared prerequisites live in `Benchmarks/OPERATIONS_GUIDE.md`; this file only covers what is specific to this paper.

## 1. Inputs

| Variable | Value | Origin | Override env var |
|---|---|---|---|
| `GSE` | `**TODO_VERIFY**` | unresolved | `WANG2026_GSE` |
| `PAIRS` | `Benchmarks/_data/wang2026_spatial_atac_hic_geo/sample.pairs.gz` | local path (must be downloaded) | `WANG2026_PAIRS` |
| `FRAGMENTS` | `Benchmarks/_data/wang2026_spatial_atac_hic_geo/fragments.tsv.gz` | local path (must be downloaded) | `WANG2026_FRAGMENTS` |

## 2. Run

```bash
bash Benchmarks/wang2026_spatial_atac_hic_geo/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark wang2026_spatial_atac_hic_geo
```

Exit 77 means a required input is missing — the message names which.

## 3. Promoting unconfirmed checks

`expected.json` ships the paper's prose-derived numbers as `"confirmed": false` with `"path": "TODO_SET_JSON_PATH"`. To turn one into a real reproducibility claim:

1. Run the chain once and open the artefact named by `primary_artefact`.
2. Find the key that holds the comparable quantity; put its dotted path in `path`.
3. Check the quoted sentence in `provenance.quote` really states that number for that quantity.
4. Set `"confirmed": true` and tighten `min`/`max` to the tolerance you are willing to defend.

## 4. Troubleshooting

See §5 of `Benchmarks/OPERATIONS_GUIDE.md`.
