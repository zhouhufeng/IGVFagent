# southard2024_comprehensive_transcription — Operations

Shared prerequisites live in `Benchmarks/OPERATIONS_GUIDE.md`; this file only covers what is specific to this paper.

## 1. Inputs

| Variable | Value | Origin | Override env var |
|---|---|---|---|
| `MODALITY` | `crispr-screen` | route default | `SOUTHARD2024_MODALITY` |

## 2. Run

```bash
bash Benchmarks/southard2024_comprehensive_transcription/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark southard2024_comprehensive_transcription
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
