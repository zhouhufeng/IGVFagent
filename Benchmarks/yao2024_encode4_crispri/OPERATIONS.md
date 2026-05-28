# OPERATIONS — Yao 2024 ENCODE4 noncoding CRISPRi screens

For shared prerequisites + the two execution paths + concordance scoring,
see `Benchmarks/OPERATIONS_GUIDE.md`. Paper context in this directory's
`README.md`.

## Quick run — online metadata pull, ~30 s

```bash
bash Benchmarks/yao2024_encode4_crispri/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark yao2024_encode4_crispri
```

## What `run.sh` does

```bash
.venv/bin/igvfagent encode retrieve \
    --assay-type "CRISPR screen" \
    --label yao2024_encode4
```

Hits `www.encodeproject.org/search/?type=Experiment&assay_title=CRISPR%20screen&limit=all`
and writes a manifest of every catalogued CRISPRi/CRISPRa screen
ENCODE has (~108 expected, per the paper).

## Where artefacts will land

`encode retrieve` uses the **flat-file** convention:

```
Docs/ENCODE/<ts>_yao2024_encode4_<purpose>.{md,csv,json}
Data/Manifests/ENCODE/<ts>_yao2024_encode4_manifest.csv
Data/<ts>_encode_yao2024_encode4_summary.json
```

The concordance scorer's `extra_search_dirs` covers
`Data/Manifests/ENCODE` and `Data` for this paper.

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | A `manifest.json` (or `*_manifest.csv`) artefact exists, non-empty | ENCODE search returned and persisted |
| 2 | `n_experiments` JSON field falls in [50, 5000] (expect ≈ 108) | Roughly the right number of screens enumerated |

If the `n_experiments` field is at a different JSON path than expected,
adjust `primary_artefact` + the `range` check's `path` in
`expected.json` after one successful run by looking at the actual
written JSON.

## Ground-truth spot-checks

| Signal | Expected |
|---|---|
| Total CRISPR screens enumerated | Approximately 108 (paper Methods) |
| Cell lines covered | K562, HepG2, Jurkat |
| Total perturbed DNA | 24.85 Mb across all screens |
| Specifically GATA1 locus screens in K562 | Should include the +24/+58 kb HS-sites region |

To verify the GATA1 locus is in the manifest, after `run.sh` finishes:

```bash
grep -i "GATA1\|chr_X\b\|K562" \
    Data/Manifests/ENCODE/*yao2024_encode4*_manifest.csv | head -5
```

## Local analytical step (optional)

`run.sh` only does the metadata pull. To actually re-run a CRISPR
screen analysis on Yao 2024's K562 GATA1 data, you'd need to:

1. Pick an accession from the manifest (e.g. `ENCSR…` for the K562
   GATA1 screen).
2. Download the guide-count matrix:
   ```bash
   .venv/bin/igvfagent encode download <accession>
   ```
3. Run the IGVFagent log-FC pipeline:
   ```bash
   .venv/bin/igvfagent crispri analyze-local --input <counts.tsv> --label yao2024_k562_gata1
   ```

This is **not** automated in `run.sh` because the ENCODE counts table
format varies per screen and IGVFagent's `crispri analyze-local`
schema is opinionated about column names. For an end-to-end CRISPR
screen reproducibility test, see `martyn2025_variant_flowfish/` —
that one wires the analysis step.

## Running through the UI

`encode_retrieve` is **not yet registered** as an LLM tool. Closest
registered alternative is `portal_search` against the IGVF portal,
which doesn't cover ENCODE. Two paths:

a) Use the shell `run.sh` (recommended).

b) Ask the LLM to call `chipatlas_search --query "CRISPR K562"` as a
   loosely-equivalent enumeration of public CRISPR-screen experiments.

```
List public CRISPR screens for K562 cells via the chipatlas_search
tool with query="CRISPR K562 GATA1", genome="hg38", limit=20.
This is a proxy for the ENCODE4 noncoding CRISPRi catalogue from
Yao 2024 (doi:10.1038/s41592-024-02216-7).

Report: total hits, top antigens/cell types, and whether any
experiments target the GATA1 locus.
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `encode retrieve` returns 0 experiments | `--assay-type` value mismatch (ENCODE renamed the field) | Check available assay titles with `igvfagent encode retrieve --help`; try `--assay "CRISPR screen"` |
| Manifest CSV exists but `n_experiments` JSON path mismatch | The summary JSON's schema differs from the placeholder | Open the actual summary JSON and update the `path` in `expected.json` |

## License + provenance

* **Paper data**: ENCODE portal — public.
* **Code**: IGVFagent Apache-2.0.
* **Citation**: Yao D et al. *Nat Methods* **21**: 1980–1992 (2024).
  doi:10.1038/s41592-024-02216-7 · PMID:38504114
