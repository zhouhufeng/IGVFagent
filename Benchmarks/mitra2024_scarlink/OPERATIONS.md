# OPERATIONS — Mitra 2024 SCARlink multiome

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Quick run — online metadata pull, ~30 s

```bash
bash Benchmarks/mitra2024_scarlink/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark mitra2024_scarlink
```

## What `run.sh` does

```bash
.venv/bin/igvfagent multiome retrieve --count 5 --label mitra2024_scarlink
```

Queries the IGVF Portal for the 5 most recent public 10x Multiome
AnalysisSets and writes a manifest. The `multiome` skill is the
clean-room IGVFagent reimplementation of the Signac/Seurat/quadbio
WNN + peak-to-gene + chromVAR + CSS + MultiVI workflow.

## Where artefacts will land

```
Data/IGVF/10xMultiome/Metadata/<ts>_*mitra2024_scarlink*_full_metadata.json
Data/Manifests/Multiome10x/<ts>_*mitra2024_scarlink*_analysis_sets.csv
Data/Manifests/Multiome10x/<ts>_*mitra2024_scarlink*_files.csv
Data/Manifests/Multiome10x/<ts>_*mitra2024_scarlink*_samples.csv
Docs/Multiome10x/<ts>_*mitra2024_scarlink*_report.md
```

The concordance scorer's `extra_search_dirs` should cover
`Data/Manifests/Multiome10x` and `Data/IGVF/10xMultiome` if you want
checks beyond the default.

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | A `manifest.json` artefact exists, non-empty | Multiome AnalysisSet enumeration succeeded |

(Current `expected.json` is minimal — extend it once you've inspected
the actual output JSON schema.)

## Ground-truth spot-checks

The paper's published expectations against IGVF data:

| Signal | Expected |
|---|---|
| At least 3 IGVF 10x Multiome AnalysisSets enumerated | The IGVF portal has dozens; `--count 5` returns 5 |
| Each set has both an scRNA and scATAC component | Inspect the `_files.csv` manifest — should list both 10x H5 + ATAC fragments |
| Tissue diversity | Should cover at least 2 tissue/cell-type contexts |

Open the manifest and verify:

```bash
.venv/bin/python -c "
import csv
m = csv.DictReader(open('Data/Manifests/Multiome10x/<ts>_..._analysis_sets.csv'))
for r in m:
    print(r.get('accession'), r.get('preferred_assay_titles'), r.get('samples'))"
```

## Local analytical step (optional but recommended)

The benchmark gets really interesting when you load an actual h5ad /
mtx bundle and run the SCARlink-style pipeline. After `run.sh` writes
the manifest, pick an AnalysisSet from it and:

```bash
# Download the 10x H5 + ATAC fragments
.venv/bin/igvfagent multiome retrieve \
    --count 1 --download-policy processed --max-download-gb 5

# Then the per-cell joint pipeline (the SCARlink analogue)
.venv/bin/igvfagent multiome process-local --input <path-to-bundle>
.venv/bin/igvfagent multiome joint-qc --label mitra2024
.venv/bin/igvfagent multiome lsi --label mitra2024
.venv/bin/igvfagent multiome wnn --label mitra2024
.venv/bin/igvfagent multiome peak2gene --label mitra2024
```

The `peak2gene` output is the direct analogue of SCARlink's per-peak
gene-link table — but using IGVFagent's correlation method
(Signac-style), not SCARlink's Poisson regression. Expect rank-order
concordance Spearman ρ ≥ 0.6 between the two on the same data.

## Ground-truth: lineage-specific enhancers to spot-check

| Lineage | Gene + enhancer | Expected behaviour |
|---|---|---|
| Erythroid | GATA1 enhancer (chrX:48,786,000-48,790,000) | Top correlate in erythroid-lineage cells |
| Monocyte | IRF8 enhancer (chr16:85,932,000-85,936,000) | Top correlate in monocyte-lineage cells |
| Lymphoid | EBF1 enhancer (chr5:158,500,000-158,710,000) | Top correlate in lymphoid-lineage cells |

After `peak2gene` runs, the output TSV will have columns
`(peak_locus, gene_name, correlation, p_value, biological_context)`.
Filter rows where `gene_name == "GATA1"` and check that the top-N rows
by correlation are in cells annotated as erythroid (column varies by
input metadata).

## Running through the UI

`multiome_retrieve` is **not yet registered** as an LLM tool. Two
paths:

a) Use the shell `run.sh` (recommended).

b) Ask the LLM to enumerate via `portal_search`:

```
Find IGVF Portal 10x Multiome AnalysisSets via portal_search:
  - type=AnalysisSet,MeasurementSet
  - field_filters="preferred_assay_titles=10x multiome"
  - limit=10
  - label="mitra2024_scarlink"

Then for each AnalysisSet, call portal_get to fetch the @id of its
associated files. Report a table of (accession, cell-class, tissue,
file count) — a manifest comparable to what mitra2024_scarlink/run.sh
would write.
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Manifest CSV has 0 rows | Portal returned no 10x multiome AnalysisSets | Check `--count` value; try `--count 50` to widen the search |
| Files manifest lists 0 files | The AnalysisSet has no released files yet | Re-run with a different AnalysisSet; the search default is "most recent" which may be embargoed |
| `multiome wnn` errors with "muon not installed" | Optional dep missing | `pip install muon` |
| `peak2gene` runs forever | The peak × cell matrix is very large | Subsample with `--max-peaks 50000 --max-cells 5000` |

## License + provenance

* **Paper data**: GSE194122 — public.
* **Code**: IGVFagent Apache-2.0 over upstream MIT (10XGenomics/analysis_guides + stuart-lab/signac) + quadbio (no LICENSE → clean-room).
* **Citation**: Mitra S et al. *Nat Genet* **56**: 627–636 (2024).
  doi:10.1038/s41588-024-01689-8 · PMC:PMC11018525
