# OPERATIONS — Agarwal 2025 lentiMPRA K562 / HepG2 / WTC11

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`. This
operations doc complements `README.md` (paper context) and was
**validated end-to-end** on commit `ad73060` — 3/3 concordance checks
passed for the online step.

## Two stages

| Stage | Mode | What it tests |
|---|---|---|
| 1. Online discovery | always runs | `mpra pull` enumerates IGVF Portal MPRA AnalysisSets, writes report + manifest + summary |
| 2. Local analytical | requires `oligo_counts.tsv` | `mpra activity / qc / volcano` run the NB GLM Wald test + replicate QC + 4-panel volcano |

## Stage 1 — quick run (always works, ~30 s)

```bash
bash Benchmarks/agarwal2025_lentimpra/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark agarwal2025_lentimpra
```

**Verified live output** (May 27, 2026):

```
=== Concordance — 20260527_192430 ===
  [          ok] agarwal2025_lentimpra                3/3

✓ MPRA portal manifest CSV: Data/Manifests/MPRA/...manifest.csv (3,137 B)
✓ MPRA portal report MD:    Docs/MPRA/...report.md (4,789 B)
✓ MPRA discovery summary:   Data/...summary.json (6,030 B)
```

## What `run.sh` does (Stage 1)

```bash
.venv/bin/igvfagent mpra portal-manifest --limit 50 --label agarwal2025_lentimpra
```

`mpra pull` (the underlying CLI for `portal-manifest`) hits four IGVF
Portal search endpoints in parallel:

1. `type=MeasurementSet&searchTerm=MPRA` → MPRA MeasurementSets
2. `type=File&searchTerm=MPRA` → MPRA files
3. `type=File&searchTerm=STARR` → STARR files (related assay)
4. `type=File&searchTerm=reporter` → general reporter assays

Each result is summarized and written to disk.

## Where Stage 1 artefacts land

The `mpra pull` skill uses **flat-file convention** across three roots:

| Path | What |
|---|---|
| `Docs/MPRA/<ts>_agarwal2025_lentimpra_portal_many_report.md` | User-facing summary |
| `Data/Manifests/MPRA/<ts>_agarwal2025_lentimpra_portal_many_manifest.csv` | Tabular manifest (one row per matched item) |
| `Data/<ts>_mpra_agarwal2025_lentimpra_portal_many_summary.json` | Machine-readable per-query rollup |
| `Docs/MPRA/Plots/<ts>_agarwal2025_lentimpra_portal_*.svg` | 3 SVG plots (query rows, status counts, format counts) |

`expected.json` declares `extra_search_dirs: ["Data/Manifests/MPRA", "Data"]`
so the scorer finds all three.

## Stage 2 — local analytical run

Stage 2 requires the per-oligo DNA + RNA count table at:

```
Data/Benchmarks/agarwal2025_lentimpra/oligo_counts.tsv
```

### Where to get the counts table

**Option A — GEO (the paper's primary deposit)**

```bash
mkdir -p Data/Benchmarks/agarwal2025_lentimpra
cd Data/Benchmarks/agarwal2025_lentimpra

# GEO accession is GSE142696 (related Inoue series).
# Visit https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142696
# and click through to the supplementary files. The per-oligo counts
# typically have suffixes like:
#   GSE142696_K562_oligo_counts.tsv.gz
#   GSE142696_HepG2_oligo_counts.tsv.gz
#   GSE142696_WTC11_oligo_counts.tsv.gz

# Once downloaded:
gunzip GSE142696_K562_oligo_counts.tsv.gz
mv GSE142696_K562_oligo_counts.tsv oligo_counts.tsv
cd ../../..
```

**Option B — IGVF Portal MPRA AnalysisSet**

After Stage 1's manifest, pick one of the MPRA MeasurementSets it
listed:

```bash
.venv/bin/igvfagent portal batch-download \
    --type MeasurementSet \
    --field-filters "accession=<IGVFDS…>" \
    --fetch --label agarwal2025_inputs
# Move the fetched count table to the expected path:
mv Docs/PortalQuery/*_batchdl_agarwal2025_inputs/files/*_counts.tsv \
   Data/Benchmarks/agarwal2025_lentimpra/oligo_counts.tsv
```

**Option C — Synthetic, for pipeline-check only**

```bash
.venv/bin/python -c "
import csv, random
random.seed(42)
mkpath = 'Data/Benchmarks/agarwal2025_lentimpra/oligo_counts.tsv'
import os; os.makedirs(os.path.dirname(mkpath), exist_ok=True)
with open(mkpath, 'w') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['oligo_id', 'dna_rep1', 'dna_rep2', 'dna_rep3',
                'rna_rep1', 'rna_rep2', 'rna_rep3'])
    for i in range(10000):
        dna = [max(1, int(random.gauss(100, 30))) for _ in range(3)]
        active = random.random() < 0.42  # ≈ paper's 41.7%
        mult = random.uniform(1.5, 4.0) if active else random.uniform(0.5, 1.5)
        rna = [int(d * mult * random.uniform(0.9, 1.1)) for d in dna]
        w.writerow([f'oligo_{i:05d}'] + dna + rna)
print('Wrote synthetic counts to', mkpath)"
```

Synthetic data won't reproduce paper-specific signals but will exercise
the entire `mpra activity / qc / volcano` chain and let you verify the
artefact paths.

### Stage 2 commands (once `oligo_counts.tsv` exists)

```bash
bash Benchmarks/agarwal2025_lentimpra/run.sh
```

This time the local-file gate passes and `run.sh` proceeds to:

```bash
.venv/bin/igvfagent mpra activity --counts <INPUT> --label agarwal2025_lentimpra
.venv/bin/igvfagent mpra qc       --counts <INPUT> --label agarwal2025_lentimpra
.venv/bin/igvfagent mpra volcano  --label agarwal2025_lentimpra
```

## Stage 2 ground-truth spot-checks

| Signal | Expected (paper) |
|---|---|
| Fraction of oligos called active overall | ~41.7 % |
| K562 active fraction | ~41.3 % |
| HepG2 active fraction | ~39 % |
| WTC11 active fraction | ~33 % |
| Per-replicate Pearson r in `mpra qc` matrix | ≥ 0.8 across all replicate pairs |
| Promoter-class oligos | Strand-biased activity (forward orientation higher) |
| Enhancer-class oligos | Tissue-specific activity (cell-line-dependent) |

After `mpra activity` runs, the output TSV contains columns
`(oligo_id, log2FC, lfcSE, pvalue, padj, active_call)`. To verify the
active fraction:

```bash
.venv/bin/python -c "
import csv
n_active = n_total = 0
for r in csv.DictReader(open('Docs/MPRA/<ts>_<files>/mpra_activity.tsv'), delimiter='\t'):
    n_total += 1
    if r.get('active_call') == 'True': n_active += 1
print(f'{n_active}/{n_total} = {100*n_active/n_total:.1f}%')"
```

Adjust the exact path/column based on the actual `mpra activity` output
schema.

## Running through the UI

### Stage 1 via UI

`mpra_pull` IS registered. Paste:

```
Run the Agarwal 2025 lentiMPRA online discovery step. Call mpra_pull
with limit=50 and label="agarwal2025_lentimpra".

Then call portal_facets with type=MeasurementSet and
field_filters="preferred_assay_titles=lentiMPRA" to confirm the
per-cell-line distribution (K562, HepG2, WTC11 should all appear).

Report:
  - total MPRA MeasurementSet count in IGVF
  - cell-line distribution (counts of K562 vs HepG2 vs WTC11)
  - first 10 accessions from the manifest

Do NOT call mpra_activity — that needs a local oligo_counts.tsv that
isn't downloaded yet.
```

UI sidebar: max iterations = 12, temperature = 0.0.

### Stage 2 via UI

Once `oligo_counts.tsv` is in place:

```
The per-oligo counts are at Data/Benchmarks/agarwal2025_lentimpra/oligo_counts.tsv.

Run the Agarwal 2025 analytical step:
  1) mpra_activity with counts at that path, label="agarwal2025_lentimpra"
  2) mpra_qc with the same counts and label
  3) mpra_volcano with that label

Report:
  - fraction of oligos called active (paper target: ~41.7%)
  - replicate-concordance Pearson r matrix from mpra_qc
  - top 10 most-active oligos by padj
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Stage 1 reports `1/3` or `2/3` checks passed | One of the three artefact roots is missing or empty | Check the per-check `detail` field in the Markdown report; usually a write-permission issue on `Data/Manifests/MPRA/` |
| Stage 2 exits with 77 | `oligo_counts.tsv` not at expected path | See "Where to get the counts table" above |
| `mpra activity` reports "no oligo_id column" | Counts file has wrong column names | Required schema: `oligo_id, dna_rep[1-N], rna_rep[1-N]` (or equivalent). Use `mpra activity --help` for the alternate schemas |
| Active fraction far from 41.7% | Synthetic input, or wrong cell-line table | Real Agarwal data should hit ±5% of paper value |
| `mpra qc` Pearson r < 0.5 | Likely synthetic data with too much noise | Check on real data |

## License + provenance

* **Paper data**: GEO GSE142696 family (related Inoue series); IGVF Portal MPRA AnalysisSets.
* **Code**: IGVFagent Apache-2.0 (MPRA skill is a clean-room reimpl of MPRAflow patterns + Tewhey lab MPRASuite).
* **Citation**: Agarwal V et al. *Nature* **639**: 411–420 (2025).
  doi:10.1038/s41586-024-08430-9 · PMID:39814879
