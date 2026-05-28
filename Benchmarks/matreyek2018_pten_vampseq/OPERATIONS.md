# OPERATIONS — Matreyek 2018 PTEN VAMP-seq (★ suite smoke-test)

> **Use this first.** It is the only MaveDB benchmark in the suite with
> a URN that's verified to resolve on the live API, and it requires no
> local data downloads. If it passes, you've proven that the MaveDB +
> Ensembl REST + codon-mapping pipeline is healthy on your machine.

For shared prerequisites + the two execution paths + concordance scoring,
see `Benchmarks/OPERATIONS_GUIDE.md`. Paper context + ground-truth
checks are in this directory's `README.md`.

## Quick run — 2 commands, ~60 seconds

```bash
# 1) Execute the IGVFagent pipeline
bash Benchmarks/matreyek2018_pten_vampseq/run.sh

# 2) Score the result against expected.json
.venv/bin/python Benchmarks/concordance.py --benchmark matreyek2018_pten_vampseq
```

Expected console output:

```
=== Concordance — <ts> ===
  [          ok] matreyek2018_pten_vampseq            4/4
```

## What `run.sh` actually does

```bash
.venv/bin/igvfagent mavedb map-scoreset \
    --urn urn:mavedb:00000013-a-1 \
    --gene PTEN \
    --label matreyek2018_pten_vampseq
```

Behind the scenes IGVFagent:

1. Downloads `urn:mavedb:00000013-a-1` from `api.mavedb.org` (~5 MB CSV,
   4,408 rows).
2. Resolves the canonical Ensembl transcript for PTEN via the REST API
   (`/xrefs/symbol/human/PTEN` → ENSG00000171862 → ENST00000371953).
3. For each variant in the CSV, queries
   `/map/translation/{ENSP_id}/{aa}..{aa}` to get the codon's genomic
   coordinates.
4. Reverse-complements where the strand is `-1`, validates against the
   genomic 3-nt substring, enumerates single-nt changes that produce
   the alt amino acid.
5. Emits per-variant TSV + VCF-4.2 + summary JSON + composite SVG/PNG.

## Where artefacts land

```
Docs/MaveDB/<ts>_matreyek2018_pten_vampseq/
├── PTEN_mapped.tsv          ~700 KB, 4,656 rows
├── PTEN_mapped.vcf          ~175 KB, 4,656 records
├── summary.json             gene metadata + type counts
├── mapping_summary.png      coverage along protein + outcome bar
├── mapping_summary.svg
└── showcase_report.md       narrative
```

The Ensembl REST cache (so reruns take ~5 s instead of ~60 s) lives at
`Data/Cache/Ensembl/*.json` — safe to delete; it just makes the next
run slower.

## Concordance interpretation

The 4 checks in `expected.json`:

| # | Check | What it verifies |
|---|---|---|
| 1 | `n_rows_in ∈ [1000, 20000]` (expect ≈ 8,000) | MaveDB download + CSV parse worked |
| 2 | `gene = "PTEN"` | Ensembl symbol resolution worked |
| 3 | `PTEN_mapped.tsv` has 1,000–50,000 rows | Codon mapping produced output |
| 4 | `PTEN_mapped.vcf` exists and is non-empty | VCF emission worked |

A `1/4` or `2/4` partial result is the standard signal of network /
Ensembl REST instability on first run — retry once.

## Ground-truth spot-checks (eyeball-level)

After scoring, open `Docs/MaveDB/<ts>_matreyek2018_pten_vampseq/PTEN_mapped.tsv`
and verify:

| Signal | Expected |
|---|---|
| Chromosome | All rows show `chr10` (PTEN locus) |
| Position range | Approximately `87,864,000-87,966,000` GRCh38 |
| `mapping_type` value distribution | `missense` ≈ 1,500–1,800 rows ; `synonymous` < 250 ; `nonsense` < 100 |
| Transcript id | `ENST00000371953` (canonical PTEN) for all rows |

Open `mapping_summary.png` — the top panel should show a roughly even
spread of variants from residue 1 to ~400 (PTEN is 403 aa); the bottom
panel breaks down missense vs synonymous vs nonsense vs error.

## Running through the UI

This benchmark's CLI uses `mavedb map-scoreset` which **is** registered
as the LLM tool `mavedb_map_scoreset`. Paste this into the chat box:

```
Run the Matreyek 2018 PTEN VAMP-seq reproducibility benchmark:
call mavedb_map_scoreset with urn="urn:mavedb:00000013-a-1",
gene="PTEN", label="matreyek2018_pten_vampseq".

Then report:
  - the row count in the resulting PTEN_mapped.tsv
  - the breakdown by mapping_type
  - the chromosome (should be chr10)
  - the canonical Ensembl transcript that was resolved
```

UI sidebar: max iterations = 5 (this is a 1-tool-call benchmark).

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `FileNotFoundError: ... urn-mavedb-00000013-a-1.csv` | MaveDB download failed silently (404 or network) | Probe `curl -s https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000013-a-1/scores | head -1` — should return CSV header `accession,hgvs_nt,...` |
| `Mapping failed: No transcript found for PTEN` | Ensembl REST timeout | Wait 30s, rerun — Ensembl REST is sometimes overloaded |
| Per-variant TSV has < 1,000 rows | Most variants hit `no_single_nt_change` (alt codons need ≥ 2 nt edits) | Look at `summary.json["type_counts"]["no_single_nt_change"]` — typically ~60% of rows for SGE-style scoresets |
| `n_rows_in` reports 0 | The Ensembl cache contains a stale 404 response | `rm Data/Cache/Ensembl/*.json` then rerun |

## License + provenance

* **Paper data** (MaveDB scoreset): MaveDB CC-BY 4.0. We fetch / link only.
* **Code**: IGVFagent Apache-2.0 over upstream MIT (`ave-dcd/dcd_mapping`).
* **Citation**: Matreyek KA et al. *Nat Genet* **50**: 874–882 (2018).
  doi:10.1038/s41588-018-0122-z · PMID:29785012
