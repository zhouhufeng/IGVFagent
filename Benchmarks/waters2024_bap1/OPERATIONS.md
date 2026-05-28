# OPERATIONS — Waters 2024 BAP1 SGE

> **URN-pending.** The scoreset URN this benchmark needs is not yet
> verified. You will get a clean exit 77 with instructions until you
> supply the right URN. Use `matreyek2018_pten_vampseq` as a smoke-test
> to confirm the surrounding plumbing works.

For shared prerequisites + the two execution paths + concordance scoring,
see `Benchmarks/OPERATIONS_GUIDE.md`. Paper context in this directory's
`README.md`.

## The URN problem

The BAP1 SGE experiment is registered on MaveDB at
`urn:mavedb:00000662-0` with 17 per-exon sub-experiments
(`urn:mavedb:00000662-a` for Exon 1 through `…-i` for Exon 11 etc.).
Each experiment has one or more **scoresets**, and the scoreset URN is
what `mavedb map-scoreset --urn …` actually needs.

The placeholder in `run.sh` is `urn:mavedb:00000662-a-1`; this is the
**conventional pattern** but isn't guaranteed to resolve. Verify before
running.

### How to find the real URNs

```bash
# List the per-exon experiments
curl -s 'https://api.mavedb.org/api/v1/experiments/' | \
    .venv/bin/python -c "
import json, sys
for exp in json.load(sys.stdin):
    title = (exp.get('title') or '')
    if 'BAP1' in title:
        print(exp.get('urn'), '-', title[:80])
"
```

For each experiment URN that comes back (e.g. `urn:mavedb:00000662-a`),
open the MaveDB web page

```
https://www.mavedb.org/#/experiments/urn:mavedb:00000662-a
```

and copy the scoreset URN(s) listed on the page (they have the form
`urn:mavedb:00000662-a-1`, `urn:mavedb:00000662-a-2`, …). Pre-flight
each one with:

```bash
curl -sI "https://api.mavedb.org/api/v1/score-sets/<URN>/scores" | head -1
```

A `HTTP/2 200` means the URN is good. A `HTTP/2 404` means the URN is wrong.

## Quick run — once URN is verified

The `run.sh` honors a `BAP1_URN` environment variable. Once you have a
verified URN:

```bash
BAP1_URN="urn:mavedb:00000662-a-1" \
    bash Benchmarks/waters2024_bap1/run.sh

.venv/bin/python Benchmarks/concordance.py --benchmark waters2024_bap1
```

If the URN is unverified, `run.sh` exits 77 cleanly (the suite-runner
treats 77 as "skipped, not failed") and prints these instructions.

## What `run.sh` does (once URN is set)

```bash
URN="${BAP1_URN:-urn:mavedb:00000662-a-1}"
.venv/bin/igvfagent mavedb map-scoreset \
    --urn "$URN" \
    --gene BAP1 \
    --label waters2024_bap1
```

Same pipeline as the PTEN smoke-test: MaveDB download → Ensembl
transcript resolution → codon-to-genomic mapping. BAP1's canonical
transcript is ENST00000460680 (ENSG00000163930) at chr3:52,401,000–52,443,000.

## Where artefacts will land

```
Docs/MaveDB/<ts>_waters2024_bap1/
├── BAP1_mapped.tsv
├── BAP1_mapped.vcf
├── summary.json
├── mapping_summary.png
└── showcase_report.md
```

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | `n_rows_in ∈ [15000, 25000]` (expect ≈ 18,108) | Full BAP1 SGE scoreset downloaded |
| 2 | `n_mapped_rows_out ∈ [5000, 80000]` (expect ≈ 17,000) | Codon-mapping yielded sensible output |
| 3 | `gene = "BAP1"` | Symbol resolved cleanly |
| 4 | `BAP1_mapped.vcf` exists, non-empty | VCF emission worked |
| 5 | `BAP1_mapped.tsv` has 5,000–80,000 rows | Plain row-count cross-check |

If you pulled only one exon's scoreset, expect ~1,000-2,000 rows, not
18,108 — adjust `expected.json` or pull all 17 exon scoresets in a
loop and concatenate.

## Ground-truth spot-checks

| Signal | Expected |
|---|---|
| Chromosome | All rows show `chr3` |
| Position range | Approximately `52,400,000-52,443,000` GRCh38 |
| Total SNVs perturbed | 18,108 across the full 17-exon scoreset (paper Table 1) |
| LOF variants | ~5,665 (paper) |
| GOF variants | ~531 (paper) |
| Most-depleted domain | UCH (exons 1–9, 15–17) |
| Pathogenicity concordance | ~99.8 % of ClinVar-pathogenic truth-set in LOF bucket |

These are paper-level signals — the spot-check requires you to slice
the resulting TSV by exon or position window and check the score
distribution. For example, `awk -F'\t' '$2 ~ /^chr3:5240[0-2]/ {print}'
BAP1_mapped.tsv` should give you Exon-1 rows.

## Running through the UI

Same as the PTEN smoke-test — `mavedb_map_scoreset` is registered.
Paste:

```
Run the Waters 2024 BAP1 SGE reproducibility benchmark. Call
mavedb_map_scoreset with urn="urn:mavedb:00000662-a-1" (or whatever
verified URN I supply via the BAP1_URN environment), gene="BAP1",
label="waters2024_bap1".

Then report:
  - the row count in BAP1_mapped.tsv
  - the chromosome (should be chr3)
  - the position range of mapped variants
  - the breakdown by mapping_type
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Exit 77 + "URN not yet verified" | The pre-flight curl check failed | Look up real URN on MaveDB; export `BAP1_URN=` and rerun |
| Score-set returned but only ~1,000 rows | You pulled one exon's scoreset, not the unified one | Loop over all 17 sub-URNs and concatenate, or pick the unified scoreset URN (if MaveDB has one) |
| Ensembl can't resolve BAP1 | Symbol-lookup race | Retry once; if persistent, force `--species human` explicitly |
| Position range looks wrong (e.g. chr2 instead of chr3) | Wrong transcript chosen | Inspect `summary.json["transcript_id"]` and `["chromosome"]`; if not chr3, the Ensembl call resolved to a paralog |

## License + provenance

* **Paper data** (MaveDB scoreset): MaveDB CC-BY 4.0.
* **Code**: IGVFagent Apache-2.0.
* **Citation**: Waters AJ et al. *Nat Genet* **56**: 1434–1445 (2024).
  doi:10.1038/s41588-024-01799-3 · PMID:38969833
