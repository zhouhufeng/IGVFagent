# Matreyek 2018 — PTEN VAMP-seq (the suite's verified smoke-test)

## Citation

Matreyek KA, ..., Fowler DM. *Nat Genet* **50**: 874–882 (2018).
DOI: [10.1038/s41588-018-0122-z](https://doi.org/10.1038/s41588-018-0122-z) · PMID: 29785012

## Data

* **MaveDB**: `urn:mavedb:00000013-a-1` — the canonical PTEN VAMP-seq scoreset. **URN verified live**.
* **GitHub**: [FowlerLab/VAMPseq](https://github.com/FowlerLab/VAMPseq)

## Why this is the suite's smoke-test

Every other MaveDB benchmark in this suite (`waters2024_bap1`, `buckley2024_vhl`) needs the user to verify the MaveDB URN before it can run end-to-end. The PTEN URN above is **verified working today** against the live MaveDB API, and IGVFagent's `mavedb showcase` defaults to it.

Use this benchmark to:

* prove your IGVFagent install is healthy
* sanity-check that the MaveDB downloader, Ensembl REST resolver, and codon-mapping pipeline all work end-to-end
* generate a reference `summary.json` you can use as the schema template when looking up real URNs for the other MaveDB benchmarks

## Headline workflow (paper)

1. VAMP-seq variant-abundance screen across all PTEN missense variants.
2. ~8,000 single-AA variants scored across 8 experimental replicates.
3. Per-residue abundance distribution; the catalytic phosphatase + C2 lipid-binding domains are most depleted.

## What IGVFagent reproduces

```bash
igvfagent mavedb map-scoreset \
    --urn urn:mavedb:00000013-a-1 \
    --gene PTEN \
    --label matreyek2018_pten_vampseq
```

* Downloads `urn:mavedb:00000013-a-1` (~5 MB CSV)
* Parses ~8 K per-variant abundance scores
* Resolves canonical PTEN Ensembl transcript ENST00000371953 (ENSG00000171862)
* Maps each variant to chr10:87864000-87966000 GRCh38 coordinates
* Emits TSV + VCF + JSON summary

## Ground truth to spot-check

* ~8,000 input rows (the per-variant CSV)
* PTEN canonical transcript resolved
* Chromosome 10 dominance in the mapped TSV
* TSV + VCF artefacts both written and non-empty
