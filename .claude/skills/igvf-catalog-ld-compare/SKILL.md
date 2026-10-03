---
name: igvf-catalog-ld-compare
description: Compare LD structure across the five 1000 Genomes superpopulations (EUR / AFR / EAS / SAS / AMR) for a single variant. Use when the user asks about cross-population LD, PRS portability, ancestry-specific tag-SNP issues, or wants to use multi-ancestry LD to narrow a fine-mapping candidate set.
argument-hint: <variant_id>
---

# Cross-ancestry LD comparison

Four-step workflow comparing a variant's LD landscape across the five 1000 Genomes superpopulations. Re-derived under Apache-2.0 from the workflow described in [IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp) (MIT) and retargeted at IGVFagent's `catalog` CLI.

`$ARGUMENTS` is a single variant — rsID, SPDI, HGVS, or CA-ID.

## Workflow

### 1. Resolve the variant

```bash
igvfagent catalog resolve-id "$ARGUMENTS"
```

Capture the canonical SPDI plus coordinates.

### 2. Query all five superpopulations in parallel

Fire these five calls concurrently — they are completely independent.

```bash
for anc in EUR AFR EAS SAS AMR; do
    igvfagent catalog find-ld "$ARGUMENTS" --r2-threshold 0.5 \
                                              --ancestry "$anc" \
                                              --limit 200 \
                                              --verbose &
done
wait
```

The skill's bucket summary (`strong / moderate / weak / negligible`) is emitted per call — you can inspect the per-ancestry JSON in `Docs/CatalogQuery/<ts>_ld_<variant>/ld.json`.

### 3. Cross-reference the proxy sets

Build three lists from the five per-ancestry proxy sets:

| Bucket | Definition |
|---|---|
| **Universal proxies** | r² ≥ 0.8 in **≥ 4** of the 5 populations. These are ancient haplotypes (often predating the out-of-Africa expansion) and are the most reliable cross-ancestry tag SNPs. |
| **Population-specific proxies** | Strong LD (r² ≥ 0.8) in **1 or 2** populations only. Useful for ancestry-aware fine-mapping and for diagnosing PRS-portability failures. |
| **AFR-narrowing proxies** | Strong in EUR/EAS/SAS but weak in AFR — diagnostic of an extended European haplotype that AFR's shorter blocks can dissect. |

### 4. Report

**Query variant.** Identifier, coordinates, reference / alternate alleles.

**Per-population summary table:**

| Population | r² > 0.5 | r² > 0.8 | Block span (kb) | Best proxy (r²) |
|---|---|---|---|---|

**Universal proxies.** Variants in strong LD across ≥ 4 populations. These are the variants you'd choose as a cross-ancestry tag for genotyping array design.

**Population-specific proxies.** Variants in strong LD in only 1 or 2 populations. List them with their per-ancestry r² values.

**Cross-population r² matrix** (top 20 proxies, sorted by mean r²):

| Variant | EUR | AFR | EAS | SAS | AMR |
|---|---|---|---|---|---|

**Interpretation, 3–5 sentences.**
- Block-size comparison: who has the longest, who has the shortest? AFR almost always has the shortest.
- PRS-portability implications: if the index is a strong EUR tag with no AFR proxies, expect EUR-trained PRS to perform poorly in AFR for any phenotype this variant tags.
- Fine-mapping value: AFR's shorter blocks can cut the EUR credible set down by an order of magnitude — call out the gain numerically (e.g. "EUR credible set has 47 variants; AFR has 6").
- Universal vs population-specific signal: the more population-specific the LD, the more confidence you have that the variant itself is causal in some populations.

## Interpretation cheat-sheet

- **AFR LD blocks are systematically shorter** because the African population is older and has had more recombination time. EUR / EAS / SAS blocks tend to be longer; AMR is highly variable due to admixture.
- A **European tag SNP** can be **completely uninformative in AFR** — this is the central driver of PRS portability failure across ancestries.
- **r² thresholds:**

| r² | Bucket |
|---|---|
| ≥ 0.8 | Strong — tagging the same haplotype |
| 0.5 – 0.8 | Moderate |
| 0.2 – 0.5 | Weak |
| < 0.2 | Negligible |

- **AFR shorter blocks = a fine-mapping advantage:** if EUR can't resolve a credible set down to a handful of variants, AFR sometimes can — provided the variant is segregating at non-trivial frequency.
- Universal proxies (strong in ≥ 4 populations) generally predate human migration out of Africa and are the safest cross-ancestry tags.

## Pairs well with

- `igvf-catalog-variant-report` — generate the per-variant deep dive first, then use this skill if the LD section showed many proxies.
- `igvf-catalog-dissect-locus` — when EUR can't fine-map a GWAS locus, fold AFR results from this skill into the credible-set step.
- `igvfagent ccre ...` — once a population-specific proxy narrows the candidate region, annotate the cCREs in it.
