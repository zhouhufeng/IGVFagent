---
name: igvf-catalog-variant-report
description: Build a complete interpretation report for a single variant — disease/phenotype associations, eQTL/sQTL regulatory effects, pharmacogenomics, coding-impact scores, and LD context. Use when the user asks "what does rs… do?", "interpret this variant", or "is this variant causal vs a tag SNP?".
argument-hint: <variant_id>
---

# Variant interpretation report

A six-section structured report built from the IGVF Catalog (Knowledge Graph). Re-derived under Apache-2.0 from the workflow described in [IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp) (MIT) and retargeted at IGVFagent's `catalog` CLI.

`$ARGUMENTS` may be an rsID (`rs429358`), SPDI (`NC_000019.10:44908683:T:C`), HGVS, or ClinGen CA-ID. IGVFagent's `get-entity` autodetects the form.

## Workflow

### 1. Normalise identifiers

```bash
igvfagent catalog resolve-id "$ARGUMENTS"
```

This produces the canonical SPDI plus rsID / HGVS / CA cross-references and the GRCh38 coordinates needed downstream.

### 2. Gather evidence in parallel

Run these four `find-associations` calls plus an LD-proxy pull. They are independent — fire them concurrently if your runtime supports it.

```bash
# Disease / phenotype associations (GWAS + ClinGen + Orphanet)
igvfagent catalog find-associations "$ARGUMENTS" --relationship genetic \
                                                    --limit 50 --verbose

# eQTL / sQTL effects across tissues
igvfagent catalog find-associations "$ARGUMENTS" --relationship regulatory \
                                                    --limit 50 --verbose

# Pharmacogenomics (PharmGKB)
igvfagent catalog find-associations "$ARGUMENTS" --relationship pharmacological \
                                                    --limit 25 --verbose

# Protein / coding impact (dbNSFP — SIFT / PolyPhen / CADD / REVEL)
igvfagent catalog find-associations "$ARGUMENTS" --relationship coding \
                                                    --limit 25 --verbose

# LD proxies at strong threshold, EUR by default
igvfagent catalog find-ld "$ARGUMENTS" --r2-threshold 0.8 --ancestry EUR --limit 100
```

For GWAS-only narrowing, add `--filters "p_value=lte:5e-8"` — IGVFagent auto-translates that into the catalog's native `log10pvalue=gte:7.301` clause.

### 3. Compile the report

Render six sections, in this order:

**Identity & coordinates.** All identifiers, GRCh38 chr/pos/ref/alt.

**Disease & phenotype associations.** Ranked by p-value. For each row: trait, source, p-value (scientific notation), effect size (OR or β), flagging `p < 5e-8` as genome-wide significant.

**Regulatory effects.** Group by target gene, then by tissue. Report effect direction (`+` / `-`) and tissue specificity (ubiquitous vs ≤2 tissues).

**Pharmacogenomics.** Group by drug. Annotation level (1A / 1B / 2A / 2B / 3 / 4) + phenotype category (Efficacy / Toxicity / Metabolism-PK / Dosage) + biogeographical group. Flag Level 1A/1B explicitly.

**Coding impact** (only if the variant has any `coding` hits). Amino-acid change + the four standard scores with thresholds:

| Score | Damaging threshold |
|---|---|
| SIFT | < 0.05 |
| PolyPhen | > 0.85 |
| CADD | > 20 |
| REVEL | > 0.5 |

**LD structure.** Number of r²>0.8 proxies (EUR). Call out any proxy that itself appears in the disease, regulatory, or coding sections above — that is evidence the index variant is a tag rather than the causal one.

### 4. Two-to-four-sentence synthesis

End with a concise interpretation:

- Is this likely causal or a tag SNP for the proxies?
- Is the mechanism regulatory or coding?
- Which gene(s) are the strongest target?
- Any actionable clinical or pharmacogenomic context?

## Interpretation cheat-sheet

- `p < 5e-8` = genome-wide significant; `p < 1e-5` = suggestive only.
- GWAS hit + eQTL in a disease-relevant tissue ⇒ strong evidence for a regulatory mechanism through the eQTL target gene.
- Few r²>0.8 proxies ⇒ the index variant itself is more likely causal (easier fine-mapping).
- Many proxies ⇒ tag-vs-causal is ambiguous; consider cross-population LD (see `igvf-catalog-ld-compare`).
- Most GWAS hits are non-coding — a strong eQTL with no coding impact is the typical regulatory pattern.
- PharmGKB level ladder: 1A (clinical guideline + strong evidence) > 1B > 2A > 2B > 3 > 4 (case report only).

## Pairs well with

- `igvf-catalog-ld-compare` — when the LD-section has many proxies, run this next to compare ancestries.
- `igvf-catalog-dissect-locus` — when the report points at a regulatory mechanism, dissect the locus to nominate the causal gene.
- `igvf-catalog-gene-dossier` — once a target gene is nominated, build its full dossier.
