---
name: igvf-catalog-dissect-locus
description: Go from a GWAS hit to candidate causal gene(s) by integrating LD block, eQTL evidence, and regulatory elements in the locus. Use when the user asks "which gene does this GWAS variant affect?", "fine-map this locus", or "is the nearest gene actually causal?".
argument-hint: <variant_or_region>
---

# GWAS-to-causal-gene locus dissection

Seven-step fine-mapping workflow built from the IGVF Catalog. Re-derived under Apache-2.0 from the workflow described in [IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp) (MIT) and retargeted at IGVFagent's `catalog` CLI.

`$ARGUMENTS` is either a lead variant (`rs429358`, SPDI, HGVS, CA-ID) or a region (`chr19:44.9M-44.92M`). If a region is given, skip steps 1-3 and start at step 4 with the region as the LD-block proxy.

## Workflow

### 1. Resolve the lead variant

```bash
igvfagent catalog resolve-id "$ARGUMENTS"
```

Capture the canonical SPDI + GRCh38 coordinates. These are the anchor for everything downstream.

### 2. Confirm the GWAS signal(s)

```bash
igvfagent catalog find-associations "$ARGUMENTS" --relationship genetic \
                                                    --filters "p_value=lte:5e-8" \
                                                    --limit 25 --verbose
```

Record the phenotype(s), p-value(s), effect size(s), and source studies. If nothing comes back at genome-wide significance, drop to `p_value=lte:1e-5` (suggestive) for an exploratory pass.

### 3. Build the credible set via LD

```bash
igvfagent catalog find-ld "$ARGUMENTS" --r2-threshold 0.5 --ancestry EUR --limit 200
```

Use the moderate threshold (0.5) so the credible set is inclusive. Note which proxies sit at r² ≥ 0.8 (strong, almost certainly tagging the same haplotype) vs 0.5 ≤ r² < 0.8 (moderate, possibly distinct).

If the locus is fine-mapped poorly in EUR (e.g. one big haplotype block > 500 kb), run

```bash
igvfagent catalog find-ld "$ARGUMENTS" --r2-threshold 0.5 --ancestry AFR --limit 200
```

— AFR usually has shorter blocks and narrows the candidate set. See also `igvf-catalog-ld-compare` for a full cross-ancestry view.

### 4. Survey the locus

Define a locus span: either the outermost LD-proxy coordinates from step 3, or lead ± 500 kb if the LD block is small.

```bash
igvfagent catalog search-region "<chrom>:<start>-<end>" \
                                  --include "genes,variants,genomic-elements" \
                                  --limit 100
```

You now have every gene, every catalogued variant, and every ENCODE SCREEN cCRE in the locus.

### 5. eQTL evidence for the lead + top proxies

Pick the lead plus the top 5 proxies by r². Run `find-associations` for each:

```bash
for v in "$LEAD" "$PROXY1" "$PROXY2" "$PROXY3" "$PROXY4" "$PROXY5"; do
    igvfagent catalog find-associations "$v" --relationship regulatory \
                                                --limit 50 --verbose &
done
wait
```

This tells you which gene(s) each variant regulates, in which tissues, and with what effect size.

### 6. Rank candidate causal genes

Bucket every gene in the locus into one of three tiers:

| Tier | Criterion |
|---|---|
| **1 — Strong** | Has an eQTL variant at r² ≥ 0.8 with the lead, in a disease-relevant tissue. |
| **2 — Moderate** | Has an eQTL at 0.5 ≤ r² < 0.8, **or** overlaps a regulatory element that harbours an LD variant. |
| **3 — Positional only** | In the locus but no functional link to the GWAS signal. |

### 7. Compile the report

Sections, in order:

**Lead variant & GWAS signal.** Identifier, coordinates, trait(s), p-value(s), effect size(s).

**LD block & credible set.** Total variants at r²≥0.5 and r²≥0.8. Block span (kb). Ancestry used. If multi-ancestry, note where blocks differ.

**Candidate causal-gene table:**

| Gene | Tier | eQTL evidence | Best LD-r² | Tissue(s) | Distance to lead |
|---|---|---|---|---|---|

For Tier 1 / 2 genes, attach which variant is the eQTL, its r² with the lead, the effect direction/size, and the tissues.

**Regulatory-element map.** ENCODE SCREEN cCREs in the locus and which overlap LD variants (these are mechanistic candidates for *how* the GWAS variant affects the gene).

**Interpretation, 3–5 sentences.**
- Which gene(s) you'd nominate, and why.
- Whether the nearest gene wins or loses against the best functional candidate.
- Likely mechanism (which enhancer / promoter / coding change).
- Caveats — block size, tissue gaps, ancestry, etc.

## Interpretation cheat-sheet

- The nearest gene is the causal one **only ~50 %** of the time at GWAS loci — always weight functional evidence over proximity.
- Tissue context dominates: an eQTL in the disease-relevant tissue (pancreatic islets for T2D, microglia for late-onset AD, liver for LDL) beats a stronger eQTL in an irrelevant tissue.
- "Two strong candidates" is a legitimate and frequent answer — report both rather than forcing one.
- Convergence is the strongest evidence: same variant being a GWAS hit *and* an eQTL *and* overlapping a tissue-matched enhancer ⇒ that's your causal regulatory variant.
- Wide LD blocks (> 500 kb) make EUR fine-mapping ambiguous; bring in AFR.

## Pairs well with

- `igvf-catalog-variant-report` — generate the per-variant deep dive for the lead before dissection.
- `igvf-catalog-ld-compare` — when EUR can't fine-map, multi-ancestry LD is the next move.
- `igvf-catalog-regulatory-landscape` — for the gene(s) you nominate, map their regulatory architecture in full.
