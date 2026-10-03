---
name: igvf-catalog-regulatory-landscape
description: Map how a gene is regulated — cis-regulatory elements in the locus, eQTL/sQTL variants by tissue, tissue-specificity patterns, and convergence between eQTLs and GWAS hits. Use when the user asks "how is X regulated?", "what controls expression of X?", or "what's the regulatory architecture of this locus?".
argument-hint: <gene> [tissue]
---

# Regulatory-architecture map for a gene

Six-section workflow that enumerates the regulatory landscape around a target gene. Re-derived under Apache-2.0 from the workflow described in [IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp) (MIT) and retargeted at IGVFagent's `catalog` CLI.

Parse `$ARGUMENTS` as `<gene> [tissue]`. The first token is the gene; any second token is an optional tissue / cell-type filter (e.g. `liver`, `pancreatic islet`, `cortex`).

## Workflow

### 1. Anchor the gene

```bash
igvfagent catalog get-entity "<gene>"
```

Capture chr / start / end and the canonical ENSG. The regulatory window is the gene span plus 500 kb on each side — far enough to include most distal enhancers.

### 2. Enumerate cis-regulatory elements

```bash
igvfagent catalog search-region "<chrom>:<start-500kb>-<end+500kb>" \
                                  --include "genomic-elements" \
                                  --limit 200
```

Classify each element by its `source_annotation` value into one of six SCREEN classes:

| Class | Meaning |
|---|---|
| PLS | Promoter-like (active TSS, high H3K4me3) |
| pELS | Proximal enhancer (≤ 2 kb from a TSS) |
| dELS | Distal enhancer (> 2 kb from any TSS) |
| CTCF | Insulator |
| CA | Open chromatin, otherwise unclassified |
| TF | TF binding-site only |

### 3. Pull eQTL / sQTL variants

```bash
igvfagent catalog find-associations "<gene>" --relationship regulatory \
                                                --limit 100 --verbose
```

If a `[tissue]` was supplied, add `--filters "biological_context=<tissue>"` to restrict to that tissue's QTL studies.

### 4. Profile tissue-specificity

Group eQTL hits by their `biological_context` (tissue / cell type). Count per tissue and rank by `-log10(P)` of the strongest hit. Then bin each variant:

- **Ubiquitous**: significant in ≥ 5 tissues — constitutive regulation.
- **Tissue-specific**: significant in ≤ 2 tissues — biologically more interesting, and a strong candidate for the mechanism behind tissue-restricted phenotypes.

If `[tissue]` was specified, compare its eQTL set to the gene's overall set and highlight what is gained / lost.

### 5. Regulatory–disease convergence

For the top 5–10 eQTL variants (highest `-log10(P)`), run in parallel:

```bash
for v in "$VAR1" "$VAR2" ... "$VARN"; do
    igvfagent catalog find-associations "$v" --relationship genetic \
                                                --filters "p_value=lte:5e-8" \
                                                --limit 25 --verbose &
done
wait
```

Any variant that is **both** an eQTL for this gene **and** a genome-wide-significant GWAS hit is a high-value finding: it directly nominates the gene as causal for the implicated trait in the implicated tissue.

### 6. Compile the report

**Gene overview.** Symbol, full name, coordinates, biotype. Regulatory window analysed.

**cis-regulatory-element census.** Total cCRE count. Breakdown by class. Note interesting patterns:
- Clusters of dELS (a super-enhancer signal).
- CTCF flanking the TAD that contains the gene.
- Presence/absence of a PLS exactly at the canonical TSS.

**eQTL / sQTL landscape.** Total regulatory variants. Top variants table: ID, p-value (scientific), β, tissue, eQTL vs sQTL.

**Tissue-specificity profile.** Tissues ranked by eQTL count. Tissue-specific (≤ 2 tissues) variants highlighted explicitly. If `[tissue]` is specified, show how it compares to "everything else".

**Regulatory–disease convergence.** Table of variants that appear in both the eQTL set and the GWAS hits, with: eQTL effect (gene + tissue + direction) and GWAS hit (trait + p-value). One-sentence per row interpreting the convergence.

**Synthesis, 3–5 sentences.** Regulatory complexity (simple vs distal enhancer-driven), key tissues, disease-relevant variants, overall architecture.

## Interpretation cheat-sheet

- A predominance of **dELS** elements ⇒ complex distal regulation likely through chromatin looping; **mostly PLS** ⇒ simpler promoter-driven.
- **Ubiquitous eQTLs** = constitutive housekeeping regulation; **tissue-specific eQTLs** = where the biology and disease relevance usually live.
- eQTL ⊕ GWAS convergence in a relevant tissue is the cleanest evidence that *this variant regulates this gene which is causal for this disease*.
- **sQTLs** (splicing QTLs) can be disease-relevant even when total expression is unchanged — don't ignore them when total-expression eQTLs are weak.
- A dELS sitting > 100 kb from the TSS can still regulate the gene; chromatin contact data (Hi-C / micro-C) is the next step to confirm.

## Pairs well with

- `igvf-catalog-dissect-locus` — when the convergence table flags a GWAS hit, dissect that locus to confirm.
- `igvf-catalog-gene-dossier` — fold this regulatory map into the wider dossier on the gene.
- `igvfagent enhancer ...` / `igvfagent ccre ...` — IGVFagent's enhancer-gene linkage + cCRE annotation skills go deeper on the cCRE side.
