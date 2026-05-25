---
name: igvf-catalog-gene-dossier
description: Compile a 360-degree dossier on a single gene — coordinates and isoforms, disease links, regulatory landscape, physical interactions, pathway memberships. Use when the user asks "tell me about gene X", "build a dossier for BRCA1", or wants a comprehensive overview before committing to a target.
argument-hint: <gene_symbol_or_id>
---

# Comprehensive gene dossier

Six-section dossier built from the IGVF Catalog. Re-derived under Apache-2.0 from the workflow described in [IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp) (MIT) and retargeted at IGVFagent's `catalog` CLI.

`$ARGUMENTS` may be a gene symbol (`APOE`), ENSG (`ENSG00000130203`), HGNC (`HGNC:613`), or Entrez (`ENTREZ:348`).

## Workflow

### 1. Resolve identity + coordinates

```bash
igvfagent catalog get-entity "$ARGUMENTS"
igvfagent catalog resolve-id "$ARGUMENTS"
```

`get-entity` returns the node JSON (location, biotype, synonyms). `resolve-id` projects the full cross-reference set — symbol ↔ ENSG ↔ HGNC ↔ Entrez ↔ synonyms — which you'll cite in the summary.

### 2. Walk every relevant edge (parallel)

These five calls cover the disease, regulatory, physical, functional, and isoform axes of the dossier. They are independent — fire them concurrently.

```bash
# Disease associations (Orphanet + ClinGen)
igvfagent catalog find-associations "$ARGUMENTS" --relationship genetic \
                                                    --limit 50 --verbose

# eQTL / sQTL variants regulating this gene
igvfagent catalog find-associations "$ARGUMENTS" --relationship regulatory \
                                                    --limit 100 --verbose

# Genetic + physical protein-protein interactions (BioGRID, IntAct)
igvfagent catalog find-associations "$ARGUMENTS" --relationship physical \
                                                    --limit 100 --verbose

# Reactome pathway memberships
igvfagent catalog find-associations "$ARGUMENTS" --relationship functional \
                                                    --limit 50 --verbose

# GENCODE transcript isoforms + protein products
igvfagent catalog find-associations "$ARGUMENTS" --relationship transcription \
                                                    --limit 50 --verbose
```

### 3. Census of regulatory elements overlapping the locus

Use the coordinates from step 1 — extend by 5 kb on each side if desired — and run:

```bash
igvfagent catalog search-region "<chr>:<start-5kb>-<end+5kb>" \
                                  --include "genomic-elements" \
                                  --limit 100
```

This surfaces ENCODE SCREEN cCREs (promoters, enhancers, CTCF sites, etc.) that physically overlap the gene body and flanks — the "regulatory landscape" section of the dossier.

### 4. Compile the dossier

Render in this order:

**Gene summary.** Symbol, full name, biotype, every cross-reference, GRCh38 coordinates (1-based) and span.

**Disease associations.** Group by the Orphanet `association_type` ranking:

```
Disease-causing germline > somatic > Major susceptibility >
Modifying > Candidate > Role in phenotype > Biomarker
```

For each row: disease term ID + name, association type, source, PMID if any.

**Regulatory landscape.** Count of cCREs in the locus, broken down by `source_annotation`:

| Class | Meaning |
|---|---|
| PLS | Promoter-like; high H3K4me3 near TSS |
| pELS | Proximal enhancer (< 2 kb from a TSS) |
| dELS | Distal enhancer (> 2 kb from any TSS) |
| CTCF | Insulator |
| CA | Open chromatin, unclassified |
| TF | TF-binding only |

Then the top 5–10 eQTL/sQTL variants ranked by `-log10(P)`, with target gene, tissue (`biological_context`), and effect size.

**Interaction network.** Genetic interactions (synthetic lethal partners first — therapeutic relevance). Then PPI partners, by source (BioGRID / IntAct). Total counts. If you have time, render the PPI subnetwork via the network-viz skill (`igvfagent kg gene X --depth 1 --viz ...`).

**Pathway context.** Reactome pathways grouped by top-level category where possible (Signal Transduction, Immune System, Metabolism, etc.).

**Transcript isoforms.** Count protein-coding vs non-coding. List up to 10 isoforms with biotype + length. Note if any are canonical / MANE Select.

### 5. Two-to-four-sentence summary

- Primary biological role.
- Strongest disease links.
- Notable interaction partners (especially synthetic-lethal in cancer contexts).
- Whether the gene is well characterised or understudied (count of citations, presence/absence in clinical guidelines).

## Interpretation cheat-sheet

- **Synthetic lethal partners** are drug targets: if gene A is mutated in cancer, inhibiting synthetic-lethal partner B selectively kills the tumour (e.g. PARPi in BRCA1/2-mutant cancers).
- **Many dELS elements distal from the TSS** ⇒ complex enhancer-driven regulation; **mostly PLS** ⇒ simpler promoter-driven.
- **Many transcript isoforms** usually indicates tissue-specific regulation worth investigating with `igvf-catalog-regulatory-landscape`.
- The dossier is most useful when followed up with locus dissection for the strongest GWAS variants in the disease section.

## Pairs well with

- `igvf-catalog-regulatory-landscape` — go deeper on the eQTL/sQTL section.
- `igvf-catalog-disease-genes` — if the disease section is rich, expand it to the full disease-gene network.
- `igvfagent network steiner --seeds <gene>` — build a context-specific subnetwork around this gene from the IGVF KG.
