---
name: igvf-catalog-disease-genes
description: Map the genetic architecture of a disease — causal genes, the pathways that converge on them, interaction-network hubs, and therapeutic-target opportunities. Use when the user asks "what genes cause X?", "genetic basis of disease X", or wants the gene network underlying a condition.
argument-hint: <disease_name_or_id>
---

# Disease-gene network and pathway architecture

Five-step workflow that goes from a disease term to a structured gene-network report. Re-derived under Apache-2.0 from the workflow described in [IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp) (MIT) and retargeted at IGVFagent's `catalog` CLI.

`$ARGUMENTS` may be a disease ontology ID (`MONDO:0008199`, `Orphanet:586`, `EFO:0000305`, `DOID:8778`) or a free-text disease name (`cystic fibrosis`, `Alzheimer's disease`, `type 2 diabetes`).

## Workflow

### 1. Find the disease and its associated genes

**If a disease ID was provided:**

```bash
igvfagent catalog get-entity "$ARGUMENTS"
igvfagent catalog find-associations "$ARGUMENTS" --relationship genetic \
                                                    --limit 100 --verbose
```

**If a disease name was provided:**

```bash
igvfagent catalog get-entity "$ARGUMENTS"
```

If this returns a single ontology term, proceed as above. Otherwise — and the disease→gene direction is frequently sparse on free-text — use a known primary gene as the entry point. For "cystic fibrosis" that's `CFTR`; for "sickle cell" that's `HBB`; etc. Then:

```bash
igvfagent catalog find-associations "<known_gene>" --relationship genetic \
                                                      --limit 50 --verbose
```

Inspect the returned disease IDs and switch to the disease-centric query once you have a confident ID.

In all cases, end this step with a working list of the **disease's gene set** (often 5–50 genes; can be larger for highly polygenic traits).

### 2. Pathway analysis — fan out across the gene set (parallel)

For each disease gene (up to ~10 to keep the workflow tractable):

```bash
for g in $DISEASE_GENES; do
    igvfagent catalog find-associations "$g" --relationship functional \
                                                --limit 50 --verbose &
done
wait
```

This pulls Reactome pathway memberships per gene.

### 3. Interaction analysis — fan out again (parallel)

```bash
for g in $DISEASE_GENES; do
    igvfagent catalog find-associations "$g" --relationship physical \
                                                --limit 100 --verbose &
done
wait
```

This pulls genetic + physical PPI partners (BioGRID / IntAct).

### 4. Convergence analysis (the core of the report)

Three counts you build by hand from the steps above:

| Metric | Definition |
|---|---|
| **Pathway convergence** | Reactome pathways that contain **≥ 2** disease genes. These are the *shared biology* of the disease. |
| **Interaction hubs** | Non-disease genes that interact with **≥ 2** disease genes. Strong candidates for novel disease genes or therapeutic targets. |
| **Direct disease-gene interactions** | Disease genes that interact with each other. Often indicates the same protein complex / pathway, and is the strongest network signal. |

For each, store the count of "supporting disease genes" — that's the rank metric.

### 5. Compile the report

**Disease overview.** Display name, all cross-reference IDs (Orphanet / MONDO / EFO / DOID), total gene count returned.

**Disease-gene table.** One row per gene, ranked by association strength (the Orphanet `association_type` ladder is the strongest single signal):

```
Disease-causing germline ≫ somatic ≫ Major susceptibility ≫
Modifying ≫ Candidate ≫ Role in phenotype ≫ Biomarker
```

Columns: gene, association type, source, PMID-count, notes.

**Pathway convergence.** Table of pathways with ≥ 2 disease genes, ranked by `n_disease_genes / pathway_size`. Show pathway name, IDs, supporting genes, pathway size.

**Interaction network.**
- Direct edges between disease genes (the tightest core).
- Top 10 hub genes (non-disease, ≥ 2 disease-gene partners), with their partners and the interaction types.

| Hub gene | Disease partners | Edge types |
|---|---|---|

**Therapeutic insights.** Surface three angles:
- Synthetic-lethal partners of disease genes (cancer-relevant).
- Druggable hubs (cross-reference the hub list against the catalog's `proteins_drugs` edge if available, otherwise note as a follow-up).
- Pathway-based intervention points (pathway hits with high disease-gene-coverage).

**Architecture summary, 3–5 sentences.**
- Monogenic vs oligogenic vs polygenic? (Count of "Disease-causing germline" rows.)
- Does the genetic signal converge on a small number of pathways or scatter?
- Hub genes that look therapeutically promising.
- Caveats — gene set incompleteness, ascertainment bias, etc.

## Interpretation cheat-sheet

- **Pathway convergence beats per-gene novelty.** Two disease genes in the same Reactome pathway tell you more about the disease's biology than ten disease genes in unrelated pathways.
- **Hub genes** with multiple disease-gene partners and no disease-gene status themselves are the strongest novel-target candidates — see PARP / BRCA-mutant cancers for the canonical example.
- **Monogenic** = 1–3 strong "Disease-causing" entries with very high effect size (Mendelian; e.g. CF/CFTR). **Polygenic** = many "Susceptibility" entries with modest individual effects (T2D, schizophrenia). The shape of the disease-gene table tells you which one you're looking at.
- The Orphanet `association_type` ladder is the single best filter for casual / causal claim strength.

## Pairs well with

- `igvf-catalog-gene-dossier` — build a deep dossier for any nominated hub or strong causal gene.
- `igvfagent network steiner --seeds <disease_gene_list>` — produces a context-specific subnetwork at MILP optimum, complementing this skill's pathway-centric view.
- `igvfagent enrich pathways --genes <disease_gene_list>` — statistical enrichment confirmation for the convergence-pathways list.
