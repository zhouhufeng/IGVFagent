# Skill: IGVF Catalog (Knowledge Graph) canonical-query layer

A faceted, ID-aware, edge-relationship-aware client for the IGVF
Catalog at ``api.catalogkg.igvf.org`` that matches the canonical query
patterns the IGVF Data Coordinating Center publishes in
[IGVF-DACC/igvf-catalog-mcp](https://github.com/IGVF-DACC/igvf-catalog-mcp)
(MIT, Cherry Lab / IGVF DACC, 2026).

Where the sibling ``portal`` skill covers the **Portal** item catalog
(data.igvf.org), this skill covers the **Knowledge Graph** — the
federated nodes-and-edges layer over genes / variants / proteins /
transcripts / ontology terms / drugs / complexes / studies / pathways
/ genomic elements, and all 30+ edges between them.

Clean-room reimplementation under Apache-2.0; stdlib-only (`urllib`).

## Commands

```bash
# Universal entity lookup — auto-detects type from the ID
igvfagent catalog get-entity APOE              # gene symbol
igvfagent catalog get-entity rs429358          # variant rsID
igvfagent catalog get-entity ENSG00000130203   # ENSG
igvfagent catalog get-entity P02649            # UniProt
igvfagent catalog get-entity MONDO:0004975     # disease ontology term
igvfagent catalog get-entity CPX-2189          # Complex Portal
igvfagent catalog get-entity R-HSA-8964038     # Reactome pathway

# Parallel fan-out over genes + variants + regulatory in a region
igvfagent catalog search-region "chr19:44,907,000-44,910,000"
igvfagent catalog search-region "19:44.9M-44.92M" --include genes,variants

# Edge query by semantic relationship
igvfagent catalog find-associations APOE --relationship pharmacological
igvfagent catalog find-associations rs429358 --relationship genetic \
                                    --filters "p_value=lte:5e-8"

# LD proxies for an index variant
igvfagent catalog find-ld rs429358 --r2-threshold 0.6 --ancestry EUR,AFR

# Translate one ID into all of its cross-references
igvfagent catalog resolve-id rs429358    # → rsid + spdi + hgvs + ca_id
igvfagent catalog resolve-id APOE        # → symbol + ENSG + HGNC + Entrez + synonyms

# Enumerate edge endpoints by semantic category
igvfagent catalog list-sources --category regulatory
igvfagent catalog list-sources --endpoint variants_genes
```

## Semantic relationships (`--relationship`)

| Category | Edges walked |
|---|---|
| `genetic` | variants_genes, variants_phenotypes, variants_diseases, variants_proteins, genes_diseases, genes_phenotypes, studies_phenotypes |
| `regulatory` | variants_genomic_elements, genes_genes, genomic_elements_genes, genomic_elements_biosamples, regulatory_regions_genes, motifs_proteins |
| `physical` | proteins_proteins, complexes_proteins |
| `functional` | genes_pathways, complexes_terms, proteins_terms, pathways_pathways, gene_products_terms (GO annotations) |
| `pharmacological` | variants_drugs, genes_drugs, proteins_drugs, drugs_diseases |
| `ld` | variants_variant_ld |
| `coding` | variants_coding_variants, coding_variants_phenotypes |
| `transcription` | genes_transcripts, genes_proteins, genes_structure |
| `all` | everything |

## ID auto-detection (the IDParser regex table)

| Pattern | Entity | Form | Example |
|---|---|---|---|
| `rs\d+` | variant | rsid | `rs429358` |
| `NC_\d+:pos:ref:alt` | variant | SPDI | `NC_000019.10:44908683:T:C` |
| `(NM|NP|ENST|ENSP):c.…` | variant | HGVS | `ENSP00000252934:p.Cys130Arg` |
| `CA\d+` | variant | ca_id | `CA388851` |
| `ENSG\d+` | gene | ensembl | `ENSG00000130203` |
| `HGNC:\d+` | gene | hgnc | `HGNC:613` |
| `ENTREZ:\d+` | gene | entrez | `ENTREZ:348` |
| `ENST\d+` | transcript | ensembl | `ENST00000252934` |
| `ENSP\d+` | protein | ensembl | `ENSP00000252934` |
| UniProt | protein | uniprot | `P02649` |
| `MONDO:\d+` / `EFO:\d+` / `DOID:\d+` / `HP:\d+` | ontology_term | (various) | `MONDO:0004975` (Alzheimer's) |
| `GO:\d+` / `UBERON:\d+` / `CL:\d+` / `CHEBI:\d+` / `OBA:\d+` | ontology_term | (various) | `GO:0006357` |
| `DB\d+` / `CHEMBL\d+` | drug | drugbank / chembl | `DB00945` |
| `CPX-\d+` | complex | complex_portal | `CPX-2189` |
| `R-HSA-\d+` | pathway | reactome | `R-HSA-8964038` |
| `GCST\d+` | study | gwas_catalog | `GCST90027158` |
| _everything else_ | gene | symbol | `APOE`, `TP53` |

## Filter DSL

Same syntax as ``portal``'s field_filters:

| Clause | Meaning |
|---|---|
| `label=eqtl` | equality |
| `method=GTEx,FANTOM5` | list (repeated params) |
| `label!=conservative` | negation |
| `log10pvalue=gte:7.30` | range op |
| `p_value=lte:5e-8` | **shortcut** — auto-translated to `log10pvalue=gte:7.301` |

Multi-clause separator: `;`.

## Region syntax

Accepts `chr19:44,907,000-44,910,000`, `19:44.9M-44.92M`, `chr1:1K-2K`,
`chrX:100000-200000`. The K/M/G suffixes follow the
**comma-friendly UCSC convention**. Coordinates are 0-based half-open
(the catalog's standard).

## P-value convention

The catalog stores `log10pvalue = -log10(P)`. To filter for
**GWAS-significant** associations:

```
--filters "p_value=lte:5e-8"       # this skill translates →
                                    # log10pvalue=gte:7.301
```

## What this skill adds over `kg` / `kg-mirror` / `igvf_client.catalog_get`

| Capability | Before | After |
|---|---|---|
| Universal `get-entity` with 20+ ID auto-detection | only gene + variant heuristics | ✓ |
| `search-region` parallel fan-out | manual region builds | ✓ |
| `find-associations` by semantic category | per-edge calls | ✓ |
| `find-ld` with r²/D'/ancestry buckets | summary endpoint only | ✓ |
| `resolve-id` cross-reference projection | none | ✓ |
| `list-sources` per-endpoint catalog | none | ✓ |
| Filter DSL with `p_value → log10pvalue` translation | none | ✓ |
| EDGE_ENDPOINTS registry (30+ edges with metadata) | scattered dicts | ✓ |
| Pagination metadata block | none | ✓ |

## License posture

Apache-2.0. Clean-room reimplementation of IGVF-DACC/igvf-catalog-mcp
(MIT) — the catalog REST contract is a documented set of facts; the
upstream IMPLEMENTATION_GUIDE is consulted as a factual reference, and
the live catalog endpoint is the wire-level source of truth. No source
code is copied. Stdlib-only — no `httpx` / `mcp` / `pydantic` runtime
deps. No GPL anywhere in this stack.
