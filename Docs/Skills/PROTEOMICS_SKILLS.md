# Skill: Proteomics & PPI knowledge graph

End-to-end proteomics skill for IGVFagent. Aggregates **BioGRID**,
**IntAct**, **HuRI**, **Reactome**, **KEGG**, **UniProt id-mapping**,
and **IGVF Portal protein assays** into a local SQLite knowledge graph
that is incrementally update-aware. Provides summary statistics,
network visualizations (degree distribution, hub bar plot, ego graph),
literature surveys for VAMP-seq family / MAVE / semi-qY2H / DUAL-IPA via
the Reference skill, and per-assay example figures from real Portal
files.

## Subcommands

### download
```
igvfagent proteomics download --source all
igvfagent proteomics download --source biogrid
igvfagent proteomics download --source intact
igvfagent proteomics download --source huri
igvfagent proteomics download --source reactome
igvfagent proteomics download --source kegg --kegg-max-pathways 400
igvfagent proteomics download --source igvf
igvfagent proteomics download --source uniprot   # idmap, optional
```
The skill stores `Data/Proteomics/_versions.json` with version, URL, sha256,
and record count per source. BioGRID is auto-resolved against the latest
`Release-Archive` listing. IntAct uses `Last-Modified` on
`psimitab/intact-micluster.txt`. HuRI is mostly static.

### versions / update
```
igvfagent proteomics versions
igvfagent proteomics update
```
`update` only re-fetches sources where the upstream version differs from
the local one (delta refresh).

### igvf-protein
```
igvfagent proteomics igvf-protein
```
Pulls all 214 protein-slim MeasurementSets, the 14 PPI-score files, the
DUAL-IPA stability file, the UniProt protein references, and the protein
language model files from the IGVF Portal. Saves metadata JSON + actual
TSV/CSV file content under `Data/Proteomics/Sources/IGVF/`.

### build-kg
```
igvfagent proteomics build-kg --sources all
igvfagent proteomics build-kg --sources biogrid,reactome,igvf
```
Ingests downloaded sources into `Data/Proteomics/KG/proteomics.sqlite`.
Schema: `interactions`, `proteins`, `pathways`, `pathway_membership`,
`igvf_evidence`, `id_map`, `versions`. Edges deduped on
`(id_a, id_b, source, source_id)`.

### kg-stats
```
igvfagent proteomics kg-stats --label initial
```
Summary report under `Docs/Proteomics/<ts>_initial/`:
total interactions, distinct proteins, pathways, IGVF evidence,
breakdowns by source / evidence type / detection method, top-30 hubs.

### kg-visualize
```
igvfagent proteomics kg-visualize --gene TP53 --label tp53
igvfagent proteomics kg-visualize --label snapshot
```
Produces `Plots/degree_distribution.png`, `Plots/top_hubs.png`,
`Plots/per_source.png`, and (when `--gene` given) an ego graph
`Plots/ego_<GENE>.png` (requires networkx).

### assay-survey
```
igvfagent proteomics assay-survey --label may2026
```
Calls the Reference skill (PubMed + Semantic Scholar + OpenAlex) to
retrieve current studies on:
- VAMP-seq (MultiSTEP)
- VAMP-seq
- MAVE
- semi-qY2H v1 / v2 / v3 (semi-quantitative yeast two-hybrid)
- DUAL-IPA
Filters to Nature / Cell / Science family journals (Nat Genet, Nat Methods,
Nat Biotechnol, Mol Cell, Mol Syst Biol, Nat Commun) and writes
`literature_survey.md` and `.json`.

### assay-figures
```
igvfagent proteomics assay-figures --label demos
```
Generates per-assay example histograms from the actual IGVF Portal files
pulled by `igvf-protein`. One PNG per assay under
`Docs/Proteomics/<ts>_demos_assay_figures/Plots/`.

### pipeline
```
igvfagent proteomics pipeline --label may2026 --gene TP53 \
  --sources biogrid,intact,huri,reactome,kegg,igvf
```
End-to-end: download → build-kg → kg-stats → kg-visualize → assay-figures
→ assay-survey. Use `--skip-download` to reuse cached files,
`--skip-literature` to skip the network-dependent literature scan,
`--max-rows N` to cap parser ingestion (useful for smoke tests).

### write-playbook
```
igvfagent proteomics write-playbook
```
Writes this document to `Docs/Skills/PROTEOMICS_SKILLS.md`.

## Storage

```
Data/Proteomics/
  _versions.json
  Sources/
    BioGRID/BIOGRID-ALL-X.Y.Z.tab3.zip
    IntAct/intact-micluster.txt
    HuRI/{HuRI.tsv, HI-union.tsv, Lit-BM.tsv}
    Reactome/{ReactomePathways.txt, UniProt2Reactome.txt,
              reactome.homo_sapiens.interactions.tab-delimited.txt}
    KEGG/{pathways_hsa.json, pathway_membership_hsa.json}
    UniProt/HUMAN_9606_idmapping_selected.tab.gz
    IGVF/{measurement_sets.json, ppi_files_meta.json,
          stability_files_meta.json, protein_reference_meta.json,
          protein_language_models.json, ppi_files/IGVFFI*}
  KG/proteomics.sqlite

Docs/Proteomics/<ts>_<label>/
  report.md
  Plots/*.png
  literature_survey.md (if assay-survey)
  per_assay_figures/Plots/*.png (if assay-figures)
```

## Cross-skill chaining

- `proteomics build-kg` → `kg gene <SYM>` to combine the proteomics PPI
  graph with the IGVF Catalog Knowledge Graph (variants, regulatory
  elements, pathways).
- `proteomics assay-survey` → `ref design --query "<assay>"` to draft
  a study design that mirrors the canonical analysis flow used by
  recent Nature/Cell/Science studies.
- `rnaseq deg` significant DEGs → query
  `proteomics kg-visualize --gene <SYM>` to see what proteins the
  upregulated genes interact with.
- `se-targets pipeline` → for each SE-target gene, look up upstream
  PPI partners with `proteomics kg-visualize` to suggest functional
  cofactors.

## Notes

- KEGG REST is throttled (default 0.4s/req) to respect their TOS;
  bulk distribution requires a license.
- IntAct full PSI-MITAB (`intact.zip`) is ~700MB — the skill defaults
  to the smaller `intact-micluster.txt`. Set `INTACT_BASE` env to
  override.
- HuRI uses Ensembl gene IDs. Run `proteomics download --source uniprot`
  to populate `id_map` for cross-resource harmonization.
