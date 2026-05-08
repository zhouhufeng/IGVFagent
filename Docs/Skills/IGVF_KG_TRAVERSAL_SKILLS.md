# Skill: IGVF Knowledge Graph traversal

Iteratively walks the IGVF Catalog Knowledge Graph (`api.catalogkg.igvf.org` + the underlying ArangoDB at `db.catalog.igvf.org`) starting from a single entity and assembles a unified evidence pack across direct neighbors, second-degree relations, and optional cross-skill enrichment (FAVOR, enhancer-gene linkage, IGVF single-cell datasets, prior literature).

Designed as the orchestrator-friendly **comprehensive context** tool: one CLI call → per-relation manifests + JSON evidence pack + markdown report.

## Subcommands

### 1. `gene <symbol>` — comprehensive gene-centric traversal

```bash
python3 Scripts/kg_traversal_skill.py gene APOE \
    --depth 2 --limit 50 \
    --max-variants 25 --subvariant-limit 10 \
    --call-favor --call-linkage --call-singlecell --call-literature \
    --literature-context Alzheimer cardiovascular \
    --label apoe_full
```

Default direct relations (each saved to its own manifest CSV):

- `variants` — variants on the gene
- `coding_variant_scores` — MutPred2 / ESM-1v predictions per coding change
- `regulatory_elements` — gene→cCRE links (CRISPRi / Perturb-seq / etc.)
- `transcripts` — gene model isoforms
- `proteins` — protein records linked to the gene
- `diseases` — gene→disease/phenotype links
- `pathways` — gene→pathway membership

At `--depth 2`, each variant additionally fans out into its own summary, QTL genes, phenotypes, biosamples (CRISPRi / MPRA), genomic-element overlaps, and prediction sets.

Optional side-calls:

- `--call-favor` — pulls FAVOR functional annotations for the gene region.
- `--call-linkage` — adds enhancer-gene linkage predictions for the gene region (rE2G / catalog regulatory-region links).
- `--call-singlecell` — searches the IGVF Portal for single-cell AnalysisSets that mention the gene; surfaces candidate datasets for downstream expression analysis with `Scripts/single_cell_data_skills.py` or `Scripts/splitseq_pipeline.py`.
- `--call-literature` — runs `Scripts/reference_skill.py validate` on the gene + your context terms.

### 2. `variant <id>` — variant-centric traversal

```bash
python3 Scripts/kg_traversal_skill.py variant rs429358 \
    --call-favor --call-literature --label apoe_e4_variant
```

Variant ID accepted as rsID, SPDI, HGVS, or chr:pos:ref:alt where the Catalog API understands it.

### 3. `region <chr:start-end>`

```bash
python3 Scripts/kg_traversal_skill.py region chr19:44903000-44912000 \
    --call-favor --label apoe_locus
```

Returns: genes overlapping the region, regulatory elements (cCREs) in the region, region-predictor enhancer-gene linkage rows, and (optional) FAVOR variant annotations.

### 4. `aql` — direct ArangoDB AQL pass-through

```bash
python3 Scripts/kg_traversal_skill.py aql \
    'FOR g IN genes FILTER g.name == "APOE" RETURN g'
```

## Outputs

Each run writes a timestamped folder under `Docs/KGTraversal/`:

  Docs/KGTraversal/<timestamp>_<label>/
    ├─ <entity>_<key>_report.md   # full markdown report
    ├─ evidence_pack.json         # complete JSON with every relation
    └─ Manifests/                  # one CSV per relation

## How this chains with other skills

- The `evidence_pack.json` is designed to be consumed by other skills — variant manifests can be fed to `Scripts/advanced_variant_analysis.py` or `Scripts/annotate_variant_list.py`; cCRE manifests can be fed to `Scripts/enhancer_gene_linkage_skills.py compare-sets`.
- The single-cell dataset hits can be drilled into with `Scripts/single_cell_data_skills.py manifest` or `Scripts/splitseq_pipeline.py manifest`.
- The literature manifest can be cross-checked with `Scripts/reference_skill.py validate`.
- The internal orchestrator (Plan → Action → Results → Evaluation) uses this skill as the primary 'comprehensive context' action: given a gene of interest from the planning step, this skill supplies all the multi-omic evidence the analysis and evaluation steps need.