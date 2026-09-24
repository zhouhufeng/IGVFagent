# Finding and retrieving data

IGVF Portal, Catalog and ENCODE access, the processed-first Portal lineage, dataset explanation, literature, GEO, Synapse, ChIP-Atlas, the Perturbation Catalogue and local reference databases.

**On this page**

- [IGVF / ENCODE / Knowledge Graph client](#igvf--encode--knowledge-graph-client)
- [Catalog / Portal / ENCODE overviews](#catalog--portal--encode-overviews)
- [Pathway databases (KEGG + Reactome + WikiPathways, integrated locally)](#pathway-databases-kegg--reactome--wikipathways-integrated-locally)
- [Start from what the IGVF Portal already computed (`processed lineage`)](#start-from-what-the-igvf-portal-already-computed-processed-lineage)
- [Data illustration and interpretation](#data-illustration-and-interpretation)
- [Reference skill (literature retrieval, validation, design)](#reference-skill-literature-retrieval-validation-design)
- [GEO retrieval](#geo-retrieval)
- [Perturbation Catalogue retrieval](#perturbation-catalogue-retrieval)
- [Human Transcription Factors database](#human-transcription-factors-database)
- [ChIP-Atlas reprocessed peak archive](#chip-atlas-reprocessed-peak-archive)
- [Synapse / Sage Bionetworks retrieval](#synapse--sage-bionetworks-retrieval)

## IGVF / ENCODE / Knowledge Graph client

The thin client behind every skill — direct HTTP calls and AQL.

```bash
python3 Scripts/igvf_client.py check
python3 Scripts/igvf_client.py catalog-api /
python3 Scripts/igvf_client.py catalog-files --limit 25
python3 Scripts/igvf_client.py gene TP53 --limit 10
python3 Scripts/igvf_client.py variant rs58658771 --limit 10
python3 Scripts/igvf_client.py encode-search --type Experiment --param assay_title=ATAC-seq
python3 Scripts/igvf_client.py aql "FOR doc IN genes LIMIT 5 RETURN doc"
```

## Catalog / Portal / ENCODE overviews

High-level inventory and smoke summaries.

```bash
python3 Scripts/igvf_data_skills.py catalog-smoke --limit 10
python3 Scripts/igvf_data_skills.py overview --limit 25
python3 Scripts/igvf_data_skills.py encode-overview --limit 25
python3 Scripts/igvf_data_skills.py encode-smoke --limit 10
python3 Scripts/igvf_data_skills.py encode-export-csv --type Experiment --param assay_title=ATAC-seq
python3 Scripts/igvf_frontpage_summary.py refresh --update-readme
```

## Pathway databases (KEGG + Reactome + WikiPathways, integrated locally)

Pulls the **current** releases, normalises every gene identifier into one
namespace, merges pathways that different databases describe differently, and
loads the result into the local knowledge graph — so pathway structure is
retrieved once and then queried offline. Each stored fact records which
databases assert it, so agreement between them is visible rather than
flattened into a duplicate. The integration method follows IntPath
(Zhou et al. 2012, cited under [References](../references.md)); the data is pulled
fresh, not taken from that release.

```bash
igvfagent pathwaydb sources                    # databases, licences, coverage
igvfagent pathwaydb pull                       # download current releases (~180 MB)
igvfagent pathwaydb pull --kgml                # + KEGG maps, for typed relations
igvfagent pathwaydb build --relations          # normalise, merge, load into the KG
igvfagent pathwaydb query --genes CLU,BIN1,PICALM,APOE,TREM2
igvfagent pathwaydb evaluate --agreement       # score the unification criteria
igvfagent pathwaydb status                     # what is cached, how old, which release
```

A pull on 2026-08-29 (KEGG 2026/08/27 · WikiPathways 20260810 · Reactome
current) yielded **211,120 gene-pathway memberships over 4,016 pathways and
14,870 genes**, plus **73,925 typed gene-gene relations** (PPrel 54,002 ·
ECrel 15,996 · GErel 3,197 · PCrel 730) — against 582 pathways in the 2012
IntPath release. Downloads are cached under `Data/References/PathwayDB/`
(gitignored) with a per-file release stamp and SHA-256, so a figure can be
traced back to the exact release behind it.

Pathways are unified across databases by IntPath's longest-common-subsequence
rule on pathway names, corroborated by word overlap and gene membership, and
constrained by Reactome's own parent/child hierarchy so curated fact outranks
string similarity. `igvfagent pathwaydb evaluate` scores each criterion on the
current data — the default reaches zero known-wrong merges where the published
rule alone reaches 3.4%. Full method, CLI reference and measurements:
**[`Docs/PATHWAYDB.md`](../../PATHWAYDB.md)**.

BioCyc/HumanCyc is not fetched — verified 2026-08-29, its download requires a
subscription "costing at least $5,000", so it cannot be pulled reproducibly.
If you hold a licence, integrate your own export with
`--extra-gmt HumanCyc=<file.gmt>`; it flows through the same normalisation,
unification and KG ingestion as the fetched sources, and nothing licensed is
redistributed.

## Start from what the IGVF Portal already computed (`processed lineage`)

Given any IGVF accession, IGVFagent now walks the links the Portal records
around it before it downloads anything. The links are the ones in the Portal's
data-model diagrams:
- `input_for` and `input_file_sets`, through intermediate, principal,
  pseudobulk, model and prediction sets;
- `related_measurement_sets` (the other modality of a multiome);
- `auxiliary_sets` (MULTI-seq, hashing, guides);
- `samples`, then `barcode_map`, then the curated barcode set;
- `seqspecs`, `onlist_files` and `documents`;
- `derived_from`, for a model's training data;
- `large_scale_loci_list`;
- the `QualityMetric` objects the uniform pipeline published.

The walk is directed. Outputs are followed downward and inputs upward, and
datasets from the same sample pool are listed but not expanded. So the "start
here" table only offers files actually built from your accession. Prediction
sets that apply a model trained on your data to other data are listed
separately from results of your data.

```bash
igvfagent processed lineage IGVFDS9875NBZW          # report.md + plan.json + lineage graph + figure
igvfagent processed fetch IGVFDS9875NBZW --max-gb 20   # the processed files + the Portal's QC, not the reads
igvfagent explain explain IGVFDS9875NBZW --download    # now fetches processed results first (--include-raw for reads)
igvfagent processed selftest                           # offline fixture Portal, 28 checks
```

On IGVFDS9875NBZW, a snRNA-seq 10x multiome with MULTI-seq, the walk finds:
- the uniform-pipeline h5ad, 3.7 GB, with its QC (74.9% pseudoaligned, 93.4%
  of reads on the barcode onlist);
- the public ATAC fragments, 5.6 GB, with 29.3% duplicates;
- cell annotations from the principal analysis;
- the MULTI-seq hashing table and the barcode-to-sample map.

That is about 18 GB of processed files, instead of 155 GB of controlled reads
across the RNA and ATAC measurement sets and the MULTI-seq auxiliary set.
Starting from the raw K562 multiome IGVFDS3910GBDQ, it reaches the scE2G
prediction sets IGVFDS5428HHMB and IGVFDS2032HBUP and their model
IGVFDS6290TXZY. Linked objects the credentials cannot see are reported as
"not visible (403)" rather than dropped. `raw-pipeline plan` and
`processed find/plan` now use the same walk.

Agent tools: `portal_lineage` (the agent's first call for any IGVF accession)
and `processed_fetch`.

The lineage report also says what the data is about and how it was built,
following the links in the lab submission diagrams (see
[IGVF Portal data model](../../Architecture/IGVF_PORTAL_DATA_MODEL.md)):

- superseded sets, and the set that replaces them;
- the construct library (guide, reporter or editing-template library) and its integrated guide or element tables;
- the sample tree: sorted fractions, time points, treatments and CRISPR modifications;
- a prediction's phenotypes, assessed genes, cell type and external training data;
- the column-definition document of each tabular output;
- the analysis step and software that made each file, and publications.

## Find IGVF data by topic (`processed discover`)

When a question names a phenotype, tissue, gene or kind of data rather than an
accession, `processed discover` finds the Portal objects about it. Your words
are matched to the terms the Portal actually uses (so "coronary artery
disease" reaches the Portal's phenotype term, and "heart" also reaches "heart
left ventricle"), and each type is searched on its own link fields:
phenotypes and assessed genes on PredictionSets, targeted genes and library
type on MeasurementSets, sample terms everywhere. Superseded sets are left
out, and a search with no match lists the nearest real values.

```bash
igvfagent processed discover --phenotype "coronary artery disease"
igvfagent processed discover --tissue heart --types MeasurementSet,AnalysisSet,PredictionSet
igvfagent processed discover --gene GATA1                       # TF binding models and predictions
igvfagent processed discover --prediction-type "element-gene links" --tissue heart
igvfagent processed discover --library-type "guide library" --tissue K562
```

Follow any accession it returns with `processed lineage`. Agent tool: `processed_discover`.

## Data illustration and interpretation

```bash
python3 Scripts/data_illustration_interpretation.py explain \
  '<igvf-portal-url>/curated-sets/IGVFDS2544COZH/'
python3 Scripts/data_illustration_interpretation.py explain \
  '<encode-portal-url>/search/?type=Annotation&searchTerm=encode-re2g&status!=archived'
python3 Scripts/data_illustration_interpretation.py explain IGVFDS2544COZH \
  --download --max-download-gb 2
python3 Scripts/data_illustration_interpretation.py write-playbook
```

## Reference skill (literature retrieval, validation, design)

The Reference Skill is the literature arm of the agent: it pulls publications
across PubMed / PMC, bioRxiv, medRxiv, arXiv, Semantic Scholar, and OpenAlex,
weights top-tier journals (Nature / Cell / Science families, NEJM, NAR,
Bioinformatics, Genome Biology / Research, eLife, PNAS), and feeds three
distinct functions used by the internal orchestrator's Evaluation Agent.

```bash
# 1) learn — what does the field do for this topic / assay / biosample?
python3 Scripts/reference_skill.py learn \
  --topic '10x multiome human putamen' --limit 30 --top 15

# 2) validate — has anyone seen these genes / variants / regulatory elements
#               before? cross-check IGVFagent discoveries against literature
python3 Scripts/reference_skill.py validate \
  --input Docs/AdvancedVariantAnalysis/<run>/<label>_summary_stats.csv \
  --context 'putamen Parkinson' --limit-per-item 5

# 3) design — given an IGVF data type, recommend a workflow + cognate studies
#             + matching IGVF Portal AnalysisSets
python3 Scripts/reference_skill.py design \
  --data-type parse_split_seq --assay-title 'Parse SPLiT-seq'

# generic multi-source search and playbook
python3 Scripts/reference_skill.py search --query 'enhancer-gene linkage rE2G ABC' --top 20
python3 Scripts/reference_skill.py write-playbook
```

All API responses are cached under `Data/Cache/References/<source>/` with a
14-day TTL; re-querying is free. Reports under
`Docs/References/<timestamp>_<subcommand>_<label>/`.

## GEO retrieval

Search NCBI GEO Series, parse SOFT metadata, list FTP file inventories,
download supplementary files, and produce a tidy sample sheet for any
downstream RNA-seq / ATAC / ChIP analysis. Pure stdlib + `requests` — no
Bioconductor dependency. Playbook:
[`Docs/Skills/GEO_RETRIEVAL_SKILLS.md`](../../Skills/GEO_RETRIEVAL_SKILLS.md).

```bash
# 1) Keyword search → top GSEs (title, organism, summary, n samples)
igvfagent geo search --query "GM12878 RNA-seq" --max-results 10

# 2) Pull the SOFT family file and parse series + sample metadata
igvfagent geo series --gse GSE9574 --label nci_breast

# 3) List supplementary files on the GEO FTP site for one Series
igvfagent geo list-files --gse GSE9574

# 4) Download chosen supplementary files (size-capped)
igvfagent geo download --gse GSE9574 --max-download-gb 5 \
  --extensions .txt.gz .tsv.gz .csv.gz

# 5) Write a tidy sample sheet (one row per GSM, condition columns)
igvfagent geo sample-sheet --gse GSE9574 --label nci_breast
```

Outputs land in `Docs/GEO/<timestamp>_<label>/` with the SOFT-derived
metadata JSON, file inventory CSV, sample sheet CSV, and a
human-readable summary. The agent runtime exposes `geo_search`,
`geo_series`, and `geo_download` as tools.

## Perturbation Catalogue retrieval

Pulls metadata and per-row perturbation effects from the public
**Perturbation Catalogue** Search API, which indexes ~1,222 datasets
across MAVE (DMS / VAMP-seq family), CRISPR screens (DepMap, Project
Score, Project Achilles, etc.), and Perturb-seq (Replogle 2022,
Nadig 2025, X-Atlas/Orion 2025). Playbook:
[`Docs/Skills/PERTURBATION_CATALOG_SKILLS.md`](../../Skills/PERTURBATION_CATALOG_SKILLS.md).

```bash
# 1) Catalogue-wide summary
igvfagent perturb-catalog summary

# 2) Global gene search (one row per perturbed gene, all modalities)
igvfagent perturb-catalog search --query BRCA1 --size 10

# 3) Modality-scoped search with filters
igvfagent perturb-catalog search-modality --modality mave \
    --perturbation-gene-name TP53 \
    --perturbation-position 100_300 \
    --effect-score-name vamp_score \
    --effect-score-value 0.5_1.0 \
    --dataset-limit 25
igvfagent perturb-catalog search-modality --modality crispr-screen \
    --perturbation-gene-name BRCA1 --dataset-limit 25
igvfagent perturb-catalog search-modality --modality perturb-seq \
    --query "lung cancer" --dataset-limit 10

# 4) Full dataset record + paginated per-perturbation rows
igvfagent perturb-catalog dataset --dataset-id replogle_2022_k562_essential_normalized
igvfagent perturb-catalog dataset-rows --modality perturb-seq \
    --dataset-id replogle_2022_k562_essential_normalized \
    --limit 500 --offset 0

# 5) Perturb-seq GSEA hallmark/pathway table
igvfagent perturb-catalog gsea --query BRCA1 --size 50

# 6) Bulk dataset download (auto-detects .csv.gz / .json / .zip)
igvfagent perturb-catalog download --modality perturb-seq \
    --dataset-id replogle_2022_k562_essential_normalized

# 7) End-to-end gene-centric pipeline (the headline command for
#    "what perturbation data exists for gene X?")
igvfagent perturb-catalog pipeline --gene BRCA1 --dataset-limit 10
```

The agent runtime registers `perturb_catalog_summary`,
`perturb_catalog_search`, `perturb_catalog_search_modality`,
`perturb_catalog_dataset`, `perturb_catalog_dataset_rows`,
`perturb_catalog_gsea`, and `perturb_catalog_pipeline` as tools.
Live smoke test on `--gene BRCA1` returns **1,191 CRISPR-screen
datasets** (509 significant), **7 Perturb-seq datasets**, and the
canonical BRCA1 GSEA signature (`HALLMARK_E2F_TARGETS`,
`HALLMARK_G2M_CHECKPOINT`, `HALLMARK_MYC_TARGETS_V1`). Outputs land
under `Docs/Perturbation/<timestamp>_<gene>/` with one JSON per
modality plus a markdown report; bulk downloads go to
`Data/Perturbation/Downloads/<modality>/`.

## Human Transcription Factors database

Local SQLite mirror of the Lambert/Jolma/Hughes **Human Transcription
Factors** database v1.01 ([humantfs.ccbr.utoronto.ca](https://humantfs.ccbr.utoronto.ca);
Lambert et al., *Cell* 2018): 2,765 assessed proteins of which **1,639
are curated TFs**, with DNA-binding-domain family, binding mode, motif
status and cross-references. Playbook:
[`Docs/Skills/HUMANTFS_SKILL.md`](../../Skills/HUMANTFS_SKILL.md).

```bash
igvfagent humantfs build-db                     # ~2 MB, seconds
igvfagent humantfs is-tf --genes "FOXP3,EOMES,ACTB,GAPDH"
igvfagent humantfs lookup --genes SATB2
igvfagent humantfs list --dbd "C2H2 ZF" --with-motif
igvfagent humantfs families
```

`is-tf` distinguishes *unassessed* from *not a TF* — absence from the
database is not evidence against. The data is downloaded at build time,
never vendored; cite Lambert 2018 for derived tables.

## ChIP-Atlas reprocessed peak archive

Browse and pull from the [ChIP-Atlas](https://chip-atlas.org) (Ohta/Oki/DBCLS)
reprocessed archive of public ChIP-seq / ATAC-seq / DNase-seq / Bisulfite-seq
peaks — a clean-room, stdlib-only client over the public HTTP surface. Ten
genomes, four −log10(q) thresholds, per-experiment BigWig/BigBed/BED, assembled
all-peaks BEDs, pre-computed Target-Genes tables, and queued TF-enrichment jobs.
Playbook: [`Docs/Skills/CHIPATLAS_SKILL.md`](../../Skills/CHIPATLAS_SKILL.md).

```bash
igvfagent chipatlas list-antigens   --genome hg38 --cell-type Blood
igvfagent chipatlas search          --genome hg38 --antigen GATA1 --cell-type Blood
igvfagent chipatlas target-genes    --antigen H3K4me3 --distance 5000
igvfagent chipatlas submit-enrichment --genes my_genes.txt --genome hg38
```

## Synapse / Sage Bionetworks retrieval

Discover and download from [Synapse](https://www.synapse.org) — a clean-room
client over the public REST API (no `synapseclient` dependency). Anonymous read
of entity metadata, child listing, recursive walks, and full-text search; PAT-
authenticated download (`SYNAPSE_AUTH_TOKEN`) for controlled-access deposits
(PsychENCODE, AMP-AD/PD, ROSMAP) that IGVF distributes off-Portal. This skill
powers the **Deng 2024 cortex lentiMPRA** benchmark. Playbook:
[`Docs/Skills/SYNAPSE_RETRIEVAL_SKILLS.md`](../../Skills/SYNAPSE_RETRIEVAL_SKILLS.md).

```bash
igvfagent synapse entity   --syn syn21392931                 # metadata (anon)
igvfagent synapse walk     --syn syn21392931 --max-depth 3   # recursive walk
igvfagent synapse search   --query "lentiMPRA cortex" --limit 20
igvfagent synapse download --syn synXXXXXXXX --out-dir Data/Input   # needs PAT
```

---

[← Documentation index](../README.md) · [Project README](../../../README.md)
