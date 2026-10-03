# Bulk genomics, RNA-seq and proteomics

The ENCODE bulk-genomics pipeline, bulk RNA-seq, and the protein-interaction knowledge graph.

**On this page**

- [ENCODE bulk-genomics pipeline](#encode-bulk-genomics-pipeline)
- [Bulk RNA-seq analysis](#bulk-rna-seq-analysis)
- [Proteomics & PPI knowledge graph](#proteomics--ppi-knowledge-graph)

## ENCODE bulk-genomics pipeline

End-to-end retrieval, description, peak QC, super-enhancer calling,
SCREEN cCRE integration, and IGV-style browser-SVG visualization for
the major ENCODE bulk assays — ChIP-seq (TF + Histone), ATAC-seq,
DNase-seq, Hi-C, capture Hi-C, ChIA-PET, plus RNA-seq / MNase-seq /
FAIRE-seq / CAGE / RAMPAGE for retrieval & description.

```bash
# 1) Discover experiments by assay × biosample × target
igvfagent encode retrieve --assay 'Histone ChIP-seq' \
  --target H3K27ac --biosample K562 --assembly GRCh38 --limit 50

# 2) Per-file inventory for one or more accessions
igvfagent encode manifest --accessions ENCSR000AKP --label k562_h3k27ac

# 3) Pull files under a size cap, filter by file format
igvfagent encode download \
  --manifest Data/Manifests/ENCODE/<files.csv> \
  --max-download-gb 5 --formats bed bigWig

# 4) Plain-language description of one experiment
igvfagent encode describe --accession ENCSR000AKP

# 5) Peak QC (count, width, score, per-chromosome) from a BED file
igvfagent encode analyze-peaks --bed peaks.bed.gz

# 6) ROSE-style super-enhancer call (H3K27ac / BRD4 / MED1 / P300)
igvfagent encode super-enhancers --bed h3k27ac_peaks.bed \
  --stitching-distance 12500 --tss-bed tss.bed --tss-distance 2000

# 7) Overlay peaks with SCREEN cCRE classes (PLS / pELS / dELS / CTCF)
igvfagent encode integrate-ccre --bed peaks.bed

# 8) IGV-style multi-track SVG browser view
igvfagent encode browser --region chr19:44903000-44912000 \
  --track 'H3K27ac peaks:peaks.bed' \
  --track 'ATAC peaks:atac_peaks.bed' \
  --with-ccre --label apoe_locus
```

### Everything a portal holds for one cell line

"Summarise the ENCODE and IGVF data for GM12878, with tables and plots" is
one command, not a search URL:

```bash
igvfagent biosample-census --biosample GM12878          # both portals
igvfagent biosample-census --biosample K562 --portal encode --status all
```

It counts every object type on both portals through their **structured
sample-term filters** (`biosample_ontology.term_name` on ENCODE,
`samples.sample_terms.term_name` on IGVF), resolves the term case-insensitively
through the shared ontology id (GM12878 is `EFO:0002784` on both), and writes
`Docs/BiosampleCensus/<run>/` with `report.md`, `totals_by_type.csv`,
`facet_counts.csv`, per-type item tables, and SVG (+PNG) figures for assays,
ChIP targets, labs, annotation types, file formats and release years. Free-text
search is deliberately not used for counting: `searchTerm=GM12878` on the IGVF
Portal matches thousands of unrelated sets. A zero-hit search, which both
portals answer with HTTP 404, is recorded as 0. The agent tool is
`biosample_portal_census`.

`explain_dataset` on a search URL that matches nothing now says *why* --
which filter field the portal does not know (legacy names such as
`biosample_term_name` are rewritten automatically) or which value is
misspelt (`gm12878` vs `GM12878`) -- instead of writing an empty report.

The agent runtime exposes `encode_retrieve`, `encode_describe`,
`encode_super_enhancers`, `encode_integrate_ccre`, and `encode_browser`
as tools, so a single `igvfagent ask` can drive the full pipeline:

```bash
igvfagent ask "Find K562 H3K27ac ChIP-seq experiments on GRCh38, pick \
  the highest-quality one, describe it, and call super-enhancers from \
  its peak BED. Then overlay the super-enhancers against SCREEN cCREs \
  and produce a browser view of the APOE locus."
```

## Bulk RNA-seq analysis

Counts QC, sample PCA + correlation heatmap, differential expression
(pyDESeq2 when available, Welch's t-test on log-CPM with BH FDR
fallback), volcano + MA + top-DEG z-scored heatmap, and DEG → cCRE
linkage via the IGVF Catalog `/api/genes/genomic-elements` endpoint.
Playbook:
[`Docs/Skills/RNASEQ_ANALYSIS_SKILLS.md`](../../Skills/RNASEQ_ANALYSIS_SKILLS.md).

```bash
# 1) Counts QC (library size, gene detection, mito %, top genes)
igvfagent rnaseq qc --counts counts.tsv --label k562_vs_gm12878

# 2) PCA + sample correlation heatmap
igvfagent rnaseq pca --counts counts.tsv --sample-sheet samples.csv \
  --label k562_vs_gm12878

# 3) Differential expression (uses pyDESeq2 if installed, else Welch + BH)
igvfagent rnaseq deg --counts counts.tsv --sample-sheet samples.csv \
  --condition-col condition --treated K562 --control GM12878 \
  --label k562_vs_gm12878

# 4) Link significant DEGs to controlling cCREs via IGVF Catalog
igvfagent rnaseq link-cre \
  --deg Docs/RNAseq/<run>/k562_vs_gm12878_deg.csv \
  --label k562_vs_gm12878 --padj-threshold 0.05 --max-genes 200

# 5) End-to-end QC → PCA → DEG → cCRE linkage in one call
igvfagent rnaseq pipeline --counts counts.tsv --sample-sheet samples.csv \
  --condition-col condition --treated K562 --control GM12878 \
  --label k562_vs_gm12878
```

Outputs include `<label>_deg.csv`, `<label>_deg_to_cre.csv`, the
volcano/MA/heatmap/PCA PNGs, and a markdown report under
`Docs/RNAseq/<timestamp>_<label>/`. Chains naturally with `geo` (pull
counts) and `se-targets` (cross-reference up-regulated genes against
SE-driven targets).

## Proteomics & PPI knowledge graph

End-to-end protein/PPI skill: pulls and version-tracks
**BioGRID**, **IntAct**, **HuRI**, **Reactome**, **KEGG**, **UniProt
id-mapping**, and **all 214 IGVF Portal protein-slim MeasurementSets**
(plus the 14 PPI-score files, the DUAL-IPA fluorescence file, and the
UniProt / protein-language-model reference files). Integrates everything
into a local SQLite knowledge graph (`Data/Proteomics/KG/proteomics.sqlite`)
with deduplicated edges keyed on `(id_a, id_b, source, source_id)`.
Provides summary stats, network visualizations, per-IGVF-assay example
figures from real Portal files, and a literature survey for the IGVF
protein assays via the Reference skill. Playbook:
[`Docs/Skills/PROTEOMICS_SKILLS.md`](../../Skills/PROTEOMICS_SKILLS.md).

```bash
# 1) Download (only sources that changed upstream)
igvfagent proteomics download --source all
# Or one at a time:
igvfagent proteomics download --source biogrid
igvfagent proteomics download --source intact
igvfagent proteomics download --source huri
igvfagent proteomics download --source reactome
igvfagent proteomics download --source kegg --kegg-max-pathways 400

# 2) Version manifest (local + upstream probe)
igvfagent proteomics versions

# 3) Pull all IGVF Portal protein assays + actual files (semi-qY2H,
#    DUAL-IPA, VAMP-seq, MAVE) from the public S3 bucket
igvfagent proteomics igvf-protein

# 4) Build the integrated PPI-KG (SQLite)
igvfagent proteomics build-kg --sources all

# 5) Summary statistics (per source / evidence type / detection method,
#    top-30 hubs)
igvfagent proteomics kg-stats --label initial

# 6) Network visualizations (degree distribution, top-30 hubs,
#    per-source breakdown, ego graph for a query gene)
igvfagent proteomics kg-visualize --gene TP53 --label tp53

# 7) Literature survey for IGVF protein assays — restricted to the
#    Nature/Cell/Science journal family
igvfagent proteomics assay-survey --label may2026

# 8) Per-assay example histograms generated from real IGVF Portal files
#    (semi-qY2H v1/v2/v3, DUAL-IPA, plus VAMP-seq family / MAVE if
#    those Portal files are pulled)
igvfagent proteomics assay-figures --label demos

# 9) End-to-end orchestrator
igvfagent proteomics pipeline --label may2026 --gene TP53 \
  --sources biogrid,intact,huri,reactome,kegg,igvf

# 10) VAMP-seq deep analysis — pull canonical MaveDB scoresets, run the
#     full Matreyek/Suiter/Clausen/Coyote-Maestas analysis pipeline,
#     and inventory the IGVF Portal raw VAMP-seq experiments
igvfagent proteomics vampseq-pull              # PTEN, TPMT, VKOR, PRKN,
                                               # CYP2C9, NUDT15 from MaveDB
igvfagent proteomics vampseq-analyze --gene PTEN --label pten_deep
igvfagent proteomics vampseq-analyze           # all 6 catalogued targets
igvfagent proteomics vampseq-inventory --label igvf_f9
```

The skill maintains a `Data/Proteomics/_versions.json` manifest with
URL, sha256, and record count per source. `update` only re-fetches
sources where the upstream version differs from the local one. KEGG
calls are throttled (default 0.4 s/req) to respect their TOS; IntAct
defaults to the smaller `intact-micluster.txt` rather than the full
~700 MB `intact.zip`. HuRI uses Ensembl gene IDs; run
`proteomics download --source uniprot` to populate the `id_map` table
for UniProt ↔ Ensembl ↔ Symbol harmonization.

The `vampseq-analyze` subcommand follows the canonical pipeline
distilled from Matreyek *Nat Genet* 2018, Suiter *eLife* 2020, Clausen
*Nat Commun* 2024, and Coyote-Maestas *Nat Commun* 2024 (MultiSTEP) —
producing six publication-grade plots per gene: score distribution,
residue × AA heatmap (the iconic VAMP-seq view), per-residue mean ±
IQR with a domain track, replicate concordance scatter with Pearson
*r*, abundance-class breakdown, and a cumulative ranked-variant curve.
Domain tracks are pre-curated for PTEN (PIP4-bind / Phosphatase / C2 /
C-tail), TPMT, VKOR, PRKN (Ubl / Linker / RING0 / RING1 / IBR / RING2),
CYP2C9, and NUDT15.

The `vampseq-inventory` subcommand decodes the alias scheme on the
IGVF Portal MeasurementSets (`<lab>:<GENE>-DMS-<antibody>-Tile<i>-Replicate<j>-Bin<k>`)
into per-gene coverage matrices: the 144 MultiSTEP sets resolve to
**F9 (Coagulation Factor IX)** across 3 tiles × 4 bins × 4 replicates ×
5 antibody readouts (Light-chain, Heavy-chain, Strep-II-tag, and two
carboxylation-sensitive Gla-domain antibodies); the 36 plain VAMP-seq
sets cover **CYP2C19** and **G6PD**.

---

[← Documentation index](../README.md) · [Project README](../../../README.md)
