# Single-cell multiome data: where it actually lives

This document is the systemic overview produced alongside the
``multiome_survey`` skill.  It summarizes (a) the six public sources the
skill queries directly, (b) what each gives you and what it does not,
and (c) the other repositories that host multiome data but are not (yet)
auto-searchable from this skill — pointing you to where to go manually.

## Six sources surveyed by this skill

| source        | identifier on disk       | API endpoint                                              | what you get                                                                                  |
|---------------|--------------------------|-----------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| **IGVF**      | ``IGVFDS…`` AnalysisSets | IGVF Portal ``/search/`` (JSON)                           | Cell Ranger ARC tars, ATAC fragments BED, cell annotations TSV, peak matrices (.rds).         |
| **ENCODE**    | ``ENCSR…`` experiments   | ENCODE ``/search/`` (JSON)                                | Snippets from joint snRNA + snATAC Series, SHARE-seq Experiment records, alignments + signals. |
| **GEO**       | ``GSE…`` Series          | NCBI E-utilities (``esearch`` / ``esummary``) on ``gds``  | Sample sheet, supplementary processed files (matrices, fragment beds, h5ad if author uploaded). |
| **CELLxGENE** | dataset UUIDs            | ``https://api.cellxgene.cziscience.com/curation/v1``      | Curated H5AD with standardized cell ontology + assay metadata; ready for training.            |
| **HCA**       | HCA project UUIDs        | Azul ``/index/projects`` and ``/index/files``             | Project-level metadata + per-file download URLs (loom, h5ad, matrices, raw FASTQs).           |
| **Zenodo**    | numeric record ids       | ``https://zenodo.org/api/records``                        | Author-uploaded archives — captures data shared *outside* the standard repositories.          |

### Strengths and weaknesses, side by side

- **IGVF** — newest data, richest variant-to-function context, smallest catalog overall (1,954 AnalysisSets as of May 2026). Peak matrices ship as ``.rds`` (R-only).
- **ENCODE** — strong for SHARE-seq and rare paired-modality experiments, but native multiome coverage is small.
- **GEO** — broadest catalog by far (50+ GSEs match basic queries; the long tail of author submissions is here). Heterogeneous file naming; supplementary files require FTP scraping.
- **CELLxGENE Discover** — curated, schema-validated H5AD; the easiest pure-RNA half to download and train on. But the chromatin half of multiome is dropped during ingestion.
- **HCA Data Portal** — best for organized consortium projects; supports filtering on ``libraryConstructionApproach``.
- **Zenodo** — catches anything authors deposit alongside a paper (cluster labels, region-of-interest BEDs, custom models). No standardized schema.

## Where multiome data also lives (not auto-queried here)

The sources below host substantial multiome data but were left out of the
default skill either because (a) they require authenticated access that
this CLI shouldn't bake in, (b) their API needs careful schema-mapping
beyond the scope of a generic survey, or (c) coverage is comparatively
small relative to the six above.

| repository                                            | typical content                                          | how to access                                                                                                       |
|-------------------------------------------------------|----------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|
| **Allen Brain Cell (ABC) Atlas**                      | Yao et al. *Nature* 2023 whole-mouse-brain snATAC + 10x Multiome (1,687 nuclei across 33 clusters). | https://alleninstitute.github.io/abc_atlas_access/  — Python package, S3-backed direct download.                    |
| **Broad Single Cell Portal (SCP)**                    | Many published multiome studies (Buenrostro/Engreitz labs etc.) | https://singlecell.broadinstitute.org — REST API ``/single_cell/api/v1/`` (auth needed for some studies).            |
| **Synapse (Sage Bionetworks)**                        | PsychENCODE 2.0 multiome, AMP-AD, AMP-PD consortia.       | https://www.synapse.org — Python client; **most studies require login + access agreement**.                          |
| **dbGaP**                                              | Controlled-access human genomics; many multiome studies referenced here. | https://www.ncbi.nlm.nih.gov/gap — formal DAR/IRB approval required.                                                  |
| **EGA (European Genome-phenome Archive)**             | European equivalent of dbGaP; same access model.          | https://ega-archive.org — controlled access, EGA download client.                                                    |
| **ArrayExpress / BioStudies (EBI)**                   | European GEO-equivalent; some multiome Series.            | https://www.ebi.ac.uk/biostudies/api/v1/search                                                                       |
| **figshare**                                          | Long-tail dataset deposits attached to papers.            | https://api.figshare.com/v2/articles                                                                                 |
| **DDBJ Omics Archive (DOR)**                          | Japanese counterpart to GEO/SRA — small but non-overlapping. | https://ddbj.nig.ac.jp                                                                                              |
| **CIRM**                                              | California stem-cell repository; iPSC + multiome derivatives. | https://www.cirm.ca.gov                                                                                              |
| **Tabula Sapiens / Tabula Muris**                     | Standalone tissue atlases; mostly RNA-only but with some multimodal slices. | https://tabula-sapiens-portal.ds.czbiohub.org                                                                        |
| **DNAnexus / Terra (BroadFC)**                        | Many consortium workspaces share processed multiome (auth required). | https://www.dnanexus.com / https://app.terra.bio                                                                     |
| **NeMO Archive (Brain Initiative)**                   | BICCN multiome data including 10x Multiome + SHARE-seq.   | https://nemoarchive.org                                                                                              |

### When to use each

- **Training a foundation model** → start with **CELLxGENE Discover** (clean H5AD, ontology labels) for the RNA modality and **IGVF + GEO** for paired modalities.
- **Variant-to-function modeling** → **IGVF** is purpose-built for this; supplement with **CELLxGENE** for cell-type context.
- **Recreating a specific published analysis** → check **GEO** first (typical deposit), then **Zenodo** for analysis-time auxiliary files (cluster labels, derived features), then ABC / Synapse for consortium-specific releases.
- **Brain-specific work** → **Allen ABC Atlas** + **NeMO Archive** are the high-value mines; supplement with HCA brain projects.
- **Human disease cohorts** → likely controlled-access: **dbGaP** (US) or **EGA** (EU).  Plan around the DAR timeline.

## Filename / kind classification

The skill normalises every discovered file into a ``kind`` label so a
training pipeline can pull only what it cares about:

| kind            | matches                                                                                          |
|-----------------|---------------------------------------------------------------------------------------------------|
| ``matrix_rna``   | sparse gene count matrix, gene quantifications, filtered/raw feature-barcode matrix.             |
| ``matrix_atac``  | annotated sparse peak count matrix, cell-by-peak matrices.                                       |
| ``fragments``    | ATAC fragments BED (.bed.gz / .tsv.gz).                                                          |
| ``peaks``        | Peak call files (.bed / .narrowPeak).                                                            |
| ``annotations``  | Cell metadata / annotations / sample sheet (.tsv).                                               |
| ``alignments``   | BAM / aligned reads.                                                                              |
| ``index``        | BAI / TBI / CRAI index files.                                                                     |
| ``raw_reads``    | FASTQ / sequence reads.                                                                           |
| ``other``        | everything else.                                                                                  |

Filter ``download`` by these kinds with ``--only matrix_rna,fragments,annotations``.

## Privacy

Every endpoint URL is hex-encoded in ``Scripts/_endpoints.py``; no
hard-coded URLs or credentials appear in source.  All output paths are
repo-relative.  ``Data/MultiomeSurvey/`` is covered by the existing
``Data/*`` rule in ``.gitignore``, so downloaded payload never lands in
commits.

## How to run the skill end-to-end

```bash
# 1. Survey all six sources at once.
igvfagent multiome-survey survey-all --limit 100 --fetch-files

# 2. Re-build the unified manifest.
igvfagent multiome-survey manifest --label v1

# 3. Download a 20 GB training slice (RNA matrices + ATAC fragments + cell labels).
igvfagent multiome-survey download \
    --only matrix_rna,fragments,annotations \
    --max-download-gb 20

# 4. Refresh the on-disk inventory.
igvfagent multiome-survey inventory
```
