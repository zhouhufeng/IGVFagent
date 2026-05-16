# IGVF Integrated Data Layer — design

How IGVFagent will go from "many independent skills pulling from many
upstream resources" to a **single centralized data warehouse** that
feeds a multi-modal **IGVF Foundation Model**.

This doc covers four things:

1. **The centralized data structure** — which database(s) host what
2. **Raw → processed pipeline** — how each existing skill becomes a
   data producer that lands QC'd outputs in the warehouse
3. **Cross-source integration & embeddings** — how variants, genes,
   proteins, regulatory elements, regulatory networks, and cell
   profiles get joined into a single feature space
4. **Foundation-model training** — concrete architecture + pre-trained
   backbones + training stages

The proposed system is laptop-runnable today (DuckDB + Parquet +
SQLite + FAISS) and scales linearly to a cluster when the dataset
size demands it (drop-in S3/MinIO + Postgres/ArangoDB).

---

## 1. Architecture — three-tier medallion

```
┌────────────────────── BRONZE ──────────────────────┐
│ Raw, immutable copies of upstream data             │
│  ├ IGVF Portal files  (FASTQ, BAM, h5ad, BED)     │
│  ├ IGVF Catalog JSON dumps (genes/variants/...)    │
│  ├ IGVF KG Arango exports (nodes/edges JSONL)      │
│  ├ ENCODE files       (bigWig, peaks, Hi-C)        │
│  ├ FAVOR snapshots    (variant annotations)        │
│  ├ MaveDB scoresets   (VAMP-seq / MAVE CSVs)       │
│  └ Perturbation Catalogue (CRISPR / Perturb-seq)   │
└──────────────────────┬─────────────────────────────┘
                       │ skill-driven extract +
                       │ schema validation
                       ▼
┌────────────────────── SILVER ──────────────────────┐
│ DuckDB analytical warehouse (canonical entities)   │
│  ├ entities/                                       │
│  │   ├ genes        (~60k human + mouse)           │
│  │   ├ variants     (~10M FAVOR cumulative)        │
│  │   ├ proteins     (~80k UniProt)                 │
│  │   ├ regulatory_elements  (~1M SCREEN cCREs)     │
│  │   ├ cells        (per-experiment cell BCs)      │
│  │   ├ samples      (~10k biosamples)              │
│  │   └ datasets     (~5k MeasurementSets +         │
│  │                   1.2k Perturbation Cat. sets)  │
│  ├ edges/                                          │
│  │   ├ variant_to_gene         (eQTL / vEP)        │
│  │   ├ gene_to_protein         (translation)       │
│  │   ├ protein_to_protein      (BioGRID + IntAct + │
│  │   │                          HuRI + Reactome)   │
│  │   ├ gene_regulated_by       (rE2G / ABC)        │
│  │   ├ variant_in_re           (FAVOR cCRE class)  │
│  │   ├ cell_to_sample          (per-experiment)    │
│  │   └ dataset_provenance      (file → entity)     │
│  ├ measurements/                                   │
│  │   ├ vampseq_scores    (~50k MAVE variants)      │
│  │   ├ perturbseq_effects (~30M cell-gene-effect)  │
│  │   ├ crispri_calls     (~100k significant)       │
│  │   ├ mpra_activities   (~1M elements)            │
│  │   └ rnaseq_deg        (per-comparison DEG)      │
│  └ provenance/                                     │
│      ├ runs       (every skill invocation)         │
│      ├ versions   (upstream version per source)    │
│      └ qc         (per-table QC pass/fail)         │
└──────────────────────┬─────────────────────────────┘
                       │ embedding extraction +
                       │ contrastive alignment
                       ▼
┌────────────────────── GOLD ────────────────────────┐
│ Feature store + vector indices + model artefacts   │
│  ├ embeddings/                                     │
│  │   ├ gene_embeddings.parquet      (Geneformer)   │
│  │   ├ variant_embeddings.parquet   (Enformer)     │
│  │   ├ protein_embeddings.parquet   (ESM-2)        │
│  │   ├ re_embeddings.parquet        (ChromBPNet)   │
│  │   └ cell_embeddings.parquet      (scGPT)        │
│  ├ vector_indices/                                 │
│  │   ├ gene.faiss          (HNSW, IP / cosine)     │
│  │   ├ variant.faiss                               │
│  │   ├ protein.faiss                               │
│  │   └ cell.faiss                                  │
│  └ models/                                         │
│      ├ ckpts/igvf-foundation-v0.pt                 │
│      ├ tokenizer.json                              │
│      └ training_runs/<ts>_<label>/                 │
└────────────────────────────────────────────────────┘
```

---

## 2. Database choice — and why

### Primary warehouse → **DuckDB**

| Criterion | DuckDB | SQLite | Postgres | BigQuery |
|---|---|---|---|---|
| Embedded (no server) | ✓ | ✓ | — | — |
| Columnar analytics | ✓ | — | (extensions) | ✓ |
| Reads Parquet natively | ✓ | — | (FDW) | ✓ |
| Writes Parquet natively | ✓ | — | — | ✓ |
| 100M-row joins on laptop | ✓ | slow | ✓ | ✓ |
| Vector functions | ✓ (1.0+) | — | (pgvector) | — |
| Python-first | ✓ | ✓ | — | — |
| Free / local | ✓ | ✓ | ✓ | $ |

DuckDB is the right choice because it is **the same shape as the
existing IGVFagent stack** (embedded, file-backed, Python-first), but
it is **orders of magnitude faster than SQLite for the analytical
joins** that dominate this workload — variants joined to genes joined
to regulatory elements joined to perturbation effects. It reads and
writes Parquet natively, which means the Gold layer is just files —
no separate vector DB, no separate analytics DB.

### Graphs → keep **SQLite** for ≤ 10M edges, optional **Kùzu / Neo4j** for bigger

The existing `Data/Proteomics/KG/proteomics.sqlite` and
`Data/KG/portal_kg.sqlite` stay. When edge counts exceed ~10M and you
need multi-hop Cypher-style traversal, switch to:

- **Kùzu** — embedded property graph DB; pairs naturally with DuckDB
  (same query model, same files-on-disk philosophy)
- **ArangoDB** — what the IGVF Catalog upstream uses; matches their schema 1:1

### Vector indices → **FAISS** (`hnswlib` if you prefer pure-Python)

For nearest-neighbour search on the Gold-layer embeddings. Index size
~ 100 MB per million 768-d float32 vectors. FAISS-CPU is enough on a
128 GB-RAM laptop; FAISS-GPU only matters at >10M vectors.

### Object storage → **filesystem now, MinIO / S3 later**

Raw Bronze files. The existing `Data/<skill>/...` convention already
plays this role. When sharing across machines: `mc mirror` to MinIO,
no code change in the skills (they all use `Path` + URL).

---

## 3. Canonical schema (Silver layer)

```sql
-- ENTITIES ---------------------------------------------------------
CREATE TABLE genes (
    gene_id           TEXT PRIMARY KEY,            -- Ensembl, ENSG…
    symbol            TEXT,                        -- HGNC
    chrom             TEXT,
    start_bp          BIGINT,
    end_bp            BIGINT,
    strand            CHAR(1),
    biotype           TEXT,
    taxon             INT,
    refseq            TEXT,
    uniprot_acs       TEXT[],                      -- canonical proteins
    source_priority   TEXT,                        -- 'igvf-catalog'/'gencode'/...
    ingested_at       TIMESTAMP DEFAULT now()
);

CREATE TABLE variants (
    variant_id        TEXT PRIMARY KEY,            -- chrN_pos_ref_alt
    chrom             TEXT,
    pos_bp            BIGINT,
    ref               TEXT,
    alt               TEXT,
    rsid              TEXT,
    type              TEXT,                        -- SNV/INS/DEL/CNV
    population_af     DOUBLE,
    cadd              DOUBLE,
    favor_score       DOUBLE,
    consequence       TEXT,                        -- VEP top
    overlaps_re       TEXT[],                      -- cCRE class names
    taxon             INT,
    ingested_at       TIMESTAMP DEFAULT now()
);

CREATE TABLE proteins (
    uniprot_ac        TEXT PRIMARY KEY,
    symbol            TEXT,
    gene_id           TEXT,
    length_aa         INT,
    organism_taxid    INT,
    sequence_md5      TEXT,
    domains_jsonb     JSON,                        -- Pfam/InterPro
    pdb_ids           TEXT[],
    ingested_at       TIMESTAMP DEFAULT now()
);

CREATE TABLE regulatory_elements (
    re_id             TEXT PRIMARY KEY,            -- SCREEN cCRE or IGVF GE
    chrom             TEXT,
    start_bp          BIGINT,
    end_bp            BIGINT,
    classification    TEXT,                        -- PLS/pELS/dELS/CTCF
    assembly          TEXT,
    score             DOUBLE,
    n_assays          INT,
    source            TEXT,
    ingested_at       TIMESTAMP DEFAULT now()
);

CREATE TABLE cells (
    cell_id           TEXT PRIMARY KEY,            -- dataset:barcode
    dataset_id        TEXT,
    sample_id         TEXT,
    n_genes           INT,
    n_counts          INT,
    pct_counts_mt     DOUBLE,
    leiden            TEXT,
    cell_type_label   TEXT,
    multiseq_call     TEXT,                        -- singlet/multiplet/negative
    ingested_at       TIMESTAMP DEFAULT now()
);

CREATE TABLE samples (
    sample_id         TEXT PRIMARY KEY,
    biosample         TEXT,
    tissue            TEXT,
    cell_line         TEXT,
    sex               TEXT,
    age               TEXT,
    treatment         TEXT,
    disease           TEXT,
    donor_id          TEXT
);

CREATE TABLE datasets (
    dataset_id        TEXT PRIMARY KEY,            -- IGVF/ENCODE/MaveDB/PCat acc.
    upstream_source   TEXT NOT NULL,               -- 'igvf-portal','encode',
                                                   --   'mavedb','perturb-cat','geo'
    accession         TEXT,
    title             TEXT,
    assay             TEXT,
    modality          TEXT,                        -- 'scRNA','scATAC','VAMP-seq',
                                                   --   'CRISPRi','MPRA','bulk',...
    biosample_id      TEXT,
    year              INT,
    n_samples         INT,
    pubmed_id         TEXT,
    license           TEXT,
    download_url      TEXT,
    bronze_path       TEXT,                        -- local Bronze location
    ingested_at       TIMESTAMP DEFAULT now()
);

-- EDGES ------------------------------------------------------------
-- All edges are typed, sourced, and provenance-tracked.
CREATE TABLE edges (
    src_type     TEXT NOT NULL,                    -- 'gene'|'variant'|'protein'|…
    src_id       TEXT NOT NULL,
    dst_type     TEXT NOT NULL,
    dst_id       TEXT NOT NULL,
    relation     TEXT NOT NULL,                    -- 'regulated_by'|'eqtl'|
                                                   --   'binds'|'overlaps'|…
    score        DOUBLE,
    pubmed_id    TEXT,
    evidence     TEXT,                             -- 'experimental'|'predicted'|…
    upstream     TEXT NOT NULL,                    -- 'biogrid','igvf-rE2G',…
    UNIQUE(src_type, src_id, dst_type, dst_id, relation, upstream)
);
CREATE INDEX idx_edges_src  ON edges(src_type, src_id);
CREATE INDEX idx_edges_dst  ON edges(dst_type, dst_id);
CREATE INDEX idx_edges_rel  ON edges(relation);

-- MEASUREMENTS (one table per modality, all reference entity ids) -
CREATE TABLE vampseq_scores (
    measurement_id   BIGINT,
    gene_id          TEXT,
    protein_ac       TEXT,
    aa_position      INT,
    aa_wt            CHAR(1),
    aa_alt           CHAR(1),
    score            DOUBLE,
    se               DOUBLE,
    abundance_class  INT,
    dataset_id       TEXT,
    upstream         TEXT
);

CREATE TABLE perturbseq_effects (
    measurement_id   BIGINT,
    perturbed_gene   TEXT,
    target_gene      TEXT,
    log2fc           DOUBLE,
    padj             DOUBLE,
    cell_type        TEXT,
    dataset_id       TEXT,
    upstream         TEXT
);

-- PROVENANCE -------------------------------------------------------
CREATE TABLE runs (
    run_id        TEXT PRIMARY KEY,
    skill         TEXT,
    subcommand    TEXT,
    args_json     JSON,
    started_at    TIMESTAMP,
    ended_at      TIMESTAMP,
    rows_emitted  BIGINT,
    target_table  TEXT,
    success       BOOLEAN
);

CREATE TABLE versions (
    source        TEXT PRIMARY KEY,
    version       TEXT,
    upstream_url  TEXT,
    sha256        TEXT,
    fetched_at    TIMESTAMP
);

CREATE TABLE qc (
    table_name    TEXT,
    metric        TEXT,
    value         DOUBLE,
    passed        BOOLEAN,
    threshold     DOUBLE,
    run_id        TEXT
);
```

---

## 4. Each skill becomes a data producer

Every existing IGVFagent skill maps to **one or more producer functions**
that land QC'd rows in Silver. No new logic; just an ingestion contract.

| Skill | Producer → Silver table |
|---|---|
| `client` / `data` / `frontpage` | `datasets`, `samples` (Portal-side) |
| `kg gene` / `kg variant` / `kg region` | `genes`, `variants`, `edges{ variant→gene, gene↔re }` |
| `portal-kg` | bulk nodes/edges (already SQLite; mirror into DuckDB) |
| `variant` / `advanced-variant` | `variants`, `edges{ variant→gene }`, FAVOR `scores` |
| `ccre` | `regulatory_elements`, `edges{ variant↔re, gene↔re }` |
| `enhancer` | `edges{ gene_regulated_by_re }` |
| `mpra` | `measurements.mpra_activities` |
| `crispri` | `measurements.crispri_calls` |
| `encode` | `regulatory_elements`, `edges{ peak_overlaps_re }`, signal Parquet |
| `se-targets` | `regulatory_elements (SE)`, `edges{ se → target_gene }` |
| `singlecell` / `sc-analyze` | `cells`, `measurements.scrna_quant` |
| `multiome` / `splitseq` | `cells`, `samples`, multi-modal counts |
| `multiseq` | `cells.multiseq_call`, `edges{ cell↔sample }` |
| `proteomics` | `proteins`, `edges{ protein↔protein }`, `vampseq_scores` |
| `perturb-catalog` | `measurements.perturbseq_effects`, `datasets` |
| `geo` / `rnaseq` | `datasets`, `measurements.rnaseq_deg` |
| `ref` | `datasets.pubmed_id` enrichment |

**The contract**: every producer writes to a single named DuckDB table
via `INSERT OR REPLACE` keyed on the canonical primary key, and emits a
row in the `runs` provenance table. Any reader can join across all
modalities with one SQL statement.

---

## 5. QC, harmonization, and identifier mapping

The Silver-layer schema is meaningless if `gene_id` from IGVF Catalog
doesn't match `gene_id` from ENCODE. The harmonization steps:

1. **Genes → Ensembl gene IDs** (HGNC symbols are unstable). Use the
   IGVF Catalog's gene-id mapping table as the authoritative source.
   The proteomics `id_map` table already does this.
2. **Variants → `chrN_pos_ref_alt` on GRCh38**. Auto-liftover when
   upstream gives hg19. Track both rsid + chr-pos-ref-alt.
3. **Proteins → canonical UniProt accession**. Joint with isoform via
   a separate `protein_isoforms` table when needed.
4. **Regulatory elements → SCREEN cCRE ID when available, otherwise
   a deterministic hash of `(chrom, start, end, source)`**.
5. **Cells → `<dataset>:<barcode>`** (uniqueness across experiments).
6. **Datasets → upstream accession** (IGVFDS…, GSM…, urn:mavedb:…).

A `qc` row is emitted at every ingestion step. Failure modes that
must block a run:

- > 1 % rows fail schema validation
- > 5 % rows have a missing canonical ID after mapping
- Any duplicate primary keys (unless explicit upsert)
- Any cross-source disagreement on a fact already in the warehouse
  (gets a `divergence_flag` row in `qc`, not a hard block)

---

## 6. Embeddings — what to compute and with what

The Gold layer pre-computes embeddings for every entity using the
best available pre-trained models. All embeddings land in Parquet
files keyed by entity id, with a parallel FAISS index for kNN search.

| Entity | Embedding source | Dimension | Compute cost |
|---|---|---|---|
| **Gene** | [Geneformer](https://huggingface.co/ctheodoris/Geneformer) (`fine-tuned-public-30M`) | 512 | ~10 min CPU for 60k genes |
| | [scGPT](https://github.com/bowang-lab/scGPT) (when cell context is available) | 512 | ~30 min GPU |
| **Variant** | [Nucleotide Transformer 2.5B](https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-2.5b-multi-species) on ±1 kb window | 2560 | ~hours CPU for 10M variants — use GPU |
| | [Enformer](https://github.com/google-deepmind/deepmind-research/tree/master/enformer) for cell-type-specific signal | 5313 tracks | precomputed predictions |
| **Protein** | [ESM-2 (650M)](https://huggingface.co/facebook/esm2_t33_650M_UR50D) | 1280 | ~30 min CPU for 80k proteins |
| | [AlphaFold structure embedding](https://alphafold.ebi.ac.uk/) | structure | precomputed |
| **Regulatory element** | [ChromBPNet](https://github.com/kundajelab/chrombpnet) or DNABERT on the cCRE sequence | 768 | ~hours CPU |
| **Cell** | [scGPT](https://github.com/bowang-lab/scGPT) over the integrated UMI matrix | 512 | ~minutes GPU per dataset |

These are all **frozen feature extractors** in this phase. The next
phase trains a fusion model that aligns them.

---

## 6b. Network-integration layer — CORNETO

Embedding-based alignment (Phase 4 below) is one half of cross-source
integration; **constrained-optimization-based network inference** is
the other. We adopt **CORNETO** (Saez-Rodriguez lab,
[*Nat Mach Intell* 2025](https://www.nature.com/articles/s42256-025-01069-9))
as the framework for the second half. CORNETO reformulates many
separate biological-network methods (CARNIVAL, COSMOS, PCSF/Steiner,
OmniPath subnetwork extraction, FBA/iMAT, shortest-path) as instances
of one MILP over a signed directed graph:

- Binary edge/vertex activation indicators + continuous flow
- Flow conservation + sign consistency constraints
- Objective = gap-to-data + L0 sparsity (edges, vertices) + L1 flow

In IGVFagent terms, the integration flow is:

```
Silver edges (BioGRID + IntAct + HuRI + Reactome PPI) ─┐
Silver edges (rE2G + ABC + SE→target)                  ├─▶ signed PKN
Silver edges (FAVOR variant → RE overlap)              ─┘
                                                         │
   Silver measurements (Perturb-seq DEG, log2FC)         ▼
   Silver measurements (VAMP-seq abundance change)   CORNETO MILP
   Silver measurements (CRISPRi hits)                    │
                                                         ▼
                              context-specific subnetwork edges
                              tagged upstream='corneto:<method>:<label>'
                              and written back to Silver `edges` table
```

The inferred subnetworks become **first-class edges in the warehouse**,
queryable like any other source. They are *also* training pairs for
the Phase-4 contrastive alignment: a CORNETO-inferred edge between two
proteins under a condition is a positive pair in the latent space for
that condition.

**License**: CORNETO is GPL-3.0; IGVFagent stays Apache-2.0 by
calling it only at runtime as a `pip`-installed dependency, never
copying source. See ``Docs/Skills/CORNETO_INTEGRATION_SKILLS.md`` for
the wrapper skill's design.

The shipped wrapper exposes three sub-commands today:
- `corneto demo`      — synthetic cascade end-to-end self-test
- `corneto pkn-from-kg` — materialise signed PKN from the proteomics SQLite KG
- `corneto carnival`  — Perturb-seq / DEG → upstream subnetwork
- `corneto steiner`   — VAMP-seq prizes / GWAS hits → connecting subnetwork

## 7. Foundation-model architecture

The IGVF Foundation Model is a **multi-modal contrastive transformer**
that learns a joint embedding space across the five entity types
above by exploiting the known cross-modal relations in the edges
table.

### Phase 1 — modality-specific encoders (frozen, off-the-shelf)

```python
gene_emb     = Geneformer(gene_id)
variant_emb  = NucleotideTransformer(seq_context)
protein_emb  = ESM2(uniprot_seq)
re_emb       = ChromBPNet(re_seq)
cell_emb     = scGPT(cell_expression)
```

### Phase 2 — projection heads + contrastive alignment

```python
g_proj  = Linear(512 → 768)(gene_emb)
v_proj  = Linear(2560 → 768)(variant_emb)
p_proj  = Linear(1280 → 768)(protein_emb)
re_proj = Linear(768 → 768)(re_emb)
c_proj  = Linear(512 → 768)(cell_emb)
```

Train with **InfoNCE / CLIP-style contrastive loss** on the
cross-modal edges in the warehouse:

```
        variant ↔ gene    (eQTL, vEP, FAVOR)
           gene ↔ protein  (canonical translation)
        protein ↔ protein  (BioGRID + IntAct + HuRI + Reactome)
           gene ↔ re       (rE2G + ABC + SE→target)
        variant ↔ re       (FAVOR cCRE overlap)
           cell ↔ sample   (multiseq + sc-analyze clusters)
           cell ↔ dataset  (per-experiment provenance)
```

Each edge in the Silver `edges` table is a positive pair. Negatives
are mined in-batch. Loss:

```
L = Σ_(x_i, y_j ∈ edges) -log( exp(sim(p(x_i), p(y_j)) / τ)
                              / Σ_k exp(sim(p(x_i), p(y_k)) / τ) )
```

This is small enough to train on a laptop GPU (one epoch over 10M
edges with 768-d projections fits in ~12 GB VRAM at batch 1024).

### Phase 3 — generative / predictive head

After the projection heads are trained:

1. **Masked-token pre-training** on the joint graph: hide a random
   subset of entities, predict their embeddings from neighbours.
2. **Downstream fine-tuning** on canonical tasks:
   - Variant effect prediction (CADD / ClinVar)
   - Disease association (OMIM / GWAS catalog)
   - Drug-target prediction (DGIdb)
   - Cell-type annotation transfer
   - Per-gene perturbation effect prediction

### Phase 4 — release

Checkpoint to `Gold/models/igvf-foundation-v0.pt`, tokenizer to
`Gold/models/tokenizer.json`. Expose via:

```bash
igvfagent warehouse embed --entity-type variant --id chr17_43044295_G_A
# -> 768-d vector

igvfagent warehouse nn --entity-type gene --id BRCA1 --k 20
# -> top-20 most similar entities across all modalities

igvfagent warehouse predict --task variant_pathogenicity --input variants.csv
```

---

## 8. Implementation roadmap

### Phase 1 — Warehouse skeleton (week 1-2)

- New skill `warehouse` with subcommands:
  - `init` — create the DuckDB schema, set up `Data/Warehouse/`
  - `ingest --source <skill>` — pull each skill's outputs into Silver
  - `query "<SQL>"` — run an SQL query against the warehouse
  - `stats` — row counts per entity / edge / measurement table
- Wire `portal-kg`, `proteomics`, `perturb-catalog`, `multiseq` first
  (highest entity-id coverage)

### Phase 2 — Harmonization + QC (week 3-4)

- Implement `id_map` ingestion from UniProt + Ensembl + HGNC
- Add `qc` checks for every ingest
- Cross-source divergence flagging

### Phase 3 — Feature extraction (week 5-8)

- Ingest the off-the-shelf pre-trained models
- Compute and cache embeddings to Parquet
- Build FAISS indices

### Phase 4 — Contrastive alignment (week 9-12)

- Implement the contrastive trainer
- Train v0 IGVF Foundation Model on the laptop GPU
- Release checkpoint + tokenizer + minimal serving layer

### Phase 5 — Predictive heads + benchmarks (ongoing)

- Variant effect prediction
- Disease association
- Cell-type transfer
- Compare to existing single-modality SOTA on each task

---

## 9. Why this approach (not the alternatives)

- **Not a single graph database** because the embedding / vector
  side of the workload is columnar, not graph-shaped. Storing 10M
  variant embeddings as graph node properties is the wrong shape.
- **Not Postgres** because we don't need concurrent writes, and
  DuckDB is 10–100× faster on the analytical joins that dominate.
- **Not a cloud data warehouse (BigQuery / Snowflake)** because the
  whole point of IGVFagent is local, auditable, reproducible. A
  cloud warehouse would break that contract.
- **Not building our own foundation model from scratch** because
  pre-trained Geneformer / scGPT / ESM-2 / Enformer / Nucleotide
  Transformer are excellent baselines and freely available. Train
  only the alignment heads — that's where the IGVF-specific value
  lies.
- **Not a separate vector DB (Pinecone, Weaviate)** because FAISS on
  Parquet is enough at this scale and stays in the same local
  filesystem as everything else.

---

## 10. References for the model components

- **Geneformer** — Theodoris et al. *Nature* 2023. PMID 37258680.
- **scGPT** — Cui et al. *Nat Methods* 2024. PMID 38409223.
- **Nucleotide Transformer** — Dalla-Torre et al. *Nat Methods* 2024.
- **Enformer** — Avsec et al. *Nat Methods* 2021. PMID 34608324.
- **ESM-2** — Lin et al. *Science* 2023. PMID 36927031.
- **AlphaFold 2 / 3** — Jumper et al. *Nature* 2021 / Abramson 2024.
- **ChromBPNet** — Kundaje lab, github.com/kundajelab/chrombpnet.
- **CLIP** — Radford et al. ICML 2021 (the contrastive recipe).
- **DuckDB** — Raasveldt & Mühleisen, SIGMOD 2019; https://duckdb.org

---

This document lives at `Docs/Architecture/INTEGRATED_DATA_LAYER.md`.
A starter skill (`igvfagent warehouse`) accompanies it as
`Scripts/warehouse_skill.py`.
