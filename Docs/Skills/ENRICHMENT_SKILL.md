# Skill: GO + Pathway enrichment (validation layer)

Validate a gene list — DEGs from a multiome / RNA-seq run, top hits from
a CRISPR / Perturb-seq screen, or the gene-side of an enhancer-gene
linkage — by asking which **Gene Ontology** terms and which **canonical
pathway databases** (Reactome, KEGG, WikiPathways, MSigDB Hallmark) are
over-represented relative to a background. Two statistical modes:

* **ORA** — hypergeometric over-representation of a discrete hit list.
* **GSEA / preranked** — Subramanian-style enrichment of a *ranked*
  list (e.g. by t-statistic), no arbitrary cutoff.

Backend: [gseapy](https://github.com/zqfang/GSEApy) (BSD-3, Fang 2023).
ORA queries the Enrichr proxy; GSEA runs offline against cached .gmt
files.

## Commands

```bash
# ORA on every default library (GO_BP/MF/CC + Reactome + KEGG +
# WikiPathways + MSigDB_Hallmark) from a one-gene-per-line file:
igvfagent enrich ora --genes Data/MyHits/my_degs.txt --label my_degs_v1

# ORA restricted to GO branches only
igvfagent enrich go --genes Data/MyHits/my_degs.txt --label my_degs_go

# ORA restricted to canonical pathway DBs
igvfagent enrich pathways --genes Data/MyHits/my_degs.txt --label my_degs_paths

# Preranked GSEA from a TSV with 'gene' + 'score' columns
igvfagent enrich gsea --ranked Data/MyHits/my_ranked.tsv --label my_gsea

# End-to-end demo (cell-cycle gene set, expected strong CC enrichment)
igvfagent enrich showcase
```

### Common arguments

| Flag | Meaning |
|---|---|
| `--libs` | `all` / `go` / `pathways` / comma-list of friendly names (GO_BP, GO_MF, GO_CC, Reactome, KEGG, WikiPathways, MSigDB_Hallmark) or raw Enrichr lib ids |
| `--organism` | `human` / `mouse` / `rat` / `fly` / `worm` / `fish` / `yeast` |
| `--background` | Optional background gene list for ORA |
| `--label` | Output subdir tag |
| `--top-k` | Top-K terms per library in the composite figure (default 8) |

### GSEA-only arguments

| Flag | Meaning |
|---|---|
| `--ranked` | TSV/CSV with 'gene' + one of 'score','stat','log2fc','rank' |
| `--min-size` / `--max-size` | Gene-set size filter (default 10–1000) |
| `--permutations` | Permutation count (default 1000) |

## Input formats

### Gene list (ORA)

Either one-gene-per-line:

```
TP53
MYC
CDKN2A
...
```

Or CSV/TSV with a `gene` (or `symbol` / `gene_symbol` / …) column:

```
gene,log2fc,padj
TP53,2.4,1.3e-9
MYC,1.8,4.5e-7
```

### Ranked gene-score (GSEA)

TSV/CSV with a `gene` column and a `score` (or `stat` / `log2fc` /
`rank` / `t` / `wald` / `signed_lp`) column. Sign-aware — positive
scores → "up in condition A", negative → "up in B". Sort order is
recomputed internally; you don't have to pre-sort.

## Output schema (ORA `enrichment.tsv`)

| Column | Meaning |
|---|---|
| library | Friendly library name (GO_BP, Reactome, …) |
| library_id | Enrichr library id (`GO_Biological_Process_2023`, …) |
| term | Term name (`Cell Cycle Mitotic R-HSA-69278`, …) |
| overlap | `n_overlap/n_term_total` (Enrichr convention) |
| p_value | Fisher's exact P |
| adjusted_p_value | Benjamini–Hochberg adj-P within the library |
| odds_ratio | Fisher odds ratio |
| combined_score | `-log(P) * z-score-of-rank` (Enrichr ranking metric) |
| genes | Semicolon-joined overlap genes |

## Output schema (GSEA `gsea.tsv`)

| Column | Meaning |
|---|---|
| library / library_id / term | Same as above |
| es | Enrichment score |
| nes | Normalized enrichment score |
| p_value | Nominal permutation P |
| fdr | FDR q-value |
| size | Gene-set size after filtering |
| lead_genes | Leading-edge subset |

## How it works

1. **ORA path**: gseapy.enrichr posts the gene list to the Enrichr API
   (`https://maayanlab.cloud/Enrichr`), reads back hypergeometric
   stats, and writes per-library TSVs. We concatenate them and apply
   our composite-figure renderer (one bar chart per library + a
   bubble overview).
2. **GSEA path**: gseapy.prerank reads the ranked scores, walks every
   gene set in the chosen libraries (.gmt downloaded once and cached),
   computes weighted Kolmogorov–Smirnov ES, builds a null via gene-
   label permutation (1000 by default), and returns NES + FDR q.

## Caching

Enrichr responses are cached by gseapy under each run's `outdir`. The
.gmt files for GSEA libraries are cached under
`Data/Cache/Enrichment/.gseapy/<lib>.gmt` (controlled by the
`GSEAPY_CACHE` env var if set).

## Showcase test

`igvfagent enrich showcase` runs ORA on a curated 47-gene cell-cycle
list (CCN*, CDK*, CDKN*, MCM*, AURK*, PLK*, BUB*, etc.) and asserts
that:

* the top Reactome term is in the **Cell Cycle** family,
* the top KEGG term is **Cell cycle**,
* the top MSigDB Hallmark term is **G2-M Checkpoint** or **E2F
  Targets**, and
* all three GO branches show a mitotic / chromosome / kinase enrichment.

This serves as a positive-control validation that the skill itself is
healthy before it's used to validate downstream results.

## License posture

Apache-2.0. Heavy deps imported lazily:

* gseapy — BSD-3 (Z. Fang, A. Liu, M. Tu, *Bioinformatics* 2023)
* pandas — BSD-3
* matplotlib — PSF-style

Gene-set libraries:

* GO terms — Ashburner 2000 (Gene Ontology Consortium), CC-BY 4.0
* Reactome — Fabregat 2018, CC-BY 4.0
* WikiPathways — Slenter 2018, CC0
* MSigDB Hallmark — Liberzon 2015, CC-BY 4.0
* KEGG — accessed via Enrichr's academic-use proxy (Kanehisa 2000)

No GPL runtime deps.
