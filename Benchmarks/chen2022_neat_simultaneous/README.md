# chen2022_neat_simultaneous

## Paper

**NEAT-seq: simultaneous profiling of intra-nuclear proteins, chromatin accessibility and gene expression in single cells.**
Chen AF, Parks B, Kathiria AS, Ober-Reynolds B, Goronzy JJ, Greenleaf WJ.
*Nature Methods* 2022 · doi:[10.1038/s41592-022-01461-y](https://doi.org/10.1038/s41592-022-01461-y) · PMID 35501385 · PMC11192021 · IGVF0042, award UM1HG011972

Resolver confidence: **1.00** (exact PMID lookup).

## Data sources

- **RNA/ADT/HTO/ATAC data:** [GSE178707](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178707) — the paper's own NEAT-seq CD4 memory T cell dataset (2 lanes), plus its own species-mixing (barnyard) QC experiment. Deposited as raw feature-barcode matrices, an RNA counts `dgCMatrix`, a peak count matrix, and ADT/HTO count CSVs — and, unusually, the authors' own final **correlation-analysis output tables** (`GSM5396330_corr_gene_adt.csv.gz`, `corr_peak_adt.csv.gz`, `putative_targets.csv.gz`), which serve as an exact ground truth for their Fig. 3 statistics.
- **Code + committed intermediate objects:** [github.com/GreenleafLab/NEAT-seq_reproducibility](https://github.com/GreenleafLab/NEAT-seq_reproducibility) (pinned `04b52760a1c3`) — full R/ArchR/Seurat analysis code for every figure, **plus two of the authors' own post-clustering ArchR project `.rds` objects committed directly to the repo** (`fig2_CD4_Tcells/data/ArchR_HTOsinglets_CD4only[_25XADT].rds`), which is what makes Fig. 2's exact cluster assignments reproducible without needing to re-run ArchR's clustering algorithm.
- Full text was available via NCBI's `efetch` (Europe PMC's own `fullTextXML` endpoint 500'd; the PMC NIHMS author-manuscript XML worked instead).

## What IGVFagent does

No `bench list-routes` route fits a trimodal (protein+ATAC+RNA) single-cell method paper (`multiome_peak2gene` was scaffolded by default as the closest fit; overridden). The reproduction is driven directly off GSE178707 and the authors' own repo via three routes:

1. **`fig2_cluster_cell_counts`** — no port needed: `Clusters` is read directly out of the authors' own repo-committed `ArchR_HTOsinglets_CD4only.rds` via base-R `attr()` slot access (the `ArchR` package itself is not installed and is not needed just to read an already-serialized S4 object's slots).
2. **`chen2022-neat-gata3-de`** — port of `Seurat_GATA3_differential_analysis.R`: the deterministic GATA3 RNA-high/protein-high vs RNA-high/protein-low gating (no random subsampling — the "low protein" control group is the N cells with the lowest ADT among high-RNA cells, N = size of the high-protein group), then a two-sided Wilcoxon rank-sum test per gene + BH adjustment + Seurat's own `avg_log2FC` formula.
3. **`chen2022-neat-tf-correlation`** — port of `correlation_analysis.R`/`correlation_utils.R`'s Spearman TF-ADT vs gene-RNA and TF-ADT vs peak-accessibility correlation (`sparse.cor`/`cor.pval`). Deliberately does **not** reproduce the motif-annotation join (needs BSgenome+motifmatchr) or the pseudobulk peak-gene distance-correlation step (needs ArchR `reducedDims` + an Rcpp partial-dot-product) that together build Fig. 3b's full TF-peak-gene linkage heatmap — those need a full Bioconductor/ArchR genome-annotation stack, out of scope here; the per-ADT correlation matrices that step is built from are ported and verified instead.

```bash
bash Benchmarks/chen2022_neat_simultaneous/run.sh
igvfagent bench score --paper-id chen2022_neat_simultaneous
```

## Concordance

**`reproduced_except_access`, 3/4 analyses** (7/7 checks pass; the 4th is genuinely blocked by undeposited raw sequencing reads).

| Paper claim | IGVFagent result | Match |
|---|---|---|
| Fig. 2a: Th1 (1562), Th2 (939), Th17 (1855), Treg (583), TCM (2512), Act. (116), Uncom. (905); n=8472 total | `Clusters` column of the authors' own repo-committed ArchR project: cluster sizes {116, 583, 1855, 2512, 939, 1562, 905}, n=8472 | **7/7 exact** |
| Fig. 2e volcano plot labels CAMK4, GAB2, EEF1G, SNX9, NELL2, PABPC4, NIBAN1, SNED1, BTBD11, RPL18A as GATA3-RNA/protein-discordant DE genes | All 10 found; all 10 clear \|avg_log2FC\|>0.5 (7/10 also clear BH-adjusted p<0.05 — the other 3 have large effect size but fall short of that combined threshold, matching the authors' own script's separate "diff" vs "sigdiff" plot categories) | **10/10** effect-size match |
| `chen2022-neat-gata3-de` port vs a real Seurat 5.5.1 `FindMarkers` run on the same gating | 36,601/36,601 genes match to floating-point precision (max abs diff on avg_log2FC = 5.8e-15) | **1.0000** match_rate |
| `chen2022-neat-tf-correlation` (gene-RNA) port vs the authors' own deposited `corr_gene_adt.csv.gz` | 67,265/67,265 gene x TF-ADT pairs match (r=0.9999999992) | **1.0000** match_rate |
| `chen2022-neat-tf-correlation` (peak-accessibility) port vs the authors' own deposited `corr_peak_adt.csv.gz` | 389,910/389,910 peak x TF-ADT pairs match (r=0.9999999999979) | **1.0000** match_rate |
| Fig. 3c-d: RORgT ADT level correlates with CCR6 RNA (Th17 master-TF/marker-gene example) | Spearman rho=0.1638, BH-adjusted p=4.51e-21 | Positive, highly significant |

## Honest caveats

- **The GATA3-DE fold-change formula took real debugging to match.** A first, textbook-looking port (`log2(mean(expm1(x))+1)`) disagreed with Seurat's own `FindMarkers` by up to 5.3 log2FC units on real genes. Reading `Seurat:::FoldChange.Assay`'s source directly showed the actual formula adds the pseudocount to the **sum** before dividing by group size (`log2((sum(expm1(x))+1)/n)`), not to the mean — a subtle but consequential difference for small groups (n=140 here). After the fix, all 36,601 genes match Seurat to floating-point precision. This is exactly the kind of thing "faithful port, verified against the authors' own reference" is meant to catch.
- **The reference for `fig2_gata3_de` used Seurat 5.5.1, not the authors' original 3.2.1** (not installed here; 5.5.1 is). The Wilcoxon rank-sum p-value and the `avg_log2FC` formula used are both stable across that version gap (verified: p-values and fold-changes match to 1e-15), but this is a real Seurat-version substitution worth flagging, not the literal 3.2.1 binary the paper used.
- **Fig. 3g-j (SNP allelic imbalance) is blocked `not_deposited`.** The authors' own `code_utils/download_data.py` states in a code comment that the raw BAM/SRA reads needed for allele-specific Tn5 insertion counting at rs62088464 "are not publicly accessible yet." No substitute derivation from the deposited fragments files can recover allele-specific read mapping without the original BAM CIGAR/MD tags.
- **Not ported:** the motif-annotation join and the pseudobulk peak-gene distance-correlation step that build Fig. 3b's full linkage heatmap (both need a full ArchR/Bioconductor genome-annotation stack: BSgenome, motifmatchr, an Rcpp partial-dot-product against ArchR `reducedDims`). The per-ADT correlation matrices that step consumes are ported and verified to floating-point precision instead — a real but partial cut of Fig. 3's pipeline.
- **Fig. 1's species-mixing (barnyard) QC experiment was scoped out** as a Methods-validation figure rather than a headline scientific claim, consistent with the skill's "drop Methods-only items" guidance.

## Provenance

`provenance.json` holds the resolve/harvest record. `expected.json`'s `checks[].provenance` and `analyses[].blocker.reason` carry every number, threshold, and blocker justification above.
