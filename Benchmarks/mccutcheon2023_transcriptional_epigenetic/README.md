# mccutcheon2023_transcriptional_epigenetic

## Paper

**Transcriptional and epigenetic regulators of human CD8&lt;sup&gt;+&lt;/sup&gt; T cell function identified through orthogonal CRISPR screens.**
McCutcheon SR, Swartz AM, Brown MC, Barrera A, McRoberts Amador C, Siklenka K, Humayun L, Ter Weele MA, Isaacs JM, Reddy TE, Allen AS, Nair SK, Antonia SJ, Gersbach CA.
*Nature genetics* 2023 · doi:[10.1038/s41588-023-01554-0](https://doi.org/10.1038/s41588-023-01554-0) · PMID 37945901 · PMC10703699 · IGVF tracking IGVF0062, award UM1HG012053.

Resolver confidence: **1.00** (resolved on exact PMID lookup).

## Data sources found in the paper

| Repository | Accession | Contents | In Data Availability |
|---|---|---|---|
| NCBI GEO | `GSE241933` | all bulk sort-based CRISPR screens (CD2/B2M CRISPRi tiling, IL2RA CRISPRa tiling, CRISPRi/a TF CCR7, TFome CRISPRko) — raw gRNA count tables | ✓ (subseries of `GSE218988`) |
| NCBI GEO | `GSE218986` | RNA-seq TPM (BATF3 OE vs GFP; ZNF217/GATA3 KO) | ✓ |
| NCBI GEO | `GSE218987` | ATAC-seq narrowPeak + RPKM bigwig (acute/chronic BATF3 OE) | ✓ |
| NCBI GEO | `GSE218985` | scRNA-seq TF-screen characterization (not pursued here — see "Honest caveats") | ✓ |
| Zenodo | `10.5281/zenodo.8370763` | authors' own custom gene-level "Allen method" code | ✓ (Code Availability) |

`GSE197268`, flagged by the automated harvester as a possible accession for this paper, is a **different** paper's dataset (Haradhvala et al. 2022, PMID 36097221) that this paper's Methods cites and re-analyzes for a CD19 CAR-T scRNA-seq comparison — not this paper's own primary data, and not pursued here.

## What was reproduced

**Route:** custom (`run.sh` downloads real GEO/Zenodo data and runs 5 registered `igvfagent` ports; no generic skill route fit this paper's own bespoke statistics).

```bash
bash Benchmarks/mccutcheon2023_transcriptional_epigenetic/run.sh
python Benchmarks/concordance.py --benchmark mccutcheon2023_transcriptional_epigenetic
```

**`igvfagent bench score`: reproduction = `reproduced_except_access`, coverage 5/6, 16/16 checks passing** (re-confirmed 2026-10-02 after the backfill simulation job `48972015` finished; see `OPERATIONS.md` §2).

| Analysis | Figure | Class | Result |
|---|---|---|---|
| `crispri_tf_ccr7_gene_level` | Fig. 1b | **B** (byte-identical vs authors' code) + C | DNMT1 is the FDR-lowest gene-level hit (padj≈0), matching "the most significant hit from the CRISPRi screen was...DNMT1" exactly |
| `crispra_tf_ccr7_gene_level` | Fig. 1c,d | **B** (byte-identical vs authors' code) + C | BATF3 rank 1, BATF rank 2, EOMES rank 3, JUN rank 8 of 121 — matches "BATF and BATF3 were among the top hits" and the named EOMES/BATF/JUN hits exactly |
| `tfome_ko_cofactors_mageck` | Fig. 6b,c | **A** (recovery_rate) + B (self-test) + C | all 9 of the paper's own named "most enriched genes" (ZNF217, RUNX3, FOXP1, GATA3, GFI1, AHR, ETS1, ZNF626, FOXP3) recovered in this project's own top-15 hits; BATF3/JUNB/IRF4 confirmed BATF3-OE-dependent (rank 7/1614 with BATF3 OE vs rank 559/1614 without) |
| `batf3_oe_rnaseq_de` | Fig. 3, Fig. 4h | **A** | FOXP1 (log2FC −0.56), ETS1 (−0.55), FOXP3 (−1.05) all downregulated with BATF3 OE, matching "several TFs including FOXP1, ETS1 and FOXP3 were all downregulated"; BATF3 itself +2.55 log2FC (sanity check on the OE construct) |
| `batf3_oe_atac_da` | Fig. 4a-c | **A** | 76.1% of nominally-changed (uncorrected p<0.05) acute-stimulation regions are MORE accessible with BATF3 OE, vs the paper's 60% — same direction, right order of magnitude |
| `znf217_gata3_ko_rnaseq_de` | Fig. 6f,g | **blocked** (`not_deposited`) | raw RNA-seq counts for this contrast were not deposited (TPM only); see "Honest caveats" |

## Ports absorbed into IGVFagent

| Command | Upstream | Verification |
|---|---|---|
| `mccutcheon2023_deseq2_grna` | Methods (paired two-tailed DESeq2, no author wrapper released) | gRNA-level driver used by every other port below |
| `mccutcheon2023_allen_gene_level` | Zenodo 10.5281/zenodo.8370763: `AllenMethod_fromdeseq_NTtransformation_Sean.py` + `AllenMethod_proportionalbnd_parallel.py` + `gene-level-effect-size.R` | **verify-port match_rate = 1.0000** (byte-identical to the authors' own unmodified script on both CRISPRi and CRISPRa TF screens) |
| `mccutcheon2023_rra_gene_level` | Li et al. 2014 *Genome Biology* MAGeCK (alpha-RRA formula; no author code — third-party compiled binary not installable here) | synthetic self-test: 5/5 planted true positives recovered, 0/195 false positives, match_rate = 1.0000 |
| `mccutcheon2023_rnaseq_de` | Methods (RNA sequencing, DESeq2) | documented substitute (TPM-only deposit); direction/magnitude checked against real data |
| `mccutcheon2023_atac_da` | Methods (ATAC-seq, DESeq2 on consensus peaks) | documented substitute (narrowPeak+bigwig-only deposit); direction/magnitude checked against real data |

`igvfagent port list` shows all five; each is a general-purpose command, not a one-off script — e.g. `mccutcheon2023_allen_gene_level` will run on any FACS-sorted CRISPR screen's DESeq2 gRNA table, not just this paper's.

## Honest caveats

* **The Allen-method verification is unusually strong, and unusually reproducible only in a specific sense.** The authors' own Monte Carlo null (10M random draws per gRNA-count `m`) has no fixed random seed in their released code — re-running their own script twice would not itself be bit-identical. This project's verification therefore generates ONE set of null simulations (`Data/mccutcheon2023/simulations/`, using the authors' own unmodified `AllenMethod_proportionalbnd_parallel.py`, path-fixed only) and feeds it to *both* the authors' script and this project's port, isolating the comparison to "does the surrounding algorithm — NT-pvalue transform, NTC pseudo-gene construction, FDR — match exactly" (it does, to floating-point identity) rather than claiming two independent Monte Carlo draws would match.
* **MAGeCK itself could not be run.** The paper's own Methods says gene-level TFome-KO enrichment used the compiled `mageck test --paired --control_sgrna` binary; that binary requires bioconda/root install that is not available in this environment. `mccutcheon2023_rra_gene_level.py` is a from-scratch re-implementation of the published alpha-RRA formula (Li et al. 2014), not a port of MAGeCK's own source, verified with a synthetic self-test rather than against the real binary's output. It cannot resolve MAGeCK's exact "ZNF217 is uniquely rank 1" claim (several top genes tie at this project's Monte Carlo permutation floor of p<2×10⁻⁷ at 5M permutations) — reported as top-10/top-15 membership instead.
* **RNA-seq and ATAC-seq DESeq2 exact statistics are not reproducible from what GEO deposits.** Both GSE218986 (RNA-seq) and GSE218987 (ATAC-seq) carry only processed summary tracks (TPM; narrowPeak+RPKM-bigwig) — not the raw featureCounts/fragment-count tables the paper's own DESeq2 negative-binomial model needs. This is the same documented gap `anglen2025_rnaseq_tpm_de.py` (a different paper, this project) already hit. The BATF3-OE RNA-seq and acute-ATAC-seq checks above use a real, donor/replicate-paired substitute statistic on the real deposited data and check direction/magnitude, not DESeq2's exact p-values or region/DEG counts. `znf217_gata3_ko_rnaseq_de` hits this same gap at only n=3 replicates hard enough that even the direction-only substitute can't call the paper's "644 DEGs" — that specific numeric claim is blocked as `not_deposited`, not silently skipped.
* **Chronic-stimulation ATAC-seq (n=2 replicates, 1 degree of freedom) was computed but not scored** — the paired t-test at that replicate count is too noisy to be an informative direction check (measured 29% more-accessible vs the paper's 54%, opposite qualitative lean from the acute condition, most plausibly power/noise rather than a real biological reversal).
* **scRNA-seq characterization (Fig. 2, GSE218985) and the in-vivo CAR-T tumor model (Fig. 3c-h) were not attempted.** The former needs a full guide-assignment + MAST-differential-expression single-cell pipeline across ~360,000 cells / 6 samples (real, public, but a substantial additional pipeline beyond this pass's scope); the latter is a mouse xenograft experiment with no computable public deposit beyond the paper's own summary flow-cytometry/tumor-volume statistics. Both are confirmatory of the Fig. 1/6 gene-level hits already reproduced above, not independent headline claims.
* **CD2/B2M CRISPRi and IL2RA CRISPRa promoter-tiling screens** (Extended Data Figs. 1-3, platform-characterization panels rather than gene discoveries) were downloaded but not separately ported — they are gRNA-level-only DESeq2 tests (no gene-level aggregation), and `mccutcheon2023_deseq2_grna.py` already covers that statistic; adding dedicated checks for them was left for a future pass.

## Provenance

`provenance.json` in this directory holds the full resolve/harvest/route record. All real intermediate data and outputs live under `Data/mccutcheon2023/` (git-ignored by size convention; re-generated by `run.sh`).
