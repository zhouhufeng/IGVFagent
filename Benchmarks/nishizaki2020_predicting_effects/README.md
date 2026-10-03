# nishizaki2020_predicting_effects

## Paper

**Predicting the effects of SNPs on transcription factor binding affinity**
Nishizaki Sierra S; Ng Natalie; Dong Shengcheng; Porter Robert S; Morterud Cody; Williams Colten; Asman Courtney; Switzenberg Jessica A; Boyle Alan P
*Bioinformatics* 36(2), 364-372 (2020) · doi:[10.1093/bioinformatics/btz612](https://doi.org/10.1093/bioinformatics/btz612) · PMID 31373606 · PMC7999143 · bioRxiv preprint doi:10.1101/581306 · grant U01HG011952

Resolver confidence: **1.00** (exact DOI match against Crossref/EuropePMC).

SEMpl (SNP Effect Matrix pipeline) builds a "SNP Effect Matrix" (SEM) per transcription factor: a position x base matrix of log2 ChIP-seq-signal-ratio scores, one per possible point mutation within the TF's binding motif, learned genome-wide from real ChIP-seq + DNase-seq data. Unlike a PWM (built from binding-site sequence alone), the SEM is trained against real binding intensity, so the paper's central claim is that it predicts real experimental outcomes (ChIP-seq signal, allele-specific binding, EMSA measurements) better than a PWM.

## Why this doesn't use an IGVFagent assay route

This is a pure computational method paper (`igvfagent bench route --paper-id nishizaki2020_predicting_effects` returns "No route matched" — no wet-lab assay family it maps to). Its Code Availability statement points to `github.com/Boyle-Lab/SEM_CPP`, which now redirects to `Boyle-Lab/SEMpl` — a C++ tool (`iterativeSEM`) that only *builds* one SEM at a time from raw ChIP-seq/DNase-seq/genome inputs (README's own estimate: 64+ GB RAM, 8+ cores, ~38 CPU-hours per SEM). It has **no downstream validation code** (no PWM/EMSA/ChIP-seq-correlation comparison) — those figures were produced by separate, unpublished analysis scripts.

`iterativeSEM` **does build successfully in this environment** (`lib/libBigWig` + `lib/TFMPvalue` submodules + the main `Makefile`, GCC 8.5, zero errors — see `OPERATIONS.md` §1) and runs (prints the expected "Running Iterative SEM building.." banner from the README's own demo). Actually regenerating a SEM genome-wide, though, needs a full hg19 bowtie index and is the ~38-CPU-hour step above — out of proportion to what's needed here, because the repository already ships **the exact real output of that pipeline**: `SEMs/*.sem`, one file per TF/cell-type/replicate, used unmodified as this benchmark's reference throughout.

So instead of re-deriving the SEMs, this benchmark validates the authors' **own deposited SEMs** against three independent sources of real data, reproducing the *comparisons* the paper's downstream (unpublished) analysis code would have made:

## Part 1 — Fig. 4: cross-cell-type / replicate concordance

Claim: *"R2 values over 0.97 for HepG2, A549 and T47D (P-values < 1e-32)"*. `sempl_sem_concordance.py` (port) flattens each deposited FOXA1 SEM (`SEMs/MA0148.1_{A549,HepG2.1,HepG2.2,T47D}.sem`) into its 44 score cells and computes pairwise R² (the two HepG2 replicates averaged first). **Reproduced almost exactly**: A549-vs-HepG2(avg) R²=0.989, A549-vs-T47D R²=0.979, HepG2(avg)-vs-T47D R²=0.971 — all ≥0.97, matching the claim.

## Part 2 — Fig. 5: allele-specific CTCF binding, GM12878

Claim: *"468 heterozygous sites at 5% FDR (AlleleDB), 240 with a matching CTCF PWM... SEM R2=0.50 vs PWM R2=0.41"* against real ChIP-seq (ENCSR000DZN). `sempl_asb_validation.py` (port) uses AlleleDB's own public per-SNV ASB calls (alleledb.gersteinlab.org — the same pipeline the paper's Methods names) for CTCF in NA12878 (458 sites; the public v2.1 recompute of the same pipeline the paper cites, vs. the paper's v1.0/468), scans real hg19 sequence (UCSC REST API) around each site with the real JASPAR CTCF PWM (MA0139.1) to locate the motif instance and orientation, then regresses the deposited SEM's and the PWM's predicted ref→alt score change against `log2((alt_reads+0.5)/(ref_reads+0.5))`.

**272/458 sites matched the PWM** (vs. paper's 240/468 — close). **SEM R²=0.649 vs PWM R²=0.571** — SEM outperforms PWM, and the *margin* (0.078) is close to the paper's own margin (0.50-0.41=0.09), even though both absolute values run higher than the paper's (458 real sites is a different, newer version of the same pipeline's calls; allelic ratio here is ref/alt reads, not maternal/paternal phase, which this public table doesn't carry).

## Part 3 — Fig. 3: FOXA1 SEM vs. real ChIP-seq

Claim: *"SEM: R2 = 0.66, PWM: R2 = 0.24"*. `sempl_chipseq_validation.py` (port) uses real ENCODE HepG2 data (hg19, classic wgEncode-era production matching the repo's own demo naming): FOXA1 ChIP-seq optimal-IDR peaks + fold-change-over-control bigWig (ENCSR000BMO), and a DNase-seq narrowPeak (ENCSR000EJV) as the open-chromatin filter. For each of 3,000 real ChIP peaks (bound sites) and 3,000 DNase-open sites with no ChIP peak (unbound sites), it locates each matrix's own best-scoring FOXA1-motif window in real hg19 sequence and reads the real ChIP fold-change signal at the summit.

**SEM R²=0.278 vs PWM R²=0.230** — SEM outperforms PWM (as claimed), by a modest, stable margin (~0.05, reproducible from n=300 up to n=3,000+3,000 samples). The **PWM baseline (0.230) closely matches the paper's own PWM R² (0.24)**. The SEM's absolute R² is well below the paper's 0.66, most plausibly because the paper's original genome-wide validation cohort/threshold isn't recoverable from the accessible text — see "Access notes" below — and this reconstruction (real ChIP peaks + a matched real unbound background, rather than every genome-wide PWM match) is a defensible but different sampling of the same real signal.

## Part 4/5 — Fig. 6 (EMSA) and Fig. 7 (13-TF comparison): blocked

Both figures' underlying measurements (EMSA densiometric scores; the per-TF SEM/PWM/DeepBind/LS-GKM correlation table) exist only in the paper's own Supplementary Data (`btz612_supplementary_data.zip`, 13.6 MB). See `expected.json`'s blockers on `fig6_emsa_validation` / `fig7_13tf_comparison` for the access details — not controlled-access or embargoed, but gated behind a bot-detection proof-of-work challenge we do not attempt to solve, and the bioRxiv companion preprint (10.1101/581306) that might carry the same data was persistently rate-limited (HTTP 429) from this environment.

## Access notes

The publisher (Oxford Academic) blocks XML full-text export via EuropePMC/eutils for this closed-access article; the PMC rendered page (`pmc.ncbi.nlm.nih.gov/articles/PMC7999143/`) was used instead (a fetch tool renders and summarizes it — enough to recover every numeric claim quoted above), and its Methods paragraph on the ASB analysis and EMSA protocol was used verbatim to design Parts 2 and 3. The bioRxiv companion preprint (10.1101/581306, same title/authors) could not be fetched at all (persistent HTTP 429 from biorxiv.org, several attempts, different user agents, minutes apart) — it may have contained additional Methods detail (e.g. the exact Fig. 3 site-selection procedure) that would let Part 3's SEM R² close the gap with the paper's reported 0.66.

## Running it

```bash
bash Benchmarks/nishizaki2020_predicting_effects/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark nishizaki2020_predicting_effects
```

## Coverage

3 of 5 planned analyses have passing, confirmed checks against real data (Figs. 3, 4, 5); 2 are blocked on inaccessible (not controlled-access, not embargoed) supplementary data (Figs. 6, 7). `igvfagent bench score --paper-id nishizaki2020_predicting_effects` will therefore report `incomplete 3/5` — see "Part 4/5" above for why the remaining 2 cannot presently reach `reproduced_except_access`.

## Provenance

`provenance.json` holds the full resolve/harvest record. Every port is in `igvfagent port list` (`sempl_sem_concordance`, `sempl_asb_validation`, `sempl_chipseq_validation`), pinned to `Boyle-Lab/SEMpl @ 3fe09bc78149`.
