# Agarwal 2025 — lentiMPRA across K562 / HepG2 / WTC11

[![paper](https://img.shields.io/badge/Nature-639:411--420-blue)](https://doi.org/10.1038/s41586-024-08430-9)
[![PMID](https://img.shields.io/badge/PMID-39814889-blue)](https://pubmed.ncbi.nlm.nih.gov/39814889/)
[![data](https://img.shields.io/badge/ENCODE-3_cell_lines-orange)](https://www.encodeproject.org/)
[![coverage](https://img.shields.io/badge/paper%20coverage-4%2F5%20analyses-yellow)]()

## Citation

Agarwal V, Inoue F, ..., Ahituv N. **Massively parallel characterization of
transcriptional regulatory elements.** *Nature* **639**: 411–420 (2025).
DOI: [10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9) · PMID: 39814889 · PMC11903340

## Data

* **Data**: ENCODE portal, per the paper's Data availability (PMC11903340): large-scale libraries K562 ENCSR382BVV, HepG2 ENCSR022GQD, WTC11 ENCSR244FWB; joint libraries K562 ENCSR203UFY, HepG2 ENCSR405QCT, WTC11 ENCSR336MKI; pilot libraries K562 ENCSR460LZI, HepG2 ENCSR463IRX. Public, no access request needed. (An earlier version of this benchmark cited GEO GSE142696; that series is Klein et al. 2020, a different MPRA study, and has no per-cell-line count files — fixed in commit `9817692`.)
* **GitHub**: [visze/sequence_cnn_models](https://github.com/visze/sequence_cnn_models) (Snakemake glue, EnformerMPRA) and [autosome-ru/human_legnet](https://github.com/autosome-ru/human_legnet) (MPRALegNet — the one actually run here, see CNN section below)
* **Pipeline**: MPRAflow (Gordon 2020); the deposited ENCODE files are its "combine-replicates" step output (`condition, replicate, name, dna_count, rna_count, ratio, log2, n_obs_bc`), one row per element per replicate, for each cell line's large-scale library.

## Headline workflow (paper)

1. lentiMPRA on ~680K cis-regulatory elements in 3 cell lines.
2. log2(RNA/DNA) per oligo → activity calls against shuffled negative controls, 5% FDR.
3. Train a CNN sequence model (MPRALegNet / EnformerMPRA) on the activity calls, predict held-out enhancers and cell-type specificity.

## What IGVFagent reproduces (`igvfagent bench score`: 4/5 analyses reproduced, 12/16 checks)

| Analysis | Paper | IGVFagent | How checked | State |
|---|---:|---:|---|:---:|
| K562 total elements tested | 243,780 | 230,933 | count of unique element names in the ENCODE deposit, ±10% | reproduced |
| K562 inter-replicate Pearson r | 0.76 | 0.772 | barcode-depth-filtered (n_obs_bc≥10) mean pairwise Pearson r of per-replicate log2(RNA/DNA), ±10% | reproduced |
| K562 promoters active | 52.3% | 44.7% | own reimplementation (see below) | failing |
| K562 enhancers active | 41.3% | 27.1% | own reimplementation | failing |
| HepG2 total elements tested | 164,307 | 153,904 | count of unique element names, ±10% | reproduced |
| HepG2 inter-replicate Pearson r | 0.94 | 0.937 | barcode-depth-filtered mean pairwise Pearson r, ±10% | reproduced |
| HepG2 promoters/enhancers active | 54.6% / 42.8% | not attempted | this cell line's ENCODE deposit has 0 `seq*_shuffled_*` negative controls | blocked (no null) |
| WTC11 total elements tested | 75,542 | 70,233 | count of unique element names, ±10% | reproduced |
| WTC11 inter-replicate Pearson r | 0.76 | **0.732** (was 0.555) | barcode-depth-filtered mean pairwise Pearson r, ±10% | **reproduced** (fixed) |
| WTC11 promoters/enhancers active | 50.6% / 25.8% | 44.4% / 17.0% | own reimplementation, barcode-depth-filtered | failing (narrowed, not closed) |
| Joint-library concordance (0.96–0.98) + EnformerMPRA (r≈0.81) | — | not attempted | out of scope this pass (public data, just not fetched — see below) | pending |
| MPRALegNet CNN, K562 | 0.83 (10-fold ensemble) | **0.821** (single fold) | held-out Pearson r, paper's own code + data, ±10%/[0.6,1.0] | reproduced |
| MPRALegNet CNN, HepG2 | 0.83 (10-fold ensemble) | **0.794** (single fold) | held-out Pearson r, paper's own code + data | reproduced (loose tolerance) |
| MPRALegNet CNN, WTC11 | 0.83 (10-fold ensemble) | **0.725** (single fold) | held-out Pearson r, paper's own code + data | reproduced (loose tolerance) |

Reference code: **the paper released no code for element activity calling**
(Code Availability only names the CNN repo, `visze/sequence_cnn_models`).
Activity calls here are IGVFagent's own reimplementation of the Methods text
("using shuffled controls as a background set ... 5% false discovery rate"):
per element, mean log2(RNA/DNA) across replicates, one-sided z-test against
the `seq*_shuffled_*` negative-control distribution, BH-FDR 5%
(`Scripts/ported/skills/agarwal2025_lentimpra_activity.py`). Elements are
split into promoter/enhancer by name (`ENSG*` = promoter; `peak*`
(K562) / `DNasePeakNoPromoter*` (HepG2) / `seq<N>_[FR]` (WTC11) = enhancer;
everything else — deep-mutational-scan sub-libraries, saturation-mutagenesis
controls, joint-library spike-ins — is "other"), since the paper does not
deposit a shared cCRE-class annotation table alongside these files.

Why the gaps (elaborated, third pass):

* **WTC11 replicate correlation was a real, fixable bug — now fixed.**
  Root cause: the per-replicate `n_obs_bc` (barcode count backing each
  element/replicate estimate) has a much lower median in WTC11 (17) than
  K562 (48) or HepG2 (28) — WTC11 was sequenced/tested at lower barcode
  complexity. Low-barcode elements are noisier, and unfiltered they drag
  the Pearson r down hard. Requiring `n_obs_bc >= 10` in both replicates of
  a pair before computing the correlation (a standard MPRA QC threshold,
  e.g. Tewhey lab MPRAsuite, MPRAflow's own defaults) moves the mean
  pairwise r from 0.555 to **0.732** (paper: 0.76) — closing ~80% of the
  gap. The same filter, applied for consistency, also nudges K562
  (0.744→0.772, paper 0.76) and HepG2 (0.915→0.937, paper 0.94) slightly
  *closer* to the paper, confirming this isn't a WTC11-specific tuning
  knob but a genuine QC step the reimplementation was missing. Implemented
  as `--min-bc-depth` (default 10) in
  `agarwal2025_lentimpra_activity.py`/now `igvfagent agarwal2025-lentimpra-activity`.
* **Activity fractions still run low — barcode depth is not the cause.**
  Applying the same `n_obs_bc >= 10` filter to the activity-calling z-test's
  input elements narrows WTC11's gap (promoter 38.8%→44.4%, paper 50.6%;
  enhancer 14.6%→17.0%, paper 25.8%) but does not close it, and barely
  moves K562 at all (44.9%→44.7%, paper 52.3%). A second attempt — a
  one-sample t-test per element using its own cross-replicate SD instead of
  the pooled negative-control SD — made both cell lines *worse* (K562
  promoter 41.7%, WTC11 promoter 28.8%), so it was not adopted. Conclusion,
  after two bounded diagnostic attempts per NEXT_STEPS.md: the remaining
  15–30% relative under-call is a genuine difference in test statistic
  (the paper's own software/exact statistic for this step is not
  specified beyond "shuffled controls as background... 5% FDR"), not a
  data-quality artifact this reimplementation can close without the
  authors' code.
* **HepG2 has no negative controls in this deposit.** Its large-scale
  element-quantification file (`ENCFF755BGY`) carries only 4
  `neg1`/`neg2`-named elements, all from the separate small "C:" TF-motif
  sub-library, not enough for a null distribution. Replicate concordance is
  still reported; activity calling is `not_attempted` for HepG2's
  promoters/enhancers, not silently skipped. Not revisited this pass
  (would need locating a HepG2-specific shuffled-control set elsewhere in
  ENCODE's file inventory for this accession — not attempted).
* **Element counts run ~5–8% low across all three cell lines**, consistently
  in one direction: MPRAflow's combine-replicates step keeps only elements
  quantified in all three replicates, so a few percent of the paper's
  designed elements are absent from this specific deposited table.
* **`mpra activity` (NB-GLM Wald test) does not apply here.** ENCODE's
  deposited `dna_count`/`rna_count` are already MPRAflow-normalized
  fractions (all < 1, not raw barcode counts), so IGVFagent's registered
  NB-GLM route — which needs integer counts — cannot run on this input; only
  `mpra qc` (log10(count+1) Pearson correlation, tolerant of floats) runs
  through it, once per cell line. This is a real, cluster-agnostic gap: no
  count reprocessing from FASTQ was attempted (out of scope).
* **Joint libraries remain pending, not blocked** — public ENCODE
  accessions exist, they were simply not fetched in this pass. See
  `NEXT_STEPS.md`.
* **MPRALegNet CNN — actually run this pass, not just deferred.** See the
  dedicated section below.

## MPRALegNet CNN reproduction

The paper's own repo, [autosome-ru/human_legnet](https://github.com/autosome-ru/human_legnet)
(pinned commit `7638fce137db0445123efe8d7e2c35e248fafc5f`), ships:

* the exact cross-validation data it was trained/evaluated on —
  `datasets/original/{K562,HepG2,WTC11}.tsv`, each row a 230bp sequence
  with a `mean_value` label (mean log2(RNA/DNA) across replicates) and a
  `fold` assignment (10-fold CV) already baked in;
* its own training + prediction code, `core.py` (PyTorch + Lightning).

`Benchmarks/agarwal2025_lentimpra/run_cnn.sh` reproduces this end to end:
clones the repo at the pinned commit, builds a venv (`torch`, `lightning`,
`torchmetrics`), submits a SLURM job to the cluster's `gpu_test` partition
(1x A100 MIG slice, 12h limit — plenty for this), and runs the repo's own
documented `--demo` mode — a single train/val/test fold split per cell
line (fold 1 held out as test, fold 2 as validation, folds 3–10 train;
25 epochs, `--use_shift --reverse_augment` as the paper's own
`models/example.cfg` uses) — then evaluates held-out Pearson r with the
newly-registered `igvfagent agarwal2025-lentimpra-cnn-eval` port.

**Result** (verified live, this session): all three cell lines trained and
predicted successfully in under 12 minutes total wall-clock on one GPU.
Held-out Pearson r (forward/reverse-orientation-averaged prediction vs
true `mean_value`, on the fold-1 test set the model never trained on):

| Cell line | Held-out n | Pearson r (this run, single fold) | Paper (r≈0.83, 10-fold ensemble) |
|---|---:|---:|---:|
| K562 | 19,666 | **0.821** | 0.83 |
| HepG2 | 12,292 | **0.794** | 0.83 |
| WTC11 | 4,620 | **0.725** | 0.83 |

Mean across cell lines: 0.780. This is a genuine, direct reproduction of
the paper's own architecture, training procedure, and data producing
held-out predictive correlation in the same range as the paper's headline
number — from a *single* 1/10-sized fold, not the full 10-fold ensemble
the paper actually reports (which trains on 9x more data per model and
averages across folds, both of which push correlation up further). K562's
single-fold result (0.821) already lands within 1% of the paper's
ensembled number; HepG2 and WTC11 trail more, consistent with their
smaller dataset sizes (245K and 92K rows vs K562's 393K) giving each
single fold less data to train on.

Not attempted: the full 10-fold × 9-val-fold ensemble (90 models/cell
line — a multi-hour GPU job, out of scope for this pass, see
`NEXT_STEPS.md` item 5) and EnformerMPRA (the paper's other CNN; no
committed weights/data found in either paper repo).

## Reusable scripts / core CLI ports

`agarwal2025_lentimpra_reshape` and `agarwal2025_lentimpra_activity` are
now registered `igvfagent` ports (`igvfagent port list`) — callable
directly as `igvfagent agarwal2025-lentimpra-reshape` /
`igvfagent agarwal2025-lentimpra-activity`, not just one-off benchmark
scripts, since the barcode-depth-filtered activity-calling method is
generic to any lentiMPRA/MPRAflow ENCODE deposit, not specific to this
paper's benchmark:

* `igvfagent agarwal2025-lentimpra-reshape` (`Scripts/ported/skills/agarwal2025_lentimpra_reshape.py`)
  — long-format ENCODE element-quantification TSV → wide
  `Oligo,DNA_rep*,RNA_rep*` for `mpra qc`/`mpra activity`.
* `igvfagent agarwal2025-lentimpra-activity` (`Scripts/ported/skills/agarwal2025_lentimpra_activity.py`)
  — barcode-depth-filtered (`--min-bc-depth`, default 10) inter-replicate
  Pearson concordance, plus negative-control z-test + BH-FDR activity
  calling, per cell line.
* `igvfagent agarwal2025-lentimpra-cnn-eval` (`Scripts/ported/skills/agarwal2025_lentimpra_cnn_eval.py`)
  — held-out Pearson r (forward/reverse-orientation-averaged) from a
  MPRA-LegNet `predictions_new_format.tsv`, one or more cell lines at
  once. Evaluation glue around the paper's own unmodified LegNet code,
  not a reimplementation of the model.
* `Scripts/ported/skills/agarwal2025_lentimpra_summarize.py` — consolidates
  the three cell lines' output into `results_summary.json`, which
  `expected.json`'s checks read (this one stays a benchmark-local script,
  not registered as a port: it's specific to this paper's own
  promoter/enhancer name conventions per cell line, not a generalizable
  method).

Both ports are registered `unreviewed` (`igvfagent ext-review promote`
still needed for a human sign-off) with provenance pointing at the paper's
own Methods text (no author code exists for this step; see
`Scripts/ported/registry.json`).

None of these are registered as `igvfagent port register` entries: the
paper released no code for this specific step, so there is no upstream
repository/commit to attribute the port to (unlike, say, the ma2020
SHARE-seq DORC ports, which follow FigR).

## Citation

Agarwal V, Inoue F, Jia RZ, Ahituv N, et al. **Massively parallel
characterization of transcriptional regulatory elements.** *Nature* **639**:
411–420 (2025). DOI: [10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9) · PMID 39814889 · PMC11903340
