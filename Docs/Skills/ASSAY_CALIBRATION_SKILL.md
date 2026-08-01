exCALIBR — functional-assay calibration to ACMG/AMP evidence (IGVFagent port).

A clean-room Python reimplementation of **exCALIBR**
(github.com/rosstewart/exCALIBR, MIT, R. Stewart, Northeastern), which
implements the gene-based assay-calibration method of

    Zeiberg D, Tejura M, McEwen AE, Fayer S, Pejaver V, Rubin AF,
    Starita LM, Fowler DM, O'Donnell-Luria A, Radivojac P.
    "Gene-based calibration of high-throughput functional assays for
    clinical variant classification." bioRxiv 2025.04.29.651326.

The method turns a raw multiplexed-assay score (MAVE / VAMP-seq / SGE /
cell-fitness …) into **calibrated ACMG/AMP evidence strengths** (PS3 / BS3
at supporting / moderate / strong / very-strong), so an assay score can be
used directly in clinical variant classification instead of as a bare number.

Chain, faithful to upstream::

  1. Label variants into four samples from the scoreset table —
     0 P/LP (ClinVar), 1 B/LB (ClinVar), 2 population (gnomAD), 3 synonymous.
  2. Fit a **constrained skew-normal mixture** by EM: K components whose
     (skew, loc, scale) are *shared* across samples, each sample carrying its
     own mixing weights. The constraint forces the density ratio of adjacent
     components to be monotone (so the resulting LR+ is monotone in score);
     it is enforced parameter-by-parameter by binary search inside every
     M-step.
  3. **Bootstrap** the whole fit (per-sample resampling; best-of-N fits per
     bootstrap chosen on held-out log-likelihood).
  4. Estimate the **prior** P(pathogenic | population) per fit by EM against
     the population sample; take the median across bootstraps.
  5. Build **LR+(score) = f_P(score) / f_B(score)** per bootstrap, take the
     conservative percentile envelope across bootstraps.
  6. Solve for **Tavtigian's C (O_PVSt)** at that prior (Tavtigian 2018
     Bayesian ACMG framework), giving evidence thresholds C^(points/8), and
     convert the LR+ curve into **score ranges per evidence point value**.
  7. Optionally choose between a 2- and 3-component model by a paired
     Wilcoxon / 5th-percentile test on bootstrap validation likelihoods.

Reimplemented from the algorithm and the papers; no upstream source lines
copied. Deviations from upstream, all deliberate:

  * **Deterministic**: every random draw (bootstrap split, k-means seed,
    method-of-moments cut points, skew-sign table index) is derived from an
    explicit seed, so a rerun reproduces the calibration bit-for-bit. This
    matches IGVFagent's cross-backend consistency contract.
  * **stdlib parallelism** (`concurrent.futures`) instead of joblib; SLURM
    generation dropped in favour of a resumable in-process run — bootstraps
    are appended to a JSON-lines ledger with a progress heartbeat, so one
    long call finishes the job and `--resume` picks up where it stopped.
  * Empty sample categories are dropped at load time (with a name→column
    map) rather than index-shifted downstream.

CLI::

    igvfagent calibrate thresholds --prior 0.1
    igvfagent calibrate prepare --pillar MSH2_Jia_2021.csv --name MSH2_Jia_2021
    igvfagent calibrate run --table scores.csv --name MSH2_Jia_2021 \
        --components 2 3 --n-bootstraps 1000 --fits-per-bootstrap 100
    igvfagent calibrate assign --calibration MSH2_2c_calibration.json \
        --scores my_variants.csv
    igvfagent calibrate selftest

License: Apache-2.0 (upstream exCALIBR is MIT; method credited above).
Heavy deps imported lazily: numpy, scipy, matplotlib, scikit-learn
(k-means only — a numpy fallback is used when absent).
