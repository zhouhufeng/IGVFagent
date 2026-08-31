# MPRA barcode QC (MPRAlib) — skill playbook

Python port of the analysis core of
[MPRAlib](https://github.com/kircherlab/MPRAlib) (Max Schubach, Berlin
Institute of Health at Charité; MIT). Cite Rosen *et al.*,
*Genome Research* (2025), doi:10.1101/gr.281462.125.

    oligo  →  ORDER  →  SEQUENCE  →  mpraflow  →  this skill
    design                          counts       barcode QC

## Why barcode outliers matter

An oligo's activity is an average over its barcodes. A few barcodes with
wild RNA counts — a PCR jackpot, an unlucky integration site — can drag
the oligo's ratio away from where its other twenty barcodes sit. The
paper defines three detectors, each catching a different failure:

| Method | Question | Rule |
|---|---|---|
| `global` | extreme for the whole replicate? | \|z\| > 3 on RNA counts across all barcodes |
| `oligo` | extreme for its own oligo? | \|z\| > 3 on RNA counts within the oligo |
| `large_expression` | runaway activity? | log2 ratio > 5 above the oligo's median |

`large_expression` is one-sided by design — it looks for runaway
expression, not silence — and yields one verdict per barcode rather than
one per replicate.

## Should you filter on these?

Usually **no**. Rosen *et al.* compared variant effects before and after
removal and found correlations above 0.95, concluding that "explicit
outlier removal is not necessary for variant effect analysis when using
BCalm". Treat a high rate as information about the experiment, not as a
cleaning step to apply by reflex.

What the rates *do* tell you (paper, Figure 4A):

| | lentiviral | episomal (plasmid) |
|---|---|---|
| global | 0.67–1.56% | 0.12–0.48% |
| oligo-specific | 0.85–1.86% | 0.85–1.86% |
| large expression | up to 0.01% | 0.05–0.53% |

## Replicate consistency

`consistency` reports what fraction of outlier barcodes are flagged in
*every* replicate. A detector firing on the same barcodes each time is
describing the library; one firing on different barcodes is describing
noise. The paper found episomal assays far more consistent (63.2–83.9%)
than lentiviral ones (5.6–49.0%), and differentiated cells more
consistent than progenitors.

## Usage

```bash
# Outlier rates, one method at a time
igvfagent mpralib outliers --barcode-file reporter_experiment.barcode.tsv.gz \
    --method global
igvfagent mpralib outliers --barcode-file reporter_experiment.barcode.tsv.gz \
    --method large_expression

# How reproducible are those calls across replicates?
igvfagent mpralib consistency \
    --barcode-file reporter_experiment.barcode.tsv.gz --method global

# Per-barcode normalised activity
igvfagent mpralib activity --barcode-file reporter_experiment.barcode.tsv.gz
```

Input is the IGVF **reporter experiment barcode** file, the same one
`mpraflow` consumes — downloadable straight from the portal.

## Deviations from upstream

* No AnnData / numpy / pandas. Upstream materialises a
  barcodes × replicates matrix; the 240K libraries have 20 M barcodes,
  where that costs several GB. This streams the file instead.
* pandas' `std()` is ddof=1 and skips NaN — both reproduced explicitly,
  since the population SD would shift every z-score.
* Upstream's `transform("std").fillna(0).replace(0, 1)` for
  oligo-specific z-scores (single-barcode or zero-variance oligos get
  std = 1) is reproduced exactly.
