# MPRA count & assignment — skill playbook

Python port of the analysis core of
[MPRAsnakeflow](https://github.com/kircherlab/MPRAsnakeflow)
(Max Schubach, Berlin Institute of Health at Charité; MIT).
Cite Rosen *et al.*, *Genome Research* (2025), doi:10.1101/gr.281462.125.

`igvfagent oligo` designs the library; this skill processes what comes
back off the sequencer.

    design → ORDER → transfect → SEQUENCE → this skill → activity

## Why the two-step structure

An MPRA measures activity indirectly. Every oligo carries random
barcodes; you sequence DNA (input) and RNA (output) and the activity is
the RNA/DNA ratio. Two things have to happen first:

1. **Assignment** — a separate sequencing run tells you which barcode
   belongs to which oligo. A barcode is trusted only when at least
   `--minimum` reads agree and at least `--fraction` of them name the
   same oligo. Barcodes whose reads disagree are *ambiguous* and dropped,
   never guessed.
2. **Normalisation** — counts are summed over an oligo's barcodes,
   divided by the number of barcodes, and scaled per million, so an
   oligo with 40 barcodes is comparable to one with 11.

## Pipeline

| Step | Subcommand | In → out |
|---|---|---|
| Count barcodes | `count-bc` | FASTQ → barcode counts |
| Assign (no aligner) | `assign-exact` | reads + design → barcode↔oligo pairs |
| Trust the assignment | `assign-filter` | pairs → assignment |
| Quantify | `merge-counts` | barcode counts + assignment → per-oligo activity |
| Library complexity | `complexity` | barcode file → Lincoln-Petersen estimate |
| Stack replicates | `master-table` | per-replicate tables → master table |
| Aggregate | `combine-replicates` | master → one row per oligo |
| Per-barcode matrix | `barcode-matrix` | per-replicate barcodes → matrix |
| Allelic skew | `variant-table` | REF/ALT pairs → log2 skew |
| Everything | `pipeline` | assignment + replicates → all of the above |

## Usage

```bash
# Assignment, starting from aligner output (barcode-sorted triples)
igvfagent mpraflow assign-filter --pairs pairs.tsv.gz \
    --minimum 3 --fraction 0.75

# Or skip the aligner entirely when inserts match the design exactly
igvfagent mpraflow assign-exact --design design.fa \
    --reads inserts.fastq.gz --barcodes index.fastq.gz --barcode-length 15

# Quantify one replicate
igvfagent mpraflow merge-counts --counts rep1_counts.tsv.gz \
    --assignment assignment.tsv.gz --min-rna-counts 1 \
    --outlier-detection ratio_mad

# Everything, three replicates, with allelic skew
igvfagent mpraflow pipeline --assignment assignment.tsv.gz \
    --replicate rep1=rep1.tsv.gz --replicate rep2=rep2.tsv.gz \
    --replicate rep3=rep3.tsv.gz --declaration variants.tsv \
    --threshold 10 --label my_experiment
```

## Outlier removal

`--outlier-detection ratio_mad` bins barcodes by RNA depth and trims
those whose DNA/RNA ratio deviates from their oligo's median by more
than `--outlier-mad-times` MADs. `rna_counts_zscore` instead drops
barcodes more than `--outlier-zscore-times` SDs from their oligo's mean
RNA count.

Two upstream quirks are **reproduced by default**, because they define
the numbers in published MPRAsnakeflow output:

| Quirk | Effect | Opt-in fix |
|---|---|---|
| MAD test is one-sided (`ratio_diff <= times*mad`) | deliberate — Rosen et al. define it as ratios "exceeding 5 log2-units *above* their oligo-specific median" | `--mad-two-sided` |
| Bin edges are `arange(0, n_bins)/n_bins`, stopping at the 95th percentile | barcodes above it get no bin, a NaN MAD, and are always dropped — the top 5% by RNA count | `--mad-include-top-bin` |

## Deviations from upstream

* No Snakemake, conda, R, pandas, numpy or pysam — the R master-table
  step and the pandas group-by/quantile logic are reimplemented on the
  standard library.
* **Alignment is not reimplemented.** Upstream offers five aligner
  backends (bwa, bbmap, pbmm2, exact, hybrid) to produce barcode↔oligo
  pairs. `assign-filter` starts from the barcode-sorted
  `barcode <tab> oligo <tab> quality` table they emit; `assign-exact`
  covers the alignment-free case directly. For bwa/bbmap runs, align
  with your existing tooling and feed the result to `assign-filter`.
* Upstream's ambiguous-barcode branch reports the quality of whichever
  oligo the loop happened to end on (an uninitialised-variable leak);
  this reports the modal oligo's quality.

## Reproducing the paper

`complexity` implements the Lincoln-Petersen estimate from Rosen *et al.*
(2025): each replicate is a capture of the barcode pool, barcodes seen in
two replicates are recaptures, and the gap between observed and estimated
barcodes is what the sequencing missed. Verified against the published
figures:

| Dataset | Metric | Paper | This skill |
|---|---|---|---|
| 8K-neurons | median assigned barcodes | 1,444,480 | 1,444,480 |
| 8K-neurons | median Lincoln index | 1,523,572 | 1,523,572 |
| 80K-neurons | median assigned barcodes | 5,459,247 | 5,459,247 |
| 80K-neurons | median Lincoln index | 6,243,618 | 6,243,618 |

The count-aggregation stage reproduces the published IGVF
`reporter experiment` artefact byte-for-byte (210,660/210,660 values). See
`Benchmarks/rosen2025_mprasnakeflow/`.

Note the portal's `reporter experiment` files were produced with
`--min-dna-counts 1`, not the tool default of 0. On the default a DNA
pseudocount of 1 per barcode is added, inflating `dna_counts` by `n_bc` per
oligo and shifting every `log2FoldChange`.

## Key outputs

Everything lands in `Docs/MPRASnakeflow/<timestamp>_<label>/`:

| File | Contents |
|---|---|
| `assignment.tsv.gz` | barcode → oligo, quality, read support |
| `counts.<rep>.tsv.gz` | per-oligo DNA/RNA counts, normalised, log2FC, n_bc |
| `master_table.tsv.gz` | all replicates stacked, filtered on n_bc |
| `combined.tsv.gz` | one row per oligo, pooled and per-replicate means |
| `barcode_matrix.tsv.gz` | per-barcode DNA/RNA across replicates |
| `variants.tsv.gz` | REF/ALT pairs with `log2FoldChange_expression` |
| `summary.json` | per-replicate and overall counts |
