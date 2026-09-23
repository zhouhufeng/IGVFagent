# MPRA and STARR-seq

Massively parallel reporter assays from library design to allelic effects, STARR-seq and scQers.

**On this page**

- [scQers — quantitative single-cell enhancer reporters](#scqers--quantitative-single-cell-enhancer-reporters)
- [MPRA / STARR / BlueSTARR](#mpra--starr--bluestarr)
- [MPRA library design → counts → barcode QC](#mpra-library-design--counts--barcode-qc)
- [STARR-seq allelic test](#starr-seq-allelic-test)

## scQers — quantitative single-cell enhancer reporters

Two barcodes per construct are the whole trick. A candidate element drives
mCherry carrying an **mBC**; the same construct constitutively expresses an
**oBC**. In single cells the oBC says *which* element a cell received and the
mBC says *how hard it is driving* — which is what makes the readout
quantitative rather than a sort-based enrichment.

```bash
igvfagent scqers extract-bc --in-r1 mBC_R1.fastq.gz --out-file mBC.txt.gz \
    --start 0 --end 15 --check-seq GCT
igvfagent scqers count-bc --in-file mBC.txt.gz --out-file high_mBC.txt \
    --threshold 300 --plot dist.png
igvfagent scqers subassembly --in-file paired_counts.txt --out-file oBC_mBC.tsv
igvfagent scqers pipeline --counts joined_counts.tsv       # activity → specificity → calls
igvfagent scqers selftest                                  # synthetic data, no inputs
```

Activity is bootstrapped against the library's `minP` and `noP` controls **at
matched sample size**, so a rare element is not beaten by a control that was
merely pooled over more cells. Specificity is a permutation of the
cell-to-cluster assignment, because the best of several clusters has a high
fold-change by construction and a fold-change alone cannot tell you otherwise.
p-values are empirical and BH-corrected within each replicate, and an element
is called only if it clears the FDR in *every* replicate it was measured in.

Method: Shendure lab, *Multiplex profiling of developmental enhancers with
quantitative, single-cell expression reporters*
([shendurelab/scQers](https://github.com/shendurelab/scQers), MIT). The
upstream repository is R and shell shared "for transparency" rather than as a
pipeline; this is the method as a runnable CLI. Barcode extraction and
counting reproduce the upstream example output exactly — all 191 barcodes
above threshold, with identical counts.

## MPRA / STARR / BlueSTARR

```bash
python3 Scripts/mpra_data_skills.py pull --source catalog --limit 25
python3 Scripts/mpra_data_skills.py portal-manifest --limit 100 --label igvf_portal_mpra_many
python3 Scripts/mpra_data_skills.py analyze-local \
  --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra
python3 Scripts/mpra_data_skills.py literature-demo \
  --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra_literature_demo
python3 Scripts/mpra_data_skills.py write-playbook
```

## MPRA library design → counts → barcode QC

Three clean-room ports of the kircherlab MPRA suite (Max Schubach / BIH, MIT),
covering the whole assay lifecycle. All standard-library-only — no Snakemake,
conda, R, bedtools, pandas or AnnData — and each verified against the upstream
implementation.

```bash
# ── 1. DESIGN the library ────────────────────────────────────────────────
# Tile regions into oligo windows (3 strategies chosen by region length)
igvfagent oligo tile --regions regions.bed.gz --oligo-length 200 --min-overlap 50

# One REF + one ALT oligo per variant x region — the contrast the assay reads
igvfagent oligo design-variants --regions regions.bed.gz \
    --variants variants.vcf.gz --reference hg38.fa --variant-edge-exclusion 20

# Filters: homopolymers, EcoRI/SbfI sites, simple repeats, TSS, CTCF motifs
igvfagent oligo filter --design design.fa --regions regions.bed.gz

# Everything, one call → synthesis-ready design.fa.gz + per-drop filter log
igvfagent oligo pipeline --regions regions.bed.gz --variants variants.vcf.gz \
    --reference hg38.fa --tile --left-adapter AGGACCGGATCAACT \
    --right-adapter CATTGCGTGAACCGA --label my_library

# ── 2. SEQUENCING → per-oligo ACTIVITY ───────────────────────────────────
# Trust the barcode→oligo assignment (min reads + majority agreement)
igvfagent mpraflow assign-filter --pairs pairs.tsv.gz --minimum 3 --fraction 0.75

# Quantify one replicate, with optional barcode outlier removal
igvfagent mpraflow merge-counts --counts rep1.tsv.gz \
    --assignment assignment.tsv.gz --outlier-detection ratio_mad

# "Would deeper sequencing help?" — Lincoln-Petersen mark-recapture
igvfagent mpraflow complexity --barcode-file reporter_experiment.barcode.tsv.gz

# Everything: per-replicate activity → master table → pooled → allelic skew
igvfagent mpraflow pipeline --assignment assignment.tsv.gz \
    --replicate 1=rep1.tsv.gz --replicate 2=rep2.tsv.gz --replicate 3=rep3.tsv.gz \
    --declaration variants.tsv --threshold 10

# ── 3. BARCODE QC + IGVF FORMAT VALIDATION ───────────────────────────────
# Three outlier detectors from Rosen et al.
igvfagent mpralib outliers --barcode-file barcodes.tsv.gz --method global
igvfagent mpralib outliers --barcode-file barcodes.tsv.gz --method large_expression

# Are those calls reproducible across replicates?
igvfagent mpralib consistency --barcode-file barcodes.tsv.gz --method global

# Validate any file against the 8 IGVF MPRA community standards
igvfagent mpralib schemas
igvfagent mpralib validate --file master_table.tsv.gz --schema reporter_experiment
```

**Verification.** `oligo` reproduces all 35 rows of MPRAOligoDesign's golden
tiling fixture. `mpraflow` is **bit-identical** to six upstream scripts run
side-by-side on the same inputs (`merge_label.py` 420/420 values,
`filterAssignmentTsv.py` 238/238 rows, `generateMasterVariantTable.py` 435/435).
Against the published IGVF artefact for AnalysisSet `IGVFDS1933XFAF` —
where the portal hosts both the input *and* the authors' own pipeline output —
it reproduces **210,660 / 210,660 values exactly**, and reproduces the paper's
stated library-complexity figures to the digit. See
[`Benchmarks/rosen2025_mprasnakeflow/`](../../../Benchmarks/rosen2025_mprasnakeflow/README.md).

Two upstream bugs are fixed and documented (a homopolymer run at a sequence's
end was never checked; `use_most_centered_region` measured distance to the
region *end*), each with a flag to restore the original behaviour for exact
reproduction.

## STARR-seq allelic test

Clean-room rewrite of
[gaochengwen/STARR-seq-Data-Analysis](https://github.com/gaochengwen/STARR-seq-Data-Analysis)
(no LICENSE — every line is paraphrased from the published `mpra::mpralm`
methods). Implements TPM counts QC + RLE + Spearman D-stat, per-(SNP,
Allele) aggregation, log activity, and the moderated allelic test with
Smyth-2004 trigamma-inversion eBayes + BH-FDR. See
[`Docs/Skills/STARRSEQ_ANALYSIS_SKILLS.md`](../../Skills/STARRSEQ_ANALYSIS_SKILLS.md).

```bash
# Discover IGVF Portal STARR-seq MeasurementSets
igvfagent starrseq pull-portal --limit 50 --label survey

# TPM QC: filter low-expression fragments, RLE matrix, outlier samples
igvfagent starrseq qc --input counts.tsv --label run1

# Collapse barcode-level counts to per-(SNP, Allele)
igvfagent starrseq aggregate --input barcode_counts.tsv --label run1

# Per-fragment log activity + per-SNP allelic skew (descriptive)
igvfagent starrseq activity --input aggregated.tsv --label run1

# mpralm-style allelic test (OLS + eBayes + BH-FDR)
igvfagent starrseq allelic-test --input aggregated.tsv --label run1
```

---

[← Documentation index](../README.md) · [Project README](../../../README.md)
