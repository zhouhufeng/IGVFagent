# Skill: STARR-seq Allelic Analysis

Use this skill when the agent needs to run end-to-end **allelic STARR-seq**
counts analytics on data from IGVF Portal STARR-seq MeasurementSets or any
other source with paired DNA + RNA replicate counts per (SNP, allele).

The math is a clean-room reimplementation of the methods in
[gaochengwen/STARR-seq-Data-Analysis](https://github.com/gaochengwen/STARR-seq-Data-Analysis)
— the awk + R pipeline that ultimately calls Bioconductor `mpra::mpralm`
for allelic significance. The reference repository has no LICENSE file,
so this skill ships a pure-Python rewrite that follows the underlying
published methods (limma + voom + eBayes, BH-FDR) only.

## Commands

```bash
# 1. Discover IGVF Portal STARR-seq MeasurementSets
python3 Scripts/starr_seq_skill.py pull-portal --limit 50 --label survey

# 2. Counts QC (TPM scaling + RLE + Spearman D-stat outliers)
python3 Scripts/starr_seq_skill.py qc --input counts.tsv --label run1

# 3. Collapse barcode-level counts to per-(SNP, Allele) wide form
python3 Scripts/starr_seq_skill.py aggregate --input barcode_counts.tsv --label run1

# 4. Per-fragment log activity + per-SNP allelic skew (descriptive)
python3 Scripts/starr_seq_skill.py activity --input aggregated.tsv --label run1

# 5. mpralm-style allelic test (voom-weighted GLS + eBayes moderation + BH-FDR)
python3 Scripts/starr_seq_skill.py allelic-test --input aggregated.tsv --label run1
```

## Input schema

`counts.tsv` columns (TSV/CSV, auto-detected):

| Column | Required | Notes |
|---|---|---|
| `SNP` | yes | SNP identifier |
| `Allele` | yes | `REF`/`ALT` or `A`/`B` |
| `Fragment` | optional | Defaults to `SNP_Allele` |
| `DNA_rep1`, `DNA_rep2`, ... | yes | Library (DNA) replicate counts |
| `RNA_rep1`, `RNA_rep2`, ... | yes | Output (RNA) replicate counts |

Replicate columns are auto-detected by the `DNA`/`RNA` prefix (case-insensitive).
For paired analysis, replicate index `1, 2, ...` in `DNA_rep*` corresponds to
the same index in `RNA_rep*`.

## Output schema — `*_allelic_test.tsv`

| Column | Meaning |
|---|---|
| `SNP` | SNP id |
| `beta_allele2` | Effect of ALT allele on log2(RNA/DNA), GLS estimate |
| `se_moderated` | eBayes-moderated standard error |
| `t_moderated` | Moderated t-statistic |
| `pvalue` | Two-sided p from moderated t |
| `padj` | BH-FDR adjusted p |
| `d0_prior_df` | Empirical-Bayes prior degrees of freedom |
| `s0_sq_prior` | Empirical-Bayes prior variance |

## Notes

- The reference repo also includes an awk + Biopython sequence-alignment
  step (`extract_alleles.py`) that QC's whether each oligo aligns to the
  intended reference window. That step is not yet absorbed here — IGVF
  Portal STARR-seq deposits already report aligned counts.
- Genome-wide STARR-seq peak callers (CRADLE, STARRPeaker, BasicSTARRseq)
  are out of scope; this skill covers the allelic-test surface only.
- IGVF Portal: 5,239 STARR-seq files and 1,066 BlueSTARR files are
  catalogued as of the last full Portal survey. Use `starrseq pull-portal`
  to refresh the list.
