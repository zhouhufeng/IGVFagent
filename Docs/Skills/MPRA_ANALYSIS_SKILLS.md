# Skill: MPRA Data Retrieval And Analysis

Use this skill when the agent needs to retrieve MPRA/STARR/BlueSTARR evidence
from IGVF Catalog, IGVF Portal, or ENCODE, or to run end-to-end MPRA
analytics (per-element activity, allelic skew, QC) on a local counts table.

The analytical methods are clean-room reimplementations of the canonical
public pipelines (algorithms only, no source copied):

- **Tewhey-lab MPRASuite** (Apache-2.0): `MPRAmodel` per-element NB GLM
  Wald test (DESeq2) of RNA vs DNA with **summit-shift** size-factor
  normalization, then paired-t allelic skew across replicates with BH-FDR.
- **WangLabTHU esMPRA**: count-based activity ratios, replicate concordance,
  barcode-per-oligo QC.

## Retrieval

```bash
python3 Scripts/mpra_data_skills.py pull --source catalog --limit 25
python3 Scripts/mpra_data_skills.py pull --source all --limit 25
python3 Scripts/mpra_data_skills.py portal-manifest --limit 100 --label igvf_portal_mpra_many
```

Use `IGVF_PORTAL_COOKIE` for unreleased Portal datasets.

## Local Result-Table Summaries

```bash
python3 Scripts/mpra_data_skills.py analyze-local --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra
python3 Scripts/mpra_data_skills.py literature-demo --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra_literature_demo
```

These take an existing MPRA result table (with columns like
`MPRA.log2FoldChange`, `MPRA.minusLog10PValue`) and emit summary stats,
SVG plots and a Markdown report.

## End-to-End Counts Analytics (NEW)

Run these on a counts table with columns `Oligo,Allele,SNP,Window,Strand,
Haplotype,DNA_rep1,DNA_rep2,...,RNA_rep1,RNA_rep2,...` (barcode-level
tables with a `Barcode` column are summed per oligo automatically):

```bash
# QC: replicate Pearson r heatmaps, BC/oligo, log10 counts/oligo
python3 Scripts/mpra_data_skills.py qc       --input counts.tsv --label myrun

# Per-oligo activity (NB GLM Wald + summit-shift) -> *_activity.out
python3 Scripts/mpra_data_skills.py activity --input counts.tsv --label myrun

# Allelic skew (paired t-test of ALT vs REF on log2(RNA/DNA)) -> *_skew.out
python3 Scripts/mpra_data_skills.py skew     --input counts.tsv --label myrun

# 4-panel volcano: activity, skew, MA-activity, MA-skew (red = padj < FDR)
python3 Scripts/mpra_data_skills.py volcano      --activity Docs/MPRA/<ts>_myrun_activity.out     --skew     Docs/MPRA/<ts>_myrun_skew.out          --label myrun --fdr 0.05
```

Dependencies (installed once into the project venv): `pandas`, `numpy`,
`scipy`, `statsmodels`, `pydeseq2` (BSD-3).

### Output schemas

- **`*_activity.out`** (TSV) -- one row per Oligo, columns:
  `baseMean, log2FoldChange, lfcSE, stat, pvalue, padj`. The
  `summit_shift_log2` value used to recenter size factors is stored in
  the JSON summary saved next to the table.
- **`*_skew.out`** (TSV) -- one row per paired element, columns:
  `Element, Log2Skew, LogSkew_SE, tstat, pvalue, padj, RefOligo, AltOligo`.
- **`*_qc_report.md`** -- per-condition replicate Pearson stats (median,
  min, max r) and links to BC/oligo + counts/oligo + Pearson heatmap SVGs.

## Literature-Informed Templates

- Cell / PubMed variant-effect MPRA: volcano plots, significant allelic-
  effect tables, integration with GWAS/eQTL/fine-mapping.
- Nature Methods MPRA design benchmarking: count QC, assay-design
  comparisons, activity dynamic range, sequence-context summaries.
- Nature Genetics single-cell MPRA: cell-type activity heatmaps and
  cluster/cell-type specificity summaries.
- Nature Biotechnology regulatory grammar MPRA: saturation-mutagenesis
  maps, position-effect heatmaps, and model-vs-observed plots.
- Nature large-scale cCRE MPRA: activity distributions, cCRE class
  enrichment, genome-browser style views.

## Reuse Rules

- Inspect metadata before downloading large MPRA files.
- For raw counts: always run `qc` first to check replicate concordance
  before trusting activity / skew calls.
- Use `summit-shift` (default in `activity`) to center the log2FC mode
  on zero before declaring "active" elements.
- Preserve variant identifiers, allele orientation, library/source,
  biosample, activity class, counts, log2FC, P/Q values, significance.
- Cross-check active or skewed elements against IGVF Catalog variant-gene,
  variant-biosample, regulatory-element, and enhancer-gene predictions.
