# Advanced Variant Analysis

This skill produces the full IGVF-style integrated variant analysis from a
generic variant CSV: per-variant feature matrix, predicted-functional
composite, logistic models against any user-supplied experimental outcome,
volcano / Miami / overlap plots, and a research-grade markdown report.

## Inputs

- `--input <csv>`: variant list with any of `rsid`, `chrom/pos/ref/alt`,
  `spdi`, `hgvs`. Extra user columns are preserved.
- `--experimental <csv>` (optional): user effect / p-value table; joined on
  `variant_id` or `rsid` automatically (override with `--join-col`).
- `--outcome <col>`: column treated as the binary outcome for logistic
  models (e.g. `BEAN_pval_lt_.05`, `mpra_sig`).
- `--predictors <comma>`: defaults to all `PFunc_*` flags.
- `--gene-list <comma>`: per-gene Miami plots near these symbols.
- `--effect-col <col>`: continuous effect for volcano (paired with `--outcome`).
- `--bottom-p <col>`: optional second p-value track for the Miami plot.

## Outputs (under `Docs/AdvancedVariantAnalysis/`)

- `<stamp>_<label>_annotated.csv`        unified per-variant matrix
- `<stamp>_<label>_summary_stats.csv`    feature counts
- `<stamp>_<label>_logistic_model.json`  marginal + joint logistic fits
- `Plots/<stamp>_<label>_logistic_model.svg`
- `Plots/<stamp>_<label>_evidence_overlap.svg`
- `Plots/<stamp>_<label>_volcano.svg` (when applicable)
- `Plots/<stamp>_<label>_miami_<gene>.svg`
- `<stamp>_<label>_analysis_report.md`   markdown research report

## Predicted-Functional axes

| Flag | Source |
|---|---|
| `PFunc_aPC` | Catalog summary aPC scores (any > 20) |
| `PFunc_MACIE` | Catalog summary MACIE scores (any > 20) |
| `PFunc_ClinVar` | ClinVar pathogenic / likely-pathogenic / drug response |
| `PFunc_QTL` | At least one QTL gene linked in Catalog |
| `PFunc_RegElement` | ENCODE cCRE class match (PLS / pELS / dELS) or Catalog element link |
| `PFunc_Prediction` | At least one prediction-set link in Catalog |

`PFunc_Any` is the OR; `PFunc_Sum` is the count of axes that fire.

## Example

```bash
python3 Scripts/advanced_variant_analysis.py run \
  --input Data/Input/VariantList/example_variants.csv \
  --label example_locus_v1
```

With user experimental data and a custom outcome:

```bash
python3 Scripts/advanced_variant_analysis.py run \
  --input Data/Input/VariantList/my_variants.csv \
  --experimental Data/Input/Experimental/my_crispri.csv \
  --join-col variant_id \
  --outcome BEAN_pval_lt_.05 \
  --gene-list LDLR,PCSK9,APOE \
  --label my_crispri_v1
```

## Reuse rules

- Provide your own variant list and (optionally) experimental table — never
  commit confidential or pre-publication variant data.
- The script caches Catalog responses under `Data/Cache/AdvancedVariantAnnotations/`;
  delete the cache to force re-fetch.
- The ENCODE cCRE BED is downloaded once from `ENCODE_CCRE_BED_URL`. Override
  the env var to use a private mirror.
- Logistic models use IRLS over numpy/scipy only; results match `glm()`
  asymptotically.
- Plots are SVG so they can be edited in any vector editor.
