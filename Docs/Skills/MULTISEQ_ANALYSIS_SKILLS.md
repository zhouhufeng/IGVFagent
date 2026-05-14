# Skill: MULTI-seq / Cell Hashing demultiplexing

Python port of **deMULTIplex2** (Zhu et al., *Nat Methods* 2024,
https://github.com/Gartner-Lab/deMULTIplex2). Fits a two-component
negative-binomial mixture per sample tag via EM and classifies every
cell as **singlet** (named by tag), **multiplet**, or **negative**.

## Algorithm in one paragraph

For each tag *j* the routine fits two NB GLMs with a log link:

  * `fit0`:  `bc_umi ~ log(tt_umi)`            (off-target background model)
  * `fit1`:  `(tt_umi - bc_umi) ~ log(tt_umi)` (background tags in true positives)

`bc_umi` is the focal-tag count and `tt_umi` is the per-cell total tag
UMI. EM is initialized by a cosine-similarity cut, then iterates
through E-/M-steps until the log-likelihood is stable (default ≤ 10
iter, tolerance 1e-3). A cell is positive for *j* when
`P(positive | data) > 0.5`. After looping over all tags:

  * 1 positive tag  →  **singlet** (named by that tag)
  * 0 positive tags →  **negative**
  * ≥2 positive tags →  **multiplet**

Randomized quantile residuals (Dunn & Smyth 1996) and Pearson
residuals are computed per tag and returned for downstream UMAP /
QC.

## Subcommands

### demultiplex
```
igvfagent multiseq demultiplex --input tag_counts.csv \
    --label demo --prob-cut 0.5 --residual-type rqr
```
Headline classifier. Writes
`classifications.csv` / `posteriors.csv` / `residuals.csv` /
`tag_coefficients.csv` + diagnostic plots under
`Docs/MultiSeq/<ts>_<label>/`.

### histogram / heatmap (standalone diagnostics)
```
igvfagent multiseq histogram --input tag_counts.csv
igvfagent multiseq heatmap   --input tag_counts.csv \
    --calls classifications.csv
```
`heatmap` will run a quick demultiplex itself when `--calls` is omitted.

### simulate (smoke test)
```
igvfagent multiseq simulate --n-cells 2000 --n-tags 6 \
    --doublet-rate 0.08 --negative-rate 0.05 --label smoke
```
Generates a synthetic cell × tag matrix with ground truth — useful for
self-tests and demos.

### pipeline (one-shot)
```
igvfagent multiseq pipeline --input tag_counts.csv \
    --ground-truth ground_truth.csv --label end_to_end
```
Load → demux → all plots → markdown report + (if `--ground-truth`
given) accuracy table vs known assignments.

### write-playbook
```
igvfagent multiseq write-playbook
```

## Input formats

| Extension | Where the tag counts live |
|---|---|
| `.h5ad`            | `adata.obsm['multiplexing']` / `obsm['mux_counts']` / `obsm['hashing']` / `obsm['HTO']` / `obsm['X_hto']` (fallback: `adata.X`) |
| 10x `.h5`          | features tagged `Multiplex Capture` / `Antibody Capture` |
| `.csv` / `.tsv`    | cells × tags (auto-oriented if rows ≪ cols) |

## Cross-skill chaining

- After `sc-analyze` clusters a 10x dataset, run `multiseq demultiplex`
  on the multiplexing matrix to assign sample identity, then write the
  joined obs back into the AnnData for downstream cluster × sample
  analysis.
- `portal_multiome` and `portal_scrna_10` discover IGVF multiome
  datasets that often include `MULTI-seq` barcoding — the
  classifications produced here drop straight into their per-pool
  manifests.
- `rnaseq deg` over multiplexed populations: use `multiseq` to subset
  each sample cleanly first.

## Dependencies

`numpy`, `scipy`, `pandas`, `statsmodels`, `matplotlib` (all in the
`analysis` extras). Optional: `anndata` / `scanpy` for `.h5ad` and 10x
`.h5` inputs.

## Notes vs the R original

* NB GLMs fit via `statsmodels.discrete.discrete_model.NegativeBinomial`
  (NB2 parameterization). Joint MLE of intercept + slope + alpha —
  matches `MASS::glm.nb`.
* EM loop, posterior thresholding, RQR residuals, and final
  singlet / multiplet / negative call follow the R source line-for-
  line.
* FASTQ → tag-count alignment (the R `readTags` / `alignTags` step) is
  *not* ported — that pipeline is better served by Cell Ranger's
  Feature Barcoding workflow or the original `deMULTIplex` aligner.
  Hand this skill a tag count matrix produced upstream.
