# Skill: MULTI-seq / Cell Hashing demultiplexing

Python port of **deMULTIplex2** (Zhu et al., *Nat Methods* 2024,
https://github.com/Gartner-Lab/deMULTIplex2), the v2 classifier for
**MULTI-seq** (McGinnis, Patterson, Winkler et al., *Nat Methods* 2019).
Fits a two-component negative-binomial mixture per sample tag via EM
and classifies every cell as **singlet** (named by tag), **multiplet**,
or **negative**.

## Background — what MULTI-seq actually does

MULTI-seq tags every cell (or nucleus) with a sample-specific **8-nt
barcode** held in the membrane by two **lipid-modified oligonucleotides
(LMOs)** that hybridize to the barcode:

| Component | Lipid | Role |
|---|---|---|
| Anchor LMO | 5'-lignoceric acid (C24) | primary hydrophobic insertion |
| Co-anchor LMO | 3'-palmitic acid (C16) or cholesterol | stabilises retention |
| Sample barcode oligo | none | unique 8-nt tag, Hamming ≥ 3 |

Labelling is a **two-step, 10-minute, 4 °C** protocol (200 nM final;
~2.5 µM for fixed / cryopreserved tissue). No permeabilisation, no
antibodies — works on live cells, isolated nuclei, and PDX/tumour
samples. Read structure on **10x Chromium 3' v2/v3**:

    R1 :  cell barcode (1-16) + UMI (17-26 v2 / 17-28 v3)
    R2 :  sample tag    (1-8)  + 30 nt poly-A capture

McGinnis et al. demonstrated **96-plex** in HMEC, snRNA-seq, MEFs,
primary human cells, and cryopreserved PDX tumours. The headline
operational win is **'Super-Loading'**: by recovering ~85 % of
inter-species doublets, you can over-load a 10x lane ~5× before the
doublet rate is uninterpretable.

## Why deMULTIplex2 (the v2 algorithm we port here)

The original 2019 paper shipped a thresholding scheme — Gaussian KDE
on each tag's log-counts, sweep candidate thresholds (q = 0.2–0.99)
to maximise singlet count, then recursively remove negatives. That
collapses when the background mode is poorly separated (few samples,
one dominant population, fixed tissue) and produces an unstable
classification across runs. **deMULTIplex2 (Zhu et al. 2024)**
replaces it with a generative two-component **negative-binomial
mixture** fit by EM per tag. This module ports that v2 algorithm
faithfully (NB2 via `statsmodels` matches `MASS::glm.nb`).

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

## References

- **McGinnis et al. 2019** — original MULTI-seq method:
  *MULTI-seq: sample multiplexing for single-cell RNA sequencing using
  lipid-tagged indices.* **Nature Methods** 16(7):619–626.
  DOI [10.1038/s41592-019-0433-8](https://www.nature.com/articles/s41592-019-0433-8)
  · PMID 31209384 · PMCID PMC6837808.
- **Zhu et al. 2024** — deMULTIplex2 (the v2 EM classifier ported here):
  *Robust sample demultiplexing for scRNA-seq.*
  GitHub: [Gartner-Lab/deMULTIplex2](https://github.com/Gartner-Lab/deMULTIplex2).
- **Sigma-Aldrich technical article** — LMO001 reagent kit + protocol:
  [MULTI-seq: Sample Multiplexing for Single-cell Analysis Sequencing](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/sequencing/multi-seq-sample-multiplexing-single-cell-analysis-sequencing).

BibTeX:

```bibtex
@article{mcginnis2019multiseq,
  title   = {{MULTI-seq}: sample multiplexing for single-cell {RNA}
             sequencing using lipid-tagged indices},
  author  = {McGinnis, Christopher S. and Patterson, David M. and
             Winkler, Juliane and Conrad, Daniel N. and Hein, Marco Y. and
             Srivastava, Vasudha and Hu, Jennifer L. and Murrow, Lyndsay M. and
             Weissman, Jonathan S. and Werb, Zena and Chow, Eric D. and
             Gartner, Zev J.},
  journal = {Nature Methods},
  volume  = {16}, number = {7}, pages = {619--626}, year = {2019},
  doi     = {10.1038/s41592-019-0433-8},
  pmid    = {31209384}, pmcid = {PMC6837808}
}

@misc{zhu2024demultiplex2,
  title        = {{deMULTIplex2}: robust sample demultiplexing for scRNA-seq},
  author       = {Zhu, Q. and {Gartner Lab}},
  year         = {2024},
  howpublished = {GitHub: Gartner-Lab/deMULTIplex2},
  url          = {https://github.com/Gartner-Lab/deMULTIplex2}
}
```
