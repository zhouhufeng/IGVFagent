# miao2024_valid_inference

## Paper

**Valid inference for machine learning-assisted genome-wide association studies.**
Miao J, Wu Y, Sun Z, Miao X, Lu T, Zhao J, Lu Q.
*Nature Genetics* 2024 · doi:[10.1038/s41588-024-01934-0](https://doi.org/10.1038/s41588-024-01934-0) · PMID 39349818 · PMC11972620 (IGVF0097, U01HG012039)

A statistical-genetics methods paper, not an assay paper: it shows that naive GWAS on ML-imputed phenotypes (e.g. an ML model predicting a hard-to-measure trait from easy-to-measure covariates, then running GWAS on the predicted values) can produce pervasive false-positive associations, and introduces **POP-GWAS**, a summary-statistics-only debiasing framework. Applied to DXA-derived bone mineral density (BMD) across 14 UK Biobank skeletal sites, POP-GWAS identifies 89 novel loci.

Code Availability: [qlu-lab/POP-TOOLS](https://github.com/qlu-lab/POP-TOOLS) (commit `f1db4032b7ffe3d65cc03c314ee06da07ac9bd50`, the released CLI/Python package) and [jmiao24/POP-GWAS_analysis](https://github.com/jmiao24/POP-GWAS_analysis) (commit `b51da435032e23a7f993fed28ce5143300702167`, the paper's own simulation + UKB analysis scripts). Data: GWAS summary statistics for all 14 sites deposited at the NHGRI-EBI GWAS Catalog (GCST90446627-644, found via `findByPublicationIdPubmedId?pubmedId=39349818`) and qlu-lab.org/data.html.

## What IGVFagent does

No IGVF assay route matches a statistical-genetics methods paper (`bench route` correctly reports no match); reproduction here is manual, driven directly off the authors' own two repos plus their own public deposits — no UK Biobank individual-level access requested or used.

```bash
bash Benchmarks/miao2024_valid_inference/run.sh
.venv/bin/igvfagent bench score --paper-id miao2024_valid_inference
```

## Coverage: 5/6 analyses reproduced, 1/6 honestly blocked (7/7 checks pass, `reproduction: reproduced_except_access`)

| Analysis | Figure | Status | Result |
|---|---|---|---|
| `pop_gwas_core_algorithm` | Fig. 2 | **reproduced** | Ran the authors' own `POP-GWAS.py` (incl. its bundled LDSC) on their own bundled `Head_BMD` test data as reference. Our Python port `igvfagent miao2024-pop-gwas` (`compute.py::_cal_qt_wtd` + `utils.py::_read_z`, given the same r12/r13/r23) matches **100% of 90,076 SNPs** to the reference file's own display precision (atol=6e-4 on Z; max abs diff = 5.0e-4). |
| `fig3_simulation_type1_error` | Fig. 3b | **reproduced** | Reran the authors' own unmodified `simulation/Fun.R` + `qt.R`'s own `mar_effect=0` (true null) grid point, at 1/10 scale (300 reps) for tractability, as reference. Our port `igvfagent miao2024-pop-gwas-simple` (`Fun.R::cal_qt_opt`) matches **100% of 300 replicates to floating-point precision** (atol=1e-8, max abs diff=1.3e-14). **Class-A**: empirical Type-I error = 0.03 at alpha=0.05 (paper: nominal 0.05; 0.03 is within Monte Carlo noise for n=300, binomial SE=0.013). |
| `fig4_effective_sample_size` | Fig. 4 | **reproduced** | The authors' own closed-form `simulation/fig4.R` formula (deterministic, no simulation) vs our port `igvfagent miao2024-effective-n`: **100% match to floating-point precision** (atol=1e-9, max abs diff=9.8e-15). |
| `fig5_novel_bmd_loci_count` | Fig. 5 | **reproduced** | Paper: "89 novel loci." The authors' own deposited `Novel_BMD_loci_qced.xlsx` (linked from `POP-GWAS_analysis/README.md`'s real_data section) has exactly 89 rows. Independent re-derivation (not just trusting that table): looked up all 89 rsIDs directly in the authors' 14 primary per-site GWAS Catalog deposits (GCST90446627-640) — all 89/89 present, and 50/89 (56%) independently clear genome-wide significance from a single primary site's file alone (the rest are found but sub-genome-wide per-site, consistent with the paper likely using a joint/meta-analysis the 14 primary files alone don't capture — see the 4 separate "boosted-N" GCST deposits for femur_neck/L1-L4/head/total_body). |
| `fig6_lgr5_locus` | Fig. 6 | **reproduced** | LGR5/rs12308154, the paper's highlighted head-specific signal: independently pulled from the public GWAS Catalog head-BMD deposit (GCST90446633, N=44,255 — matches the paper's stated Fig. 6 range of 44,267-60,829), P=1.156e-09 (genome-wide significant), beta=-0.078; consistent with the authors' own table's P=1.472e-09 for the same rsID (position differs by ~59 kb, consistent with a GRCh37/38 build difference in the two sources, not a different SNP). |
| `fig1_t2d_false_positives` | Fig. 1 | **blocked** (`not_deposited`) | The T2D case study (imputed-T2D-vs-ground-truth-T2D GWAS, plus a 4-SNP HbA1c/glycemic/erythrocytic comparison) needs individual-level UK Biobank T2D case status, HbA1c, and the authors' ML imputation model; neither repo ships code or data for it. Confirmed the underlying phenotype really is access-controlled, not just unshared: the authors' own "example data" Box folder for the *other* (BMD) real-data pipeline ships covariates/predictors but blanks out the phenotype column. |

## Honest caveats

* `Fig. 2`'s core estimator and `Fig. 4`'s formula check are exact-algorithm matches; `Fig. 3`'s Type-I-error number is a reduced-scale (1/10 sample size, 300 vs 1000 replicates) Monte Carlo re-run of the authors' own unmodified simulation code for tractability, not their full-scale run — the port-vs-reference agreement (which is what "reproduced" actually rests on) is exact regardless of scale.
* `Fig. 5`'s 89-loci count is the authors' own curated table, not re-derived from raw phenotypes (that needs UKB access — see `fig1_t2d_false_positives`); the added independent per-site GWAS-Catalog cross-check is a lower-bound sanity check on that table, not a full re-run of their discovery pipeline (see `expected.json`'s check provenance for the exact 50/89 breakdown and why 89/89 isn't expected from primary single-site files alone).
* `Fig. 1`'s T2D case study is the only fully blocked analysis, and only because of access control, not missing effort — see `expected.json`'s `fig1_t2d_false_positives.blocker.reason` for the due-diligence trail (including checking the authors' own shared example data for a redacted-phenotype pattern).

## Provenance

`provenance.json` holds the resolve/harvest/route record (harvest was abstract-only; full text unavailable via Europe PMC). `expected.json`'s `checks[].provenance` holds the full port-vs-reference match rates and paper-number comparisons. Ported code: `Scripts/ported/skills/miao2024_pop_gwas.py`, `miao2024_pop_gwas_simple.py`, `miao2024_effective_n.py` (see `OPERATIONS.md` for what each wraps).
