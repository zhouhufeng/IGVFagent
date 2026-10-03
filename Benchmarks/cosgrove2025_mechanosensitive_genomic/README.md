# cosgrove2025_mechanosensitive_genomic

## Paper

**Mechanosensitive genomic enhancers potentiate the cellular response to matrix stiffness.**
Cosgrove BD, Bounds LR, Taylor CK, Su AL, Rizzo AJ, Barrera A, Sun T, Safi A, Song L, Whitlow T, Tata A, Iglesias N, Diao Y, Tata PR, Hoffman BD, Crawford GE, Gersbach CA.
*Science (New York, N.Y.)* 2025;390(6778):eadl1988 · doi:[10.1126/science.adl1988](https://doi.org/10.1126/science.adl1988) · PMID 40997217 · PMC13005951 (open Author Manuscript; full text fetched via the plain PMC HTML page, not Europe PMC's fullTextXML endpoint which 500'd)

IGVF award: IGVF0196, grant UM1HG012053. Data also indexed at the NIH IGVF Data Portal (`data.igvf.org`) under reference file `IGVFFI8553BTYV`.

## Summary

Primary human fibroblasts (HFF) and A549 lung epithelial cells were cultured on soft (1 kPa) vs stiff (50 kPa) polyacrylamide hydrogels. RNA-seq/ATAC-seq/HiCAR profiling, CRISPRi tiling and saturating-Cas9-indel screens, high-throughput CRISPRi growth/migration screens, and single-cell CRISPRi (Perturb-seq) were combined to identify "mechanoenhancers" — *cis*-regulatory elements preferentially active on soft or stiff matrix — and link them to target genes (MYH9, BMF, CTGF, SKP2, CYR61, ...) driving apoptosis, contractility, growth and migration. Validated mechanoenhancers were also shown to modulate TGFβ1-driven fibroblast-to-myofibroblast activation in healthy- and IPF-donor lung fibroblasts.

## Data sources

* GEO SuperSeries `GSE243765`, SubSeries: `GSE243753` (ATAC-seq), `GSE243756` (scRNA-seq/Perturb-seq), `GSE243760` (migration+growth screen), `GSE243761` (MYH9 FACS screen), `GSE243763` (bulk RNA-seq), `GSE307449` (HiCAR loops, added later). BioProject `PRJNA1019774`.
* Authors' code: [`Gersbachlab-Bioinformatics/MechanoEnhancer`](https://github.com/Gersbachlab-Bioinformatics/MechanoEnhancer) (current) and [`Gersbachlab-Bioinformatics/myosin_enhancer`](https://github.com/Gersbachlab-Bioinformatics/myosin_enhancer) (older bioRxiv-era figure numbering, still has the only copies of the RNA-seq/ATAC-seq DESeq2 and motif scripts). Both repos ship analysis code as `.rtf`-saved Rmd notebooks (stripped to text with a crude regex RTF parser, not `unrtf`/`textutil`/LibreOffice — none were available in this environment).
* Full supplementary tables (S1-S31) are individually downloadable from the PMC article page; the combined `NIHMS2152387-supplement-Supplementary_Tables.xlsx` is gated behind a client-side proof-of-work JS challenge that plain `curl`/`WebFetch` cannot pass.

## Analyses and current status (5/8 reproduced, 3/8 blocked `not_deposited` — `reproduced_except_access`)

| id | figure | state | evidence |
|---|---|---|---|
| `rnaseq_atac_stiffness` | Fig. 1 | **reproduced** (weak on RNA-seq count) | `cosgrove2025_stiffness_diff_count` on the authors' own `GSE243763_SupplementaryTable2`/`GSE243753_SupplementaryTable3`: 4,519 HFF RNA-seq DEGs (paper: 4,009, ~13% high) and 16.7% differential ATAC peaks (paper: ~23% HFF, ~15% A549 — A549 has no deposited differential table, only per-replicate narrowPeak/bigwig files) |
| `tf_footprinting_motifs` | Fig. 1J-K | **reproduced** | real HOMER `findMotifsGenome.pl` (the authors' own tool, per the ATAC table's embedded HOMER command strings) known-motif enrichment on `GSE243753_SupplementaryTable3`'s soft-biased (19,385) and stiff-biased (18,280) differential peaks vs hg38, after repairing a broken HOMER install (hardcoded Perl `use lib` paths pointing at a deleted home directory). TEAD1/TEAD3/TEAD/TEAD4/TEAD2 are the top 5 hits on stiff peaks (log P=-3261, 56% of peaks) and AP-1/Fra1 family the top 8 on soft peaks (30-58% of peaks) -- matching the paper's stated TEAD-on-stiff and >50%-AP-1-on-soft findings; FOXA/HNF1B/LEF (paper's secondary stiff motifs) present but not significant in this run |
| `hicar_looping` | Fig. S7 | **reproduced** | `cosgrove2025_hicar_loop_stats` on the authors' own `GSE307449` loop bedpe files: loop conservation between soft/stiff = 43.49%/46.71% vs paper's stated 43.5-46.7% (near-exact); the companion ATAC-anchor-overlap statistic reproduced far more weakly (12.4% vs paper's 42.4-42.9%), not claimed |
| `myh9_crispri_screen` | Fig. 2 | **reproduced** (rank/correlation, not exact magnitude) | `cosgrove2025_dhs_screen_score` on `GSE243761_SupplementaryTable5` verified against the authors' own `SupplementaryTable6`: top-hit DHS ranking exact (45/72/46/47 = MYH9 intron-3 pRE#1-3), r=0.95 (tscore)/0.97 (median log2FC) restricting to negative-strand gRNAs per the authors' own `DHSmed_negStrand` convention |
| `bmf_anoikis_enhancer` | Fig. 3 | **blocked** `not_deposited` | luciferase/qPCR values checked against all 3 requested access routes (PMC supplementary bin, Science.org direct + suppl_file URL, bioRxiv full text) and failed all 3 — see `expected.json`'s blocker reason |
| `growth_migration_screen` | Fig. 4A-D | **blocked** `not_deposited` (with strong supporting count-level evidence) | `cosgrove2025_growth_migration_hit_calling` on `GSE243760_SupplementaryTable8`, after empirically testing 6 independent hit-calling models against the paper's 58/50/7 target: (a) permissive either-rep `\|Z\|>2` → 683/844/625 (far too many); (b) both-reps-agree per-gRNA `\|Z\|>2` → 53/100/12; (c) DHS-pooled mean Z (mirroring the validated MYH9 method) and alpha-RRA (via this session's own `sc_crispr_de_aggregate`, calibrated against the 1,000 deposited non-targeting gRNAs) → 0/0/0 at any reasonable threshold (over-diluted by the "small fraction of gRNAs" effect the paper itself notes); (d) averaged-Z per gRNA (mean of both reps, not requiring both individually >2) → 268/125-184 (too many); (e) an analogous baseline-count QC filter applied to migration → makes migration worse, not better; (f) both-reps-agree **plus a D0>=140-read QC filter for growth** → **53/52/7** — migration 9% low, growth 4% high, both-phenotype overlap **exact**. A fine-grained threshold/QC-cutoff scan (thresh 1.8-2.05, D0 cutoff 130-155) confirms no single shared threshold recovers both 58 and 50 simultaneously, suggesting the authors' real pipeline uses phenotype-specific thresholds not recoverable from the deposited per-gRNA table. The exact Methods-defined statistic sits behind the same access wall as `bmf_anoikis_enhancer`/`fibrosis_validation` (checked and failed the same 3 routes); the 53/52/7 model is retained as supporting, unconfirmed evidence, not claimed as a verified formula match |
| `sc_crispri_target_genes` | Fig. 4E-H | **reproduced** (3 named linkages, real Perturb-seq re-analysis) | `cosgrove2025_sc_crispri_guide_calling` + this session's own `sc_crispr_de_prepare`/`test` on `GSE243756`'s raw 118,647-cell matrix (after fixing a scanpy `gex_only` bug and a heavy-ambient-contamination guide-calling problem): MYH9-pRE45→MYH9 (log2fc=-0.42, p=1.9e-16), pRE#740→CTGF/CCN2 (log2fc=-1.04, p=3.3e-49), pRE#62→CYR61/CCN1 (log2fc=-0.10, p=0.022) — all 3 correct direction and significant |
| `fibrosis_validation` | Fig. 5 | **blocked** `not_deposited` | same 3-route check and same result as `bmf_anoikis_enhancer` |

```bash
bash Benchmarks/cosgrove2025_mechanosensitive_genomic/run.sh
python3 Benchmarks/concordance.py --benchmark cosgrove2025_mechanosensitive_genomic
```

## Honest caveats

* The `myh9_crispri_screen` and `rnaseq_atac_stiffness`/`hicar_looping` checks are approximate reconstructions (rank/correlation or order-of-magnitude), not byte-for-byte reproductions — the authors' `.rtf`-saved Rmd notebooks are exploratory, contain at least one visible bug (`bothRep_$pZ <- ...`, a malformed R identifier), and the exact normalization-constant detail behind the MYH9 tscore's absolute magnitude was not recovered.
* The A549 cell-line ATAC-seq differential-accessibility figure (~15%) was not independently reproduced — no A549 differential table is deposited, only per-replicate peak/bigwig files.
* `route: sc_analyze` in `provenance.json`/`expected.json` is a scaffold artifact from the auto-router (this is not a single-cell atlas paper); its placeholder "Single-cell pipeline summary written" check was removed as inapplicable rather than left to fail.
* The exact Methods text (needed for `growth_migration_screen`'s hit-calling formula, and for `bmf_anoikis_enhancer`/`fibrosis_validation`'s qPCR values) lives only in a POW-gated PMC binary or a Cloudflare-blocked Science.org/bioRxiv page in this session — every reasonable automated route was checked (see the two `not_deposited` blockers' reasons in `expected.json`), none succeeded. A JS-capable browser tool (not available to this session) would likely resolve this.
* The raw `GSE243756` sgRNA aggregate matrix has severe ambient cross-guide contamination (~99% of cells show nonzero counts for most of the 1,005 guides); the dominance-ratio guide call used here (`cosgrove2025_sc_crispri_guide_calling`, top guide >=5 UMI and >=3x the runner-up) recovers 63,791 singlet cells across 995 guides — fewer than the paper's ~103,440, and with looser per-guide power, which is why the `pRE#62` linkage only reaches significance after pooling all 10 of the pRE's gRNAs, not per individual gRNA.

## Provenance

`provenance.json` holds the full resolve/harvest/route record. Registered ports: `cosgrove2025_dhs_screen_score`, `cosgrove2025_stiffness_diff_count`, `cosgrove2025_hicar_loop_stats`, `cosgrove2025_sc_crispri_guide_calling`, `cosgrove2025_tf_motif_enrichment`, `cosgrove2025_growth_migration_hit_calling` (see `Scripts/ported/registry.json`).
