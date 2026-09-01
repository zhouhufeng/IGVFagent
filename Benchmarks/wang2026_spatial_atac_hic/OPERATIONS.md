# OPERATIONS — Wang 2026 Spatial-ATAC-Hi-C (`spatial-hic`)

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## Quick run (~20 s; no network, no downloads)

```bash
bash Benchmarks/wang2026_spatial_atac_hic/run.sh
python3 Benchmarks/concordance.py --benchmark wang2026_spatial_atac_hic
```

Expected: `18/18 checks PASSED` (+3 unconfirmed, by design).

Requires only `pip install 'igvfagent[analysis]'` — numpy, scipy,
matplotlib. Nothing else, and no optional Hi-C stack.

To additionally confirm the real GEO series is reachable (this leg *does*
need network):

```bash
bash Benchmarks/wang2026_spatial_atac_hic/run.sh --with-geo
```

## Steps `run.sh` performs

1. **`prep_input.py`** — synthesises `Benchmarks/_data/wang2026_spatial_atac_hic/`:
   a barcoded `sample.pairs.gz` (144,000 pairs over a 6×6 = 36-pixel grid),
   `fragments.tsv.gz` (120,132 ATAC fragments), two 8 bp barcode
   whitelists, `chrom.sizes`, a 4-gene GTF, a 2-loop BEDPE and clone
   labels. Deterministic (`SEED = 20260901`). Then **re-measures** the
   realised cis fraction, long-range ratio, CN fold and TSS density ratio
   off the emitted files and writes them to `truth.json`.
2. **`spatial-hic pixel-demux`** — splits the barcoded pairs into the
   spatial grid using `--layout BA` (barcode-B then barcode-A, matching
   upstream `bcsplit.py`).
3. **`spatial-hic qc`** — per-pixel cis / trans / long-range contacts plus
   ArchR-style TSS enrichment.
4. **`spatial-hic gas` / `gad`** — gene activity score from ATAC
   fragments (keyed by pixel, via the barcode whitelists) and
   gene-associated domain score from Hi-C pair ends.
5. **`spatial-hic compartment`** — A/B PC1 on chr2 at 200 kb, oriented by
   gene density.
6. **`spatial-hic cnv`** — pseudobulk CN at 500 kb with HMM segmentation,
   plus the per-pixel matrix with MAGIC smoothing.
7. **`spatial-hic loops`** — per-pixel loop strength, one-way ANOVA across
   clones, APA pileup.
8. **`spatial-hic viz`** — renders GENEX activity onto the 50×50 grid.
9. **`make_figures.py`** — scores recovery against `truth.json`, writes
   `concordance_metrics.json` into the QC run dir and a copy plus the
   figure into `figures/`.

## Where artefacts land

```
Docs/SpatialATACHiC/<ts>_wang2026_spatial_atac_hic_demux/
  pixels/                    36 per-pixel .pairs files
  demux_summary.json         assigned_fraction, pixels_with_data
Docs/SpatialATACHiC/<ts>_wang2026_spatial_atac_hic_qc/
  pixel_qc.tsv               per-pixel contact statistics
  qc_summary.json            medians + TSS enrichment
  concordance_metrics.json   recovery vs planted truth   ← scored artefact
  Plots/benchmark_recovery.png
Docs/SpatialATACHiC/<ts>_wang2026_spatial_atac_hic_{gas,gad}/
  gene_activity_score.tsv  |  gene_associated_domain_score.tsv
Docs/SpatialATACHiC/<ts>_wang2026_spatial_atac_hic_compartment/
  compartments.tsv           chrom/start/end/pc1/compartment
Docs/SpatialATACHiC/<ts>_wang2026_spatial_atac_hic_cnv/
  cnv_pseudobulk.tsv         cn_ratio + cn_segment per bin
  cnv_per_pixel.tsv          pixel × bin CN matrix (MAGIC-smoothed)
Docs/SpatialATACHiC/<ts>_wang2026_spatial_atac_hic_loops/
  cluster_specific_loops.tsv F, p, specific, top_cluster
  apa_matrix.tsv
```

## Troubleshooting

**`no run directory Docs/SpatialATACHiC/*_..._demux`** — `make_figures.py`
was run before `run.sh`. Run the whole script.

**Assigned fraction near 0 in a real dataset** — the barcode layout is
reversed. Retry `pixel-demux` with `--layout AB`. The benchmark always
uses `BA` because its generator writes that order.

**APA reports "corner background is zero"** — expected here. With only
two loops, the pileup has no background to normalise against, so no
enrichment ratio is defined and none is asserted.

**Different numbers after editing `prep_input.py`** — `truth.json` is
regenerated on every run, so the recovery ratios stay meaningful, but the
absolute values in `expected.json` (`expected` fields and the
paper-band checks) assume the committed generator. Re-tune the ranges if
you change the generator.
