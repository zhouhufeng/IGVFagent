# Agarwal 2025 — lentiMPRA across K562 / HepG2 / WTC11

## Citation

Agarwal V, Inoue F, ..., Ahituv N. *Nature* **639**: 411–420 (2025).
DOI: [10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9) · PMID: 39814889 · PMC11903340

## Data

* **Data**: ENCODE portal, per the paper's Data availability (PMC11903340): large-scale libraries K562 ENCSR382BVV, HepG2 ENCSR022GQD, WTC11 ENCSR244FWB; joint libraries K562 ENCSR203UFY, HepG2 ENCSR405QCT, WTC11 ENCSR336MKI; pilot libraries K562 ENCSR460LZI, HepG2 ENCSR463IRX. Public, no access request needed. (An earlier version cited GEO GSE142696; that series is Klein et al. 2020, a different MPRA study, and has no per-cell-line count files.)
* **GitHub**: [visze/sequence_cnn_models](https://github.com/visze/sequence_cnn_models) — CNN code
* **Pipeline**: MPRAflow v1.0.20 (Gordon 2020)

## Headline workflow (paper)

1. lentiMPRA on ~680K cis-regulatory elements in 3 cell lines.
2. log2(RNA/DNA) per oligo → MPRAflow standard activity calls.
3. Train a CNN sequence model on the activity calls, predict held-out enhancers.

## What IGVFagent reproduces

* **Online step**: `mpra portal-manifest` discovers IGVF Portal MPRA / STARR / reporter AnalysisSets.
* **Local step** (requires DNA + RNA per-oligo count table at `Data/Benchmarks/agarwal2025_lentimpra/oligo_counts.tsv`):
  * `mpra activity` → per-oligo NB GLM Wald test
  * `mpra qc` → replicate-concordance + barcodes/oligo
  * `mpra volcano` → 4-panel volcano

## Ground truth to spot-check

* ~41 % of oligos called active overall (paper: 41.7 %)
* Per-cell-type fractions: K562 ~41.3 %, HepG2 ~39 %, WTC11 ~33 %
* Manifest should enumerate at least a handful of IGVF MPRA AnalysisSets
