# Reproducibility benchmarks

Published analyses that IGVFagent reproduces from public data. The up-to-date dashboard is [`Benchmarks/README.md`](../../Benchmarks/README.md).

IGVFagent ships a **reproducibility benchmark suite** in
[`Benchmarks/`](../../Benchmarks/README.md): recent Nature / Cell / Science /
Nat Genet / Nat Methods / Genome Biol papers whose published analyses IGVFagent
reproduces **directly from public data**. Each benchmark is a self-contained
directory — data sources, a deterministic `run.sh`, machine-readable
`expected.json` checks, regenerable figures, and a paper-vs-IGVFagent
`README.md` with a *Concordance / Verdict / Honest caveats* structure. A
stdlib-only scorer (`concordance.py`) turns each run into pass/fail checks.

The table below is a snapshot of the first 22 benchmarks. Newer ones (Tabula Sapiens 2.0, Spatial-ATAC-Hi-C, the scE2G K562 crowdsourced features, the E2G CRISPR benchmark, CRISPR-BEAN and others) are listed with their current results in [`Benchmarks/README.md`](../../Benchmarks/README.md).

| # | Paper | Skill exercised | Headline result |
|---|---|---|---|
| ⭐ | **Matreyek 2018** PTEN VAMP-seq *Nat Genet* | `mavedb` | 4/4 concordance checks pass (8,000 variants) |
| 1 | **Waters 2024** BAP1 SGE *Nat Genet* | `mavedb` (SGE path) | LOF +1.1 %, GOF +7.7 % vs paper |
| 2 | **Buckley 2024** VHL SGE *Nat Genet* | `mavedb` (SGE path) | 2,268 / 2,268 variants recovered |
| 3 | **Zou 2024** ChIP-Atlas 3.0 *Nucleic Acids Res* | `chipatlas` | 815 TFs catalogued; GATA1 rank #7 |
| 4 | **Agarwal 2025** lentiMPRA *Nature* | `mpra` | 3/3 discovery artefacts written |
| 5 | **Yao 2024** ENCODE4 CRISPRi *Nat Methods* | `encode` (FCE path) | 368 CRISPR-screen FCEs enumerated |
| 6 | **Mitra 2024** SCARlink multiome *Nat Genet* | `multiome` | 505 multiome AnalysisSets; 4/4 content types |
| 7 | **Weinstock 2024** CD4 CRISPR *Cell Genomics* | `perturb-catalog` + `geo` | 1,197 CRISPR-screen datasets; 99.7 % CRISPRn |
| 8 | **Zheng 2024** in-vivo Perturb-seq *Cell* | `geo` + `sc-analyze` | GSE249416 metadata + 9-file inventory |
| 9 | **Martyn 2025** Variant-FlowFISH *Cell* | `flowfish` | end-to-end chain: 20 elements → 7 Significant |
| 10 | **Joung 2025** TF Perturb-seq *Nat Genet* | `perturb-catalog` | modality scale confirmed (15 datasets) |
| 11 | **Deng 2024** cortex lentiMPRA *Science* | `mpra` + `synapse` | 166-node Synapse walk; 12/12 annotations recovered |
| 12 | **Gschwind 2023 / Sheth 2024** ENCODE-rE2G & scE2G *bioRxiv* | `portal` + `synapse` | ENCODE-rE2G base AUPRC 0.634 reproduces the published 0.634 exactly |
| 13 | **Liu 2025** kidney multiome scorecard *Science* | `open4gene` + `figshare` | port matches R `pscl::hurdle` β to r=1.0; 125,699 unique peaks (exact) |
| 14 | **Travaglini 2020** lung atlas *Nature* | `sc-analyze` (**full local repro**) | 49 clusters vs 46 author types, AMI 0.81; 9/10 markers — 8/8 |
| 15 | **Trevino 2021** cortex multiome *Cell* | `multiome peak2gene` (**full local repro**) | all 26 lineage genes linked; 93 % positive cis links — 4/4 |
| 16 | **deMULTIplex2 / Stoeckius 2018** cell hashing *Genome Biol* | `multiseq` (**full local repro**) | all 8 donor HTOs; 83 % singlets — 7/7 |
| 17 | **Rosenberg 2018** SPLiT-seq CNS *Science* | `splitseq` (**full local repro**) | 156,049-nucleus atlas (exact); all 8 CNS lineages — 4/4 |
| 18 | **Ma 2020** SHARE-seq mouse skin *Cell* | `share` (**full local repro**) | 34,774-cell skin set (exact); all 23 author cell types; Leiden AMI 0.63 — 5/5 |
| 19 | **Wang 2025** developing neocortex multiome *Nature* | `sc-analyze` (**full local repro**) | 232,328-nucleus atlas (exact); 29 author cell types; AMI 0.65 — 7/7 |
| 20 | **Zou/Shi 2026** scEPS GWAS × single-cell *medRxiv* | `sceps` port (**full local repro**) | microglia significantly AD-associated (d=4.76e-5, **Z=4.99, P=6e-7**); d>0 fraction 58 % vs paper 58 % |
| 21 | **Rosen 2025** MPRAsnakeflow *Genome Res / bioRxiv* | `mpraflow` + `mpralib` (**full local repro**) | **210,660 / 210,660 values byte-identical** to the authors' published IGVF artefact; 3/3 stated complexity figures reproduced to the digit — 8/8 |

Benchmarks 14–21 are **full local reproductions** — they download the public
data and run IGVFagent's real single-cell / multiome / GWAS-integration chain,
then score concordance against the authors' own results. The newest, **scEPS**
(Zou/Shi 2026), reimplements the full estimate→cluster→aggregate pipeline
(validated vs upstream, corr 1.0) and reproduces the paper's microglia AD-
association from public SEA-AD + Bellenguez-AD MAGMA scores. Most other
benchmarks run end-to-end as pure online calls; a few (e.g. **Zheng 2024**,
**Deng 2024**) need a user-fetched local file to complete their full chains.

```bash
# Verify the suite works (~60 s)
bash Benchmarks/matreyek2018_pten_vampseq/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark matreyek2018_pten_vampseq

# Run all online benchmarks and score them
bash Benchmarks/run_all.sh --online-only
.venv/bin/python Benchmarks/concordance.py --all
```

---

[← Documentation index](README.md) · [Project README](../../README.md)
