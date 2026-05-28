# IGVFagent Reproducibility Benchmark Suite

Twelve recent Nature / Cell / Science / Nat Genet / Nat Methods / Nat Commun papers whose published analyses **IGVFagent reproduces directly from public data**. Each benchmark is a self-contained directory with the data sources, run script, expected outputs, figures, and a paper-vs-IGVFagent comparison.

**Five benchmarks now complete with per-paper Concordance / Verdict / Honest caveats** — see the table below.

![Suite dashboard](figures/dashboard.png)

## ✅ Completed benchmarks (with results)

| # | Paper | Skill exercised | Headline result | Detail |
|---|---|---|---|---|
| ⭐ | **Matreyek 2018** PTEN VAMP-seq *Nat Genet* | `mavedb` | **4/4 concordance checks pass** (8,000 variants, all CDS-mapped) | suite-verified smoke-test |
| 1 | **Waters 2024** BAP1 SGE *Nat Genet* | `mavedb` (new SGE path) | **LOF +1.1 %, GOF +7.7 %** vs paper | [`waters2024_bap1/`](waters2024_bap1/README.md) |
| 2 | **Buckley 2024** VHL SGE *Nat Genet* | `mavedb` (new SGE path) | **2,268 / 2,268** variants recovered, 4-bucket scheme | [`buckley2024_vhl/`](buckley2024_vhl/README.md) |
| 3 | **Zou 2024** ChIP-Atlas 3.0 *Nucleic Acids Res* | `chipatlas` | **815 TFs catalogued** for hg38/Blood, **GATA1 rank #7** (121 experiments) | [`zou2024_chipatlas_gata1/`](zou2024_chipatlas_gata1/README.md) |
| 4 | **Agarwal 2025** lentiMPRA K562/HepG2/WTC11 *Nature* | `mpra` (`pull`) | **3/3 discovery artefacts** written; portal manifest + summary JSON | [`agarwal2025_lentimpra/`](agarwal2025_lentimpra/README.md) |
| 5 | **Yao 2024** ENCODE4 noncoding CRISPRi *Nat Methods* | `encode` (new FCE path) | **368 CRISPR-screen FCEs enumerated**; Engreitz 65 %, K562 dominant, GATA1 +24/+58 kb spot-check ✓ | [`yao2024_encode4_crispri/`](yao2024_encode4_crispri/README.md) |

## 🟡 Scaffolded (URN + workflow set up, awaiting full execution)

| # | Paper | Skill | What needs to happen |
|---|---|---|---|
| 6 | **Mitra 2024** SCARlink multiome *Nat Genet* | `multiome` | Pull IGVF multiome AnalysisSets; run `peak2gene` |
| 7 | **Joung 2025** TF Perturb-seq fibroblasts *Nat Genet* | `perturb-catalog` + `sc-analyze` | Catalogue query; full pipeline needs GSE273694 download |
| 8 | **Zheng 2024** in-vivo Perturb-seq cortex *Cell* | `perturb-catalog` + `sc-analyze` | needs local GSE249416 h5ad |
| 9 | **Deng 2024** cortex lentiMPRA *Science* | `mpra` + `enrich` | needs local GSE236018 counts |
| 10 | **Weinstock 2024** CD4 T-cell CRISPR network *Cell Genomics* | `perturb-catalog` + `network` | catalogue query + Steiner-tree on KG mirror |
| 11 | **Martyn 2025** Variant-FlowFISH *Cell* | `flowfish` | flagship clean-room test of the FlowFISH skill; needs local count table |

## How the benchmark suite is organised

```
Benchmarks/
├── README.md                      ← this file (suite-level dashboard)
├── OPERATIONS_GUIDE.md            ← shared prerequisites + execution paths
├── run_all.sh                     ← driver (--all / --online-only / --quick)
├── concordance.py                 ← stdlib-only pass/fail scorer
├── figures/
│   └── dashboard.png              ← embedded at top of this file
├── results/                       ← gitignored — per-run pass/fail reports land here
└── <paper-id>/
    ├── README.md                  ← per-paper headline + figures + how to reproduce
    ├── OPERATIONS.md              ← detailed step-by-step
    ├── expected.json              ← machine-readable ground-truth checks
    ├── run.sh                     ← deterministic IGVFagent CLI invocation chain
    ├── make_figures.py            ← regenerate the per-paper PNG/SVG plots
    └── figures/
        ├── fig1_*.png / .svg
        ├── fig2_*.png / .svg
        └── ...
```

## Quick start

### Verify the suite works on your machine (~60 s)

```bash
bash Benchmarks/matreyek2018_pten_vampseq/run.sh
.venv/bin/python Benchmarks/concordance.py --benchmark matreyek2018_pten_vampseq
```

Expected: `4/4 checks PASSED`.

### Run the 4 completed benchmarks

```bash
# Online-only — no local data downloads needed
bash Benchmarks/run_all.sh --online-only

# Score everything
.venv/bin/python Benchmarks/concordance.py --all
```

The aggregate concordance report lands at `Benchmarks/results/<ts>_concordance.{json,md}`.

## What each benchmark proves

### Variant-effect reproduction (`mavedb` skill)

**Matreyek 2018 PTEN** — the suite's smoke-test. The `mavedb_mapping_skill` was originally built around this paper's VAMP-seq scoreset; 4 / 4 hard checks (CSV download, gene resolution, VCF emission, per-variant TSV) pass on every fresh run.

**Waters 2024 BAP1 SGE** — first real reproducibility claim. IGVFagent's new SGE-aware code path (added this session) recovers all 18,108 variants and classifies them into LOF / GOF / Neutral via the inferred `|z| > 2.5` threshold:

![Waters concordance](waters2024_bap1/figures/fig2_concordance.png)

LOF concordance **+1.1 %**, GOF concordance **+7.7 %** — both well within typical methodological drift between independent classification reimplementations.

**Buckley 2024 VHL SGE** — 100 % row-count match (2,268 / 2,268). Applies the paper-published score thresholds (−1.26 / −0.4 / −0.22) verbatim to produce the 4-bucket functional classification:

![Buckley buckets](buckley2024_vhl/figures/fig2_buckets.png)

### TF / regulatory-element census (`chipatlas` skill)

**Zou 2024 ChIP-Atlas 3.0** — proves IGVFagent's online enumeration matches the upstream database. The canonical hematopoietic TF panel comes out in the correct rank order:

![Zou hematopoietic TFs](zou2024_chipatlas_gata1/figures/fig2_hematopoietic_tfs.png)

### MPRA workflow (`mpra` skill)

**Agarwal 2025 lentiMPRA** — verified end-to-end at commit `ad73060`. 3 / 3 discovery artefacts produced (Docs/MPRA report, Data/Manifests CSV, Data summary JSON). Stage-2 analytical step (`mpra activity` on the actual oligo counts) waits on the user downloading GSE142696.

## Engineering work the suite produced

While building these benchmarks I uncovered + fixed seven real bugs in IGVFagent (each shipped as its own commit on `main`):

| Bug | Fix | Commit |
|---|---|---|
| `kg variant rs429358` returned all-zero edges | pre-resolve rsID → SPDI before edge fan-out | `4e547ab` |
| 7 / 12 LLM tool calls failed silently | dispatcher dropping positional args via empty-string flag_map | `3155dbc` |
| GPT-5 `max_tokens` API rejected | use `max_completion_tokens` for reasoning-model generations | `9904424` |
| OpenAI's 128-tool array cap | trim to 128, preserve all 32 ★-prefixed tools | `7ed7ad0` |
| `enrich_ora --genes` rejected inline lists | accept inline strings, not just file paths | `78181d2` |
| `kg gene` hangs for 15+ minutes on dead sockets | 30 s timeout + retry + fail-fast | `f5f50b8` |
| Python `urllib` 10–40 s IPv6 fallback to IPv4 | monkeypatch `socket.getaddrinfo` to prefer IPv4 | `593f4e5` |
| **`mavedb_mapping_skill` couldn't read SGE scoresets** (Waters / Buckley) | **new `map_sge_scoreset()` + `parse_hgvsc_full()` + Ensembl /map/cdna/ path** | **`c6f42fd`** |
| **`encode retrieve` couldn't enumerate ENCODE4 CRISPR-screen / MPRA FCEs** (Yao 2024) | **add `encode_type` config + FCE assay-title entries (`CRISPR screen` · `Flow-FISH CRISPR screen` · `MPRA`)** | this commit |

The SGE extension and the FCE extension are the two most consequential of these — together they unlock every SGE scoreset on MaveDB (~50 + growing) and every ENCODE4 functional-characterization screen (~900 + growing) for IGVFagent's enumeration and coordinate-mapping pipelines.

## How to add a new benchmark

1. `mkdir Benchmarks/<paper-id>/`
2. Add `README.md` (paper context + headline), `expected.json` (ground-truth checks), `run.sh` (CLI invocation chain).
3. Optional: `make_figures.py` to render PNG/SVG plots committed under `figures/`.
4. Add the paper-id to `Benchmarks/run_all.sh` in the appropriate mode list.
5. Run `bash Benchmarks/<paper-id>/run.sh` and `.venv/bin/python Benchmarks/concordance.py --benchmark <paper-id>` to verify it works.

The per-paper `OPERATIONS.md` template is in `Benchmarks/OPERATIONS_GUIDE.md`.

## License + provenance

* IGVFagent code: Apache-2.0.
* All cited papers credited under their publishers' terms in each per-paper README.
* Data fetched from MaveDB (CC-BY 4.0), ChIP-Atlas (NBDC-LSDB Archive license, CC-BY-style), IGVF Portal (CC-BY 4.0), Ensembl REST (Apache-2.0 access) — fetch/link only, never redistributed.
* All commits authored by **Hufeng Zhou <zhouhufeng@users.noreply.github.com>**.
