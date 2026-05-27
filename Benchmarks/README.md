# IGVFagent Reproducibility Benchmark Suite

A curated set of **11 recent Nature / Cell / Science / Nat Genet / Nat Methods / Nat Commun papers** whose published analyses IGVFagent can reproduce — or at least produce comparable summary statistics for — using only its own skill surface (no upstream pipeline installs required).

The suite is the artefact that turns "IGVFagent has 141 tools and a clean Apache-2.0 codebase" into a defensible scientific reproducibility claim: *given the same inputs as paper X, our skill chain produces concordant outputs.*

## Layout

```
Benchmarks/
├── README.md                            ← this file (suite index)
├── run_all.sh                           ← driver that executes every online-only benchmark
├── concordance.py                       ← reads each run's outputs, scores vs expected.json
├── results/                             ← timestamped concordance JSON + MD reports land here
└── <paper-id>/
    ├── README.md                        ← citation, accessions, workflow, IGVFagent skill mapping
    ├── expected.json                    ← ground-truth metrics + tolerance ranges
    └── run.sh                           ← the verified IGVFagent CLI invocations
```

Every `<paper-id>/run.sh` writes its output under `Docs/<skill>/<ts>_<label>/` (IGVFagent's standard per-run output convention). `concordance.py` then walks the latest run directory matching the label, compares against `expected.json`, and emits a row in the suite report.

## The 12 candidates

| # | Paper | Journal · Year | Skill | Mode |
|---|---|---|---|---|
| **0** | **Matreyek et al.** — PTEN VAMP-seq (suite smoke-test, **URN verified live**) | *Nat Genet* 2018 | `mavedb` | **online** |
| 1 | Waters et al. — BAP1 saturation genome editing *(URN needs verification)* | *Nat Genet* 2024 | `mavedb` | **online** |
| 2 | Buckley et al. — VHL saturation genome editing | *Nat Genet* 2024 | `mavedb` + `catalog` | **online** |
| 3 | Yao et al. — ENCODE4 noncoding CRISPRi screens | *Nat Methods* 2024 | `encode` + `crispri` | online + local |
| 4 | Mitra et al. — SCARlink multiome regression | *Nat Genet* 2024 | `multiome` | online + local |
| 5 | Agarwal et al. — lentiMPRA K562/HepG2/WTC11 | *Nature* 2025 | `mpra` | local |
| 6 | Joung et al. — TF Perturb-seq fibroblasts | *Nat Genet* 2025 | `perturb-catalog` + `sc-analyze` | **online** + local |
| 7 | Zheng et al. — in-vivo AAV Perturb-seq cortex | *Cell* 2024 | `perturb-catalog` + `sc-analyze` | local |
| 8 | Deng et al. — cortex lentiMPRA + psychiatric variants | *Science* 2024 | `mpra` + `enrich` | local |
| 9 | Weinstock et al. — CD4+ T cell CRISPR network | *Cell Genomics* 2024 | `perturb-catalog` + `network` | **online** |
| 10 | Martyn et al. — Variant-FlowFISH | *Cell* 2025 | `flowfish` | online + local |
| 11 | Zou et al. — ChIP-Atlas 3.0 GATA1 hematopoietic | *Nucleic Acids Res* 2024 | `chipatlas` + `enrich` | **online** |

**online** = exercises only IGVFagent + public REST APIs (no local inputs).
**local** = requires the user to have first downloaded the counts table / h5ad / BED.

The **cleanest 3 to run first** (per the in-chat analysis) are `waters2024_bap1`, `martyn2025_variant_flowfish`, and `mitra2024_scarlink`. The first is one anonymous MaveDB call with a crisp ground-truth table; the second is the most direct test of an IGVFagent clean-room absorption (the `flowfish` skill *is* Variant-FlowFISH); the third probes the multiome peak-to-gene pipeline against well-known PBMC enhancers.

## Running the suite

### One benchmark at a time

```bash
# Run a single paper's IGVFagent workflow
bash Benchmarks/waters2024_bap1/run.sh

# Score the result against expected.json
.venv/bin/python Benchmarks/concordance.py --benchmark waters2024_bap1
```

### Whole online-only suite

```bash
# Run every benchmark that doesn't require local data downloads
bash Benchmarks/run_all.sh --online-only

# Score everything
.venv/bin/python Benchmarks/concordance.py --all
```

### Including local-data benchmarks

For the four "local-input-required" benchmarks (Agarwal lentiMPRA, Zheng in-vivo Perturb-seq, Deng cortex MPRA, Martyn Variant-FlowFISH analyse-step), the per-paper `README.md` lists which file you need to download and where to place it. Each `run.sh` checks for the required input at a canonical path under `Data/Benchmarks/<paper-id>/` and aborts cleanly with instructions if it's missing.

## Concordance scoring

Each `expected.json` declares two kinds of ground truth:

* **Hard metrics** with tolerance ranges, e.g.:
  ```json
  "n_lof_variants": { "min": 5000, "max": 6300, "expected": 5665 }
  ```
* **Qualitative checks**, e.g.:
  ```json
  "top_term_must_be_in": ["Cell Cycle", "G2-M Checkpoint", "E2F Targets"]
  ```

`concordance.py` walks the latest output directory under `Docs/<skill>/2*_<label>/` matching the benchmark's `label`, opens the canonical artefact (usually a `summary.json` or `*_mapped.tsv`), and compares each declared metric. A pass-rate (`n_passed / n_total`) is reported per paper plus across the suite.

## Provenance

The full forensic trail of every benchmark run lands under three locations:

| Location | Content |
|---|---|
| `Docs/<skill>/<ts>_<label>/` | The skill's own per-run output (TSV / JSON / SVG / PNG) |
| `Docs/Agent/<ts>_<query>/transcript.json` | Per-agent-run tool-call transcript (when run through `igvfagent ask` / the Streamlit UI) |
| `Benchmarks/results/<ts>_concordance.{json,md}` | The pass/fail concordance table |

This is the auditable artefact that distinguishes IGVFagent from agent stacks that only return a chat-style PDF — every reproducibility claim is backed by a directory of timestamped artefacts.

## License

All eleven referenced papers and their associated public data are cited under their original publisher's terms. IGVFagent itself remains Apache-2.0. The benchmark plan files (this directory) are part of IGVFagent under the same Apache-2.0 license.
