# OPERATIONS — [CORRECTED] TF Perturb-seq fibroblasts

> **Correction notice.** This benchmark was originally scaffolded under a
> fabricated citation ("Joung J...Zhang F., *Nat Genet* 57:828–838") with
> fabricated data accessions (`GSE237056`, `SCP2169`). The real paper behind
> this DOI is Southard et al. 2025 (*Nat Genet* 57:2323–2334,
> doi:10.1038/s41588-025-02284-1); see `README.md` for the full correction
> and real identifiers. **For the real reproduction, including a genuine
> data-verified confirmation of the paper's guide-library counts, use
> [`Benchmarks/southard2024_comprehensive_transcription/OPERATIONS.md`](../../southard2024_comprehensive_transcription/OPERATIONS.md).**

For shared prerequisites, see `Benchmarks/OPERATIONS_GUIDE.md`.

## What this directory still does

`run.sh` performs only the generic, paper-agnostic Perturbation Catalogue
census check — this part was never fabricated, only the paper metadata
attached to it was:

```bash
bash Benchmarks/southard2024_comprehensive_transcription/perturb_catalog/run.sh
# part of Benchmarks/southard2024_comprehensive_transcription; not scored on its own
```

```bash
.venv/bin/igvfagent perturb-catalog summary
.venv/bin/igvfagent perturb-catalog search-modality \
    --modality perturb-seq --query KLF4 --dataset-limit 20 || true
```

Hits the Perturbation Catalogue API and asks: "show me every Perturb-seq
dataset where KLF4 was perturbed." This is a generic, live discovery step;
it makes no paper-specific claim.

## Where artefacts land

```
Docs/Perturbation/<ts>_summary/
└── summary.json          ← scored
Docs/Perturbation/<ts>_search_perturb-seq_KLF4/
└── search.json           ← informational; upstream often times out
```

## Concordance interpretation

| # | Check | What it verifies |
|---|---|---|
| 1 | `summary.json` exists, non-empty | Perturbation Catalogue API returned its census (structural check only — asserted by the route, not by the paper). Scored on the summary rather than the KLF4 search because that endpoint times out intermittently and the catalogue CLI takes no `--label`. |

## Real data and the real reproduction

The paper's own analysis code, real accessions, and a genuine
data-verified confirmation of two headline numbers (10,979 guides; 1,836
TFs targeted, both measured directly from the authors' Zenodo-deposited
guide×droplet UMI matrix rather than trusted from prose) live at:

```
Benchmarks/southard2024_comprehensive_transcription/
├── README.md              — corrected citation + what was verified
├── OPERATIONS.md           — how to reproduce
├── expected.json           — paper claims, 2 now confirmed against real data
├── run.sh
└── verify_guide_library.py — the actual verification script
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `.venv/bin/igvfagent: bad interpreter` | Repo-local `.venv` was built on a different machine (its shebang is baked to that machine's path) | Use `$(command -v igvfagent)` instead, or run `bash Benchmarks/southard2024_comprehensive_transcription/run.sh`, which already falls back to `$PATH` |
| Perturbation Catalogue API 503/504 | Upstream Cloud Run cold start | Retry once; `run.sh` already wraps the modality search in `|| true` |

## License + provenance

* **Paper code** (real): [norman-lab-msk/TFs_CRISPRa](https://github.com/norman-lab-msk/TFs_CRISPRa) — pinned by `igvfagent paper-code fetch` under `Data/PaperCode/norman-lab-msk__TFs_CRISPRa/`.
* **IGVFagent code**: Apache-2.0; `Scripts/perturbation_catalog_skill.py`.
* **Figure-generation script**: `make_figures.py` in this directory (unchanged; generic catalogue-summary figures only).
