# liu2025_kidney_multiome — Operations

Shared prerequisites live in `Benchmarks/OPERATIONS_GUIDE.md`; this file only covers what is specific to this paper.

## 1. Inputs

| Input | Source | Where it lands |
|---|---|---|
| eGFRcrea GWAS (multi-ancestry, 2.2M), RASQUAL ASE (tubule, glomeruli), bASA, snASA, Open4Gene significant links + summary statistics, snATAC peak BED, Kidney Disease Genetic Scorecard | Figshare 26299093 (CC BY), md5-checked by `run.sh` | `Data/liu2025/figshare/` (~2.4 GB) |
| Authors' code | `hbliu/Kidney_Epi_Pri@552c461`, `hbliu/Open4Gene@6e7f36a` | `Data/liu2025/repo_*` |
| LD reference | plink 1.9; 1000 Genomes phase 3 EUR (MAGMA `g1000_eur`, GRCh37) re-keyed `CHR:POS`; HapMap-II GRCh37 genetic map (Eagle) | `Data/liu2025/ref/` (~1 GB) |
| R reference runs | R 4.3.3 (`module load R/4.3.3-fasrc01` on FASRC) with `pscl`, `progress`, `GenomicRanges`, `Matrix`, `dplyr`, `tidyr`; `pscl`/`progress` are installed into `Data/liu2025/Rlib` if missing | — |

Individual-level kidney data (CMDGA, AMP consortium sign-in) and per-cohort GWAS are controlled access and not used.

## 2. Run

On a compute node (clumping 102k variants against 1000G EUR and parsing the 700 MB GWAS take ~15 min and ~20 GB RAM):

```bash
sbatch -p <partition> -c 4 --mem 32G -t 2:00:00 --wrap "bash Benchmarks/liu2025_kidney_multiome/run.sh"
igvfagent bench score --paper-id liu2025_kidney_multiome
```

`run.sh` writes `Docs/Benchmark/<ts>_liu2025_kidney_multiome/summary.json`; port verifications land in `Data/liu2025/verify/<tag>/validation_vs_reference.json` and the paper-number recomputation in `Data/liu2025/verify/derived/derived_metrics.json`.

## 3. Pieces

| File | Role |
|---|---|
| `reference/authors_gwas_step12.sh` | authors' awk Steps 1-2 (P < 5e-8, MHC removal) — reference for `liu2025-gwas-loci` |
| `reference/authors_gwas_step34.sh`, `reference/authors_gwas_step4.R` | authors' Step 3 formatting + Step 4 cM merge — reference for `liu2025-gwas-loci` |
| `reference/authors_open4gene.R` | Open4Gene R package on its bundled 62,278-cell test slice — reference for `liu2025-open4gene` |
| `verify_ports.py` | builds the comparison tables and runs `igvfagent bench verify-port` for both ports |
| `verify_derived_tables.py` | recomputes the paper's Fig. 2-6 numbers from the deposited tables |

## 4. Troubleshooting

* `verify-port` names its output directory by the second; `verify_ports.py` sleeps 1 s between calls and copies each result to a stable path. Do not run two copies at once.
* The port `liu2025-open4gene` reproduces pscl's numerics (glm.fit starts, R's `vmmin` BFGS, `optimhess`, R 4.3 `dnbinom`). Count-component estimates for poorly identified fits (few expressing cells, quasi-separation, theta -> infinity) still stop at slightly different points than R; see README.
