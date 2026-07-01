#!/usr/bin/env bash
# Open4Gene benchmark — Liu et al. Science 2025 (kidney multiome).
# (1) validate the IGVFagent port vs the R pscl::hurdle reference on the package
#     test data; (2) pull the paper's published Open4Gene links from Figshare and
#     check headline consistency.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
EXP="Data/Papers/open4gene/export"
PAPER="Data/Papers/open4gene/paper_results"
mkdir -p "$PAPER"

# 1) Pull the paper's Open4Gene significant links (public Figshare, md5-verified)
python3 Scripts/figshare_skill.py download --id 26299093 --file-id 49324033 --out-dir "$PAPER"

# 2) Export the Open4Gene R-package test data + R reference (needs R + pscl):
#    Rscript export_testdata_and_reference.R   (see Data/Papers/open4gene)

# 3) Run the IGVFagent port on the test data
python3 Scripts/open4gene_skill.py link \
  --rna "$EXP/RNA.mtx" --rna-genes "$EXP/RNA_genes.txt" \
  --atac "$EXP/ATAC.mtx" --atac-peaks "$EXP/ATAC_peaks.txt" \
  --meta "$EXP/Meta.csv" --pairs "$EXP/PeakGene.csv" \
  --covariates lognCount_RNA,percent.mt --celltype All --min-cells 5 \
  --out "$EXP/port_links.tsv"

# 4) Concordance + figures (see the inline analysis in the benchmark build)
echo "Open4Gene benchmark: see results/concordance.json + figures/"
