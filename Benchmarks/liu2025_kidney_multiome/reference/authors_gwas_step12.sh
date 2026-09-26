#!/usr/bin/env bash
# Reference for liu2025_gwas_loci steps 1-2: the authors' awk from
# hbliu/Kidney_Epi_Pri@552c461 eGFR_GWAS/Independent.Loci.eGFR.GWAS.sh (Step1,
# Step2) verbatim, except that P.value is field $11 in the deposited Figshare
# file (it was $10 in the 2022 input the script was written for).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
D="$ROOT/Data/liu2025"; R="$D/ref_authors"; mkdir -p "$R"
G="$D/figshare/eGFRcrea_GWAS_Multi.txt.gz"
zcat "$G" | awk '{if($11 == "P.value" || $11 < 5e-8) print $0}' > "$R/Sig5e8.txt"
head -n 1 "$R/Sig5e8.txt" > "$R/Sig5e8.NoMHC.txt"
cat "$R/Sig5e8.txt" | awk 'NR >1 {if (!(($2 == 6) && ($3 > 25000000) && ($3 < 35000000))) print $0}' >> "$R/Sig5e8.NoMHC.txt"
awk 'BEGIN{OFS="\t"; print "ID","P"} NR>1{print $1"_"$2":"$3, $11}' "$R/Sig5e8.NoMHC.txt" > "$R/ref_ID_P.tsv"
