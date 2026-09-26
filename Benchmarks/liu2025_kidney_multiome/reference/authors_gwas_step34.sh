#!/usr/bin/env bash
# Reference for liu2025_gwas_loci steps 3-4: the authors' Step3 formatting awk
# and Step4 (plink --cm-map + R merge loop) from hbliu/Kidney_Epi_Pri@552c461
# eGFR_GWAS/Independent.Loci.eGFR.GWAS.sh, verbatim except: the input .clumped
# is the port's plink run (identical flags), --bfile replaces --file (binary
# copy of the same 1000G EUR panel), cM maps are HapMap-II GRCh37 per
# chromosome, and the R loop bound 2:1310 (hard-coded for 2022) is nrow().
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
D="$ROOT/Data/liu2025"; L="$D/ref"; W="$D/ref_authors/step34"; mkdir -p "$W"
C="$D/out/gwas_loci_cli/plink/PlinkClumping.p5e8.r20.1.clumped"
cd "$W"
echo CHR$'\t'POS$'\t'SNP$'\t'P$'\t'Left$'\t'Right > clumped.formatted
awk 'NR>1{print $3"\t"$5"\t"$12}' $C | sed 's/(1)//g' | sed 's/NONE//g' | awk -F',' '{print $1"\t"$NF}' | sed 's/:/\t/g' |\
awk '{if($1 >0) print $1"\t"$2"\t"$1":"$2"\t"$3"\t"$5"\t"$7}' >> clumped.formatted
awk 'NR>1{print $3}' clumped.formatted > clumped.formatted.SNP.txt
"$L/plink" --bfile "$L/eur_chrpos" --extract clumped.formatted.SNP.txt --recode --cm-map "$L/gmap/genetic_map_chr@_combined_b37.txt" --out clumpedSNP > /dev/null
Rscript "$ROOT/Benchmarks/liu2025_kidney_multiome/reference/authors_gwas_step4.R"
