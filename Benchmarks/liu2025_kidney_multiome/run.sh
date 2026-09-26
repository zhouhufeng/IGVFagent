#!/usr/bin/env bash
# liu2025_kidney_multiome — Kidney multiome-based genetic scorecard (Liu et al.,
# Science 2025, doi:10.1126/science.adp4753)
#
# What this reproduces, and from what:
#   * the authors' deposited derived data (Figshare 26299093, CC BY): multi-ancestry
#     eGFRcrea GWAS, RASQUAL ASE / bASA tables, snASA, Open4Gene links and summary
#     statistics, the snATAC peak set and the Kidney Disease Genetic Scorecard;
#   * two methods absorbed into IGVFagent as ports and verified against the
#     authors' own code (igvfagent port list):
#       liu2025-gwas-loci   <- hbliu/Kidney_Epi_Pri@552c461 eGFR_GWAS/Independent.Loci.eGFR.GWAS.sh
#       liu2025-open4gene   <- hbliu/Open4Gene@6e7f36a R/Open4Gene.R
# Individual-level kidney data (CMDGA, AMP sign-in) and the per-cohort GWAS are
# controlled access and are not used (see expected.json, analysis
# rerun_from_primary_data).
#
# Heavy steps (700 MB GWAS, 1000G EUR clumping) take ~15 min and ~20 GB RAM:
# run on a compute node (sbatch / salloc), not a login node. Re-runs reuse the
# downloads under Data/liu2025/.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
export IGVF_ROOT="$ROOT"
LABEL="liu2025_kidney_multiome"
IGVF=".venv/bin/igvfagent"
"$IGVF" --help >/dev/null 2>&1 || IGVF="$(command -v igvfagent)"
PY="$(dirname "$(readlink -f "$(command -v "$IGVF")")")/python"
[ -x "$PY" ] || PY="$(command -v python3)"
D="Data/liu2025"
mkdir -p "$D/figshare" "$D/ref" "$D/out"

# ---- 1. Authors' deposited data (Figshare 26299093), md5-checked
fetch() {  # id name md5
  local f="$D/figshare/$2"
  if [ ! -s "$f" ] || [ "$(md5sum "$f" | cut -d' ' -f1)" != "$3" ]; then
    curl -sfL -o "$f" "https://ndownloader.figshare.com/files/$1"
    [ "$(md5sum "$f" | cut -d' ' -f1)" = "$3" ] || { echo "md5 mismatch: $2" >&2; exit 1; }
  fi
}
fetch 47668987 eGFRcrea_GWAS_Multi.txt.gz   9a53ad9b896052b0ac6413c423a74d62
fetch 49324036 ASE_Glomeruli.txt.gz         90f4e4c6365d5d68d0ce338006a0bf69
fetch 49324039 ASE_Tubule.txt.gz            4c545db7492bd2813c4618d84bf729b7
fetch 49324030 bASA.txt.gz                  92b26bf1f121fdf74af2a4e082a97110
fetch 49324033 Open4Gene_sig.txt.gz         f39c585f83c2a486ac46b1915770f95a
fetch 49324063 Open4Gene_all.txt.gz         77bf478f7fe3d264eb9c996cdba5735a
fetch 49324027 snASA.txt.gz                 1f595acdd0a89573eb0308ed91ed949b
fetch 49817118 snATAC_peaks.bed.gz          3f784cb5e94ec255758013937af68c16
fetch 57066494 Scorecard.xlsx               1612a96bad0e11b17169552d64f80f8f

# ---- 2. Authors' code, pinned
[ -d "$D/repo_Kidney_Epi_Pri" ] || git clone -q https://github.com/hbliu/Kidney_Epi_Pri.git "$D/repo_Kidney_Epi_Pri"
git -C "$D/repo_Kidney_Epi_Pri" checkout -q 552c461f85b91f0fe888cf4d8448f6906bea5042
[ -d "$D/repo_Open4Gene" ] || git clone -q https://github.com/hbliu/Open4Gene.git "$D/repo_Open4Gene"
git -C "$D/repo_Open4Gene" checkout -q 6e7f36aa80e6ff778737fc25978337342394cc68

# ---- 3. LD reference: plink 1.9, 1000 Genomes phase 3 EUR (MAGMA build, GRCh37)
#         re-keyed CHR:POS as in the authors' script, HapMap-II GRCh37 genetic map
R_="$D/ref"
[ -x "$R_/plink" ] || { curl -sfL -o "$R_/plink.zip" https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && unzip -oq "$R_/plink.zip" plink -d "$R_"; }
[ -s "$R_/g1000_eur.bed" ] || { curl -sfL -o "$R_/g1000_eur.zip" "https://vu.data.surf.nl/index.php/s/VZNByNwpD8qqINe/download" && unzip -oq "$R_/g1000_eur.zip" -d "$R_"; }
if [ ! -s "$R_/eur_chrpos.bed" ]; then
  awk 'BEGIN{OFS="\t"}{$2=$1":"$4; print}' "$R_/g1000_eur.bim" > "$R_/eur_tmp.bim"
  ln -sf g1000_eur.bed "$R_/eur_tmp.bed"; ln -sf g1000_eur.fam "$R_/eur_tmp.fam"
  "$R_/plink" --bfile "$R_/eur_tmp" --make-bed --out "$R_/eur_chrpos" > /dev/null
fi
[ -s "$R_/genetic_map_hg19_withX.txt.gz" ] || curl -sfL -o "$R_/genetic_map_hg19_withX.txt.gz" \
  https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz
if [ ! -s "$R_/gmap/genetic_map_chr1_combined_b37.txt" ]; then
  mkdir -p "$R_/gmap"
  "$PY" -c "
import pandas as pd
m = pd.read_csv('$R_/genetic_map_hg19_withX.txt.gz', sep=' ')
for c, g in m.groupby('chr'):
    g[['position', 'COMBINED_rate(cM/Mb)', 'Genetic_Map(cM)']].to_csv(f'$R_/gmap/genetic_map_chr{c}_combined_b37.txt', sep=' ', index=False)"
fi

# ---- 4. Fig. 1A-B: independent loci with the absorbed port, then the authors' code as reference
"$IGVF" liu2025-gwas-loci --gwas "$D/figshare/eGFRcrea_GWAS_Multi.txt.gz" \
    --bfile "$R_/eur_chrpos" --plink "$R_/plink" \
    --genetic-map "$R_/genetic_map_hg19_withX.txt.gz" --out "$D/out/gwas_loci_cli"
module load R/4.3.3-fasrc01 2>/dev/null || true
bash "$HERE/reference/authors_gwas_step12.sh"
bash "$HERE/reference/authors_gwas_step34.sh"

# ---- 5. Fig. 5A: Open4Gene — authors' R package (reference) and the port, same input
Rscript -e 'lib <- file.path(Sys.getenv("IGVF_ROOT"), "Data/liu2025/Rlib"); dir.create(lib, FALSE, TRUE);
  .libPaths(c(lib, .libPaths())); for (p in c("pscl", "progress")) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, lib = lib, repos = "https://cloud.r-project.org", quiet = TRUE)'
Rscript "$HERE/reference/authors_open4gene.R" > /dev/null
O4G="$D/ref_authors/open4gene"
mkdir -p "$D/out/open4gene"
for CT in All Each; do
  "$IGVF" liu2025-open4gene --rna "$O4G/rna.mtx" --rna-genes "$O4G/rna_genes.txt" \
      --atac "$O4G/atac.mtx" --atac-peaks "$O4G/atac_peaks.txt" --meta "$O4G/meta.csv" \
      --pairs "$O4G/pairs.csv" --covariates lognCount_RNA,percent.mt --celltype "$CT" \
      --min-cells 5 --out "$D/out/open4gene/port_$CT.res.txt" > /dev/null
done

# ---- 6. Port-vs-authors verification (class B) and paper numbers from the deposited tables (class A)
"$PY" "$HERE/verify_ports.py"
"$PY" "$HERE/verify_derived_tables.py" > /dev/null

# ---- 7. Run directory for concordance.py
RUN="Docs/PaperReproduction/$(date +%Y%m%d_%H%M%S)_${LABEL}"
mkdir -p "$RUN"
"$PY" - "$RUN" <<'EOF'
import json, sys
from pathlib import Path
run = Path(sys.argv[1])
s = {"paper_id": "liu2025_kidney_multiome",
     "gwas_loci": json.loads(Path("Data/liu2025/out/gwas_loci_cli/summary.json").read_text()),
     "derived": json.loads(Path("Data/liu2025/verify/derived/derived_metrics.json").read_text()),
     "port_verification": {p.parent.name: json.loads(p.read_text()).get("validation_vs_reference", {})
                           for p in sorted(Path("Data/liu2025/verify").glob("*/validation_vs_reference.json"))}}
(run / "summary.json").write_text(json.dumps(s, indent=2))
print(f"wrote {run}/summary.json")
EOF

echo ""
echo "== liu2025_kidney_multiome benchmark complete =="
echo "Score with: $IGVF bench score --paper-id $LABEL"
