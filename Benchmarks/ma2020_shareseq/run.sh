#!/usr/bin/env bash
# Ma 2020 SHARE-seq (Cell 183:1103) — reproduction from GEO GSE140203.
#
# Steps 1-3 are the original RNA cell-type check; steps 4-9 reproduce the
# paper's analyses with the ported commands (see README "What is ported"):
#   shareseq-dorc              peak-gene associations + DORCs (port of FigR runGenePeakcorr)
#   shareseq-species-mix       human/mouse species-mixing calls (Fig. 1B-D)
#   shareseq-chromatin-potential, shareseq-dorc-residuals  (Fig. 5H, Fig. 4C; from STAR Methods)
#   shareseq-gene-activity     Seurat v3 peak gene activity for computational pairing (Fig. S2N-S)
# The skin peak-gene step tests ~1M peak-gene pairs x 101 correlations over
# 34,774 cells: run it through sbatch (IGVF_SBATCH=1) on a cluster.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"; cd "$ROOT"
LABEL="ma2020_shareseq"
D="Benchmarks/_data/$LABEL"
IGVF="${IGVF:-$(command -v igvfagent)}"
PY="${PY:-$(dirname "$IGVF")/python}"
mkdir -p "$D"

# 1) Per-sample GEO supplementary files (not the 7.5 GB RAW tar).
GEO=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4156nnn
fetch() { local gsm=$1 f=$2; [ -s "$D/$f" ] || curl -sfL --retry 5 -C - -o "$D/$f" "$GEO/$gsm/suppl/$f"; }
fetch GSM4156597 GSM4156597_skin.late.anagen.counts.txt.gz
fetch GSM4156597 GSM4156597_skin.late.anagen.barcodes.txt.gz
fetch GSM4156597 GSM4156597_skin.late.anagen.peaks.bed.gz
fetch GSM4156597 GSM4156597_skin_celltype.txt.gz
fetch GSM4156608 GSM4156608_skin.late.anagen.rna.counts.txt.gz
fetch GSM4156593 GSM4156593_GM12878.3T3.hg19.counts.txt.gz
fetch GSM4156593 GSM4156593_GM12878.3T3.hg19.barcodes.txt.gz
fetch GSM4156594 GSM4156594_GM12878.3T3.mm10.counts.txt.gz
fetch GSM4156594 GSM4156594_GM12878.3T3.mm10.barcodes.txt.gz
fetch GSM4156604 GSM4156604_GM12878.3T3.rna.hg19.counts.txt.gz
fetch GSM4156605 GSM4156605_GM12878.3T3.rna.mm10.counts.txt.gz

# 2-3) Original RNA cell-type recovery check (writes the run dir).
$PY "$HERE/prep_input.py"
$IGVF share rna-qc --h5ad "$D/shareseq_skin_rna.h5ad" --label "$LABEL"
$PY "$HERE/make_figures.py"
RUN="$(ls -d Docs/SHAREseq/2*_${LABEL} | sort | tail -1)"

# 4) Species mixing (Fig. 1B-D).
$IGVF shareseq-species-mix \
  --atac-human-mtx "$D/GSM4156593_GM12878.3T3.hg19.counts.txt.gz" --atac-human-barcodes "$D/GSM4156593_GM12878.3T3.hg19.barcodes.txt.gz" \
  --atac-mouse-mtx "$D/GSM4156594_GM12878.3T3.mm10.counts.txt.gz" --atac-mouse-barcodes "$D/GSM4156594_GM12878.3T3.mm10.barcodes.txt.gz" \
  --rna-human "$D/GSM4156604_GM12878.3T3.rna.hg19.counts.txt.gz" --rna-mouse "$D/GSM4156605_GM12878.3T3.rna.mm10.counts.txt.gz" \
  --out "$RUN/species_mix"

# 5) Skin peak-gene associations and DORCs (Fig. 3; +/-50 kb, p < 0.05, >10 peaks).
#    Background peaks + TSS table come from the FigR reference run
#    (Data/ma2020/ref/run_figr_peakgene.R), so the port uses the same
#    chromVAR background set.
REF="${MA2020_REF:-Data/ma2020/ref/skin_full}"
$IGVF shareseq-dorc peakgene \
  --atac-mtx "$D/GSM4156597_skin.late.anagen.counts.txt.gz" --atac-barcodes "$D/GSM4156597_skin.late.anagen.barcodes.txt.gz" \
  --peaks "$D/GSM4156597_skin.late.anagen.peaks.bed.gz" --rna "$D/GSM4156608_skin.late.anagen.rna.counts.txt.gz" \
  --tss "$REF/mm10_TSS_figr.tsv" --bg-matrix "$REF/bg_peaks_chromvar.tsv" \
  --window 50000 --dorc-cutoff 10 --out "$RUN/peakgene"
$IGVF shareseq-dorc dorc-scores \
  --atac-mtx "$D/GSM4156597_skin.late.anagen.counts.txt.gz" --atac-barcodes "$D/GSM4156597_skin.late.anagen.barcodes.txt.gz" \
  --peaks "$D/GSM4156597_skin.late.anagen.peaks.bed.gz" --gene-peak "$RUN/peakgene/gene_peak_cor.tsv.gz" \
  --out "$RUN/peakgene/dorc_scores.tsv.gz"

# 6-7) Chromatin potential (Fig. 5H) and DORC residuals over pseudotime (Fig. 4C)
#      on the hair-follicle lineage cells.
awk 'NR>1 && $2>10 {print $1}' "$RUN/peakgene/dorc_rank.tsv" > "$RUN/peakgene/dorc_genes.txt"
HF="TAC-1,TAC-2,IRS,Medulla,Hair Shaft-cuticle.cortex"
$IGVF shareseq-chromatin-potential \
  --atac-mtx "$D/GSM4156597_skin.late.anagen.counts.txt.gz" --atac-barcodes "$D/GSM4156597_skin.late.anagen.barcodes.txt.gz" \
  --dorc-scores "$RUN/peakgene/dorc_scores.tsv.gz" --dorc-genes "$RUN/peakgene/dorc_genes.txt" \
  --rna "$D/GSM4156608_skin.late.anagen.rna.counts.txt.gz" --celltypes "$D/GSM4156597_skin_celltype.txt.gz" \
  --types "$HF" --progenitor "TAC-1,TAC-2" --differentiated "IRS,Medulla,Hair Shaft-cuticle.cortex" \
  --out "$RUN/chromatin_potential"
$IGVF shareseq-dorc-residuals \
  --atac-mtx "$D/GSM4156597_skin.late.anagen.counts.txt.gz" --atac-barcodes "$D/GSM4156597_skin.late.anagen.barcodes.txt.gz" \
  --dorc-scores "$RUN/peakgene/dorc_scores.tsv.gz" --dorc-genes "$RUN/peakgene/dorc_genes.txt" \
  --rna "$D/GSM4156608_skin.late.anagen.rna.counts.txt.gz" --celltypes "$D/GSM4156597_skin_celltype.txt.gz" \
  --progenitor "TAC-1,TAC-2" --terminals "IRS,Medulla,Hair Shaft-cuticle.cortex" \
  --out "$RUN/residuals"

# 8) Computational pairing (Fig. S2N-S): gene activity here, CCA label
#    transfer in Seurat (Data/ma2020/ref/run_seurat_label_transfer.R).
$IGVF shareseq-gene-activity \
  --atac-mtx "$D/GSM4156597_skin.late.anagen.counts.txt.gz" --atac-barcodes "$D/GSM4156597_skin.late.anagen.barcodes.txt.gz" \
  --peaks "$D/GSM4156597_skin.late.anagen.peaks.bed.gz" --gtf "${MM10_GTF:-Data/ma2020/annot/Mus_musculus.GRCm38.102.gtf.gz}" \
  --out "$RUN/gene_activity/skin_activity.mtx.gz"

echo "Run dir: $RUN"
echo "Score with: igvfagent bench score --paper-id $LABEL"
