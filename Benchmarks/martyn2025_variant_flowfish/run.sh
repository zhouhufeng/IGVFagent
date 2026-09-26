#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
LABEL="martyn2025_variant_flowfish"
SIM_DIR="Data/Benchmarks/$LABEL"
mkdir -p "$SIM_DIR"

# ---- Step 1. Online: pull the IGVF Portal MeasurementSet manifest ----
.venv/bin/igvfagent flowfish pull-portal --limit 500 --label "$LABEL"

# ---- Step 2. Analytical-chain mechanics demo (synthetic 20-element screen)
# ---- NOT a reproduction of this paper: Martyn 2025's own data (Step 3
# below) is already scored by IGVF's uniform pipeline, so this chain is
# never actually run on it. This only demonstrates that IGVFagent's
# `flowfish` skill (a generic flow-sort + CRISPR-screen effect-scoring
# pipeline, unrelated in name to this paper's own "Variant-EFFECTS" method)
# can take raw guide×bin counts through MLE -> real-space -> per-element
# significance, the same general shape of analysis this paper's Methods
# describe, on synthetic input calibrated to a known truth.
.venv/bin/igvfagent flowfish simulate \
    --out-dir "$SIM_DIR" \
    --n-elements 20 --guides-per-element 5 \
    --knockdown-frac 0.5 --cells-per-guide 200 --seed 42

# Step 2a: per-guide MLE on (guide × bin) counts
.venv/bin/igvfagent flowfish estimate-effects \
    --counts "$SIM_DIR/counts.tsv" \
    --sortparams "$SIM_DIR/sortparams.tsv" \
    --label "${LABEL}_pipeline"

# Step 2b: real-space rescaling against negative controls
RAW=$(ls -t Docs/FlowFISH/*${LABEL}_pipeline_raw_effects.tsv | head -1)
.venv/bin/igvfagent flowfish real-space --input "$RAW" --label "${LABEL}_pipeline"

# Step 2c: per-element Mann-Whitney + Welch + BH-FDR
EFFECTS=$(ls -t Docs/FlowFISH/*${LABEL}_pipeline*real_space*.tsv | head -1)
.venv/bin/igvfagent flowfish score-elements --effects "$EFFECTS" --label "${LABEL}_pipeline"

# ---- Step 3. Real paper data: Martyn 2025's own published Variant-EFFECTS
# variant-effect tables, pulled directly from the IGVF Portal (uniform-
# pipeline output the paper's own Data Availability statement points at —
# not synthetic, not re-derived from raw reads). Accessions are fixed
# per-paper facts, not discovered generically:
#   IGVFDS5056OAGR -> IGVFFI4057VSBO  PPIF promoter   (GRCh38 variant effects)
#   IGVFDS5031MNRR -> IGVFFI4333XLOF  PPIF enhancer   (GRCh38 variant effects)
#   IGVFDS1824XDMU -> IGVFFI4854DWEG  IL2RA promoter  (GRCh38 variant effects)
TS="$(date +%Y%m%d_%H%M%S)"
REAL_DIR="$SIM_DIR/real_data/${TS}_${LABEL}"
mkdir -p "$REAL_DIR"
curl -sL "https://api.data.igvf.org/tabular-files/IGVFFI4057VSBO/@@download/IGVFFI4057VSBO.tsv.gz" \
    | gunzip -c > "$REAL_DIR/PPIF_promoter_GRCh38.tsv"
curl -sL "https://api.data.igvf.org/tabular-files/IGVFFI4333XLOF/@@download/IGVFFI4333XLOF.tsv.gz" \
    | gunzip -c > "$REAL_DIR/PPIF_enhancer_GRCh38.tsv"
curl -sL "https://api.data.igvf.org/tabular-files/IGVFFI4854DWEG/@@download/IGVFFI4854DWEG.tsv.gz" \
    | gunzip -c > "$REAL_DIR/IL2RA_promoter_GRCh38.tsv"
python3 Benchmarks/martyn2025_variant_flowfish/analyze_real_data.py "$REAL_DIR"

# ---- Step 4. Follow-up real paper data: PPIF splice-site edits (3
# untransfected-control replicate files, each a subset of the paper's own
# "three edits disrupting the splice donor motif" proof-of-concept, Fig.
# 1e/f) and the lentiMPRA (episomal) parallel arm at the PPIF promoter,
# compared against the endogenous PPIF-promoter numbers above (paper's own
# "Pearson's r=0.54" cross-assay claim):
#   IGVFDS8174BPPS -> IGVFFI0524YUIL   PPIF splice site, rep A (GRCh38)
#   IGVFDS8267JJHO -> IGVFFI2542METL   PPIF splice site, rep B (GRCh38)
#   IGVFDS7090INDM -> IGVFFI5097SDKA   PPIF splice site, rep C (GRCh38)
#   IGVFDS1376WOXJ -> IGVFFI2620KDMB   lentiMPRA reporter variants, PPIF promoter (GRCh38)
# (IGVFDS1003XTAF is the lentiMPRA library's barcode-to-element mapping only
# -- no scored variant effects -- and is not fetched here.)
curl -sL "https://api.data.igvf.org/tabular-files/IGVFFI0524YUIL/@@download/IGVFFI0524YUIL.tsv.gz" \
    | gunzip -c > "$REAL_DIR/PPIF_splice_repA_GRCh38.tsv"
curl -sL "https://api.data.igvf.org/tabular-files/IGVFFI2542METL/@@download/IGVFFI2542METL.tsv.gz" \
    | gunzip -c > "$REAL_DIR/PPIF_splice_repB_GRCh38.tsv"
curl -sL "https://api.data.igvf.org/tabular-files/IGVFFI5097SDKA/@@download/IGVFFI5097SDKA.tsv.gz" \
    | gunzip -c > "$REAL_DIR/PPIF_splice_repC_GRCh38.tsv"
curl -sL "https://api.data.igvf.org/tabular-files/IGVFFI2620KDMB/@@download/IGVFFI2620KDMB.tsv.gz" \
    | gunzip -c > "$REAL_DIR/PPIF_promoter_lentiMPRA_GRCh38.tsv"
python3 Benchmarks/martyn2025_variant_flowfish/analyze_followups.py "$REAL_DIR"

echo ""
echo "== Martyn 2025 Variant-EFFECTS benchmark complete =="
echo "Real-data summary: $REAL_DIR/summary.json"
echo "Generate figures with:"
echo "  .venv/bin/python Benchmarks/martyn2025_variant_flowfish/make_figures.py"
