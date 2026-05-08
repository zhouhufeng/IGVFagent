# Skill: ENCODE bulk-genomics pipeline

End-to-end retrieval, description, peak QC, super-enhancer calling, cCRE integration, and browser-style SVG visualization for the major ENCODE bulk assays:

`ChIP-seq`, `Histone ChIP-seq`, `ATAC-seq`, `DNase-seq`, `Hi-C`, `capture Hi-C`, `ChIA-PET`, `RNA-seq`, `MNase-seq`, `FAIRE-seq`, `CAGE`, `RAMPAGE`

## Subcommands

### `retrieve` — find experiments by assay / biosample / target

```bash
igvfagent encode retrieve --assay 'Histone ChIP-seq' \
    --target H3K27ac --biosample K562 --assembly GRCh38 \
    --limit 50 --label k562_h3k27ac
igvfagent encode retrieve --assay ATAC-seq --biosample 'liver'
igvfagent encode retrieve --assay Hi-C --biosample GM12878
```

### `manifest` — per-file inventory for one or more accessions

```bash
igvfagent encode manifest --accessions ENCSR000DUB --label demo
```

### `download` — pull files under a size cap

```bash
igvfagent encode download \
    --manifest Data/Manifests/ENCODE/<files.csv> \
    --max-download-gb 5 --formats bed bigWig
```

### `describe` — plain-language report for one experiment

```bash
igvfagent encode describe --accession ENCSR000DUB
```

### `analyze-peaks` — peak QC stats from a BED file

```bash
igvfagent encode analyze-peaks --bed peaks.bed.gz
```

Outputs: peak count, width distribution, score distribution, per-chromosome counts, plus a 3-panel QC PNG/SVG.

### `super-enhancers` — ROSE-style super-enhancer call

```bash
igvfagent encode super-enhancers --bed h3k27ac_peaks.bed \
    --stitching-distance 12500 --tss-bed tss.bed --tss-distance 2000 \
    --label k562_h3k27ac
```

Stitches enhancer-mark peaks within `--stitching-distance` (default 12.5 kb), optionally excludes TSS-proximal peaks given a TSS BED, ranks the stitched regions by summed signal, and finds the geometric inflection point. Emits two BED files (super- and typical-enhancers) plus a hockey-stick PNG/SVG.

Reliable for H3K27ac, BRD4, MED1, MED12, P300; not meaningful for polycomb / heterochromatin marks.

### `integrate-ccre` — overlay peaks with SCREEN cCREs

```bash
igvfagent encode integrate-ccre --bed peaks.bed
igvfagent encode integrate-ccre --bed peaks.bed \
    --ccre-bed Data/Cache/cCRELinkage/<custom.bed>
```

Auto-downloads the SCREEN V4 cCRE BED on first use (cached in `Data/Cache/ENCODE/`). Annotates each peak with its cCRE class (PLS / pELS / dELS / CTCF-only / DNase-H3K4me3) and emits a stacked-bar overview.

### `browser` — IGV-style multi-track SVG

```bash
igvfagent encode browser \
    --region chr19:44903000-44912000 \
    --track 'H3K27ac peaks:peaks.bed' \
    --track 'ATAC peaks:atac_peaks.bed' \
    --with-ccre \
    --label apoe_locus
```

Each `--track LABEL:PATH` adds a horizontal track strip; `--with-ccre` overlays SCREEN cCREs colored by class.

### `bigwig-frip` — fraction-of-signal-in-peaks (needs `[hic]` extras)

```bash
pip install 'igvfagent[hic]'
igvfagent encode bigwig-frip --bigwig signal.bw --bed peaks.bed
```

Returns the global FRiP plus per-chromosome breakdown. ENCODE's FRiP minimum is 0.01 for ChIP-seq; > 0.05 is good. Histone ChIP-seq tends to score lower because the signal is broader.

### `bigwig-tss-heatmap` — anchor-centered signal heatmap (needs `[hic]`)

```bash
igvfagent encode bigwig-tss-heatmap \
    --bigwig signal.bw --anchor-bed tss.bed \
    --window 4000 --bins 200 --max-anchors 5000 --label tss_h3k27ac
```

Aggregates bigWig signal over `--window` bp centered on each anchor (TSS BED, peak summits, etc.), sorts rows by total signal, and renders a heatmap + meta-profile in one figure.

### `hic-matrix` — Hi-C contact heatmap (needs `[hic]`)

```bash
pip install 'igvfagent[hic]'
igvfagent encode hic-matrix \
    --input ENCFF000XYZ.mcool --region chr19:44900000-45100000 \
    --resolution 10000 --balance --label apoe_hic
igvfagent encode hic-matrix \
    --input ENCFF111ABC.hic --region chr19:44900000-45100000 \
    --resolution 10000
```

Loads the contact matrix for a region from `.mcool` (cooler) or `.hic` (hic-straw), saves the raw matrix as `.npy`, and renders a `log1p`-scaled contact heatmap PNG/SVG.

### `hic-insulation` — Crane-style TAD-boundary score (needs `[hic]`)

```bash
igvfagent encode hic-insulation \
    --input ENCFF000XYZ.mcool --region chr19:44000000-46000000 \
    --resolution 10000 --window 200000 --balance \
    --boundary-threshold -0.3 --label apoe_ins
```

Computes the insulation score along the diagonal, identifies local minima below `--boundary-threshold` as candidate TAD boundaries, and emits a BED of boundaries plus a 1-D track plot. For loop-level calls, run a dedicated tool (HiCCUPS / Mustache / Peakachu) and feed the resulting bedpe into `loops-analyze`.

### `loops-analyze` — Hi-C / ChIA-PET / capture Hi-C loops

```bash
igvfagent encode loops-analyze --bedpe loops.bedpe \
    --peaks 'CTCF peaks:ctcf.bed' --peaks 'H3K27ac peaks:h3k27ac.bed' \
    --label gm12878_loops
```

Parses a `.bedpe` and reports total / intra / inter / median loop length plus optional anchor ↔ peak overlap stats per input peak set. Works on outputs from any loop caller (Mustache, HiCCUPS, MAPS, Fit-Hi-C, ChIA-PET clusters, etc.).

### `motif-enrichment` — TF motif enrichment in peak sequences

```bash
pip install 'igvfagent[motif]'
# Download a genome FASTA from UCSC (one-time):
#   curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O
igvfagent encode motif-enrichment \
    --bed peaks.bed --genome hg38.fa.gz \
    --top 2000 --score-cutoff 8 --label k562_h3k27ac_motifs
```

Scans peak sequences against a curated, bundled JASPAR-derived PFM set (CTCF, AP-1, GATA1, ETS, NFkB, STAT1, FOXA1, TP53, MYC, SP1/KLF) and reports log2 odds enrichment vs mononucleotide-shuffled background. For full TF coverage use HOMER `findMotifsGenome.pl` or MEME-ChIP.

## How this chains with other skills

- After `retrieve` / `manifest` / `download`, hand peak BEDs to `analyze-peaks` and `integrate-ccre`. For H3K27ac (or any enhancer mark), pipe to `super-enhancers` next.
- `browser` accepts any combination of BED tracks — including outputs from `splitseq_pipeline` / `multiome_10x_pipeline` ATAC peak BEDs, advanced-variant-analysis cCRE overlaps, etc.
- The agent runtime exposes `encode_retrieve`, `encode_describe`, `encode_super_enhancers`, `encode_integrate_ccre`, and `encode_browser` as tools so a single `igvfagent ask` can drive the whole pipeline end-to-end.