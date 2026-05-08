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

## How this chains with other skills

- After `retrieve` / `manifest` / `download`, hand peak BEDs to `analyze-peaks` and `integrate-ccre`. For H3K27ac (or any enhancer mark), pipe to `super-enhancers` next.
- `browser` accepts any combination of BED tracks — including outputs from `splitseq_pipeline` / `multiome_10x_pipeline` ATAC peak BEDs, advanced-variant-analysis cCRE overlaps, etc.
- The agent runtime exposes `encode_retrieve`, `encode_describe`, `encode_super_enhancers`, `encode_integrate_ccre`, and `encode_browser` as tools so a single `igvfagent ask` can drive the whole pipeline end-to-end.