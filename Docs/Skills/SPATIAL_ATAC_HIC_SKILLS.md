# Spatial-ATAC-Hi-C skill

Spatially resolved co-profiling of 3D genome organisation and chromatin
accessibility, after Wang P., Wang J., Wang Q., Youngblood M. W. et al.,
*Nature Methods* 2026 (doi:10.1038/s41592-026-03217-4), GEO **GSE307620**.

The assay puts a 50 x 50 microfluidic barcode grid over a tissue
section: 2,500 spatial pixels, each with a Hi-C contact set *and* an
ATAC fragment set from the same molecules.

## What this skill consumes

Processed deposits only — the same boundary `share_seq` draws:

| Input | Flag | Notes |
|---|---|---|
| Contact table | `--pairs-dir` | 4DN `.pairs`, **or** GEO's tabix-indexed `*.hic.fragments.sorted.header.tsv.gz`; one file or a directory of per-pixel ones |
| ATAC fragments | `--fragments` | `chrom start end barcode count` |
| Spatial positions | `--positions` | AtlasXBrowser `*_tissue_positions_list.csv` |
| Gene model | `--gene-model` | GTF/GFF or BED6 |
| Chrom sizes | `--chrom-sizes` | optional if the pairs header has `#chromsize:` |

**Contact-file shapes.** A 4DN header (`#columns:` / `#chromsize:`), a
bare column-header row, and a headerless TSV all parse identically, and
column-name synonyms (`chr1`/`chrom1`, `position1`/`pos1`, ...) are
normalised. Pixel ids are taken from an `AAxBB` token anywhere in the
filename, so `12x47.hic.fragments.sorted.header.tsv.gz` keys as `12x47`
and joins with the cluster labels and the tissue grid. Use
`--pairs-glob` if a directory holds TSVs you do *not* want swept in.

What GSE307620 actually deposits, per sample (inside `GSE307620_RAW.tar`,
9.9 GB):

```
GSM9228169_MBR6merge.fragments.tsv.gz                    <- --fragments
GSM9228169_MBR6merge.hic.fragments.sorted.header.tsv.gz  <- --pairs-dir
GSM9228169_MouseBrainR6_tissue_positions_list.csv.gz     <- --positions
GSM9228169_MouseBrainR6_scalefactors_json.json.gz
GSM9228169_MouseBrainR6_tissue_{hires,lowres}_image.png.gz
```

Read alignment is **not** reimplemented. It is the one upstream step with
no tractable clean-room form, and the tools involved (Cell Ranger ATAC;
runHiC, GPL-3.0) are excluded from this Apache-2 codebase by policy.
Bring aligned pairs/fragments — GSE307620 publishes them.

## Typical run

```bash
# 0. See what the series deposits, then fetch the parts you want.
igvfagent spatial-hic pull-geo --gse GSE307620 \
    --download 'pairs|fragments|positions'

# 1. If you have one barcoded pairs file, split it into the pixel grid.
igvfagent spatial-hic pixel-demux --pairs sample.pairs.gz \
    --barcode-a barcodes_A.txt --barcode-b barcodes_B.txt \
    --layout BA --label mouse_R6

# 2. Per-pixel QC. The paper reports medians of 25,343-58,403 total
#    contacts, 88.1-90.3% cis, 24-33.3% long-range (>=10 kb over cis).
igvfagent spatial-hic qc --pairs-dir <run>/pixels \
    --fragments fragments.tsv.gz --gene-model gencode.gtf.gz \
    --label mouse_R6

# 3. Gene-level scores. GAS from ATAC (promoter + body), GAD from Hi-C
#    pair ends over the gene body.
igvfagent spatial-hic gas --fragments fragments.tsv.gz \
    --gene-model gencode.gtf.gz --label mouse_R6
igvfagent spatial-hic gad --pairs-dir <run>/pixels \
    --gene-model gencode.gtf.gz --label mouse_R6

# 4. Imputation, compartments, copy number.
igvfagent spatial-hic impute --pairs-dir <run>/pixels --chrom chr2 \
    --resolution 100000 --pad 1 --per-pixel --label mouse_R6
igvfagent spatial-hic compartment --pairs-dir <run>/pixels \
    --resolution 100000 --gene-model gencode.gtf.gz --label mouse_R6
igvfagent spatial-hic cnv --pairs-dir <run>/pixels \
    --resolution 5000000 --per-pixel --smooth --label tumour

# 5. Loops: quantify at known anchors, test across clusters, pile up.
igvfagent spatial-hic loops --pairs-dir <run>/pixels --bedpe loops.bedpe \
    --clusters clusters.tsv --apa --label mouse_R6

# 6. Put any per-pixel value back into tissue space.
igvfagent spatial-hic viz --table <run>/gene_activity_score.tsv \
    --column Satb2 --positions tissue_positions.csv --label mouse_R6
```

## Method notes

**Pixel demux.** The barcode is the two microfluidic rounds
concatenated. Upstream `bcsplit.py` writes barcode-B then barcode-A
(Read2 offsets 22:30 and 60:68), which is the `--layout BA` default.
Pass `--layout AB` if yours is the other way round; the `assigned`
fraction in `demux_summary.json` tells you immediately if you guessed
wrong.

**QC.** `long_range_ratio` uses cis as the denominator, matching the
paper's Fig. 1f ("fraction of long-range contacts over the unique
intra-chromosomal contacts"). `long_range_over_total` is reported
alongside so neither reading is ambiguous.

**TSS enrichment.** ArchR's definition, as the paper specifies:
mean insertion density in a 50 bp window on the TSS over the density in
the +/-1,900-2,000 bp flanks. Densities, not raw counts — the windows
differ in width by 4x.

**Imputation.** scHiCluster's convolution + random-walk-with-restart.
Note the restart term puts `rp` of the walk's mass on the diagonal
(over half of it at the default `rp=0.5`); pass `--zero-diagonal` when
whatever reads the matrix next cares about off-diagonal structure.
Imputation earns its keep at low depth: on a decay toy at 0.3x depth,
off-diagonal correlation with truth goes 0.43 -> 0.85.

**Compartments.** Observed/expected, Pearson correlation, leading
eigenvector. The eigenvector's sign is arbitrary, so pass
`--gene-model` to orient A toward gene-dense bins — otherwise A and B
can flip between chromosomes.

**Copy number.** Binned coverage scaled so the genome median is
`--ploidy` (2 by default; the paper treats every sample as diploid and
reports CN/2). GC and mappability corrections are a linear fit, not
NeoLoopFinder's Poisson GLM — concordant, not bit-identical. `--smooth`
applies MAGIC diffusion, which is what the paper uses before plotting
single-pixel CN in tissue space.

**Loops.** This quantifies loops at anchors you supply and tests them
across clusters (one-way ANOVA, `p < 0.05`, as in the paper). It does
**not** call loops de novo — Peakachu is a trained model, so bring its
BEDPE (or any other caller's).

## Outputs

Everything lands under `Docs/SpatialATACHiC/<timestamp>_<label>/`:
per-pixel TSVs, `.npz` matrix bundles, a `*_summary.json` per command,
and PNG+SVG figures under `Plots/`.

## Provenance and licensing

Apache-2.0. Algorithms are reimplemented from published descriptions;
no upstream source is copied, imported or vendored.

| Capability | Reference implementation | License | Approach |
|---|---|---|---|
| Assay + preprocessing | [wangjuan001/Spatial-ATAC-Hi-C](https://github.com/wangjuan001/Spatial-ATAC-Hi-C) | MIT | clean-room; barcode offsets and linker literals follow the published protocol |
| Contact-matrix build | [XiaoTaoWang/HiC_pipeline](https://github.com/XiaoTaoWang/HiC_pipeline) (runHiC) | **GPL-3.0** | external tool only — never imported. This skill starts from its `.pairs` output |
| Adapter trimming | [FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore) | **GPL-3.0** | external tool only — upstream of this skill's inputs |
| Imputation | scHiCluster (Zhou 2019 PNAS) | MIT | clean-room convolution + RWR |
| Compartments | cooltools `eigs-cis` | MIT | clean-room O/E + correlation eigenvector |
| CNV / segmentation | NeoLoopFinder `calculate-cnv`, `segment-cnv` | MIT | clean-room; linear bias fit in place of a Poisson GLM |
| Gene activity | SnapATAC2 `make_gene_matrix` | MIT | clean-room promoter+body insertion counting |
| Spatial smoothing | MAGIC (van Dijk 2018) | GPL-2.0 | clean-room diffusion — not imported |
| Loop calling | Peakachu (Salameh 2020) | MIT | **not** reimplemented; bring a BEDPE |

Runtime dependencies: `numpy`, `scipy`, `matplotlib` (all in the
`analysis` extra). No GPL runtime dependencies.
