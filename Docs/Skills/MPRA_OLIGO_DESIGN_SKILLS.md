# MPRA Oligo Design — skill playbook

Python port of [MPRAOligoDesign](https://github.com/kircherlab/MPRAOligoDesign)
(Max Schubach, Berlin Institute of Health / kircherlab; MIT).
Cite Rosen *et al.*, *Genome Research* (2025), doi:10.1101/gr.281462.125.

The upstream Snakemake workflow needs conda, bedtools, pysam, pyranges,
vcfpy, pyfaidx and BioPython. This port is standard-library only, so a
design runs wherever IGVFagent runs.

## The four decisions in a design

| Stage | Question | Subcommand |
|---|---|---|
| Tiling | Where does the oligo window sit? | `tile` |
| Sequence | What do we write on the array? | `design-regions` / `design-variants` |
| Filtering | What breaks synthesis or the readout? | `filter` |
| Adapters | What constant flanks get bolted on? | `adapters` |

## Tiling strategies

Chosen automatically by region length:

| Region length | Strategy | Result |
|---|---|---|
| `<= --centering-max` | no tiling | one oligo, region centred and padded |
| `<= two_tiles_max` | two tiles | two oligos pinned to the region's ends |
| longer | centred tiling | centre oligo, then outwards with `>= --min-overlap` |

`two_tiles_max = min(--two-tiles-max, 2*oligo_length - 2*min_overlap - edge)`,
where `edge = 2*--variant-edge-exclusion` when `--include-variant-edge` is set.

## Filters

| Filter | Applies to | Default |
|---|---|---|
| Homopolymer run | every sequence | `> 10` fails |
| EcoRI / SbfI site | every sequence | any occurrence fails |
| Simple repeat | region-derived only | `> 25%` of the window fails |
| TSS overlap | region-derived only | any overlap fails |
| CTCF motif | region-derived only | any overlap fails |

The last three need annotation BEDs under `Data/Reference/OligoDesign/`
(`simpleRepeat.bed.gz`, `TSS_pos.bed.gz`,
`CTCF-MA0139-1_intCTCF_fp25.hg38.bed.gz`) or explicit `--simple-repeats`
/ `--tss-positions` / `--ctcf-motifs` paths. A missing file disables that
filter with a warning rather than failing the run.

If a REF oligo fails, its ALT partners are dropped with it — an ALT with
no REF to compare against is not interpretable.

## Usage

```bash
# Tile regions into oligo windows
igvfagent oligo tile --regions regions.bed.gz --oligo-length 200 \
    --min-overlap 50 --centering-max 200

# Region-only design
igvfagent oligo design-regions --regions regions.tiles.bed.gz \
    --reference hg38.fa --oligo-length 200

# Variant design: one REF + one ALT oligo per variant x region
igvfagent oligo design-variants --regions regions.bed.gz \
    --variants variants.vcf.gz --reference hg38.fa \
    --variant-edge-exclusion 20

# Filter, then add adapters
igvfagent oligo filter --design design.fa --regions regions.bed.gz \
    --map region_map.tsv.gz --max-homopolymer-length 10
igvfagent oligo adapters --design design_filtered.fa \
    --left AGGACCGGATCAACT --right CATTGCGTGAACCGA

# Everything in one run directory
igvfagent oligo pipeline --regions regions.bed.gz --variants variants.vcf.gz \
    --reference hg38.fa --oligo-length 200 --tile --label my_library
```

## Deviations from upstream

* **Homopolymer bug fixed.** Upstream's `nucleotideruns` never closes the
  final run, so `"ACCCC"` scores 1 and `"AAAA"` scores 0 — homopolymers at
  a sequence's end slip through. Fixed here.
* **Most-centred region fixed.** Upstream computes
  `abs(pos - 1 - start - length)`, the distance to the region *end*. The
  corrected centre distance is the default; `--legacy-centre-metric`
  restores upstream behaviour for bit-identical reproduction.
* No Snakemake / conda / bedtools / pysam / pyranges / vcfpy / pyfaidx /
  BioPython. Tabix lookups are replaced by an in-memory interval index.

## Outputs

Everything lands in `Docs/MPRAOligoDesign/<timestamp>_<label>/`:

| File | Contents |
|---|---|
| `design.fa.gz` | final oligos, ids namespaced `<sample>:<id>` |
| `regions.tiles.bed.gz` | oligo windows after tiling |
| `variant_region_map.tsv.gz` | Variant → Region → REF_ID → ALT_ID |
| `region_map.tsv.gz` | Region → sequence ID |
| `variants.vcf.gz` | placed variants, INFO annotated with Region/REF_ID/ALT_ID |
| `variants.removed.vcf.gz` | variants no oligo could carry |
| `filter.log.tsv` | one row per dropped sequence, with the reason |
| `summary.json` | counts for every stage |
