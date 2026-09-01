# Enhancer / regulatory-region annotation — skill playbook

Annotates a list of **regions** against the ENCODE SCREEN candidate
cis-regulatory element (cCRE) registry. This is interval overlap, not
point lookup: `ccre annotate-variants` answers "what overlaps this
coordinate", while this answers "how much of this enhancer is covered by
a cCRE, of which classes, and if none — how far is the nearest".

## Data source

The registry is served from the Weng lab for
[SCREEN](https://screen.encodeproject.org):

| Registry | File |
|---|---|
| `V3` (default) | `Registry-V3/GRCh38-cCREs.bed` |
| `V4` | `Registry-V4/GRCh38-cCREs.bed` |

Downloaded once into `Data/Reference/cCRE/` and indexed into
`Data/Reference/cCRE/ccre.sqlite`, so annotation runs offline afterwards
and is reproducible against a pinned registry version. The endpoint is
resolved at runtime like every other archive — no URL is embedded in
source.

## Usage

```bash
# once — download and index the registry (~64 MB, ~2.3 M elements)
igvfagent enhancer build-db --registry V3
igvfagent enhancer db-stats

# look before you leap: which sheet holds coordinates?
igvfagent enhancer inspect --input Data/Input/EnhancerList/library.xlsx

# annotate
igvfagent enhancer annotate \
    --input Data/Input/EnhancerList/library.xlsx \
    --sheet High.sgRNA0820 --label high_enhancers
```

## Input formats

`.xlsx` / `.xlsm`, BED, TSV, CSV, optionally gzipped. Coordinate columns
are detected by name (`Enhancer_chr` / `chrom` / `chr`, `Enhancer_start` /
`start` / `chromStart`, and so on) and a headerless BED is recognised by
shape. Override with `--chrom-col` / `--start-col` / `--end-col` /
`--name-col` when detection is wrong.

`.xlsx` is read with `zipfile` + `xml.etree`, streaming, so no third-party
spreadsheet dependency is required and million-row sheets do not have to
be held in memory.

A CRISPR library lists one row per sgRNA, so the same enhancer repeats
several times. Rows are deduplicated on coordinates by default; pass
`--keep-duplicates` to annotate every row.

## Outputs

Written to `Docs/EnhancerAnnotation/<timestamp>_<label>/`:

| File | Contents |
|---|---|
| `annotated_regions.csv` | one row per region with the columns below |
| `regions_without_ccre.csv` | the subset with no overlap, for follow-up |
| `annotated_regions.bed` | BED with overlap fraction as score, class as name |
| `summary.json` | counts, class distribution, coverage statistics |

| Column | Meaning |
|---|---|
| `n_ccre` | number of cCREs overlapping the region |
| `top_class` | one representative class, ranked PLS > pELS > dELS > … |
| `ccre_classes` | every class seen across overlapping elements |
| `ccre_overlap_bp` | bases covered, **union** — overlapping elements are not double-counted |
| `ccre_overlap_frac` | that coverage as a fraction of the region |
| `nearest_ccre` / `nearest_distance_bp` | filled only when nothing overlaps |

## Interpreting the result

A region with no overlapping cCRE is not necessarily inert — it may be a
true negative, or the registry may simply lack coverage in that cell
type. `nearest_distance_bp` is what separates "sits in a cCRE desert"
from "just missed one", so both are reported rather than collapsing the
no-overlap set into a single count.

`top_class` is a convenience for summarising; a region overlapping both a
promoter-like and a distal-enhancer-like element is genuinely both, and
`ccre_classes` keeps that.
