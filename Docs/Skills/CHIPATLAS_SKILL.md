# Skill: ChIP-Atlas canonical-query layer

A polite anonymous client for the [ChIP-Atlas](https://chip-atlas.org)
public REST + bulk archive — reprocessed ChIP-seq / ATAC-seq /
DNase-seq / Bisulfite-seq peak calls and metadata for hundreds of
thousands of public SRA experiments — plus the WABI Enrichment / Diff
Analysis job queue at NIG/DDBJ.

Clean-room reimplementation under Apache-2.0; stdlib-only. Modelled on
[inutano/chip-atlas](https://github.com/inutano/chip-atlas) (MIT,
Ohta/Oki/DBCLS 2017–2026) — we re-derive the wire contract; no source
copied. The skill addresses three hosts indirected through
``_endpoints.py``:

```
chipatlas_api   → https://chip-atlas.org             JSON browse / search / POST-download
chipatlas_data  → https://chip-atlas.dbcls.jp/data   bulk static archive
chipatlas_wabi  → https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas  Enrichment / Diff job queue
```

## Commands

```bash
# Browse — anonymous, JSON
igvfagent chipatlas list-genomes
igvfagent chipatlas list-qvalues
igvfagent chipatlas list-experiment-types --genome hg38
igvfagent chipatlas list-antigens --genome hg38 --ag-class Histone \
                                    --cell-class 'Pluripotent stem cell'
igvfagent chipatlas list-cell-types --genome hg38 --ag-class Histone

# Search + fetch metadata
igvfagent chipatlas search --query 'CTCF K562' --genome hg38 --limit 25
igvfagent chipatlas get-experiment SRX150531

# Per-experiment files (BigBed, BigWig, BED at q-thresholds 05/10/20)
igvfagent chipatlas download-experiment SRX150531 --genome hg38 \
                                                     --kinds bw,bb05,bb10
igvfagent chipatlas download-experiment SRX150531 --genome hg38 \
                                                     --kinds bw,bb05 --urls-only

# Assemble an all-peaks BED for a (ag × cellClass × qval) tuple
igvfagent chipatlas assemble-bed --genome hg38 \
                                   --ag-class 'TFs and others' --antigen CTCF \
                                   --cell-class 'Blood' --qval 05
# add --fetch to also stream the assembled BED to disk

# Bulk all-peaks (HEAD only by default — these are 1–30 GB files!)
igvfagent chipatlas download-all-peaks --genome hg38 --qval 05 --head-only
igvfagent chipatlas download-all-peaks --genome hg38 --qval 05 \
                                          --max-bytes 100000000   # 100 MB sample

# Pre-computed Target-Genes: discover which antigens have tables, then fetch
igvfagent chipatlas target-genes --genome hg38 --list
igvfagent chipatlas target-genes --genome hg38 --antigen H3K4me3 --distance 5000

# WABI Enrichment Analysis (over-representation of TF peaks at query regions)
igvfagent chipatlas submit-enrichment --mode genes --genome hg38 \
                                         --query gene_list.txt --background bg_list.txt \
                                         --label my_gene_enrichment_v1
igvfagent chipatlas poll-enrichment <JOB_ID> --interval 15 --timeout 1800

# One-command showcase
igvfagent chipatlas showcase --genome hg38 --cell-class 'Pluripotent stem cell'
```

## Conventions

| Concept | Allowed values |
|---|---|
| `genome` | `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` |
| `agClass` | `Histone`, `TFs and others`, `RNA polymerase`, `Input control`, `ATAC-Seq`, `DNase-seq`, `Bisulfite-Seq`, `Annotation tracks` |
| `qval` | `05` / `10` / `20` / `50` — peak-call `-log10(q)` threshold; lower digit → looser peaks |
| `kinds` (per-experiment) | `bw` (BigWig), `bb` (BigBed all-peaks), `bb05`/`bb10`/`bb20` (BigBed per-qval), `bed05`/`bed10`/`bed20` (BED per-qval) |
| Cell-type tree | 2-level: `clClass` → `clSubClass`, discoverable via `list-cell-types` |

## What this skill adds over IGVFagent's ENCODE/cCRE/enhancer skills

| Capability | Before | After |
|---|---|---|
| Public SRA/GEO ChIP-seq beyond ENCODE | ❌ | ✓ (hundreds of thousands of SRX) |
| Per-experiment BigBed at fixed q-thresholds | partial (ENCODE only) | ✓ |
| (ag × cellClass × qval) assembled BED | ❌ | ✓ |
| Bulk all-peaks_light archive | ❌ | ✓ (HEAD-probe by default) |
| Pre-computed Target-Genes tables | ❌ | ✓ |
| WABI Enrichment / Diff jobs | ❌ | ✓ |
| Cell-type-aware antigen browser | partial | ✓ (counts per slice) |

## Politeness defaults

The client paces itself at **1 request per second** by default to be
gentle with a small group's servers. Override via `CHIPATLAS_MIN_SECONDS`
(env var, accepts fractional seconds, set `0` only for genuinely-needed
bulk pulls).

## License posture

* **Code**: Apache-2.0 IGVFagent ⊃ MIT upstream (`inutano/chip-atlas`).
* **Data**: NBDC/DBCLS LSDB Archive license (typically CC-BY 4.0). We
  only **fetch** / **link**; never redistribute. Read
  `https://dbarchive.biosciencedbc.jp/en/chip-atlas/lic.html` before
  re-publishing downloaded BED/BigBed bytes.
* **Citation**: Zou A, Ohta T, Oki S. *Nucleic Acids Res.* 2024;
  doi:[10.1093/nar/gkae358](https://doi.org/10.1093/nar/gkae358).
  Oki S et al. *EMBO Rep.* 2018;
  doi:[10.15252/embr.201846255](https://doi.org/10.15252/embr.201846255).

## Pairs well with

- `igvfagent enhancer ...` — once a ChIP-Atlas BED slice is downloaded,
  feed it to IGVFagent's ABC enhancer–gene linkage skill.
- `igvfagent ccre ...` — annotate ChIP-Atlas peaks against ENCODE SCREEN
  cCREs to harmonise with the IGVF-Catalog convention.
- `igvfagent enrich ora` — independent statistical check of any
  Target-Genes table.
- `igvfagent catalog find-associations` — cross-reference TF identity
  → IGVF Catalog edges for downstream regulatory analysis.
