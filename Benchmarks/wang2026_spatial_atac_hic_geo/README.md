# Wang et al. 2026 — Spatial-ATAC-Hi-C, on the paper's own GEO data

*Spatial chromatin architecture and accessibility co-profiling of mammalian
tissues*, Nature Methods (2026).
DOI `10.1038/s41592-026-03217-4` · PMID 42680843 · data **GSE307620**
Authors' code: https://github.com/wangjuan001/Spatial-ATAC-Hi-C

## Two benchmarks, on purpose

| | `wang2026_spatial_atac_hic` | `wang2026_spatial_atac_hic_geo` (this one) |
|---|---|---|
| Input | synthetic, planted signal, 6×6 grid | the paper's own GSE307620 deposit |
| Network | none | yes |
| Runtime | ~20 s | download-bound |
| Asks | "does the math recover known truth?" | "do we get the paper's numbers?" |

The synthetic one is the unit test and should stay fast and green. This one
is the reproduction, and it is allowed to be slow and to fail honestly.

## Chain, mapped to the paper's figures

| Paper | IGVFagent |
|---|---|
| 50×50 barcoded pixel grid (Fig. 1a,b) | `spatial-hic pixel-demux` |
| contacts / cis / long-range ≥10 kb, TSS enrichment (Fig. 1e–g) | `spatial-hic qc` |
| gene activity score, ATAC (Fig. 2b) | `spatial-hic gas` |
| gene-associated domain score, Hi-C (Fig. 2a) | `spatial-hic gad` |
| A/B compartment PC1 at 100 kb (Fig. 3a) | `spatial-hic compartment` |
| per-pixel CNV at 5 Mb, diploid baseline (Figs. 4, 5) | `spatial-hic cnv` |
| cell-type-specific loops, one-way ANOVA p<0.05 (Fig. 3c–j) | `spatial-hic loops` † |

† `loops` does **not** call loops de novo. The paper used Peakachu at 10 kb;
bring that BEDPE and a pixel→cluster table.

## Scored: 12 checks

The three that actually test reproduction are the QC bands the paper reports
across all samples:

- median total contacts per pixel ∈ [25,343 – 58,403] (Fig. 1e)
- median cis fraction ∈ [88.1% – 90.3%] (Extended Data Fig. 4c)
- median long-range ≥10 kb ratio ∈ [24% – 33.3%] (Fig. 1f)

plus structural invariants (2,500-pixel grid, 100 kb compartments, 5 Mb
diploid CNV) and the artefacts each leg must write.

## Recorded but not scored: 10 claims

Each carries its source sentence and why it cannot be scored yet.
`concordance.py` reports these under `unconfirmed` and never counts them, so
this benchmark cannot pass by grading its own text extraction. To close them:

- **Seurat clustering + differential testing** → 915 and 431 marker genes,
  1,006 and 1,731 specific open chromatin regions, 1,683 differential
  compartments (765 at PCC > 0.3)
- **scATAC label transfer vs BICCN GSE246791** → 16 cell-type identities
- **Peakachu de novo loop calling** → 268/556/41 (R6), 1,688/243/17 (R8),
  6,492 (62%) loops shared with in situ Hi-C
- **EagleC + Delly/Lumpy WGS validation** → 53/73 and 7/12 SVs supported
- **matched bulk ATAC control library** → SCC 0.895

## A caveat about harvesting

This paper is closed-access, so `bench harvest` sees only title + abstract
and finds **0 accessions**. GSE307620 is therefore set by hand in
`expected.json` from the Data availability statement. The routing still works
(the assay terms fire on the abstract), but do not expect the scaffold to
discover the deposit on its own for paywalled articles.

## Blocked (2026-09-26): the deposited Hi-C file isn't pairs data

Staged real data from GSM9228169 (`GSM9228169_MBR6merge.hic.fragments.sorted.header.tsv.gz`,
2.47 GB) and GSM9228170 (`GSM9228170_MBR8merge...`, same shape) and hit a
genuine data-integrity problem in the series, not a gap on our end:

* The file GEO labels as the Hi-C contact file has **5 columns** — `chrom,
  start, end, barcode, count` — structurally identical to a plain ATAC
  fragments file. `end - start` is ~150 bp across 200k+ rows sampled (Tn5
  insert size), never kb-to-Mb, so there is **no second locus** (no
  `chrom2`/`pos2`) anywhere in the file. A real 4DN pairs file needs both.
* Confirmed on **both** mouse samples, so it isn't a one-off corrupt
  download.
* GSM9228169's file (named `...MBR6merge...`) carries an internal
  `# id=MouseBrainR8` header — the wrong sample's id. Likely fallout from
  the series' own `[UPDATE on July 1, 2026] processed data files have been
  replaced with the merged version` note.
* IGVFagent's `spatial-hic` code (correctly) expects this filename shape to
  carry a bare column-header row per an earlier version of this deposit;
  that assumption no longer holds after the July 2026 replacement. Forcing
  the existing `--pairs` code path through this file would **silently
  misassign columns** (e.g. reading a 16bp barcode string as `chrom2`)
  rather than failing loudly — `run.sh` now checks the column count and
  exits 77 with a clear message instead.
* Raw reads exist — BioProject `PRJNA1322418`, SRX30428684 → `SRR35327743`
  for the R6 sample, 641M paired reads / ~189 GB bases (~69 GB gzipped).
  R8's raw reads aren't on ENA/SRA yet. But the authors' own repo
  ([wangjuan001/Spatial-ATAC-Hi-C](https://github.com/wangjuan001/Spatial-ATAC-Hi-C))
  publishes only the **ATAC** alignment step (`03.atac_alignment.sh`) plus
  two downstream tools that consume an already-built contact matrix
  (Higashi, scHiCluster) — there is no published Hi-C alignment/pairing
  recipe. Rebuilding pairs from raw reads would mean guessing at an
  undocumented protocol (this looks like a Tn5-tagmentation-based contact
  assay, not classic restriction-digest Hi-C), with no way to verify the
  reconstruction matches what the authors actually ran.

**Net effect:** the metadata check (#1) still scores; the 11 checks that
need real contact pairs (pixel-demux grid/assigned-fraction, QC bands, GAS,
GAD, compartments, CNV) cannot be computed from what GSE307620 currently
serves. This is the same kind of "not_deposited" outcome as several
`demultiplex2_stoeckius` analyses — documented rather than forced.

## Run

```bash
igvfagent bench pipeline --query "10.1038/s41592-026-03217-4" --execute
# or directly:
bash Benchmarks/wang2026_spatial_atac_hic_geo/run.sh
python3 Benchmarks/concordance.py --benchmark wang2026_spatial_atac_hic_geo
```

`bench pipeline --force` regenerates the scaffold and **overwrites the curated
`expected.json`** with the route's six structural checks. The paper-value
checks here were read out of the article by hand; re-apply them after any
forced rescaffold.
