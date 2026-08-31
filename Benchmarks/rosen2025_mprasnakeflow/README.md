# rosen2025_mprasnakeflow

Reproduction benchmark for:

> Rosen JD, Vasanthakumari AD, Salomon K, de Lange N, Dash PM, Keukeleire P,
> Hassan A, Barrera A, Kircher M, Love MI, Schubach M.
> **Uniform processing and analysis of IGVF massively parallel reporter assay
> data with MPRAsnakeflow.** bioRxiv 2025.09.25.678548 (PMC12621732);
> published in *Genome Research* (2026), doi:10.1101/gr.281462.125.

## What this actually tests

Most benchmarks here check that a skill runs and writes plausible artefacts.
This one is stronger: it reproduces a **published numeric artefact**.

The IGVF portal hosts both sides of a pipeline step for AnalysisSet
`IGVFDS1933XFAF` (lentiMPRA, Ahituv lab, processed on the Kircher lab's
cluster — the paper's own MPRAsnakeflow run):

| File | IGVF content type | Role here |
|---|---|---|
| `IGVFFI4856GPLO` | reporter experiment barcode | **input** — per-barcode DNA/RNA counts, 3 replicates, 1,900,181 barcodes |
| `IGVFFI3491STGH` | reporter experiment | **reference** — the authors' aggregated output, 35,110 rows |

`run.sh` feeds the input to `igvfagent mpraflow pipeline` and compares the
result to the reference, value by value.

## Result

**Byte-for-byte reproduction.** All 210,660 compared values
(35,110 rows × 6 columns) match exactly:

| Column | max abs. difference | rows differing |
|---|---|---|
| `dna_counts` | 0 | 0 |
| `rna_counts` | 0 | 0 |
| `n_bc` | 0 | 0 |
| `dna_normalized` | 0 | 0 |
| `rna_normalized` | 0 | 0 |
| `log2FoldChange` | 0 | 0 |

The `(replicate, oligo_name)` key sets are identical — no missing or extra
rows in either direction.

## One parameter had to be recovered

The published file was produced with `minDNACounts=1`, not the tool default
of `0`. This matters because `merge_label.py` adds a DNA pseudocount of 1
per barcode when `minDNACounts == 0`, which inflates `dna_counts` by exactly
`n_bc` per oligo and shifts every `log2FoldChange`.

The benchmark recovered this from the data rather than from documentation:
with defaults, `rna_counts` and `n_bc` matched exactly while `dna_counts` was
high by precisely the barcode count — which identifies the pseudocount as
the only discrepancy. `run.sh` pins `--min-dna-counts 1 --min-rna-counts 1`.

This is worth knowing when processing IGVF MPRA data: the portal's
`reporter experiment` files are not on the tool's default parameters.

## Running it

```bash
bash Benchmarks/rosen2025_mprasnakeflow/run.sh
python3 Benchmarks/concordance.py --paper-id rosen2025_mprasnakeflow
```

Inputs download from the IGVF portal automatically (public, no auth). The
barcode file is ~30 MB compressed; the run takes a few minutes and needs
roughly 2 GB of RAM for the 1.9 M-barcode join.

## Published claims reproduced

Beyond the artefact comparison, the paper's *Results* section states library
complexity figures ("Complexity and sequencing depth estimation"). Each is
recomputed by `igvfagent mpraflow complexity` from the published IGVF
barcode file and checked by `claims.py`:

| Dataset | Metric | Paper | IGVFagent |
|---|---|---|---|
| 8K-neurons (`IGVFFI0032GUOI`) | median assigned barcodes | 1,444,480 | **1,444,480** |
| 8K-neurons | median pairwise Lincoln index | 1,523,572 | **1,523,572** |
| 8K-neurons | library missing | "around 5%" | 5.2% |
| 80K-neurons (`IGVFFI8345QIJJ`) | median assigned barcodes | 5,459,247 | **5,459,247** |
| 80K-neurons | median pairwise Lincoln index | 6,243,618 | **6,243,618** |
| 80K-neurons | library missing | "approximately 12%" | 12.6% |

Both to the digit.

### Which file version the paper used

80K-neurons has two published barcode files. Only
`IGVFFI8345QIJJ` (`results/defaultNGN2`) reproduces the paper; the newer
`IGVFFI0128CWLH` (`results.0.5.4/NGN2`) gives 5,459,214 / 6,243,624 —
off by 33 barcodes in 5.5 M. The benchmark pins the version that matches.

### One claim that does not reproduce exactly

For **240K-HepG2** the paper reports a median of 13,038,694 assigned
barcodes and a Lincoln index of 20,408,916. Neither published file gives
those numbers:

| File | median assigned | median Lincoln |
|---|---|---|
| `IGVFFI1492BOLP` (`results/`) | 13,038,538 | 20,410,542 |
| paper | **13,038,694** | **20,408,916** |
| `IGVFFI0347YDFI` (`results.0.4.6/`) | 13,041,356 | 20,411,802 |

The paper's value falls *between* the two published versions, so the figure
was computed from a third, intermediate pipeline run that is not on the
portal. The gap is 0.001% (156 barcodes in 13 M) and the derived quantity —
"approximately 36% of missing barcodes" — reproduces either way (36.1%).
This is a data-provenance difference, not a disagreement in method: the
same code reproduces 8K-neurons and 80K-neurons to the digit.

## Scope

This validates the **count-aggregation and normalisation** stage and the
**library-complexity estimation**, i.e. the parts of MPRAsnakeflow that
`igvfagent mpraflow` ports.

Not covered: barcode assignment from raw FASTQs (the portal files are
already assigned; the reads are in the ConstructLibrarySets if wanted), the
mapping-strategy comparison (needs BBMap/BWA-MEM), and the outlier-detection
and GC-bias analyses — those live in MPRAlib, which is the next port.
