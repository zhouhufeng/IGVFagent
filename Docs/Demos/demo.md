# IGVFagent — demo questions

Question text is verbatim: copy-paste it straight into the chat box at
<https://igvfagent.genohub.org>.

Every dataset below was looked up on the Portal rather than assumed, so the
assay, route and **download size** are measured. Size is the thing that
decides what is safe to run live: the agent fetches raw reads, and a 72 GB
dataset does not finish inside a talk.

## Read this before presenting

**Run the big ones ahead of time.** Counts and matrices are cached, so a
dataset analysed once returns in seconds afterwards. The LDLR screen went
from 15 minutes to 9.6 seconds on its second run. Anything over ~2 GB in the
tables below should be pre-warmed the day before, not discovered on stage.

**One bin of a screen is not a dataset.** For the CRISPR screens, the
measurement IS the comparison between sorted bins, so the agent pulls the
screen's siblings and analyses them together. Asking about one accession
alone yields library composition and no biology, and the tools say so rather
than quantifying it anyway.

**Heavy jobs are not yet rate-limited.** `raw_pipeline_run`,
`crispr_screen_analyze`, `sge_analyze`, `mct_analyze` and
`gradient_screen_analyze` can all run concurrently, and several large ones at
once can exhaust the container's 22 GB and drop every session. With an
audience able to type their own queries, pre-warm the demo datasets and keep
the unbounded ones off the slide.

---

## Act 1 · Single-cell and multiome (the transcript route)

These go FASTQ -> count matrix (kallisto|bustools) -> QC -> clustering ->
plots, or reuse a published matrix from a derived AnalysisSet when one
exists, which is faster and is what the pipeline prefers.

| # | Accession | Assay | Size | Live-safe? |
|---|---|---|---|---|
| 1 | IGVFDS7013XXYV | 10x multiome | 1.35 GB | yes |
| 2 | IGVFDS5498IKCV | Perturb-seq (CRISPRi, MHC) | 21.9 GB | pre-warm |
| 3 | IGVFDS5414UFNC | 10x multiome | 26.9 GB | pre-warm |
| 4 | IGVFDS2149HIIR | Parse SPLiT-seq | 34.4 GB | pre-warm |
| 5 | IGVFDS6179NVYH | 10x multiome | 41.0 GB | pre-warm |
| 6 | IGVFDS3532MONX | TAP-seq | 45.6 GB | pre-warm |
| 7 | IGVFDS9875NBZW | 10x multiome with MULTI-seq | 72.4 GB | pre-warm |
| 8 | IGVFDS2597SSBK | MORF-SHARE-seq | 16.7 GB | pre-warm |
| 9 | IGVFDS5674PNNL | SHARE-seq | 30.4 GB | pre-warm |

```
can you explain and analyze the raw data of IGVFDS7013XXYV and visualize the analysis results?
```
- **Start here.** Smallest transcript dataset, so it is the one that finishes
  while you talk.
- **Watch for:** the route decision and download size printed before any
  work begins; then QC violins, a PCA scree plot, UMAP and marker heatmap
  rendered inline.

```
can you perform QC and then in-depth data analysis of the raw data in IGVFDS5498IKCV and then visualize the QC and analysis results?
```
- **Why this one is interesting:** it is a targeted-panel Perturb-seq, where
  cells legitimately express few genes (median 68 genes, ~2,500 UMIs per
  cell). A default `min_genes=200` collapses 340,599 barcodes to 188 cells.
  Good moment to show the agent reporting a QC threshold as the cause rather
  than presenting 188 cells as the answer.
- **Route:** `matrix_derived` — it reuses the published matrix from
  AnalysisSet IGVFDS3959LESA instead of re-aligning 21.9 GB.

```
can you analyze this dataset? IGVFDS9875NBZW
can you help me understand analyze and process this data? IGVFDS3532MONX
can you analyze the raw data of IGVFDS9875NBZW and visualize the analysis results?
can you explain and analyze the raw data of IGVFDS9875NBZW and visualize the analysis results?
can you analyze IGVFDS6179NVYH and visualize the analysis results?
can you analyze IGVFDS2149HIIR and visualize the analysis results?
can you perform QC and then in-depth data analysis of the raw data in IGVFDS5414UFNC and then visualize the QC and analysis results?
can you perform QC and then in-depth data analysis of the raw data in IGVFDS2597SSBK and then visualize the QC and analysis results?
can you perform QC and then in-depth data analysis of the raw data in IGVFDS5674PNNL and then visualize the QC and analysis results?
```
- The phrasings differ deliberately: "analyze this dataset", "explain and
  analyze", "perform QC and then in-depth analysis". All should reach the
  same route. Useful for showing the request does not have to be phrased in
  the tool's language.

---

## Act 2 · CRISPR screens (two different shapes)

The important point for an IGVF audience: these look similar and are not.
A **tail sort** compares bottom20% against top20%. A **lettered-bin
gradient** has bins A-F and no tails, so the statistic is each construct's
frequency-weighted mean bin. 554 Portal MeasurementSets have the gradient
shape.

| # | Accession | Shape | Tool | Size |
|---|---|---|---|---|
| 10 | IGVFDS6464SOVZ | tail sort, 4 reps x 4 bins | `crispr-screen` | 0.02 GB |
| 11 | IGVFDS5997IVEM | tail sort, prime editing | `crispr-screen` | 0.32 GB |
| 12 | IGVFDS8710ZSOZ | gradient, bins A-F | `gradient-screen` | 0.06 GB |
| 13 | IGVFDS3899ANMJ | gradient, bins A-F | `gradient-screen` | 0.06 GB |

All four are small. **These are the best live demos in the set.**

```
can you analyze IGVFDS6464SOVZ and visualize the analysis results?
can you analyze IGVFDS6464SOVZ and visualize the analysis results? it needs to pull related data and analyze together
```
- **Watch for:** one accession in, 16 libraries out — the agent finds the
  screen's other bins and replicates (`18loci_uptake`, 4 replicates x 4
  bins) and scores 1,656 targets. The second phrasing makes the sibling
  pulling explicit; the first should do it anyway.

```
can you perform QC and then in-depth data analysis of the raw data in IGVFDS5997IVEM and then visualize the QC and analysis results?
```
- **The set-piece.** This is Rep2 top20% of the Sherwood `LDLR137-219`
  screen, and its library is **prime editing**: 1,741 pegRNAs share just 52
  spacers, because the spacer only sets the nick site and the variant lives
  in the RT template. Counting it by spacer would return effect sizes and
  FDRs for variants whose reads were never told apart.
- **Watch for:** the key-calibration table. The tool tests each candidate
  column against real reads and picks `rt_template_sequence` (62% of reads,
  0% ambiguous) over `spacer` (3% of constructs separated) and
  `peg_sequence` (perfectly unique, 0% of reads — its entries are longer
  than the 128 bp read). Then: 20 libraries, 1,732 variants scored, and
  **zero** significant after correcting for 1,732 tests. That last part is
  the honest answer for this screen, not a failure.

```
can you perform QC and then in-depth data analysis of the raw data in IGVFDS8710ZSOZ and then visualize the QC and analysis results?
```
- **The other set-piece.** KITLG CRISPRi FlowFISH, 24 sets over 4 flow
  replicates, bins A-F. Counts by `spacer` at 81% of reads.
- **Watch for:** 65 hits, all negative — perturbing KITLG regulatory
  elements lowers KITLG. And the caveat the run prints: its four
  "replicates" are flow re-sorts of one biological sample, so their spread
  is instrument precision. The FDR is calibrated against the scatter of the
  library's 434 non-targeting controls instead, giving a 0.46% control
  false-positive rate. Using replicate spread called 27% of those controls
  hits.

```
can you perform QC and then in-depth data analysis of the raw data in IGVFDS3899ANMJ and then visualize the QC and analysis results?
```
- PPIF promoter prime editing, 66 sets over 11 sorts, read out by
  sequencing the **endogenous** locus — the RT template is written into the
  genome, so it is found in the amplicon (65% of R2 reads; 0% of R1, because
  the edit sits at one end of an amplicon longer than 2x150 bp).
- Unlike KITLG this screen has real biological replicates, so replicate
  spread is a legitimate null here.

---

## Act 3 · Variant-effect assays

| # | Accession | Assay | Tool | Size |
|---|---|---|---|---|
| 14 | IGVFDS4629JYPY | SGE | `sge` | 0.14 GB |
| 15 | IGVFDS0865HLQG | VAMP-seq (MultiSTEP) | `sge` / MAVE | 1.64 GB |
| 16 | IGVFDS4826YNLK | snMCT-seq | `mct` | 50.6 GB |

```
can you analyze this dataset IGVFDS4629JYPY and then visualize the analysis results
can you perform QC and then in-depth data analysis of the raw data in IGVFDS0865HLQG and then visualize the QC and analysis results?
```
- Saturation genome editing and a multiplexed assay of variant effect. Both
  small enough to run live. `raw_pipeline` refuses these for a stated
  reason: the reads are a fixed amplicon, so quantifying them against a
  transcriptome would restate the amplicon design.

```
can you perform QC and then in-depth data analysis of the raw data in IGVFDS4826YNLK and then visualize the QC and analysis results?
```
- **Two measurements from the same nuclei**: RNA and DNA methylation.
  Pre-warm it (50.6 GB).
- **Watch for:** both halves analysed and cross-compared. Routing this as a
  transcript assay analyses the RNA and silently drops the methylation --
  the half the assay exists for. 933 Portal datasets are snMCT-seq. Also
  worth saying out loud: methylation ratios are NOT log-normalised, which is
  the mistake that makes a methylome look like an expression matrix.

---

## Act 4 · Knowledge Graph and mechanism

No raw data, so these are fast and safe to take from the audience.

```
whether there are enhancer-gene predictions that overlap the variant rs1250566 and show me the visualization also
```
- **Watch for:** the distinction between "no evidence" and "no such edge in
  the graph". There is no direct variant->element edge for this variant, so
  the honest route is to query by genomic region. A tool that returns zero
  rows and reports "no predictions" is indistinguishable from absence of
  evidence, and that is a wrong answer, not a null one.

```
Which enhancers regulate PCSK9 in hepatocytes, which variants alter their activity, and what IGVF evidence confirms the mechanism?
```
- A full mechanism chain: element -> gene -> variant -> functional evidence.
  Good closing question because the answer is a narrative with citations to
  specific IGVF datasets.

```
BRCA1 in the Perturbation Catalogue - MAVE, CRISPR screens, Perturb-seq, and top GSEA hallmarks.
```

```
Everything in the IGVF KG about APOE - variants, transcripts, proteins, regulatory elements, diseases, and pathways.
```
- One tool call, multi-hop evidence pack. Open the generated `report.md`
  from the artefact pane.

---

## Act 5 · What the agent declines, and why

Worth a slot of its own in front of a methods audience: the useful answer is
sometimes "not this way, and here is the reason".

| # | Accession | Assay | Route | Size |
|---|---|---|---|---|
| 17 | IGVFDS9961QTUB | ATAC-seq | `chromatin` | 6.35 GB |

```
can you perform QC and then in-depth data analysis of the raw data in IGVFDS9961QTUB and then visualize the QC and analysis results?
```
- **Expect no plots, and that is the correct outcome.** ATAC-seq reads are
  genomic, and the deployment ships only a transcriptome aligner
  (kallisto|bustools). The agent reports:

  > ROUTE: assay_mismatch — this is ATAC-seq (route: chromatin), whose reads
  > are not a transcript library. Quantifying them against a transcriptome
  > would report which gene the amplicon covers -- the assay design restated,
  > not a result. The readout you want is peak or contact analysis against
  > the genome -- NOT supported end-to-end; no aligner for genomic reads.

- **The point to make:** only 59.5% of the Portal's 11,070 MeasurementSets
  are transcript assays. A pipeline that defaults to "quantify it as RNA"
  produces a plausible-looking answer for the other 4,486, and a UMAP built
  that way is indistinguishable from a real one. Declining is the feature.
- `force_align=true` overrides it if someone insists, which is also worth
  showing: the guard is a default, not a wall.

---

## If something goes wrong on stage

- **"Agent run ended before completion"** — the model did not use the tool
  protocol. It retries once automatically; if you see this, re-run the query
  or pick a different model in the sidebar. It is not a data problem.
- **A run reports "Done" having made no tool calls** — same cause, and the
  answer text may contain an invented explanation such as tools being
  unavailable. Treat the tool-call count in the footer as the ground truth.
- **Plots described in prose but not shown** — the answer names a directory
  the renderer could not expand. The artefact pane at the bottom is the
  fallback; every figure is also written to disk under `Docs/`.
- **A long job with no visible progress** — check the sidebar for the
  running-jobs panel and the heavy-slot line.
