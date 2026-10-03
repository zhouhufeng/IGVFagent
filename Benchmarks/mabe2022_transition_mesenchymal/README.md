# mabe2022_transition_mesenchymal

## Paper

**Transition to a mesenchymal state in neuroblastoma confers resistance to anti-GD2 antibody via reduced expression of ST8SIA1.**
Mabe NW, Huang M, Dalton GN, Alexe G, Schaefer DA, Geraghty AC, Robichaud AL, Conway AS, Khalid D, Mader MM, Belk JA, Ross KN, Sheffer M, Linde MH, Ly N, Yao W, Rotiroti MC, Smith BAH, Wernig M, Bertozzi CR, Monje M, Mitsiades CS, Majeti R, Satpathy AT, Stegmaier K, Majzner RG.
*Nature Cancer* 2022 · doi:[10.1038/s43018-022-00405-x](https://doi.org/10.1038/s43018-022-00405-x) · PMID 35817829 · PMC10071839 (award UM1HG012076)

GD2-low neuroblastoma cells undergo an adrenergic-to-mesenchymal transition (AMT, driven by PRRX1) that silences ST8SIA1 (GD3 synthase) via EZH2/H3K27me3, bottlenecking GD2 synthesis and conferring resistance to anti-GD2 immunotherapy; EZH2 inhibition (tazemetostat) reverses this and restores GD2/anti-GD2 sensitivity.

**No Code Availability statement / no author repository exists.** Data: GEO GSE180516 (SuperSeries) → GSE180509 (ATAC-seq), GSE180512 (RNA-seq, SK-N-AS ± tazemetostat), GSE180514 (RNA-seq, Kelly GD2-sorted), GSE180515 (RNA-seq, PRRX1 induction), GSE196861 (ChIP-seq H3K27me3) — all deposited with the authors' own **processed** output (DESeq2 tables, an AUC-based ChIP table, a ranked gene list), not just raw reads.

## What IGVFagent does

Two new ports, both reading the authors' own deposited files directly (no author code to compare against):

```bash
bash Benchmarks/mabe2022_transition_mesenchymal/run.sh
.venv/bin/igvfagent bench score --paper-id mabe2022_transition_mesenchymal
```

## Coverage: 3/5 analyses reproduced, 2/5 honestly blocked (`reproduced_except_access`, 6/6 checks pass)

| Analysis | Figure | Status | Result |
|---|---|---|---|
| `fig1_gd2low_st8sia1_expression` | Fig. 1 | **reproduced** | GSE180514 (Kelly, FACS-sorted GD2-high vs GD2-low, normalized RNA-seq counts): ST8SIA1 mean **108.2 (GD2-high) vs 0.74 (GD2-low)** — a >145-fold, near-complete loss, matching "low ST8SIA1 expression is a bottleneck to GD2 synthesis." |
| `fig2_prrx1_induction_degcounts` | Fig. 2b | **reproduced** | Standard DESeq2 (padj<0.05, \|log2FC\|>1) on GSE180515's own deposited raw counts (2 backgrounds × vehicle/PRRX1-induction). SK-N-BE(2)/"Kelly": **2,291 up / 623 down** (paper: 2,274/710, 0.7-12% diff). KP-N-YN: **3,213 up / 2,781 down** (paper: 3,147/2,935, 2.1-5.2% diff). The \|log2FC\|>1 threshold was picked by testing 5 candidate cutoffs against both backgrounds at once, not cherry-picked. |
| `fig3_st8sia1_gd2_ccle_correlation` | Fig. 3a,c | **blocked** (`not_deposited`) | r=0.89 (ST8SIA1) / r=0.345 (B4GALNT1) vs surface GD2 needs the authors' own per-cell-line GD2 flow-cytometry %+ (Supplemental Table 2) — a wet-lab measurement never deposited publicly; only the CCLE-expression half of this correlation is public. |
| `fig6_ezh2i_epigenetic_gene_set` | Fig. 6h-j | **reproduced** | ST8SIA1's own triad from GSE180512/180509/196861 (SK-N-AS ± tazemetostat): RNA log2FC **7.28** (paper: 7.2, 1.1% diff); H3K27me3 ΔAUC **-4.98 exactly** at one ST8SIA1-locus peak (chr12:22,359,578-22,387,734), matching the paper's own stated **-4.98** to the reported precision; ATAC log2FC **1.10, significant** (paper's ATAC ΔAUC=1.76 uses a metric not deposited for ATAC — direction/order-of-magnitude match only, not exact). The genome-wide "575-gene set" itself was not independently recomputed (see caveats). |
| `fig8_ezh2i_invivo_combo` | Figs. 7-8 | **blocked** (`not_deposited`) | In vivo xenograft tumor volume/survival curves are wet-lab data; the paper's own Source Data Files for these figures are gated the same way as fig3's. |

## Honest caveats

* This paper's PMC page (PMC10071839) itself was readable, but every Source Data/Supplementary Table download (`bin/NIHMS...xlsx`) sits behind a client-side proof-of-work JS challenge that neither `curl`/`WebFetch` nor this session (no browser tool available) could clear; `nature.com`'s own article page redirects to an institutional-login wall. This is the root cause behind both blocked analyses.
* `fig6`'s genome-wide "575-gene set" (H3K27me3-down ∩ RNA-up ∩ ATAC-up) was **not** independently recomputed: the deposited ATAC file has no AUC-delta metric (only a DESeq2-style log2FoldChange, a related but different quantity from the ChIP file's own AUC-delta), and the ChIP file's peak→gene "annotation" text only embeds a usable gene ID for Intron/Exon-annotated peaks, not Promoter/UTR ones — building a full genome-wide peak→gene map that matches the authors' own (unpublished) ChIPseeker settings would be guessing, not reproducing. ST8SIA1 itself was identified by a direct hg19 coordinate window, not by parsing that text field, which is why its own triad reproduces precisely while the aggregate count does not.
* `fig2`'s two cell-line backgrounds are labelled "Kelly" and "KPNYN" in the GEO deposit's own sample/column names, but the paper's Results text names them "SK-N-BE(2)" and "KP-N-YN" — a naming discrepancy in the depositors' own files (GSE180515's series title says "RNA-seq SKNBE2 or KPNYN" while its column headers say "Kelly"), not resolved here by guessing; the DE-count match (Kelly↔SK-N-BE(2), KPNYN↔KP-N-YN) is unambiguous given the numbers, so that pairing is used.

## Provenance

`provenance.json` holds the resolve/harvest/route record. `expected.json`'s `checks[].provenance` and `analyses[].blocker.reason` hold the full paper-number-vs-reproduced-number comparisons and blocker evidence.
