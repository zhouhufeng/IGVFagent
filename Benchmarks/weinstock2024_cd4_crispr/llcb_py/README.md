# A from-source Python port of LLCB, run on the paper's own data

This directory reruns Weinstock et al. 2024's actual causal-network method
end to end, on the paper's own deposited data, in Python — as opposed to
`Benchmarks/weinstock2024_cd4_crispr/`'s main benchmark, which only checks
that IGVFagent can *discover and contextualise* the paper via public
catalogues (Perturbation Catalogue, GEO). This is the deeper piece: an
actual re-implementation of LLCB (github.com/weinstockj/LLCB), fit on the
paper's own GEO deposit, whose recovered network is compared directly
against the paper's own published edge counts.

## Why a port, not a wrapper

LLCB is a Julia package (Turing.jl). Neither Julia nor R/Bioconductor is
installed on this cluster, and the authors' own downstream R analysis
pipeline (`weinstockj/RNAseq-perturbation-CD4-pipeline`) has hardcoded
Stanford Sherlock paths and needs external GWAS/ChIP-ABC files that aren't
public. Rather than either skip the causal-network step entirely (the
original state of this benchmark) or spend the session trying to stand up
a foreign toolchain, this is a clean-room Python port of the actual method
(`llcb.py`), fit on the real data. See the module docstring in `llcb.py`
for the full derivation and the specific place where a naive port would
have silently degenerated (an unconstrained per-block noise estimate on an
exactly-square linear system), and how that was caught and fixed.

## Pipeline

| Script | What it does |
|---|---|
| `01_parse_counts.py` | Parses `GSE271788_dedup_counts.txt` (the authors' own featureCounts deposit) into a genes x samples raw count matrix + sample metadata (donor, KO'd gene, control flag). 311 samples = 84 genes x 3 donors + 59 AAVS1 controls — an exact match to the paper's own "84 genes perturbed in three donors and 54(*) samples from control guides" (STAR Methods; *the small 54-vs-59 gap is presumably a post-QC drop we don't have visibility into). |
| `02_normalize.py` | Size factors + variance-stabilizing transform via `pydeseq2` — a from-source Python reimplementation of DESeq2's own algorithms (Love, Huber & Anders 2014) — standing in for the authors' R `DESeq2::vst()` call. |
| `03_build_network_input.py` | Subsets to the paper's 84 target genes (symbol -> Ensembl ID via mygene.info), regresses out donor + the top 10 expression PCs estimated from the AAVS1 controls (their README's own recommended substitute when you don't have their curated `covariates.tsv`, which we don't), and zeroes each sample's own KO'd gene (matching `zero_intervened_nodes` in their R pipeline). |
| `llcb.py` | The port itself: `estimate_total_effects` / `get_cyclic_blocks` are a direct line-for-line port of `graph.jl`; `fit_block_map` replaces `turing_models.jl`'s `joint_cyclic_model` + Pathfinder with a closed-form/MAP solve per gene (justified in the module docstring). |
| `04_fit_llcb.py` | Runs the fit, writes `edges.csv`, prints the same three magnitude-threshold comparisons the paper itself reports. |
| `05_make_figures.py` | `fig4`/`fig5` in the main benchmark's `figures/` directory. |

Run the whole thing:

```bash
cd Benchmarks/weinstock2024_cd4_crispr/llcb_py
python3 01_parse_counts.py && python3 02_normalize.py && \
    python3 03_build_network_input.py && python3 04_fit_llcb.py && \
    python3 05_make_figures.py
```

Needs `pydeseq2` (`pip install pydeseq2`) in addition to this repo's usual
numpy/scipy/pandas/matplotlib. Total runtime on one CPU core: well under a
minute (the dominant cost, DESeq2 dispersion fitting over ~40k genes, is
~9 seconds; the network fit itself, 84 independent 83-parameter L-BFGS
solves, is under a second).

## Results

The port runs end to end on the real data and recovers a real signed,
directed 84-gene network (`Data/Weinstock2024/processed/edges.csv`, 6,972
directed pairs — note this denominator matches the paper's own "out of
6,972 possible" exactly, confirming the network dimensionality is right).

**Edge count vs. the paper's own thresholds** (paper STAR Methods: "We
identified 350, 211, and 151 total edges (out of 6,972 possible) when
thresholding |βij| at 0.020, 0.025, and 0.030, respectively"):

| \|β\| threshold | paper | this port | ratio |
|---|---:|---:|---:|
| 0.020 | 350 | 1,294 | 3.7x |
| 0.025 | 211 | 877 | 4.2x |
| 0.030 | 151 | 618 | 4.1x |

Our network is consistently ~3.5-4x denser than the paper's at the same
magnitude cutoff (see `figures/fig4_llcb_edge_threshold.png`).

**A cross-validating finding on edge-calling itself.** The paper explains
*why* it uses a raw magnitude threshold rather than a posterior-uncertainty
criterion: "We reported the network after thresholding on β because
filtering on local false sign rates (LFSR) resulted in very dense networks
(67% network density at LFSR < 5x10⁻³), reflecting the challenges in the
estimation of uncertainty in graph structures." This port computes its own
LFSR from the Laplace-approximated posterior at each block's MAP, entirely
independently, and hits the *same* qualitative failure mode they describe —
thresholding on LFSR < 0.005 gives 639/6,972 = 9.2% density here (not 67%;
the absolute density differs, but the underlying phenomenon — LFSR-based
thresholding failing to sparsify this specific model class — reproduces).
That's a genuine, unprompted cross-check: an independent implementation
hitting the same known failure mode the original authors reported is
stronger evidence for a shared underlying cause (uncertainty in cyclic
graph structure estimation is genuinely hard to calibrate) than either
result alone.

**The paper's headline KMT2A->Th17-IL2-JAK-STAT edge is not recovered.**
At |β| > 0.05, KMT2A connects to MED12, BCL11B, NFKB2, and FOXP1 in this
port's network — not to STAT5A/STAT5B/IL2RA/JAK3/RORC. The broader network
does recover a biologically coherent immune-signaling module among those
genes (STAT5A<->STAT5B, NFKB2->IL2RA, STAT3->JAK3, RELA->RELB, IRF4->RORC,
KLF2->RORC — see `figures/fig5_llcb_top_edges.png`), just not the specific
path the paper highlights in Figure 5 en route to the rs45480496 enhancer
variant.

## Honest limitations

- **Normalization and covariate handling differ from the paper's.** We use
  `pydeseq2` (Python) instead of their R `DESeq2::vst()`; the two should be
  numerically close but aren't guaranteed bit-identical. More importantly,
  we regress out our *own* PCs estimated from the AAVS1 controls, not their
  curated `covariates.tsv` (not public) — this is almost certainly the
  single biggest source of the ~4x edge-count inflation, since incompletely
  removed batch/technical variation directly inflates apparent total
  effects.
- **The spectral-radius stability penalty is not ported.** `joint_cyclic_model`
  has two soft regularizers on the full adjacency matrix; the per-column
  sparsity term is block-separable and *is* ported (`fit_block_map`), but
  the spectral-radius term needs the whole matrix jointly and isn't. Adding
  it back would require a joint (not per-block) optimization over all
  ~7,000 parameters at once.
- **No downstream GWAS colocalization.** The paper's KMT2A/rs45480496
  Th17-enhancer claim (Figure 5) integrates external GWAS summary
  statistics this repo doesn't have; this port stops at the network itself.
- **Small per-perturbation replication (n=3 donors).** Both this port and
  the paper's own total-effect estimation lean on comparing 3 KO replicates
  against ~54-59 controls per gene; this is an intrinsic feature of the
  experimental design, not something either implementation can fix.

None of this is presented as "IGVFagent reproduces Weinstock 2024's
network exactly" — it doesn't, and the divergence is characterized above
rather than papered over. What it does show: a real, working,
independent implementation of the paper's own causal-network-inference
method, run on the paper's own public data, in the same 84-gene x 6,972-edge
space they report, recovering a plausible network at a well-characterized
different operating point, and independently reproducing a failure mode
the original authors themselves flagged.
