# Skill: Super-enhancer → target-gene pipeline

End-to-end workflow that takes any ENCODE-supported cell line / tissue (GM12878, K562, HepG2, liver, hippocampus, cortex, …), discovers H3K27ac (or BRD4 / MED1 / P300) ChIP-seq plus optional 3D-chromatin (Hi-C / ChIA-PET / capture Hi-C) experiments, downloads the peak files, calls super-enhancers ROSE-style, and links each super-enhancer to candidate target genes.

Linkage combines four evidence streams per (SE, gene) pair:

1. **3D loops** — Hi-C / ChIA-PET / capture-Hi-C anchors that overlap the SE; the other anchor's region implies a target gene.
2. **IGVF Catalog rE2G / ABC predictions** — region-level regulatory-region → gene maps in the same biosample.
3. **Proximity** — nearest gene TSS within ±N kb of the SE midpoint.
4. **SCREEN cCRE classes** — per-SE composition (PLS / pELS / dELS / CTCF) as quality / specificity evidence (not direct linkage but useful for filtering).

Each (SE, gene) pair is scored by how many of those streams support it. The top pairs are written to a ranked CSV and rendered as a bipartite SE↔gene network.

## Subcommands

### `pipeline` — one-command workflow

```bash
# Default: GM12878 H3K27ac, no 3D
igvfagent se-targets pipeline --biosample GM12878

# Include Hi-C / ChIA-PET 3D chromatin layer + focus on APOE
igvfagent se-targets pipeline \
    --biosample GM12878 --target H3K27ac --assembly GRCh38 \
    --include-3d --gene APOE --label gm12878_apoe_run

# Different cell line, different mark
igvfagent se-targets pipeline --biosample K562 --target BRD4
igvfagent se-targets pipeline --biosample 'liver' --target H3K27ac
```

Outputs land under `Docs/SETargets/<timestamp>_<label>/`:

  - `<label>_super_enhancers.bed`
  - `<label>_se_targets.csv`           # aggregated, sorted by support
  - `<label>_se_targets_loops.csv`     # loop-anchor evidence rows
  - `<label>_se_targets_catalog.csv`   # rE2G / ABC predictions
  - `<label>_se_targets_proximity.csv` # nearest-TSS rows
  - `<label>_se_ccre_composition.csv`  # PLS/pELS/dELS counts per SE
  - `Plots/<label>_se_target_network.png`
  - `<label>_report.md`                # human-readable summary

### `discover` — list candidate experiments

```bash
igvfagent se-targets discover --biosample GM12878 --target H3K27ac
```

Writes a manifest of candidate ChIP-seq + DNase/ATAC + 3D-chromatin experiments for the biosample (no downloads, no SE calling).

### `link-targets` — standalone linkage step

```bash
igvfagent se-targets link-targets \
    --se-bed path/to/super_enhancers.bed \
    --loops-bedpe path/to/hic_loops.bedpe
```

Useful when you've called SEs from a custom upstream pipeline and just want the IGVFagent linkage layer on top.

## How this chains with other skills

- After `pipeline`, drill into any candidate gene with `igvfagent kg gene <SYMBOL> --depth 2 --call-literature`.
- For locus-level visualization, `igvfagent encode browser --region chrN:start-end --track 'SE:<bed>' --with-ccre`.
- For literature corroboration of the high-support gene list, `igvfagent ref validate --input <gene_list.csv>`.
- The internal Plan → Action → Results → Evaluation orchestrator exposes `se_targets_pipeline` as a single tool so a natural-language query like *'For GM12878, call super-enhancers from H3K27ac and tell me which genes they target via Hi-C'* fires the whole workflow in one shot.