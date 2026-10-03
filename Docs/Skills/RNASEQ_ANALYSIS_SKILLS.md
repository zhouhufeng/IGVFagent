# Skill: bulk RNA-seq analysis + DEG → controlling cCRE

Picks up where the GEO retrieval skill leaves off: takes a counts matrix + sample sheet, runs QC, PCA, differential expression, and (optionally) links each significant DEG back to the regulatory elements that control it via the IGVF Catalog.

DEG engine: pyDESeq2 if installed (preferred); otherwise Welch's t-test on log-CPM with Benjamini-Hochberg FDR. The fallback is fast and dependency-free; for publication-quality calls, install pyDESeq2 or use proper tools (DESeq2/edgeR/limma in R).

## Subcommands

### `pipeline` — one-command end-to-end

```bash
igvfagent rnaseq pipeline \
    --counts counts.tsv --sample-sheet samples.csv \
    --condition-col condition --group-a control --group-b treated \
    --label demo_run
```

Writes QC CSV + PCA + volcano + MA + heatmap + DEG CSV + (with Catalog reachable) DEG → controlling cCRE table + a summary markdown report under `Docs/RNAseq/<ts>_<label>/`.

### `qc`, `pca`, `deg`, `link-cre` — individual steps

```bash
igvfagent rnaseq qc --counts counts.tsv
igvfagent rnaseq pca --counts counts.tsv --sample-sheet samples.csv \
    --condition-col condition
igvfagent rnaseq deg --counts counts.tsv --sample-sheet samples.csv \
    --condition-col condition --group-a control --group-b treated
igvfagent rnaseq link-cre --deg <deg.csv> --padj-cut 0.05 --fc-cut 1.0
```

## How this chains with other skills

1. `igvfagent geo series --gse GSE9574 --full-samples` →
2. `igvfagent geo download --gse GSE9574 --only matrix` →
3. `igvfagent geo sample-sheet --gse GSE9574` →
4. `igvfagent rnaseq pipeline --counts <matrix> --sample-sheet <sheet>` →
5. `igvfagent kg gene <SYMBOL> --depth 2 --call-literature` for any DEG you'd like to drill into.

## Outputs

- Per-sample QC: `Docs/RNAseq/<ts>_<label>/<label>_qc.csv`
- DEG table:     `<label>_deg.csv`
- Plots:         volcano / MA / heatmap / PCA (.png + .svg)
- DEG → cCRE:    `<label>_deg_to_cre.csv` (when Catalog reachable)
- Report:        `<label>_report.md`