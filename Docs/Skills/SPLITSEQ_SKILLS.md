# Skill: Parse SPLiT-seq pipeline (IGVF)

Reads, processes, and analyzes Parse-Biosciences-style combinatorial-barcoding single-nucleus RNA-seq datasets from the IGVF Portal. Designed around SPLiT-seq's structural quirks: multiplexed sub-pools, tens of donors per pool, mtx tarballs, and the analytical value living in **per-strain / per-donor / per-pool** comparisons.

## Subcommands

### `retrieve` — discover SPLiT-seq AnalysisSets

```bash
python3 Scripts/splitseq_pipeline.py retrieve \
    --limit 100 --fetch-file-details --label mortazavi_8cube
python3 Scripts/splitseq_pipeline.py retrieve \
    --limit 50 --lab 'Ali Mortazavi, UCI' --taxa 'Mus musculus'
python3 Scripts/splitseq_pipeline.py retrieve \
    --limit 20 --sample-type gonad
```

### `manifest` — per-pool / per-donor manifest for one or more accessions

```bash
python3 Scripts/splitseq_pipeline.py manifest \
    --accessions IGVFDS3222WCZH,IGVFDS6290WNNH \
    --label 8cube_mouse_demo
```

### `download` — pull files for a manifest under a size cap

```bash
python3 Scripts/splitseq_pipeline.py download \
    --manifest Data/Manifests/SPLiTseq/<files_per_pool>.csv \
    --max-download-gb 5 \
    --only 'sparse gene count matrix' 'cell annotations'
```

### `process-local` — load downloaded mtx bundles into AnnData

```bash
python3 Scripts/splitseq_pipeline.py process-local \
    --label 8cube_mouse_demo
```

Concatenates every downloaded AnalysisSet under `Data/IGVF/SPLiTseq/Downloads/` into a single `.h5ad` under `Data/Cache/SPLiTseq/`. Each cell carries a `dataset` and `sample_id` column for downstream batch correction.

### `analyze` — full pipeline (QC → norm → integrate → UMAP → cluster → label)

```bash
python3 Scripts/splitseq_pipeline.py analyze \
    --input Data/Cache/SPLiTseq/8cube_mouse_demo.h5ad \
    --tissue gonad \
    --resolution 0.6 \
    --label 8cube_mouse_demo
```

Tunable QC thresholds via `--min-genes`, `--max-genes`, `--min-counts`, `--max-counts`, `--max-pct-mito`. Pool integration uses Harmony if `scanpy.external` is available, else falls back to uncorrected PCA.

Auto-annotation panels currently shipped: `gonad`, `adrenal`, `brain`, `liver`, `heart`, `kidney`, `muscle`. Pass `--tissue` to use one. If the AnnData already carries a `cell_name` obs column from an IGVF cell-annotation TSV (as in some principal AnalysisSets), that label is used directly.

### `plot` — publication-style multi-panel figure

```bash
python3 Scripts/splitseq_pipeline.py plot \
    --input Data/Cache/SPLiTseq/8cube_mouse_demo_processed.h5ad \
    --tissue gonad --label 8cube_mouse_demo
```

Panels: UMAP by cell type, UMAP by pool/batch, UMAP by strain (if the `strain` obs column is populated), marker dot plot, stacked composition bar.

### `compare-strains` — per-cell-type DEG across the 8 founders

```bash
python3 Scripts/splitseq_pipeline.py compare-strains \
    --input Data/Cache/SPLiTseq/8cube_mouse_demo_processed.h5ad \
    --label 8cube_strain_DEG --padj-cutoff 0.05 --min-cells 50
```

Requires a `strain` (or `strain_background` / `DonorStrain`) column on `obs`. Populate it after running donor demultiplexing.

### `demux-script` — emit a runnable demux pipeline

```bash
python3 Scripts/splitseq_pipeline.py demux-script --tool souporcell --n-donors 8
python3 Scripts/splitseq_pipeline.py demux-script --tool vireo --n-donors 8
```

The skill does NOT run the demultiplexer — these tools need BAMs and (for vireo) per-donor VCFs that live outside the agent. The emitted script wires up the common-case command lines so you can run them locally and feed the resulting `donor_ids.tsv` back as a `strain` obs column.

## Outputs

- `Data/IGVF/SPLiTseq/Metadata/<accession>.json` — hydrated portal JSON
- `Data/IGVF/SPLiTseq/Downloads/<accession>/` — per-AnalysisSet bundles
- `Data/Manifests/SPLiTseq/` — dataset + per-pool manifests, download log
- `Data/Cache/SPLiTseq/<label>.h5ad` — processed AnnData
- `Docs/SPLiTseq/<timestamp>_<label>/` — reports + plots

## How this chains with other skills

- After `retrieve`, hand the manifest to `download` and then `process-local`.
- After `analyze` produces a cell-type-labeled `.h5ad`, hand a gene list (e.g. top markers) to `Scripts/reference_skill.py validate` to check prior literature.
- For workflow scaffolding before retrieval, run `Scripts/reference_skill.py design --data-type parse_split_seq` to see the curated workflow plus cognate published studies.
- The bundled mouse marker panels (gonad / adrenal / brain / liver / heart / kidney / muscle) match the tissues covered by the Mortazavi 8-cube founder atlas, so a one-shot run through `retrieve → download → process-local → analyze → plot` produces publication-style figures with no manual annotation.