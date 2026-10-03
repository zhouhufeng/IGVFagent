# Skill: Reference (literature retrieval, validation, study design)

Use this skill when an IGVF agent run needs to consult prior literature: to scope a problem, to validate discoveries, or to design a study workflow that mirrors successful published practice.

## Sources

- PubMed / PMC (NCBI E-utilities)
- bioRxiv + medRxiv (via Crossref)
- arXiv (Atom API)
- Semantic Scholar (paper search + recommendation)
- OpenAlex (open scholarly graph)

Top-tier journals are weighted in ranking: Nature/Cell/Science families, NEJM, Nucleic Acids Research, Bioinformatics, Genome Biology / Genome Research, eLife, PNAS.

## Subcommands

### 1. `learn` — what does the field do?

```bash
python3 Scripts/reference_skill.py learn \
    --topic '10x multiome human putamen' --limit 30 --top 15
```

Returns a manifest of top-ranked papers and a report extracting:
- recurring methodology phrasing,
- visualization vocabulary used in abstracts (UMAP / dot-plot / volcano / track / heatmap …),
- a consensus figure recipe.

### 2. `validate` — has anyone seen this before?

```bash
python3 Scripts/reference_skill.py validate \
    --input Docs/AdvancedVariantAnalysis/<run>/<label>_summary_stats.csv \
    --context 'putamen Parkinson' --limit-per-item 5
```

Input is any CSV with `gene` / `variant` / `rsid` / `element` / `term` columns (case-insensitive). Each row is searched against PubMed and Semantic Scholar; the manifest links each discovery to the prior-evidence papers, and the report ranks the strength of prior support.

### 3. `design` — what should our study look like?

```bash
python3 Scripts/reference_skill.py design \
    --data-type parse_split_seq --assay-title 'Parse SPLiT-seq'
```

Known data types: `10x_multiome`, `parse_split_seq`, `mpra`, `crispri`, `enhancer_gene`.

Outputs: a recommended pipeline (phases + plots/tables), a literature manifest of cognate published studies, and a manifest of IGVF Portal AnalysisSets that match the assay.

### `search` — generic multi-source search

```bash
python3 Scripts/reference_skill.py search \
    --query 'enhancer-gene linkage rE2G ABC' --top 20
```

## Caching

All API responses are cached under `Data/Cache/References/<source>/` with a 14-day TTL. Re-running the same query is free.

## Outputs

- Reports: `Docs/References/<timestamp>_<subcommand>_<label>/*.md`
- Manifests: `Data/Manifests/References/`
- Logs: `Docs/Logs/reference_skill_*.log`

## How this skill chains with other IGVF agent skills

1. After `data_illustration_interpretation.py` produces a dataset summary, call `learn` with the assay/biosample to get the typical analysis recipe.
2. After `advanced_variant_analysis.py` writes a discovery table, feed the `summary_stats.csv` (or any gene/variant CSV) into `validate`.
3. Before starting a new study, call `design` with the planned IGVF data type to seed the workflow and surface matching IGVF datasets.