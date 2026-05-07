# Skill: MPRA Data Retrieval And Analysis

Use this skill when the agent needs to retrieve MPRA/STARR/BlueSTARR evidence from IGVF Catalog, IGVF Portal, or ENCODE, or analyze local MPRA result tables.

## Retrieval

```bash
python3 Scripts/mpra_data_skills.py pull --source catalog --limit 25
python3 Scripts/mpra_data_skills.py pull --source all --limit 25
python3 Scripts/mpra_data_skills.py portal-manifest --limit 100 --label igvf_portal_mpra_many
```

Use `IGVF_PORTAL_COOKIE` for unreleased Portal datasets.

## Local Analysis

```bash
python3 Scripts/mpra_data_skills.py analyze-local --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra
python3 Scripts/mpra_data_skills.py literature-demo --input Data/Input/VariantList/example_variants.csv --label my_locus_mpra_literature_demo
```

The analysis writes summary stats, SVG plots, and a Markdown report.

## Literature-Informed Templates

- Cell / PubMed variant-effect MPRA: use volcano plots, significant allelic-effect tables, and integration with GWAS/eQTL/fine-mapping.
- Nature Methods MPRA design benchmarking: use count QC, assay-design comparisons, activity dynamic range, and sequence-context summaries.
- Nature Genetics single-cell MPRA: use cell-type activity heatmaps and cluster/cell-type specificity summaries.
- Nature Biotechnology regulatory grammar MPRA: use saturation-mutagenesis maps, position-effect heatmaps, and model-vs-observed plots.
- Nature large-scale cCRE MPRA: use activity distributions, cCRE class enrichment, and genome-browser style views.

## Reuse Rules

- Inspect metadata before downloading large MPRA files.
- Preserve variant identifiers, allele orientation, library/source, biosample, activity class, input/output counts, log2 fold-change, P/Q values, and significance calls.
- Check input balance and low-count variants before interpreting effect sizes.
- Compare MPRA effects with IGVF Catalog variant-gene, variant-biosample, regulatory-element, and enhancer-gene prediction evidence.
- Use `portal-manifest` to pull many IGVF Portal MPRA/STARR/reporter datasets before choosing files to download.
- Use `literature-demo` to create paper-style interpretation plots and a research-use report from a local MPRA table.
