# Skill: Enhancer-Gene Linkage Retrieval And Overview

Use this skill when the agent needs enhancer-gene, regulatory element-gene, variant-gene, or enhancer-supporting metadata from IGVF Catalog, IGVF Portal, or ENCODE.

## Evidence Classes

- Experimental linkage: CRISPRi/CRISPRa/Perturb-seq, MPRA/STARR/BlueSTARR, and other tested element-gene or variant-biosample assays.
- eQTL/QTL linkage: variant-gene expression or splicing links from IGVF Catalog `/api/variants/genes` and `/api/variants/genes/summary`.
- Computational linkage: enhancer-gene predictions, variant-element-gene predictions, co-accessibility, distance, or model-based links.
- Context evidence: ENCODE/IGVF accessibility, peaks, chromatin interaction, expression, bigWig signal, and BED interval files.

## Commands

```bash
python3 Scripts/enhancer_gene_linkage_skills.py overview --source catalog --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py overview --source encode --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py overview --source both --limit 10
python3 Scripts/enhancer_gene_linkage_skills.py pull-sets --region chr1:903900-904900 --gene SAMD11
python3 Scripts/enhancer_gene_linkage_skills.py compare-sets --include-local-catalog --demo-if-empty
python3 Scripts/enhancer_gene_linkage_skills.py compare-sets --inputs 'Data/Linkages/*.bed.gz' 'Data/Linkages/*.tsv.gz' --include-local-catalog
python3 Scripts/enhancer_gene_linkage_skills.py write-playbook
```

Use `--source portal` or `--source both` with `IGVF_PORTAL_COOKIE` for unreleased Portal datasets.

## Workflow

1. Start with Catalog linkage endpoints for direct enhancer-gene or variant-gene evidence.
2. Add ENCODE and IGVF Portal metadata for supporting tracks, peaks, expression, and chromatin interaction context.
3. Build a manifest before downloading large files.
4. For a variant list, intersect variants with enhancer/peak intervals, then attach linked genes from Catalog evidence.
5. Prioritize links with convergent support across experimental, eQTL/QTL, and computational methods.
6. Use `compare-sets` to normalize multiple linkage tables into a common element-gene schema, compute pairwise consistency, find links supported by multiple evidence sets, and make bar, heatmap, and arc-view SVG plots.

## Comparison Outputs

- Normalized enhancer-gene row CSV with evidence set, evidence class, method, source, element, gene, context, and score.
- Pair support table showing which element-gene bins are supported by which datasets.
- Evidence set and evidence class bar plots.
- Pairwise Jaccard heatmap for consistency across linkage datasets.
- Arc-style enhancer-gene visualization colored by evidence class.
