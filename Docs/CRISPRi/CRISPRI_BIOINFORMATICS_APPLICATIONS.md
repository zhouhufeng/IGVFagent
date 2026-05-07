# CRISPRi Bioinformatics Research Applications

Generated: 2026-05-07

CRISPRi data are the agent's strongest evidence class for moving from predicted regulatory elements to experimentally perturbed regulatory function. In IGVFdataAgent, CRISPRi should be interpreted together with cCRE, rE2G, single-cell linkage, MPRA, QTL, accessibility, and variant annotations.

## Application Patterns

1. Functional cCRE discovery: identify which SCREEN/ENCODE cCREs have endogenous regulatory effects in a specific cell context.
2. Enhancer-to-gene mapping: connect perturbed regulatory elements to target genes using CRISPRi-FlowFISH, CRISPRi growth screens, TAP-seq, Perturb-seq, or related readouts.
3. Variant-to-gene prioritization: upgrade variants that overlap CRISPRi-supported elements, especially when MPRA, rE2G, QTL, or single-cell linkage evidence points to the same gene.
4. Gene regulatory network inference: use Perturb-seq CRISPRi knockdowns as causal perturbations for downstream transcriptome programs and TF/gene modules.
5. Multiome perturbation interpretation: use paired RNA and ATAC perturbation readouts to connect accessibility changes with expression programs.
6. Screen QC and method comparison: evaluate guide specificity, guide efficiency, effect-size distributions, strand/orientation bias, replicate concordance, and analysis-method sensitivity.

## Implementation In This Repo

```bash
python3 Scripts/crispri_data_skills.py pull --source all --limit 25
python3 Scripts/crispri_data_skills.py analyze-local --input Data/Input/VariantList/example_variants.csv --label my_locus_crispri
python3 Scripts/ccre_linkage_annotation_skills.py linkage-manifest --source all --limit 100 --hydrate-limit 50
python3 Scripts/ccre_linkage_annotation_skills.py browser-demo --region chr19:44850000-44910000
```

## Literature Anchors

- Nature Methods 2024 multicenter noncoding CRISPRi screen analysis: 108 screens, >540,000 perturbations, 24.85 Mb tested sequence, 332 K562 CRE-gene links, and benchmarking of analysis tools.
- Nature Methods 2020 TAP-seq: targeted Perturb-seq for high-throughput enhancer-target maps, including 1,778 enhancers.
- Nature Biotechnology 2024 compressed Perturb-seq: lower-cost single-cell perturbation screens for regulatory circuits and genetic interactions.
- Cell Systems 2025 Multiome Perturb-seq: CRISPRi with paired transcriptome and chromatin accessibility readout.
- Genome Biology 2024 Enhlink: computational enhancer-promoter linkage from scATAC-seq, useful as a comparison/partner evidence class for CRISPRi validation.

## Output Expectations

- CRISPRi metadata manifest in `Data/Manifests/CRISPRi/`.
- Local CRISPRi annotation table and summary stats in `Docs/CRISPRi/`.
- Plots for CRISPRi effect distributions, functional annotation score, CRISPRi-vs-MPRA agreement, and evidence overlap counts.
- IGV-like browser SVGs for prioritized loci through `ccre_linkage_annotation_skills.py browser-demo`.
