# Skill: CRISPRi Data Access, Download Planning, And Functional Annotation Integration

Use this skill when the agent needs CRISPRi, CRISPR FACS screen, or Perturb-seq regulatory perturbation evidence from IGVF Portal, IGVF Catalog, or ENCODE.

## Bioinformatics Research Use Cases

CRISPRi data are most useful when the question needs endogenous perturbation evidence rather than correlation-only regulatory annotation.

1. **Functional cCRE discovery and benchmarking.** Noncoding CRISPRi screens can directly test candidate cis-regulatory elements and compare hits against cCRE, ATAC/DNase, H3K27ac, TF-binding, MPRA, and rE2G annotations. The ENCODE multicenter analysis assembled 108 noncoding CRISPR screens, >540,000 perturbations, and 24.85 Mb of tested sequence, then benchmarked analysis tools and screen design rules.
2. **Enhancer-to-gene mapping.** CRISPRi-FlowFISH, CRISPRi growth screens, and TAP-seq/Perturb-seq can connect perturbed enhancers to target genes. TAP-seq demonstrated perturbation-based enhancer-target maps for 1,778 enhancers.
3. **Variant-to-gene interpretation.** Variants inside CRISPRi-supported regulatory elements should be prioritized when the same locus also has cCRE overlap, accessibility, MPRA allelic effect, rE2G/single-cell linkage, eQTL/QTL, or disease fine-mapping evidence.
4. **Single-cell perturbation programs and GRNs.** Perturb-seq and compressed Perturb-seq use CRISPRi/CRISPR perturbations plus single-cell RNA-seq readouts to infer regulatory circuits, cell states, pathway programs, and genetic interactions.
5. **Multiome perturbation interpretation.** Multiome Perturb-seq extends CRISPRi screens to paired gene expression and chromatin accessibility, enabling analysis of how perturbations alter both transcription and regulatory DNA state.
6. **Screen QC and method comparison.** CRISPRi screen analysis should inspect guide efficiency, low-specificity guides, strand/orientation artifacts, replicate concordance, negative controls, local chromatin context, and whether effects are expression, growth, reporter, or cell-state readouts.

## Literature Anchors

- Multicenter integrated analysis of noncoding CRISPRi screens, Nature Methods 2024: https://www.nature.com/articles/s41592-024-02216-7
- Targeted Perturb-seq enables genome-scale genetic screens in single cells, Nature Methods 2020: https://www.nature.com/articles/s41592-020-0837-5
- Scalable genetic screening for regulatory circuits using compressed Perturb-seq, Nature Biotechnology 2024: https://www.nature.com/articles/s41587-023-01964-9
- Multiome Perturb-seq unlocks scalable discovery of integrated perturbation effects on the transcriptome and epigenome, Cell Systems 2025 / PubMed: https://pubmed.ncbi.nlm.nih.gov/39091800/
- Enhlink infers distal and context-specific enhancer-promoter linkages, Genome Biology 2024: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03374-9

## Retrieval

```bash
python3 Scripts/crispri_data_skills.py pull --source catalog --limit 25
python3 Scripts/crispri_data_skills.py pull --source portal --limit 25
python3 Scripts/crispri_data_skills.py pull --source all --limit 25
```

Use `IGVF_PORTAL_COOKIE` for unreleased IGVF Portal datasets.

## Functional Annotation Integration

```bash
python3 Scripts/crispri_data_skills.py analyze-local --input Data/Input/VariantList/example_variants.csv --label my_locus_crispri
```

The integration score combines CRISPRi significance/effect size, MPRA significance, HepG2 ATAC peak overlap, cCRE annotation, FAVOR predicted function, and base-editing target evidence.

## Suggested IGVFdataAgent Workflow

1. Pull CRISPRi/Perturb-seq metadata from Catalog, Portal, and ENCODE.
2. Build a download manifest before fetching large count, guide, or processed result files.
3. Normalize guide-level or element-level output to stable columns: element interval, guide ID, target gene, biosample, method, effect size, p/q value, direction, source file, and screen readout.
4. Join with cCRE classes, ENCODE-rE2G links, single-cell linkage files, MPRA effects, eQTL/QTL, ATAC/DNase/H3K27ac peaks, and FAVOR/IGVF variant annotations.
5. Summarize hit rates by cCRE class, cell type, target gene, distance-to-TSS, and evidence overlap.
6. Plot effect-size distributions, guide concordance, CRISPRi-vs-MPRA agreement, evidence overlap counts, and IGV-like browser views around prioritized loci.
7. Report research interpretation: candidate causal CREs, likely target genes, cell contexts, phenotype direction, and follow-up experiments.

## Reuse Rules

- Preserve guide/element, target gene, biosample, assay, source fileset, effect size, significance, and direction fields.
- Check whether CRISPRi direction has been synchronized to MPRA or gene loss-of-function assumptions before combining evidence.
- Integrate CRISPRi with MPRA, ATAC/cCRE, eQTL/QTL, enhancer-gene prediction, and base-editing evidence.
- Use manifests before downloading large files.
- Treat CRISPRi evidence as context-specific: a negative result in one cell type does not rule out regulatory function in another context.
