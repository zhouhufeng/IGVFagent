# Buckley 2024 — VHL saturation genome editing

## Citation

Buckley M, ..., Findlay GM. *Nat Genet* **56**: 1446–1455 (2024).
DOI: [10.1038/s41588-024-01800-z](https://doi.org/10.1038/s41588-024-01800-z) · PMID: 38969834

## Data

* **MaveDB**: `urn:mavedb:00001183-a-1` (VHL SGE, RCC10)
* **Explorer**: https://vhl-board.onrender.com
* **GitHub**: [FindlayLab/VHL_SGE](https://github.com/FindlayLab/VHL_SGE)

## Headline workflow (paper)

1. SGE of 2,268 VHL SNVs in RCC10 renal cell carcinoma cells.
2. Sort cells over HIF2A stabilisation phenotype, recover variant fitness.
3. Classify each variant into 4 functional classes: LOF1 (severe), LOF2 (mild), Intermediate, Neutral.

## What IGVFagent reproduces

`mavedb map-scoreset --urn urn:mavedb:00001183-a-1 --gene VHL` pulls the MaveDB CSV, resolves the canonical Ensembl transcript for VHL (ENST00000256474 / ENSG00000134086), translates each scored variant into chr / pos / ref / alt, and emits TSV + VCF + summary JSON. Then `catalog get-entity VHL` and `catalog find-associations VHL --relationship genetic` cross-reference the variant set against IGVF Catalog disease + phenotype annotations.

Online-only.

## Ground truth to spot-check

* 2,268 SNVs perturbed (paper)
* Type-2A pheochromocytoma variants (e.g. p.Tyr98His at position 188) ∈ functional bucket [–0.4, –0.22]
* Type-1 truncating LOF1 variants score below –1.26
* Mapping should land on chr3 (VHL is at chr3:10141778-10153667)
