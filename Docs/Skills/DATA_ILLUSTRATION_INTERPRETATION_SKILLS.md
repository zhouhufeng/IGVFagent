# Skill: Data Illustration And Interpretation For IGVF/ENCODE

Use this skill when a user pastes an IGVF Portal or ENCODE accession, object URL, or search URL and needs to understand what the data actually are, which files matter, how to download them, and how to apply them in research.

## Commands

```bash
python3 Scripts/data_illustration_interpretation.py explain '<igvf-portal-url>/curated-sets/IGVFDS2544COZH/'
python3 Scripts/data_illustration_interpretation.py explain '<encode-portal-url>/search/?type=Annotation&searchTerm=encode-re2g&status!=archived'
python3 Scripts/data_illustration_interpretation.py explain IGVFDS2544COZH --download --max-download-gb 2
python3 Scripts/data_illustration_interpretation.py explain '<search-url>' --hydrate-limit 50
```

## Workflow

1. Resolve the pasted ID or URL to JSON metadata using the configured IGVF API base or ENCODE `format=json`.
2. Preserve raw JSON under `Data/Interpreted/Metadata/`.
3. Extract directly linked files and write a download manifest under `Data/Manifests/DataIllustration/`.
4. Generate SVG plots for file formats, content types, and statuses under `Docs/DataIllustration/Plots/`.
5. Write a plain-language report explaining what the object/search appears to represent and how to use it.
6. Only download payloads when explicitly requested with `--download` and a size cap.

## Interpretation Rules

- Translate vague metadata into concrete analysis questions: expression, accessibility, enhancer-gene linkage, variant effect, binding effect, SGE, annotation, or curated collection.
- Prefer processed files for first-pass research use; raw reads are for reprocessing or method development.
- Always check assembly, biosample, assay, status, controlled access, file size, and provenance before joining data sources.
- Suggest IGVF Catalog/KG and ENCODE integrations when the data can support variant-to-gene, enhancer-to-gene, regulatory element, or disease-locus interpretation.
