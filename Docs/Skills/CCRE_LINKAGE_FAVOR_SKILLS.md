# Skill: cCRE, Linkage, FAVOR, And COSMIC-Aware Variant Annotation

Use this skill when users need cCRE downloads from SCREEN, rE2G/single-cell linkage files, FAVOR functional annotations, IGVF Catalog regulatory evidence, or COSMIC-like variant filtering through FAVOR outputs.

## Main Commands

```bash
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
python3 Scripts/ccre_linkage_annotation_skills.py screen-download --only PLS --download --max-download-gb 1
python3 Scripts/ccre_linkage_annotation_skills.py linkage-manifest --source all --limit all --hydrate-limit -1
python3 Scripts/ccre_linkage_annotation_skills.py linkage-download --manifest Data/Manifests/cCRELinkage/<manifest.csv> --download
python3 Scripts/ccre_linkage_annotation_skills.py annotate-variants --input Data/Input/VariantList/example_variants.csv --max-rows 10
python3 Scripts/ccre_linkage_annotation_skills.py cosmic-from-favor --region chr19:44851820-44908922
python3 Scripts/ccre_linkage_annotation_skills.py summarize-local --ccre-bed Data/cCRE/GRCh38-cCREs.bed --linkage-files 'Data/Linkages/*'
python3 Scripts/ccre_linkage_annotation_skills.py browser-demo --region chr19:44850000-44910000
```

## Workflow

1. Run `screen-manifest` first. Review sizes and categories before downloading.
2. Download cCRE class files or cCRE-gene association files with `screen-download`.
3. Run `linkage-manifest --source all --limit all --hydrate-limit -1` to discover ENCODE-rE2G, single-cell linkage, IGVF Portal, and IGVF Catalog linkage evidence. Use a smaller `--hydrate-limit` for smoke tests.
4. Download selected linkage files only after checking the manifest.
5. Annotate user variants with `annotate-variants`, which joins FAVOR fields and IGVF Catalog evidence.
6. Use `cosmic-from-favor` for region-level FAVOR retrieval and COSMIC text filtering. FAVOR docs do not expose a COSMIC-only endpoint, so the skill filters returned variant records.
7. Use `summarize-local` to make cCRE class and linkage evidence plots from downloaded files.
8. Use `browser-demo` to generate an IGV-like static SVG with coordinate axis, cCRE tracks, rE2G/cCRE-gene arcs, and gene models.

## Interpretation

- cCRE classes answer what regulatory element class a variant or linkage lands in.
- rE2G links provide computational enhancer-to-gene predictions across ENCODE biosamples.
- Single-cell linkage files support cell-type-specific enhancer-gene interpretation when available.
- IGVF Catalog edges add experimental, QTL, and prediction evidence around variants, genes, and genomic elements.
- FAVOR adds broad functional annotation, including coding class, nearby genes, population frequency, epigenetic scores, GeneHancer/CAGE/super-enhancer fields, and disease fields.
- IGV-like browser SVGs are for static reporting and quick interpretation; use IGV/JBrowse/UCSC for interactive inspection of very large tracks.
