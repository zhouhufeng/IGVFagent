# IGVF Catalog Smoke Analysis

Generated: 2026-05-06 13:12:59 EDT

## Purpose

These public Catalog API checks exercise data classes that are derived from or linked back to IGVF Portal submissions: file/fileset inventory, measured noncoding activity, coding variant scores, and regulatory element-gene evidence.

## catalog_files_filesets

Catalog file/fileset inventory used to connect portal filesets to KG-facing rows.

Endpoint: `/api/files-filesets`
Parameters: `{'limit': '3', 'offset': '0'}`
HTTP status: 200
Returned rows: 3
Saved response: `<repo>/Data/20260506_131258_catalog_files_filesets.json`

## catalog_mpra_variants

Measured noncoding variant activity evidence.

Endpoint: `/api/variants/biosamples`
Parameters: `{'method': 'MPRA', 'limit': '3', 'page': '0'}`
HTTP status: 200
Returned rows: 3
Saved response: `<repo>/Data/20260506_131258_catalog_mpra_variants.json`

## catalog_coding_variant_scores

Coding variant score evidence for one well-known gene.

Endpoint: `/api/genes/coding-variants/scores`
Parameters: `{'gene_name': 'TP53', 'limit': '3', 'page': '0'}`
HTTP status: 200
Returned rows: 25
Saved response: `<repo>/Data/20260506_131258_catalog_coding_variant_scores.json`

## catalog_enhancer_gene_predictions

Regulatory element-gene evidence rows for a small genomic region.

Endpoint: `/api/genomic-elements/genes`
Parameters: `{'region': 'chr1:903900-904900', 'limit': '3', 'page': '0'}`
HTTP status: 200
Returned rows: 3
Saved response: `<repo>/Data/20260506_131259_catalog_enhancer_gene_predictions.json`
