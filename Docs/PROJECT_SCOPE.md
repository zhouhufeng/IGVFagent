# IGVF Agent Scope

The IGVF agent will provide a local, auditable workflow for discovering and reading data from the IGVF ecosystem.

## Initial Targets

1. IGVF Portal at `https://data.igvf.org/`
   - Requires authenticated access for protected data.
   - Authentication should use a local, user-approved Google/OAuth/browser-session workflow.
   - The repository must not store Gmail passwords or raw session exports.

2. IGVF Catalog docs at `https://docs.catalog.igvf.org/introduction`
   - Public documentation describing Catalog pages, data entities, and query features.
   - The docs identify `https://docs.catalog.igvf.org/llms.txt` as the documentation index.

3. IGVF Catalog API at `https://api.catalogkg.igvf.org/`
   - Preferred programmatic interface for Knowledge Graph data.
   - The Catalog UI's X-ray button can reveal the exact API call behind supported tables.
   - Table-oriented endpoints should be captured as reusable client methods as we discover them.

4. IGVF Knowledge Graph at `https://db.catalog.igvf.org/_db/igvf/`
   - Backed by ArangoDB.
   - Guest access is configured locally with `IGVF_ARANGO_USER` and `IGVF_ARANGO_PASSWORD`.
   - Use direct ArangoDB access when the REST API is insufficient or when we need schema-level inspection.

5. ENCODE Portal at `https://www.encodeproject.org/`
   - Public source for ENCODE experiments, files, biosamples, assays, genome references, schemas, pipelines, and related functional genomics data.
   - The agent should learn to query ENCODE search endpoints, inspect metadata JSON, and download selected files with provenance.

## Agent Responsibilities

- Fetch and cache public documentation.
- Query KG data through the Catalog API.
- Explore and interpret genes, variants, gene-variant associations, coding variant scores, QTL evidence, phenotype links, regulatory element links, and LD summaries.
- Query available Knowledge Graph collections and run controlled AQL queries when needed.
- Query and interpret ENCODE metadata and download selected ENCODE files.
- Build reusable metadata-first workflows for single-cell RNA-seq, single-cell ATAC-seq, and Perturb-seq from IGVF Portal and ENCODE datasets.
- Retrieve and summarize enhancer-gene linkage evidence from IGVF Catalog, IGVF Portal, and ENCODE, including experimental, eQTL/QTL, and computational evidence.
- Retrieve MPRA/STARR/BlueSTARR metadata and evidence, and analyze local MPRA result tables with summary statistics and plots.
- Retrieve CRISPRi/CRISPR FACS/Perturb-seq regulatory perturbation evidence and integrate CRISPRi results with functional annotation evidence.
- Read portal endpoints when authenticated session material is provided locally.
- Save raw outputs in `Data/`.
- Save run logs in `Docs/Logs/`.
- Keep scripts and agent code in `Scripts/`.

## Security Notes

- Never commit `.env`, cookies, OAuth tokens, passwords, or browser profiles.
- Prefer read-only credentials for automation.
- Log request URLs and status codes, but avoid logging credential-bearing headers.
