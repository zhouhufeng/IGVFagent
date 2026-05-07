# IGVF Agent

A local, auditable agent for discovering, retrieving, and analyzing data from
the [IGVF](https://igvf.org/) ecosystem (Portal, Catalog, Knowledge Graph) and
related public resources (ENCODE, FAVOR).

This repository contains command-line skills for:

- IGVF Portal, IGVF Catalog API, IGVF Knowledge Graph (ArangoDB) access
- ENCODE Portal metadata search and file download
- Variant annotation against IGVF Catalog evidence (CADD, QTL, phenotypes,
  regulatory elements, predictions)
- **Advanced variant analysis**: integrates Catalog + ENCODE cCRE evidence
  with optional user experimental tables (CRISPRi / MPRA / GWAS), fits
  logistic models, and produces volcano / Miami / evidence-overlap plots
  plus a research-grade markdown report
- Single-cell RNA-seq, single-cell ATAC-seq, Perturb-seq, and 10x Multiome
  metadata-first workflows
- Enhancer-gene linkage retrieval and comparison (ABC, rE2G, ENCODE-rE2G,
  catalog-based predictions, eQTL-based linkage)
- MPRA / STARR / BlueSTARR retrieval, summary statistics, and plotting
- CRISPRi / CRISPR-FACS / Perturb-seq evidence integration with functional
  annotation
- cCRE (SCREEN) discovery and FAVOR-based variant annotation, plus IGV-like
  browser views
- Data illustration and interpretation across IGVF and ENCODE search URLs

The agent is CLI-first: every skill exposes shell-runnable subcommands so it
can be driven by Codex, Claude, Ollama, or any orchestration layer that can
invoke `python3 Scripts/...`.

## Quick start

```bash
git clone <your-fork-url> IGVFagent
cd IGVFagent
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt   # optional analysis extras
cp .env.example .env                         # then edit locally
```

Smoke test:

```bash
python3 Scripts/igvf_client.py check
python3 Scripts/igvf_data_skills.py overview --limit 5
python3 Scripts/igvf_data_skills.py encode-overview --limit 5
```

See `Scripts/README.md` for the full skill catalog and `Docs/DEPLOYMENT.md`
for environment configuration, LLM-backend wiring, and recommended workflows.

## Repository layout

```
Scripts/   CLI skills and shared client code
Data/      Cached API responses and user-supplied input lists (gitignored)
  Input/VariantList/   User-provided variant CSVs (NOT committed)
Docs/      Project scope, skill documentation, smoke analyses
  Logs/    Runtime logs (gitignored)
  Skills/  Per-skill reference docs
```

## Configuration

The agent reads endpoints and credentials from environment variables:

| Variable                | Purpose                                       |
|-------------------------|-----------------------------------------------|
| `IGVF_CATALOG_API_BASE` | IGVF Catalog API root                         |
| `IGVF_CATALOG_DOCS_BASE`| IGVF Catalog docs root                        |
| `IGVF_PORTAL_BASE`      | IGVF Portal root                              |
| `IGVF_PORTAL_COOKIE`    | Optional: authenticated portal session cookie |
| `IGVF_ARANGO_BASE`      | IGVF Knowledge Graph ArangoDB endpoint        |
| `IGVF_ARANGO_USER`      | KG user (default `guest`)                     |
| `IGVF_ARANGO_PASSWORD`  | KG password (set locally; never commit)       |
| `ENCODE_BASE`           | ENCODE Portal root                            |
| `FAVOR_API_BASE`        | FAVOR API root                                |

See `.env.example`. Never commit `.env`, cookies, OAuth tokens, or browser
session exports.

## Variant lists

The variant-annotation skills accept any user-provided CSV via `--input`. A
tiny illustrative example is shipped at
`Data/Input/VariantList/example_variants.csv`; replace with your own list.
**Do not commit confidential or pre-publication variant data.**

## Security

- The repository ships with an aggressive `.gitignore` that excludes runtime
  outputs, caches, logs, and any `.env` / cookie files.
- All credentials are read from environment variables; nothing is hardcoded.
- Logs record request URLs and HTTP status codes only, not credential headers.

## License

See `LICENSE` (add a license file appropriate to your distribution).
