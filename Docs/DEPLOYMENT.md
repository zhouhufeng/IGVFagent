# IGVFdataAgent Deployment Guide

This guide describes how another user can run IGVFdataAgent on a workstation with Codex API, Claude API, or a local Ollama model such as Qwen.

## 1. Clone And Prepare

```bash
git clone <your-repo-url> IGVFdataAgent
cd IGVFdataAgent
mkdir -p Data Scripts Docs/Logs
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
```

The current scripts use the Python standard library where possible. Optional downstream single-cell analysis can add packages such as `pandas`, `scanpy`, `anndata`, `numpy`, `scipy`, `matplotlib`, and `pyarrow`.

## 2. Configure Data Access

Released IGVF Portal and ENCODE data can be accessed without login. Unreleased IGVF Portal data require the user's own authorized browser or cookie/session workflow.

Environment variables:

```bash
export IGVF_CATALOG_API_BASE="https://api.catalogkg.igvf.org"
export ENCODE_BASE="https://www.encodeproject.org"
export IGVF_PORTAL_BASE="https://data.igvf.org"
export FAVOR_API_BASE="https://api.genohub.org"
```

For IGVF Knowledge Graph / ArangoDB access:

```bash
export IGVF_ARANGO_USER="guest"
export IGVF_ARANGO_PASSWORD="<set locally>"
```

Do not commit cookies, API keys, passwords, or `.env` files.

## 3. LLM Backends

IGVFdataAgent scripts are CLI-first, so they can run with any orchestration layer that can call shell commands and read files.

### Codex API

Use Codex as the coding/data agent layer and give it this repository as the workspace root. The agent should be instructed to keep all reads, scripts, data, and logs inside the repository.

Recommended runtime instruction:

```text
Use /path/to/IGVFdataAgent as the project root. Put scripts in Scripts, data in Data, and logs/reports in Docs/Logs or Docs. Use the local CLI skills before writing new one-off code.
```

### Claude API

Use Claude with a tool runner that exposes shell commands in the repo. The same workspace rule applies: restrict file access to the IGVFdataAgent folder, and call the scripts as tools.

Example tool call targets:

```bash
python3 Scripts/igvf_client.py check
python3 Scripts/annotate_variant_list.py --max-rows 10
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
```

### Ollama Local Models

Install Ollama and pull a model:

```bash
ollama pull qwen3
ollama pull llama3.1
```

Run a local model server:

```bash
ollama serve
```

Then connect your preferred local agent runner to `http://localhost:11434`. For workstation use, Qwen-class coding models are useful for command planning and report drafting, while the Python scripts do the deterministic data access and annotation.

## 4. Smoke Test

Run these commands from the repo root:

```bash
python3 Scripts/igvf_client.py catalog-api /
python3 Scripts/igvf_data_skills.py overview --limit 5
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
python3 Scripts/ccre_linkage_annotation_skills.py linkage-manifest --source catalog
```

Expected outputs:

- Runtime logs in `Docs/Logs/`.
- Cached API responses in `Data/`.
- Manifests in `Data/Manifests/`.
- Reports in `Docs/`.

## 5. Typical User Workflows

Annotate a variant list:

```bash
python3 Scripts/annotate_variant_list.py --input Data/Input/VariantList/my_variants.csv
python3 Scripts/ccre_linkage_annotation_skills.py annotate-variants --input Data/Input/VariantList/my_variants.csv
```

Discover and download cCRE resources:

```bash
python3 Scripts/ccre_linkage_annotation_skills.py screen-manifest
python3 Scripts/ccre_linkage_annotation_skills.py screen-download --only PLS --download --max-download-gb 1
```

Discover rE2G and single-cell linkage data:

```bash
python3 Scripts/ccre_linkage_annotation_skills.py linkage-manifest --source all --limit all --hydrate-limit -1
```

Download selected linkage files only after inspecting the generated manifest:

```bash
python3 Scripts/ccre_linkage_annotation_skills.py linkage-download --manifest Data/Manifests/cCRELinkage/<manifest.csv> --only rE2G --download
```

## 6. Workstation Notes

- Keep large downloads on a disk with enough space. Full cCRE, rE2G, and single-cell linkage data can be many gigabytes.
- Prefer manifest commands before download commands.
- Use `--max-rows` for smoke tests.
- Commit scripts and docs, but do not commit large downloaded data unless the project policy explicitly allows it.
- For reproducibility, preserve generated manifest CSVs and Markdown reports.
