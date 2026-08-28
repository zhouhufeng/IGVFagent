"""Tool registry for the IGVFagent ReAct runtime.

Curated, hand-described surface of the highest-value subcommands across
every IGVFagent skill. Each tool entry carries:

  * ``name``        — snake_case, stable identifier the LLM will use
  * ``description`` — one-paragraph natural-language description for the
                      model (kept short and concrete)
  * ``parameters``  — JSON Schema for the tool's arguments
  * ``cli``         — list of CLI argv tokens the tool runs under the
                      hood (this skill always shells out to ``igvfagent``
                      so the audit trail and side-effects are identical
                      to the human-driven flow)
  * ``arg_map``     — how to translate parameter dict → CLI args

The registry is intentionally smaller than the union of every skill
subcommand. Surfacing too many low-level tools degrades model planning;
this curated set covers the comprehensive-context, dataset-discovery,
literature, ETL, and pipeline workflows in one or two tool calls each.

``execute(tool_call)`` runs the wrapped CLI as a subprocess, captures
stdout / stderr / exit code, and parses any obvious "Report:" /
"Manifest:" lines so the model can reference downstream artifacts in
its next turn.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)


# --------------------------- Tool dataclass --------------------------------


@dataclasses.dataclass
class Tool:
    name: str
    description: str
    parameters: dict
    cli: "list[str]"               # baseline argv (e.g. ["kg", "gene"])
    positional: "list[str]"        # parameter names that go positional
    flag_map: "dict[str, str]"     # parameter name -> CLI flag name
    flag_repeat: "set[str]"        # parameter names that are list -> repeat flag
    bool_flags: "set[str]"         # parameter names that are bare flags
    # User-extension tools only: full argv of an arbitrary executable to
    # run instead of the ``igvfagent`` console script. Empty for every
    # built-in tool.
    command: "list[str]" = dataclasses.field(default_factory=list)

    def to_dict(self) -> dict:
        return {"name": self.name, "description": self.description,
                "parameters": self.parameters}


# --------------------------- Tool definitions -------------------------------

def _T(name, description, parameters, cli, *,
        positional=(), flag_map=None, flag_repeat=(), bool_flags=()) -> Tool:
    return Tool(
        name=name, description=description, parameters=parameters,
        cli=list(cli), positional=list(positional),
        flag_map=dict(flag_map or {}), flag_repeat=set(flag_repeat),
        bool_flags=set(bool_flags),
    )


_S_OBJECT  = {"type": "object"}
_S_STRING  = {"type": "string"}
_S_INTEGER = {"type": "integer"}
_S_BOOLEAN = {"type": "boolean"}
_S_ARRAY_S = {"type": "array", "items": {"type": "string"}}


_TOOLS: "list[Tool]" = [

    # ──────────────────────────────────────────────────────────────────
    # Self-extension. These are what let the agent add capability instead
    # of only consuming it: it writes a manifest or a Python module into a
    # user-extension directory and the loader absorbs it with no restart.
    # Every write refuses unless IGVF_ALLOW_AGENT_AUTHORING=1, because an
    # authored skill is Python this host later executes as a subprocess.
    # ──────────────────────────────────────────────────────────────────
    _T(
        "annotate_variant_list",
        "★ ANNOTATE A LIST OF VARIANTS PASTED BY THE USER ★ and add them to "
        "the local knowledge graph as variant vertices. Pass the list "
        "VERBATIM in `variants` — any notation and any separator works "
        "(2-21001846-G-A, chr2:21001846:G:A, rs763341676, SPDI, HGVS, VCF "
        "lines; spaces, commas or newlines). THIS IS THE RIGHT TOOL whenever "
        "someone pastes variants and asks to annotate, characterise, or "
        "add them to the KG — do NOT ask them for a CSV file first. Queries "
        "FAVOR (CADD, ClinVar, GENCODE, enhancer/promoter overlap) and the "
        "IGVF Catalog, writes a report plus CSV/JSON, and upserts variant "
        "nodes with edges to genes and diseases.",
        {
            "type": "object",
            "properties": {
                "variants":  {**_S_STRING, "description":
                              "The variant list verbatim, any notation."},
                "input":     {**_S_STRING, "description":
                              "Alternative: path to a file of variants."},
                "label":     {**_S_STRING, "description":
                              "Short label for the run directory."},
                "sources":   {**_S_STRING, "description":
                              "Comma list: favor, catalog (default both)."},
                "max_rows":  {**_S_INTEGER, "description": "0 = all."},
                "no_kg":     {**_S_BOOLEAN, "description":
                              "Annotate without writing to the KG."},
            },
        },
        ["variant-list", "annotate"],
        flag_map={"variants": "--variants", "input": "--input",
                   "label": "--label", "sources": "--sources",
                   "max_rows": "--max-rows", "no_kg": "--no-kg"},
        bool_flags=("no_kg",),
    ),

    _T(
        "document_plan",
        "★ READ AN UPLOADED PAPER AND PLAN ITS REPRODUCTION ★. Give it the "
        "path of a PDF/DOCX/text manuscript (uploads land in "
        "`Data/Uploads/`). Extracts the text, finds every repository "
        "accession (GEO, IGVF, ENCODE, MaveDB, Synapse, dbGaP…) and assay "
        "family, and returns a concrete `igvfagent bench` chain to reproduce "
        "the analysis. THIS IS THE RIGHT TOOL whenever someone uploads a "
        "paper and asks to reproduce, replicate or analyse it.",
        {
            "type": "object",
            "properties": {
                "path":  {**_S_STRING, "description":
                          "Path to the manuscript, e.g. Data/Uploads/paper.pdf"},
                "label": {**_S_STRING, "description": "Short run label."},
            },
            "required": ["path"],
        },
        ["document", "plan"],
        flag_map={"path": "--path", "label": "--label"},
    ),

    _T(
        "document_read",
        "★ Extract the text of an uploaded document ★ (PDF/DOCX/text) so it "
        "can be quoted or searched. Use document_plan instead when the goal "
        "is to reproduce the paper's analysis. Reports honestly when a PDF "
        "has no text layer rather than returning an empty result.",
        {
            "type": "object",
            "properties": {
                "path": {**_S_STRING, "description": "Path to the document."},
                "head": {**_S_INTEGER, "description":
                         "Print only the first N characters."},
            },
            "required": ["path"],
        },
        ["document", "read"],
        flag_map={"path": "--path", "head": "--head"},
    ),

    _T(
        "read_artifact",
        "★ READ A FILE THIS AGENT PRODUCED ★ — reports, manifests, JSON, CSV, "
        "TSV, logs. Skills announce outputs as `Report: <path>`; use this to "
        "OPEN that path and quote what is actually inside. ALWAYS call this "
        "before summarising a report: without it you only know the counts a "
        "skill happened to print, so you would say '18 diseases' instead of "
        "naming them. Use --head/--tail to page a large file.",
        {
            "type": "object",
            "properties": {
                "path":      {**_S_STRING, "description":
                              "Workspace-relative or absolute path to read."},
                "head":      {**_S_INTEGER, "description":
                              "Return only the first N lines."},
                "tail":      {**_S_INTEGER, "description":
                              "Return only the last N lines."},
                "max_bytes": {**_S_INTEGER, "description":
                              "Byte cap (default 200000)."},
            },
            "required": ["path"],
        },
        ["artifact", "read"],
        flag_map={"path": "--path", "head": "--head", "tail": "--tail",
                   "max_bytes": "--max-bytes"},
    ),

    _T(
        "grep_artifacts",
        "★ SEARCH INSIDE produced artefacts ★ for a regex — find which report "
        "mentions a gene, disease, or accession without reading each one. "
        "Returns file, line number, and the matching line.",
        {
            "type": "object",
            "properties": {
                "pattern":  {**_S_STRING, "description": "Regex to search for."},
                "path":     {**_S_STRING, "description":
                             "Directory or file to search (default Docs)."},
                "max_hits": {**_S_INTEGER, "description": "Cap on hits (default 50)."},
            },
            "required": ["pattern"],
        },
        ["artifact", "grep"],
        flag_map={"pattern": "--pattern", "path": "--path",
                   "max_hits": "--max-hits"},
    ),

    _T(
        "list_artifacts",
        "★ LIST a run directory ★ to discover what a skill actually wrote "
        "before reading it. Use on the run dir a skill reported.",
        {
            "type": "object",
            "properties": {
                "path":  {**_S_STRING, "description": "Directory (default Docs)."},
                "limit": {**_S_INTEGER, "description": "Max entries (default 200)."},
            },
        },
        ["artifact", "ls"],
        flag_map={"path": "--path", "limit": "--limit"},
    ),

    _T(
        "ext_author_tool",
        "★ WRAP AN EXISTING COMMAND as a new tool ★. Use ONLY when the "
        "command already exists — an igvfagent subcommand (--cli) or a real "
        "executable like `cat` or `sort` (--command). It is REJECTED if the "
        "--cli subcommand does not exist, because such a tool registers "
        "cleanly and then fails on every call. **If the capability does not "
        "exist yet, use ext_author_skill instead** — it writes the "
        "implementation AND registers the tool in one step. Give `parameters` "
        "as a JSON Schema object string, or the model cannot pass it input.",
        {
            "type": "object",
            "properties": {
                "name":        {**_S_STRING, "description":
                                "snake_case tool name, 3-49 chars."},
                "description": {**_S_STRING, "description":
                                "What it does and when to call it — this is "
                                "all a model sees when choosing it."},
                "cli":         {**_S_STRING, "description":
                                'igvfagent subcommand tail, e.g. "kg gene".'},
                "command":     {**_S_STRING, "description":
                                "argv for any executable (alternative to cli)."},
                "parameters":  {**_S_STRING, "description":
                                'JSON Schema object, e.g. {"type":"object",'
                                '"properties":{"gene":{"type":"string"}}}'},
                "positional":  {**_S_STRING, "description":
                                "space-separated params passed positionally."},
            },
            "required": ["name", "description"],
        },
        ["extauthor", "write-tool"],
        flag_map={"name": "--name", "description": "--description",
                   "cli": "--cli", "command": "--command",
                   "parameters": "--parameters", "positional": "--positional"},
    ),

    _T(
        "ext_author_skill",
        "★ WRITE A NEW CAPABILITY FROM SCRATCH ★ — the right tool whenever "
        "the user asks for something IGVFagent cannot currently do. Writes a "
        "Python module AND registers a matching callable tool in one step, so "
        "you can invoke it immediately afterwards. Pass the COMPLETE module "
        "source in `source` (must define a top-level main(); it is "
        "syntax-checked before writing) and declare `tool_parameters` as a "
        "JSON Schema object so the new tool can receive arguments. Prefer "
        "this over ext_author_tool for any new logic — parsing, analysis, a "
        "new API client. Then call ext_validate to confirm it loaded.",
        {
            "type": "object",
            "properties": {
                "name":        {**_S_STRING, "description":
                                "snake_case skill name; registers as "
                                "`igvfagent <name-with-hyphens>`."},
                "description": {**_S_STRING, "description": "What it does."},
                "source":      {**_S_STRING, "description":
                                "Full Python source, must define main()."},
                "source_file": {**_S_STRING, "description":
                                "Alternative: read source from this path."},
                "tool_parameters": {**_S_STRING, "description":
                                "JSON Schema object for the auto-registered "
                                'tool, e.g. {"type":"object","properties":'
                                '{"variants":{"type":"string"}}}. Without it '
                                "the new tool takes no arguments."},
            },
            "required": ["name", "description"],
        },
        ["extauthor", "write-skill"],
        flag_map={"name": "--name", "description": "--description",
                   "source": "--source", "source_file": "--source-file",
                   "tool_parameters": "--tool-parameters"},
    ),

    _T(
        "ext_validate",
        "★ Check an authored extension actually registered ★ and, for a "
        "skill, that it imports. Call this right after ext_author_skill — a "
        "module that fails to import is skipped silently by the loader.",
        {
            "type": "object",
            "properties": {
                "name": {**_S_STRING, "description": "Extension name."},
            },
            "required": ["name"],
        },
        ["extauthor", "validate"],
        flag_map={"name": "--name"},
    ),

    _T(
        "ext_list",
        "★ List every user-authored tool and skill ★ currently discovered, "
        "the directories searched, whether authoring is enabled, and any "
        "manifests that failed to load. Call this before authoring to avoid "
        "duplicating something that already exists.",
        {"type": "object", "properties": {}},
        ["extauthor", "list"],
    ),

    _T(
        "kg_gene",
        "Comprehensive gene context from the IGVF Catalog KG: variants, "
        "transcripts, proteins, regulatory elements, diseases, pathways. "
        "Default 'tell me about gene X' tool. Optional flags pull FAVOR, "
        "enhancer-gene linkage, single-cell datasets, prior literature.",
        {
            "type": "object",
            "properties": {
                "symbol":   {**_S_STRING, "description": "Gene symbol, e.g. APOE"},
                "depth":    {**_S_INTEGER, "default": 1,
                              "description": "1 = direct relations only; 2 = also fan out per variant."},
                "limit":    {**_S_INTEGER, "default": 25},
                "call_favor":      {**_S_BOOLEAN, "default": False},
                "call_linkage":    {**_S_BOOLEAN, "default": False},
                "call_singlecell": {**_S_BOOLEAN, "default": False},
                "call_literature": {**_S_BOOLEAN, "default": False},
                "literature_context": {**_S_ARRAY_S,
                    "description": "Free-text disease / tissue terms for the literature side-call."},
                "label": {**_S_STRING, "description": "Run label for output dir."},
            },
            "required": ["symbol"],
        },
        cli=["kg", "gene"],
        positional=["symbol"],
        flag_map={
            "depth": "--depth", "limit": "--limit",
            "literature_context": "--literature-context",
            "label": "--label",
        },
        flag_repeat={"literature_context"},
        bool_flags={"call_favor", "call_linkage",
                     "call_singlecell", "call_literature"},
    ),

    _T(
        "kg_variant",
        "Variant-centric KG: linked genes, regulatory elements, "
        "phenotypes, biosamples, predictions. Accepts rsID/SPDI/HGVS. "
        "Spans several Catalog collections, each with its own IGVF "
        "'method' vocabulary: variant->gene effects are 'Variant-EFFECTS', "
        "variant->phenotype functional calls are 'cV2F', and "
        "variant->biosample reporter records are 'STARR-seq' / "
        "'BlueSTARR'. For variant effects on protein binding "
        "specifically (SEMVAR / ADASTRA allele-specific binding), use "
        "grn_protein_variants — that collection is not traversed here.",
        {
            "type": "object",
            "properties": {
                "variant": {**_S_STRING,
                    "description": "rsID (rs429358), SPDI (NC_000019.10:44908821:C:T), or HGVS"},
                "limit":   {**_S_INTEGER, "default": 25},
                "call_favor":      {**_S_BOOLEAN, "default": False},
                "call_literature": {**_S_BOOLEAN, "default": False},
                "literature_context": {**_S_ARRAY_S},
                "label": {**_S_STRING},
            },
            "required": ["variant"],
        },
        cli=["kg", "variant"],
        positional=["variant"],
        flag_map={"limit": "--limit",
                   "literature_context": "--literature-context",
                   "label": "--label"},
        flag_repeat={"literature_context"},
        bool_flags={"call_favor", "call_literature"},
    ),

    _T(
        "kg_region",
        "Region-centric KG: genes + cCREs + enhancer-gene linkage in "
        "window. Format chr19:44903000-44912000.",
        {
            "type": "object",
            "properties": {
                "region": {**_S_STRING},
                "limit":  {**_S_INTEGER, "default": 50},
                "call_favor": {**_S_BOOLEAN, "default": False},
                "label":  {**_S_STRING},
            },
            "required": ["region"],
        },
        cli=["kg", "region"],
        positional=["region"],
        flag_map={"limit": "--limit", "label": "--label"},
        bool_flags={"call_favor"},
    ),

    _T(
        "explain_dataset",
        "Plain-language explainer for an IGVF/ENCODE accession or URL: "
        "metadata, file inventory, SVG overview plots, how-to-use report.",
        {
            "type": "object",
            "properties": {
                "accession_or_url": {**_S_STRING},
                "download": {**_S_BOOLEAN, "default": False},
                "max_download_gb": {"type": "number", "default": 0.3},
            },
            "required": ["accession_or_url"],
        },
        cli=["explain", "explain"],
        positional=["accession_or_url"],
        flag_map={"max_download_gb": "--max-download-gb"},
        bool_flags={"download"},
    ),

    _T(
        "annotate_variants",
        "Annotate a variant CSV against IGVF Catalog evidence "
        "(CADD/QTL/phenotypes/regulatory). CSV needs rsid/hgvs/spdi or "
        "chr/pos/ref/alt.",
        {
            "type": "object",
            "properties": {
                "input": {**_S_STRING, "description": "Path to variant CSV."},
                "max_rows": {**_S_INTEGER, "default": 25},
            },
            "required": ["input"],
        },
        cli=["variant"],
        flag_map={"input": "--input", "max_rows": "--max-rows"},
    ),

    _T(
        "advanced_variant_analysis",
        "Integrated variant pipeline: catalog evidence + cCRE overlay + "
        "Predicted_Functional composite + optional experimental join + "
        "logistic model + markdown report with volcano/Miami/overlap plots.",
        {
            "type": "object",
            "properties": {
                "input":         {**_S_STRING},
                "experimental":  {**_S_STRING},
                "outcome":       {**_S_STRING},
                "gene_list":     {**_S_STRING},
                "label":         {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["advanced-variant", "run"],
        flag_map={"input": "--input", "experimental": "--experimental",
                   "outcome": "--outcome", "gene_list": "--gene-list",
                   "label": "--label"},
    ),

    _T(
        "splitseq_retrieve",
        "★ Discover IGVF Parse SPLiT-seq AnalysisSets by lab/tissue/taxa. "
        "Writes a manifest with donor + founder-strain metadata.",
        {
            "type": "object",
            "properties": {
                "limit":       {**_S_INTEGER, "default": 50},
                "lab":         {**_S_STRING},
                "tissue":      {**_S_STRING},
                "taxa":        {**_S_STRING},
                "label":       {**_S_STRING},
                "fetch_file_details": {**_S_BOOLEAN, "default": False},
            },
        },
        cli=["splitseq", "retrieve"],
        flag_map={"limit": "--limit", "lab": "--lab",
                   "tissue": "--sample-type", "taxa": "--taxa",
                   "label": "--label"},
        bool_flags={"fetch_file_details"},
    ),

    _T(
        "splitseq_manifest",
        "Per-pool/per-donor file manifest for one or more SPLiT-seq "
        "AnalysisSet accessions.",
        {
            "type": "object",
            "properties": {
                "accessions": {**_S_STRING,
                    "description": "Comma-separated AnalysisSet accessions."},
                "label": {**_S_STRING},
            },
            "required": ["accessions"],
        },
        cli=["splitseq", "manifest"],
        flag_map={"accessions": "--accessions", "label": "--label"},
    ),

    _T(
        "multiome_retrieve",
        "Discover IGVF 10x Multiome AnalysisSets; writes file/sample/donor "
        "manifest.",
        {
            "type": "object",
            "properties": {
                "count": {**_S_INTEGER, "default": 25},
                "fetch_file_details": {**_S_BOOLEAN, "default": False},
            },
        },
        cli=["multiome", "retrieve"],
        flag_map={"count": "--count"},
        bool_flags={"fetch_file_details"},
    ),

    _T(
        "enhancer_gene_overview",
        "Pull enhancer-gene linkage (ABC/rE2G/eQTL) from Catalog or ENCODE.",
        {
            "type": "object",
            "properties": {
                "source": {**_S_STRING, "enum": ["catalog", "encode", "all"],
                            "default": "catalog"},
                "limit":  {**_S_INTEGER, "default": 25},
            },
        },
        cli=["enhancer", "overview"],
        flag_map={"source": "--source", "limit": "--limit"},
    ),

    _T(
        "mpra_pull",
        "Pull MPRA/STARR/BlueSTARR metadata from Catalog or Portal. "
        "Backed by /api/variants/biosamples, where the reporter assay is "
        "carried by the Catalog 'method' field: 'STARR-seq' for the "
        "measured assay and 'BlueSTARR' for the neural-network predictions "
        "over it (source='IGVF'). Those are distinct record sets — say "
        "which one the user means rather than conflating them. This "
        "collection is region-queryable (e.g. chr4:155600-155770).",
        {
            "type": "object",
            "properties": {
                "source": {**_S_STRING, "enum": ["catalog", "portal"],
                            "default": "catalog"},
                "limit":  {**_S_INTEGER, "default": 25},
            },
        },
        cli=["mpra", "pull"],
        flag_map={"source": "--source", "limit": "--limit"},
    ),

    _T(
        "mpra_activity",
        "Per-oligo MPRA activity: negative-binomial GLM Wald test (DESeq2 "
        "via pydeseq2) of RNA vs DNA counts, followed by summit-shift "
        "normalization that moves the mode of the log2-fold-change density "
        "to zero (Tewhey-lab MPRAmodel convention). Emits a .out table with "
        "baseMean, log2FoldChange, lfcSE, stat, pvalue, padj per Oligo. "
        "Input: a counts table with an 'Oligo' column and DNA_*/RNA_* "
        "replicate columns; barcode-level tables (with a 'Barcode' column) "
        "are summed per oligo automatically.",
        {
            "type": "object",
            "properties": {
                "input": {**_S_STRING, "description":
                          "Path to counts table (CSV/TSV)."},
                "label": {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["mpra", "activity"],
        flag_map={"input": "--input", "label": "--label"},
    ),

    _T(
        "mpra_skew",
        "Allelic-skew analysis: pairs ref/alt oligos by "
        "SNP_window_strand_haplotype, computes per-replicate "
        "log2((RNA+1)/(mean DNA+1)) for each allele, runs a paired t-test "
        "of alt vs ref across replicates, and applies BH-FDR. Emits a .out "
        "table with Log2Skew, LogSkew_SE, tstat, pvalue, padj per element. "
        "Mirrors the t-test path of Tewhey-lab MPRAmodel runSkew.",
        {
            "type": "object",
            "properties": {
                "input": {**_S_STRING, "description":
                          "Counts table with Allele plus at least one of "
                          "SNP/Window/Strand/Haplotype."},
                "label": {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["mpra", "skew"],
        flag_map={"input": "--input", "label": "--label"},
    ),

    _T(
        "mpra_qc",
        "MPRA QC suite: replicate concordance (pairwise Pearson r on log10 "
        "counts, separately for DNA and RNA, with an N x N SVG heatmap), "
        "unique-barcodes-per-oligo histogram (when a 'Barcode' column is "
        "present), and total-counts-per-oligo histogram (log10). Writes a "
        "markdown QC report alongside the SVGs.",
        {
            "type": "object",
            "properties": {
                "input": {**_S_STRING, "description":
                          "Counts table (barcode-level supported)."},
                "label": {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["mpra", "qc"],
        flag_map={"input": "--input", "label": "--label"},
    ),

    _T(
        "mpra_volcano",
        "Render a 4-panel volcano figure (activity volcano, allelic-skew "
        "volcano, activity MA plot, skew MA plot) from the .out tables "
        "produced by `mpra activity` and `mpra skew`. Points with padj "
        "below the FDR threshold are highlighted in red.",
        {
            "type": "object",
            "properties": {
                "activity": {**_S_STRING, "description":
                              "Activity .out table from `mpra activity`."},
                "skew":     {**_S_STRING, "description":
                              "Skew .out table from `mpra skew`."},
                "label":    {**_S_STRING},
                "title":    {**_S_STRING},
                "fdr":      {"type": "number", "default": 0.05},
            },
        },
        cli=["mpra", "volcano"],
        flag_map={"activity": "--activity", "skew": "--skew",
                   "label": "--label", "title": "--title", "fdr": "--fdr"},
    ),

    # ──────────────────────────────────────────────────────────────────
    # 10x Multiome analytics — clean-room reimplementations of methods in
    # 10XGenomics/analysis_guides + Stuart-lab Signac (MIT) extensions.
    # ──────────────────────────────────────────────────────────────────
    _T(
        "multiome_qc_atac",
        "★ Per-barcode ATAC QC from a fragments TSV/BED ★. Computes "
        "fragments per barcode, TSS enrichment (Signac TSSEnrichment "
        "formula: center reads / max(0.2, flank reads), ±100 bp center "
        "/ ±900-1000 bp flank), nucleosome signal (mono-nucleosome / "
        "NFR fragment-length ratio), and FRIP (reads in peaks). Use "
        "this BEFORE joint-qc.",
        {
            "type": "object",
            "properties": {
                "fragments": {**_S_STRING, "description":
                              "Path to fragments TSV/BED (.gz OK)."},
                "tss_bed":   {**_S_STRING, "description":
                              "Optional TSS BED for TSS enrichment."},
                "peaks_bed": {**_S_STRING, "description":
                              "Optional peaks BED for FRIP."},
                "label":     {**_S_STRING},
            },
            "required": ["fragments"],
        },
        cli=["multiome", "qc-atac"],
        flag_map={"fragments": "--fragments", "tss_bed": "--tss-bed",
                   "peaks_bed": "--peaks-bed", "label": "--label"},
    ),
    _T(
        "multiome_joint_qc",
        "Merge ATAC + RNA per-barcode QC tables and apply Signac-"
        "convention joint thresholds (RNA: 1k-25k UMIs, ≥200 genes, "
        "≤20%% MT; ATAC: 1.8k-100k frags, TSS>1, nucleosome<2, "
        "FRIP>0.15). Outputs a per-barcode QC TSV with 'qc' label "
        "(both / RNA only / ATAC only / neither).",
        {
            "type": "object",
            "properties": {
                "rna_qc":   {**_S_STRING},
                "atac_qc":  {**_S_STRING},
                "min_umis": {**_S_INTEGER, "default": 1000},
                "max_umis": {**_S_INTEGER, "default": 25000},
                "min_genes": {**_S_INTEGER, "default": 200},
                "max_pct_mt": {"type": "number", "default": 0.20},
                "min_frags": {**_S_INTEGER, "default": 1800},
                "max_frags": {**_S_INTEGER, "default": 100000},
                "min_tss":  {"type": "number", "default": 1.0},
                "max_nuc":  {"type": "number", "default": 2.0},
                "min_frip": {"type": "number", "default": 0.15},
                "label":    {**_S_STRING},
            },
            "required": ["rna_qc", "atac_qc"],
        },
        cli=["multiome", "joint-qc"],
        flag_map={"rna_qc": "--rna-qc", "atac_qc": "--atac-qc",
                   "min_umis": "--min-umis", "max_umis": "--max-umis",
                   "min_genes": "--min-genes", "max_pct_mt": "--max-pct-mt",
                   "min_frags": "--min-frags", "max_frags": "--max-frags",
                   "min_tss": "--min-tss", "max_nuc": "--max-nuc",
                   "min_frip": "--min-frip", "label": "--label"},
    ),
    _T(
        "multiome_lsi",
        "TF-IDF normalization + truncated SVD on a peak × cell matrix. "
        "Drops dim 1 (correlates with sequencing depth) per Signac "
        "RunSVD + DepthCor convention. Writes the embedding to "
        "obsm['X_lsi'] of the input h5ad and saves a new h5ad.",
        {
            "type": "object",
            "properties": {
                "input": {**_S_STRING},
                "n_components": {**_S_INTEGER, "default": 50},
                "seed":  {**_S_INTEGER, "default": 7},
                "label": {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["multiome", "lsi"],
        flag_map={"input": "--input", "n_components": "--n-components",
                   "seed": "--seed", "label": "--label"},
    ),
    _T(
        "multiome_wnn",
        "Joint weighted-nearest-neighbor embedding via muon (BSD-3). "
        "Equivalent to Seurat 5 FindMultiModalNeighbors. Requires both "
        "AnnData inputs to share per-cell barcodes; ATAC must have "
        "obsm['X_lsi'] (run `multiome lsi` first). Optionally clusters "
        "the joint graph with Leiden.",
        {
            "type": "object",
            "properties": {
                "rna_h5ad":  {**_S_STRING},
                "atac_h5ad": {**_S_STRING},
                "n_pca":     {**_S_INTEGER, "default": 50},
                "n_neighbors": {**_S_INTEGER, "default": 20},
                "min_dist":  {"type": "number", "default": 0.3},
                "cluster":   {**_S_BOOLEAN, "default": False},
                "resolution": {"type": "number", "default": 1.0},
                "seed":      {**_S_INTEGER, "default": 7},
                "label":     {**_S_STRING},
            },
            "required": ["rna_h5ad", "atac_h5ad"],
        },
        cli=["multiome", "wnn"],
        flag_map={"rna_h5ad": "--rna-h5ad", "atac_h5ad": "--atac-h5ad",
                   "n_pca": "--n-pca", "n_neighbors": "--n-neighbors",
                   "min_dist": "--min-dist", "cluster": "--cluster",
                   "resolution": "--resolution", "seed": "--seed",
                   "label": "--label"},
    ),
    _T(
        "multiome_peak2gene",
        "Peak-to-gene correlation. For each peak, identify genes whose "
        "TSS is within --window bp (default 500 kb) and compute "
        "Pearson/Spearman correlation between peak accessibility and "
        "gene expression across cells. BH-FDR adjusted p-values. Sign "
        "of correlation indicates enhancer-like (positive) vs "
        "repressor-like (negative) association.",
        {
            "type": "object",
            "properties": {
                "rna_h5ad":  {**_S_STRING},
                "atac_h5ad": {**_S_STRING},
                "tss_bed":   {**_S_STRING},
                "window":    {**_S_INTEGER, "default": 500000},
                "method":    {**_S_STRING, "default": "pearson",
                               "enum": ["pearson", "spearman"]},
                "max_pairs": {**_S_INTEGER, "default": 200000},
                "label":     {**_S_STRING},
            },
            "required": ["rna_h5ad", "atac_h5ad", "tss_bed"],
        },
        cli=["multiome", "peak2gene"],
        flag_map={"rna_h5ad": "--rna-h5ad", "atac_h5ad": "--atac-h5ad",
                   "tss_bed": "--tss-bed", "window": "--window",
                   "method": "--method", "max_pairs": "--max-pairs",
                   "label": "--label"},
    ),
    _T(
        "multiome_da_peaks",
        "Differential accessibility on TF-IDF-normalized peaks via "
        "Wilcoxon rank-genes. Pairs naturally with peak2gene — peak2gene "
        "gives enhancer→gene candidates; da-peaks tells you WHICH peaks "
        "are differentially open between user-defined clusters / "
        "cell-types.",
        {
            "type": "object",
            "properties": {
                "input":       {**_S_STRING, "description":
                                 "ATAC h5ad (peaks will be TF-IDF normalized)."},
                "cluster_key": {**_S_STRING, "default": "leiden_wnn"},
                "top_n":       {**_S_INTEGER, "default": 50},
                "label":       {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["multiome", "da-peaks"],
        flag_map={"input": "--input", "cluster_key": "--cluster-key",
                   "top_n": "--top-n", "label": "--label"},
    ),
    _T(
        "multiome_atac_spectral",
        "Jaccard-Laplacian spectral embedding alternative to LSI for "
        "ATAC peak matrices (snapATAC2-style). Robust to depth "
        "differences without explicit depth-correction drop. Good for "
        "datasets where TF-IDF + LSI gives a strong depth axis.",
        {
            "type": "object",
            "properties": {
                "input":          {**_S_STRING},
                "n_components":   {**_S_INTEGER, "default": 30},
                "n_neighbors":    {**_S_INTEGER, "default": 20},
                "max_cells":      {**_S_INTEGER, "default": 5000},
                "seed":           {**_S_INTEGER, "default": 7},
                "label":          {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["multiome", "atac-spectral"],
        flag_map={"input": "--input", "n_components": "--n-components",
                   "n_neighbors": "--n-neighbors", "max_cells": "--max-cells",
                   "seed": "--seed", "label": "--label"},
    ),
    _T(
        "multiome_chromvar",
        "Clean-room chromVAR-style TF motif activity per cell. For each "
        "(cell, motif) pair, computes raw deviations from expected "
        "accessibility, then bias-corrects via K=50 GC-content + log-mean-"
        "accessibility-matched background motif sets, yielding a per-cell "
        "z-score for every motif. Output: cells × motifs z-score TSV.",
        {
            "type": "object",
            "properties": {
                "input":       {**_S_STRING, "description":
                                 "ATAC h5ad (peak × cell)."},
                "motif_hits":  {**_S_STRING, "description":
                                 "Peak × motif binary TSV."},
                "gc_content":  {**_S_STRING, "description":
                                 "Optional per-peak GC content TSV."},
                "k_background": {**_S_INTEGER, "default": 50},
                "seed":        {**_S_INTEGER, "default": 7},
                "label":       {**_S_STRING},
            },
            "required": ["input", "motif_hits"],
        },
        cli=["multiome", "chromvar"],
        flag_map={"input": "--input", "motif_hits": "--motif-hits",
                   "gc_content": "--gc-content",
                   "k_background": "--k-background",
                   "seed": "--seed", "label": "--label"},
    ),
    _T(
        "multiome_css",
        "Cluster Similarity Spectrum batch correction (He 2020 Genome "
        "Biol). Per-batch HVG → per-batch Leiden → per-batch cluster "
        "centroids → represent each cell as its vector of correlations "
        "to all batch×cluster centroids. Apache/BSD-friendly alternative "
        "to GPL Harmony.",
        {
            "type": "object",
            "properties": {
                "input":       {**_S_STRING},
                "batch_key":   {**_S_STRING},
                "n_hvg":       {**_S_INTEGER, "default": 2000},
                "n_pca":       {**_S_INTEGER, "default": 50},
                "n_neighbors": {**_S_INTEGER, "default": 20},
                "resolution":  {"type": "number", "default": 1.0},
                "label":       {**_S_STRING},
            },
            "required": ["input", "batch_key"],
        },
        cli=["multiome", "css"],
        flag_map={"input": "--input", "batch_key": "--batch-key",
                   "n_hvg": "--n-hvg", "n_pca": "--n-pca",
                   "n_neighbors": "--n-neighbors",
                   "resolution": "--resolution", "label": "--label"},
    ),
    _T(
        "multiome_multivi",
        "MultiVI deep joint VAE (Ashuach 2023) via scvi-tools. "
        "Optional dep — install with `pip install scvi-tools` to enable. "
        "Joint generative model over RNA (ZINB) + ATAC (Bernoulli on "
        "binarized peaks) yielding a shared latent z, with batch-key-"
        "conditioned encoder/decoder.",
        {
            "type": "object",
            "properties": {
                "rna_h5ad":  {**_S_STRING},
                "atac_h5ad": {**_S_STRING},
                "batch_key": {**_S_STRING, "default": "batch"},
                "epochs":    {**_S_INTEGER, "default": 50},
                "label":     {**_S_STRING},
            },
            "required": ["rna_h5ad", "atac_h5ad"],
        },
        cli=["multiome", "multivi"],
        flag_map={"rna_h5ad": "--rna-h5ad", "atac_h5ad": "--atac-h5ad",
                   "batch_key": "--batch-key", "epochs": "--epochs",
                   "label": "--label"},
    ),

    _T(
        "multiome_showcase",
        "★ ONE-COMMAND 10x MULTIOME QC SHOWCASE ★. Runs qc-atac + "
        "(optional) joint-qc + builds a 6-panel composite figure "
        "(fragments/cell, TSS enrichment, nucleosome signal, FRIP, "
        "reads in TSS, reads in peaks — each with Signac thresholds "
        "as red dashed lines) + writes a narrative report. THIS IS "
        "THE RIGHT TOOL FOR ANY 10x MULTIOME QC DEMO QUESTION.",
        {
            "type": "object",
            "properties": {
                "fragments": {**_S_STRING},
                "tss_bed":   {**_S_STRING},
                "peaks_bed": {**_S_STRING},
                "rna_qc":    {**_S_STRING, "description":
                               "Optional RNA QC TSV (output of share rna-qc)."},
                "label":     {**_S_STRING},
            },
            "required": ["fragments"],
        },
        cli=["multiome", "showcase"],
        flag_map={"fragments": "--fragments", "tss_bed": "--tss-bed",
                   "peaks_bed": "--peaks-bed", "rna_qc": "--rna-qc",
                   "label": "--label"},
    ),

    _T(
        "crispri_pull",
        "Pull CRISPRi/CRISPR-FACS/Perturb-seq evidence from the Catalog. "
        "Backed by /api/genes/genomic-elements, where the IGVF assay that "
        "produced an element->gene link is carried by the Catalog 'method' "
        "field — 'Perturb-seq' and 'CRISPR screen' are the two values on "
        "this collection (source='IGVF'). For per-edge effect sizes and "
        "significance from the same perturbation experiments, use "
        "grn_network instead; this tool returns the element/gene evidence "
        "records, not the dEx statistics.",
        {
            "type": "object",
            "properties": {
                "source": {**_S_STRING, "enum": ["catalog", "portal"],
                            "default": "catalog"},
                "limit":  {**_S_INTEGER, "default": 25},
            },
        },
        cli=["crispri", "pull"],
        flag_map={"source": "--source", "limit": "--limit"},
    ),

    _T(
        "grn_network",
        "Differential-expression gene regulatory network (dEx GRN) edges "
        "from the IGVF Catalog. Call this for 'what genes does TF X "
        "regulate?', 'what regulates gene Y?', 'is there a Perturb-seq "
        "self-edge for TF Z?', or any question about regulator->target "
        "effect size or significance. Returns one row per element->gene "
        "edge with log2FC, neg_log10_pvalue, significant, crispr_modality "
        "and the perturbed element's coordinates. Query from either end: "
        "'regulator' gives that TF's targets, 'response' gives that gene's "
        "regulators, and passing BOTH with the same symbol tests the "
        "self-edge (many TFs legitimately have none). This is the only "
        "tool that reaches /api/gene-regulatory-network — kg_gene and "
        "catalog_find_associations do not cover it.",
        {
            "type": "object",
            "properties": {
                "regulator": {**_S_STRING, "description":
                    "Regulator (TF) gene symbol; returns its targets. "
                    "e.g. EOMES"},
                "response": {**_S_STRING, "description":
                    "Response gene symbol; returns its regulators. "
                    "e.g. ARHGEF3"},
                "method": {**_S_STRING, "description":
                    "Assay method: 'Perturb-seq' or 'CRISPR screen'. Omit "
                    "for all methods."},
                "p_value": {"type": "number", "description":
                    "Keep edges with p <= this (typical: 0.05). Omit to "
                    "return all edges regardless of significance."},
                "host": {**_S_STRING, "enum": ["prod", "dev"],
                          "default": "prod",
                          "description":
                    "'dev' queries the pre-release demo Catalog; use only "
                    "when the user explicitly asks for dev data."},
                "max_results": {**_S_INTEGER, "default": 10000},
                "label": {**_S_STRING},
            },
        },
        cli=["grn", "network"],
        flag_map={"regulator": "--regulator", "response": "--response",
                   "method": "--method", "p_value": "--p-value",
                   "host": "--host", "max_results": "--max-results",
                   "label": "--label"},
    ),

    _T(
        "grn_protein_variants",
        "Sequence-variant effects on proteins from the IGVF Catalog: "
        "allele-specific binding and motif-disruption calls (SEMVAR, "
        "ADASTRA). Use for 'which variants alter binding of TF X?' or "
        "'what is the allele-specific binding evidence for this protein?'. "
        "Returns sequence_variant, protein_complex, biosample_term, the "
        "effect label (e.g. 'allele-specific binding' / 'binding modulated "
        "by'), method and source_url. This is the only tool that reaches "
        "/api/proteins/variants.",
        {
            "type": "object",
            "properties": {
                "protein": {**_S_STRING, "description":
                    "Protein / gene symbol, e.g. ELF2, TP53."},
                "method": {**_S_STRING, "description":
                    "Scoring method, e.g. 'SEMVAR' or 'ADASTRA'."},
                "source": {**_S_STRING, "description":
                    "Originating source, e.g. IGVF, ADASTRA."},
                "host": {**_S_STRING, "enum": ["prod", "dev"],
                          "default": "prod"},
                "max_results": {**_S_INTEGER, "default": 10000},
                "label": {**_S_STRING},
            },
        },
        cli=["grn", "protein-variants"],
        flag_map={"protein": "--protein", "method": "--method",
                   "source": "--source", "host": "--host",
                   "max_results": "--max-results", "label": "--label"},
    ),

    _T(
        "ccre_screen_manifest",
        "Build a manifest of SCREEN cCRE downloads (PLS/pELS/dELS/CTCF).",
        {"type": "object", "properties": {}},
        cli=["ccre", "screen-manifest"],
    ),

    _T(
        "ccre_favor",
        "Annotate variants in a region with FAVOR (CADD/GERP/conservation).",
        {
            "type": "object",
            "properties": {
                "region": {**_S_STRING},
            },
            "required": ["region"],
        },
        cli=["ccre", "cosmic-from-favor"],
        flag_map={"region": "--region"},
    ),

    _T(
        "ref_learn",
        "Multi-source literature search (PubMed/bioRxiv/arXiv/SemanticScholar/"
        "OpenAlex). Ranks by journal+citations; extracts methods + plot "
        "vocabulary; emits a consensus figure recipe.",
        {
            "type": "object",
            "properties": {
                "topic": {**_S_STRING},
                "limit": {**_S_INTEGER, "default": 25},
                "top":   {**_S_INTEGER, "default": 12},
                "label": {**_S_STRING},
            },
            "required": ["topic"],
        },
        cli=["ref", "learn"],
        flag_map={"topic": "--topic", "limit": "--limit",
                   "top": "--top", "label": "--label"},
    ),

    _T(
        "ref_validate",
        "Cross-check a CSV of genes/variants/regulatory elements against "
        "prior literature.",
        {
            "type": "object",
            "properties": {
                "input":   {**_S_STRING},
                "context": {**_S_ARRAY_S,
                    "description": "Disease / tissue / phenotype context terms."},
                "limit_per_item": {**_S_INTEGER, "default": 5},
                "label":   {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["ref", "validate"],
        flag_map={"input": "--input", "context": "--context",
                   "limit_per_item": "--limit-per-item",
                   "label": "--label"},
        flag_repeat={"context"},
    ),

    _T(
        "ref_design",
        "★ Recommend a study workflow + cognate published studies + matching "
        "IGVF Portal AnalysisSets for an assay type.",
        {
            "type": "object",
            "properties": {
                "data_type":   {**_S_STRING,
                    "enum": ["10x_multiome", "parse_split_seq", "mpra",
                              "crispri", "enhancer_gene"]},
                "assay_title": {**_S_STRING},
                "label":       {**_S_STRING},
            },
            "required": ["data_type"],
        },
        cli=["ref", "design"],
        flag_map={"data_type": "--data-type",
                   "assay_title": "--assay-title",
                   "label": "--label"},
    ),

    _T(
        "portal_kg_pull",
        "ETL: pull Portal entities by tissue/gene/assay/lab into the local "
        "SQLite KG; expands linked samples/donors/files at depth 1.",
        {
            "type": "object",
            "properties": {
                "type":   {**_S_STRING,
                    "enum": ["AnalysisSet", "MeasurementSet", "Sample",
                              "Donor", "File", "FileSet", "Document"],
                    "default": "AnalysisSet"},
                "tissue": {**_S_STRING},
                "gene":   {**_S_STRING},
                "assay":  {**_S_STRING},
                "lab":    {**_S_STRING},
                "limit":  {**_S_INTEGER, "default": 50},
                "depth":  {**_S_INTEGER, "default": 1},
            },
        },
        cli=["portal-kg", "pull"],
        flag_map={"type": "--type", "tissue": "--tissue",
                   "gene": "--gene", "assay": "--assay", "lab": "--lab",
                   "limit": "--limit", "depth": "--depth"},
    ),

    _T(
        "portal_kg_annotate",
        "Mine local-KG Portal node descriptions for gene/variant mentions, "
        "confirm via Catalog, add mentions_gene / mentions_variant edges.",
        {"type": "object", "properties": {}},
        cli=["portal-kg", "annotate"],
    ),

    _T(
        "portal_kg_enrich",
        "For each Gene node in the local KG, hydrate Catalog evidence "
        "(variants, transcripts, cCREs, diseases, pathways).",
        {
            "type": "object",
            "properties": {
                "symbols": {**_S_STRING,
                    "description": "Comma-separated gene symbols. Default = all."},
                "limit":   {**_S_INTEGER, "default": 25},
            },
        },
        cli=["portal-kg", "enrich"],
        flag_map={"symbols": "--symbols", "limit": "--limit"},
    ),

    _T(
        "sce2g_kg_pull",
        "★ BULK-ingest scE2G element→gene regulatory linkages from the IGVF "
        "Catalog into the local KG as `regulates` edges. Deterministic, "
        "resumable, runs to completion over the whole genome regardless of "
        "size (adaptive region tiling around the API's 500-row cap) with a "
        "progress heartbeat. USE THIS for 'download/integrate all scE2G "
        "linkages' — one call, not a per-gene loop. Scope with `region` "
        "(default 'all' = whole genome) or `chromosomes`.",
        {
            "type": "object",
            "properties": {
                "region":      {**_S_STRING, "description":
                    "'all' (whole genome, default), a chromosome ('chr19'), "
                    "or a locus ('chr19:44900000-45000000')."},
                "chromosomes": {**_S_STRING, "description":
                    "Comma-separated chroms to restrict 'all' (e.g. '19,20,X')."},
                "min_window":  {**_S_INTEGER, "default": 20000},
                "heartbeat":   {**_S_INTEGER, "default": 25},
                "max_windows": {**_S_INTEGER, "default": 0,
                    "description": "Stop after N windows (0 = unlimited)."},
            },
        },
        cli=["sce2g-kg", "pull"],
        flag_map={"region": "--region", "chromosomes": "--chromosomes",
                   "min_window": "--min-window", "heartbeat": "--heartbeat",
                   "max_windows": "--max-windows"},
    ),

    _T(
        "portal_kg_query",
        "Query the local KG by gene / tissue / node-id.",
        {
            "type": "object",
            "properties": {
                "gene":    {**_S_STRING},
                "tissue":  {**_S_STRING},
                "node_id": {**_S_STRING},
                "limit":   {**_S_INTEGER, "default": 50},
            },
        },
        cli=["portal-kg", "query"],
        flag_map={"gene": "--gene", "tissue": "--tissue",
                   "node_id": "--node-id", "limit": "--limit"},
    ),

    _T(
        "portal_kg_stats",
        "Local KG counts (node/edge/source) + recent run history.",
        {"type": "object", "properties": {}},
        cli=["portal-kg", "stats"],
    ),

    _T(
        "frontpage_summary",
        "Refresh IGVF Portal + KG front-page summary stats.",
        {
            "type": "object",
            "properties": {
                "update_readme": {**_S_BOOLEAN, "default": False},
            },
        },
        cli=["frontpage", "refresh"],
        bool_flags={"update_readme"},
    ),

    # ---- ENCODE pipeline tools (step 6) ----

    _T(
        "encode_retrieve",
        "Search ENCODE for experiments by assay (ChIP-seq / Histone / "
        "ATAC-seq / DNase / Hi-C / capture Hi-C / ChIA-PET / RNA-seq / "
        "MNase / FAIRE / CAGE / RAMPAGE), biosample, target, assembly.",
        {
            "type": "object",
            "properties": {
                "assay":       {**_S_STRING,
                    "description": "Assay group; e.g. 'Histone ChIP-seq', "
                                    "'ATAC-seq', 'Hi-C'."},
                "biosample":   {**_S_STRING,
                    "description": "Biosample term name (K562, GM12878, "
                                    "liver, hippocampus, etc.)."},
                "target":      {**_S_STRING,
                    "description": "ChIP-seq / Histone target (H3K27ac, "
                                    "CTCF, etc.). Ignored for non-ChIP."},
                "assembly":    {**_S_STRING},
                "limit":       {**_S_INTEGER, "default": 50},
                "fetch_file_details": {**_S_BOOLEAN, "default": False},
                "label":       {**_S_STRING},
            },
            "required": ["assay"],
        },
        cli=["encode", "retrieve"],
        flag_map={"assay": "--assay", "biosample": "--biosample",
                   "target": "--target", "assembly": "--assembly",
                   "limit": "--limit", "label": "--label"},
        bool_flags={"fetch_file_details"},
    ),

    _T(
        "encode_describe",
        "Plain-language report for one ENCODE experiment: assay, biosample, "
        "target, replicates, file inventory.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING,
                    "description": "ENCODE accession (e.g. ENCSR000DUB)."},
            },
            "required": ["accession"],
        },
        cli=["encode", "describe"],
        flag_map={"accession": "--accession"},
    ),

    _T(
        "encode_super_enhancers",
        "ROSE-style super-enhancer call from an enhancer-mark peak BED "
        "(H3K27ac/BRD4/MED1/P300). Stitches + ranks + inflection point.",
        {
            "type": "object",
            "properties": {
                "bed":                 {**_S_STRING},
                "stitching_distance":  {**_S_INTEGER, "default": 12500},
                "tss_bed":             {**_S_STRING,
                    "description": "Optional TSS BED for excluding "
                                    "promoter-proximal peaks."},
                "tss_distance":        {**_S_INTEGER, "default": 2000},
                "label":               {**_S_STRING},
            },
            "required": ["bed"],
        },
        cli=["encode", "super-enhancers"],
        flag_map={"bed": "--bed",
                   "stitching_distance": "--stitching-distance",
                   "tss_bed": "--tss-bed", "tss_distance": "--tss-distance",
                   "label": "--label"},
    ),

    _T(
        "encode_integrate_ccre",
        "Overlay a peak BED with SCREEN cCRE classes (PLS/pELS/dELS/CTCF). "
        "Auto-downloads the cCRE registry on first use.",
        {
            "type": "object",
            "properties": {
                "bed":         {**_S_STRING},
                "ccre_bed":    {**_S_STRING,
                    "description": "Optional local cCRE BED override."},
                "ccre_class":  {**_S_STRING,
                    "enum": ["PLS", "pELS", "dELS", "CTCF"]},
                "label":       {**_S_STRING},
            },
            "required": ["bed"],
        },
        cli=["encode", "integrate-ccre"],
        flag_map={"bed": "--bed", "ccre_bed": "--ccre-bed",
                   "ccre_class": "--ccre-class", "label": "--label"},
    ),

    _T(
        "encode_browser",
        "IGV-style multi-track SVG for a region. LABEL:PATH BED tracks + "
        "optional SCREEN cCRE overlay. Set `re2g='auto'` (and optionally a "
        "gene filter) to also draw enhancer→gene linkage arcs from the "
        "IGVF Catalog (rE2G / ENCODE-rE2G / ABC). Use this for any user "
        "request that mentions 'rE2G links', 'enhancer-gene arcs', or "
        "'browser view of gene X with its regulatory elements'.",
        {
            "type": "object",
            "properties": {
                "region": {**_S_STRING,
                    "description": "chr19:44903000-44912000"},
                "tracks": {**_S_ARRAY_S,
                    "description": "List of `LABEL:PATH` strings for "
                                    "each BED track to render."},
                "with_ccre": {**_S_BOOLEAN, "default": False},
                "ccre_bed":  {**_S_STRING},
                "width":     {**_S_INTEGER, "default": 1000},
                "label":     {**_S_STRING},
                "re2g":      {**_S_STRING,
                    "description": "Pass 'auto' to pull enhancer→gene "
                                    "links from the IGVF Catalog, or a "
                                    "path to a BEDPE / CSV linkage file."},
                "re2g_gene_filter": {**_S_STRING,
                    "description": "Comma-separated gene symbols to "
                                    "restrict arcs (e.g. 'BRCA1' or "
                                    "'APOE,TOMM40')."},
                "re2g_score_cut": {"type": "number", "default": 0.0,
                    "description": "Minimum rE2G/ABC score to render."},
                "re2g_limit":     {**_S_INTEGER, "default": 200},
                "re2g_arcs_height": {**_S_INTEGER, "default": 140},
            },
            "required": ["region"],
        },
        cli=["encode", "browser"],
        flag_map={"region": "--region", "tracks": "--track",
                   "ccre_bed": "--ccre-bed", "width": "--width",
                   "label": "--label",
                   "re2g": "--re2g",
                   "re2g_gene_filter": "--re2g-gene-filter",
                   "re2g_score_cut": "--re2g-score-cut",
                   "re2g_limit": "--re2g-limit",
                   "re2g_arcs_height": "--re2g-arcs-height"},
        flag_repeat={"tracks"},
        bool_flags={"with_ccre"},
    ),

    _T(
        "encode_bigwig_frip",
        "FRiP (fraction-of-signal-in-peaks) from bigWig + peak BED. "
        "Requires [hic] extras.",
        {
            "type": "object",
            "properties": {
                "bigwig": {**_S_STRING},
                "bed":    {**_S_STRING},
                "label":  {**_S_STRING},
            },
            "required": ["bigwig", "bed"],
        },
        cli=["encode", "bigwig-frip"],
        flag_map={"bigwig": "--bigwig", "bed": "--bed", "label": "--label"},
    ),

    _T(
        "encode_bigwig_tss_heatmap",
        "bigWig signal heatmap centered on anchors (TSS BED / peak "
        "summits) + meta-profile. Requires [hic] extras.",
        {
            "type": "object",
            "properties": {
                "bigwig":      {**_S_STRING},
                "anchor_bed":  {**_S_STRING},
                "window":      {**_S_INTEGER, "default": 4000},
                "bins":        {**_S_INTEGER, "default": 200},
                "max_anchors": {**_S_INTEGER, "default": 5000},
                "label":       {**_S_STRING},
            },
            "required": ["bigwig", "anchor_bed"],
        },
        cli=["encode", "bigwig-tss-heatmap"],
        flag_map={"bigwig": "--bigwig", "anchor_bed": "--anchor-bed",
                   "window": "--window", "bins": "--bins",
                   "max_anchors": "--max-anchors", "label": "--label"},
    ),

    _T(
        "encode_hic_matrix",
        "Hi-C contact heatmap for a region from .mcool or .hic. "
        "Requires [hic] extras.",
        {
            "type": "object",
            "properties": {
                "input":      {**_S_STRING,
                    "description": "Path to .mcool / .hic file."},
                "region":     {**_S_STRING,
                    "description": "chr19:44900000-45100000"},
                "resolution": {**_S_INTEGER, "default": 10000},
                "balance":    {**_S_BOOLEAN, "default": False},
                "label":      {**_S_STRING},
            },
            "required": ["input", "region"],
        },
        cli=["encode", "hic-matrix"],
        flag_map={"input": "--input", "region": "--region",
                   "resolution": "--resolution", "label": "--label"},
        bool_flags={"balance"},
    ),

    _T(
        "encode_hic_insulation",
        "Crane-style insulation score for TAD-boundary calls from "
        ".mcool/.hic. Requires [hic] extras.",
        {
            "type": "object",
            "properties": {
                "input":      {**_S_STRING},
                "region":     {**_S_STRING},
                "resolution": {**_S_INTEGER, "default": 10000},
                "window":     {**_S_INTEGER, "default": 200000},
                "balance":    {**_S_BOOLEAN, "default": False},
                "boundary_threshold": {"type": "number", "default": -0.3},
                "label":      {**_S_STRING},
            },
            "required": ["input", "region"],
        },
        cli=["encode", "hic-insulation"],
        flag_map={"input": "--input", "region": "--region",
                   "resolution": "--resolution", "window": "--window",
                   "boundary_threshold": "--boundary-threshold",
                   "label": "--label"},
        bool_flags={"balance"},
    ),

    _T(
        "encode_loops_analyze",
        "Loop QC from .bedpe (Hi-C / ChIA-PET / capture Hi-C): length "
        "distribution, intra/inter split, optional anchor-peak overlap.",
        {
            "type": "object",
            "properties": {
                "bedpe":  {**_S_STRING},
                "peaks":  {**_S_ARRAY_S,
                    "description": "Optional `LABEL:PATH` peak BEDs to "
                                    "intersect against the loop anchors."},
                "label":  {**_S_STRING},
            },
            "required": ["bedpe"],
        },
        cli=["encode", "loops-analyze"],
        flag_map={"bedpe": "--bedpe", "peaks": "--peaks",
                   "label": "--label"},
        flag_repeat={"peaks"},
    ),

    _T(
        "encode_motif_enrichment",
        "TF motif enrichment in peak sequences (CTCF/AP-1/GATA1/ETS/NFkB/"
        "STAT1/FOXA1/TP53/MYC/SP1) vs shuffled background. Requires "
        "--genome FASTA and [motif] extras.",
        {
            "type": "object",
            "properties": {
                "bed":          {**_S_STRING},
                "genome":       {**_S_STRING,
                    "description": "Path to indexed genome FASTA, e.g. "
                                    "hg38.fa or hg38.fa.gz."},
                "top":          {**_S_INTEGER, "default": 2000},
                "score_cutoff": {"type": "number", "default": 8.0},
                "label":        {**_S_STRING},
            },
            "required": ["bed", "genome"],
        },
        cli=["encode", "motif-enrichment"],
        flag_map={"bed": "--bed", "genome": "--genome", "top": "--top",
                   "score_cutoff": "--score-cutoff", "label": "--label"},
    ),

    _T(
        "se_targets_pipeline",
        "★ End-to-end super-enhancer → target-gene pipeline. For a chosen "
        "biosample (GM12878 / K562 / HepG2 / liver / brain / ...): "
        "discovers H3K27ac (or BRD4/MED1/P300) ChIP-seq + optional "
        "Hi-C/ChIA-PET 3D experiments, downloads peak BED, calls "
        "super-enhancers ROSE-style, and links each SE to candidate "
        "target genes via four streams (3D loops, IGVF Catalog rE2G/ABC "
        "predictions, proximity, SCREEN cCRE composition). Outputs a "
        "ranked SE↔gene table, network plot, and report.",
        {
            "type": "object",
            "properties": {
                "biosample":  {**_S_STRING,
                    "description": "ENCODE biosample term (GM12878, K562, "
                                    "liver, hippocampus, ...)."},
                "target":     {**_S_STRING, "default": "H3K27ac"},
                "assembly":   {**_S_STRING, "default": "GRCh38"},
                "include_3d": {**_S_BOOLEAN, "default": False},
                "gene":       {**_S_STRING,
                    "description": "Optional gene to focus a per-locus "
                                    "browser view on."},
                "label":      {**_S_STRING},
            },
            "required": ["biosample"],
        },
        cli=["se-targets", "pipeline"],
        flag_map={"biosample": "--biosample", "target": "--target",
                   "assembly": "--assembly", "gene": "--gene",
                   "label": "--label"},
        bool_flags={"include_3d"},
    ),

    _T(
        "se_targets_discover",
        "List candidate ChIP-seq + DNase/ATAC + Hi-C/ChIA-PET "
        "experiments for a biosample without downloading or analyzing "
        "anything. Use this to scope what's available before running "
        "the full pipeline.",
        {
            "type": "object",
            "properties": {
                "biosample": {**_S_STRING},
                "target":    {**_S_STRING, "default": "H3K27ac"},
                "assembly":  {**_S_STRING, "default": "GRCh38"},
                "limit":     {**_S_INTEGER, "default": 25},
                "label":     {**_S_STRING},
            },
            "required": ["biosample"],
        },
        cli=["se-targets", "discover"],
        flag_map={"biosample": "--biosample", "target": "--target",
                   "assembly": "--assembly", "limit": "--limit",
                   "label": "--label"},
    ),

    _T(
        "geo_search",
        "Search NCBI GEO Series by keyword / organism / platform "
        "(e.g. 'GM12878 RNA-seq' organism 'Homo sapiens').",
        {
            "type": "object",
            "properties": {
                "query":    {**_S_STRING},
                "organism": {**_S_STRING},
                "platform": {**_S_STRING},
                "study_type": {**_S_STRING},
                "limit":    {**_S_INTEGER, "default": 25},
                "label":    {**_S_STRING},
            },
            "required": ["query"],
        },
        cli=["geo", "search"],
        flag_map={"query": "--query", "organism": "--organism",
                   "platform": "--platform", "study_type": "--study-type",
                   "limit": "--limit", "label": "--label"},
    ),

    _T(
        "geo_series",
        "Pull metadata + sample sheet for one GEO Series accession "
        "(e.g. GSE9574). Writes a markdown report and a per-GSM CSV "
        "ready to feed to rnaseq_pipeline.",
        {
            "type": "object",
            "properties": {
                "gse":           {**_S_STRING,
                    "description": "GSE accession, e.g. GSE9574."},
                "full_samples":  {**_S_BOOLEAN, "default": False},
                "label":         {**_S_STRING},
            },
            "required": ["gse"],
        },
        cli=["geo", "series"],
        flag_map={"gse": "--gse", "label": "--label"},
        bool_flags={"full_samples"},
    ),

    _T(
        "geo_download",
        "Download supplementary / matrix files for a GEO Series. "
        "Filter by category (matrix / suppl / soft) and a regex over "
        "filenames; cap by --max-download-gb.",
        {
            "type": "object",
            "properties": {
                "gse":     {**_S_STRING},
                "only":    {**_S_ARRAY_S,
                    "description": "Subset of [matrix, suppl, soft]."},
                "pattern": {**_S_STRING,
                    "description": "Case-insensitive regex over filenames."},
                "max_download_gb": {"type": "number", "default": 1.0},
            },
            "required": ["gse"],
        },
        cli=["geo", "download"],
        flag_map={"gse": "--gse", "only": "--only",
                   "pattern": "--pattern",
                   "max_download_gb": "--max-download-gb"},
        flag_repeat={"only"},
    ),

    _T(
        "rnaseq_pipeline",
        "★ End-to-end bulk RNA-seq: QC + PCA + differential expression "
        "(pyDESeq2 if installed, Welch's t-test + BH FDR fallback) + "
        "volcano/MA/heatmap plots + DEG → controlling cCRE linkage "
        "via the IGVF Catalog. Inputs: counts matrix + sample sheet + "
        "two group labels in --condition-col.",
        {
            "type": "object",
            "properties": {
                "counts":         {**_S_STRING},
                "sample_sheet":   {**_S_STRING},
                "condition_col":  {**_S_STRING, "default": "condition"},
                "group_a":        {**_S_STRING,
                    "description": "Control / reference group label."},
                "group_b":        {**_S_STRING,
                    "description": "Treated / test group label."},
                "padj_cut":       {"type": "number", "default": 0.05},
                "fc_cut":         {"type": "number", "default": 1.0},
                "skip_link_cre":  {**_S_BOOLEAN, "default": False},
                "label":          {**_S_STRING},
            },
            "required": ["counts", "sample_sheet", "group_a", "group_b"],
        },
        cli=["rnaseq", "pipeline"],
        flag_map={"counts": "--counts", "sample_sheet": "--sample-sheet",
                   "condition_col": "--condition-col",
                   "group_a": "--group-a", "group_b": "--group-b",
                   "padj_cut": "--padj-cut", "fc_cut": "--fc-cut",
                   "label": "--label"},
        bool_flags={"skip_link_cre"},
    ),

    _T(
        "rnaseq_link_cre",
        "Given a DEG CSV (with `gene`, `log2FC`, `padj` columns), query "
        "the IGVF Catalog for the regulatory elements that control "
        "each significant gene. Useful as a follow-up to an existing "
        "DEG analysis.",
        {
            "type": "object",
            "properties": {
                "deg":             {**_S_STRING},
                "padj_cut":        {"type": "number", "default": 0.05},
                "fc_cut":          {"type": "number", "default": 1.0},
                "limit_per_gene":  {**_S_INTEGER, "default": 10},
                "max_genes":       {**_S_INTEGER, "default": 50},
                "label":           {**_S_STRING},
            },
            "required": ["deg"],
        },
        cli=["rnaseq", "link-cre"],
        flag_map={"deg": "--deg", "padj_cut": "--padj-cut",
                   "fc_cut": "--fc-cut",
                   "limit_per_gene": "--limit-per-gene",
                   "max_genes": "--max-genes", "label": "--label"},
    ),

    _T(
        "proteomics_download",
        "Download latest PPI/pathway sources (BioGRID, IntAct, HuRI, "
        "Reactome, KEGG, UniProt idmap) and IGVF Portal protein assays. "
        "Maintains _versions.json and only re-fetches changed releases.",
        {
            "type": "object",
            "properties": {
                "source": {**_S_STRING,
                            "description": "biogrid|intact|huri|reactome|kegg|"
                                           "igvf|uniprot|all (comma-separated)"},
                "biogrid_version": {**_S_STRING,
                            "description": "Pin a BioGRID release (e.g. 4.4.244)."},
                "kegg_max_pathways": {**_S_INTEGER, "default": 400},
            },
            "required": ["source"],
        },
        cli=["proteomics", "download"],
        flag_map={"source": "--source",
                   "biogrid_version": "--biogrid-version",
                   "kegg_max_pathways": "--kegg-max-pathways"},
    ),

    _T(
        "proteomics_versions",
        "Show locally-installed PPI/pathway source versions and probe "
        "upstream-latest. Use before deciding to update.",
        {"type": "object", "properties": {}},
        cli=["proteomics", "versions"],
    ),

    _T(
        "proteomics_igvf_protein",
        "Pull all IGVF Portal protein-assay metadata + actual PPI / "
        "stability files (semi-qY2H, DUAL-IPA, VAMP-seq) into "
        "Data/Proteomics/Sources/IGVF/.",
        {"type": "object", "properties": {}},
        cli=["proteomics", "igvf-protein"],
    ),

    _T(
        "proteomics_build_kg",
        "Build the local SQLite proteomics knowledge graph from "
        "previously-downloaded sources. Edges deduped on "
        "(id_a, id_b, source, source_id).",
        {
            "type": "object",
            "properties": {
                "sources":  {**_S_STRING, "default": "all",
                              "description": "biogrid,intact,huri,reactome,"
                                             "kegg,igvf,uniprot or 'all'."},
                "max_rows": {**_S_INTEGER, "default": 0,
                              "description": "Cap rows per source (0 = no cap)."},
            },
        },
        cli=["proteomics", "build-kg"],
        flag_map={"sources": "--sources", "max_rows": "--max-rows"},
    ),

    _T(
        "proteomics_kg_stats",
        "Summary statistics on the integrated proteomics PPI-KG: total "
        "interactions, distinct proteins, per-source / per-evidence-type / "
        "per-detection-method breakdowns, top hubs.",
        {"type": "object",
         "properties": {"label": {**_S_STRING}}},
        cli=["proteomics", "kg-stats"],
        flag_map={"label": "--label"},
    ),

    _T(
        "proteomics_kg_visualize",
        "Generate degree distribution, top-hubs, per-source breakdown, "
        "and (when --gene given) an ego graph PNG. Saves under "
        "Docs/Proteomics/<ts>_<label>/Plots/.",
        {
            "type": "object",
            "properties": {
                "gene":           {**_S_STRING,
                                    "description": "Symbol/UniProt for ego graph."},
                "max_neighbors":  {**_S_INTEGER, "default": 60},
                "label":          {**_S_STRING},
            },
        },
        cli=["proteomics", "kg-visualize"],
        flag_map={"gene": "--gene",
                   "max_neighbors": "--max-neighbors",
                   "label": "--label"},
    ),

    _T(
        "proteomics_assay_survey",
        "Use the Reference skill to retrieve recent Nature/Cell/Science "
        "studies on VAMP-seq (MultiSTEP), VAMP-seq, MAVE, semi-qY2H, and "
        "DUAL-IPA. Writes literature_survey.md/.json.",
        {
            "type": "object",
            "properties": {
                "label":          {**_S_STRING},
                "max_per_assay":  {**_S_INTEGER, "default": 20},
            },
        },
        cli=["proteomics", "assay-survey"],
        flag_map={"label": "--label",
                   "max_per_assay": "--max-per-assay"},
    ),

    _T(
        "proteomics_assay_figures",
        "Generate per-assay example histograms from the IGVF Portal "
        "files for VAMP-seq family, MAVE, semi-qY2H v1/v2/v3, DUAL-IPA. "
        "Requires `proteomics igvf-protein` first.",
        {"type": "object",
         "properties": {"label": {**_S_STRING}}},
        cli=["proteomics", "assay-figures"],
        flag_map={"label": "--label"},
    ),

    _T(
        "proteomics_pipeline",
        "End-to-end: download sources → build KG → stats → visualize → "
        "per-assay figures → literature survey. Use --skip-download to "
        "reuse cached files; --gene to add an ego graph.",
        {
            "type": "object",
            "properties": {
                "sources":  {**_S_STRING, "default": "all"},
                "label":    {**_S_STRING, "default": "pipeline"},
                "gene":     {**_S_STRING},
                "max_rows": {**_S_INTEGER, "default": 0},
                "skip_download":       {**_S_BOOLEAN, "default": False},
                "skip_literature":     {**_S_BOOLEAN, "default": False},
                "skip_assay_figures":  {**_S_BOOLEAN, "default": False},
            },
        },
        cli=["proteomics", "pipeline"],
        flag_map={"sources": "--sources", "label": "--label",
                   "gene": "--gene", "max_rows": "--max-rows"},
        bool_flags={"skip_download", "skip_literature", "skip_assay_figures"},
    ),

    _T(
        "proteomics_vampseq_pull",
        "★ Download canonical published VAMP-seq scoresets from MaveDB "
        "(PTEN, TPMT, VKOR, PRKN, CYP2C9, NUDT15). Pulls per-replicate "
        "score CSVs for downstream analysis.",
        {
            "type": "object",
            "properties": {
                "gene": {**_S_STRING,
                          "description": "Gene symbol; omit to pull all curated "
                                         "VAMP-seq targets."},
            },
        },
        cli=["proteomics", "vampseq-pull"],
        flag_map={"gene": "--gene"},
    ),

    _T(
        "proteomics_vampseq_analyze",
        "Deep VAMP-seq analysis on a MaveDB scoreset. Produces the "
        "Fowler-lab-style suite of plots: score-density distribution, "
        "residue×AA heatmap, per-residue mean (± IQR) AND per-residue "
        "median + 3-residue moving average, replicate concordance "
        "(rep-1 vs rep-2) AND N×N replicate matrix when ≥ 3 reps, "
        "abundance-class bar, cumulative ranked variants WITH 95 % CI "
        "band, nonsense-by-position QC scatter, biophysical-feature "
        "Spearman ρ panel (when RSA / B-factor / hydrophobicity / "
        "Grantham / PSIC columns are present), and an optional "
        "PyMOL .pml export when `pdb_id` is supplied.",
        {
            "type": "object",
            "properties": {
                "gene":   {**_S_STRING,
                            "description": "PTEN | TPMT | VKOR | PRKN | CYP2C9 | "
                                           "NUDT15. Omit for all targets."},
                "label":  {**_S_STRING},
                "pdb_id": {**_S_STRING,
                            "description": "Optional PDB id (e.g. 1d5r for "
                                           "PTEN, 2bzg for TPMT). When set, "
                                           "emits a PyMOL .pml that overlays "
                                           "the abundance scale on the "
                                           "structure."},
            },
        },
        cli=["proteomics", "vampseq-analyze"],
        flag_map={"gene": "--gene", "label": "--label",
                   "pdb_id": "--pdb-id"},
    ),

    _T(
        "proteomics_vampseq_showcase",
        "★ ONE-COMMAND COMPREHENSIVE VAMP-SEQ DEMO ★. Auto-downloads the "
        "MaveDB scoreset for a target gene (PTEN by default; also TPMT, "
        "VKOR, PRKN, CYP2C9, NUDT15), runs the full Fowler-lab 10-plot "
        "suite (score distribution, residue×AA heatmap, per-residue "
        "mean+IQR + median+moving-avg, replicate concordance + N×N "
        "matrix, abundance-class bar, nonsense-by-position QC, "
        "cumulative ranked variants with 95% CI, biophysical-feature "
        "correlations) AND a publication-grade 9-panel composite figure "
        "AND a deep narrative report explaining how to read every plot. "
        "USE THIS TOOL for any 'show me VAMP-seq', 'analyze VAMP-seq for "
        "<gene>', 'demonstrate variant abundance' question — do NOT call "
        "proteomics_vampseq_pull + proteomics_vampseq_analyze separately.",
        {
            "type": "object",
            "properties": {
                "gene":  {**_S_STRING,
                          "description": "PTEN | TPMT | VKOR | PRKN | CYP2C9 | NUDT15"},
                "label": {**_S_STRING},
                "pdb_id": {**_S_STRING},
            },
        },
        cli=["proteomics", "vampseq-showcase"],
        flag_map={"gene": "--gene", "label": "--label", "pdb_id": "--pdb-id"},
    ),

    _T(
        "proteomics_vampseq_inventory",
        "Decode IGVF Portal raw VAMP-seq MeasurementSets into a "
        "tile × bin × replicate × antibody coverage matrix per target gene "
        "(currently F9 from the MultiSTEP deposit, plus CYP2C19 / G6PD).",
        {
            "type": "object",
            "properties": {
                "label": {**_S_STRING},
            },
        },
        cli=["proteomics", "vampseq-inventory"],
        flag_map={"label": "--label"},
    ),

    _T(
        "sc_pipeline",
        "★ Full single-cell analysis pipeline: QC + filter, log-normalize, "
        "HVG selection, PCA, k-NN graph, UMAP, t-SNE, Leiden clustering, "
        "and marker-gene DE. Accepts .h5ad / 10x .h5 / .mtx / .csv / .tsv. "
        "Saves processed.h5ad, markers.csv, and publication PNGs (UMAP "
        "colored by cluster / sample / top markers). The headline tool "
        "when a user asks for a UMAP / t-SNE / cluster analysis.",
        {
            "type": "object",
            "properties": {
                "input":            {**_S_STRING,
                                      "description": "Counts file path "
                                                     "(.h5ad/.h5/.mtx/.csv)."},
                "label":            {**_S_STRING},
                "min_genes":        {**_S_INTEGER, "default": 200},
                "min_cells":        {**_S_INTEGER, "default": 3},
                "max_mito":         {"type": "number", "default": 20.0},
                "mito_prefix":      {**_S_STRING, "default": "MT-"},
                "n_hvg":            {**_S_INTEGER, "default": 2000},
                "n_pcs":            {**_S_INTEGER, "default": 50},
                "n_neighbors":      {**_S_INTEGER, "default": 15},
                "resolution":       {"type": "number", "default": 1.0},
                "n_markers":        {**_S_INTEGER, "default": 25},
                "sample_col":       {**_S_STRING,
                                      "description": "Obs column to color "
                                                     "UMAP by (sample, batch)."},
                "highlight_genes":  {**_S_STRING,
                                      "description": "Comma-separated gene "
                                                     "symbols to overlay."},
                "skip_tsne":        {**_S_BOOLEAN, "default": False},
                "transpose":        {**_S_BOOLEAN, "default": False},
            },
            "required": ["input"],
        },
        cli=["sc-analyze", "pipeline"],
        flag_map={"input": "--input", "label": "--label",
                   "min_genes": "--min-genes", "min_cells": "--min-cells",
                   "max_mito": "--max-mito", "mito_prefix": "--mito-prefix",
                   "n_hvg": "--n-hvg", "n_pcs": "--n-pcs",
                   "n_neighbors": "--n-neighbors",
                   "resolution": "--resolution",
                   "n_markers": "--n-markers",
                   "sample_col": "--sample-col",
                   "highlight_genes": "--highlight-genes"},
        bool_flags={"skip_tsne", "transpose"},
    ),

    _T(
        "sc_qc",
        "Standalone QC pass: load counts, filter cells (min-genes), filter "
        "genes (min-cells), drop high-mito cells, write QC violins, save a "
        "cleaned processed.h5ad. Use when the user wants QC numbers before "
        "committing to a full pipeline run.",
        {
            "type": "object",
            "properties": {
                "input":       {**_S_STRING},
                "label":       {**_S_STRING},
                "min_genes":   {**_S_INTEGER, "default": 200},
                "min_cells":   {**_S_INTEGER, "default": 3},
                "max_mito":    {"type": "number", "default": 20.0},
                "mito_prefix": {**_S_STRING, "default": "MT-"},
                "transpose":   {**_S_BOOLEAN, "default": False},
            },
            "required": ["input"],
        },
        cli=["sc-analyze", "qc"],
        flag_map={"input": "--input", "label": "--label",
                   "min_genes": "--min-genes", "min_cells": "--min-cells",
                   "max_mito": "--max-mito",
                   "mito_prefix": "--mito-prefix"},
        bool_flags={"transpose"},
    ),

    _T(
        "sc_umap",
        "Run k-NN + UMAP on an already-PCA'd anndata, write a UMAP figure "
        "colored by clusters (if present) or QC. Skips PCA if X_pca exists.",
        {
            "type": "object",
            "properties": {
                "input":       {**_S_STRING},
                "label":       {**_S_STRING},
                "n_pcs":       {**_S_INTEGER, "default": 40},
                "n_neighbors": {**_S_INTEGER, "default": 15},
            },
            "required": ["input"],
        },
        cli=["sc-analyze", "umap"],
        flag_map={"input": "--input", "label": "--label",
                   "n_pcs": "--n-pcs", "n_neighbors": "--n-neighbors"},
    ),

    _T(
        "sc_cluster",
        "Leiden clustering on the k-NN graph (falls back to Louvain if "
        "leidenalg missing). Saves UMAP colored by leiden cluster.",
        {
            "type": "object",
            "properties": {
                "input":       {**_S_STRING},
                "label":       {**_S_STRING},
                "resolution":  {"type": "number", "default": 1.0},
                "n_pcs":       {**_S_INTEGER, "default": 40},
                "n_neighbors": {**_S_INTEGER, "default": 15},
            },
            "required": ["input"],
        },
        cli=["sc-analyze", "cluster"],
        flag_map={"input": "--input", "label": "--label",
                   "resolution": "--resolution",
                   "n_pcs": "--n-pcs", "n_neighbors": "--n-neighbors"},
    ),

    _T(
        "perturb_catalog_summary",
        "Landing-page stats for the Perturbation Catalogue: total "
        "datasets / experiments, top tissues / cell types / cell lines "
        "/ diseases / perturbation types. Use to orient before deeper "
        "queries.",
        {"type": "object", "properties": {}},
        cli=["perturb-catalog", "summary"],
    ),

    _T(
        "perturb_catalog_search",
        "Global gene/term search across the Perturbation Catalogue. "
        "Returns one row per perturbed gene with counts across MAVE / "
        "CRISPR-screen / Perturb-seq and the top GSEA terms. "
        "Headline tool when the user asks 'what perturbation data "
        "exists for gene X?'.",
        {
            "type": "object",
            "properties": {
                "query": {**_S_STRING,
                            "description": "Gene symbol or free-text query."},
                "size":  {**_S_INTEGER, "default": 25},
                "page":  {**_S_INTEGER, "default": 0},
                "facets": {**_S_STRING,
                            "description": "Comma-separated facet names."},
            },
            "required": ["query"],
        },
        cli=["perturb-catalog", "search"],
        flag_map={"query": "--query", "size": "--size", "page": "--page",
                   "facets": "--facets"},
    ),

    _T(
        "perturb_catalog_search_modality",
        "Modality-scoped search of the Perturbation Catalogue. Use "
        "`modality='mave'` for VAMP-seq / DMS data, `'crispr-screen'` "
        "for pooled CRISPR screens, `'perturb-seq'` for single-cell. "
        "Supports filters on gene name, position range, score "
        "name/value, tissue, cell line, disease, study year.",
        {
            "type": "object",
            "properties": {
                "modality": {**_S_STRING,
                              "description": "mave | crispr-screen | perturb-seq"},
                "query":    {**_S_STRING},
                "perturbation_gene_name": {**_S_STRING,
                              "description": "Filter to one perturbed gene."},
                "perturbation_position":  {**_S_STRING,
                              "description": "Position or range '100_300'."},
                "effect_score_name":  {**_S_STRING},
                "effect_score_value": {**_S_STRING,
                              "description": "Score range e.g. '0.5_1.0'."},
                "dataset_limit":  {**_S_INTEGER, "default": 25},
                "dataset_offset": {**_S_INTEGER, "default": 0},
                "rows_per_dataset_limit": {**_S_INTEGER, "default": 5},
                "sort": {**_S_STRING},
            },
            "required": ["modality"],
        },
        cli=["perturb-catalog", "search-modality"],
        flag_map={
            "modality": "--modality", "query": "--query",
            "perturbation_gene_name": "--perturbation-gene-name",
            "perturbation_position": "--perturbation-position",
            "effect_score_name": "--effect-score-name",
            "effect_score_value": "--effect-score-value",
            "dataset_limit": "--dataset-limit",
            "dataset_offset": "--dataset-offset",
            "rows_per_dataset_limit": "--rows-per-dataset-limit",
            "sort": "--sort",
        },
    ),

    _T(
        "perturb_catalog_dataset",
        "Fetch the full record for one Perturbation Catalogue dataset "
        "by id.",
        {
            "type": "object",
            "properties": {"dataset_id": {**_S_STRING}},
            "required": ["dataset_id"],
        },
        cli=["perturb-catalog", "dataset"],
        flag_map={"dataset_id": "--dataset-id"},
    ),

    _T(
        "perturb_catalog_dataset_rows",
        "Paginate the per-perturbation rows inside one dataset (variant- "
        "or gRNA-level effect scores).",
        {
            "type": "object",
            "properties": {
                "modality":   {**_S_STRING,
                                "description": "mave|crispr-screen|perturb-seq"},
                "dataset_id": {**_S_STRING},
                "limit":      {**_S_INTEGER, "default": 100},
                "offset":     {**_S_INTEGER, "default": 0},
            },
            "required": ["modality", "dataset_id"],
        },
        cli=["perturb-catalog", "dataset-rows"],
        flag_map={"modality": "--modality",
                   "dataset_id": "--dataset-id",
                   "limit": "--limit", "offset": "--offset"},
    ),

    _T(
        "perturb_catalog_gsea",
        "Perturb-seq GSEA hallmark/pathway enrichment table. Pass a "
        "gene `query` and/or a specific `dataset_id`.",
        {
            "type": "object",
            "properties": {
                "query":      {**_S_STRING},
                "dataset_id": {**_S_STRING},
                "page":       {**_S_INTEGER, "default": 0},
                "size":       {**_S_INTEGER, "default": 50},
            },
        },
        cli=["perturb-catalog", "gsea"],
        flag_map={"query": "--query", "dataset_id": "--dataset-id",
                   "page": "--page", "size": "--size"},
    ),

    _T(
        "perturb_catalog_pipeline",
        "End-to-end gene-centric pull from the Perturbation Catalogue: "
        "summary + global gene search + modality-scoped searches for "
        "MAVE / CRISPR-screen / Perturb-seq, writing a markdown report. "
        "The default tool when the user asks 'show me all perturbation "
        "data for gene X'.",
        {
            "type": "object",
            "properties": {
                "gene":         {**_S_STRING},
                "label":        {**_S_STRING},
                "dataset_limit": {**_S_INTEGER, "default": 10},
                "rows_per_dataset_limit": {**_S_INTEGER, "default": 3},
            },
            "required": ["gene"],
        },
        cli=["perturb-catalog", "pipeline"],
        flag_map={"gene": "--gene", "label": "--label",
                   "dataset_limit": "--dataset-limit",
                   "rows_per_dataset_limit": "--rows-per-dataset-limit"},
    ),

    _T(
        "multiseq_demultiplex",
        "Demultiplex a MULTI-seq / Cell Hashing tag-count matrix into "
        "singlet / multiplet / negative calls. Python port of "
        "deMULTIplex2 (Zhu et al. Nat Methods 2024). Fits a "
        "two-component negative-binomial mixture per tag via EM and "
        "writes classifications.csv + posteriors + residuals + "
        "diagnostic plots. The headline tool when the user asks 'who "
        "is each cell in a multiplexed run?'.",
        {
            "type": "object",
            "properties": {
                "input":         {**_S_STRING,
                                    "description": "Tag counts: .h5ad / 10x .h5 / "
                                                   ".csv / .tsv."},
                "label":         {**_S_STRING},
                "obsm_key":      {**_S_STRING,
                                    "description": ".h5ad obsm key with the "
                                                   "multiplexing matrix."},
                "init_cos_cut":  {"type": "number", "default": 0.5},
                "max_iter":      {**_S_INTEGER, "default": 10},
                "prob_cut":      {"type": "number", "default": 0.5},
                "residual_type": {**_S_STRING, "default": "rqr",
                                    "description": "rqr | pearson"},
                "seed":          {**_S_INTEGER, "default": 1},
                "transpose":     {**_S_BOOLEAN, "default": False},
                "skip_diagnostics": {**_S_BOOLEAN, "default": False},
            },
            "required": ["input"],
        },
        cli=["multiseq", "demultiplex"],
        flag_map={"input": "--input", "label": "--label",
                   "obsm_key": "--obsm-key",
                   "init_cos_cut": "--init-cos-cut",
                   "max_iter": "--max-iter",
                   "prob_cut": "--prob-cut",
                   "residual_type": "--residual-type",
                   "seed": "--seed"},
        bool_flags={"transpose", "skip_diagnostics"},
    ),

    _T(
        "multiseq_pipeline",
        "End-to-end MULTI-seq workflow: load tag counts, run EM "
        "demultiplexing, generate per-tag histograms + call heatmap + "
        "per-tag 4-panel diagnostics + markdown report. Accepts either "
        "a local file via `input` OR an IGVF Portal accession via "
        "`igvf_accession` (auto-pulls then runs). Optional ground-truth "
        "CSV produces an accuracy table for benchmarking.",
        {
            "type": "object",
            "properties": {
                "input":          {**_S_STRING,
                                    "description": "Local tag counts "
                                                   "(.h5ad / 10x .h5 / "
                                                   ".csv / .tsv)."},
                "igvf_accession": {**_S_STRING,
                                    "description": "IGVF Portal file "
                                                   "accession (auto-pull). "
                                                   "Mutually exclusive "
                                                   "with `input`."},
                "label":          {**_S_STRING},
                "obsm_key":       {**_S_STRING},
                "ground_truth":   {**_S_STRING},
                "init_cos_cut":   {"type": "number", "default": 0.5},
                "max_iter":       {**_S_INTEGER, "default": 10},
                "prob_cut":       {"type": "number", "default": 0.5},
                "residual_type":  {**_S_STRING, "default": "rqr"},
                "seed":           {**_S_INTEGER, "default": 1},
                "transpose":      {**_S_BOOLEAN, "default": False},
            },
        },
        cli=["multiseq", "pipeline"],
        flag_map={"input": "--input", "label": "--label",
                   "igvf_accession": "--igvf-accession",
                   "obsm_key": "--obsm-key",
                   "ground_truth": "--ground-truth",
                   "init_cos_cut": "--init-cos-cut",
                   "max_iter": "--max-iter",
                   "prob_cut": "--prob-cut",
                   "residual_type": "--residual-type",
                   "seed": "--seed"},
        bool_flags={"transpose"},
    ),

    _T(
        "multiseq_showcase",
        "★ ONE-COMMAND COMPREHENSIVE MULTI-SEQ DEMO ★. Generates synthetic "
        "tag matrix by default (2000 cells × 6 tags), OR pulls an IGVF "
        "Portal tag-count file via --igvf-accession, OR consumes a local "
        "--input matrix. Runs the deMULTIplex2 EM classifier "
        "(NB-GLM-based singlet/doublet/negative call), emits per-tag UMI "
        "histograms + call-group heatmap + per-tag posterior diagnostics "
        "(8+ plots), builds a publication composite figure, and writes a "
        "deep narrative report (including accuracy vs ground truth when "
        "available — typical synthetic run hits ≥85% singlet recovery). "
        "USE THIS TOOL for any 'show me MULTI-seq', 'demonstrate cell "
        "hashing', 'demultiplex MULTI-seq tags' question — do NOT call "
        "multiseq_pipeline + multiseq_histogram + multiseq_heatmap "
        "separately.",
        {
            "type": "object",
            "properties": {
                "input":           {**_S_STRING},
                "igvf_accession":  {**_S_STRING},
                "label":           {**_S_STRING},
                "n_cells":         {**_S_INTEGER, "default": 2000},
                "n_tags":          {**_S_INTEGER, "default": 6},
                "doublet_rate":    {"type": "number", "default": 0.08},
                "negative_rate":   {"type": "number", "default": 0.05},
                "seed":            {**_S_INTEGER, "default": 1},
            },
        },
        cli=["multiseq", "showcase"],
        flag_map={
            "input": "--input", "igvf_accession": "--igvf-accession",
            "label": "--label", "n_cells": "--n-cells", "n_tags": "--n-tags",
            "doublet_rate": "--doublet-rate", "negative_rate": "--negative-rate",
            "seed": "--seed",
        },
    ),

    _T(
        "multiseq_discover",
        "List MULTI-seq tag-count files available on the IGVF Portal. "
        "Defaults to content_type='cell hashing barcodes'. The output "
        "is a ranked manifest (largest files first) the user / agent "
        "can pick from to feed into `multiseq_pipeline`.",
        {
            "type": "object",
            "properties": {
                "content_type": {**_S_STRING,
                                  "default": "cell hashing barcodes"},
                "assay_title":  {**_S_STRING,
                                  "description": "Restrict to one "
                                                 "preferred_assay_titles."},
                "limit":        {**_S_INTEGER, "default": 20},
                "label":        {**_S_STRING},
            },
        },
        cli=["multiseq", "discover"],
        flag_map={"content_type": "--content-type",
                   "assay_title": "--assay-title",
                   "limit": "--limit", "label": "--label"},
    ),

    _T(
        "multiseq_pull_igvf",
        "Download a single MULTI-seq tag-count file from the IGVF "
        "Portal by accession into Data/MultiSeq/. Returns the local "
        "path which can be passed straight to `multiseq_pipeline`'s "
        "`input`.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING,
                                "description": "IGVF file accession, e.g. "
                                               "IGVFFI7138DMIL."},
            },
            "required": ["accession"],
        },
        cli=["multiseq", "pull-igvf"],
        flag_map={"accession": "--accession"},
    ),

    _T(
        "multiseq_simulate",
        "Generate a synthetic MULTI-seq cell × tag UMI matrix with "
        "ground truth. Useful for smoke-testing the demultiplexer or "
        "for the user to walk through the workflow without real data.",
        {
            "type": "object",
            "properties": {
                "n_cells":        {**_S_INTEGER, "default": 1000},
                "n_tags":         {**_S_INTEGER, "default": 4},
                "doublet_rate":   {"type": "number", "default": 0.05},
                "negative_rate":  {"type": "number", "default": 0.05},
                "pos_mean":       {"type": "number", "default": 1000.0},
                "bg_mean":        {"type": "number", "default": 20.0},
                "seed":           {**_S_INTEGER, "default": 7},
                "label":          {**_S_STRING, "default": "sim"},
            },
        },
        cli=["multiseq", "simulate"],
        flag_map={"n_cells": "--n-cells", "n_tags": "--n-tags",
                   "doublet_rate": "--doublet-rate",
                   "negative_rate": "--negative-rate",
                   "pos_mean": "--pos-mean", "bg_mean": "--bg-mean",
                   "seed": "--seed", "label": "--label"},
    ),

    _T(
        "sc_plot_embedding",
        "Re-render UMAP or t-SNE coloured by any combination of obs columns "
        "and/or gene symbols. Use when the user asks 'show me APOE on the "
        "UMAP we just made' or wants a sample/batch overlay.",
        {
            "type": "object",
            "properties": {
                "input":      {**_S_STRING},
                "label":      {**_S_STRING},
                "embedding":  {**_S_STRING, "default": "umap",
                                "description": "umap | tsne | both"},
                "color":      {**_S_STRING,
                                "description": "Comma-separated obs cols + genes."},
                "ncols":      {**_S_INTEGER, "default": 2},
            },
            "required": ["input", "color"],
        },
        cli=["sc-analyze", "plot-embedding"],
        flag_map={"input": "--input", "label": "--label",
                   "embedding": "--embedding", "color": "--color",
                   "ncols": "--ncols"},
    ),

    _T(
        "network_carnival",
        "CARNIVAL — given a signed perturbation set + signed measurement "
        "set (e.g. perturbed genes + DEG log2FCs from Perturb-seq), "
        "infer the minimum-cost upstream subnetwork in a signed PPI that "
        "explains the perturbations → measurements. Clean-room cvxpy "
        "MILP; no GPL dependencies. Appends selected edges to the "
        "warehouse with upstream='network:carnival:<label>'. Math "
        "reference: Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md.",
        {
            "type": "object",
            "properties": {
                "perturbations":  {**_S_STRING,
                                    "description": "CSV: gene,sign in {-1,+1}."},
                "measurements":   {**_S_STRING,
                                    "description": "CSV: gene,score "
                                                   "(signed log2FC)."},
                "pkn":            {**_S_STRING,
                                    "description": "SIF file path; if "
                                                   "omitted, build from "
                                                   "the proteomics KG."},
                "pkn_limit":      {**_S_INTEGER},
                "taxon":          {**_S_INTEGER, "default": 9606},
                "beta":           {"type": "number", "default": 0.2,
                                    "description": "L0 edge-sparsity vs "
                                                   "data-fit trade-off."},
                "lambda_v":       {"type": "number", "default": 0.0,
                                    "description": "L0 vertex sparsity."},
                "solver":         {**_S_STRING, "default": "SCIP"},
                "label":          {**_S_STRING, "default": "run"},
            },
            "required": ["perturbations", "measurements"],
        },
        cli=["network", "carnival"],
        flag_map={"perturbations": "--perturbations",
                   "measurements": "--measurements",
                   "pkn": "--pkn", "pkn_limit": "--pkn-limit",
                   "taxon": "--taxon", "beta": "--beta",
                   "lambda_v": "--lambda-v",
                   "solver": "--solver", "label": "--label"},
    ),

    _T(
        "network_viz",
        "★ PUBLICATION-GRADE NETWORK VISUALIZATION ★. Given a signed-SIF "
        "subnetwork file (from `network carnival`, `network steiner`, "
        "`network demo`, or any external CARNIVAL-style output), produces "
        "a force-directed graph (PNG + SVG) with signed edges (green = "
        "activation, red = inhibition), node coloring by role "
        "(perturbation / measurement / inferred up/down), node sizing by "
        "prize, plus a pathway-enrichment bar chart, degree distribution, "
        "edge-sign breakdown, and a 4-panel publication composite figure. "
        "Also emits an interactive vis.js HTML (with --html) and a per-"
        "node summary CSV. USE THIS TOOL whenever the user asks to "
        "'visualize', 'plot', 'render', or 'show' a network result.",
        {
            "type": "object",
            "properties": {
                "sif":      {**_S_STRING, "description":
                              "Path to signed-SIF input."},
                "prizes":   {**_S_STRING, "description":
                              "Optional per-node CSV (node, prize, sign, role)."},
                "pathways": {**_S_STRING, "description":
                              "Optional node->pathway CSV for enrichment."},
                "layout":   {**_S_STRING, "default": "spring",
                              "enum": ["spring", "kamada", "circular", "shell"]},
                "html":     {**_S_BOOLEAN, "default": False},
                "label":    {**_S_STRING},
                "title":    {**_S_STRING},
            },
            "required": ["sif"],
        },
        cli=["network", "viz"],
        flag_map={
            "sif": "--sif", "prizes": "--prizes", "pathways": "--pathways",
            "layout": "--layout", "html": "--html",
            "label": "--label", "title": "--title",
        },
    ),

    _T(
        "network_steiner",
        "Prize-collecting Steiner tree — given per-gene prizes (e.g. "
        "VAMP-seq abundance change, GWAS hit strength) and a PPI prior, "
        "find the connected subnetwork that maximises (prizes − costs). "
        "Clean-room cvxpy MILP. Appends selected edges to the warehouse.",
        {
            "type": "object",
            "properties": {
                "terminals":  {**_S_STRING,
                                "description": "CSV: gene,prize."},
                "pkn":        {**_S_STRING},
                "pkn_limit":  {**_S_INTEGER},
                "taxon":      {**_S_INTEGER, "default": 9606},
                "edge_cost":  {"type": "number", "default": 1.0},
                "solver":     {**_S_STRING, "default": "SCIP"},
                "label":      {**_S_STRING, "default": "run"},
            },
            "required": ["terminals"],
        },
        cli=["network", "steiner"],
        flag_map={"terminals": "--terminals", "pkn": "--pkn",
                   "pkn_limit": "--pkn-limit", "taxon": "--taxon",
                   "edge_cost": "--edge-cost", "solver": "--solver",
                   "label": "--label"},
    ),

    _T(
        "network_demo",
        "Self-test: synthetic EGFR → MYC cascade. Proves the CARNIVAL "
        "MILP recovers all 6 cascade edges and that the warehouse picks "
        "up the inferred subnetwork.",
        {
            "type": "object",
            "properties": {
                "beta":   {"type": "number", "default": 0.05},
                "solver": {**_S_STRING, "default": "SCIP"},
                "label":  {**_S_STRING, "default": "demo"},
            },
        },
        cli=["network", "demo"],
        flag_map={"beta": "--beta", "solver": "--solver", "label": "--label"},
    ),

    # ------------------------------------------------------------------
    # STARR-seq allelic analysis (clean-room rewrite of mpralm)
    # Ref: gaochengwen/STARR-seq-Data-Analysis (no LICENSE; clean-room)
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # Local IGVF Knowledge Graph mirror (Arango -> Parquet + DuckDB)
    # ------------------------------------------------------------------
    # ──────────────────────────────────────────────────────────────────
    # MaveDB → genomic coordinates mapping (clean-room reimpl of
    # ave-dcd/dcd_mapping, MIT)
    # ──────────────────────────────────────────────────────────────────
    _T(
        "mavedb_map_scoreset",
        "★ MAP MAVEDB SCORESET VARIANTS TO GENOMIC COORDINATES ★. Takes "
        "either a MaveDB URN or a gene symbol (looked up in the curated "
        "VAMP-seq catalog: PTEN/TPMT/VKOR/PRKN/CYP2C9/NUDT15), parses "
        "every variant's HGVSp, resolves chr/pos/ref/alt via the public "
        "Ensembl REST API, and emits TSV + VCF + summary JSON. Clean-room "
        "reimplementation of ave-dcd/dcd_mapping (MIT) without UTA / "
        "SeqRepo / BLAT — uses only Ensembl REST + an on-disk JSON cache. "
        "For amino acids with multiple alt codons, emits one row per "
        "candidate single-nt change with `candidate_idx` + `n_candidates`. "
        "USE THIS to make MAVE / VAMP-seq scores cross-referenceable "
        "with ClinVar, gnomAD, GWAS catalogues.",
        {
            "type": "object",
            "properties": {
                "urn":     {**_S_STRING, "description":
                            "MaveDB URN (e.g. urn:mavedb:00000013-a-1)."},
                "gene":    {**_S_STRING, "description":
                            "HGNC gene symbol (PTEN/TPMT/VKOR/PRKN/CYP2C9/NUDT15)."},
                "species": {**_S_STRING, "default": "human"},
                "label":   {**_S_STRING},
            },
        },
        cli=["mavedb", "map-scoreset"],
        flag_map={"urn": "--urn", "gene": "--gene", "species": "--species",
                   "label": "--label"},
    ),
    _T(
        "mavedb_showcase",
        "★ ONE-COMMAND MAVEDB DEMO ★. Downloads the canonical "
        "MaveDB scoreset for a gene, maps every variant to genomic "
        "coordinates, writes TSV + VCF + summary JSON + a 2-panel "
        "composite figure (per-protein-position variant coverage + "
        "mapping-outcome bar) + a narrative report. USE THIS for any "
        "'show me MAVE mapping for <gene>' demo question.",
        {
            "type": "object",
            "properties": {
                "gene":    {**_S_STRING, "default": "PTEN"},
                "species": {**_S_STRING, "default": "human"},
                "label":   {**_S_STRING},
            },
        },
        cli=["mavedb", "showcase"],
        flag_map={"gene": "--gene", "species": "--species", "label": "--label"},
    ),

    _T(
        "kg_mirror_inventory",
        "List Arango collections in the IGVF Catalog KG with per-collection "
        "document counts and on-disk byte sizes. Writes a CSV inventory.",
        {"type": "object", "properties": {}},
        cli=["kg-mirror", "inventory"],
    ),
    _T(
        "kg_mirror_pull",
        "Mirror a single Arango collection locally: stream via AQL cursor "
        "and persist as zstd-compressed Parquet shards under "
        "Data/Warehouse/KG/<collection>/. Resumable — re-run to continue.",
        {"type": "object", "properties": {
            "collection": {**_S_STRING},
            "batch_size": {**_S_INTEGER, "default": 5000},
            "max_rows":   {**_S_INTEGER},
            "restart":    {**_S_BOOLEAN, "default": False},
        }, "required": ["collection"]},
        cli=["kg-mirror", "pull"],
        flag_map={"collection": "--collection", "batch_size": "--batch-size",
                   "max_rows": "--max-rows", "restart": "--restart"},
    ),
    _T(
        "kg_mirror_pull_all",
        "Mirror every Arango collection except the skip list (default skip: "
        "variants and variants_variants — together ~1.5 TB). Small "
        "collections first, then medium, then large. Resumable.",
        {"type": "object", "properties": {
            "skip":      {**_S_STRING, "default": "variants,variants_variants"},
            "only":      {**_S_STRING},
            "include_giants": {**_S_BOOLEAN, "default": False},
            "max_collection_bytes": {**_S_INTEGER},
            "batch_size": {**_S_INTEGER, "default": 5000},
            "max_rows":   {**_S_INTEGER},
            "restart":    {**_S_BOOLEAN, "default": False},
        }},
        cli=["kg-mirror", "pull-all"],
        flag_map={"skip": "--skip", "only": "--only",
                   "include_giants": "--include-giants",
                   "max_collection_bytes": "--max-collection-bytes",
                   "batch_size": "--batch-size", "max_rows": "--max-rows",
                   "restart": "--restart"},
    ),
    _T(
        "kg_mirror_register",
        "Register the on-disk Parquet shards as DuckDB views in "
        "Data/Warehouse/igvf_kg_mirror.duckdb (one view per collection, "
        "named kg_<collection>).",
        {"type": "object", "properties": {}},
        cli=["kg-mirror", "register"],
    ),
    _T(
        "kg_mirror_verify",
        "Print row counts for every kg_* view in the local DuckDB warehouse.",
        {"type": "object", "properties": {}},
        cli=["kg-mirror", "verify"],
    ),

    _T(
        "starr_pull_portal",
        "★ Discover IGVF Portal STARR-seq MeasurementSets. Writes a TSV "
        "manifest with accession, assay_titles, n_files, status.",
        {"type": "object", "properties": {
            "limit": {**_S_INTEGER, "default": 50},
            "label": {**_S_STRING},
        }},
        cli=["starrseq", "pull-portal"],
        flag_map={"limit": "--limit", "label": "--label"},
    ),
    _T(
        "starr_qc",
        "STARR-seq counts QC: TPM-style scaling, low-expression filter, "
        "RLE matrix, per-sample Spearman D-statistic outlier flag (top "
        "5%). Mirrors count_qc.R from gaochengwen/STARR-seq-Data-Analysis.",
        {"type": "object", "properties": {
            "input": {**_S_STRING, "description": "Counts TSV/CSV."},
            "label": {**_S_STRING},
        }, "required": ["input"]},
        cli=["starrseq", "qc"],
        flag_map={"input": "--input", "label": "--label"},
    ),
    _T(
        "starr_aggregate",
        "Collapse barcode-level STARR-seq counts to a per-(SNP, Allele) "
        "wide table (long form -> wide).",
        {"type": "object", "properties": {
            "input": {**_S_STRING},
            "label": {**_S_STRING},
        }, "required": ["input"]},
        cli=["starrseq", "aggregate"],
        flag_map={"input": "--input", "label": "--label"},
    ),
    _T(
        "starr_activity",
        "Per-fragment STARR-seq log activity = log(RNA TPM / DNA TPM) "
        "across paired replicates, plus per-SNP allelic skew "
        "(log activity ALT - REF).",
        {"type": "object", "properties": {
            "input": {**_S_STRING},
            "label": {**_S_STRING},
        }, "required": ["input"]},
        cli=["starrseq", "activity"],
        flag_map={"input": "--input", "label": "--label"},
    ),
    _T(
        "starr_allelic_test",
        "STARR-seq allelic significance test: per-(SNP, Allele) OLS on "
        "log2((RNA+0.5)/(DNA+0.5)) per replicate, eBayes variance "
        "moderation (Smyth 2004 trigamma inversion) with graceful "
        "fallback when no variance excess is detected, BH-FDR. "
        "Clean-room rewrite of Bioconductor mpra::mpralm.",
        {"type": "object", "properties": {
            "input": {**_S_STRING},
            "label": {**_S_STRING},
        }, "required": ["input"]},
        cli=["starrseq", "allelic-test"],
        flag_map={"input": "--input", "label": "--label"},
    ),

    # ------------------------------------------------------------------
    # SHARE-seq joint scATAC + scRNA QC
    # Ref: broadinstitute/epi-SHARE-seq-pipeline (MIT)
    #      Ma et al. Cell 2020 doi:10.1016/j.cell.2020.09.056
    # ------------------------------------------------------------------
    _T(
        "share_pull_portal",
        "★ Discover IGVF Portal SHARE-seq AnalysisSets and MeasurementSets "
        "(`preferred_assay_titles=SHARE-seq`). Writes a TSV manifest "
        "covering both processed AnalysisSets (h5ad + fragments BED) and "
        "raw MeasurementSets (FASTQ + seqspec).",
        {"type": "object", "properties": {
            "limit": {**_S_INTEGER, "default": 50},
            "label": {**_S_STRING},
        }},
        cli=["share", "pull-portal"],
        flag_map={"limit": "--limit", "label": "--label"},
    ),
    _T(
        "share_demultiplex_bcs",
        "SHARE-seq round-1/2/3 24-mer barcode demultiplex on a gz FASTQ. "
        "Exact match + 1-Hamming-mismatch lookup + optional +/-1 bp shift "
        "correction. Pure stdlib; no pysam/dnaio. Mirrors correct_fastq.py "
        "from broadinstitute/epi-SHARE-seq-pipeline.",
        {"type": "object", "properties": {
            "fastq": {**_S_STRING},
            "whitelist": {**_S_STRING, "description":
                          "24-mer whitelist (one barcode per line)."},
            "out": {**_S_STRING},
            "r1_offset": {**_S_INTEGER, "default": 14},
            "r2_offset": {**_S_INTEGER, "default": 52},
            "r3_offset": {**_S_INTEGER, "default": 90},
            "shift_correct": {**_S_BOOLEAN, "default": False},
            "max_reads": {**_S_INTEGER, "default": 0},
            "label": {**_S_STRING},
        }, "required": ["fastq", "whitelist", "out"]},
        cli=["share", "demultiplex-bcs"],
        flag_map={"fastq": "--fastq", "whitelist": "--whitelist",
                  "out": "--out", "r1_offset": "--r1-offset",
                  "r2_offset": "--r2-offset", "r3_offset": "--r3-offset",
                  "shift_correct": "--shift-correct",
                  "max_reads": "--max-reads", "label": "--label"},
    ),
    _T(
        "share_fragment_qc",
        "Per-barcode ATAC QC from a SHARE-seq fragments BED: total "
        "fragments, reads in TSS +/-2kb window, reads in flanking 100bp "
        "regions, TSS enrichment (Ma 2020 formula with 0.2 floor), reads "
        "in peaks, FRIP. Mirrors qc_atac_compute_tss_enrichment.py + "
        "qc_atac_compute_reads_in_peaks.py.",
        {"type": "object", "properties": {
            "fragments": {**_S_STRING},
            "tss_bed": {**_S_STRING, "description": "Optional TSS BED."},
            "peaks_bed": {**_S_STRING, "description": "Optional peaks BED."},
            "label": {**_S_STRING},
        }, "required": ["fragments"]},
        cli=["share", "fragment-qc"],
        flag_map={"fragments": "--fragments", "tss_bed": "--tss-bed",
                  "peaks_bed": "--peaks-bed", "label": "--label"},
    ),
    _T(
        "share_rna_qc",
        "Per-barcode RNA QC from a SHARE-seq h5ad sparse gene-count "
        "matrix: total UMIs, expressed genes, percent mitochondrial "
        "(genes auto-detected by MT-/mt- prefix).",
        {"type": "object", "properties": {
            "h5ad": {**_S_STRING},
            "label": {**_S_STRING},
        }, "required": ["h5ad"]},
        cli=["share", "rna-qc"],
        flag_map={"h5ad": "--h5ad", "label": "--label"},
    ),
    _T(
        "share_joint_qc",
        "Merge SHARE-seq RNA + ATAC per-barcode tables and apply the "
        "Ma 2020 joint cell-calling thresholds (UMIs>=min_umis AND "
        "genes>=min_genes AND TSS>=min_tss AND fragments>=min_frags). "
        "Mirrors joint_cell_plotting.py.",
        {"type": "object", "properties": {
            "rna_qc": {**_S_STRING},
            "atac_qc": {**_S_STRING},
            "min_umis": {**_S_INTEGER, "default": 100},
            "min_genes": {**_S_INTEGER, "default": 200},
            "min_tss": {"type": "number", "default": 4.0},
            "min_frags": {**_S_INTEGER, "default": 100},
            "label": {**_S_STRING},
        }, "required": ["rna_qc", "atac_qc"]},
        cli=["share", "joint-qc"],
        flag_map={"rna_qc": "--rna-qc", "atac_qc": "--atac-qc",
                  "min_umis": "--min-umis", "min_genes": "--min-genes",
                  "min_tss": "--min-tss", "min_frags": "--min-frags",
                  "label": "--label"},
    ),
    _T(
        "share_multiplet_detect",
        "Pairwise Jaccard multiplet detection on a SHARE-seq fragments "
        "BED. Builds a per-barcode set of (chrom, start) coordinates, "
        "samples a null Jaccard distribution from 10000 random pairs, "
        "and flags any pair whose Jaccard exceeds the 99-th percentile "
        "of null. Mirrors detect_multiplets.py.",
        {"type": "object", "properties": {
            "fragments": {**_S_STRING},
            "min_fragments": {**_S_INTEGER, "default": 1000},
            "max_pairs": {**_S_INTEGER, "default": 200000},
            "null_samples": {**_S_INTEGER, "default": 10000},
            "seed": {**_S_INTEGER, "default": 7},
            "label": {**_S_STRING},
        }, "required": ["fragments"]},
        cli=["share", "multiplet-detect"],
        flag_map={"fragments": "--fragments",
                  "min_fragments": "--min-fragments",
                  "max_pairs": "--max-pairs",
                  "null_samples": "--null-samples",
                  "seed": "--seed", "label": "--label"},
    ),

    # ------------------------------------------------------------------
    # CRISPRi Flow-FISH screen analysis
    # Ref: EngreitzLab/CRISPRi-FlowFISH-pipeline (MIT)
    #      Fulco 2019 Nat Genet doi:10.1038/s41588-019-0538-0
    #      Nasser 2021 Nature  doi:10.1038/s41586-021-03446-x
    # ------------------------------------------------------------------
    _T(
        "flowfish_pull_portal",
        "Discover IGVF Portal CRISPRi-FlowFISH MeasurementSets.",
        {"type": "object", "properties": {
            "limit": {**_S_INTEGER, "default": 50},
            "label": {**_S_STRING},
        }},
        cli=["flowfish", "pull-portal"],
        flag_map={"limit": "--limit", "label": "--label"},
    ),
    _T(
        "flowfish_estimate_effects",
        "Per-guide log-normal MLE on a guide x FACS-bin counts matrix. "
        "Fits (logMean, logSD) by minimizing the bin-multinomial NLL with "
        "L-BFGS-B and a trailing 'outside' bin EM imputation, then "
        "carries metadata through to the raw_effects.tsv output. Clean-"
        "room implementation of estimate_effect_sizes.R from "
        "EngreitzLab/CRISPRi-FlowFISH-pipeline.",
        {"type": "object", "properties": {
            "counts": {**_S_STRING},
            "sortparams": {**_S_STRING, "description":
                            "TSV with Bin, LowBound, HighBound."},
            "label": {**_S_STRING},
        }, "required": ["counts", "sortparams"]},
        cli=["flowfish", "estimate-effects"],
        flag_map={"counts": "--counts", "sortparams": "--sortparams",
                  "label": "--label"},
    ),
    _T(
        "flowfish_real_space",
        "Convert per-guide (logMean, logSD) to log-normal mean expression "
        "(mleAvg = exp(mu*ln10 + (sigma*ln10)^2/2)), divide by the "
        "median negative-control mleAvg to get a fold-change vs null, "
        "clamp at +/-clamp, and rescale so null centers at 1. Clean-room "
        "implementation of convert_to_real_space.py.",
        {"type": "object", "properties": {
            "input": {**_S_STRING},
            "target_col": {**_S_STRING, "default": "target"},
            "negative_label": {**_S_STRING, "default": "negative_control"},
            "clamp": {"type": "number", "default": 5.0},
            "label": {**_S_STRING},
        }, "required": ["input"]},
        cli=["flowfish", "real-space"],
        flag_map={"input": "--input", "target_col": "--target-col",
                  "negative_label": "--negative-label",
                  "clamp": "--clamp", "label": "--label"},
    ),
    _T(
        "flowfish_score_elements",
        "Per-element collapse + significance: two tests vs negative-"
        "control distribution (Mann-Whitney U + Welch t-test) with BH-"
        "FDR. Flags Significant = (FDR<thr) AND (mean<1) AND (n>=min) "
        "and Regulated = Significant AND |1-mean|>=min_effect. Clean-room "
        "implementation of ScoreEnhancers.py (Fulco 2019).",
        {"type": "object", "properties": {
            "effects": {**_S_STRING},
            "target_col": {**_S_STRING, "default": "target"},
            "element_col": {**_S_STRING, "default": "ElementName"},
            "negative_label": {**_S_STRING, "default": "negative_control"},
            "min_guides": {**_S_INTEGER, "default": 5},
            "min_negative": {**_S_INTEGER, "default": 10},
            "fdr": {"type": "number", "default": 0.05},
            "min_effect": {"type": "number", "default": 0.10},
            "label": {**_S_STRING},
        }, "required": ["effects"]},
        cli=["flowfish", "score-elements"],
        flag_map={"effects": "--effects", "target_col": "--target-col",
                  "element_col": "--element-col",
                  "negative_label": "--negative-label",
                  "min_guides": "--min-guides",
                  "min_negative": "--min-negative",
                  "fdr": "--fdr", "min_effect": "--min-effect",
                  "label": "--label"},
    ),
    _T(
        "flowfish_simulate",
        "Generate a synthetic guide x bin counts table for smoke testing "
        "the Flow-FISH skill: 100 negative controls + N elements each "
        "with G guides, fraction-knockdown elements pushed to a lower "
        "log-normal mean.",
        {"type": "object", "properties": {
            "out_dir": {**_S_STRING},
            "n_elements": {**_S_INTEGER, "default": 40},
            "guides_per_element": {**_S_INTEGER, "default": 8},
            "knockdown_frac": {"type": "number", "default": 0.30},
            "cells_per_guide": {**_S_INTEGER, "default": 800},
            "seed": {**_S_INTEGER, "default": 7},
        }, "required": ["out_dir"]},
        cli=["flowfish", "simulate"],
        flag_map={"out_dir": "--out-dir", "n_elements": "--n-elements",
                  "guides_per_element": "--guides-per-element",
                  "knockdown_frac": "--knockdown-frac",
                  "cells_per_guide": "--cells-per-guide", "seed": "--seed"},
    ),

    # ──────────────────────────────────────────────────────────────────
    # GO + Pathway enrichment (validation layer)
    # ──────────────────────────────────────────────────────────────────
    _T(
        "enrich_ora",
        "★ OVER-REPRESENTATION ENRICHMENT (ORA) FOR A GENE LIST ★. Runs "
        "hypergeometric / Fisher enrichment of a discrete gene list "
        "against GO_BP, GO_MF, GO_CC, Reactome 2022, KEGG 2021 Human, "
        "WikiPathways 2024 Human, and MSigDB Hallmark 2020 via gseapy "
        "→ Enrichr. Emits a per-library TSV, a composite figure "
        "(top-K bar charts + bubble overview) and a JSON summary with "
        "top-10 terms. USE THIS to validate a DEG list, CRISPR-screen "
        "hits, or the gene-side of an enhancer-gene linkage by asking "
        "which biological processes / pathways are over-represented.",
        {
            "type": "object",
            "properties": {
                "genes":      {**_S_STRING, "description":
                                "Path to a gene list (one per line, or "
                                "CSV/TSV with a 'gene' column)."},
                "background": {**_S_STRING, "description":
                                "Optional background gene list."},
                "libs":       {**_S_STRING, "default": "all",
                                "description":
                                "'all' / 'go' / 'pathways' or comma-list of "
                                "friendly lib names (GO_BP, GO_MF, GO_CC, "
                                "Reactome, KEGG, WikiPathways, MSigDB_Hallmark) "
                                "or raw Enrichr lib ids."},
                "organism":   {**_S_STRING, "default": "human"},
                "label":      {**_S_STRING},
                "top_k":      {**_S_INTEGER, "default": 8},
            },
            "required": ["genes"],
        },
        cli=["enrich", "ora"],
        flag_map={"genes": "--genes", "background": "--background",
                   "libs": "--libs", "organism": "--organism",
                   "label": "--label", "top_k": "--top-k"},
    ),
    _T(
        "enrich_gsea",
        "★ PRERANKED GSEA FOR A RANKED GENE-SCORE TABLE ★. Subramanian-"
        "style enrichment (no arbitrary cutoff) of a ranked gene list "
        "against GO + Reactome + KEGG + WikiPathways + MSigDB Hallmark "
        "via gseapy.prerank. Input: TSV/CSV with a 'gene' column and a "
        "'score' (or 'stat' / 'log2fc' / 'rank' / 't' / 'wald') column. "
        "Emits NES + FDR table, composite NES bubble figure, JSON "
        "summary. USE THIS when you have continuous statistics for "
        "every gene (e.g. limma t-stats, DESeq2 Wald) instead of a "
        "discrete hit list.",
        {
            "type": "object",
            "properties": {
                "ranked":       {**_S_STRING, "description":
                                  "Path to a TSV/CSV with 'gene' + 'score' "
                                  "columns."},
                "libs":         {**_S_STRING, "default": "all"},
                "organism":     {**_S_STRING, "default": "human"},
                "min_size":     {**_S_INTEGER, "default": 10},
                "max_size":     {**_S_INTEGER, "default": 1000},
                "permutations": {**_S_INTEGER, "default": 1000},
                "label":        {**_S_STRING},
                "top_k":        {**_S_INTEGER, "default": 8},
            },
            "required": ["ranked"],
        },
        cli=["enrich", "gsea"],
        flag_map={"ranked": "--ranked", "libs": "--libs",
                   "organism": "--organism", "min_size": "--min-size",
                   "max_size": "--max-size",
                   "permutations": "--permutations",
                   "label": "--label", "top_k": "--top-k"},
    ),
    _T(
        "enrich_go",
        "★ GENE-ONTOLOGY ORA (BP + MF + CC) ★. Convenience wrapper around "
        "enrich_ora restricted to the three Gene Ontology branches "
        "(GO_Biological_Process_2023, GO_Molecular_Function_2023, "
        "GO_Cellular_Component_2023). USE THIS when the validation "
        "question is specifically 'which GO terms are over-represented' "
        "rather than 'which canonical pathways'.",
        {
            "type": "object",
            "properties": {
                "genes":      {**_S_STRING},
                "background": {**_S_STRING},
                "organism":   {**_S_STRING, "default": "human"},
                "label":      {**_S_STRING},
                "top_k":      {**_S_INTEGER, "default": 8},
            },
            "required": ["genes"],
        },
        cli=["enrich", "go"],
        flag_map={"genes": "--genes", "background": "--background",
                   "organism": "--organism", "label": "--label",
                   "top_k": "--top-k"},
    ),
    _T(
        "enrich_pathways",
        "★ CANONICAL-PATHWAY ORA (Reactome + KEGG + WikiPathways + "
        "MSigDB Hallmark) ★. Convenience wrapper around enrich_ora "
        "restricted to the four canonical pathway databases. USE THIS "
        "when the validation question is 'which signalling / metabolic "
        "pathways are enriched' rather than 'which GO terms'.",
        {
            "type": "object",
            "properties": {
                "genes":      {**_S_STRING},
                "background": {**_S_STRING},
                "organism":   {**_S_STRING, "default": "human"},
                "label":      {**_S_STRING},
                "top_k":      {**_S_INTEGER, "default": 8},
            },
            "required": ["genes"],
        },
        cli=["enrich", "pathways"],
        flag_map={"genes": "--genes", "background": "--background",
                   "organism": "--organism", "label": "--label",
                   "top_k": "--top-k"},
    ),
    # ──────────────────────────────────────────────────────────────────
    # IGVF Portal canonical-query layer (clean-room reimpl of
    # IGVF-DACC/igvf-portal-mcp, MIT)
    # ──────────────────────────────────────────────────────────────────
    _T(
        "portal_search",
        "★ IGVF PORTAL FACETED SEARCH (canonical DACC pattern) ★. Searches "
        "the IGVF Portal `/search/` endpoint with a typed ItemType filter "
        "(MeasurementSet / AnalysisSet / SequenceFile / HumanDonor / etc.) "
        "and an optional free-text query, returning matching `@graph` "
        "items + total count + saved JSON. Supports the DACC field-filter "
        "DSL: dotted embedded fields (`lab.@id=/labs/x`), negation "
        "(`field!=val`), range ops (`gte:`/`lte:`/`gt:`/`lt:`), and list "
        "values (`field=a,b,c`). Clauses joined by `;`. USE THIS for any "
        "'find me IGVF Portal items where ...' question. Works "
        "anonymously over `api.data.igvf.org` for released items; HTTP "
        "Basic auth via IGVF_ACCESS_KEY/SECRET unlocks restricted ones.",
        {
            "type": "object",
            "properties": {
                "type":          {**_S_STRING, "description":
                                  "Comma-list of CamelCase ItemTypes "
                                  "(MeasurementSet, AnalysisSet, SequenceFile, "
                                  "HumanDonor, Tissue, Gene, etc.)."},
                "query":         {**_S_STRING, "description": "Free-text search."},
                "field_filters": {**_S_STRING, "description":
                                   "DSL: 'lab.@id=/labs/x;file_format=bam,bed'."},
                "limit":         {**_S_STRING, "default": "25",
                                   "description": "Integer or 'all'."},
                "sort":          {**_S_STRING},
                "frame":         {**_S_STRING},
                "label":         {**_S_STRING},
            },
        },
        cli=["portal", "search"],
        flag_map={"type": "--type", "query": "--query",
                   "field_filters": "--field-filters",
                   "limit": "--limit", "sort": "--sort",
                   "frame": "--frame", "label": "--label"},
    ),
    _T(
        "portal_get",
        "★ FETCH ONE IGVF PORTAL ITEM ★. Resolves a single item by "
        "@id ('/measurement-sets/IGVFDS3909HJKS/'), accession "
        "('IGVFDS3909HJKS'), or UUID, returning the full embedded JSON "
        "+ @type breakdown + saved file. USE THIS to inspect a specific "
        "MeasurementSet / AnalysisSet / Donor / etc. that you know the "
        "id of.",
        {
            "type": "object",
            "properties": {
                "resource_id": {**_S_STRING, "description":
                                  "@id, accession (IGVFFI...), or UUID."},
            },
            "required": ["resource_id"],
        },
        cli=["portal", "get"],
        flag_map={"resource_id": ""},   # positional
    ),
    _T(
        "portal_schema",
        "★ IGVF JSON SCHEMA INTROSPECTION ★. Fetches the canonical JSON "
        "schema (`/profiles/<Type>.json`) for any CamelCase ItemType "
        "(MeasurementSet, SequenceFile, HumanDonor, etc.) and reports "
        "title + description + all 60+ documented properties. USE THIS "
        "to discover what fields exist on a type before constructing a "
        "search / report query.",
        {
            "type": "object",
            "properties": {
                "item_type": {**_S_STRING, "description":
                                "CamelCase ItemType (e.g. MeasurementSet)."},
            },
            "required": ["item_type"],
        },
        cli=["portal", "schema"],
        flag_map={"item_type": ""},
    ),
    _T(
        "portal_facets",
        "★ IGVF PORTAL FACETS-ONLY CALL ★. Calls `/search/` with "
        "`limit=0` to retrieve only the `facets[]` aggregation block — "
        "no items, just per-category value counts (assay titles, lab "
        "ids, taxa, sample terms, perturbation modality, file formats, "
        "etc.). USE THIS to summarise what IGVF has across a slice of "
        "data (e.g. 'how many of each assay class for human putamen "
        "tissue?'). Set `field` to the name of ONE facet (e.g. "
        "'content_type') to dump EVERY value of that facet with counts "
        "instead of the top-5 preview — the reliable way to discover the "
        "exact filter term for a rare value in a single call, rather than "
        "guessing field_filters and getting 404s.",
        {
            "type": "object",
            "properties": {
                "type":          {**_S_STRING},
                "query":         {**_S_STRING},
                "field_filters": {**_S_STRING},
                "field":         {**_S_STRING, "description":
                    "List every value of this one facet field (e.g. "
                    "'content_type', 'assay_titles'); lists available "
                    "fields if the name does not match."},
                "label":         {**_S_STRING},
            },
        },
        cli=["portal", "facets"],
        flag_map={"type": "--type", "query": "--query",
                   "field_filters": "--field-filters", "field": "--field",
                   "label": "--label"},
    ),
    _T(
        "portal_report",
        "★ IGVF PORTAL TSV REPORT EXPORT ★. Streams the canonical "
        "`/report.tsv` export for an ItemType (paginated server-side) "
        "to disk. Returns row + column count + first-8-column preview "
        "+ saved TSV path. USE THIS when you need a spreadsheet-shaped "
        "dump of every matching item (lab inventories, file manifests, "
        "donor cohorts) rather than the per-item JSON.",
        {
            "type": "object",
            "properties": {
                "type":          {**_S_STRING},
                "query":         {**_S_STRING},
                "field_filters": {**_S_STRING},
                "limit":         {**_S_STRING, "default": "all"},
                "label":         {**_S_STRING},
            },
            "required": ["type"],
        },
        cli=["portal", "report"],
        flag_map={"type": "--type", "query": "--query",
                   "field_filters": "--field-filters",
                   "limit": "--limit", "label": "--label"},
    ),
    _T(
        "portal_batch_download",
        "★ IGVF PORTAL BATCH-DOWNLOAD MANIFEST ★. Hits "
        "`/batch-download/` for a FileSet type (MeasurementSet / "
        "AnalysisSet / etc.) and returns a manifest of pre-signed S3 "
        "URLs for every File hanging off the selected sets. Optional "
        "`fetch=true` actually pulls each file to disk. USE THIS to "
        "grab raw data files for a slice of IGVF in one round-trip.",
        {
            "type": "object",
            "properties": {
                "type":          {**_S_STRING},
                "query":         {**_S_STRING},
                "field_filters": {**_S_STRING},
                "limit":         {**_S_STRING, "default": "all"},
                "fetch":         {**_S_BOOLEAN, "default": False,
                                   "description":
                                   "Also download each file in the manifest."},
                "label":         {**_S_STRING},
            },
            "required": ["type"],
        },
        cli=["portal", "batch-download"],
        flag_map={"type": "--type", "query": "--query",
                   "field_filters": "--field-filters",
                   "limit": "--limit", "fetch": "--fetch",
                   "label": "--label"},
    ),
    _T(
        "portal_endpoint_params",
        "★ IGVF PORTAL ENDPOINT-PARAM INTROSPECTION (the DACC UX trick) ★. "
        "For any portal collection (`measurement-sets`, `analysis-sets`, "
        "`sequence-files`, ...), returns the documented filter set as a "
        "table mapping the Python-friendly snake_case `agent_param` "
        "(file_set_id, lab_title, samples_taxa, ...) to the dotted "
        "`search_field` the search index understands (file_set.@id, "
        "lab.title, samples.taxa, ...). Lets an LLM pivot between "
        "tooled and raw search without re-learning fields.",
        {
            "type": "object",
            "properties": {
                "collection": {**_S_STRING, "description":
                                "Snake-case collection (e.g. measurement-sets)."},
            },
            "required": ["collection"],
        },
        cli=["portal", "endpoint-params"],
        flag_map={"collection": ""},
    ),
    # ──────────────────────────────────────────────────────────────────
    # IGVF Catalog (Knowledge Graph) canonical-query layer (clean-room
    # reimpl of IGVF-DACC/igvf-catalog-mcp, MIT)
    # ──────────────────────────────────────────────────────────────────
    _T(
        "catalog_get_entity",
        "★ IGVF CATALOG UNIVERSAL ENTITY LOOKUP ★. Resolves any IGVF "
        "Catalog ID — gene symbol (APOE), ENSG, HGNC:n, ENTREZ:n, "
        "rsID (rs429358), SPDI (NC_000019.10:...), HGVS, CA-id, ENSP, "
        "UniProt (P02649), MONDO:n / EFO:n / GO:n / HPO / DOID / UBERON "
        "/ CL / CHEBI / OBA, drugbank (DB...), CHEMBL..., CPX-n "
        "(Complex Portal), R-HSA-n (Reactome), GCST... (GWAS Catalog) "
        "— and returns the full node JSON. Auto-detects entity type "
        "from the ID format; --hint overrides.",
        {
            "type": "object",
            "properties": {
                "id":    {**_S_STRING, "description":
                          "Any IGVF Catalog ID (auto-detected)."},
                "hint":  {**_S_STRING, "description":
                          "Override entity-type detection."},
                "limit": {**_S_INTEGER, "default": 1},
            },
            "required": ["id"],
        },
        cli=["catalog", "get-entity"],
        flag_map={"id": "", "hint": "--hint", "limit": "--limit"},
    ),
    _T(
        "catalog_search_region",
        "★ IGVF CATALOG REGION FAN-OUT SEARCH ★. Parallel query of "
        "genes + variants + genomic-elements within a region. Accepts "
        "'chr19:44,907,000-44,910,000', '19:44.9M-44.92M', "
        "'chr1:1K-2K'. Returns three lists (one per type) with "
        "per-type pagination metadata. USE THIS to enumerate "
        "everything the Catalog knows about a locus.",
        {
            "type": "object",
            "properties": {
                "region":   {**_S_STRING, "description":
                              "chr1:1000-2000, 1:1K-2K, 19:44.9M-44.92M."},
                "include":  {**_S_STRING, "default":
                              "genes,variants,genomic-elements"},
                "organism": {**_S_STRING, "default": "Homo sapiens"},
                "limit":    {**_S_INTEGER, "default": 25},
                "page":     {**_S_INTEGER, "default": 0},
            },
            "required": ["region"],
        },
        cli=["catalog", "search-region"],
        flag_map={"region": "", "include": "--include",
                   "organism": "--organism", "limit": "--limit",
                   "page": "--page"},
    ),
    _T(
        "catalog_find_associations",
        "★ IGVF CATALOG EDGE QUERY BY SEMANTIC RELATIONSHIP ★. Walks "
        "every edge in a semantic category (genetic / regulatory / "
        "physical / functional / pharmacological / ld / coding / "
        "transcription / all) for an entity, aggregating hits across "
        "multiple edge endpoints in one call. Accepts the catalog "
        "filter DSL (label=eqtl;method=GTEx,FANTOM5;p_value=lte:5e-8) "
        "with automatic p_value → log10pvalue conversion. IGVF 'method' "
        "values worth filtering on (all with source=IGVF): 'Perturb-seq' "
        "and 'CRISPR screen' on element->gene edges, 'STARR-seq' and "
        "'BlueSTARR' on variant->biosample, 'Variant-EFFECTS' on "
        "variant->gene, 'cV2F' on variant->phenotype, 'SEMVAR' on "
        "protein->variant. Note this tool does NOT cover the "
        "gene-regulatory-network collection — use grn_network for dEx "
        "regulator/target edges.",
        {
            "type": "object",
            "properties": {
                "entity_id":     {**_S_STRING, "description":
                                  "Any IGVF Catalog ID (auto-detected)."},
                "relationship":  {**_S_STRING, "description":
                                  "genetic / regulatory / physical / "
                                  "functional / pharmacological / ld / "
                                  "coding / transcription / all"},
                "filters":       {**_S_STRING, "description":
                                  "DSL e.g. 'label=eqtl;p_value=lte:5e-8'."},
                "limit":         {**_S_INTEGER, "default": 25},
                "page":          {**_S_INTEGER, "default": 0},
                "verbose":       {**_S_BOOLEAN, "default": False},
            },
            "required": ["entity_id", "relationship"],
        },
        cli=["catalog", "find-associations"],
        flag_map={"entity_id": "", "relationship": "--relationship",
                   "filters": "--filters", "limit": "--limit",
                   "page": "--page", "verbose": "--verbose"},
    ),
    _T(
        "catalog_find_ld",
        "★ IGVF CATALOG LD PROXIES ★. Dedicated query of `/api/variants/"
        "variant-ld` for an index variant, with r² / D' / ancestry "
        "filters and a strong/moderate/weak/negligible bucket summary "
        "plus per-ancestry breakdown.",
        {
            "type": "object",
            "properties": {
                "variant_id":         {**_S_STRING, "description":
                                        "rsID / SPDI / HGVS / CA-ID."},
                "r2_threshold":       {"type": "number"},
                "d_prime_threshold":  {"type": "number"},
                "ancestry":           {**_S_STRING, "description":
                                        "Comma list AFR,AMR,EAS,EUR,SAS."},
                "limit":              {**_S_INTEGER, "default": 100},
                "verbose":            {**_S_BOOLEAN, "default": False},
            },
            "required": ["variant_id"],
        },
        cli=["catalog", "find-ld"],
        flag_map={"variant_id": "", "r2_threshold": "--r2-threshold",
                   "d_prime_threshold": "--d-prime-threshold",
                   "ancestry": "--ancestry", "limit": "--limit",
                   "verbose": "--verbose"},
    ),
    _T(
        "catalog_resolve_id",
        "★ IGVF CATALOG ID CROSS-REFERENCE PROJECTION ★. Translates one "
        "ID into all of its cross-references — rsID ↔ SPDI ↔ HGVS ↔ "
        "CA-ID for variants; symbol ↔ ENSG ↔ HGNC ↔ Entrez ↔ synonyms "
        "for genes; UniProt ↔ ENSP for proteins; etc. USE THIS as a "
        "preflight when an analysis needs identifiers in a specific "
        "namespace.",
        {
            "type": "object",
            "properties": {
                "id": {**_S_STRING, "description":
                        "Any IGVF Catalog ID."},
            },
            "required": ["id"],
        },
        cli=["catalog", "resolve-id"],
        flag_map={"id": ""},
    ),
    # ──────────────────────────────────────────────────────────────────
    # ChIP-Atlas (Ohta/Oki — chip-atlas.org)
    # ──────────────────────────────────────────────────────────────────
    _T(
        "chipatlas_list_antigens",
        "★ CHIP-ATLAS ANTIGEN BROWSER ★. Lists every antigen (histone "
        "mark / TF / ATAC-Seq / DNase-Seq / Bisulfite-Seq class) that "
        "ChIP-Atlas has reprocessed experiments for, scoped to a "
        "(genome × agClass × cellClass) slice, with per-antigen "
        "experiment counts. USE THIS to answer 'does ChIP-Atlas have X "
        "ChIP-seq in cell-type Y, and how many?' before pulling files.",
        {
            "type": "object",
            "properties": {
                "genome":     {**_S_STRING, "description":
                                "hg38 / hg19 / mm10 / mm9 / rn6 / dm6 / dm3 / "
                                "ce11 / ce10 / sacCer3"},
                "ag_class":   {**_S_STRING, "description":
                                "Histone / 'TFs and others' / 'RNA polymerase' / "
                                "ATAC-Seq / DNase-seq / Bisulfite-Seq / etc."},
                "cell_class": {**_S_STRING, "default": "All cell types"},
                "limit":      {**_S_INTEGER, "default": 40},
            },
            "required": ["genome", "ag_class"],
        },
        cli=["chipatlas", "list-antigens"],
        flag_map={"genome": "--genome", "ag_class": "--ag-class",
                   "cell_class": "--cell-class", "limit": "--limit"},
    ),
    _T(
        "chipatlas_search",
        "★ CHIP-ATLAS FREE-TEXT EXPERIMENT SEARCH ★. Search across "
        "hundreds of thousands of reprocessed public ChIP-seq / "
        "ATAC-seq / DNase-seq / Bisulfite-seq SRX experiments by any "
        "term (TF name, cell line, condition, GSE/PRJNA). Returns SRX "
        "accessions + antigen + cell-type + title. USE THIS when the "
        "user wants 'find me CTCF ChIP-seq in K562' or any similar "
        "name/keyword query, before drilling into per-experiment files.",
        {
            "type": "object",
            "properties": {
                "query":  {**_S_STRING, "description":
                            "Free-text query, e.g. 'CTCF K562'."},
                "genome": {**_S_STRING, "description": "Optional genome filter."},
                "limit":  {**_S_INTEGER, "default": 25},
            },
            "required": ["query"],
        },
        cli=["chipatlas", "search"],
        flag_map={"query": "--query", "genome": "--genome", "limit": "--limit"},
    ),
    _T(
        "chipatlas_get_experiment",
        "★ CHIP-ATLAS EXPERIMENT METADATA ★. Full metadata for one "
        "SRX/DRX/ERX accession — antigen, cell type, taxon, processing "
        "stats. USE THIS to inspect a specific experiment that "
        "chipatlas_search returned.",
        {
            "type": "object",
            "properties": {
                "experiment_id": {**_S_STRING,
                                    "description": "SRX/DRX/ERX accession."},
            },
            "required": ["experiment_id"],
        },
        cli=["chipatlas", "get-experiment"],
        flag_map={"experiment_id": ""},
    ),
    _T(
        "chipatlas_download_experiment",
        "★ CHIP-ATLAS PER-EXPERIMENT FILE PULL ★. For a given SRX + "
        "genome, download (or just enumerate the URLs of) the "
        "per-experiment files at the requested kinds. Kinds: bw "
        "(BigWig signal), bb (all-peaks BigBed), bb05/bb10/bb20 "
        "(BigBed at -log10(q)=5/10/20 thresholds), bed05/bed10/bed20 "
        "(plain BED equivalents). Pass urls_only=true to skip the "
        "actual download.",
        {
            "type": "object",
            "properties": {
                "experiment_id": {**_S_STRING},
                "genome":        {**_S_STRING},
                "kinds":         {**_S_STRING, "default": "bw,bb05",
                                   "description":
                                   "Comma-list (bw / bb / bb05 / bb10 / bb20 / "
                                   "bed05 / bed10 / bed20)."},
                "urls_only":     {**_S_BOOLEAN, "default": False},
            },
            "required": ["experiment_id", "genome"],
        },
        cli=["chipatlas", "download-experiment"],
        flag_map={"experiment_id": "", "genome": "--genome",
                   "kinds": "--kinds", "urls_only": "--urls-only"},
    ),
    _T(
        "chipatlas_assemble_bed",
        "★ CHIP-ATLAS ASSEMBLED ALL-PEAKS BED ★. POST a "
        "(genome × agClass × antigen × cellClass × cellSubclass × qval) "
        "tuple to ChIP-Atlas's /download endpoint to get the URL of an "
        "assembled all-peaks BED that unions every reprocessed peak "
        "call matching the slice. Pass fetch=true to also stream the "
        "BED to disk. USE THIS to get a single TF-in-cellClass peakset "
        "without manually concatenating per-experiment files.",
        {
            "type": "object",
            "properties": {
                "genome":        {**_S_STRING},
                "ag_class":      {**_S_STRING},
                "antigen":       {**_S_STRING},
                "cell_class":    {**_S_STRING, "default": "All cell types"},
                "cell_subclass": {**_S_STRING},
                "qval":          {**_S_STRING, "default": "05"},
                "fetch":         {**_S_BOOLEAN, "default": False},
                "max_bytes":     {**_S_INTEGER},
            },
            "required": ["genome", "ag_class"],
        },
        cli=["chipatlas", "assemble-bed"],
        flag_map={"genome": "--genome", "ag_class": "--ag-class",
                   "antigen": "--antigen", "cell_class": "--cell-class",
                   "cell_subclass": "--cell-subclass", "qval": "--qval",
                   "fetch": "--fetch", "max_bytes": "--max-bytes"},
    ),
    _T(
        "chipatlas_target_genes",
        "★ CHIP-ATLAS PRE-COMPUTED TARGET-GENES TABLES ★. Discover "
        "(list=true) which antigens have pre-computed Target-Genes "
        "tables for a genome, or fetch one for a specific antigen at a "
        "given TSS-proximity distance (default 5000 bp). Each row "
        "contains a target gene with its mean peak score across all "
        "experiments. USE THIS to get TF→target-gene scoring without "
        "running your own peak overlap.",
        {
            "type": "object",
            "properties": {
                "genome":   {**_S_STRING},
                "list":     {**_S_BOOLEAN, "default": False},
                "antigen":  {**_S_STRING},
                "distance": {**_S_INTEGER, "default": 5000},
                "limit":    {**_S_INTEGER, "default": 50},
            },
            "required": ["genome"],
        },
        cli=["chipatlas", "target-genes"],
        flag_map={"genome": "--genome", "list": "--list",
                   "antigen": "--antigen", "distance": "--distance",
                   "limit": "--limit"},
    ),
    _T(
        "chipatlas_showcase",
        "★ CHIP-ATLAS ONE-COMMAND DEMO ★. End-to-end probe of the "
        "ChIP-Atlas surface: genomes, top histone antigens for a "
        "cell-class slice, bulk allPeaks_light HEAD probe, per-"
        "experiment BigBed HEAD probe, and the count of antigens with "
        "Target-Genes tables. All HEAD/JSON only — no GB-scale downloads.",
        {
            "type": "object",
            "properties": {
                "genome":         {**_S_STRING, "default": "hg38"},
                "cell_class":     {**_S_STRING, "default": "Pluripotent stem cell"},
                "canonical_srx":  {**_S_STRING, "default": "SRX150531"},
            },
        },
        cli=["chipatlas", "showcase"],
        flag_map={"genome": "--genome", "cell_class": "--cell-class",
                   "canonical_srx": "--canonical-srx"},
    ),
    _T(
        "chipatlas_submit_enrichment",
        "★ CHIP-ATLAS WABI ENRICHMENT JOB ★. Submit a gene-list or "
        "BED-region over-representation job to the NIG/DDBJ WABI queue. "
        "mode='genes' asks 'which TFs are enriched at the regulatory "
        "regions of my gene list?'; mode='regions' asks 'which TFs are "
        "enriched at my BED regions?'. Returns a job id you then pass "
        "to chipatlas_poll_enrichment.",
        {
            "type": "object",
            "properties": {
                "mode":       {**_S_STRING, "description": "'genes' or 'regions'"},
                "genome":     {**_S_STRING},
                "query":      {**_S_STRING,
                                "description": "File path OR literal text content."},
                "background": {**_S_STRING},
                "ag_class":   {**_S_STRING, "default": "TFs and others"},
                "cell_class": {**_S_STRING, "default": "All cell types"},
                "qval":       {**_S_STRING, "default": "05"},
                "distance":   {**_S_INTEGER, "default": 5000},
                "label":      {**_S_STRING},
            },
            "required": ["mode", "genome", "query"],
        },
        cli=["chipatlas", "submit-enrichment"],
        flag_map={"mode": "--mode", "genome": "--genome", "query": "--query",
                   "background": "--background", "ag_class": "--ag-class",
                   "cell_class": "--cell-class", "qval": "--qval",
                   "distance": "--distance", "label": "--label"},
    ),

    _T(
        "catalog_list_sources",
        "Enumerate IGVF Catalog edge endpoints by semantic category, "
        "or probe one endpoint for its observed sources / methods.",
        {
            "type": "object",
            "properties": {
                "category": {**_S_STRING},
                "endpoint": {**_S_STRING},
            },
        },
        cli=["catalog", "list-sources"],
        flag_map={"category": "--category", "endpoint": "--endpoint"},
    ),

    _T(
        "portal_list_types",
        "List the canonical IGVF CamelCase ItemTypes and their "
        "snake-case collection paths (e.g. MeasurementSet ↔ "
        "/measurement-sets/). Useful as a quick lookup before calling "
        "portal_search or portal_schema.",
        {"type": "object", "properties": {}},
        cli=["portal", "list-types"],
    ),

    _T(
        "enrich_showcase",
        "★ ONE-COMMAND ENRICHMENT DEMO ★. Runs ORA on a curated 47-gene "
        "cell-cycle / G2-M-checkpoint list (CCN*, CDK*, CDKN*, MCM*, "
        "AURK*, PLK*, BUB*, etc.) and writes a composite figure + "
        "narrative report. Positive-control validation that the skill "
        "is healthy — expected strong 'Cell Cycle' / 'G2-M Checkpoint' "
        "enrichment across Reactome, KEGG, MSigDB Hallmark, and all "
        "three GO branches.",
        {
            "type": "object",
            "properties": {
                "libs":     {**_S_STRING, "default": "all"},
                "organism": {**_S_STRING, "default": "human"},
                "label":    {**_S_STRING},
                "top_k":    {**_S_INTEGER, "default": 8},
            },
        },
        cli=["enrich", "showcase"],
        flag_map={"libs": "--libs", "organism": "--organism",
                   "label": "--label", "top_k": "--top-k"},
    ),

    _T(
        "calibrate_thresholds",
        "★ ACMG/AMP EVIDENCE THRESHOLDS FROM A PRIOR ★. Solves Tavtigian's "
        "constant C (O_PVSt) for a given prior probability of pathogenicity "
        "and prints the likelihood-ratio (LR+) threshold for every evidence "
        "strength — supporting / moderate / strong / very strong, on both the "
        "PS3 (pathogenic) and BS3 (benign) sides. Instant, no data or fitting "
        "needed. USE THIS to answer 'how strong must an assay LR+ be to count "
        "as PS3 moderate?' or to sanity-check a published calibration.",
        {
            "type": "object",
            "properties": {
                "prior": {"type": "number", "description":
                    "P(pathogenic) in the reference population, e.g. 0.1 for a "
                    "well-studied disease gene, 0.01 for a low-prior gene."},
                "point_values": {**_S_STRING, "description":
                    "Comma list of evidence point values (default 1..8)."},
            },
            "required": ["prior"],
        },
        cli=["calibrate", "thresholds"],
        flag_map={"prior": "--prior", "point_values": "--point-values"},
    ),

    _T(
        "calibrate_run",
        "★ CALIBRATE A FUNCTIONAL ASSAY TO ACMG/AMP EVIDENCE ★. Full "
        "exCALIBR chain (Zeiberg et al. 2025) on a variant-effect scoreset: "
        "bootstrap constrained skew-normal mixture EM over the P/LP, B/LB, "
        "gnomAD-population and synonymous samples → EM prior → LR+(score) → "
        "Tavtigian C → the score window that earns each evidence strength "
        "(PS3 / BS3 supporting → very strong), plus 2c-vs-3c model selection, "
        "a calibration JSON and a figure. Input is either an IGVF / "
        "Pillar-format scoreset CSV (`pillar`) or a score/sample table from "
        "calibrate prepare (`table`). LONG JOB: runs to completion in one "
        "call with a progress heartbeat and a resumable ledger — lower "
        "`n_bootstraps` / `fits_per_bootstrap` for a quick look (defaults "
        "1000 x 100 are the paper's settings).",
        {
            "type": "object",
            "properties": {
                "pillar": {**_S_STRING, "description":
                    "Path to an IGVF / Pillar-format scoreset CSV (carries "
                    "auth_reported_score + ClinVar + gnomAD + SpliceAI)."},
                "table": {**_S_STRING, "description":
                    "Path to a score/sample CSV (from calibrate prepare)."},
                "name": {**_S_STRING, "description":
                    "Dataset name; also selects one Dataset from a "
                    "multi-dataset Pillar table."},
                "n_bootstraps": {**_S_INTEGER, "default": 1000},
                "fits_per_bootstrap": {**_S_INTEGER, "default": 100},
                "benign_method": {**_S_STRING, "description":
                    "'avg' (benign+synonymous, default), 'benign', "
                    "'synonymous'."},
                "jobs": {**_S_INTEGER, "default": -1,
                          "description": "-1 = all CPUs."},
            },
        },
        cli=["calibrate", "run"],
        flag_map={"pillar": "--pillar", "table": "--table", "name": "--name",
                   "n_bootstraps": "--n-bootstraps",
                   "fits_per_bootstrap": "--fits-per-bootstrap",
                   "benign_method": "--benign-method", "jobs": "--jobs"},
    ),

    _T(
        "calibrate_assign",
        "Apply an existing calibration to a table of assay scores: each "
        "score gets its evidence point value and the ACMG/AMP criterion + "
        "strength it supports (e.g. 'PS3 moderate', 'BS3 strong', 'no "
        "evidence'). USE THIS after calibrate_run, or with a published "
        "calibration JSON, to interpret new variants.",
        {
            "type": "object",
            "properties": {
                "calibration": {**_S_STRING, "description":
                    "Path to a *_calibration.json from calibrate run."},
                "scores": {**_S_STRING, "description":
                    "CSV/TSV with a score column (plus any ID columns)."},
            },
            "required": ["calibration", "scores"],
        },
        cli=["calibrate", "assign"],
        flag_map={"calibration": "--calibration", "scores": "--scores"},
    ),

    _T(
        "calibrate_prepare",
        "Build the score/sample input table a calibration needs. Labels "
        "variants into the four calibration samples — 0 P/LP and 1 B/LB from "
        "ClinVar (with a review-star quality gate), 2 population from gnomAD "
        "membership, 3 synonymous — from an IGVF / Pillar-format scoreset "
        "(`pillar`), a plain score/sample CSV (`table`), or a `mavedb "
        "map-scoreset` output joined to ClinVar (`mapped` + `clinvar_tsv`). "
        "Reports the per-sample counts so you can see whether the scoreset "
        "has enough labelled variants to calibrate at all.",
        {
            "type": "object",
            "properties": {
                "pillar": {**_S_STRING},
                "table": {**_S_STRING},
                "mapped": {**_S_STRING, "description":
                    "mavedb map-scoreset output with chr/pos/ref/alt + score."},
                "clinvar_tsv": {**_S_STRING, "description":
                    "ClinVar variant_summary.txt.gz (used with `mapped`)."},
                "name": {**_S_STRING},
                "min_clinvar_star": {**_S_INTEGER, "default": 1},
            },
        },
        cli=["calibrate", "prepare"],
        flag_map={"pillar": "--pillar", "table": "--table",
                   "mapped": "--mapped", "clinvar_tsv": "--clinvar-tsv",
                   "name": "--name",
                   "min_clinvar_star": "--min-clinvar-star"},
    ),

]


_BY_NAME = {t.name: t for t in _TOOLS}


# --------------------------- User-extension tools ----------------------------


def _merge_user_tools() -> None:
    """Absorb user-defined tools into the registry at import time.

    Users drop YAML / JSON manifests under ``~/.igvfagent/tools/`` or
    ``<root>/UserExtensions/tools/`` (see ``_userext``) and they show
    up in ``igvfagent tools``, the ``ask`` agent loop, and the UI tool
    picker exactly like built-ins. Defensive by design: a broken
    manifest is skipped with a warning and can never take down the
    built-in registry.
    """
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        try:
            from igvfagent import _userext
        except ImportError:
            import _userext  # type: ignore[no-redef]
        specs = _userext.discover_tools()
    except Exception as exc:
        logger.warning("user-extension tool discovery failed: %s", exc)
        return
    for spec in specs:
        if spec["name"] in _BY_NAME:
            logger.warning("user tool `%s` (%s) shadows an existing tool; "
                           "skipped", spec["name"], spec.get("source"))
            continue
        tool = Tool(
            name=spec["name"], description=spec["description"],
            parameters=spec["parameters"], cli=list(spec["cli"]),
            positional=list(spec["positional"]),
            flag_map=dict(spec["flag_map"]),
            flag_repeat=set(spec["flag_repeat"]),
            bool_flags=set(spec["bool_flags"]),
            command=list(spec["command"]),
        )
        _TOOLS.append(tool)
        _BY_NAME[tool.name] = tool


_merge_user_tools()


def refresh_user_tools() -> int:
    """Re-scan the extension directories and absorb any new user tools.

    For long-lived processes (the Streamlit UI) where the import-time
    merge already happened: call after a manifest is added on disk.
    Existing names are never redefined. Returns the number of tools added.
    """
    before = len(_TOOLS)
    _merge_user_tools()
    return len(_TOOLS) - before


# --------------------------- Public registry API ----------------------------

def list_tools() -> "list[Tool]":
    """Return all registered tools."""
    return list(_TOOLS)


def get_tool(name: str) -> Optional[Tool]:
    return _BY_NAME.get(name)


def _llm_module():
    """Import ``_llm`` whether we're running as a package submodule or as
    a top-level script (when ``Scripts/`` is on ``sys.path``)."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        from igvfagent import _llm as mod  # installed package
    except Exception:
        import _llm as mod  # type: ignore[no-redef]
    return mod


def to_anthropic_schema() -> "list[dict]":
    return _llm_module().to_anthropic_tools(t.to_dict() for t in _TOOLS)


def to_openai_schema() -> "list[dict]":
    return _llm_module().to_openai_tools(t.to_dict() for t in _TOOLS)


# --------------------------- Argument materialization -----------------------


def _coerce_value(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


def _build_argv(tool: Tool, arguments: dict) -> "list[str]":
    """Translate the parameter dict the LLM provided into ``igvfagent``
    argv tokens (or the user tool's own ``command`` argv)."""
    argv = [*(tool.command or ["igvfagent"]), *tool.cli]
    args = dict(arguments or {})

    # Positional first
    for p in tool.positional:
        if p in args and args[p] is not None and args[p] != "":
            argv.append(str(args.pop(p)))

    for name, value in args.items():
        if value is None or value == "":
            continue
        if name in tool.bool_flags:
            if bool(value):
                flag = tool.flag_map.get(name, "--" + name.replace("_", "-"))
                argv.append(flag)
            continue
        if name in tool.flag_repeat:
            flag = tool.flag_map.get(name, "--" + name.replace("_", "-"))
            for v in (value if isinstance(value, (list, tuple)) else [value]):
                argv.extend([flag, _coerce_value(v)])
            continue
        flag = tool.flag_map.get(name, "--" + name.replace("_", "-"))
        # Convention: ``flag_map={name: ""}`` means the argument is
        # positional. Skip emitting an empty-string flag token before
        # the value — that would make argparse choke with exit_code=2.
        if flag == "":
            argv.append(_coerce_value(value))
        else:
            argv.extend([flag, _coerce_value(value)])
    return argv


# --------------------------- Execution --------------------------------------

# Lines that announce an artefact path. Two patterns:
#   1. A curated allow-list of explicit prefixes ("Report:", "Manifest:", …)
#   2. A generic fallback: any line matching ``<Label>: <path>`` where the
#      path ends in a known viewable extension. This catches one-off
#      announcements like ``Browser SVG:``, ``Output:``, ``Saved:`` that
#      individual skills print without us having to enumerate every one.
_REPORT_RE   = re.compile(
    r"^(?:Report|Lit manifest|IGVF manifest|Manifest|Evidence pack|"
    r"Local KG|Plot|Plots|Wrote|Wrote report|Wrote: |Playbook|"
    r"Browser SVG|SVG|Output|Saved|Downloaded|Pulled|"
    r"rE2G linkage table|CSV|TSV|JSON|Figure|Figures):\s*(.+?)\s*$",
    re.M,
)
_REPORT_BY_EXT_RE = re.compile(
    r"^\s*[A-Z][\w \-/]{0,40}:\s*"
    r"(\S+?\.(?:csv|tsv|json|jsonl|md|svg|png|jpg|jpeg|gif|pdf|html|h5ad))"
    r"\s*$", re.M,
)


def _resolve_igvfagent() -> "list[str]":
    """Find an ``igvfagent`` invocation that works in this environment.

    Prefers the installed console script (so the audit trail mirrors the
    end-user shell). Falls back to ``python -m igvfagent.cli`` when the
    package is importable but no console script is on PATH (e.g. running
    out of a checkout without ``pip install -e .``).
    """
    import shutil
    binary = shutil.which("igvfagent")
    if binary:
        return [binary]
    return [sys.executable, "-m", "igvfagent.cli"]


def execute(name: str, arguments: dict, *, timeout: Optional[float] = None,
             cwd: Optional[str] = None,
             extra_env: Optional[dict] = None) -> dict:
    """Run the tool's wrapped CLI as a subprocess.

    Returns a dict with ``exit_code``, ``stdout``, ``stderr``, and any
    ``Report:`` / ``Manifest:`` / ``Wrote:`` artefact paths the wrapped
    skill announced on stdout. Raises ``KeyError`` if ``name`` is unknown.
    """
    tool = get_tool(name)
    if not tool:
        raise KeyError(f"Unknown tool: {name}")

    # Evaluation hook: force named tools to fail so recovery becomes
    # observable (Scripts/eval_tiers_skill.py, Tier 2). A run that cannot
    # proceed when an archive errors has a planning defect, not a skill
    # defect, and nothing else in the suite exercises that path. Inert
    # unless IGVF_EVAL_FAIL_TOOL is set, and it is never set in normal use.
    _fail = os.environ.get("IGVF_EVAL_FAIL_TOOL", "")
    if _fail and name in {t.strip() for t in _fail.split(",") if t.strip()}:
        logger.warning("tool=%s failed by IGVF_EVAL_FAIL_TOOL (evaluation)", name)
        return {"name": name, "argv": [], "exit_code": 1, "stdout": "",
                "stderr": f"injected failure for evaluation: {name} "
                          f"unavailable (simulated upstream error)",
                "timed_out": False, "artifacts": {}}

    argv = _build_argv(tool, arguments or {})
    # Built-in tools carry a leading "igvfagent" placeholder — replace it
    # with whatever runner works here (binary or `python -m igvfagent.cli`).
    # User tools with their own `command` run that argv verbatim.
    if not tool.command:
        argv = _resolve_igvfagent() + argv[1:]

    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)
    logger.info("tool=%s argv=%s", name, shlex.join(argv))
    try:
        proc = subprocess.run(
            argv, cwd=cwd, env=env, capture_output=True,
            text=True, timeout=timeout, check=False,
        )
    except subprocess.TimeoutExpired as e:
        return {"name": name, "argv": argv, "exit_code": 124,
                "stdout": (e.stdout or ""), "stderr": (e.stderr or ""),
                "timed_out": True, "artifacts": {}}

    artefacts = _parse_artefacts(proc.stdout)
    return {
        "name":      name,
        "argv":      argv,
        "exit_code": proc.returncode,
        "stdout":    proc.stdout,
        "stderr":    proc.stderr,
        "timed_out": False,
        "artifacts": artefacts,
    }


def _parse_artefacts(stdout: str) -> "dict[str, list[str]]":
    artefacts: "dict[str, list[str]]" = {}
    seen: "set[str]" = set()

    def _add(key: str, value: str) -> None:
        if value in seen:
            return
        seen.add(value)
        artefacts.setdefault(key.strip().lower().replace(" ", "_"),
                              []).append(value)

    for m in _REPORT_RE.finditer(stdout or ""):
        line = m.group(0)
        key, _, value = line.partition(":")
        _add(key, value.strip())
    # Generic fallback: catch any "Label: <path-with-known-ext>" line that
    # the explicit prefix list didn't cover (e.g. ``Browser SVG``,
    # ``Output``, ``Saved``).
    for m in _REPORT_BY_EXT_RE.finditer(stdout or ""):
        path = m.group(1).strip()
        # Recover the label so the artefact bucket is human-readable.
        line = m.group(0)
        label, _, _ = line.partition(":")
        _add(label, path)
    return artefacts


# --------------------------- Pretty render ---------------------------------

def render_tool_summary(tool: Tool) -> str:
    runner = (shlex.join(tool.command + tool.cli) if tool.command
              else "igvfagent " + " ".join(tool.cli))
    return (f"{tool.name}\n  cli: {runner}\n"
            f"  desc: {tool.description}\n"
            f"  params: {json.dumps(tool.parameters.get('properties', {}), default=str)[:200]}")


__all__ = [
    "Tool", "list_tools", "get_tool",
    "to_anthropic_schema", "to_openai_schema",
    "execute", "refresh_user_tools",

]
