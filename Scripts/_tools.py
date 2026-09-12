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
_S_NUMBER  = {"type": "number"}


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
        "pathway_db_query",
        "Pathways for a gene list, from the locally integrated KEGG + "
        "Reactome + WikiPathways release. Faster and broader than querying a "
        "single pathway API, and reports WHICH databases assert each "
        "pathway, so agreement between them is visible. Use this to ask "
        "which pathways a gene set shares. If nothing is cached yet, call "
        "pathway_db_refresh first.",
        {
            "type": "object",
            "properties": {
                "genes":   {**_S_STRING, "description":
                            "Gene symbols, comma-separated, e.g. "
                            "'CLU,BIN1,PICALM'."},
                "sources": {**_S_STRING, "description":
                            "Comma list: kegg, reactome, wikipathways. "
                            "Default all three."},
                "top":     {**_S_INTEGER, "description":
                            "How many pathways to show (default 25)."},
                "verbose": {**_S_BOOLEAN, "description":
                            "Also list which query genes hit each pathway."},
            },
            "required": ["genes"],
        },
        ["pathwaydb", "query"],
        flag_map={"genes": "--genes", "sources": "--sources",
                   "top": "--top", "verbose": "--verbose"},
        bool_flags=("verbose",),
    ),

    _T(
        "pathway_db_refresh",
        "Download the CURRENT KEGG, Reactome and WikiPathways releases, "
        "normalise gene identifiers through NCBI gene_info, merge pathways "
        "that different databases describe differently, and load the result "
        "into the local knowledge graph. Use when pathway data is missing or "
        "stale, or when the user asks to pull/refresh pathway databases. "
        "Downloads ~180 MB the first time; later runs reuse the cache.",
        {
            "type": "object",
            "properties": {
                "sources":   {**_S_STRING, "description":
                              "Comma list: kegg, reactome, wikipathways. "
                              "Default all three."},
                "relations": {**_S_BOOLEAN, "description":
                              "Also build typed gene-gene relations from "
                              "KEGG KGML and Reactome interactions. Adds a "
                              "few minutes on the first run."},
                "min_sources": {**_S_INTEGER, "description":
                                "Only store facts asserted by at least N "
                                "databases (default 1)."},
                "label":     {**_S_STRING, "description": "Run label."},
                "no_kg":     {**_S_BOOLEAN, "description":
                              "Write files without touching the KG."},
            },
            "required": [],
        },
        ["pathwaydb", "build"],
        flag_map={"sources": "--sources", "relations": "--relations",
                   "min_sources": "--min-sources", "label": "--label",
                   "no_kg": "--no-kg"},
        bool_flags=("relations", "no_kg"),
    ),

    _T(
        "pathway_network_viz",
        "★ VISUALISE PATHWAY / INTERACTION NETWORKS FOR A GENE LIST ★. Pass "
        "the gene symbols INLINE in `genes` (comma-separated) — no file "
        "staging is needed and no file-writing tool exists. Draws a network "
        "from STRING protein interactions, Reactome and KEGG pathway "
        "membership, and the local knowledge graph, and downloads KEGG's own "
        "rendered pathway diagrams. THIS IS THE RIGHT TOOL whenever someone "
        "asks to visualise pathways, networks or interactions for a set of "
        "genes. Produces PNG + SVG, an edge-list CSV, and enrichment.",
        {
            "type": "object",
            "properties": {
                "genes":    {**_S_STRING, "description":
                             "Gene symbols, comma-separated, e.g. "
                             "'CLU,BIN1,PICALM'. A file path also works."},
                "sources":  {**_S_STRING, "description":
                             "Comma list: string, reactome, kegg, local. "
                             "Default all four."},
                "min_score": {**_S_STRING, "description":
                              "STRING confidence cutoff, default 0.4."},
                "layout":   {**_S_STRING, "description":
                             "spring | circular | kamada"},
                "label":    {**_S_STRING, "description": "Run label."},
            },
            "required": ["genes"],
        },
        ["pathway-viz", "network"],
        flag_map={"genes": "--genes", "sources": "--sources",
                   "min_score": "--min-score", "layout": "--layout",
                   "label": "--label"},
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
        "audit_manifest",
        "★ AUDIT AN UPLOADED OR GENERATED MANIFEST — USE THIS INSTEAD OF "
        "WRITING A PASS/FAIL TABLE YOURSELF ★ Checks row count, column "
        "uniqueness, one-to-one pairing between two columns, and "
        "cross-tabulated coverage. Its point is that it reports THREE "
        "outcomes, not two: pass, fail, and **limitation**. A `limitation` is "
        "a VALID manifest truthfully recording something the upstream source "
        "does not provide, or a study design that cannot support a given "
        "comparison — it is NOT a validation failure and must never be "
        "presented as one. Reporting a correctly-recorded absence as a failed "
        "check tells a user their file is broken when it is not: a GSE213151 "
        "audit showed `raw_rna_matrix_listed=false` and "
        "`atac_peak_matrix_listed=false` as FAILURES when both values were "
        "correct, because GEO genuinely does not supply those files. Pass "
        "such columns as `absent_ok`; use `require_true` only for columns "
        "where false really is a defect. Nothing is written or downloaded.",
        {
            "type": "object",
            "properties": {
                "path": {**_S_STRING, "description":
                          "Workspace-relative path to the CSV/TSV manifest."},
                "unique": {**_S_STRING, "description":
                            "Comma-separated columns that must be unique, "
                            "e.g. sample_id,rna_gsm,atac_gsm"},
                "pair": {**_S_STRING, "description":
                          "'a:b' — require a one-to-one pairing, e.g. "
                          "rna_gsm:atac_gsm"},
                "group": {**_S_STRING, "description":
                           "Comma-separated columns to cross-tabulate for "
                           "coverage, e.g. cell_line,diff_day. Levels present "
                           "for only some groups are reported as limitations."},
                "absent_ok": {**_S_STRING, "description":
                               "Boolean columns whose false values are a DATA "
                               "LIMITATION rather than a failure."},
                "require_true": {**_S_STRING, "description":
                                  "Boolean columns where false IS a failure."},
            },
            "required": ["path"],
        },
        cli=["artifact", "audit"],
        flag_map={"path": "--path", "unique": "--unique", "pair": "--pair",
                   "group": "--group", "absent_ok": "--absent-ok",
                   "require_true": "--require-true"},
    ),
    _T(
        "rank_artifact",
        "★ USE THIS FOR ANY 'TOP', 'HIGHEST', 'BEST' OR 'RANKED' CLAIM ★ — "
        "sorts a CSV/TSV artefact by a numeric column across the ENTIRE file "
        "and returns the true top N. **Never use grep_artifacts or "
        "read_artifact to answer a ranking question**: both return a bounded "
        "sample, and rows drawn from a sample are NOT the highest-scoring "
        "rows. Doing so has produced answers that named a 'top 3' whose "
        "scores were 0.99 when records at 0.9999999981 existed in the same "
        "file — the rows were real and correctly attributed, but the ranking "
        "claim was false. This tool reports n_scanned and n_ranked so the "
        "answer can state how many records the ranking actually covered. Use "
        "`where` (comma-separated, keeps a row matching ANY term) and "
        "`exclude` to filter BEFORE ranking: 'highest-scoring kidney record' "
        "is where='kidney,renal' exclude='adrenal', because 'kidney' alone "
        "misses 'renal cortical epithelial cell' and 'renal' alone also "
        "matches 'adrenal gland'.",
        {
            "type": "object",
            "properties": {
                "path": {**_S_STRING, "description":
                          "Workspace-relative path to the CSV/TSV artefact."},
                "column": {**_S_STRING, "description":
                            "Numeric column to rank by, e.g. score."},
                "n": {**_S_INTEGER, "description": "How many rows. Default 10."},
                "where": {**_S_STRING, "description":
                           "Comma-separated include terms; a row matching ANY "
                           "is kept."},
                "where_column": {**_S_STRING, "description":
                                  "Restrict the include/exclude match to one "
                                  "column, e.g. biosample."},
                "exclude": {**_S_STRING, "description":
                             "Comma-separated terms; drop rows matching any."},
                "ascending": {**_S_BOOLEAN, "description":
                               "Rank lowest-first instead."},
            },
            "required": ["path", "column"],
        },
        cli=["artifact", "top"],
        flag_map={"path": "--path", "column": "--column", "n": "--n",
                   "where": "--where", "where_column": "--where-column",
                   "exclude": "--exclude"},
        bool_flags=("ascending",),
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
        "kg_genes_batch",
        "★ SEVERAL GENES' REGULATORY EVIDENCE IN ONE CALL — USE THIS FOR ANY "
        "MULTI-GENE QUESTION ★ Do NOT loop kg_gene over a gene list. A "
        "six-gene question driven gene-by-gene took over 25 minutes and did "
        "not finish, because each gene costs many agent iterations at 5-12 "
        "seconds each; this does all six in ONE call in about 37 seconds. "
        "For every gene it resolves the target gene id, pulls "
        "element-to-gene edges exhaustively, keeps only rows whose OWN "
        "target gene matches (the endpoint returns edges for every gene near "
        "the locus, and most are other genes), separates observed "
        "measurements from predictions, applies the tissue filter, and "
        "writes one combined manifest. Pass `tissue` as a comma-separated "
        "list: 'kidney' alone misses 'renal cortical epithelial cell', and "
        "'renal' alone also matches 'adrenal gland', so kidney questions "
        "want tissue='kidney,renal' with exclude_tissue='adrenal'. The "
        "result reports `catalog_retrieval` per gene, which says whether the "
        "CATALOG traversal was exhaustive — it is NOT a statement about how "
        "much of the manifest a later read covered. To report highest-scoring "
        "rows, rank the combined manifest with rank_artifact; never quote a "
        "ranking from an excerpt.",
        {
            "type": "object",
            "properties": {
                "symbols": {**_S_STRING, "description":
                             "Comma-separated gene symbols, e.g. "
                             "PAX2,LHX1,WT1,HNF4A,GATA3,SOX9"},
                "tissue": {**_S_STRING, "description":
                            "Comma-separated include terms; a row matching "
                            "ANY is kept, e.g. kidney,renal"},
                "exclude_tissue": {**_S_STRING, "description":
                                    "Comma-separated terms to drop. Defaults "
                                    "to adrenal."},
                "max_pages": {**_S_INTEGER},
                "label": {**_S_STRING},
            },
            "required": ["symbols"],
        },
        cli=["kg", "genes"],
        positional=("symbols",),
        flag_map={"tissue": "--tissue", "exclude_tissue": "--exclude-tissue",
                   "max_pages": "--max-pages", "label": "--label"},
    ),
    _T(
        "kg_gene",
        "Comprehensive gene context from the IGVF Catalog KG: variants, "
        "transcripts, proteins, regulatory elements, diseases, pathways. "
        "Default 'tell me about gene X' tool. The Catalog is queried first "
        "and FAVOR then supplements its variants with CADD / GERP / "
        "conservation / ClinVar; optional flags add enhancer-gene linkage, "
        "single-cell datasets, prior literature.",
        {
            "type": "object",
            "properties": {
                "symbol":   {**_S_STRING, "description": "Gene symbol, e.g. APOE"},
                "depth":    {**_S_INTEGER, "default": 1,
                              "description": "1 = direct relations only; 2 = also fan out per variant."},
                "limit":    {**_S_INTEGER, "default": 25},
                "no_favor":        {**_S_BOOLEAN, "default": False,
                    "description": "Skip the FAVOR supplement. FAVOR runs by "
                        "default AFTER the IGVF Catalog and adds CADD / GERP / "
                        "conservation / ClinVar to the variants the Catalog "
                        "returned; set true to use Catalog data only."},
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
            "no_favor": "--no-call-favor",
            "literature_context": "--literature-context",
            "label": "--label",
        },
        flag_repeat={"literature_context"},
        bool_flags={"no_favor", "call_linkage",
                     "call_singlecell", "call_literature"},
    ),

    _T(
        "kg_variant",
        "Variant-centric KG: linked genes, regulatory elements, "
        "phenotypes, biosamples, predictions, plus a FAVOR supplement "
        "(CADD / GERP / conservation / ClinVar) layered on the Catalog "
        "record. Accepts rsID/SPDI/HGVS. "
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
                "no_favor":        {**_S_BOOLEAN, "default": False,
                    "description": "Skip the FAVOR supplement. FAVOR runs by "
                        "default AFTER the IGVF Catalog and adds CADD / GERP / "
                        "conservation / ClinVar to the variants the Catalog "
                        "returned; set true to use Catalog data only."},
                "call_literature": {**_S_BOOLEAN, "default": False},
                "literature_context": {**_S_ARRAY_S},
                "label": {**_S_STRING},
            },
            "required": ["variant"],
        },
        cli=["kg", "variant"],
        positional=["variant"],
        flag_map={"limit": "--limit", "no_favor": "--no-call-favor",
                   "literature_context": "--literature-context",
                   "label": "--label"},
        flag_repeat={"literature_context"},
        bool_flags={"no_favor", "call_literature"},
    ),

    _T(
        "kg_region",
        "Region-centric KG: genes + cCREs + enhancer-gene linkage in "
        "window, plus FAVOR annotations for the Catalog's variants in "
        "that window. Format chr19:44903000-44912000.",
        {
            "type": "object",
            "properties": {
                "region": {**_S_STRING},
                "limit":  {**_S_INTEGER, "default": 50},
                "no_favor":        {**_S_BOOLEAN, "default": False,
                    "description": "Skip the FAVOR supplement. FAVOR runs by "
                        "default AFTER the IGVF Catalog and adds CADD / GERP / "
                        "conservation / ClinVar to the variants the Catalog "
                        "returned; set true to use Catalog data only."},
                "label":  {**_S_STRING},
            },
            "required": ["region"],
        },
        cli=["kg", "region"],
        positional=["region"],
        flag_map={"limit": "--limit", "label": "--label",
                   "no_favor": "--no-call-favor"},
        bool_flags={"no_favor"},
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
    # MPRA oligo LIBRARY DESIGN (the wet-lab step before `mpra activity`).
    # Python port of kircherlab/MPRAOligoDesign (Max Schubach, BIH; MIT),
    # reimplemented on the standard library — no Snakemake/conda/bedtools.
    # ──────────────────────────────────────────────────────────────────
    _T(
        "oligo_tile",
        "★ MPRA TILING: REGIONS → OLIGO WINDOWS ★. Turns a region BED into "
        "oligo-sized windows using the three upstream strategies, picked "
        "automatically by region length: a region shorter than "
        "--centering-max becomes ONE centred, padded oligo; a mid-length "
        "region becomes TWO tiles pinned to its ends; anything longer is "
        "tiled from the centre outwards with at least --min-overlap bp of "
        "overlap so no base falls in a gap. Emits regions.tiles.bed.gz with "
        "tile names '<region>_<fwd|rev|none>_tile<i>-<n>'. USE THIS before "
        "designing sequences whenever regions are not already oligo-sized.",
        {
            "type": "object",
            "properties": {
                "regions":       {**_S_STRING, "description":
                                   "Region BED (.bed or .bed.gz)."},
                "oligo_length":  {**_S_INTEGER, "default": 200},
                "min_overlap":   {**_S_INTEGER, "default": 50},
                "centering_max": {**_S_INTEGER, "description":
                                   "Regions <= this get one centred oligo."},
                "two_tiles_max": {**_S_INTEGER},
                "variant_edge_exclusion": {**_S_INTEGER, "default": 20},
                "include_variant_edge":   {**_S_BOOLEAN, "default": False},
                "label":         {**_S_STRING},
            },
            "required": ["regions"],
        },
        cli=["oligo", "tile"],
        flag_map={"regions": "--regions", "oligo_length": "--oligo-length",
                   "min_overlap": "--min-overlap",
                   "centering_max": "--centering-max",
                   "two_tiles_max": "--two-tiles-max",
                   "variant_edge_exclusion": "--variant-edge-exclusion",
                   "label": "--label"},
        bool_flags=("include_variant_edge",),
    ),
    _T(
        "oligo_design_variants",
        "★ MPRA REF/ALT VARIANT OLIGO DESIGN ★. For every variant in a VCF, "
        "finds the region(s) that contain it and writes a matched pair of "
        "oligos differing at exactly that position — REF_<region> and "
        "ALT_<region>_<variant> — which is the contrast an allelic MPRA "
        "actually measures. Handles SNVs, MNVs and indels (a deletion's REF "
        "window is extended so both oligos stay the same length), and "
        "reverse-complements '-' strand regions. --variant-edge-exclusion "
        "keeps variants away from the oligo edge where flanking context is "
        "truncated; variants that fit no region are written to "
        "variants.removed.vcf.gz rather than silently dropped. Emits "
        "design.fa plus Variant→Region→REF_ID→ALT_ID maps.",
        {
            "type": "object",
            "properties": {
                "regions":   {**_S_STRING},
                "variants":  {**_S_STRING, "description": "VCF (.vcf/.vcf.gz)."},
                "reference": {**_S_STRING, "description":
                               "Genome FASTA (uncompressed; .fai used if present)."},
                "variant_edge_exclusion": {**_S_INTEGER, "default": 20},
                "use_all_regions": {**_S_BOOLEAN, "default": False,
                                     "description":
                                     "Design against every matching region "
                                     "(default: only the most centred)."},
                "remove_regions_without_variants": {**_S_BOOLEAN,
                                                     "default": False},
                "label":     {**_S_STRING},
            },
            "required": ["regions", "variants", "reference"],
        },
        cli=["oligo", "design-variants"],
        flag_map={"regions": "--regions", "variants": "--variants",
                   "reference": "--reference",
                   "variant_edge_exclusion": "--variant-edge-exclusion",
                   "label": "--label"},
        bool_flags=("use_all_regions", "remove_regions_without_variants"),
    ),
    _T(
        "oligo_design_regions",
        "★ MPRA REGION OLIGO DESIGN ★. Extracts the oligo sequence for every "
        "region in a BED from the reference genome, reverse-complementing "
        "'-' strand regions. Fails loudly if any region is not exactly "
        "--oligo-length bp (run oligo_tile first, or pass no_length_check).",
        {
            "type": "object",
            "properties": {
                "regions":         {**_S_STRING},
                "reference":       {**_S_STRING},
                "oligo_length":    {**_S_INTEGER, "default": 200},
                "no_length_check": {**_S_BOOLEAN, "default": False},
                "label":           {**_S_STRING},
            },
            "required": ["regions", "reference"],
        },
        cli=["oligo", "design-regions"],
        flag_map={"regions": "--regions", "reference": "--reference",
                   "oligo_length": "--oligo-length", "label": "--label"},
        bool_flags=("no_length_check",),
    ),
    _T(
        "oligo_filter",
        "★ MPRA DESIGN FILTERS ★. Drops oligos that would break synthesis or "
        "confound the readout, and writes filter.log.tsv naming the reason "
        "for every drop. Sequence-level: homopolymer runs longer than "
        "--max-homopolymer-length, and any EcoRI (G^AATTC) or SbfI "
        "(CCTGCA^GG) site — those are the cloning sites, so an insert "
        "containing one cannot be cut back out. Coordinate-level (needs "
        "--regions): simple repeats covering more than "
        "--max-simple-repeat-fraction of the window, TSS overlap, and CTCF "
        "motif overlap. If a REF oligo fails its ALT partners are dropped "
        "too, since an ALT with no REF is not interpretable. Annotation BEDs "
        "are looked up in Data/References/OligoDesign/; a missing one "
        "disables just that filter with a warning.",
        {
            "type": "object",
            "properties": {
                "design":      {**_S_STRING, "description": "Designed oligo FASTA."},
                "regions":     {**_S_STRING, "description":
                                 "Region BED — enables repeat/TSS/CTCF filters."},
                "map":         {**_S_STRING},
                "variant_map": {**_S_STRING},
                "max_homopolymer_length":     {**_S_INTEGER, "default": 10},
                "max_simple_repeat_fraction": {"type": "number", "default": 0.25},
                "simple_repeats": {**_S_STRING},
                "tss_positions":  {**_S_STRING},
                "ctcf_motifs":    {**_S_STRING},
                "label":       {**_S_STRING},
            },
            "required": ["design"],
        },
        cli=["oligo", "filter"],
        flag_map={"design": "--design", "regions": "--regions",
                   "map": "--map", "variant_map": "--variant-map",
                   "max_homopolymer_length": "--max-homopolymer-length",
                   "max_simple_repeat_fraction": "--max-simple-repeat-fraction",
                   "simple_repeats": "--simple-repeats",
                   "tss_positions": "--tss-positions",
                   "ctcf_motifs": "--ctcf-motifs", "label": "--label"},
    ),
    _T(
        "oligo_pipeline",
        "★ MPRA OLIGO LIBRARY, END TO END ★. tile → design → filter → "
        "adapters in one run directory: takes regions (and optionally a VCF) "
        "plus a reference genome and returns a synthesis-ready design.fa.gz "
        "with ids namespaced '<sample>:<id>', every intermediate map, a "
        "per-drop filter log, and summary.json with counts for each stage. "
        "USE THIS to go from a candidate region/variant list to an order "
        "form in one call; use the individual tools only to inspect or "
        "re-run a single stage.",
        {
            "type": "object",
            "properties": {
                "regions":       {**_S_STRING},
                "variants":      {**_S_STRING, "description":
                                   "Optional VCF; omit for a region-only design."},
                "reference":     {**_S_STRING},
                "tile":          {**_S_BOOLEAN, "default": False,
                                   "description":
                                   "Tile regions first (skip if oligo-sized)."},
                "oligo_length":  {**_S_INTEGER, "default": 200},
                "min_overlap":   {**_S_INTEGER, "default": 50},
                "centering_max": {**_S_INTEGER},
                "variant_edge_exclusion": {**_S_INTEGER, "default": 20},
                "left_adapter":  {**_S_STRING},
                "right_adapter": {**_S_STRING},
                "max_homopolymer_length": {**_S_INTEGER, "default": 10},
                "use_all_regions": {**_S_BOOLEAN, "default": False},
                "remove_regions_without_variants": {**_S_BOOLEAN, "default": False},
                "label":         {**_S_STRING},
            },
            "required": ["regions", "reference"],
        },
        cli=["oligo", "pipeline"],
        flag_map={"regions": "--regions", "variants": "--variants",
                   "reference": "--reference", "oligo_length": "--oligo-length",
                   "min_overlap": "--min-overlap",
                   "centering_max": "--centering-max",
                   "variant_edge_exclusion": "--variant-edge-exclusion",
                   "left_adapter": "--left-adapter",
                   "right_adapter": "--right-adapter",
                   "max_homopolymer_length": "--max-homopolymer-length",
                   "label": "--label"},
        bool_flags=("tile", "use_all_regions",
                     "remove_regions_without_variants"),
    ),

    # ──────────────────────────────────────────────────────────────────
    # MPRA COUNT + ASSIGNMENT (what comes back off the sequencer, after
    # `oligo` designed the library). Python port of the analysis core of
    # kircherlab/MPRAsnakeflow (Max Schubach, BIH; MIT), validated
    # bit-identical against the upstream scripts.
    # ──────────────────────────────────────────────────────────────────
    _T(
        "mpraflow_assign_filter",
        "★ MPRA BARCODE→OLIGO ASSIGNMENT ★. Turns barcode-sorted "
        "'barcode<TAB>oligo<TAB>quality' triples (what any aligner "
        "backend emits) into a trusted assignment. A barcode is kept only "
        "when at least --minimum reads support it AND at least --fraction "
        "of them name the same oligo; barcodes whose reads disagree are "
        "reported as ambiguous and DROPPED, never guessed. Emits the "
        "assignment plus barcodes-per-oligo coverage stats — the headline "
        "assignment QC, since oligos with few barcodes give noisy ratios.",
        {
            "type": "object",
            "properties": {
                "pairs":    {**_S_STRING, "description":
                              "Barcode-sorted barcode/oligo/quality TSV."},
                "minimum":  {**_S_INTEGER, "default": 3},
                "fraction": {"type": "number", "default": 0.75,
                              "description": "Must exceed 0.5."},
                "report_other":     {**_S_BOOLEAN, "default": False},
                "report_ambiguous": {**_S_BOOLEAN, "default": False},
                "label":    {**_S_STRING},
            },
            "required": ["pairs"],
        },
        cli=["mpraflow", "assign-filter"],
        flag_map={"pairs": "--pairs", "minimum": "--minimum",
                   "fraction": "--fraction", "label": "--label"},
        bool_flags=("report_other", "report_ambiguous"),
    ),
    _T(
        "mpraflow_merge_counts",
        "★ MPRA PER-OLIGO ACTIVITY FROM BARCODE COUNTS ★. Joins per-barcode "
        "DNA/RNA counts to the assignment and produces the per-oligo "
        "activity table: counts summed over an oligo's barcodes, divided "
        "by the barcode count, scaled per million, then ratio = RNA/DNA "
        "and log2FoldChange. Optional outlier removal drops barcodes whose "
        "DNA/RNA ratio deviates from their oligo's median ('ratio_mad') or "
        "whose RNA count is extreme for their oligo ('rna_counts_zscore'). "
        "Reports how many barcodes matched the assignment — a low match "
        "rate means the assignment and the counts came from different "
        "libraries.",
        {
            "type": "object",
            "properties": {
                "counts":     {**_S_STRING, "description":
                                "barcode/dna_count/rna_count TSV."},
                "assignment": {**_S_STRING},
                "min_dna_counts": {**_S_INTEGER, "default": 0},
                "min_rna_counts": {**_S_INTEGER, "default": 1},
                "outlier_detection": {**_S_STRING, "enum":
                                       ["ratio_mad", "rna_counts_zscore"]},
                "label":      {**_S_STRING},
            },
            "required": ["counts", "assignment"],
        },
        cli=["mpraflow", "merge-counts"],
        flag_map={"counts": "--counts", "assignment": "--assignment",
                   "min_dna_counts": "--min-dna-counts",
                   "min_rna_counts": "--min-rna-counts",
                   "outlier_detection": "--outlier-detection",
                   "label": "--label"},
    ),
    _T(
        "mpraflow_variant_table",
        "★ MPRA ALLELIC SKEW ★. Given a per-oligo count table and a "
        "declaration naming which oligo is REF and which is ALT for each "
        "variant, computes log2FoldChange_expression = "
        "log2(ratio_ALT / ratio_REF) — how much the alternate allele "
        "changes reporter activity. This is the number an allelic MPRA "
        "exists to produce.",
        {
            "type": "object",
            "properties": {
                "counts":      {**_S_STRING},
                "declaration": {**_S_STRING, "description":
                                 "TSV with columns ID, REF, ALT."},
                "label":       {**_S_STRING},
            },
            "required": ["counts", "declaration"],
        },
        cli=["mpraflow", "variant-table"],
        flag_map={"counts": "--counts", "declaration": "--declaration",
                   "label": "--label"},
    ),
    # ──────────────────────────────────────────────────────────────────
    # Enhancer / regulatory-REGION annotation against the ENCODE SCREEN
    # cCRE registry. Distinct from ccre_* which is variant/point-level:
    # this is interval overlap over a user-supplied region list.
    # ──────────────────────────────────────────────────────────────────
    _T(
        "enhancer_annotate",
        "★ ENHANCER / REGULATORY-REGION cCRE ANNOTATION ★. Takes a "
        "user-supplied list of REGIONS (an enhancer library, CRISPRi "
        "target set, peak list) and annotates each against the ENCODE "
        "SCREEN candidate cis-regulatory element registry by INTERVAL "
        "OVERLAP — not point lookup. Per region returns: how many cCREs "
        "overlap, a representative class (ranked PLS > pELS > dELS > "
        "DNase-H3K4me3 > CA-CTCF > CA-TF > CTCF-only), every class seen, "
        "bases covered as a UNION (overlapping elements are not "
        "double-counted), coverage as a fraction of the region, and — for "
        "regions with no overlap — the nearest cCRE and its distance, "
        "which is what separates 'in a cCRE desert' from 'just missed "
        "one'. Accepts .xlsx/.xlsm (any sheet), BED, TSV, CSV, gzipped; "
        "coordinate columns are auto-detected and a headerless BED is "
        "recognised by shape. A CRISPR library repeats each enhancer once "
        "per sgRNA, so rows are deduplicated on coordinates by default. "
        "Requires the local registry: run enhancer_build_db once first.",
        {
            "type": "object",
            "properties": {
                "input":  {**_S_STRING, "description":
                            "Region list — xlsx / bed / tsv / csv (.gz ok)."},
                "sheet":  {**_S_STRING, "description":
                            "Worksheet name for xlsx input; see enhancer_inspect."},
                "registry": {**_S_STRING, "enum": ["V3","V4"], "default": "V3"},
                "chrom_col": {**_S_STRING, "description":
                               "Override chromosome column name."},
                "start_col": {**_S_STRING}, "end_col": {**_S_STRING},
                "name_col":  {**_S_STRING},
                "keep_duplicates": {**_S_BOOLEAN, "default": False,
                                     "description":
                                     "Annotate every row instead of unique coordinates."},
                "label":  {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["enhancer-annot", "annotate"],
        flag_map={"input": "--input", "sheet": "--sheet",
                   "registry": "--registry", "chrom_col": "--chrom-col",
                   "start_col": "--start-col", "end_col": "--end-col",
                   "name_col": "--name-col", "label": "--label"},
        bool_flags=("keep_duplicates",),
    ),
    _T(
        "enhancer_inspect",
        "★ INSPECT A REGION LIST BEFORE ANNOTATING ★. Shows an xlsx file's "
        "worksheets with their column headers and flags which sheets "
        "actually carry coordinate columns — a library workbook typically "
        "mixes region sheets with sgRNA-control and comparison sheets that "
        "have none. USE THIS FIRST when the user uploads a spreadsheet and "
        "the right sheet is not obvious, rather than guessing a sheet name.",
        {
            "type": "object",
            "properties": {"input": {**_S_STRING}},
            "required": ["input"],
        },
        cli=["enhancer-annot", "inspect"],
        flag_map={"input": "--input"},
    ),
    _T(
        "enhancer_build_db",
        "★ BUILD THE LOCAL cCRE REGISTRY DATABASE ★. Downloads the ENCODE "
        "SCREEN candidate cis-regulatory element registry (GRCh38, ~64 MB, "
        "~1.06 M elements) and indexes it into a local SQLite database, so "
        "region annotation afterwards runs entirely offline and is "
        "reproducible against a pinned registry version. Run once before "
        "enhancer_annotate; the download is cached.",
        {
            "type": "object",
            "properties": {
                "registry": {**_S_STRING, "enum": ["V3","V4"], "default": "V3"},
                "force": {**_S_BOOLEAN, "default": False,
                           "description": "Re-download even if cached."},
            },
        },
        cli=["enhancer-annot", "build-db"],
        flag_map={"registry": "--registry"},
        bool_flags=("force",),
    ),
    _T(
        "enhancer_db_stats",
        "Report what the local cCRE registry database contains — registry "
        "version, element count, build time, source URL and the most "
        "common class combinations. Use to confirm the database exists and "
        "which registry version an annotation would run against.",
        {"type": "object", "properties": {}},
        cli=["enhancer-annot", "db-stats"],
        flag_map={},
    ),
    _T(
        "mpralib_outliers",
        "★ MPRA BARCODE OUTLIER DETECTION ★. An oligo's activity is an "
        "average over its barcodes, so a few barcodes with wild RNA counts "
        "(a PCR jackpot, an unlucky integration site) can drag its ratio "
        "away from where the rest sit. Three detectors from Rosen et al. "
        "(2025), each catching a different failure: 'global' flags a "
        "barcode extreme for the whole replicate (|z|>3 on RNA counts); "
        "'oligo' flags one extreme for its OWN oligo (|z|>3 within the "
        "oligo); 'large_expression' flags activity running more than 5 "
        "log2-units ABOVE the oligo median (one-sided by design — runaway "
        "expression, not silence). NOTE the paper found removal changes "
        "variant calls very little (r>0.95), so report these as "
        "diagnostics rather than filtering by reflex.",
        {
            "type": "object",
            "properties": {
                "barcode_file": {**_S_STRING, "description":
                                  "IGVF 'reporter experiment barcode' TSV."},
                "method": {**_S_STRING,
                            "enum": ["global", "oligo", "large_expression"],
                            "default": "global"},
                "times_zscore":   {"type": "number", "default": 3.0},
                "times_activity": {"type": "number", "default": 5.0},
                "label":          {**_S_STRING},
            },
            "required": ["barcode_file"],
        },
        cli=["mpralib", "outliers"],
        flag_map={"barcode_file": "--barcode-file", "method": "--method",
                   "times_zscore": "--times-zscore",
                   "times_activity": "--times-activity", "label": "--label"},
    ),
    _T(
        "mpralib_validate",
        "★ IGVF MPRA FORMAT VALIDATION ★. Checks a file against one of the "
        "eight community file standards agreed by the IGVF MPRA focus "
        "group (Rosen et al. 2025, Supplementary Note S1): "
        "reporter_sequence_design, reporter_barcode_to_element_mapping, "
        "reporter_experiment_barcode, reporter_experiment, "
        "reporter_element, reporter_variant, reporter_genomic_element, "
        "reporter_genomic_variant. Reports the offending column and reason "
        "per bad row, not just pass/fail. USE THIS before submitting files "
        "to the IGVF portal, and to check that a file someone sent you is "
        "the format it claims to be. Every file `mpraflow` writes passes.",
        {
            "type": "object",
            "properties": {
                "file":   {**_S_STRING, "description": "TSV/BED (.gz ok)."},
                "schema": {**_S_STRING, "description":
                            "Standard name; see mpralib_schemas."},
                "max_errors": {**_S_INTEGER, "default": 20},
                "max_rows":   {**_S_INTEGER, "default": 0,
                                "description": "0 validates the whole file."},
                "label":  {**_S_STRING},
            },
            "required": ["file", "schema"],
        },
        cli=["mpralib", "validate"],
        flag_map={"file": "--file", "schema": "--schema",
                   "max_errors": "--max-errors", "max_rows": "--max-rows",
                   "label": "--label"},
    ),
    _T(
        "mpralib_schemas",
        "List the eight IGVF MPRA standard file formats with their "
        "required columns and which ones are headerless/positional. Use to "
        "pick the right --schema for mpralib_validate, or to answer 'what "
        "columns does an IGVF reporter experiment file need'.",
        {"type": "object", "properties": {}},
        cli=["mpralib", "schemas"],
        flag_map={},
    ),
    _T(
        "mpralib_consistency",
        "★ MPRA OUTLIER REPRODUCIBILITY ACROSS REPLICATES ★. Of the "
        "barcodes a replicate flags as outliers, what fraction is flagged "
        "in EVERY replicate. A detector firing on the same barcodes each "
        "time is describing the library; one firing on different barcodes "
        "each time is describing noise. Rosen et al. found episomal assays "
        "far more consistent (63-84%) than lentiviral ones (6-49%), and "
        "differentiated cells more consistent than progenitors. Reported "
        "as flagged-in-all / mean-flagged-per-replicate, the paper's "
        "definition.",
        {
            "type": "object",
            "properties": {
                "barcode_file": {**_S_STRING},
                "method": {**_S_STRING, "enum": ["global", "oligo"],
                            "default": "global"},
                "label":  {**_S_STRING},
            },
            "required": ["barcode_file"],
        },
        cli=["mpralib", "consistency"],
        flag_map={"barcode_file": "--barcode-file", "method": "--method",
                   "label": "--label"},
    ),
    _T(
        "mpraflow_complexity",
        "★ MPRA LIBRARY COMPLEXITY (LINCOLN-PETERSEN) ★. Answers 'would "
        "deeper sequencing help?'. Treats each replicate as a capture of "
        "the barcode pool and barcodes seen in two replicates as "
        "recaptures, giving a mark-recapture estimate of the TRUE library "
        "size; the gap between observed and estimated barcodes is the "
        "fraction of the library the sequencing missed. Reads the IGVF "
        "'reporter experiment barcode' wide file directly (one row per "
        "barcode, empty cells where a replicate did not see it) in a "
        "single streaming pass, so it scales to the 13M-barcode 240K "
        "libraries. Reproduces the figures in Rosen et al. (2025) exactly "
        "for 8K-neurons and 80K-neurons.",
        {
            "type": "object",
            "properties": {
                "barcode_file": {**_S_STRING, "description":
                                  "IGVF 'reporter experiment barcode' TSV."},
                "label":        {**_S_STRING},
            },
            "required": ["barcode_file"],
        },
        cli=["mpraflow", "complexity"],
        flag_map={"barcode_file": "--barcode-file", "label": "--label"},
    ),
    _T(
        "mpraflow_pipeline",
        "★ MPRA EXPERIMENT, END TO END ★. assignment + per-replicate "
        "barcode counts → per-oligo activity per replicate → master table "
        "filtered on barcodes-per-oligo → pooled per-oligo table → "
        "per-replicate and pooled allelic variant tables. Replicates are "
        "passed as repeated NAME=PATH pairs. USE THIS to take a finished "
        "MPRA experiment from counts to activity and allelic skew in one "
        "call; use the individual tools to inspect or re-run one stage.",
        {
            "type": "object",
            "properties": {
                "assignment":  {**_S_STRING},
                "replicate":   {**_S_ARRAY_S, "description":
                                 "Repeated NAME=PATH, e.g. rep1=rep1.tsv.gz."},
                "declaration": {**_S_STRING, "description":
                                 "Optional ID/REF/ALT TSV for allelic skew."},
                "labels":      {**_S_STRING},
                "threshold":   {**_S_INTEGER, "default": 10,
                                 "description":
                                 "Minimum barcodes per oligo per replicate."},
                "outlier_detection": {**_S_STRING, "enum":
                                       ["ratio_mad", "rna_counts_zscore"]},
                "min_dna_counts": {**_S_INTEGER, "default": 0},
                "min_rna_counts": {**_S_INTEGER, "default": 1},
                "label":       {**_S_STRING},
            },
            "required": ["assignment", "replicate"],
        },
        cli=["mpraflow", "pipeline"],
        flag_map={"assignment": "--assignment", "replicate": "--replicate",
                   "declaration": "--declaration", "labels": "--labels",
                   "threshold": "--threshold",
                   "outlier_detection": "--outlier-detection",
                   "min_dna_counts": "--min-dna-counts",
                   "min_rna_counts": "--min-rna-counts", "label": "--label"},
        flag_repeat=("replicate",),
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

    # ──────────────────────────────────────────────────────────────────
    # Spatial-ATAC-Hi-C — spatially resolved 3D genome + accessibility.
    # Ref: wangjuan001/Spatial-ATAC-Hi-C (MIT), Wang et al. Nat Methods
    # 2026, GSE307620. Clean-room; runHiC / Trim Galore are GPL-3.0 and
    # are external tools only, never runtime deps.
    # ──────────────────────────────────────────────────────────────────
    # ──────────────────────────────────────────────────────────────────
    # Human Transcription Factors database (Lambert 2018, Cell).
    # Local SQLite mirror; data downloaded at build time, never vendored.
    # ──────────────────────────────────────────────────────────────────
    # ──────────────────────────────────────────────────────────────────
    # Tabula Sapiens 2.0 human cell atlas (Cell 2026).
    # Retrieval + figure reproduction. Raw reads are DTA-gated; the open
    # routes are figshare, GEO GSE306755 and CELLxGENE.
    # ──────────────────────────────────────────────────────────────────
    _T(
        "tabula_status",
        "★ WHAT TABULA SAPIENS DATA IS LOCAL ★ and how it compares to the "
        "paper's own numbers (1,136,218 cells, 24 donors, 28 tissues, 182 "
        "fine / 38 broad cell types, 701 populations). Run this FIRST for "
        "any Tabula Sapiens question — it says whether the atlas is "
        "downloaded and whether the Human TF database is built, so you "
        "know which other commands can run.",
        {"type": "object", "properties": {}},
        cli=["tabula", "status"],
    ),

    _T(
        "tabula_pull_figshare",
        "Fetch the Tabula Sapiens 2.0 processed per-tissue h5ads from "
        "figshare (28 files, 57 GB total; each carries raw_counts, "
        "decontXcounts, log_normalized layers plus scVI embeddings and "
        "expert cell-type annotations). Pass `tissues` to fetch only what "
        "you need — most questions do not need all 57 GB. Resumable: "
        "already-complete files are skipped.",
        {
            "type": "object",
            "properties": {
                "tissues": {**_S_STRING, "description":
                             "Comma list e.g. 'Lung,Heart'; default all 28."},
            },
        },
        cli=["tabula", "pull-figshare"],
        flag_map={"tissues": "--tissues"},
        bool_flags=("list_only",),
    ),

    _T(
        "tabula_pull_geo",
        "Fetch the Tabula Sapiens GEO deposit (GSE306755). By default "
        "grabs the full 1.1M-cell metadata CSV (41 MB) — donor, tissue, "
        "cell type, age, sex per cell — which is all that Figure 1 and "
        "any demographic question needs, without the 57 GB of matrices.",
        {
            "type": "object",
            "properties": {
                "gse":      {**_S_STRING, "default": "GSE306755"},
                "download": {**_S_STRING, "description":
                              "Regex of files to fetch; default 'metadata'."},
            },
        },
        cli=["tabula", "pull-geo"],
        flag_map={"gse": "--gse", "download": "--download"},
        bool_flags=("list_only",),
    ),

    _T(
        "tabula_s3_manifest",
        "Enumerate and size the raw Tabula Sapiens AWS bucket WITHOUT "
        "downloading. IMPORTANT: the bucket is listable but every "
        "GetObject returns 403 — raw reads need a signed data transfer "
        "agreement for donor genetic privacy. Total is 103 TB (43 TB BAM, "
        "32 TB FASTQ, 28 TB STAR intermediates, 0.3 TB matrices). Use "
        "this to answer 'what raw data exists' and to plan a DTA request; "
        "never promise a raw-read download.",
        {
            "type": "object",
            "properties": {
                "donors": {**_S_STRING, "description":
                            "Comma list e.g. 'TSP1,TSP2'; default all 24."},
            },
        },
        cli=["tabula", "s3-manifest"],
        flag_map={"donors": "--donors", "label": "--label"},
    ),

    _T(
        "tabula_overview",
        "★ TABULA SAPIENS DATASET OVERVIEW (paper Figure 1) ★ — donor "
        "demographics (age range, sex split, age groups, ethnicity), "
        "cells and cell types per tissue, and every tissue x cell-type "
        "population. Answers 'what is in Tabula Sapiens', 'how many "
        "donors/tissues/cells', 'which donors contributed most organs'. "
        "Needs only the GEO metadata, not the 57 GB of matrices.",
        {
            "type": "object",
            "properties": {
                "metadata": {**_S_STRING, "description":
                              "Metadata CSV path; default the local GEO copy."},
                "label":    {**_S_STRING, "default": "overview"},
            },
        },
        cli=["tabula", "overview"],
        flag_map={"metadata": "--metadata", "label": "--label"},
        bool_flags=("no_figure",),
    ),

    _T(
        "tabula_donors",
        "Query Tabula Sapiens donors — the offline core of the upstream "
        "ChatTS tool. Filter by donor id, tissue, sex or age range and "
        "get each donor's age, sex, ethnicity, tissues contributed, cells "
        "and cell types. Use for 'which donors are over 60', 'who "
        "contributed lung', 'how many cells did TSP25 give'.",
        {
            "type": "object",
            "properties": {
                "donor":   {**_S_STRING, "description": "Comma list, e.g. TSP25."},
                "tissue":  {**_S_STRING},
                "sex":     {**_S_STRING},
                "min_age": {"type": "number"},
                "max_age": {"type": "number"},
                "limit":   {**_S_INTEGER, "default": 30},
            },
        },
        cli=["tabula", "donors"],
        flag_map={"donor": "--donor", "tissue": "--tissue", "sex": "--sex",
                   "min_age": "--min-age", "max_age": "--max-age",
                   "limit": "--limit", "label": "--label"},
    ),

    _T(
        "tabula_tf_specificity",
        "★ WHICH TRANSCRIPTION FACTORS ARE CELL-TYPE SPECIFIC ★ (paper "
        "Figure 2). Computes the tau specificity statistic for every "
        "human TF across cell types and splits specific (tau > 0.85) from "
        "ubiquitous, reporting each TF's peak cell type and DNA-binding "
        "domain. Answers 'is TF X cell-type specific', 'which TFs mark "
        "cell type Y', 'which TFs are universally expressed'. Needs the "
        "atlas and the Human TF database. Pass `matrix` to reuse a "
        "previously computed mean-expression .npz instead of recomputing.",
        {
            "type": "object",
            "properties": {
                "matrix":    {**_S_STRING, "description":
                               "Reuse a tf-matrix .npz."},
                "tissues":   {**_S_STRING},
                "threshold": {"type": "number", "default": 0.85},
                "label":     {**_S_STRING, "default": "tf_specificity"},
            },
        },
        cli=["tabula", "tf-specificity"],
        flag_map={"matrix": "--matrix", "tissues": "--tissues",
                   "threshold": "--threshold", "label": "--label"},
        bool_flags=("no_figure",),
    ),

    _T(
        "tabula_tf_enrichment",
        "GO enrichment of the NON-cell-type-specific transcription "
        "factors (paper Figure 3) — what the ubiquitous TFs do, grouped "
        "into gene expression/regulation, tissue maintenance, response to "
        "stimulus and metabolism. Takes the tf_tau.tsv written by "
        "tabula_tf_specificity.",
        {
            "type": "object",
            "properties": {
                "tau_table": {**_S_STRING, "description":
                               "tf_tau.tsv from tf-specificity."},
                "threshold": {"type": "number", "default": 0.85},
                "max_padj":  {"type": "number", "default": 0.02},
            },
            "required": ["tau_table"],
        },
        cli=["tabula", "tf-enrichment"],
        flag_map={"tau_table": "--tau-table", "threshold": "--threshold",
                   "max_padj": "--max-padj", "label": "--label"},
        bool_flags=("tf_background", "no_figure"),
    ),

    _T(
        "tabula_tf_regulons",
        "SCENIC-style TF regulon inference and activity scoring: "
        "co-expression (GRNBoost2-style), motif support, then AUCell "
        "activity per cell type. Answers 'is TF X actually ACTIVE here', "
        "as opposed to merely transcribed — pair it with "
        "tabula_tf_specificity, which measures expression only. "
        "Clean-room: pySCENIC is GPL-3 and is not imported.",
        {
            "type": "object",
            "properties": {
                "tissues":           {**_S_STRING},
                "group_key":         {**_S_STRING, "default": "cell_ontology_class"},
                "min_regulon_size":  {**_S_INTEGER, "default": 10},
                "label":             {**_S_STRING, "default": "tf_regulons"},
            },
        },
        cli=["tabula", "tf-regulons"],
        flag_map={"tissues": "--tissues", "group_key": "--group-key",
                   "min_regulon_size": "--min-regulon-size", "label": "--label"},
        bool_flags=("all_genes",),
    ),

    _T(
        "tabula_senescence",
        "★ SENESCENT-CELL BURDEN ACROSS HUMAN TISSUES ★ (paper Figure 4). "
        "Identifies CDKN2A+ MKI67- cells and reports their fraction by "
        "tissue, donor age and cell type. The paper finds ~48,114 such "
        "cells (4.4%), highest in eye/bladder/tongue and lowest in "
        "heart/muscle/ovary, rising modestly with donor age. Answers "
        "'which tissues accumulate senescent cells', 'does senescence "
        "increase with age', 'what fraction of cells are senescent'.",
        {
            "type": "object",
            "properties": {
                "tissues":       {**_S_STRING},
                "marker":        {**_S_STRING, "default": "CDKN2A"},
                "proliferation": {**_S_STRING, "default": "MKI67"},
                "label":         {**_S_STRING, "default": "senescence"},
            },
        },
        cli=["tabula", "senescence"],
        flag_map={"tissues": "--tissues", "marker": "--marker",
                   "proliferation": "--proliferation", "label": "--label"},
        bool_flags=("hallmarks", "no_figure"),
    ),

    _T(
        "tabula_sex_de",
        "Pseudobulk differential expression between sexes per tissue "
        "(paper Figure 5). Pseudobulks by donor first — cells within a "
        "donor are not independent replicates — then runs a "
        "negative-binomial test. Sex-specific organs (uterus, ovary, "
        "prostate, testis) are excluded by default since they carry no "
        "sex contrast. XIST and Y-linked genes are the positive control.",
        {
            "type": "object",
            "properties": {
                "tissues":     {**_S_STRING},
                "min_donors":  {**_S_INTEGER, "default": 2},
                "max_padj":    {"type": "number", "default": 0.05},
                "min_log2fc":  {"type": "number", "default": 1.0},
            },
        },
        cli=["tabula", "sex-de"],
        flag_map={"tissues": "--tissues", "min_donors": "--min-donors",
                   "max_padj": "--max-padj", "min_log2fc": "--min-log2fc",
                   "label": "--label"},
        bool_flags=("include_sex_specific",),
    ),

    _T(
        "humantfs_build_db",
        "Download the Human Transcription Factors database (Lambert 2018, "
        "humantfs.ccbr.utoronto.ca v1.01) and build the local SQLite "
        "mirror at Data/HumanTFs/human_tfs.sqlite. ~2 MB, seconds. Run "
        "this once before any TF question; every other humantfs tool "
        "needs it. Warns if the curated-TF count is not 1,639.",
        {
            "type": "object",
            "properties": {
                "with_pwms": {**_S_BOOLEAN, "description":
                               "Also fetch position weight matrices."},
            },
        },
        cli=["humantfs", "build-db"],
        bool_flags=("force", "with_pwms"),
    ),

    _T(
        "humantfs_is_tf",
        "★ IS THIS GENE A TRANSCRIPTION FACTOR? ★ Partitions a pasted "
        "gene list into curated TFs, non-TFs, and UNASSESSED (absent "
        "from the database entirely — which is not the same as 'not a "
        "TF'). Returns each TF's DNA-binding-domain family. Use this "
        "whenever someone asks which genes in a list are TFs, or wants "
        "to filter a marker/DE list down to regulators.",
        {
            "type": "object",
            "properties": {
                "genes": {**_S_STRING, "description":
                           "Gene symbols or Ensembl IDs, any separator."},
                "input": {**_S_STRING, "description": "File of genes."},
                "out":   {**_S_STRING, "description": "Write a TSV here."},
            },
        },
        cli=["humantfs", "is-tf"],
        flag_map={"genes": "--genes", "input": "--input", "out": "--out"},
    ),

    _T(
        "humantfs_lookup",
        "Everything the Human TF database knows about one or more genes: "
        "DNA-binding domain, TF assessment, binding mode, motif status, "
        "InterPro/PDB cross-references and the count of CIS-BP motifs. "
        "Use for a single gene in depth; use humantfs_is_tf to classify "
        "a whole list.",
        {
            "type": "object",
            "properties": {
                "genes": {**_S_STRING},
                "input": {**_S_STRING},
            },
        },
        cli=["humantfs", "lookup"],
        flag_map={"genes": "--genes", "input": "--input"},
    ),

    _T(
        "humantfs_list",
        "List curated TFs, optionally filtered by DNA-binding-domain "
        "family (e.g. 'C2H2 ZF', 'Homeodomain', 'bHLH', 'Forkhead'), "
        "binding mode, or whether a motif is known. Good for 'how many "
        "zinc-finger TFs are there' and for building a candidate set.",
        {
            "type": "object",
            "properties": {
                "dbd":          {**_S_STRING, "description": "DBD substring."},
                "binding_mode": {**_S_STRING},
                "limit":        {**_S_INTEGER, "default": 40},
                "out":          {**_S_STRING},
            },
        },
        cli=["humantfs", "list"],
        flag_map={"dbd": "--dbd", "binding_mode": "--binding-mode",
                   "limit": "--limit", "out": "--out"},
        bool_flags=("with_motif",),
    ),

    _T(
        "humantfs_families",
        "Census of DNA-binding-domain families across all 1,639 curated "
        "TFs, largest first. C2H2 zinc fingers dominate by a wide margin "
        "and many are computationally predicted with no validated motif.",
        {
            "type": "object",
            "properties": {
                "limit": {**_S_INTEGER, "default": 40},
                "out":   {**_S_STRING},
            },
        },
        cli=["humantfs", "families"],
        flag_map={"limit": "--limit", "out": "--out"},
    ),

    _T(
        "humantfs_motifs",
        "CIS-BP motif records for one TF: motif IDs, source, evidence "
        "type, and which are flagged best. Motif status is evidence that "
        "a motif EXISTS, not that the TF is active in any given cell "
        "type — pair it with regulon activity before claiming function.",
        {
            "type": "object",
            "properties": {
                "gene":  {**_S_STRING},
                "limit": {**_S_INTEGER, "default": 30},
            },
            "required": ["gene"],
        },
        cli=["humantfs", "motifs"],
        flag_map={"gene": "--gene", "limit": "--limit"},
    ),

    _T(
        "spatial_hic_pull_geo",
        "★ Fetch Spatial-ATAC-Hi-C data from GEO ★ — list (and optionally "
        "download) the supplementary deposits of a spatial Hi-C + ATAC "
        "series. Defaults to GSE307620, the Wang 2026 Nature Methods "
        "series with mouse brain, human cerebellum and glioblastoma "
        "samples. Use this FIRST when someone asks to work with, "
        "reproduce, or explore Spatial-ATAC-Hi-C data.",
        {
            "type": "object",
            "properties": {
                "gse":      {**_S_STRING, "default": "GSE307620",
                              "description": "GEO series accession."},
                "download": {**_S_STRING, "description":
                              "Regex; matching supplementary files are "
                              "downloaded, e.g. 'pairs|fragments|positions'."},
                "dest":     {**_S_STRING},
                "label":    {**_S_STRING, "default": "spatial_hic_geo"},
            },
        },
        cli=["spatial-hic", "pull-geo"],
        flag_map={"gse": "--gse", "download": "--download",
                   "dest": "--dest", "label": "--label"},
    ),

    _T(
        "spatial_hic_pixel_demux",
        "Split one barcoded Spatial-ATAC-Hi-C pairs file into the 50x50 "
        "spatial pixel grid (2,500 tissue pixels). Needs the two "
        "microfluidic barcode whitelists. Check `assigned_fraction` in "
        "the summary — a near-zero value means the barcode layout is the "
        "other way round, so retry with layout 'AB'.",
        {
            "type": "object",
            "properties": {
                "pairs":     {**_S_STRING, "description": "Barcoded .pairs[.gz]."},
                "barcode_a": {**_S_STRING, "description": "Barcode A whitelist."},
                "barcode_b": {**_S_STRING, "description": "Barcode B whitelist."},
                "layout":    {**_S_STRING, "default": "BA",
                               "description": "Concatenation order: BA or AB."},
                "barcode_field": {**_S_INTEGER, "description":
                               "1-based pairs column holding the barcode."},
                "max_pairs": {**_S_INTEGER, "description": "0 = all."},
                "label":     {**_S_STRING, "default": "spatial_hic_demux"},
            },
            "required": ["pairs", "barcode_a", "barcode_b"],
        },
        cli=["spatial-hic", "pixel-demux"],
        flag_map={"pairs": "--pairs", "barcode_a": "--barcode-a",
                   "barcode_b": "--barcode-b", "layout": "--layout",
                   "barcode_field": "--barcode-field",
                   "max_pairs": "--max-pairs", "label": "--label"},
    ),

    _T(
        "spatial_hic_qc",
        "Per-pixel Spatial-ATAC-Hi-C quality control: total / cis / trans "
        "and long-range (>=10 kb) contacts per tissue pixel, plus ArchR-"
        "style TSS enrichment when fragments and a gene model are given. "
        "The Wang 2026 paper reports medians of 25,343-58,403 total "
        "contacts, 88.1-90.3% cis and 24-33.3% long-range — use those as "
        "the reference band when judging a new dataset.",
        {
            "type": "object",
            "properties": {
                "pairs_dir":  {**_S_STRING, "description":
                                "A .pairs file or a directory of per-pixel ones."},
                "fragments":  {**_S_STRING, "description":
                                "fragments.tsv.gz, enables TSS enrichment."},
                "gene_model": {**_S_STRING, "description": "GTF/GFF or BED6."},
                "chrom_sizes": {**_S_STRING},
                "min_contacts": {**_S_INTEGER, "default": 1000},
                "label":      {**_S_STRING, "default": "spatial_hic_qc"},
            },
            "required": ["pairs_dir"],
        },
        cli=["spatial-hic", "qc"],
        flag_map={"pairs_dir": "--pairs-dir", "fragments": "--fragments",
                   "gene_model": "--gene-model", "chrom_sizes": "--chrom-sizes",
                   "min_contacts": "--min-contacts", "label": "--label"},
    ),

    _T(
        "spatial_hic_gas",
        "Gene activity score matrix from Spatial-ATAC-Hi-C ATAC "
        "fragments — Tn5 insertions over each gene's promoter (2 kb "
        "upstream of the TSS) plus gene body, SnapATAC2-style. Pass both "
        "barcode whitelists so the columns are AAxBB pixel ids; without "
        "them the columns are raw barcodes and will not join with the "
        "GAD matrix, the QC table or the spatial renderer.",
        {
            "type": "object",
            "properties": {
                "fragments":  {**_S_STRING, "description": "fragments.tsv.gz."},
                "gene_model": {**_S_STRING, "description": "GTF/GFF or BED6."},
                "barcode_a":  {**_S_STRING},
                "barcode_b":  {**_S_STRING},
                "layout":     {**_S_STRING, "default": "BA"},
                "upstream":   {**_S_INTEGER, "default": 2000},
                "min_cells":  {**_S_INTEGER, "default": 5},
                "label":      {**_S_STRING, "default": "spatial_hic_gas"},
            },
            "required": ["fragments", "gene_model"],
        },
        cli=["spatial-hic", "gas"],
        flag_map={"fragments": "--fragments", "gene_model": "--gene-model",
                   "barcode_a": "--barcode-a", "barcode_b": "--barcode-b",
                   "layout": "--layout", "upstream": "--upstream",
                   "min_cells": "--min-cells", "label": "--label"},
    ),

    _T(
        "spatial_hic_gad",
        "Gene-associated domain score matrix from Spatial-ATAC-Hi-C "
        "contacts — Hi-C pair ends counted over each gene body (scGAD, "
        "Shen 2022). This is the Hi-C counterpart of the ATAC gene "
        "activity score, and the paper clusters on it to recover the same "
        "spatial domains the ATAC modality gives.",
        {
            "type": "object",
            "properties": {
                "pairs_dir":  {**_S_STRING, "description":
                                "Directory of per-pixel .pairs files."},
                "gene_model": {**_S_STRING, "description": "GTF/GFF or BED6."},
                "min_cells":  {**_S_INTEGER, "default": 5},
                "label":      {**_S_STRING, "default": "spatial_hic_gad"},
            },
            "required": ["pairs_dir", "gene_model"],
        },
        cli=["spatial-hic", "gad"],
        flag_map={"pairs_dir": "--pairs-dir", "gene_model": "--gene-model",
                   "min_cells": "--min-cells", "label": "--label"},
    ),

    _T(
        "spatial_hic_compartment",
        "A/B compartment PC1 from Spatial-ATAC-Hi-C contacts, at 100 kb "
        "by default. Observed/expected, Pearson correlation, leading "
        "eigenvector. ALWAYS pass `gene_model` — without it the A/B sign "
        "is arbitrary and can flip between chromosomes.",
        {
            "type": "object",
            "properties": {
                "pairs_dir":  {**_S_STRING},
                "chrom_sizes": {**_S_STRING},
                "resolution": {**_S_INTEGER, "default": 100000},
                "chroms":     {**_S_STRING, "description": "Comma list."},
                "gene_model": {**_S_STRING, "description":
                                "Orients A toward gene-dense bins."},
                "label":      {**_S_STRING, "default": "spatial_hic_compartment"},
            },
            "required": ["pairs_dir"],
        },
        cli=["spatial-hic", "compartment"],
        flag_map={"pairs_dir": "--pairs-dir", "chrom_sizes": "--chrom-sizes",
                   "resolution": "--resolution", "chroms": "--chroms",
                   "gene_model": "--gene-model", "label": "--label"},
        bool_flags=("impute",),
    ),

    _T(
        "spatial_hic_cnv",
        "★ Spatial copy-number and tumour clone structure ★ from "
        "Spatial-ATAC-Hi-C contacts. Binned coverage scaled to a diploid "
        "baseline, at 5 Mb per pixel (the paper's setting for resolving "
        "clones) or finer for pseudobulk. Set `per_pixel` to get the "
        "pixel-by-bin matrix that reveals spatially segregated clones, "
        "and `smooth` for MAGIC diffusion before plotting. This is the "
        "right tool for questions about tumour heterogeneity, CNV, "
        "amplifications (EGFR / CDK4 / MDM2) or clonal structure in "
        "spatial Hi-C data.",
        {
            "type": "object",
            "properties": {
                "pairs_dir":  {**_S_STRING},
                "chrom_sizes": {**_S_STRING},
                "resolution": {**_S_INTEGER, "default": 5000000},
                "chroms":     {**_S_STRING},
                "gc":         {**_S_STRING, "description": "bedGraph of GC."},
                "ploidy":     {"type": "number", "default": 2.0},
                "label":      {**_S_STRING, "default": "spatial_hic_cnv"},
            },
            "required": ["pairs_dir"],
        },
        cli=["spatial-hic", "cnv"],
        flag_map={"pairs_dir": "--pairs-dir", "chrom_sizes": "--chrom-sizes",
                   "resolution": "--resolution", "chroms": "--chroms",
                   "gc": "--gc", "ploidy": "--ploidy", "label": "--label"},
        bool_flags=("per_pixel", "segment", "smooth"),
    ),

    _T(
        "spatial_hic_loops",
        "Quantify chromatin loops per spatial pixel and test which are "
        "cell-type / cluster specific (one-way ANOVA, p<0.05, as in the "
        "paper), with an optional APA pileup. Needs a BEDPE of anchors — "
        "this does NOT call loops de novo, so bring Peakachu output or "
        "another caller's. Pass `clusters` (pixel<TAB>cluster) to get the "
        "differential test.",
        {
            "type": "object",
            "properties": {
                "pairs_dir":  {**_S_STRING},
                "bedpe":      {**_S_STRING, "description": "Loop anchors."},
                "chrom_sizes": {**_S_STRING},
                "clusters":   {**_S_STRING, "description":
                                "pixel<TAB>cluster; enables the ANOVA."},
                "resolution": {**_S_INTEGER, "default": 10000},
                "alpha":      {"type": "number", "default": 0.05},
                "label":      {**_S_STRING, "default": "spatial_hic_loops"},
            },
            "required": ["pairs_dir", "bedpe"],
        },
        cli=["spatial-hic", "loops"],
        flag_map={"pairs_dir": "--pairs-dir", "bedpe": "--bedpe",
                   "chrom_sizes": "--chrom-sizes", "clusters": "--clusters",
                   "resolution": "--resolution", "alpha": "--alpha",
                   "label": "--label"},
        bool_flags=("apa",),
    ),

    _T(
        "spatial_hic_viz",
        "Render any per-pixel Spatial-ATAC-Hi-C value back into tissue "
        "space as a 50x50 heatmap (PNG + SVG). `column` accepts either a "
        "column of a per-pixel table (e.g. long_range_ratio from qc) or a "
        "gene name that is a ROW of a GAS/GAD score matrix — use this to "
        "show where a marker gene is active across the tissue.",
        {
            "type": "object",
            "properties": {
                "table":     {**_S_STRING, "description":
                               "A per-pixel TSV or a gene x pixel matrix."},
                "column":    {**_S_STRING, "description":
                               "Column name, or a gene row in a score matrix."},
                "positions": {**_S_STRING, "description":
                               "AtlasXBrowser tissue_positions CSV."},
                "cmap":      {**_S_STRING, "default": "magma"},
                "title":     {**_S_STRING},
                "label":     {**_S_STRING, "default": "spatial_hic_viz"},
            },
            "required": ["table", "column"],
        },
        cli=["spatial-hic", "viz"],
        flag_map={"table": "--table", "column": "--column",
                   "positions": "--positions", "cmap": "--cmap",
                   "title": "--title", "label": "--label"},
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
        flag_map={"limit": "--limit", "label": "--label",
                   "no_favor": "--no-call-favor"},
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
        flag_map={"limit": "--limit", "label": "--label",
                   "no_favor": "--no-call-favor"},
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
        flag_map={"limit": "--limit", "label": "--label",
                   "no_favor": "--no-call-favor"},
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
        "catalog_variant_enhancers",
        "★ ARE THERE ENHANCER-GENE PREDICTIONS OVERLAPPING THIS VARIANT ★ — "
        "the tool for 'which genes are predicted to be regulated by "
        "enhancers overlapping rs...', i.e. the ENCODE-rE2G / scE2G "
        "predictions shown under Gene Regulation on the Catalog's variant "
        "page. USE THIS rather than a variant-edge lookup: the Catalog has "
        "NO variant->genomic-element edge, so asking that way returns 0 and "
        "looks like an absence of evidence. The predictions hang off the "
        "ELEMENT and the link to the variant is positional overlap, so the "
        "query must be by region. On rs1250566 the variant-edge route "
        "reported 0 while this returns hundreds of ENCODE-rE2G rows "
        "targeting PPIF, LINC00595, ZCCHC24, SFTPA2 and KCNMA1. Writes a "
        "TSV of every prediction plus a per-gene score plot.",
        {
            "type": "object",
            "properties": {
                "variant": {**_S_STRING, "description":
                             "rsID (rs1250566) or SPDI "
                             "(NC_000010.11:79286695:G:A)."},
                "window":  {**_S_INTEGER, "default": 1,
                             "description":
                             "Overlap window in bp. 1 = enhancers that "
                             "actually contain the variant; widen only to "
                             "look at the neighbourhood, and say so if you "
                             "do."},
            },
            "required": ["variant"],
        },
        cli=["catalog", "variant-enhancers"],
        positional=("variant",),
        flag_map={"window": "--window"},
    ),
    _T(
        "catalog_variant_evidence",
        "★ WHICH VARIANTS HAVE EVIDENCE FROM MANY DIFFERENT ASSAYS ★ — "
        "ranks variants by how many DISTINCT assay types support them "
        "(eQTL, spliceQTL, caQTL, pQTL, GWAS, ADASTRA, MPRA, CRISPR, ...) "
        "and filters with `min_assays`. THIS IS THE RIGHT TOOL whenever "
        "someone asks for 'variants with >N assays', 'multi-assay "
        "variants', 'variants with the most functional evidence', or "
        "wants to prioritise variants by breadth of support — those are "
        "aggregations across every edge a variant has, which the "
        "per-variant and per-edge queries cannot express. Reports "
        "experimental assays separately from computational predictions "
        "(SEMVAR, cV2F) and curated resources (PharmGKB), because "
        "counting them together overstates experimental support. LD "
        "edges and the in-silico coding-variant predictor scores are "
        "excluded by default. IMPORTANT: there is no genome-wide mode — "
        "the Catalog rejects an unfiltered variant scan and the variants "
        "collection (~944 GB) is not mirrored. Supply `variants`, or "
        "point `from_traversal` at kg-traversal run directories to go "
        "from a gene panel to per-variant assay counts.",
        {
            "type": "object",
            "properties": {
                "variants":       {**_S_STRING, "description":
                                    "Comma/space separated variant ids "
                                    "(rsIDs or SPDI)."},
                "variants_file":  {**_S_STRING, "description":
                                    "File of variant ids, one per line."},
                "from_traversal": {**_S_STRING, "description":
                                    "Comma-separated kg-traversal run dirs "
                                    "(each holding evidence_pack.json). "
                                    "Offline; the gene-panel route."},
                "min_assays":     {**_S_INTEGER, "default": 3},
                "min_experimental": {**_S_INTEGER, "default": 0,
                                    "description":
                                    "Also require this many distinct "
                                    "EXPERIMENTAL assays. Set it when the "
                                    "asker means wet-lab assays (MPRA, "
                                    "CRISPR, eQTL) and should not be given "
                                    "variants padded over the threshold by "
                                    "predictors or curated resources."},
                "include_ld":     {**_S_BOOLEAN, "description":
                                    "Count 'linkage disequilibrum' as an "
                                    "assay. Off by default: it attaches to "
                                    "nearly every common variant."},
                "include_predictions": {**_S_BOOLEAN, "description":
                                    "Count the in-silico coding-variant "
                                    "predictors (SIFT/PolyPhen2/CADD/...) "
                                    "as assays. Off by default."},
                "label":          {**_S_STRING, "default": "run"},
            },
        },
        cli=["catalog", "variant-evidence"],
        flag_map={"variants": "--variants",
                   "variants_file": "--variants-file",
                   "from_traversal": "--from-traversal",
                   "min_assays": "--min-assays",
                   "min_experimental": "--min-experimental",
                   "include_ld": "--include-ld",
                   "include_predictions": "--include-predictions",
                   "label": "--label"},
        bool_flags=("include_ld", "include_predictions"),
    ),

    _T(
        "mct_analyze",
        "★ ANALYSE snMCT-seq / snm3C-seq — BOTH MODALITIES ★ — these assays "
        "measure methylation AND the transcriptome from the SAME nucleus, "
        "published as separate matrices. The generic single-cell route "
        "analyses whichever matrix sorts first and silently drops the "
        "other, so it reports an RNA clustering and never mentions the "
        "methylation the assay exists to measure. This analyses each half "
        "with the normalisation it needs -- RNA counts are log-normalised, "
        "methylation values are ALREADY a ratio centred on 1 and must not "
        "be -- then cross-tabulates the two clusterings on the cells they "
        "share and reports an adjusted Rand index. It also interrogates the "
        "methylation clusters against the per-cell QC report and says "
        "whether they track a methylation level (biology) or read depth "
        "(library artefact).",
        {
            "type": "object",
            "properties": {
                "accession":  {**_S_STRING, "description":
                                "snMCT-seq MeasurementSet (IGVFDS...)."},
                "resolution": {**_S_NUMBER, "default": 1.0,
                                "description": "Leiden resolution for the "
                                               "methylation clustering."},
                "n_pcs":      {**_S_INTEGER, "default": 30},
                "skip_rna":   {**_S_BOOLEAN, "description":
                                "Methylation only."},
                "label":      {**_S_STRING},
            },
            "required": ["accession"],
        },
        cli=["mct", "analyze"],
        positional=("accession",),
        flag_map={"resolution": "--resolution", "n_pcs": "--n-pcs",
                   "label": "--label"},
        bool_flags=("skip_rna",),
    ),
    _T(
        "mct_discover",
        "★ WHICH snMCT-seq MODALITIES ARE PUBLISHED ★ — lists the "
        "cell-by-gene, cell-by-bin methylation, cell-by-position "
        "methylation and per-cell QC files for an snMCT-seq dataset, with "
        "sizes, and says plainly when a half is missing. Call it before "
        "mct_analyze, or to explain what a dual-modality dataset contains.",
        {
            "type": "object",
            "properties": {"accession": {**_S_STRING}},
            "required": ["accession"],
        },
        cli=["mct", "discover"],
        positional=("accession",),
    ),
    _T(
        "crispr_screen_analyze",
        "★ REPROCESS A WHOLE CRISPR FACS SCREEN FROM RAW READS ★ — NOT for "
        "a BASE-EDITING screen: if the library's guide names carry ABE or "
        "CBE, use base_editing_screen_analyze instead, because exact guide "
        "matching discards the self-edited reads (36.7% vs 62.5% assigned on "
        "IGVFDS6464SOVZ). This tool detects that and refuses rather than "
        "undercounting. Otherwise: give it "
        "ANY sorted-bin accession and it finds the screen's other bins and "
        "replicates, counts constructs in the libraries the comparison "
        "needs, compares the low tail against the high tail per replicate, "
        "aggregates constructs onto the variant each installs, and reports "
        "per-variant effects with FDR plus a volcano and "
        "replicate-agreement plot. USE THIS rather than analysing one "
        "accession: a single bin is one tail of one replicate and the "
        "measurement IS the comparison between bins, so counting one "
        "library alone yields composition and no biology. It handles both "
        "bin naming conventions in use (bottom20/top20 and Bot20/Top20) and "
        "the unsorted Bulk bins some screens add. On IGVFDS6464SOVZ it "
        "pulls the 18loci_uptake screen (4 replicates x 4 bins, 0.39 GB) "
        "and scores 1,656 targets; on IGVFDS5542IBUS it finds all 20 "
        "libraries of the LDLR137-219 screen. It also chooses the counting "
        "key by testing candidates against real reads, so a prime-editing "
        "library whose 1,741 pegRNAs share 52 spacers is counted by RT "
        "template rather than collapsed onto the shared spacers -- the "
        "output states which key was used and what fraction of reads it "
        "assigned. Positive score = enriched in the LOW tail = the variant "
        "reduces the sorted phenotype.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING, "description":
                               "Any sorted-bin MeasurementSet of the screen."},
                "tail":      {**_S_INTEGER, "default": 20,
                               "description":
                               "Which tail pair to compare, e.g. 20 for "
                               "bottom20% vs top20%."},
                "min_count": {**_S_INTEGER, "default": 10},
                "max_reads": {**_S_INTEGER},
                "all_bins":  {**_S_BOOLEAN, "description":
                               "Also count bins the tail comparison does "
                               "not use (the other tail, and Bulk). Slower; "
                               "for QC of library complexity."},
                "label":     {**_S_STRING},
            },
            "required": ["accession"],
        },
        cli=["crispr-screen", "analyze"],
        positional=("accession",),
        bool_flags={"all_bins"},
        flag_map={"tail": "--tail", "min_count": "--min-count",
                   "max_reads": "--max-reads", "label": "--label"},
    ),
    _T(
        "crispr_screen_discover",
        "★ WHAT ELSE BELONGS TO THIS CRISPR SCREEN ★ — lists every sorted "
        "bin and replicate of the screen a given accession belongs to, read "
        "from the submitter alias (e.g. 18loci_uptake_Rep1_bottom20). Call "
        "it to show a user why one accession is not analysable alone, or to "
        "check the screen was identified correctly before spending compute.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING},
            },
            "required": ["accession"],
        },
        cli=["crispr-screen", "discover"],
        positional=("accession",),
    ),
    _T(
        "base_editing_screen_analyze",
        "★ ANALYSE A BASE-EDITING (ABE/CBE) SCREEN THE WAY BEAN DOES ★ — use "
        "this INSTEAD OF crispr_screen_analyze whenever the library is a base "
        "editor, because a base editor edits the guide's own locus too and the "
        "protospacer sequenced back carries A>G (ABE) or C>T (CBE) changes of "
        "its own. Exact matching throws those reads away: measured on "
        "IGVFDS6464SOVZ, an 8,192-guide ABE screen, exact matching assigns "
        "36.7% of reads and BEAN-style masked matching assigns 62.5%. It "
        "follows crispr-bean's method (mask the edited base on both sides "
        "before comparing, rather than allowing free mismatches), detects the "
        "editor from the library's guide names AND by measuring which masking "
        "actually recovers reads, and estimates per-guide EDITING ACTIVITY "
        "from self-editing so a weakly-editing guide is reported as "
        "underpowered rather than as having no effect. It reports what it "
        "does NOT do: BEAN's Bayesian variant/tiling model, reporter-allele "
        "and bystander analysis (the IGVF library's reporter column is "
        "empty), and the bcmatch/semimatch split (no guide barcode is "
        "published, so masked collisions stay ambiguous instead of being "
        "assigned). **To actually run crispr-bean, pass run_bean=true** -- the "
        "real `bean` binary is installed and this tool drives it end to end "
        "(create-screen + run sorting variant) and plots its posteriors. "
        "Without run_bean it only does the frequentist per-guide scoring, so "
        "a request for BEAN, a Bayesian model, posteriors or credible "
        "intervals REQUIRES run_bean=true. Note BEAN's activity-normalised "
        "MixtureNormal model is unavailable on IGVF data (it needs an "
        "X_bcmatch layer built from a guide barcode IGVF does not publish), "
        "so the posteriors come from BEAN's Normal model and the output "
        "labels them as not activity-normalised.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING, "description":
                               "Any sorted-bin MeasurementSet of the screen."},
                "editor":    {**_S_STRING, "description":
                               "ABE or CBE. Omit to detect it from the "
                               "library names and the reads."},
                "tail":      {**_S_INTEGER, "default": 20},
                "min_count": {**_S_INTEGER, "default": 10},
                "max_reads": {**_S_INTEGER},
                "label":     {**_S_STRING},
                "run_bean":  {**_S_BOOLEAN, "description":
                               "★ SET THIS WHENEVER THE USER ASKS FOR "
                               "CRISPR-BEAN ★ Hands the base-edit-aware counts "
                               "to the REAL `bean` binary: writes BEAN's four "
                               "input tables, runs `bean create-screen`, then "
                               "`bean run sorting variant`, and reports its "
                               "per-target Bayesian posteriors (mu, mu_sd, "
                               "mu_z, n_guides) plus a bean_posteriors.png. "
                               "Without this flag only the per-guide "
                               "frequentist scoring runs and NO BEAN output "
                               "is produced."},
                "bean_iter": {**_S_INTEGER, "description":
                               "Override BEAN's --n-iter. 0 = BEAN's default."},
                "bean_mode": {**_S_STRING, "description":
                               "variant (default) or tiling."},
            },
            "required": ["accession"],
        },
        cli=["bean", "analyze"],
        positional=("accession",),
        flag_map={"editor": "--editor", "tail": "--tail",
                   "min_count": "--min-count", "max_reads": "--max-reads",
                   "label": "--label", "bean_iter": "--bean-iter",
                   "bean_mode": "--bean-mode"},
        bool_flags=("run_bean",),
    ),
    _T(
        "crispr_bean_analyze",
        "★ RUN THE REAL CRISPR-BEAN BAYESIAN MODEL ★ — THIS is the tool for "
        "any request naming crispr-bean, BEAN, a Bayesian variant-effect "
        "model, posterior distributions or credible intervals on a "
        "base-editing screen. It does the base-edit-aware guide assignment, "
        "writes BEAN's four input tables, then invokes the real `bean` "
        "binary: `bean create-screen` followed by `bean run sorting "
        "variant`. It returns BEAN's per-target posteriors (mu, mu_sd, "
        "mu_z, n_guides), the bean_element_result CSV, the bean_screen.h5ad, "
        "and a bean_posteriors.png. Do NOT use base_editing_screen_analyze "
        "for a BEAN request -- that runs only the frequentist per-guide "
        "scoring and produces no BEAN output. Expect the output to say "
        "'Normal model, NOT activity-normalised': BEAN's activity-normalised "
        "MixtureNormal needs an X_bcmatch layer built from a guide barcode "
        "IGVF does not publish, so the tool falls back to BEAN's Normal "
        "model and labels it. That is correct behaviour, not a failure.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING, "description":
                               "Any sorted-bin MeasurementSet of the screen, "
                               "e.g. IGVFDS6464SOVZ."},
                "editor":    {**_S_STRING, "description":
                               "ABE or CBE. Omit to detect it."},
                "tail":      {**_S_INTEGER, "default": 20},
                "max_reads": {**_S_INTEGER, "description":
                               "Reads per library. 400000 is a good default; "
                               "counts are cached per library."},
                "bean_iter": {**_S_INTEGER},
                "bean_mode": {**_S_STRING, "description":
                               "variant (default) or tiling."},
                "label":     {**_S_STRING},
            },
            "required": ["accession"],
        },
        cli=["bean", "bayesian"],
        positional=("accession",),
        flag_map={"editor": "--editor", "tail": "--tail",
                   "max_reads": "--max-reads", "bean_iter": "--bean-iter",
                   "bean_mode": "--bean-mode", "label": "--label"},
    ),
    _T(
        "bean_paper_benchmark",
        "★ BENCHMARK IGVFAGENT AGAINST THE CRISPR-BEAN PAPER ★ — Ryu et al., "
        "Nat Genet 56:925-937 (2024). Runs on the AUTHORS' OWN deposited "
        "screen objects (Zenodo 10.5281/zenodo.10139794), which are already "
        "downloaded, and compares 18 published claims to what we measure. "
        "**Use this, not an IGVF screen, for any question about reproducing "
        "or benchmarking against that paper.** The IGVF-deposited Sherwood "
        "screens are NOT the paper's screens: IGVFDS6464SOVZ has 8,192 "
        "guides and no labelled non-targeting controls, the paper's LDL-C "
        "GWAS library has 3,455 and 100, and only LDLR and HNF4A overlap. "
        "CRITICALLY: the paper's deposit DOES carry the reporter, the guide "
        "barcode (layer X_bcmatch) and per-guide edit rates that IGVF does "
        "not publish — so BEAN's full activity-normalised MixtureNormal and "
        "accessibility models CAN be fitted here. Any statement that those "
        "are untestable applies to IGVF data only, and is wrong about this "
        "benchmark. Raw SRA data (PRJNA1042659) is NOT needed. Subcommands: "
        "`describe` (what the deposit holds — run first), `measure` "
        "(library composition, replicate agreement, editing rates), `run` "
        "(fit BEAN / BEAN-Reporter / BEAN-Uniform and score AUPRC), "
        "`report` (measured vs published, claim by claim).",
        {
            "type": "object",
            "properties": {
                "subcommand": {**_S_STRING, "description":
                                "describe | measure | run | report | fetch"},
                "screen": {**_S_STRING, "description":
                            "ldlvar (LDL-C GWAS library) or ldlrcds (LDLR CDS "
                            "tiling library). Required for describe/run."},
                "model": {**_S_STRING, "description":
                           "For `run`: bean (MixtureNormal + accessibility), "
                           "reporter, uniform, or all."},
            },
            "required": ["subcommand"],
        },
        cli=["bean-benchmark"],
        positional=("subcommand", "screen"),
        flag_map={"model": "--model"},
    ),
    # ── IGVF Portal submission ────────────────────────────────────────
    # These exist because the submission loop is monthly: submit, wait for
    # DACC audits, learn at the next meeting what was missing. The same
    # properties are missing every time, so each rule cites the meeting that
    # produced it and the checks run in seconds against the live Portal.
    _T(
        "igvf_submit_audit",
        "★ WILL THIS SUBMISSION PASS THE DACC AUDITS? ★ — run the recurring-"
        "findings checklist against an object ALREADY ON the IGVF Portal, "
        "before the monthly submission meeting does it for you. Catches the "
        "properties that are missing every month: `description` on a file "
        "set, `input_file_sets` on a PredictionSet, `derived_from` on a "
        "derived file, `reference_files`, `file_format_specifications`, "
        "`analysis_step_version`, and inputs whose status has become "
        "revoked/archived/replaced (which silently blocks release of "
        "everything downstream, and the Portal names the replacement). Use "
        "--recurse to include a file set's member files. Exits non-zero when "
        "release is blocked. Call this whenever the user asks whether their "
        "submission is ready, why an audit is failing, what the DACC "
        "flagged, or what is missing from a prediction set / curated set / "
        "tabular file.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING, "description":
                              "IGVF accession, e.g. IGVFDS2581EDPS or "
                              "IGVFFI0856UJHP."},
                "recurse": {**_S_BOOLEAN, "description":
                            "Also audit the file set's member files."},
                "mode": {**_S_STRING, "description":
                         "prod (default) or staging."},
            },
            "required": ["accession"],
        },
        cli=["submit", "audit"],
        positional=("accession",),
        flag_map={"mode": "-m"},
        bool_flags=("recurse",),
    ),
    _T(
        "igvf_submit_inputs",
        "★ ARE ANY OF MY INPUTS DEAD? ★ — lists every input linked to a "
        "Portal object with its status, flagging revoked, archived, replaced "
        "and deleted ones and naming the replacement the Portal points to. "
        "This is the check for 'my file was fine last month and now release "
        "is blocked': an input was corrected upstream. Use it before "
        "igvf_submit_crosscheck, which decides whether the correction "
        "actually touched your data.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING},
                "mode": {**_S_STRING},
            },
            "required": ["accession"],
        },
        cli=["submit", "inputs"],
        positional=("accession",),
        flag_map={"mode": "-m"},
    ),
    _T(
        "igvf_submit_crosscheck",
        "★ AN INPUT WAS REVOKED AND CORRECTED — DO I HAVE TO REUPLOAD? ★ — "
        "compares the old and corrected input against YOUR derived file and "
        "answers the only question that matters: did any of the revised "
        "variants actually reach your file? If none did, repoint "
        "`derived_from` and release; if some did, the file needs regenerating. "
        "Keys on (chrom, pos) -> {(ref, alt)}, because a reference-allele "
        "mismatch is corrected IN PLACE -- the coordinate does not move while "
        "ref/alt change, so comparing positions alone reports exactly this "
        "case as 'no change'. Refuses to give a verdict unless unrevised "
        "input variants are also found in your file, since otherwise the two "
        "files' conventions differ and a clean result would be false.",
        {
            "type": "object",
            "properties": {
                "old": {**_S_STRING, "description":
                        "The revoked input accession, or a local path."},
                "new": {**_S_STRING, "description":
                        "The corrected input accession, or a local path."},
                "mine": {**_S_STRING, "description":
                         "Your derived file, accession or local path."},
                "mode": {**_S_STRING},
            },
            "required": ["old", "new", "mine"],
        },
        cli=["submit", "crosscheck"],
        flag_map={"old": "--old", "new": "--new", "mine": "--mine",
                   "mode": "-m"},
    ),
    _T(
        "igvf_submit_validate",
        "★ CHECK A FILE BEFORE UPLOADING IT ★ — local, no Portal call: gzip "
        "(the DACC asks for gzipped files every month), ragged rows, md5 and "
        "content md5, emptiness, and with --bed3 that a BED file is 0-based "
        "half-open rather than 1-based. Run this before upload, not after the "
        "audit.",
        {
            "type": "object",
            "properties": {
                "path": {**_S_STRING, "description": "Local file to check."},
                "bed3": {**_S_BOOLEAN, "description":
                         "Also apply BED coordinate checks."},
            },
            "required": ["path"],
        },
        cli=["submit", "validate"],
        positional=("path",),
        bool_flags=("bed3",),
    ),
    _T(
        "igvf_submit_fileset_files",
        "Member file accessions of a file set, as a paste-ready JSON array — "
        "the list the DACC repeatedly asks submitters to collect by hand and "
        "paste into `derived_from` or `input_file_sets` in the Portal's JSON "
        "editing view.",
        {
            "type": "object",
            "properties": {"accession": {**_S_STRING},
                            "mode": {**_S_STRING}},
            "required": ["accession"],
        },
        cli=["submit", "fileset-files"],
        positional=("accession",),
        flag_map={"mode": "-m"},
    ),
    _T(
        "igvf_submit_template",
        "The required properties of a Portal profile as a TSV header, and "
        "`enum` for a property's allowed values — both read from the LIVE "
        "schema at api.data.igvf.org/profiles, so they cannot drift from what "
        "the Portal will accept. Use when the user asks what fields a "
        "prediction set / curated set / tabular file needs, or what values "
        "`file_set_type` or `content_type` allows.",
        {
            "type": "object",
            "properties": {
                "profile": {**_S_STRING, "description":
                            "e.g. prediction_set, curated_set, tabular_file."},
                "mode": {**_S_STRING},
            },
            "required": ["profile"],
        },
        cli=["submit", "template"],
        positional=("profile",),
        flag_map={"mode": "-m"},
    ),
    _T(
        "igvf_portal_qc_audits",
        "★ WHAT IS WRONG ACROSS THE PORTAL, OR ACROSS MY LAB'S DATA? ★ — "
        "reads the Portal's OWN audit facets (audit.ERROR.category, "
        "NOT_COMPLIANT, WARNING, INTERNAL_ACTION) so a complete picture "
        "arrives in one request instead of 150,000. Scope to a lab or an "
        "item type. Use for 'what does my lab still need to fix', not for a "
        "single object -- that is igvf_submit_audit.",
        {
            "type": "object",
            "properties": {
                "item_type": {**_S_STRING, "description":
                              "e.g. TabularFile, PredictionSet, File."},
                "lab": {**_S_STRING, "description": "Lab name to scope to."},
                "mode": {**_S_STRING},
            },
        },
        cli=["portal-qc", "audits"],
        flag_map={"item_type": "--type", "lab": "--lab", "mode": "-m"},
    ),
    _T(
        "igvf_portal_qc_provenance",
        "Provenance integrity across the Portal: files deriving from inputs "
        "that are revoked or archived, with the replacement named. The Portal "
        "audits most of this at INTERNAL_ACTION severity, mixed in with "
        "thousands of other status mismatches; and it has a measured blind "
        "spot -- of 55 live files derived from an unusable file, 30 carried no "
        "audit at all, every one of them status=released deriving from "
        "status=archived. --dead-status picks the policy.",
        {
            "type": "object",
            "properties": {
                "lab": {**_S_STRING},
                "dead_status": {**_S_STRING, "description":
                                "Comma-separated, e.g. revoked,archived."},
                "mode": {**_S_STRING},
            },
        },
        cli=["portal-qc", "provenance"],
        flag_map={"lab": "--lab", "dead_status": "--dead-status",
                   "mode": "-m"},
    ),
    _T(
        "base_editing_screen_discover",
        "★ IS THIS A BASE-EDITING SCREEN, AND WHAT DOES BEAN NEED THAT IGVF "
        "PUBLISHES? ★ — names the screen and its bins, the guide library and "
        "its control classes, which base editor the library's guide names "
        "state, and which BEAN inputs are actually present in the IGVF "
        "library table (guide barcode, reporter allele). Call it before "
        "analysing a CRISPR screen whose guides mention ABE or CBE, to show "
        "why exact guide matching would undercount, and it prints the real "
        "BEAN pipeline command for the screen.",
        {
            "type": "object",
            "properties": {"accession": {**_S_STRING}},
            "required": ["accession"],
        },
        cli=["bean", "discover"],
        positional=("accession",),
    ),
    _T(
        "gradient_screen_analyze",
        "★ ANALYSE A FACS SCREEN SORTED INTO LETTERED EXPRESSION BINS ★ — "
        "for a screen whose bins are BinA..BinF rather than a bottom/top "
        "tail pair. 554 MeasurementSets on the Portal have this shape and no "
        "other tool fits them: crispr_screen_analyze compares two tails, "
        "while this computes each construct's frequency-weighted MEAN BIN "
        "across all six, which is the construct's centre of mass along the "
        "expression axis. Give it ANY bin of the screen and it finds the "
        "rest through the shared construct library plus the alias series, "
        "counts every bin, centres each sort on its non-targeting controls "
        "so re-drawn gates do not read as biology, and reports per-construct "
        "delta-mean-bin with FDR plus QC and plots. Covers both families: "
        "Variant-EFFECTS / prime-editing read out by allelic sequencing "
        "(e.g. IGVFDS3899ANMJ, PPIF promoter, 66 sets over 11 sorts) and "
        "CRISPR FlowFISH read out by gRNA sequencing (e.g. IGVFDS8710ZSOZ, "
        "KITLG, 24 sets over 4 flow replicates). Positive delta_bins = "
        "sorted into HIGHER expression bins. USE THIS rather than "
        "raw_pipeline_run, which correctly refuses these reads, and rather "
        "than sge_analyze, which wants an SGE editing-template design these "
        "libraries do not publish.",
        {
            "type": "object",
            "properties": {
                "accession":  {**_S_STRING, "description":
                                "Any lettered-bin MeasurementSet of the screen."},
                "min_count":  {**_S_INTEGER, "default": 20,
                                "description":
                                "Minimum total reads for a construct in a sort."},
                "max_reads":  {**_S_INTEGER, "description":
                                "Cap reads scanned per library."},
                "max_sorts":  {**_S_INTEGER, "description":
                                "Score only the first N sorts, for a quick look."},
                "label":      {**_S_STRING},
            },
            "required": ["accession"],
        },
        cli=["gradient-screen", "analyze"],
        positional=("accession",),
        flag_map={"min_count": "--min-count", "max_reads": "--max-reads",
                   "max_sorts": "--max-sorts", "label": "--label"},
    ),
    _T(
        "gradient_screen_discover",
        "★ WHAT ELSE BELONGS TO THIS LETTERED-BIN SCREEN ★ — lists every bin "
        "and sort of the screen a given accession belongs to, reading the "
        "submitter alias (e.g. CRUDO_KITLG-Auxin6hrs-FF2-BinC, or "
        "PPIF_promoter-BioRep3-FFrep3-BinF). Call it to show a user why one "
        "bin is not analysable alone, or to check the screen and its "
        "replicate structure before spending compute.",
        {
            "type": "object",
            "properties": {"accession": {**_S_STRING}},
            "required": ["accession"],
        },
        cli=["gradient-screen", "discover"],
        positional=("accession",),
    ),
    _T(
        "sge_analyze",
        "★ ANALYSE SATURATION GENOME EDITING (SGE) DATA PROPERLY ★ — the "
        "tool for an SGE / saturation-mutagenesis MeasurementSet, which "
        "raw_pipeline REFUSES for good reason: SGE reads are a fixed "
        "amplicon, not a transcript library, so quantifying them against a "
        "transcriptome only restates which gene the amplicon covers. What "
        "SGE measures is the fate of each programmed variant -- damaging "
        "variants are DEPLETED from the population over time -- so this "
        "reads the editing-template design, calls variants against the "
        "reference amplicon, finds the matching EARLY timepoint on the same "
        "target and replicate, and reports log2(late/early) per variant "
        "with plots. Negative score = depleted = damaging. Give it the LATE "
        "(e.g. day12) accession; the early mate is discovered automatically, "
        "or pass `early` explicitly.",
        {
            "type": "object",
            "properties": {
                "accession":  {**_S_STRING, "description":
                                "SGE MeasurementSet, the LATE timepoint "
                                "(IGVFDS...)."},
                "early":      {**_S_STRING, "description":
                                "Early-timepoint accession. Auto-discovered "
                                "from the sample alias if omitted."},
                "target":     {**_S_STRING, "description":
                                "Design target, e.g. PALB2_X7A. Inferred "
                                "from the set alias if omitted."},
                "max_reads":  {**_S_INTEGER, "description":
                                "Cap reads per FASTQ for a quick look; omit "
                                "to use all of them."},
                "min_count":  {**_S_INTEGER, "default": 5},
                "label":      {**_S_STRING},
            },
            "required": ["accession"],
        },
        cli=["sge", "analyze"],
        positional=("accession",),
        flag_map={"early": "--early", "target": "--target",
                   "max_reads": "--max-reads", "min_count": "--min-count",
                   "label": "--label"},
    ),
    _T(
        "sge_design",
        "★ WHAT DOES AN SGE LIBRARY TARGET ★ — lists the editing-template "
        "targets behind an SGE MeasurementSet: each amplicon's genomic "
        "coordinates, its length, and the fixed 'required edits' every "
        "template carries. Use it to see which exon tile a dataset covers, "
        "or to pick a `target` for sge_analyze when a library spans many "
        "(PALB2's spans 37).",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING, "description":
                               "SGE MeasurementSet or construct library set."},
            },
            "required": ["accession"],
        },
        cli=["sge", "design"],
        positional=("accession",),
    ),
    _T(
        "sge_count",
        "★ VARIANT COUNTS FOR ONE SGE SAMPLE ★ — calls variants against the "
        "reference amplicon and reports per-variant counts, depth and "
        "frequency for a SINGLE timepoint. Counts are NOT functional "
        "scores: SGE measures depletion between timepoints, so use "
        "sge_analyze unless you specifically want one sample's composition "
        "(library QC, coverage checks).",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING},
                "target":    {**_S_STRING},
                "max_reads": {**_S_INTEGER},
                "min_count": {**_S_INTEGER, "default": 5},
                "label":     {**_S_STRING},
            },
            "required": ["accession"],
        },
        cli=["sge", "count"],
        positional=("accession",),
        flag_map={"target": "--target", "max_reads": "--max-reads",
                   "min_count": "--min-count", "label": "--label"},
    ),
    _T(
        "raw_pipeline_plan",
        "★ WHAT WOULD IT TAKE TO ANALYSE THIS DATASET'S RAW DATA ★ — "
        "resolves an IGVF FileSet accession, inventories its files, pairs "
        "the FASTQ reads, infers the sequencing chemistry from the seqspec, "
        "and reports the ROUTE it would take plus the download size. "
        "Downloads nothing except the small seqspec. CALL THIS FIRST "
        "whenever someone asks to analyse, process, or run a pipeline on a "
        "dataset accession, so the cost and the route are known before any "
        "large transfer. It also names which files are controlled-access. "
        "The route comes back as `assay_mismatch` for assays whose reads "
        "are NOT a transcript library -- SGE, MPRA, STARR-seq, protein "
        "scanning -- because quantifying those against a transcriptome "
        "restates the amplicon design instead of measuring anything; it "
        "then names the analysis that does apply. Do not force past that "
        "unless the user asked for gene counts knowing what they mean. "
        "This tool DESCRIBES the work; it does not do it. When the request "
        "was to analyse or process the data, follow it with "
        "raw_pipeline_run (detach=true for anything over a few GB) rather "
        "than reporting the plan as though it were the result.",
        {
            "type": "object",
            "properties": {
                "accession":   {**_S_STRING, "description":
                                 "IGVF FileSet accession (IGVFDS...)."},
                "technology":  {**_S_STRING, "description":
                                 "kb technology override, e.g. 10XV3. "
                                 "Inferred from the seqspec if omitted."},
                "force_align": {**_S_BOOLEAN, "description":
                                 "Plan the alignment route even if a "
                                 "published matrix already exists."},
            },
            "required": ["accession"],
        },
        cli=["raw-pipeline", "plan"],
        positional=("accession",),
        flag_map={"technology": "--technology"},
        bool_flags=("force_align",),
    ),
    _T(
        "raw_pipeline_run",
        "★ ACTUALLY ANALYSE A DATASET'S RAW DATA, END TO END ★ — the tool "
        "that DOES the work rather than describing it. Takes an IGVF "
        "accession and: prefers an already-published count matrix (on the "
        "set, or on an AnalysisSet derived from it); otherwise downloads "
        "the FASTQs and quantifies them with kallisto|bustools "
        "(`kb count`), inferring the chemistry from the seqspec; then runs "
        "the single-cell pipeline on the resulting matrix (QC, UMAP, "
        "Leiden, markers). USE THIS when the user asks to analyse or "
        "process raw data -- do NOT stop at explain_dataset and suggest "
        "commands for them to run. Set workflow='kite' to assign CRISPR "
        "sgRNA feature barcodes. Alignment can move tens of GB, so call "
        "raw_pipeline_plan first, or pass dry_run=true to see every "
        "command and transfer without performing any. A PLAN IS NOT AN "
        "ANSWER: if the user asked you to analyse or process the data, "
        "planning and then describing what would happen has not done what "
        "they asked. Start the run -- with detach=true when the plan "
        "reports more than a few GB -- and report the job id.",
        {
            "type": "object",
            "properties": {
                "accession":       {**_S_STRING, "description":
                                     "IGVF FileSet accession (IGVFDS...)."},
                "technology":      {**_S_STRING, "description":
                                     "kb technology override (e.g. 10XV3, "
                                     "BULK). Inferred from the seqspec if "
                                     "omitted."},
                "workflow":        {**_S_STRING, "description":
                                     "kb workflow: standard | nac | kite | "
                                     "kite:10xFB. 'kite' assigns feature "
                                     "barcodes such as CRISPR sgRNAs."},
                "reference":       {**_S_STRING, "default": "human",
                                     "description":
                                     "Prebuilt kb reference (human, mouse)."},
                "max_download_gb": {**_S_NUMBER, "default": 100,
                                     "description":
                                     "Refuse to transfer more than this."},
                "label":           {**_S_STRING},
                "dry_run":         {**_S_BOOLEAN, "description":
                                     "Print every command and download "
                                     "nothing."},
                "detach":          {**_S_BOOLEAN, "description":
                                     "★ USE THIS FOR ANYTHING LARGE ★ Start "
                                     "the run in its own process and return "
                                     "a job id immediately, instead of "
                                     "blocking until it finishes. A 45 GB "
                                     "dataset takes hours; a synchronous "
                                     "call holds the whole conversation open "
                                     "and the result is lost when the "
                                     "session ends. Poll with "
                                     "raw_pipeline_status."},
                "force_align":     {**_S_BOOLEAN, "description":
                                     "Align even if a published matrix "
                                     "exists, AND override the "
                                     "assay_mismatch refusal for a "
                                     "non-transcript assay (SGE, MPRA, "
                                     "STARR-seq, protein scanning). Only "
                                     "when the user has asked for gene "
                                     "counts knowing they do not answer "
                                     "that assay's question."},
                "skip_analysis":   {**_S_BOOLEAN, "description":
                                     "Stop at the matrix; skip the "
                                     "single-cell pipeline."},
                "no_reuse":        {**_S_BOOLEAN, "description":
                                     "Recompute from the reads even if this "
                                     "exact analysis was already completed "
                                     "on this machine. By default a repeat "
                                     "returns the existing matrix in "
                                     "seconds, which is what makes entering "
                                     "the same accession twice fast."},
            },
            "required": ["accession"],
        },
        cli=["raw-pipeline", "run"],
        positional=("accession",),
        flag_map={"technology": "--technology", "workflow": "--workflow",
                   "reference": "--reference", "label": "--label",
                   "max_download_gb": "--max-download-gb"},
        bool_flags=("dry_run", "force_align", "skip_analysis", "detach",
                     "no_reuse"),
    ),
    _T(
        "raw_pipeline_status",
        "★ HOW IS THE ALIGNMENT GOING ★ — reports on runs started with "
        "raw_pipeline_run(detach=true): whether each is still running, and "
        "the tail of its log. Call this when a user asks about a job "
        "already under way, or right after starting one to confirm it took. "
        "Jobs survive the conversation, so this also answers 'did that "
        "analysis I asked for yesterday finish'.",
        {
            "type": "object",
            "properties": {
                "job":   {**_S_STRING, "description":
                           "Job id from raw_pipeline_run. Omit for all."},
                "tail":  {**_S_INTEGER, "default": 12},
                "limit": {**_S_INTEGER, "default": 5},
            },
        },
        cli=["raw-pipeline", "status"],
        positional=("job",),
        flag_map={"tail": "--tail", "limit": "--limit"},
    ),
    _T(
        "raw_pipeline_assay_coverage",
        "★ WHICH IGVF ASSAYS CAN IGVFAGENT ACTUALLY ANALYSE ★ — reports "
        "every assay type on the Portal, the analysis route each needs, how "
        "many datasets that covers, and whether IGVFagent supports it. Use "
        "it to answer 'can you analyse assay X' honestly, or when a user "
        "asks what the system can do. Key fact it makes concrete: only "
        "59.5% of the 11,070 MeasurementSets are transcript assays, so "
        "quantifying reads against a transcriptome is the WRONG analysis "
        "for 40% of IGVF and would look plausible every time.",
        {"type": "object", "properties": {}},
        cli=["raw-pipeline", "assay-coverage"],
    ),
    _T(
        "raw_pipeline_guide_count",
        "★ ANALYSE A CRISPR SCREEN WHOSE READOUT IS gRNA SEQUENCING ★ — for "
        "a screen where `crispr_screen_readout` is 'gRNA sequencing', the "
        "reads ARE the guide library, not transcripts, and raw_pipeline "
        "refuses to quantify them against a transcriptome (it gave 1.4% "
        "pseudoalignment on IGVFDS6464SOVZ -- the measurement saying it was "
        "the wrong question). This counts how often each DESIGNED guide "
        "appears, matching spacers in both orientations against the "
        "library's own guide table, and reports per-guide counts, "
        "frequencies, the assignment rate and library coverage. On that "
        "dataset: 31% of reads assigned, 7,444 of 8,192 guides detected. "
        "One population gives library composition; a FACS screen scores "
        "guides by comparing sorted against unsorted, so run both and "
        "compare the freq columns.",
        {
            "type": "object",
            "properties": {
                "accession":       {**_S_STRING, "description":
                                     "MeasurementSet with gRNA-sequencing "
                                     "reads (IGVFDS...)."},
                "max_reads":       {**_S_INTEGER, "description":
                                     "Cap reads for a quick look."},
                "max_download_gb": {**_S_NUMBER, "default": 20},
                "label":           {**_S_STRING},
            },
            "required": ["accession"],
        },
        cli=["raw-pipeline", "guide-count"],
        positional=("accession",),
        flag_map={"max_reads": "--max-reads", "label": "--label",
                   "max_download_gb": "--max-download-gb"},
    ),
    _T(
        "raw_pipeline_guide_library",
        "★ IS THE sgRNA LIBRARY FOR THIS CRISPR SCREEN KNOWN YET ★ — a "
        "screen cannot have its guides assigned without the protospacer "
        "table, and IGVF records that link on the AnalysisSet "
        "(construct_library_sets -> integrated_content_files), NOT on the "
        "MeasurementSet. This resolves it properly and reports the guide "
        "count and spacer lengths, or says exactly why it cannot -- most "
        "often because the dataset is still 'in progress' and has no "
        "AnalysisSet, so the link does not exist yet. Call it before "
        "attempting workflow='kite', and to re-check a dataset that was "
        "not ready earlier. NEVER guess the library from its name: "
        "assigning guides against the wrong library yields confident, "
        "entirely wrong perturbation calls.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING, "description":
                               "IGVF FileSet accession (IGVFDS...)."},
            },
            "required": ["accession"],
        },
        cli=["raw-pipeline", "guide-library"],
        positional=("accession",),
    ),

    _T(
        "warehouse_query",
        "★ RUN SQL OVER THE LOCAL WAREHOUSE ★ (DuckDB, read-only). The "
        "only way to ASK A QUESTION of data already mirrored locally "
        "rather than re-fetching it: counts, GROUP BY, joins across "
        "ingested tables, top-N. Use `warehouse_stats` first to see "
        "which tables exist and how many rows each holds — then write "
        "SQL against those names. Read-only: SELECT works, writes are "
        "refused by the connection.",
        {
            "type": "object",
            "properties": {
                "sql":   {**_S_STRING, "description":
                           "A SELECT statement. Table names come from "
                           "warehouse_stats."},
                "limit": {**_S_INTEGER, "default": 50,
                           "description": "Max rows printed."},
            },
            "required": ["sql"],
        },
        cli=["warehouse", "query"],
        positional=("sql",),
        flag_map={"limit": "--limit"},
    ),

    _T(
        "warehouse_stats",
        "List every table in the local DuckDB warehouse with its row "
        "count. Run this before `warehouse_query` so the SQL references "
        "tables that actually exist.",
        {"type": "object", "properties": {}},
        cli=["warehouse", "stats"],
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
