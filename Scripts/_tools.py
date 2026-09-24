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
                "force":       {**_S_BOOLEAN, "description":
                                "Author even if the duplication guard says a "
                                "core tool already covers it."},
            },
            "required": ["name", "description"],
        },
        ["extauthor", "write-tool"],
        flag_map={"name": "--name", "description": "--description",
                   "cli": "--cli", "command": "--command",
                   "parameters": "--parameters", "positional": "--positional"},
        bool_flags={"force"},
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
                "force":       {**_S_BOOLEAN, "description":
                                "Author even if the duplication guard says a "
                                "core tool already covers it. Re-authoring "
                                "your OWN extension under the same name never "
                                "needs this: that is an update and is allowed."},
            },
            "required": ["name", "description"],
        },
        ["extauthor", "write-skill"],
        flag_map={"name": "--name", "description": "--description",
                   "source": "--source", "source_file": "--source-file",
                   "tool_parameters": "--tool-parameters"},
        bool_flags={"force"},
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
        "★ EXPLAIN **AND DOWNLOAD** AN IGVF / ENCODE ACCESSION OR URL ★. "
        "Metadata, file inventory, SVG overview plots, how-to-use report "
        "-- and with `download=true`, THE FILES THEMSELVES. **This is the "
        "tool for \"download ENCODE file ENCFFxxxxxxx\"**, or any bigWig / "
        "BED / BAM by accession. portal_get returns IGVF Portal JSON, not "
        "ENCODE and not file payloads; geo_download only accepts a GSE "
        "series -- neither can fetch an ENCFF file, and NEVER author a "
        "script to do it. Files land under Data/Interpreted/Downloads/ and "
        "each one's download_status is reported. For \"all the data for cell "
        "line X\" use biosample_portal_census instead of a search URL: a "
        "search that matches nothing answers 404 here.",
        {
            "type": "object",
            "properties": {
                "accession_or_url": {**_S_STRING},
                "download": {**_S_BOOLEAN, "default": False},
                "max_download_gb": {"type": "number", "default": 200.0,
                    "description": "Transfer ceiling in GB. Applies ONLY "
                        "when download=true; nothing is fetched without "
                        "it. A file over the ceiling is SKIPPED, not "
                        "truncated."},
                "include_raw": {**_S_BOOLEAN, "default": False, "description":
                    "IGVF only: with download=true, also fetch raw reads when the "
                    "Portal already holds processed results (default: fetch the "
                    "processed files instead)."},
            },
            "required": ["accession_or_url"],
        },
        cli=["explain", "explain"],
        positional=["accession_or_url"],
        flag_map={"max_download_gb": "--max-download-gb"},
        bool_flags={"download", "include_raw"},
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
        bool_flags=("cluster",),
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
        "biosample_portal_census",
        "★ SYSTEMATIC SUMMARY OF EVERYTHING ENCODE **AND** IGVF HOLD FOR ONE "
        "BIOSAMPLE / CELL LINE ★ (GM12878, K562, HepG2, WTC11, liver, ...). "
        "THE tool for \"summarise all ENCODE and IGVF data for X with tables "
        "and plots\": one call counts every object type on both portals via "
        "their structured sample-term filters (never free text, which on the "
        "IGVF Portal matches thousands of unrelated sets), tabulates assays, "
        "ChIP targets, labs, annotation types, file formats and release years, "
        "and writes report.md + CSV tables + SVG/PNG figures. Do NOT use "
        "explain_dataset on a search URL, encode_retrieve per assay, or author "
        "a new skill for this question. A zero-hit 404 from a portal is "
        "recorded as 0, not treated as an error.",
        {
            "type": "object",
            "properties": {
                "biosample": {**_S_STRING,
                    "description": "Cell line or tissue term as the portals "
                                    "spell it (case-insensitive; resolved to "
                                    "the ontology term)."},
                "label":     {**_S_STRING},
                "portal":    {**_S_STRING,
                    "description": "both (default), encode, or igvf."},
                "status":    {**_S_STRING,
                    "description": "Status filter for item tables: released "
                                    "(default) or all."},
                "max_items": {**_S_INTEGER, "default": 5000},
                "top":       {**_S_INTEGER, "default": 20,
                    "description": "Top-N terms per table and figure."},
                "no_plots":  {**_S_BOOLEAN, "default": False},
            },
            "required": ["biosample"],
        },
        cli=["biosample-census"],
        flag_map={"biosample": "--biosample", "label": "--label",
                   "portal": "--portal", "status": "--status",
                   "max_items": "--max-items", "top": "--top"},
        bool_flags={"no_plots"},
    ),

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
                "max_download_gb": {"type": "number", "default": 200.0,
                    "description": "Transfer ceiling in GB (guard, not a "
                        "sampling limit)."},
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
        bool_flags=("html",),
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

    # ── Single-cell CRISPR differential expression ────────────────────
    _T(
        "tf_perturb_seq_analyze",
        "★ IGVF TF PERTURB-SEQ: CALIBRATED RESULTS -> DISEASE / GWAS OVERLAY ★. "
        "Port of the IGVF tf_perturb_seq Working Group 3 jamboree notebook. "
        "Takes the IGVF CRISPR pipeline's CALIBRATED inference tables "
        "(<prefix>direct_target_results.tsv / cis_results.tsv / "
        "trans_results.tsv, the 22M-row trans table is streamed and cached) "
        "and produces: significant direct / cis / trans effects, top trans "
        "regulators with their up/down targets, GWAS Catalog SNPs near TF "
        "elements and near trans target genes with Fisher trait enrichment, "
        "scE2G enhancer-gene links for perturbed TFs and targets, GWAS SNPs "
        "inside E2G elements (incl. non-nearest-gene assignments), optional "
        "ChIP-seq bigWig support -- plus report.md, CSV tables and PNG "
        "figures. Use for 'which TFs regulate disease/GWAS genes in this "
        "lineage' on a TF Perturb-seq run. NOT for raw counts: that is "
        "sc_crispr_de_* / crispr_pipeline.",
        {
            "type": "object",
            "properties": {
                "calibrated_prefix": {**_S_STRING, "description":
                    "Path prefix of the calibrated TSVs, e.g. "
                    "data/<run>_calibrated_ (the three files are appended)."},
                "mudata":     {**_S_STRING, "description":
                    "inference_mudata.h5mu (guide overview + gene coordinates)."},
                "gene_coords": {**_S_STRING, "description":
                    "TSV with gene_id, symbol, chr, start, end when no mudata."},
                "gwas":       {**_S_STRING, "description":
                    "GWAS Catalog associations TSV (full download)."},
                "e2g":        {**_S_ARRAY_S, "description":
                    "scE2G TSVs as LABEL=PATH, e.g. ['ESC=h7.e2g.tsv', 'DE=de.e2g.tsv']."},
                "bigwig":     {**_S_ARRAY_S, "description":
                    "ChIP-seq bigWig paths for the support step (needs pyBigWig)."},
                "tf_list":    {**_S_STRING},
                "label":      {**_S_STRING},
                "padj":       {**_S_NUMBER, "default": 0.05},
                "lfc_cis":    {**_S_NUMBER, "default": 0.2},
                "lfc_trans":  {**_S_NUMBER, "default": 1.0},
                "top_regulators": {**_S_INTEGER, "default": 20},
                "gwas_window": {**_S_INTEGER, "default": 50000},
                "e2g_score_min": {**_S_NUMBER, "default": 0.177},
                "min_hits":   {**_S_INTEGER, "default": 5},
                "signal_frac": {**_S_NUMBER, "default": 0.1},
                "no_plots":   {**_S_BOOLEAN, "default": False},
                "no_cache":   {**_S_BOOLEAN, "default": False},
            },
            "required": ["calibrated_prefix"],
        },
        cli=["tf-perturb", "run"],
        flag_map={"calibrated_prefix": "--calibrated-prefix", "mudata": "--mudata",
                   "gene_coords": "--gene-coords", "gwas": "--gwas", "e2g": "--e2g",
                   "bigwig": "--bigwig", "tf_list": "--tf-list", "label": "--label",
                   "padj": "--padj", "lfc_cis": "--lfc-cis", "lfc_trans": "--lfc-trans",
                   "top_regulators": "--top-regulators", "gwas_window": "--gwas-window",
                   "e2g_score_min": "--e2g-score-min", "min_hits": "--min-hits",
                   "signal_frac": "--signal-frac"},
        flag_repeat={"e2g", "bigwig"},
        bool_flags={"no_plots", "no_cache"},
    ),

    _T(
        "tf_perturb_qc_gene",
        "★ TF PERTURB-SEQ STAGE 3 QC: GENE MAPPING ★ on inference_mudata.h5mu "
        "(IGVF CRISPR pipeline output): per-cell UMIs, genes detected and mito % "
        "(median / mean / sd / quartiles), overall and per batch; histograms and knee "
        "plot. Port of tf_perturb_seq/qc/mapping_gene.py.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING},
            "batch_col": {**_S_STRING, "default": 'batch'},
            "prefix": {**_S_STRING, "default": 'gene'},
            "label": {**_S_STRING},
            "no_plots": {**_S_BOOLEAN, "default": False}},
         "required": ['mudata']},
        cli=['tf-perturb', 'qc-gene'],
        flag_map={"mudata": "--mudata", "batch_col": "--batch-col", "prefix": "--prefix", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "tf_perturb_qc_guide",
        "★ TF PERTURB-SEQ STAGE 3 QC: GUIDE MAPPING ★: guide UMIs per cell, guides "
        "assigned per cell, cells per guide, fraction of cells with a guide, and a "
        "per-guide capture table (detected / assigned cells, UMI stats). Port of "
        "tf_perturb_seq/qc/mapping_guide.py.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING},
            "batch_col": {**_S_STRING, "default": 'batch'},
            "prefix": {**_S_STRING, "default": 'guide'},
            "label": {**_S_STRING},
            "no_plots": {**_S_BOOLEAN, "default": False}},
         "required": ['mudata']},
        cli=['tf-perturb', 'qc-guide'],
        flag_map={"mudata": "--mudata", "batch_col": "--batch-col", "prefix": "--prefix", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "tf_perturb_qc_target",
        "★ TF PERTURB-SEQ STAGE 3 QC: INTENDED-TARGET KNOCKDOWN ★ from the pipeline's "
        "trans_per_guide_results: strong knockdowns (fold change <= 0.4), significant "
        "tests, median log2FC, and the AUROC / AUPRC of targeting-vs-non-targeting "
        "guides on a balanced evaluation table (CRISPR_Pipeline evaluate_controls "
        "convention). Volcano + ROC/PR plots. Port of qc/intended_target.py.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING},
            "results_key": {**_S_STRING, "default": 'trans_per_guide_results'},
            "fc_threshold": {**_S_NUMBER, "default": 0.4},
            "pval_threshold": {**_S_NUMBER, "default": 0.05},
            "prefix": {**_S_STRING, "default": 'intended_target'},
            "label": {**_S_STRING},
            "no_plots": {**_S_BOOLEAN, "default": False}},
         "required": ['mudata']},
        cli=['tf-perturb', 'qc-target'],
        flag_map={"mudata": "--mudata", "results_key": "--results-key", "fc_threshold": "--fc-threshold", "pval_threshold": "--pval-threshold", "prefix": "--prefix", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "tf_perturb_calibrate",
        "★ TF PERTURB-SEQ DEG CALIBRATION ★: empirical p-values for PerTurbo per-element "
        "results (perturbo_trans_per_element_output.tsv) against the non-targeting "
        "null -- z = log2FC / SE, eCDF or t-fit -- BH on targeting tests only, cis "
        "(same chromosome within --cis-window) and direct-target annotation, then the "
        "four calibrated tables (all / direct_target / cis / trans) that "
        "tf_perturb_seq_analyze reads. Port of tf_perturb_seq/inference/calibrate.py.",
        {"type": "object", "properties": {
            "trans_results": {**_S_STRING, "description": 'perturbo_trans_per_element_output.tsv(.gz)'},
            "mudata": {**_S_STRING},
            "prefix": {**_S_STRING, "description": 'output prefix, e.g. <dataset>_<run>'},
            "method": {**_S_STRING, "default": 't-fit', "description": 'ecdf or t-fit'},
            "cis_window": {**_S_INTEGER, "default": 100000},
            "label": {**_S_STRING}},
         "required": ['trans_results', 'mudata', 'prefix']},
        cli=['tf-perturb', 'calibrate'],
        flag_map={"trans_results": "--trans-results", "mudata": "--mudata", "prefix": "--prefix", "method": "--method", "cis_window": "--cis-window", "label": "--label"},
    ),

    _T(
        "tf_perturb_pathways",
        "★ TF PERTURB-SEQ PATHWAYS ★: per perturbed element, Fisher over-representation of "
        "its significant DEGs (calibrated FDR < --fdr) in GMT gene sets (MSigDB, GO, "
        "KEGG), BH per element. Port of tf_perturb_seq/inference/pathways.py.",
        {"type": "object", "properties": {
            "calibrated": {**_S_STRING, "description": '<prefix>_calibrated_trans_results.tsv or all_results'},
            "gmt": {**_S_STRING},
            "prefix": {**_S_STRING},
            "fdr": {**_S_NUMBER, "default": 0.1},
            "min_genes": {**_S_INTEGER, "default": 3},
            "label": {**_S_STRING}},
         "required": ['calibrated', 'gmt', 'prefix']},
        cli=['tf-perturb', 'pathways'],
        flag_map={"calibrated": "--calibrated", "gmt": "--gmt", "prefix": "--prefix", "fdr": "--fdr", "min_genes": "--min-genes", "label": "--label"},
    ),

    _T(
        "tf_perturb_filter_cells",
        "Drop cells carrying more than N assigned guides from an inference_mudata.h5mu "
        "(reads the assignment layer directly from HDF5). Port of "
        "scripts/filter_mudata_by_sgrna_count.py.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING},
            "out": {**_S_STRING, "description": 'output .h5mu'},
            "max_guides": {**_S_INTEGER, "default": 15}},
         "required": ['mudata', 'out']},
        cli=['tf-perturb', 'filter-cells'],
        flag_map={"mudata": "--mudata", "out": "--out", "max_guides": "--max-guides"},
    ),

    _T(
        "tf_perturb_filter_guides",
        "Drop outlier guides (BH on the energy-distance outlier tables' pval_outlier) and "
        "every cell that received one, recomputing per-cell guide counts. Port of "
        "scripts/filter_outlier_guides.py.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING},
            "out": {**_S_STRING},
            "targeting_outliers": {**_S_STRING},
            "non_targeting_outliers": {**_S_STRING},
            "fdr": {**_S_NUMBER, "default": 0.05}},
         "required": ['mudata', 'out', 'targeting_outliers', 'non_targeting_outliers']},
        cli=['tf-perturb', 'filter-guides'],
        flag_map={"mudata": "--mudata", "out": "--out", "targeting_outliers": "--targeting-outliers", "non_targeting_outliers": "--non-targeting-outliers", "fdr": "--fdr"},
    ),

    _T(
        "tf_perturb_edist_prep",
        "★ TF PERTURB-SEQ ENERGY DISTANCE, STEP 0 ★: normalise / log1p / scale / PCA(50) "
        "on the gene modality, guide -> cells dictionary, annotation table with "
        "<ENSG>|chr:start-end target labels. Writes the folder edist-filter and edist read.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING},
            "out_dir": {**_S_STRING},
            "n_comps": {**_S_INTEGER, "default": 50},
            "label": {**_S_STRING}},
         "required": ['mudata']},
        cli=['tf-perturb', 'edist-prep'],
        flag_map={"mudata": "--mudata", "out_dir": "--out-dir", "n_comps": "--n-comps", "label": "--label"},
    ),

    _T(
        "tf_perturb_edist_filter",
        "★ ENERGY DISTANCE STEP 1: OUTLIER GUIDES ★: DISCO permutation test among sibling "
        "guides of each target, hypergeometric outlier ranking on intra-target energy "
        "distances, K-means (k=2) outliers among non-targeting guides. Writes "
        "targeting_outlier_table.csv / non_targeting_outlier_table.csv (pval_outlier). "
        "Clean-room CPU rewrite of Chikara-Takeuchi/energy_dist_pipeline step 1.",
        {"type": "object", "properties": {
            "prep_dir": {**_S_STRING},
            "out_dir": {**_S_STRING},
            "min_cells": {**_S_INTEGER, "default": 20},
            "disco_permutations": {**_S_INTEGER, "default": 200},
            "fdr": {**_S_NUMBER, "default": 0.05},
            "seed": {**_S_INTEGER, "default": 0}},
         "required": ['prep_dir']},
        cli=['tf-perturb', 'edist-filter'],
        flag_map={"prep_dir": "--prep-dir", "out_dir": "--out-dir", "min_cells": "--min-cells", "disco_permutations": "--disco-permutations", "fdr": "--fdr", "seed": "--seed"},
    ),

    _T(
        "tf_perturb_edist",
        "★ ENERGY DISTANCE STEP 2: PER-TARGET E-DISTANCE VS NON-TARGETING ★ with "
        "permutation p-values: num_bg random backgrounds of non_target_pick NT cells, "
        "permutations label shuffles each; distance_i / pval_i per background, "
        "distance_mean, pval_mean, pval_mean_log, distance_mean_log, cell_count, type "
        "-> pval_edist_full.csv (the upstream schema). Outlier guides from edist-filter "
        "are excluded. Clean-room CPU rewrite; use the upstream GPU container for "
        "hundreds of thousands of cells.",
        {"type": "object", "properties": {
            "prep_dir": {**_S_STRING},
            "out_dir": {**_S_STRING},
            "num_bg": {**_S_INTEGER, "default": 20},
            "permutations": {**_S_INTEGER, "default": 1000},
            "non_target_pick": {**_S_INTEGER, "default": 2000},
            "target_cell_max": {**_S_INTEGER, "default": 2000},
            "min_cells": {**_S_INTEGER, "default": 20},
            "keep_outliers": {**_S_BOOLEAN, "default": False},
            "seed": {**_S_INTEGER, "default": 0}},
         "required": ['prep_dir']},
        cli=['tf-perturb', 'edist'],
        flag_map={"prep_dir": "--prep-dir", "out_dir": "--out-dir", "num_bg": "--num-bg", "permutations": "--permutations", "non_target_pick": "--non-target-pick", "target_cell_max": "--target-cell-max", "min_cells": "--min-cells", "seed": "--seed"},
        bool_flags={'keep_outliers'},
    ),

    _T(
        "tf_perturb_edist_summary",
        "Cross-dataset roll-up of energy-distance results: per dataset the target counts by "
        "type, median distance_mean by type, targets above the negative-control maximum "
        "(the calibration-robust significance proxy), pval_mean == 0 and < 0.05 counts. "
        "Port of energy_dist/cross_dataset_edistance_summary.py, on local files.",
        {"type": "object", "properties": {
            "inputs": {**_S_ARRAY_S, "description": 'NAME=dir-or-pval_edist_full.csv entries'},
            "label": {**_S_STRING}},
         "required": ['inputs']},
        cli=['tf-perturb', 'edist-summary'],
        flag_map={"inputs": "--inputs", "label": "--label"},
        flag_repeat={'inputs'},
    ),

    _T(
        "tf_perturb_edist_validate",
        "Schema and value-range checks on an energy-distance output folder "
        "(pval_edist_full.csv required columns, per-background columns, cell_count > 0, "
        "pval_mean in [0,1], known type values; outlier tables' pval_outlier).",
        {"type": "object", "properties": {
            "dir": {**_S_STRING}},
         "required": ['dir']},
        cli=['tf-perturb', 'edist-validate'],
        flag_map={"dir": "--dir"},
    ),

    _T(
        "tf_perturb_cnmf_export",
        "★ TF PERTURB-SEQ STAGE 5 INPUT ★: inference_mudata.h5mu -> PerturbNMF / torch-cNMF "
        "single-modality .h5ad (gene symbols as var_names, guide_assignment in obsm, "
        "guide_names / guide_targets in uns, cells with zero counts in the HVG set "
        "dropped). Port of cnmf/h5mu_to_perturbnmf_h5ad.py.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING},
            "out_h5ad": {**_S_STRING},
            "num_highvar_genes": {**_S_INTEGER, "default": 2000},
            "no_hvg_filter": {**_S_BOOLEAN, "default": False}},
         "required": ['mudata', 'out_h5ad']},
        cli=['tf-perturb', 'cnmf-export'],
        flag_map={"mudata": "--mudata", "out_h5ad": "--out-h5ad", "num_highvar_genes": "--num-highvar-genes"},
        bool_flags={'no_hvg_filter'},
    ),

    _T(
        "tf_perturb_cnmf_validate",
        "File / shape / value checks on a torch-cNMF run folder for the selected k: "
        "gene_spectra_score is k x genes, usages is cells x k and non-negative, no NaNs; "
        "lists the k values present. Port of cnmf/validate_cnmf_outputs.py layers 1-3.",
        {"type": "object", "properties": {
            "dir": {**_S_STRING},
            "selected_k": {**_S_INTEGER}},
         "required": ['dir', 'selected_k']},
        cli=['tf-perturb', 'cnmf-validate'],
        flag_map={"dir": "--dir", "selected_k": "--selected-k"},
    ),

    _T(
        "tf_perturb_pipeline_summary",
        "One row per dataset from the Stage 3 QC metric tables (cells, median UMIs, "
        "mito %, guide UMIs, guides per cell, fraction with guide, intended-target "
        "significance and AUROC / AUPRC). Port of "
        "crispr_pipeline/cross_dataset_pipeline_summary.py, on local folders.",
        {"type": "object", "properties": {
            "datasets": {**_S_ARRAY_S, "description": 'NAME=folder entries (folders holding *_metrics.tsv)'},
            "label": {**_S_STRING}},
         "required": ['datasets']},
        cli=['tf-perturb', 'pipeline-summary'],
        flag_map={"datasets": "--datasets", "label": "--label"},
        flag_repeat={'datasets'},
    ),

    _T(
        "tf_perturb_seq_selftest",
        "Self-test of tf_perturb_seq_analyze on synthetic calibrated tables, "
        "a synthetic GWAS catalog and E2G files with planted signals; checks "
        "the planted regulators, trait and enhancer-SNP are recovered. Run "
        "this to confirm the environment (pandas/scipy/matplotlib) before a "
        "real run.",
        {"type": "object", "properties": {
            "no_plots": {**_S_BOOLEAN, "default": False}}},
        cli=["tf-perturb", "selftest"],
        bool_flags={"no_plots"},
    ),

    _T(
        "sc_crispr_de_prepare",
        "★ SINGLE-CELL CRISPR SCREEN: STEP 1, PREPARE ★. 10x gene-expression "
        "+ sgRNA matrices -> a filtered checkpoint with each cell's guide "
        "call. Cells carrying exactly one guide become the test groups; cells "
        "with NO detected guide are the shared background every guide is "
        "tested against. Use for Perturb-seq / CRISPRi single-cell screens "
        "where you want per-guide differential expression.",
        {
            "type": "object",
            "properties": {
                "gex":    {**_S_STRING, "description":
                            "10x filtered_feature_bc_matrix dir, or .h5ad."},
                "sgrna":  {**_S_STRING, "description":
                            "sgRNA matrix dir/.h5ad. Omit if guide features "
                            "are in `gex` (CRISPR Guide Capture)."},
                "min_genes":  {**_S_INTEGER, "default": 200},
                "min_counts": {**_S_INTEGER, "default": 500},
                "label":  {**_S_STRING},
            },
            "required": ["gex"],
        },
        cli=["sc-crispr-de", "prepare"],
        flag_map={"gex": "--gex", "sgrna": "--sgrna",
                   "min_genes": "--min-genes", "min_counts": "--min-counts",
                   "label": "--label"},
    ),

    _T(
        "sc_crispr_de_test",
        "★ SINGLE-CELL CRISPR SCREEN: STEP 2, PER-GUIDE DE ★. Negative-"
        "binomial GLM per guide per gene, guide-bearing cells vs no-guide "
        "cells, with library size as a covariate. This is the compute-heavy "
        "step: upstream runs one SLURM array task per guide, here guides run "
        "across a process pool. Restrict with `gene_whitelist` or `guides` "
        "when you only care about a locus.",
        {
            "type": "object",
            "properties": {
                "checkpoint": {**_S_STRING, "description":
                                "Checkpoint .h5ad from sc_crispr_de_prepare."},
                "latentvar":  {**_S_STRING, "default": "nCount_RNA",
                                "description": "Comma list of .obs covariates."},
                "gene_whitelist": {**_S_STRING, "description":
                                    "Path to a plain-text gene list."},
                "guides":     {**_S_STRING, "description": "Comma list."},
                "max_guides": {**_S_INTEGER},
                "workers":    {**_S_INTEGER},
                "label":      {**_S_STRING},
            },
            "required": ["checkpoint"],
        },
        cli=["sc-crispr-de", "test"],
        flag_map={"checkpoint": "--checkpoint", "latentvar": "--latentvar",
                   "gene_whitelist": "--gene-whitelist", "guides": "--guides",
                   "max_guides": "--max-guides", "workers": "--workers",
                   "label": "--label"},
    ),

    _T(
        "sc_crispr_de_aggregate",
        "★ SINGLE-CELL CRISPR SCREEN: STEP 3, GENE LEVEL ★. Aggregates the "
        "per-guide p-values into one score per target-gene pair using "
        "alpha-RRA, calibrated against the non-targeting guides when enough "
        "are present. Emits rra_rho / rra_pval / rra_fdr / rra_effect_size. "
        "NOTE these are alpha-RRA, NOT FRACTEL: the same family of method as "
        "the upstream pipeline, deliberately named differently because the "
        "numbers are not interchangeable.",
        {
            "type": "object",
            "properties": {
                "per_guide":    {**_S_STRING, "description":
                                  "TSV from sc_crispr_de_test."},
                "guide_map":    {**_S_STRING, "description":
                                  "guide<TAB>target lines; default splits "
                                  "<gene>_<n>."},
                "nontargeting": {**_S_STRING, "description":
                                  "Comma list of NTC guides — without these "
                                  "the null is assumed, not measured."},
                "alpha":        {**_S_NUMBER, "default": 0.25},
                "label":        {**_S_STRING},
            },
            "required": ["per_guide"],
        },
        cli=["sc-crispr-de", "aggregate"],
        flag_map={"per_guide": "--per-guide", "guide_map": "--guide-map",
                   "nontargeting": "--nontargeting", "alpha": "--alpha",
                   "label": "--label"},
    ),

    _T(
        "spatial_hic_matrix",
        "Build a Hi-C contact matrix from Spatial-ATAC-Hi-C pairs, for one "
        "chromosome or a region. This is what the paper's Hi-C MAP FIGURES "
        "are drawn from (Fig. 1b,c; Fig. 3c-h; Fig. 5d-g). Point `pairs_dir` "
        "at a whole pairs file for pseudobulk, or at a per-pixel directory "
        "from `spatial_hic_pixel_demux`. Pair it with spatial_hic_impute when "
        "a single pixel is too sparse to show structure.",
        {
            "type": "object",
            "properties": {
                "pairs_dir":   {**_S_STRING, "description":
                                 "A contact file, or a directory of per-pixel "
                                 "ones (4DN .pairs or GEO contact TSVs)."},
                "chrom_sizes": {**_S_STRING},
                "pairs_glob":  {**_S_STRING},
                "chrom":       {**_S_STRING, "description": "e.g. chr2"},
                "resolution":  {**_S_INTEGER, "description":
                                 "Bin size in bp. The paper renders at 25000."},
                "start":       {**_S_INTEGER},
                "end":         {**_S_INTEGER},
                "label":       {**_S_STRING},
            },
            "required": ["pairs_dir"],
        },
        cli=["spatial-hic", "matrix"],
        flag_map={"pairs_dir": "--pairs-dir", "chrom_sizes": "--chrom-sizes",
                   "pairs_glob": "--pairs-glob", "chrom": "--chrom",
                   "resolution": "--resolution", "start": "--start",
                   "end": "--end", "label": "--label"},
    ),

    _T(
        "spatial_hic_impute",
        "★ IMPUTE SPARSE SINGLE-PIXEL Hi-C CONTACT MAPS ★ (scHiCluster-style "
        "linear convolution + random walk with restart). A single 50x50 "
        "tissue pixel holds only tens of thousands of contacts, far too few "
        "to read structure off directly, so the paper imputes before "
        "compartment calling and before drawing single-pixel maps "
        "(Methods, 'Imputation of Spatial-ATAC-Hi-C data'; Fig. 1b). Its "
        "resolutions: 100000 for compartments, 25000 for visualisation, "
        "10000 for fine structure. Use `per_pixel` to impute every pixel "
        "rather than the pseudobulk.",
        {
            "type": "object",
            "properties": {
                "pairs_dir":    {**_S_STRING, "description":
                                  "A contact file, or a directory of per-pixel "
                                  "ones from spatial_hic_pixel_demux."},
                "chrom_sizes":  {**_S_STRING},
                "pairs_glob":   {**_S_STRING},
                "chrom":        {**_S_STRING, "description": "e.g. chr2"},
                "resolution":   {**_S_INTEGER, "description":
                                  "100000 compartments / 25000 visualisation "
                                  "/ 10000 fine structure, as in the paper."},
                "start":        {**_S_INTEGER},
                "end":          {**_S_INTEGER},
                "pad":          {**_S_INTEGER, "description":
                                  "Convolution half-width (scHiCluster pad)."},
                "restart":      {**_S_NUMBER, "description":
                                  "Random-walk restart probability."},
                "tol":          {**_S_NUMBER},
                "per_pixel":    {**_S_BOOLEAN, "description":
                                  "Impute each pixel separately."},
                "zero_diagonal": {**_S_BOOLEAN},
                "min_contacts": {**_S_INTEGER, "description":
                                  "Skip pixels below this contact count."},
                "label":        {**_S_STRING},
            },
            "required": ["pairs_dir"],
        },
        cli=["spatial-hic", "impute"],
        flag_map={"pairs_dir": "--pairs-dir", "chrom_sizes": "--chrom-sizes",
                   "pairs_glob": "--pairs-glob", "chrom": "--chrom",
                   "resolution": "--resolution", "start": "--start",
                   "end": "--end", "pad": "--pad", "restart": "--restart",
                   "tol": "--tol", "min_contacts": "--min-contacts",
                   "label": "--label"},
        bool_flags={"per_pixel", "zero_diagonal"},
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
        bool_flags=("restart",),
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
        bool_flags=("include_giants", "restart",),
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
        bool_flags=("shift_correct",),
    ),
    _T(
        "share_align_atac",
        "★ SHARE-seq ATAC FASTQs → FRAGMENTS ★ via chromap. The alignment "
        "stage SHARE-seq-alignmentV2 performs, run on the tools installed "
        "here: that pipeline is GPL-3.0 and wraps STAR/bowtie2/fastp/"
        "umi_tools/samtools/Picard/featureCounts/bedtools, none of which are "
        "present. chromap IS (it is half the IGVF uniform pipeline). The "
        "SHARE-seq part — three 24-mer combinatorial barcodes at fixed "
        "offsets in the barcode read — is expressed through chromap's "
        "read-format rather than reimplemented. A DIFFERENT aligner from "
        "upstream's: expect concordant fragments, not identical ones. Feed "
        "the output to share_fragment_qc.",
        {
            "type": "object",
            "properties": {
                "read1":     {**_S_STRING}, "read2": {**_S_STRING},
                "barcode":   {**_S_STRING, "description": "FASTQ with the barcode read."},
                "index":     {**_S_STRING, "description": "chromap index."},
                "ref":       {**_S_STRING, "description": "Reference FASTA."},
                "whitelist": {**_S_STRING},
                "r1_offset": {**_S_INTEGER, "default": 14},
                "r2_offset": {**_S_INTEGER, "default": 52},
                "r3_offset": {**_S_INTEGER, "default": 90},
                "threads":   {**_S_INTEGER, "default": 4},
                "label":     {**_S_STRING},
            },
            "required": ["read1", "read2", "barcode", "index", "ref"],
        },
        cli=["share", "align-atac"],
        flag_map={"read1": "--read1", "read2": "--read2", "barcode": "--barcode",
                   "index": "--index", "ref": "--ref", "whitelist": "--whitelist",
                   "r1_offset": "--r1-offset", "r2_offset": "--r2-offset",
                   "r3_offset": "--r3-offset", "threads": "--threads",
                   "label": "--label"},
    ),

    _T(
        "share_align_rna",
        "★ SHARE-seq RNA FASTQs → COUNT MATRIX ★ via kb-python. The RNA half "
        "of the SHARE-seq alignment stage, using kallisto|bustools where "
        "upstream uses STAR + featureCounts. The barcode geometry is passed "
        "as a kb technology string, derived from the SHARE-seq offsets when "
        "not given. Feed the output to share_rna_qc, then share_joint_qc to "
        "pair it with the ATAC half.",
        {
            "type": "object",
            "properties": {
                "fastqs":     {**_S_STRING, "description":
                                "Space-separated FASTQ paths, in kb order."},
                "index":      {**_S_STRING, "description": "kallisto index."},
                "t2g":        {**_S_STRING},
                "technology": {**_S_STRING, "description":
                                "kb -x string; derived when omitted."},
                "threads":    {**_S_INTEGER, "default": 4},
                "label":      {**_S_STRING},
            },
            "required": ["fastqs", "index", "t2g"],
        },
        cli=["share", "align-rna"],
        positional=("fastqs",),
        flag_map={"index": "--index", "t2g": "--t2g",
                   "technology": "--technology", "threads": "--threads",
                   "label": "--label"},
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
        bool_flags=("fetch",),
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
        "paper_benchmark",
        "★ REPRODUCE AND BENCHMARK ANY PUBLISHED PAPER ★ — THE tool for "
        "\"reproduce this paper\", \"can IGVFagent replicate X\", "
        "\"benchmark against this study\", or a pasted title / DOI / PMID / "
        "URL. It resolves the publication, harvests the accessions and "
        "numeric claims out of its own text, ROUTES the deposit onto the "
        "right IGVFagent analysis chain (Spatial-ATAC-Hi-C, MPRA, Flow-FISH, "
        "Perturb-seq, SHARE-seq, MaveDB, ...), scaffolds "
        "`Benchmarks/<paper-id>/`, runs it, and scores what we measured "
        "against what the paper published. USE THIS FIRST for any "
        "reproduction request — do NOT reach for a paper-specific benchmark "
        "tool unless the question names that exact paper, and do NOT adopt a "
        "document found on disk when the conversation already named a paper. "
        "Subcommands: `resolve` (identifier -> one paper), `harvest` (text -> "
        "accessions + claims), `route` (-> analysis chain; `list-routes` "
        "shows the table), `scaffold`, `run`, `score`, `report`, and "
        "`pipeline` (resolve -> harvest -> route -> scaffold; add "
        "execute=true to also run -> score -> report).",
        {
            "type": "object",
            "properties": {
                "subcommand": {**_S_STRING, "description":
                                "pipeline | resolve | harvest | route | "
                                "scaffold | run | score | report | "
                                "list-routes | selftest"},
                "query":    {**_S_STRING, "description":
                              "Title, URL, DOI, PMID or free text naming the "
                              "paper. The usual entry point."},
                "doi":      {**_S_STRING},
                "pmid":     {**_S_STRING},
                "url":      {**_S_STRING},
                "paper_id": {**_S_STRING, "description":
                              "Benchmark id, e.g. wang2026_spatial_atac_hic. "
                              "REQUIRED for run / score / report / scaffold."},
                "route":    {**_S_STRING, "description":
                              "Force a route instead of the ranked pick, e.g. "
                              "spatial_atac_hic."},
                "execute":  {**_S_BOOLEAN, "description":
                              "For `pipeline`: also run, score and report."},
                "force":    {**_S_BOOLEAN, "description":
                              "Overwrite an existing Benchmarks/<paper-id>/."},
            },
            "required": ["subcommand"],
        },
        cli=["bench"],
        positional=("subcommand",),
        flag_map={"query": "--query", "doi": "--doi", "pmid": "--pmid",
                   "url": "--url", "paper_id": "--paper-id",
                   "route": "--route", "execute": "--execute",
                   "force": "--force"},
        bool_flags={"execute", "force"},
    ),

    _T(
        "pgboost_train",
        "★ LEARN TO COMBINE PEAK→GENE LINK EVIDENCE ★ (pgBoost approach). "
        "Several methods score a peak-gene link and they disagree — Signac, "
        "SCENT, Cicero, distance, and this repo's sce2g_predict. Picking one "
        "is arbitrary and averaging ignores that they are not equally "
        "trustworthy. This trains a gradient-boosted classifier on links with "
        "known labels (e.g. fine-mapped eQTLs) and reports honest metrics via "
        "LEAVE-ONE-CHROMOSOME-OUT CV — random CV inflates the score, because "
        "links on one chromosome share LD structure and near-duplicates leak "
        "into the test split. Needs a chromosome column.",
        {
            "type": "object",
            "properties": {
                "training_file": {**_S_STRING, "description":
                                   "Labelled links, TSV/CSV, with a chrom "
                                   "column and numeric predictors."},
                "label_col":     {**_S_STRING, "default": "label"},
                "predictors":    {**_S_STRING, "description":
                                   "Comma list; numeric columns if omitted."},
                "label":         {**_S_STRING},
            },
            "required": ["training_file"],
        },
        cli=["pgboost", "train"],
        flag_map={"training_file": "--training-file", "label_col": "--label-col",
                   "predictors": "--predictors", "label": "--label"},
    ),

    _T(
        "pgboost_predict",
        "★ SCORE CANDIDATE PEAK→GENE LINKS ★ with a model from pgboost_train. "
        "Emits pgBoost_probability and pgBoost_percentile alongside every "
        "input column. Refuses, naming the missing columns, when the input "
        "lacks a predictor the model was trained on — scoring against absent "
        "features would return confident numbers from nothing.",
        {
            "type": "object",
            "properties": {
                "data_file": {**_S_STRING, "description": "Candidate links."},
                "model":     {**_S_STRING, "description":
                               "pgboost_model.pkl from pgboost_train."},
                "label":     {**_S_STRING},
            },
            "required": ["data_file", "model"],
        },
        cli=["pgboost", "predict"],
        flag_map={"data_file": "--data-file", "model": "--model",
                   "label": "--label"},
    ),

    _T(
        "scqers_activity",
        "★ IS THIS ENHANCER ACTIVE, AND HOW STRONGLY ★ from a scQers "
        "single-cell reporter experiment. Bootstraps each element's own "
        "integrations against the minP and noP controls AT MATCHED SAMPLE "
        "SIZE, so a rare element is not beaten by a control that was simply "
        "pooled over more cells. Reports median bootstrap activity, the "
        "control median, the best-expressing cluster, and an empirical "
        "p-value. Input is a joined count table with one row per integration "
        "per cell (cell_bc, CRE_id, UMIs_mBC, gex_UMI, cluster_id, biol_rep) "
        "— the shape of GSE217686 / GSE217689.",
        {
            "type": "object",
            "properties": {
                "counts":     {**_S_STRING, "description":
                                "Joined count table, one row per integration "
                                "per cell."},
                "bootstraps": {**_S_INTEGER, "default": 200},
                "cre":        {**_S_STRING, "description":
                                "Restrict to a single element."},
                "cluster_groups": {**_S_STRING, "description":
                                    "JSON mapping lineage -> cluster ids, to "
                                    "analyse at coarse lineage level."},
            },
            "required": ["counts"],
        },
        cli=["scqers", "activity"],
        flag_map={"counts": "--counts", "bootstraps": "--bootstraps",
                   "cre": "--cre", "cluster_groups": "--cluster-groups"},
    ),

    _T(
        "scqers_specificity",
        "★ IS THIS ENHANCER CELL-TYPE SPECIFIC ★ — a fold-change between "
        "clusters cannot answer this on its own, because the BEST of several "
        "clusters is high by construction. This permutes the cell-to-cluster "
        "assignment and recomputes the best cluster's fold-change, giving the "
        "null for 'how specific does this look when cell type carries no "
        "information'. Cells are permuted rather than rows, so two "
        "integrations in one cell keep the same label. Reports observed "
        "fold-change against the permuted 95th percentile plus an empirical "
        "p-value.",
        {
            "type": "object",
            "properties": {
                "counts":       {**_S_STRING},
                "permutations": {**_S_INTEGER, "default": 200},
                "cre":          {**_S_STRING},
                "cluster_groups": {**_S_STRING},
            },
            "required": ["counts"],
        },
        cli=["scqers", "specificity"],
        flag_map={"counts": "--counts", "permutations": "--permutations",
                   "cre": "--cre", "cluster_groups": "--cluster-groups"},
    ),

    _T(
        "scqers_pipeline",
        "★ FULL scQers ANALYSIS ★ — bootstrap activity, permutation "
        "specificity, then BH-corrected calls in one run. An element counts "
        "as reproducibly active or specific only if it clears the FDR in "
        "EVERY replicate it was measured in, which is the criterion the "
        "method is built around: one replicate is an observation, agreement "
        "across replicates is the claim. Emits activity.tsv, specificity.tsv "
        "and calls.tsv.",
        {
            "type": "object",
            "properties": {
                "counts":       {**_S_STRING},
                "bootstraps":   {**_S_INTEGER, "default": 200},
                "permutations": {**_S_INTEGER, "default": 200},
                "fdr":          {**_S_NUMBER, "default": 0.05},
                "min_fc":       {**_S_NUMBER, "default": 1.5,
                                  "description": "Fold-change floor for "
                                                 "calling an element specific."},
                "min_reps":     {**_S_INTEGER, "default": 1},
                "cluster_groups": {**_S_STRING},
            },
            "required": ["counts"],
        },
        cli=["scqers", "pipeline"],
        flag_map={"counts": "--counts", "bootstraps": "--bootstraps",
                   "permutations": "--permutations", "fdr": "--fdr",
                   "min_fc": "--min-fc", "min_reps": "--min-reps",
                   "cluster_groups": "--cluster-groups"},
    ),

    _T(
        "scqers_extract_bc",
        "★ PULL BARCODES OUT OF scQers / MPRA FASTQ ★ at fixed positions, "
        "with a constant-region sequence check. The check is the point: a "
        "read with an upstream indel still yields something at those "
        "positions, and nothing downstream can tell that apart from a real "
        "barcode. Handles one read or a paired barcode amplicon (oBC+mBC). "
        "Follow with scqers_count_bc to threshold.",
        {
            "type": "object",
            "properties": {
                "in_r1":     {**_S_STRING, "description": "R1 FASTQ (.gz ok)."},
                "in_r2":     {**_S_STRING, "description":
                               "R2, for paired barcode amplicons."},
                "out_file":  {**_S_STRING},
                "start":     {**_S_INTEGER, "default": 0},
                "end":       {**_S_INTEGER, "default": 15},
                "check_seq": {**_S_STRING, "description":
                               "Constant sequence expected right after the "
                               "barcode, e.g. GCT."},
            },
            "required": ["in_r1", "out_file"],
        },
        cli=["scqers", "extract-bc"],
        flag_map={"in_r1": "--in-r1", "in_r2": "--in-r2",
                   "out_file": "--out-file", "start": "--start",
                   "end": "--end", "check_seq": "--check-seq"},
    ),

    _T(
        "scqers_subassembly",
        "★ BUILD THE oBC <-> mBC BARCODE DICTIONARY ★ from a counted pair "
        "table, keeping only unambiguous pairings. A barcode mapping to two "
        "partners is a chimera or a collision, and keeping it attributes one "
        "element's reporter counts to another — the most damaging error "
        "available in this assay because it is invisible downstream. Requires "
        "the best partner to beat the runner-up by --min-ratio rather than "
        "just a majority.",
        {
            "type": "object",
            "properties": {
                "in_file":   {**_S_STRING},
                "out_file":  {**_S_STRING},
                "col1":      {**_S_STRING, "default": "oBC"},
                "col2":      {**_S_STRING, "default": "mBC"},
                "min_count": {**_S_INTEGER, "default": 3},
                "min_ratio": {**_S_NUMBER, "default": 5.0},
            },
            "required": ["in_file", "out_file"],
        },
        cli=["scqers", "subassembly"],
        flag_map={"in_file": "--in-file", "out_file": "--out-file",
                   "col1": "--col1", "col2": "--col2",
                   "min_count": "--min-count", "min_ratio": "--min-ratio"},
    ),

    _T(
        "upstream_check",
        "★ HAS AN UPSTREAM PROJECT CHANGED SINCE WE BUILT AGAINST IT ★ — many "
        "skills here reimplement a published method whose reference "
        "implementation lives in someone else's repository. This compares the "
        "revision each was pinned to against that project's current head and "
        "latest release, and reports the drift with a compare link. Use it "
        "when asked whether a method is current, when a result disagrees with "
        "a published one, or before trusting a reimplementation for new work.",
        {"type": "object", "properties": {}, "required": []},
        cli=["upstream", "check"],
    ),

    _T(
        "upstream_list",
        "★ WHERE DID THIS METHOD COME FROM ★ — the upstream project behind "
        "each skill, the exact revision it was built against, its licence, "
        "and whether it was reimplemented clean-room, ported, vendored, or "
        "invoked as an installed binary. Answers provenance and licensing "
        "questions precisely instead of from a docstring's prose.",
        {
            "type": "object",
            "properties": {
                "unpinned": {**_S_BOOLEAN, "description":
                              "Only those with no recorded revision."},
                "verbose":  {**_S_BOOLEAN},
            },
            "required": [],
        },
        cli=["upstream", "list"],
        flag_map={"unpinned": "--unpinned", "verbose": "--verbose"},
        bool_flags=("unpinned", "verbose"),
    ),

    _T(
        "processed_find",
        "★ HAS IGVF ALREADY PROCESSED THIS? ★ CALL THIS BEFORE PLANNING ANY "
        "ALIGNMENT OR REPROCESSING. A MeasurementSet's own files are raw "
        "reads, but `input_for` points at AnalysisSets, and one whose "
        "`uniform_pipeline_status` is `completed` already holds the matrices, "
        "fragments, peaks and alignments IGVF derived from it. On "
        "IGVFDS9875NBZW that is a 3.7 GB h5ad and a 5.6 GB PUBLIC fragments "
        "file, versus 72 GB of CONTROLLED FASTQ and hours of compute to "
        "recreate them. Access matters as much as size: raw reads often need "
        "credentials the caller does not have, while a derived file may be "
        "public. Lists every processed output with its size and access.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING, "description":
                               "MeasurementSet, AnalysisSet or sample accession."},
            },
            "required": ["accession"],
        },
        cli=["processed", "find"],
        positional=("accession",),
    ),

    _T(
        "job_start",
        "★ START A LONG-RUNNING AGENT JOB ★ for work that cannot finish in one "
        "reply: reproducing a paper, a full pipeline from raw reads, many "
        "downloads, anything that needs a plan of several stages or hours of "
        "compute. The job runs in its own background process (it survives the "
        "browser closing), plans its stages on disk, has every stage checked "
        "by the harness, retries failures, and is reviewed by an independent "
        "verifier before it reports done. Returns the job id; tell the user "
        "it is running and that progress appears in the Jobs panel. Use "
        "orchestrator=claude_code to drive it with Claude Code on the "
        "Anthropic API instead of IGVFagent's own loop.",
        {"type": "object", "properties": {
            "query": {**_S_STRING, "description": "The full task, with every detail the job needs (paper, accessions, what to reproduce)."},
            "orchestrator": {**_S_STRING, "description": "internal (default) or claude_code."},
            "budget_minutes": {**_S_NUMBER, "description": "Wall-clock budget (default 240)."},
            "max_rounds": {**_S_INTEGER, "description": "Maximum agent rounds (default 12)."}},
         "required": ["query"]},
        cli=["job", "start"], positional=["query"],
        flag_map={"budget_minutes": "--budget-minutes", "max_rounds": "--max-rounds"},
    ),

    _T(
        "job_status",
        "Status of a long-running agent job: state, rounds, cost, the "
        "harness-verified plan (done / failed / blocked stages), the "
        "verifier's verdict and the latest answer or final report.",
        {"type": "object", "properties": {"job": {**_S_STRING, "description": "Job id (J...)."}},
         "required": ["job"]},
        cli=["job", "status"], positional=["job"],
    ),

    _T(
        "plan_set",
        "INSIDE A JOB: record the staged plan before starting multi-step work. "
        "stages_json is a JSON list of stages, each {id, title, goal, success, "
        "evidence: [paths the stage will write], check}. check is optional and "
        "run by the harness: {kind: files, paths} | {kind: json, path, key, op, "
        "value} | {kind: rows, path, min} | {kind: concordance, benchmark}. "
        "Calling it again revises the plan and keeps finished stages.",
        {"type": "object", "properties": {"stages_json": {**_S_STRING, "description": "JSON list of stages."}},
         "required": ["stages_json"]},
        cli=["job", "plan-set"], flag_map={"stages_json": "--stages-json"},
    ),

    _T(
        "plan_update",
        "INSIDE A JOB: update a stage. status=done is verified by the harness "
        "(the evidence files must exist and the stage's check must pass); if not, "
        "the stage becomes failed and the reason is returned. Use failed to "
        "record an attempt, and blocked (with a note giving the reason) after 3 "
        "failed attempts or when data or credentials are unavailable.",
        {"type": "object", "properties": {
            "stage": {**_S_STRING}, "status": {**_S_STRING, "description": "pending, running, done, failed or blocked."},
            "evidence": {**_S_ARRAY_S, "description": "Paths this stage produced."},
            "note": {**_S_STRING, "description": "What happened; required for blocked."}},
         "required": ["stage", "status"]},
        cli=["job", "plan-update"], flag_repeat={"evidence"},
    ),

    _T(
        "plan_show",
        "INSIDE A JOB: the current plan with harness-verified stage states.",
        {"type": "object", "properties": {}},
        cli=["job", "plan-show"],
    ),

    _T(
        "delegate_tasks",
        "INSIDE A JOB: run up to 4 focused sub-agents in parallel, each in a "
        "fresh context, and get back each one's concise result and "
        "artefacts. Use for independent sub-questions (look up several "
        "datasets, check several claims, run several independent analyses). "
        "Each task text must be self-contained.",
        {"type": "object", "properties": {"tasks": {**_S_ARRAY_S, "description": "1-4 self-contained task descriptions."}},
         "required": ["tasks"]},
        cli=["job", "delegate"], flag_map={"tasks": "--task"}, flag_repeat={"tasks"},
    ),

    _T(
        "job_wait",
        "INSIDE A JOB: wait (up to timeout_min, default 240) for a detached "
        "raw-pipeline job to finish or for an output file to appear and stop "
        "growing, instead of ending the round. Returns the job's final status "
        "or the file's presence.",
        {"type": "object", "properties": {
            "raw_job": {**_S_STRING, "description": "Detached raw-pipeline job id."},
            "path": {**_S_STRING, "description": "Workspace path to wait for."},
            "timeout_min": {**_S_NUMBER}}},
        cli=["job", "wait"], flag_map={"raw_job": "--raw-job", "timeout_min": "--timeout-min"},
    ),

    # Paper reproduction from the authors' own code (paper_code_skill.py),
    # following Paper2Agent's code-first route.
    _T(
        "paper_code_find",
        "Which GitHub repository holds a paper's analysis code. Give the "
        "harvest.json written by paper_benchmark/bench harvest (the Code "
        "Availability statement is ranked first), a repository, or search text.",
        {"type": "object", "properties": {
            "harvest": {**_S_STRING, "description": "Path to a bench harvest.json."},
            "repo": {**_S_STRING, "description": "owner/name or GitHub URL."},
            "query": {**_S_STRING, "description": "GitHub search text (title, first author)."}}},
        cli=["paper-code", "find"],
    ),

    _T(
        "paper_code_inventory",
        "Clone (pinned) and inventory a paper's code repository: the main "
        "analysis (.Rmd/.qmd/.ipynb/.R/.py), the inputs it reads and whether "
        "they are in the repository, packages it loads, output directories, the "
        "authors' own rendered output (for comparison), the package versions "
        "the authors ran, and the analysis sections.",
        {"type": "object", "properties": {
            "repo": {**_S_STRING, "description": "owner/name or GitHub URL."},
            "entry": {**_S_STRING, "description": "Analysis file (default: the main Rmd/notebook)."}},
         "required": ["repo"]},
        cli=["paper-code", "inventory"], positional=["repo"],
    ),

    _T(
        "paper_code_reproduce",
        "★ REPRODUCE A PAPER BY RUNNING THE AUTHORS' CODE ★ when a paper's Code "
        "Availability names a repository. Clones it at a pinned commit, builds "
        "the R (micromamba) or Python (uv) environment it needs, runs the "
        "analysis UNMODIFIED in a work copy (one failing chunk does not hide the "
        "rest), compares printed values and figures with the authors' own "
        "rendered output, and writes a per-section report with every figure. "
        "Runs in the background: wait with job_wait on <run_dir>/done.json (or "
        "paper_code_status), then read summary.json. Use pins to repair an "
        "environment (e.g. r-reshape2=1.4.4 or cran:reshape2@1.4.4).",
        {"type": "object", "properties": {
            "repo": {**_S_STRING, "description": "owner/name or GitHub URL."},
            "harvest": {**_S_STRING, "description": "Or: pick the repository from this harvest.json."},
            "ref": {**_S_STRING, "description": "Commit or tag to pin (default: HEAD, recorded)."},
            "entry": {**_S_STRING, "description": "Analysis file inside the repository."},
            "pin": {**_S_ARRAY_S, "description": "conda/pip pins or cran:pkg@version."},
            "strict": {**_S_BOOLEAN, "description": "Stop at the first failing chunk."},
            "replay": {**_S_BOOLEAN, "description": "Execute twice and check the outputs are identical."},
            "input": {**_S_ARRAY_S, "description": "Apply the authors' code to new data: NAME=PATH replaces an input "
                                                  "the analysis reads (names from paper_code_inventory)."},
            "paper": {**_S_STRING, "description": "Paper title/DOI for the report."}}},
        cli=["paper-code", "pipeline", "--detach"], positional=["repo"], flag_repeat={"pin", "input"},
        bool_flags={"strict", "replay"},
    ),

    _T(
        "write_text_file",
        "Write (or append) a UTF-8 text file inside the workspace: section "
        "verdicts, notes, reports, small JSON/TSV evidence. Paths under "
        "Docs/, Data/ or Benchmarks/; never secrets, source code or .git; "
        "text files only, up to 5 MB. Returns the path to use as evidence.",
        {"type": "object", "properties": {
            "path": {**_S_STRING, "description": "Workspace-relative path, e.g. Docs/PaperCode/<run>/verdicts.md."},
            "content": {**_S_STRING, "description": "The full text to write."},
            "append": {**_S_BOOLEAN, "description": "Append instead of replacing."}},
         "required": ["path", "content"]},
        cli=["files", "write"], flag_map={"path": "--path", "content": "--content"}, bool_flags={"append"},
    ),

    _T(
        "read_text_file",
        "Read lines of a text file inside the workspace (never secrets).",
        {"type": "object", "properties": {
            "path": {**_S_STRING}, "offset": {**_S_INTEGER, "description": "First line, 0-based."},
            "limit": {**_S_INTEGER, "description": "Lines to return (default 400)."}},
         "required": ["path"]},
        cli=["files", "read"], flag_map={"path": "--path", "offset": "--offset", "limit": "--limit"},
    ),

    _T(
        "reproductions_search",
        "Search the reproduction records: one per paper, kept across attempts, "
        "with the outcome (reproduced / partial / not reproduced), route "
        "(authors' code or public data), agreement with the paper and the "
        "authors' outputs, the verifier's verdict and the report path. Use it to "
        "answer which papers have been reproduced and how well.",
        {"type": "object", "properties": {"query": {**_S_STRING, "description": "Title, DOI, repository, gene, assay or outcome."}}},
        cli=["repro", "list"], flag_map={"query": "--query"},
    ),

    _T(
        "paper_code_status",
        "State of a paper_code_reproduce run directory (running stage, or the "
        "final result with report, summary and figure counts).",
        {"type": "object", "properties": {"run_dir": {**_S_STRING}}, "required": ["run_dir"]},
        cli=["paper-code", "status"], positional=["run_dir"],
    ),

    _T(
        "processed_discover",
        "★ FIND IGVF PORTAL DATA BY TOPIC ★ when the user names a phenotype, "
        "tissue, gene or kind of data rather than an accession (e.g. 'IGVF "
        "data on coronary artery disease / heart', 'element-gene links in "
        "liver', 'screens targeting GATA1'). Matches the words to the terms "
        "the Portal actually uses (phenotype terms, sample terms, gene "
        "symbols, prediction / library / curated set types), then searches "
        "PredictionSets (associated_phenotypes, assessed_genes, virtual "
        "samples), ModelSets, MeasurementSets (samples, targeted genes, "
        "construct library type, assay), AnalysisSets, CuratedSets and "
        "ConstructLibrarySets. Superseded sets are left out. An empty match "
        "is reported with the nearest real values. Writes report.md and "
        "discover.json; follow each accession with portal_lineage.",
        {"type": "object", "properties": {
            "phenotype": {**_S_STRING, "description": "Phenotype or disease, e.g. 'coronary artery disease'."},
            "tissue": {**_S_STRING, "description": "Tissue / cell type / cell line sample term, e.g. heart, liver, K562."},
            "gene": {**_S_STRING, "description": "Gene symbol assessed by predictions or targeted by screens."},
            "prediction_type": {**_S_STRING, "description": "PredictionSet / ModelSet type, e.g. 'element-gene links', 'non-coding variant effects', 'coding variant effects', 'disease associations', 'variant TF binding effects'."},
            "library_type": {**_S_STRING, "description": "Construct library type: 'guide library', 'reporter library', 'editing template library', 'expression vector library'."},
            "assay": {**_S_STRING, "description": "Assay title, e.g. 'MPRA', '10x multiome', 'CRISPR FlowFISH'."},
            "curated_type": {**_S_STRING, "description": "CuratedSet type, e.g. variants, elements, 'guide RNAs', 'training data for predictive models'."},
            "types": {**_S_STRING, "description": "Comma-separated Portal types to search (default: prediction, model, measurement, analysis, curated and construct library sets)."},
            "limit": {**_S_INTEGER, "description": "Results per type (default 25)."},
            "include_superseded": {**_S_BOOLEAN, "description": "Also list superseded sets."}}},
        cli=["processed", "discover"],
        flag_map={"prediction_type": "--prediction-type", "library_type": "--library-type",
                  "curated_type": "--curated-type", "include_superseded": "--include-superseded"},
        bool_flags={"include_superseded"},
    ),

    _T(
        "portal_lineage",
        "★ CALL THIS FIRST FOR ANY IGVF ACCESSION (IGVFDS / IGVFFI / IGVFSM) "
        "BEFORE DOWNLOADING OR REPROCESSING ANYTHING ★ Walks every link the "
        "IGVF Portal records around the accession, in both directions, and "
        "says what has ALREADY been computed. For a MeasurementSet it "
        "follows input_for to intermediate analyses (IGVF uniform pipeline "
        "or lab), the principal analysis built on them (cell annotations), "
        "pseudobulk sets and the prediction sets built on those (e.g. scE2G "
        "element-gene links). It also follows the multiome partner "
        "(related_measurement_sets), auxiliary sets (MULTI-seq / hashing / "
        "guides), the multiplexed sample and its barcode to sample map, and "
        "seqspec / onlist files. For a PredictionSet it climbs to its "
        "ModelSet, the model's training data and the input data. It fetches "
        "the QualityMetric objects the pipeline published, so QC numbers "
        "(pseudoalignment %, reads on onlist, duplicates) come without "
        "processing reads. Returns a START-HERE table with the best "
        "processed file per product (RNA matrix, ATAC matrix, fragments, "
        "peaks, cell annotations, demultiplexing, perturbation effects, "
        "predictions, model, alignments). Each row gives size, access and "
        "source set, set against the raw-read volume. It lists linked "
        "objects the credentials cannot see (403 = unreleased or controlled) "
        "and writes report.md, plan.json, lineage_graph.json and a lineage "
        "figure. On IGVFDS9875NBZW it finds a 3.7 GB h5ad, 5.6 GB public "
        "fragments, cell annotations and a hashing table, instead of 155 GB "
        "of controlled reads. It also reports what the data is about and "
        "how it was designed: superseded sets and their replacements, the "
        "construct library (guide / reporter / editing-template library and "
        "its integrated guide or element tables), the sample tree (sorted "
        "fractions, time points, treatments, CRISPR modifications), a "
        "prediction's phenotypes, assessed genes, cell type and external "
        "training data, the column-definition documents of tabular outputs, "
        "the analysis step and software that made each file, and "
        "publications. Follow with processed_fetch.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING, "description":
                              "Any IGVF accession or @id: measurement / analysis / auxiliary / prediction / model set, file, sample."},
                "want": {**_S_STRING, "description":
                         "Optional comma-separated products to list, e.g. rna_matrix,fragments,cell_annotations,predictions."},
                "depth": {**_S_INTEGER, "default": 4},
                "max_nodes": {**_S_INTEGER, "default": 300},
                "no_plots": {**_S_BOOLEAN, "default": False},
            },
            "required": ["accession"],
        },
        cli=["processed", "lineage"],
        positional=("accession",),
        flag_map={"want": "--want", "depth": "--depth", "max_nodes": "--max-nodes"},
        bool_flags={"no_plots"},
    ),

    _T(
        "processed_fetch",
        "★ DOWNLOAD WHAT IGVF ALREADY COMPUTED, NOT THE RAW READS ★ After "
        "portal_lineage: downloads the best processed file per requested "
        "product (default RNA matrix, ATAC matrix, fragments, peaks, cell "
        "annotations, demultiplexing, predictions) under a total GB budget. "
        "It also saves the Portal's QC metric objects and their attachments "
        "(kb_info.json, barcode summaries). Controlled-access files are "
        "skipped with a reason when no IGVF credentials are configured. "
        "Files land in Data/Processed/<accession>/ with fetch_manifest.json. "
        "Use dry_run=true to see the transfer first.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING},
                "want": {**_S_STRING, "default": "rna_matrix,atac_matrix,fragments,peaks,cell_annotations,demultiplexing,predictions"},
                "max_gb": {**_S_NUMBER, "default": 20.0},
                "dry_run": {**_S_BOOLEAN, "default": False},
                "no_qc": {**_S_BOOLEAN, "default": False},
            },
            "required": ["accession"],
        },
        cli=["processed", "fetch"],
        positional=("accession",),
        flag_map={"want": "--want", "max_gb": "--max-gb"},
        bool_flags={"dry_run", "no_qc"},
    ),

    _T(
        "processed_plan",
        "★ DOWNLOAD THE EXISTING RESULT, OR RUN THE PIPELINE? ★ Answers that "
        "for one product. `want` is matrix, fragments, peaks, alignments or "
        "index. Returns the specific file to fetch when IGVF has already "
        "produced it — preferring output of a COMPLETED uniform pipeline, "
        "which is the version IGVF stands behind — and says to run the "
        "pipeline only when nothing suitable exists. Use this to decide the "
        "route before committing to a long job.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING},
                "want":      {**_S_STRING, "description":
                               "matrix | fragments | peaks | alignments | index",
                               "default": "matrix"},
            },
            "required": ["accession"],
        },
        cli=["processed", "plan"],
        positional=("accession",),
        flag_map={"want": "--want"},
    ),

    _T(
        "matrix_summary",
        "★ WHAT IS ACTUALLY IN THIS h5ad ★ — shape, obs/var columns, layers, "
        "raw, and TRUE whole-matrix statistics: min, max, where the max is, "
        "non-zero count, whether values are integer-valued. The full "
        "reduction is the DEFAULT and walks the sparse matrix in row blocks, "
        "so it never builds a dense copy. Use this instead of writing your "
        "own inspector: a hand-rolled one that peeked at the first 50 rows "
        "reported a maximum of 14 on a matrix whose true maximum was 3,724, "
        "and nothing in its output revealed that it had sampled. If you "
        "genuinely want a quick peek, pass sample_rows — the result is then "
        "labelled `sampled` and every field is named sample_*, so it cannot "
        "be quoted as a whole-matrix figure by accident.",
        {
            "type": "object",
            "properties": {
                "h5ad":        {**_S_STRING, "description": "Path to the .h5ad."},
                "sample_rows": {**_S_INTEGER, "description":
                                 "Examine only the first N rows; the output "
                                 "then says so. Omit for the true statistics."},
                "no_x_stats":  {**_S_BOOLEAN, "description":
                                 "Structure only; do not read matrix values."},
                "backed":      {**_S_BOOLEAN, "description":
                                 "Open backed on disk for a very large file."},
            },
            "required": ["h5ad"],
        },
        cli=["matrix-qc", "summary"],
        flag_map={"h5ad": "--h5ad", "sample_rows": "--sample-rows",
                   "no_x_stats": "--no-x-stats", "backed": "--backed"},
        bool_flags=("no_x_stats", "backed"),
    ),

    _T(
        "table_threshold",
        "★ COUNT ROWS MEETING NUMERIC CONDITIONS, WITH MEDIANS ★ — the "
        "'now filter this QC table and tell me how many pass' step. Reads a "
        "plain TSV/CSV and needs NOTHING initialised: not a warehouse, not a "
        "database, not a new tool. Use it whenever a question asks how many "
        "rows clear some cut-offs and what the medians are before and after. "
        "Rule syntax: \"total_counts>=1000,n_genes_by_counts>=200\". Reports "
        "rows passing as QC-PASSING, never as validated cells — calling cells "
        "takes a cell-calling method, not a threshold.",
        {
            "type": "object",
            "properties": {
                "table":  {**_S_STRING, "description": "TSV or CSV path."},
                "where":  {**_S_STRING, "description":
                            'Comma-separated conditions, e.g. '
                            '"total_counts>=1000,n_genes_by_counts>=200".'},
                "median": {**_S_STRING, "description":
                            "Comma-separated columns to take medians of."},
            },
            "required": ["table", "where"],
        },
        cli=["matrix-qc", "threshold"],
        flag_map={"table": "--table", "where": "--where", "median": "--median"},
    ),

    _T(
        "history_flag",
        "★ MARK A PAST ANSWER WRONG ★ — call this the moment you establish "
        "that a recalled prior result is incorrect: it contradicts the "
        "portal, the literature, a fresh computation, or plain arithmetic. "
        "Prior answers are served automatically to everyone who asks about "
        "that dataset, so leaving a wrong one in place means the next person "
        "inherits it as though it were settled. Flagging stops it being "
        "recalled; it is never deleted and stays searchable. Also takes "
        "'correct' when you have verified one holds up. Give a reason — the "
        "next reader needs to know WHAT was wrong.",
        {
            "type": "object",
            "properties": {
                "ref":     {**_S_STRING, "description":
                             "Run directory, e.g. Docs/Agent/20260914_150552_…"},
                "verdict": {**_S_STRING, "description":
                             "wrong | correct | unsure", "default": "wrong"},
                "reason":  {**_S_STRING, "description":
                             "What is wrong with it, specifically."},
            },
            "required": ["ref", "reason"],
        },
        cli=["project", "flag"],
        positional=("ref",),
        flag_map={"verdict": "--verdict", "reason": "--reason"},
    ),

    _T(
        "history_recall",
        "★ WHAT HAVE WE ALREADY PRODUCED ABOUT THIS ACCESSION ★ — CALL THIS "
        "FIRST whenever a question names an IGVF or ENCODE accession. Returns "
        "every past agent session, skill run and download that ever touched "
        "it, newest first, with the run directory holding the figures, tables "
        "and report. A dataset that was analysed last week does not need "
        "analysing again: read the recorded result and say where it came "
        "from. Returns nothing when the accession is genuinely new, which is "
        "the signal to go and do the work.",
        {
            "type": "object",
            "properties": {
                "accession": {**_S_STRING, "description":
                               "IGVFDS… / IGVFFI… / ENCSR… — one accession."},
                "limit":     {**_S_INTEGER, "default": 50},
            },
            "required": ["accession"],
        },
        cli=["project", "recall"],
        positional=("accession",),
        flag_map={"limit": "--limit"},
    ),

    _T(
        "history_search",
        "★ SEARCH EVERY PAST ANALYSIS BY WORDS ★ — full-text over the "
        "questions asked, the answers given, the skill runs executed and the "
        "project notes written. Use it when the user refers to earlier work "
        "without an accession (\"the spatial ATAC paper we did\", \"that "
        "CRISPR screen\"), or before starting something that may already have "
        "been done. Each hit carries a run directory that history_show "
        "replays in full.",
        {
            "type": "object",
            "properties": {
                "query":   {**_S_STRING, "description": "Free text keywords."},
                "kind":    {**_S_STRING, "description":
                             "Restrict to session | analysis | download | "
                             "project | item."},
                "project": {**_S_STRING, "description":
                             "Restrict to one project by name or id."},
                "limit":   {**_S_INTEGER, "default": 20},
            },
            "required": ["query"],
        },
        cli=["project", "search"],
        positional=("query",),
        flag_map={"kind": "--kind", "project": "--project", "limit": "--limit"},
    ),

    _T(
        "history_show",
        "★ REPLAY ONE RECORDED SESSION ★ — the original question, the answer "
        "as given, the accessions involved and every artefact path it "
        "produced. Takes a run directory from history_recall or "
        "history_search. This is how a past result is quoted accurately "
        "instead of paraphrased from memory.",
        {
            "type": "object",
            "properties": {
                "ref":   {**_S_STRING, "description":
                           "Run directory, e.g. Docs/Agent/20260914_150552_…"},
                "brief": {**_S_BOOLEAN, "description":
                           "Omit the answer text; metadata only."},
            },
            "required": ["ref"],
        },
        cli=["project", "show"],
        positional=("ref",),
        flag_map={"brief": "--brief"},
        bool_flags={"brief"},
    ),

    _T(
        "project_manage",
        "★ CREATE / LIST / RENAME A PROJECT, OR FILE WORK INTO ONE ★ — a "
        "project is a permanent named container for analyses. Subcommands: "
        "create (needs name), list, use (make active — later runs are filed "
        "automatically), rename (needs name + to; old references still "
        "resolve), add (needs kind + ref), items, describe, archive, stats. "
        "Use it when the user says they are starting a study, or asks to keep "
        "or group results. Nothing a project holds is ever deleted.",
        {
            "type": "object",
            "properties": {
                "subcommand": {**_S_STRING, "description":
                                "create | list | use | rename | add | items | "
                                "describe | archive | stats | recent"},
                "name":        {**_S_STRING, "description":
                                 "Project name or id, for every subcommand "
                                 "that names one."},
                "to":          {**_S_STRING, "description": "New name, for rename."},
                "description": {**_S_STRING},
                "kind":        {**_S_STRING, "description":
                                 "For add: session | analysis | artifact | "
                                 "dataset | figure | paper | note."},
                "ref":         {**_S_STRING, "description":
                                 "For add: run directory, accession or path."},
                "title":       {**_S_STRING},
                "project":     {**_S_STRING, "description":
                                 "For add/items: target project; defaults to "
                                 "the active one."},
            },
            "required": ["subcommand"],
        },
        cli=["project"],
        positional=("subcommand", "name"),
        flag_map={"to": "--to", "description": "--description", "kind": "--kind",
                   "ref": "--ref", "title": "--title", "project": "--project"},
    ),

    _T(
        "regulome_annotate",
        "★ IS THIS NON-CODING VARIANT IN A REGULATORY ELEMENT, AND IN WHICH "
        "TISSUES ★ — RegulomeDB rank (1a best .. 7 no evidence), probability, "
        "the evidence behind it (ChIP-seq, open chromatin, footprint, motif, "
        "QTL), and a PER-TISSUE score across organs. Answers a DIFFERENT "
        "question from the other annotators: FAVOR/ccre_favor give "
        "conservation and deleteriousness (CADD, GERP), ClinVar gives "
        "clinical significance — none say whether the base sits in something "
        "regulatory, or where. Accepts chr1:39492462 or chr1-39492462-A-G. "
        "The API returns ~1.3 MB per variant; this extracts the score and "
        "tissue table and COUNTS the supporting experiments rather than "
        "returning them.",
        {
            "type": "object",
            "properties": {
                "variants": {**_S_STRING, "description":
                              "Comma/space separated, e.g. chr1:39492462."},
                "input":    {**_S_STRING, "description":
                              "Or a file with one variant per line."},
                "assembly": {**_S_STRING, "default": "GRCh38"},
                "label":    {**_S_STRING},
            },
        },
        cli=["regulome", "annotate"],
        flag_map={"variants": "--variants", "input": "--input",
                   "assembly": "--assembly", "label": "--label"},
    ),

    _T(
        "bcalm_activity",
        "★ MPRA ACTIVITY FROM BARCODE-LEVEL COUNTS ★ (BCalm approach). Use "
        "INSTEAD OF mpra_activity when you have per-barcode counts and care "
        "about confidence, not just effect size. mpra_activity SUMS an "
        "element's barcodes before testing, which throws away the spread "
        "among them — the best evidence about how noisy that measurement "
        "really is. Two elements with the same mean activity, one whose "
        "barcodes agree and one whose barcodes disagree, come out of a summed "
        "analysis with the SAME confidence; here they do not. Per-element "
        "variances are moderated by empirical Bayes, so an element with few "
        "barcodes is shrunk toward the trend rather than trusted on its own "
        "noisy estimate. Needs DNA*/RNA* count columns per replicate.",
        {
            "type": "object",
            "properties": {
                "counts":     {**_S_STRING, "description":
                                "Barcode-level table: oligo, barcode, "
                                "DNA_rep1/RNA_rep1, ... (from "
                                "mpraflow barcode-matrix or equivalent)."},
                "oligo_col":  {**_S_STRING},
                "barcode_col": {**_S_STRING},
                "min_obs":    {**_S_INTEGER, "default": 3},
                "fdr":        {**_S_NUMBER, "default": 0.05},
                "label":      {**_S_STRING},
            },
            "required": ["counts"],
        },
        cli=["bcalm", "activity"],
        flag_map={"counts": "--counts", "oligo_col": "--oligo-col",
                   "barcode_col": "--barcode-col", "min_obs": "--min-obs",
                   "fdr": "--fdr", "label": "--label"},
    ),

    _T(
        "scnt_seq_count",
        "★ NEW vs OLD RNA FROM 4sU METABOLIC LABELLING ★ (scNT-seq). Counts "
        "T>C conversions per read in an aligned BAM and splits new from old "
        "RNA per cell and per gene. NOT the same as kb's `nac` workflow, "
        "which IGVFagent also has: `nac` calls a transcript nascent from "
        "INTRON content, this calls it new from CHEMICAL LABELLING — an "
        "intronless transcript made an hour ago is new here and mature there. "
        "Needs a BAM with MD tags (conversions are read from them) and cell "
        "barcodes. Reports the rate of every OTHER substitution type as "
        "background, and warns when labelling conversions are not clearly "
        "above it — a low ratio means the 'new' calls are sequencing error.",
        {
            "type": "object",
            "properties": {
                "bam":             {**_S_STRING, "description":
                                     "Aligned BAM with MD tags."},
                "cell_tag":        {**_S_STRING, "default": "CB"},
                "gene_tag":        {**_S_STRING, "default": "GX"},
                "min_conversions": {**_S_INTEGER, "default": 2,
                                     "description": "Conversions needed to "
                                     "call a read new. 1 is weak — SNPs and "
                                     "sequencing error also make T>C."},
                "max_reads":       {**_S_INTEGER},
                "label":           {**_S_STRING},
            },
            "required": ["bam"],
        },
        cli=["scnt-seq", "count"],
        flag_map={"bam": "--bam", "cell_tag": "--cell-tag",
                   "gene_tag": "--gene-tag",
                   "min_conversions": "--min-conversions",
                   "max_reads": "--max-reads", "label": "--label"},
    ),

    # ---- scE2G workbench (EngreitzLab/scE2G training + CRISPR_comparison) ----

    _T(
        "sce2g_setup",
        "★ scE2G TRAINING, STEP 1: PREPARE THE CHECKOUT ★. Clones "
        "EngreitzLab/scE2G (branch fix/dag-staleness-integration) with the "
        "ENCODE-rE2G submodule, or takes an existing checkout, and applies the "
        "crowdsourced-features walkthrough's patches: drop `conda: \"mamba\"` "
        "from both Snakefile_training files, add SCRIPTS_DIR, add "
        "RNA_matrix_filtered / max_cell_count to config_training.yaml, put the "
        "missing resources/feature_tables/multiome_arc_n6.tsv in place. "
        "Idempotent. `check=true` only reports.",
        {"type": "object", "properties": {
            "repo_dir": {**_S_STRING, "description": "Path of the scE2G checkout (created if absent)."},
            "branch":   {**_S_STRING, "default": "fix/dag-staleness-integration"},
            "check":    {**_S_BOOLEAN, "default": False},
            "keep_branch": {**_S_BOOLEAN, "default": False},
            "offline":  {**_S_BOOLEAN, "default": False}},
         "required": ["repo_dir"]},
        cli=["sce2g", "setup"],
        flag_map={"repo_dir": "--repo-dir", "branch": "--branch"},
        bool_flags={"check", "keep_branch", "offline"},
    ),

    _T(
        "sce2g_features",
        "★ scE2G TRAINING, STEP 2: ADD CROWDSOURCED FEATURES ★. Converts an "
        "E2G feature table (ElementChr/ElementStart/ElementEnd/GeneSymbol + "
        "feature columns, as shared on Synapse, e.g. syn73717888 for K562) "
        "into the scE2G source_file (chr/start/end/TargetGene, spaces -> "
        "underscores, .tsv.gz), writes config/external_features_config_<name>"
        ".tsv (input_col, source_col, aggregate_function=mean, join_by=overlap, "
        "source_file) and resources/feature_tables/multiome_arc_n6_<name>.tsv "
        "extending the base table with one row per feature (max, fill 0, "
        "nice_name).",
        {"type": "object", "properties": {
            "features": {**_S_STRING, "description": "Feature table TSV."},
            "name":     {**_S_STRING, "description": "Feature-set name used in every file name."},
            "repo_dir": {**_S_STRING},
            "out_dir":  {**_S_STRING},
            "select":   {**_S_STRING, "description": "Comma list of feature columns to keep."},
            "merge_aggregate": {**_S_STRING, "default": "mean"},
            "benchmark_aggregate": {**_S_STRING, "default": "max"},
            "fill_value": {**_S_STRING, "default": "0"}},
         "required": ["features", "name"]},
        cli=["sce2g", "features"],
        flag_map={"features": "--features", "name": "--name", "repo_dir": "--repo-dir",
                   "out_dir": "--out-dir", "select": "--select",
                   "merge_aggregate": "--merge-aggregate",
                   "benchmark_aggregate": "--benchmark-aggregate", "fill_value": "--fill-value"},
    ),

    _T(
        "sce2g_configure",
        "★ scE2G TRAINING, STEP 3: WRITE THE CONFIG ROWS ★. Adds or replaces "
        "the cluster row in config/config_cell_clusters.tsv (rna_matrix_file, "
        "atac_frag_file, model_dir=models/multiome_powerlaw_v3, plus the "
        "external_features_config column) and the model row in "
        "config/config_models.tsv (dataset == cluster, ABC_directory blank, "
        "polynomial False, feature_table from step 2), and writes "
        "run_training_<model>.sh with the Slurm-profile snakemake command and "
        "the memory/runtime overrides from the walkthrough.",
        {"type": "object", "properties": {
            "repo_dir": {**_S_STRING}, "cluster": {**_S_STRING},
            "rna": {**_S_STRING, "description": "rna_matrix_file (.csv.gz)."},
            "atac_frag": {**_S_STRING, "description": "atac_frag_file (.tsv.gz)."},
            "model_dir": {**_S_STRING, "default": "models/multiome_powerlaw_v3"},
            "name": {**_S_STRING, "description": "Feature-set name from sce2g_features."},
            "model": {**_S_STRING}, "feature_table": {**_S_STRING},
            "profile": {**_S_STRING, "description": "Snakemake profile dir, e.g. profiles/slurm."},
            "jobs": {**_S_INTEGER, "default": 4},
            "hic": {**_S_STRING}, "hic_type": {**_S_STRING}, "hic_resolution": {**_S_STRING}},
         "required": ["repo_dir", "cluster", "atac_frag"]},
        cli=["sce2g", "configure"],
        flag_map={"repo_dir": "--repo-dir", "cluster": "--cluster", "rna": "--rna",
                   "atac_frag": "--atac-frag", "model_dir": "--model-dir", "name": "--name",
                   "model": "--model", "feature_table": "--feature-table", "profile": "--profile",
                   "jobs": "--jobs", "hic": "--hic", "hic_type": "--hic-type",
                   "hic_resolution": "--hic-resolution"},
    ),

    _T(
        "sce2g_check",
        "★ scE2G TRAINING, STEP 4: CHECK BEFORE SPENDING COMPUTE ★. Verifies "
        "every path the configs reference exists, cluster == dataset, the "
        "feature table lists every external input_col, the source_file has "
        "chr/start/end/TargetGene and the source_cols, no feature name has a "
        "space, the Snakefile patches are in place, and prints the snakemake "
        "command. Exit 1 with the list of problems.",
        {"type": "object", "properties": {
            "repo_dir": {**_S_STRING}, "model": {**_S_STRING},
            "profile": {**_S_STRING}, "jobs": {**_S_INTEGER, "default": 4}},
         "required": ["repo_dir"]},
        cli=["sce2g", "check"],
        flag_map={"repo_dir": "--repo-dir", "model": "--model", "profile": "--profile", "jobs": "--jobs"},
    ),

    _T(
        "sce2g_run",
        "★ scE2G TRAINING, STEP 5: RUN SNAKEMAKE ★ (dry-run unless "
        "execute=true). Runs sce2g_check first and refuses on problems. Needs "
        "snakemake on PATH; training takes hours and is normally run with a "
        "Slurm profile on a cluster, not on the hosted instance.",
        {"type": "object", "properties": {
            "repo_dir": {**_S_STRING}, "model": {**_S_STRING}, "profile": {**_S_STRING},
            "jobs": {**_S_INTEGER, "default": 4},
            "execute": {**_S_BOOLEAN, "default": False},
            "force": {**_S_BOOLEAN, "default": False}},
         "required": ["repo_dir"]},
        cli=["sce2g", "run"],
        flag_map={"repo_dir": "--repo-dir", "model": "--model", "profile": "--profile", "jobs": "--jobs"},
        bool_flags={"execute", "force"},
    ),

    _T(
        "sce2g_predictions",
        "★ DESCRIBE / COMPARE scE2G PREDICTION TABLES ★ (*.e2g.tsv, or any "
        "element-gene table with ElementChr/Start/End + GeneSymbol or chr/start/"
        "end + TargetGene). Links, genes, elements, score quantiles, links above "
        "the threshold, element classes, self-promoter links; with two or more "
        "tables also shared pairs, Jaccard and Spearman of scores on shared "
        "pairs. Score-distribution figure.",
        {"type": "object", "properties": {
            "predictions": {**_S_ARRAY_S, "description": "LABEL=PATH entries."},
            "score_col": {**_S_STRING, "default": "E2G.Score.qnorm"},
            "threshold": {**_S_NUMBER, "default": 0.177},
            "label": {**_S_STRING}, "no_plots": {**_S_BOOLEAN, "default": False}},
         "required": ["predictions"]},
        cli=["sce2g", "predictions"],
        flag_map={"predictions": "--predictions", "score_col": "--score-col",
                   "threshold": "--threshold", "label": "--label"},
        flag_repeat={"predictions"}, bool_flags={"no_plots"},
    ),

    _T(
        "sce2g_benchmark",
        "★ CRISPR BENCHMARK OF E2G PREDICTIONS ★ in the manner of "
        "EngreitzLab/CRISPR_comparison: each CRISPR-tested element-gene pair "
        "(EPCrisprBenchmark: chrom, chromStart, chromEnd, measuredGeneSymbol, "
        "Regulated) gets the aggregate (max/mean/sum) of overlapping predicted "
        "elements for the same gene, fill_value when none, inverse_predictor "
        "and boolean semantics from a pred_config.txt; AUPRC (step rule) plus "
        "`auprc_crispr_comparison` (the upstream pipeline's exact definition, "
        "validated to 4 decimals against it), bootstrap 95% interval, "
        "precision at 70% recall, PR-curve figure. Also writes "
        "pred_config.txt + config.yml so the upstream Snakemake pipeline can "
        "be run on the same inputs. Use to compare a newly trained scE2G "
        "model against the base model or ABC.",
        {"type": "object", "properties": {
            "predictions": {**_S_ARRAY_S, "description": "LABEL=PATH entries."},
            "crispr": {**_S_STRING, "description": "EPCrisprBenchmark TSV."},
            "pred_config": {**_S_STRING},
            "score_col": {**_S_STRING, "default": "E2G.Score.qnorm"},
            "cell_type": {**_S_STRING},
            "bootstrap": {**_S_INTEGER, "default": 200},
            "all_features": {**_S_BOOLEAN, "default": False, "description":
                "Benchmark EVERY feature column of each table as its own "
                "predictor -- the crowdsourced-feature benchmark. Tables are "
                "streamed and restricted to CRISPR-tested genes, so 11M-row "
                "Synapse tables are fine."},
            "gene_universe_filter": {**_S_BOOLEAN, "default": False, "description":
                "CRISPR_comparison semantics: drop tested pairs whose gene is "
                "absent from the prediction table instead of scoring them as "
                "fill_value (raises scE2G K562 AUPRC from 0.531 to 0.551)."},
            "tss_bed": {**_S_STRING, "description": "TSS universe BED (e.g. CollapsedGeneBounds.hg38.TSS500bp.bed); needed for baselines and distance bins."},
            "gene_bed": {**_S_STRING, "description": "Gene-body BED for distToGene / nearestGene / within100kbGene."},
            "expressed_genes": {**_S_STRING, "description": "TSV cell_type, gene, expressed for the *Expr* baselines."},
            "baselines": {**_S_ARRAY_S, "description": "CRISPR_comparison baseline predictors: distToTSS, distToGene, nearestTSS, nearestGene, within100kbTSS, within100kbGene, nearestExprTSS, nearestExprGene, within100kbExprTSS, within100kbExprGene (validated: distToTSS K562 AUPRC 0.4359 = upstream)."},
            "filter_pred_tss": {**_S_STRING, "description": "TSS BED: drop predicted elements overlapping a gene TSS (upstream default filter_pred_tss: True)."},
            "dist_bins_kb": {"type": "array", "items": {"type": "number"}, "description": "Distance-to-TSS bin edges in kb for AUPRC by bin, e.g. [0, 20, 100, 2500]."},
            "delta": {**_S_ARRAY_S, "description": "PRED1,PRED2 pairs: bootstrap delta AUPRC with CI and p-value on shared pairs."},
            "label": {**_S_STRING}, "no_plots": {**_S_BOOLEAN, "default": False}},
         "required": ["predictions", "crispr"]},
        cli=["sce2g", "benchmark"],
        flag_map={"predictions": "--predictions", "crispr": "--crispr", "pred_config": "--pred-config",
                   "score_col": "--score-col", "cell_type": "--cell-type", "bootstrap": "--bootstrap",
                   "tss_bed": "--tss-bed", "gene_bed": "--gene-bed", "expressed_genes": "--expressed-genes",
                   "baselines": "--baselines", "filter_pred_tss": "--filter-pred-tss", "dist_bins_kb": "--dist-bins-kb",
                   "delta": "--delta", "label": "--label"},
        flag_repeat={"predictions", "baselines", "dist_bins_kb", "delta"}, bool_flags={"no_plots", "all_features", "gene_universe_filter"},
    ),

    _T(
        "sce2g_inventory",
        "★ INVENTORY A FOLDER OF CROWDSOURCED E2G FEATURE TABLES ★ (e.g. the "
        "Synapse K562 folder syn73717888 after `synapse download`): rows, "
        "columns, feature names, missing fraction per feature, and whether "
        "each table shares the first table's element-gene universe. Run "
        "before sce2g_features / sce2g_benchmark to see what you have.",
        {"type": "object", "properties": {
            "dir": {**_S_STRING, "description": "Folder of *.tsv.gz feature tables."},
            "max_rows": {**_S_INTEGER, "default": 0, "description": "Rows per file to scan (0 = all)."},
            "label": {**_S_STRING}},
         "required": ["dir"]},
        cli=["sce2g", "inventory"],
        flag_map={"dir": "--dir", "max_rows": "--max-rows", "label": "--label"},
    ),

    _T(
        "sce2g_benchmark_merge",
        "Merge several sce2g_benchmark run directories into one ranked table "
        "(AUPRC, bootstrap interval, precision at 70% recall, AUPRC over "
        "random, best orientation) and one ranked-bar figure. Use after "
        "benchmarking feature tables in batches.",
        {"type": "object", "properties": {
            "runs": {**_S_ARRAY_S, "description": "Benchmark run directories."},
            "top": {**_S_INTEGER, "default": 30},
            "label": {**_S_STRING}, "no_plots": {**_S_BOOLEAN, "default": False}},
         "required": ["runs"]},
        cli=["sce2g", "merge"],
        flag_map={"runs": "--runs", "top": "--top", "label": "--label"},
        flag_repeat={"runs"}, bool_flags={"no_plots"},
    ),

    _T(
        "eqtl_enrichment_setup",
        "Fetch the EngreitzLab/eQTLEnrichment resource files (hg38 partition, TSS "
        "reference / gene universe, gene bounds, chromosome sizes) into "
        "Data/eQTLEnrichment/resources; with synapse=true also the 1000G "
        "background SNPs (syn52264319) and the GTEx SuSiE fine-mapping release "
        "(syn52264297) via the synapse skill (needs SYNAPSE_AUTH_TOKEN).",
        {"type": "object", "properties": {
            "synapse": {**_S_BOOLEAN, "default": False},
            "force": {**_S_BOOLEAN, "default": False}}},
        cli=["eqtl-enrich", "setup"], bool_flags={"synapse", "force"},
    ),

    _T(
        "eqtl_enrichment_prepare_gtex",
        "Raw GTEx fine-mapping table (19 columns: GTEx_30tissues_release1.tsv.gz) -> "
        "the 7-column eQTL input for eqtl_enrichment_run: SUSIE rows in a credible "
        "set, Ensembl -> HGNC through the gene bounds table, and (with expression) "
        "eGene median TPM > tpm in its tissue from the GTEx .gct.",
        {"type": "object", "properties": {
            "raw": _S_STRING, "out": _S_STRING, "method": {**_S_STRING, "default": "SUSIE"},
            "expression": {**_S_STRING, "description": "GTEx median-TPM .gct(.gz)"},
            "tpm": {**_S_NUMBER, "default": 1.0}, "keep_all_cs": {**_S_BOOLEAN, "default": False}},
         "required": ["raw", "out"]},
        cli=["eqtl-enrich", "prepare-gtex"],
        flag_map={"raw": "--raw", "out": "--out", "method": "--method", "expression": "--expression", "tpm": "--tpm"},
        bool_flags={"keep_all_cs"},
    ),

    _T(
        "eqtl_enrichment_variants",
        "Filter fine-mapped eQTL variants (PIP >= threshold, distal noncoding via the "
        "hg38 partition, eGene in the gene universe, eVariant-eGene TSS distance bins) "
        "and the background SNPs (distal noncoding) once, writing a variants dir that "
        "eqtl_enrichment_run reuses with variants_dir.",
        {"type": "object", "properties": {
            "eqtl": {**_S_STRING, "description": "chr,start,end,varID_hg38,gene_hgnc,tissue,pip"},
            "bg_variants": {**_S_STRING, "description": "1000G SNP bed: chr,start,end,rsid"},
            "threshold_pip": {**_S_NUMBER, "default": 0.5},
            "distances": {"type": "array", "items": _S_INTEGER, "description": "distance bin edges in kb (default 10 100 250 1000)"},
            "label": _S_STRING},
         "required": ["eqtl", "bg_variants"]},
        cli=["eqtl-enrich", "variants"],
        flag_map={"eqtl": "--eqtl", "bg_variants": "--bg-variants", "threshold_pip": "--threshold-pip",
                  "distances": "--distances", "label": "--label"},
        flag_repeat={"distances"},
    ),

    _T(
        "eqtl_enrichment_run",
        "★ eQTL ENRICHMENT BENCHMARK OF ENHANCER-GENE PREDICTIONS ★ (port of "
        "EngreitzLab/eQTLEnrichment, the ENCODE-rE2G / scE2G papers' eQTL "
        "benchmark). Inputs: a methods table (method, boolean, inverse_predictor, "
        "pred_name_long, threshold, score_col, color) and a predictions table "
        "(biosample, one column per method with the prediction file, GTExTissue = "
        "comma-separated matched eQTL tissues); prediction files need chr, start, end, "
        "TargetGene and the score column. Computes, per method and (eQTL tissue x "
        "prediction biosample): enrichment of PIP-filtered distal-noncoding eQTL "
        "variants in predicted enhancers vs 1000G SNPs (log-RR CI, hypergeometric "
        "p, Bonferroni), recall (total) and recall (linking to the eGene) across a "
        "quantile threshold span and by eVariant-eGene distance bin, aggregated "
        "all_matches rows, enrichment-recall curves, enrichment at target recalls "
        "with pairwise z-tests, enhancer set sizes, heatmap matrices, figures, "
        "report.md and summary.json. Two upstream quirks (SE formula; unthresholded "
        "by-distance background counts) are corrected unless upstream_compat.",
        {"type": "object", "properties": {
            "methods_table": _S_STRING, "predictions_table": _S_STRING,
            "methods": {**_S_ARRAY_S, "description": "subset of methods to run"},
            "eqtl": _S_STRING, "bg_variants": _S_STRING,
            "variants_dir": {**_S_STRING, "description": "output dir of eqtl_enrichment_variants (skips variant filtering)"},
            "threshold_pip": {**_S_NUMBER, "default": 0.5},
            "distances": {"type": "array", "items": _S_INTEGER},
            "recalls": {"type": "array", "items": _S_NUMBER, "description": "target recall (linking) values, default 0.03 0.05 0.1"},
            "threshold_pval": {**_S_NUMBER, "default": 0.05},
            "n_threshold_steps": {**_S_INTEGER, "default": 50},
            "upstream_compat": {**_S_BOOLEAN, "default": False},
            "no_plots": {**_S_BOOLEAN, "default": False},
            "label": _S_STRING},
         "required": ["methods_table", "predictions_table"]},
        cli=["eqtl-enrich", "run"],
        flag_map={"methods_table": "--methods-table", "predictions_table": "--predictions-table", "methods": "--methods",
                  "eqtl": "--eqtl", "bg_variants": "--bg-variants", "variants_dir": "--variants-dir",
                  "threshold_pip": "--threshold-pip", "distances": "--distances", "recalls": "--recalls",
                  "threshold_pval": "--threshold-pval", "n_threshold_steps": "--n-threshold-steps", "label": "--label"},
        flag_repeat={"methods", "distances", "recalls"}, bool_flags={"upstream_compat", "no_plots"},
    ),

    _T(
        "eqtl_enrichment_selftest",
        "Self-test of the eQTL enrichment benchmark on a synthetic genome with planted "
        "enhancers, eQTLs and background SNPs: overlap engine, statistics, prepare-gtex, "
        "variants and a four-predictor run (good, random, inverse distance, binary) "
        "all asserted.",
        {"type": "object", "properties": {"no_plots": {**_S_BOOLEAN, "default": False}}},
        cli=["eqtl-enrich", "selftest"], bool_flags={"no_plots"},
    ),

    _T(
        "igvf_crispr_parse_seqspec",
        "Parses a seqspec YAML into the IGVF CRISPR pipeline's <modality>_parsed_seqSpec.txt: the reads of each modality, the kallisto-bustools technology string (the `seqspec index -t kb` equivalent, e.g. 0,0,16:0,16,26:1,23,43 for a CROP-seq guide library) and the barcode onlist path. Reads the YAML directly in Python, or uses the seqspec binary when installed. Writes the TSV, report.md and summary.json with per-read region coordinates.",
        {"type": "object", "properties": {
            "yaml": {**_S_STRING, "description": "seqspec YAML path"},
            "modality": {**_S_ARRAY_S, "description": "modalities to parse, e.g. rna, guide, multiseq"},
            "directory": {**_S_STRING, "description": "seqspec directory joined to the onlist filename (default: the YAML's directory)"},
            "engine": {**_S_STRING, "default": "auto", "description": "auto | python | seqspec"},
            "strict_region_types": {**_S_BOOLEAN, "default": False, "description": "only region_type 'barcode' counts as a barcode (seqspec binary behaviour)"},
            "label": {**_S_STRING, "default": "seqspec", "description": "run label"}},
         "required": ["yaml", "modality"]},
        cli=["igvf-crispr-pipeline", "parse-seqspec"],
        flag_map={"yaml": "--yaml", "modality": "--modality", "directory": "--directory", "engine": "--engine", "strict_region_types": "--strict-region-types", "label": "--label"},
        flag_repeat={'modality'},
        bool_flags={'strict_region_types'},
    ),

    _T(
        "igvf_crispr_guide_features",
        "Converts a guide library table (xlsx, tsv or csv with sgRNA_ID and sgRNA_sequences columns) into guide_features.txt (sequence<TAB>ID, no header), the input of kb ref --workflow kite. Reports duplicate sequences/IDs and guide lengths.",
        {"type": "object", "properties": {
            "guide_table": {**_S_STRING, "description": "guide table path"},
            "out": {**_S_STRING, "description": "explicit output path (skips the run dir/report)"},
            "label": {**_S_STRING, "default": "guide_features", "description": "run label"}},
         "required": ["guide_table"]},
        cli=["igvf-crispr-pipeline", "guide-features"],
        flag_map={"guide_table": "--guide-table", "out": "--out", "label": "--label"},
    ),

    _T(
        "igvf_crispr_kb_ref",
        "Builds (or, without kb on PATH, prints) the kb ref command of the IGVF CRISPR pipeline: a prebuilt transcriptome index (kb ref -d human) or a kite guide index from guide_features.txt. Exits 0 with the exact command when kb is not installed.",
        {"type": "object", "properties": {
            "mode": {**_S_STRING, "default": "transcriptome", "description": "transcriptome | kite"},
            "species": {**_S_STRING, "default": "human", "description": "kb ref -d species"},
            "guide_features": {**_S_STRING, "description": "guide_features.txt (kite)"},
            "guide_table": {**_S_STRING, "description": "guide table; guide_features.txt is written from it (kite)"},
            "genome": {**_S_STRING, "description": "genome FASTA passed as -f1 only with upstream_compat"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "pass the genome as -f1 like the upstream"},
            "dry_run": {**_S_BOOLEAN, "default": False, "description": "print the command only"},
            "label": {**_S_STRING, "default": "kb_ref", "description": "run label"}}},
        cli=["igvf-crispr-pipeline", "kb-ref"],
        flag_map={"mode": "--mode", "species": "--species", "guide_features": "--guide-features", "guide_table": "--guide-table", "genome": "--genome", "upstream_compat": "--upstream-compat", "dry_run": "--dry-run", "label": "--label"},
        bool_flags={'dry_run', 'upstream_compat'},
    ),

    _T(
        "igvf_crispr_kb_count",
        "Runs (or prints) kb count for RNA or guide reads with the -x technology string and -w whitelist taken from a parsed seqspec, exactly as the upstream mappingscRNA / mappingGuide processes (6 threads, -m 48G). Output ks_transcripts_out or ks_guide_out with counts_unfiltered/adata.h5ad.",
        {"type": "object", "properties": {
            "modality": {**_S_STRING, "default": "rna", "description": "rna | guide"},
            "index": {**_S_STRING, "description": "kallisto index"},
            "t2g": {**_S_STRING, "description": "transcript-to-gene / t2guide file"},
            "parsed_seqspec": {**_S_STRING, "description": "<modality>_parsed_seqSpec.txt"},
            "technology": {**_S_STRING, "description": "explicit kb -x string"},
            "whitelist": {**_S_STRING, "description": "barcode onlist"},
            "fastqs": {**_S_ARRAY_S, "description": "FASTQs in kb order R1 R2 ..."},
            "fastq_dir": {**_S_STRING, "description": "directory for reads1/reads2 lists"},
            "reads1": {**_S_STRING, "description": "';'-separated read-1 files"},
            "reads2": {**_S_STRING, "description": "';'-separated read-2 files"},
            "threads": {**_S_INTEGER, "default": 6, "description": "threads"},
            "memory": {**_S_STRING, "default": "48G", "description": "kb -m"},
            "dry_run": {**_S_BOOLEAN, "default": False, "description": "print the command only"},
            "label": {**_S_STRING, "default": "kb_count", "description": "run label"}},
         "required": ["index", "t2g"]},
        cli=["igvf-crispr-pipeline", "kb-count"],
        flag_map={"modality": "--modality", "index": "--index", "t2g": "--t2g", "parsed_seqspec": "--parsed-seqspec", "technology": "--technology", "whitelist": "--whitelist", "fastqs": "--fastqs", "fastq_dir": "--fastq-dir", "reads1": "--reads1", "reads2": "--reads2", "threads": "--threads", "memory": "--memory", "dry_run": "--dry-run", "label": "--label"},
        flag_repeat={'fastqs'},
        bool_flags={'dry_run'},
    ),

    _T(
        "igvf_crispr_count_features",
        "Counts guides per cell barcode from paired FASTQs in pure Python (a stand-in for kb count --workflow kite): splits reads by the seqspec technology string, matches the protospacer exactly or by a unique 1-mismatch neighbour, corrects barcodes against the onlist, and counts distinct UMIs. Writes counts_unfiltered/adata.h5ad, feature_totals.tsv, a knee plot, report.md and summary.json.",
        {"type": "object", "properties": {
            "fastqs": {**_S_ARRAY_S, "description": "FASTQs in technology order (R1 R2 ...), .gz allowed"},
            "guide_table": {**_S_STRING, "description": "guide table with sgRNA_ID/sgRNA_sequences"},
            "guide_features": {**_S_STRING, "description": "guide_features.txt instead of a table"},
            "seqspec": {**_S_STRING, "description": "guide seqspec YAML (gives -x and onlist)"},
            "technology": {**_S_STRING, "description": "explicit kb -x string"},
            "parsed_seqspec": {**_S_STRING, "description": "parsed seqspec TSV"},
            "whitelist": {**_S_STRING, "description": "barcode onlist (plain or .gz)"},
            "no_whitelist": {**_S_BOOLEAN, "default": False, "description": "keep barcodes uncorrected"},
            "exact": {**_S_BOOLEAN, "default": False, "description": "exact guide matches only"},
            "both_strands": {**_S_BOOLEAN, "default": False, "description": "also try the reverse complement"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "count_features", "description": "run label"}},
         "required": ["fastqs"]},
        cli=["igvf-crispr-pipeline", "count-features"],
        flag_map={"fastqs": "--fastqs", "guide_table": "--guide-table", "guide_features": "--guide-features", "seqspec": "--seqspec", "technology": "--technology", "parsed_seqspec": "--parsed-seqspec", "whitelist": "--whitelist", "no_whitelist": "--no-whitelist", "exact": "--exact", "both_strands": "--both-strands", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'fastqs'},
        bool_flags={'no_plots', 'both_strands', 'no_whitelist', 'exact'},
    ),

    _T(
        "igvf_crispr_preprocess",
        "Runs the IGVF CRISPR pipeline's RNA AnnData QC: applies gene names, knee plot, filter_cells(min_genes=100), filter_genes(min_cells=3), flags MT-/Mt- and RPS/RPL genes and computes scanpy calculate_qc_metrics (log1p). Writes filtered_anndata.h5ad, cell_qc.tsv, knee/violin/scatter figures, report.md and summary.json.",
        {"type": "object", "properties": {
            "adata_rna": {**_S_STRING, "description": "RNA counts h5ad"},
            "gene_names": {**_S_STRING, "description": "cells_x_genes.genes.names.txt (one name per gene)"},
            "min_genes": {**_S_INTEGER, "default": 100, "description": "min genes per cell"},
            "min_cells": {**_S_INTEGER, "default": 3, "description": "min cells per gene"},
            "reference": {**_S_STRING, "default": "human", "description": "'human' -> MT- prefix, otherwise Mt-"},
            "mt_prefix": {**_S_STRING, "description": "override the mitochondrial prefix"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "preprocess", "description": "run label"}},
         "required": ["adata_rna"]},
        cli=["igvf-crispr-pipeline", "preprocess"],
        flag_map={"adata_rna": "--adata-rna", "gene_names": "--gene-names", "min_genes": "--min-genes", "min_cells": "--min-cells", "reference": "--reference", "mt_prefix": "--mt-prefix", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "igvf_crispr_create_mdata",
        "Assembles the IGVF CRISPR pipeline MuData: names guide features sgRNA_ID|sgRNA_sequences, adds number_of_nonzero_guides and batch_number, draws the guide knee plot, intersects RNA and guide barcodes and writes mudata.h5mu with modalities 'transcripts' and 'guides', plus guide_cell_summary.tsv, report.md and summary.json.",
        {"type": "object", "properties": {
            "adata_rna": {**_S_STRING, "description": "filtered_anndata.h5ad"},
            "adata_guide": {**_S_STRING, "description": "guide counts h5ad"},
            "guide_metadata": {**_S_STRING, "description": "guide table (xlsx/tsv) with sgRNA_ID and sgRNA_sequences"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "positional guide naming and unordered barcode intersection as upstream"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "mudata", "description": "run label"}},
         "required": ["adata_rna", "adata_guide", "guide_metadata"]},
        cli=["igvf-crispr-pipeline", "create-mdata"],
        flag_map={"adata_rna": "--adata-rna", "adata_guide": "--adata-guide", "guide_metadata": "--guide-metadata", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'upstream_compat', 'no_plots'},
    ),

    _T(
        "igvf_crispr_run",
        "Runs the IGVF CRISPR pipeline chain on whatever inputs are given: seqspec parsing (RNA and guide), guide_features.txt, guide counting (kb when installed, else the pure-Python counter), RNA QC and MuData assembly. Writes every intermediate plus report.md and summary.json in one run directory.",
        {"type": "object", "properties": {
            "rna_seqspec": {**_S_STRING, "description": "RNA seqspec YAML"},
            "guide_seqspec": {**_S_STRING, "description": "guide seqspec YAML"},
            "guide_table": {**_S_STRING, "description": "guide table"},
            "adata_guide": {**_S_STRING, "description": "existing guide counts h5ad (skips counting)"},
            "guide_fastqs": {**_S_ARRAY_S, "description": "guide FASTQs R1 R2 ..."},
            "technology": {**_S_STRING, "description": "guide kb -x string when no guide seqspec"},
            "whitelist": {**_S_STRING, "description": "barcode onlist"},
            "no_whitelist": {**_S_BOOLEAN, "default": False, "description": "keep barcodes uncorrected"},
            "python_counts": {**_S_BOOLEAN, "default": False, "description": "use the Python counter even if kb is installed"},
            "adata_rna": {**_S_STRING, "description": "RNA counts h5ad"},
            "gene_names": {**_S_STRING, "description": "gene names file"},
            "min_genes": {**_S_INTEGER, "default": 100, "description": "min genes per cell"},
            "min_cells": {**_S_INTEGER, "default": 3, "description": "min cells per gene"},
            "reference": {**_S_STRING, "default": "human", "description": "species"},
            "mt_prefix": {**_S_STRING, "description": "override mitochondrial prefix"},
            "engine": {**_S_STRING, "default": "auto", "description": "auto | python | seqspec"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "upstream naming/ordering quirks"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "igvf_crispr_run", "description": "run label"}}},
        cli=["igvf-crispr-pipeline", "run"],
        flag_map={"rna_seqspec": "--rna-seqspec", "guide_seqspec": "--guide-seqspec", "guide_table": "--guide-table", "adata_guide": "--adata-guide", "guide_fastqs": "--guide-fastqs", "technology": "--technology", "whitelist": "--whitelist", "no_whitelist": "--no-whitelist", "python_counts": "--python-counts", "adata_rna": "--adata-rna", "gene_names": "--gene-names", "min_genes": "--min-genes", "min_cells": "--min-cells", "reference": "--reference", "mt_prefix": "--mt-prefix", "engine": "--engine", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'guide_fastqs'},
        bool_flags={'no_whitelist', 'python_counts', 'upstream_compat', 'no_plots'},
    ),

    _T(
        "sceptre_igvf_convert",
        "Converts an IGVF CRISPR MuData into a SCEPTRE object (response and gRNA matrices, gRNA-target table, cell covariates incl. response_n_umis/n_nonzero/p_mito and gRNA totals) with the sceptreIGVF multicollinearity check. Input: .h5mu. Output: run directory with matrices, covariate_data_frame.tsv, the auto-constructed formula, report.md and summary.json.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu) with mod/gene, mod/guide (intended_target_name, targeting, uns moi) and uns pairs_to_test"},
            "remove_collinear_covariates": {**_S_BOOLEAN, "default": False, "description": "drop extra covariates if their design matrix is rank deficient"},
            "collinear_mode": {**_S_STRING, "default": "all-or-nothing", "description": "all-or-nothing (sceptreIGVF) or greedy (IGVF CRISPR_Pipeline)"},
            "label": {**_S_STRING, "default": "convert", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["sceptre-igvf", "convert"],
        flag_map={"mudata": "--mudata", "label": "--label", "remove_collinear_covariates": "--remove-collinear-covariates", "collinear_mode": "--collinear-mode"},
        bool_flags={'remove_collinear_covariates'},
    ),

    _T(
        "sceptre_igvf_assign_guides",
        "Assigns gRNAs to cells as sceptreIGVF assign_grnas_sceptre does: SCEPTRE's Poisson-mixture EM per gRNA (or thresholding / maximum). Input: .h5mu with guide counts. Output: a MuData with the guide_assignment layer in mod/guide, guide_assignment.mtx (guides x cells), a per-gRNA table, report and summary.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu) with mod/gene, mod/guide (intended_target_name, targeting, uns moi) and uns pairs_to_test"},
            "method": {**_S_STRING, "default": "mixture", "description": "mixture (default), thresholding or maximum"},
            "threshold": {**_S_NUMBER, "description": "UMI cutoff for thresholding (default 5)"},
            "probability_threshold": {**_S_NUMBER, "default": 0.8, "description": "mixture posterior cutoff"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "assign_guides", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["sceptre-igvf", "assign-guides"],
        flag_map={"mudata": "--mudata", "label": "--label", "method": "--method", "threshold": "--threshold", "probability_threshold": "--probability-threshold", "no_plots": "--no-plots"},
        bool_flags={'no_plots'},
    ),

    _T(
        "sceptre_igvf_qc",
        "Runs SCEPTRE cell-wise QC (1%-99% response UMI / n_nonzero range, mitochondrial fraction, low-MOI zero-or-multiple-gRNA removal) after gRNA assignment and reports how many cells each filter removes. Input: .h5mu. Output: cells_in_use.tsv, report and summary.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu) with mod/gene, mod/guide (intended_target_name, targeting, uns moi) and uns pairs_to_test"},
            "method": {**_S_STRING, "default": "default", "description": "assignment method (default: sceptre default for the MOI)"},
            "p_mito_threshold": {**_S_NUMBER, "default": 0.2, "description": "max mitochondrial UMI fraction"},
            "label": {**_S_STRING, "default": "qc", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["sceptre-igvf", "qc"],
        flag_map={"mudata": "--mudata", "label": "--label", "method": "--method", "p_mito_threshold": "--p-mito-threshold"},
    ),

    _T(
        "sceptre_igvf_inference",
        "Tests the pairs_to_test of an IGVF MuData with SCEPTRE (NB GLM score statistic with conditional-resampling or permutation null and skew-normal tail), as sceptreIGVF inference_sceptre: guide assignment by 1-UMI threshold, relaxed QC. Writes uns test_results (pairs_to_test columns + p_value, log2_fc), or in pipeline mode per_element_results / per_guide_results (optionally cis_/trans_ prefixed). Output: annotated .h5mu, TSVs, QQ/volcano figures, report and summary.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu) with mod/gene, mod/guide (intended_target_name, targeting, uns moi) and uns pairs_to_test"},
            "side": {**_S_STRING, "default": "both", "description": "test sidedness: both, left (knockdown) or right"},
            "control_group": {**_S_STRING, "default": "default", "description": "default, complement (all other cells) or nt_cells (non-targeting cells)"},
            "resampling_mechanism": {**_S_STRING, "default": "default", "description": "default, crt (conditional resampling on covariates) or permutations"},
            "formula": {**_S_STRING, "description": "covariate formula, e.g. '~ log(response_n_umis) + batch'; 'default' = upstream default"},
            "grna_integration_strategy": {**_S_STRING, "default": "union", "description": "union (per element), singleton (per guide) or bonferroni"},
            "alpha": {**_S_NUMBER, "default": 0.1, "description": "BH FDR level for significance"},
            "mode": {**_S_STRING, "default": "sceptreIGVF", "description": "sceptreIGVF (test_results) or pipeline (per-element + per-guide results)"},
            "scope": {**_S_STRING, "description": "pipeline mode: cis or trans prefix for the uns keys"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "inference", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["sceptre-igvf", "inference"],
        flag_map={"mudata": "--mudata", "label": "--label", "side": "--side", "control_group": "--control-group", "resampling_mechanism": "--resampling-mechanism", "formula": "--formula", "grna_integration_strategy": "--grna-integration-strategy", "alpha": "--alpha", "mode": "--mode", "scope": "--scope", "no_plots": "--no-plots"},
        bool_flags={'no_plots'},
    ),

    _T(
        "sceptre_igvf_calibration_check",
        "Runs SCEPTRE's calibration check: undercover non-targeting gRNA groups paired with random genes, tested exactly like discovery pairs; a calibrated analysis makes ~0 discoveries. Input: .h5mu with non-targeting guides. Output: calibration_check_results.tsv, number of false discoveries, KS uniformity p, QQ plot, report and summary.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu) with mod/gene, mod/guide (intended_target_name, targeting, uns moi) and uns pairs_to_test"},
            "side": {**_S_STRING, "default": "both", "description": "test sidedness: both, left (knockdown) or right"},
            "control_group": {**_S_STRING, "default": "default", "description": "default, complement (all other cells) or nt_cells (non-targeting cells)"},
            "resampling_mechanism": {**_S_STRING, "default": "default", "description": "default, crt (conditional resampling on covariates) or permutations"},
            "formula": {**_S_STRING, "description": "covariate formula, e.g. '~ log(response_n_umis) + batch'; 'default' = upstream default"},
            "grna_integration_strategy": {**_S_STRING, "default": "union", "description": "union (per element), singleton (per guide) or bonferroni"},
            "alpha": {**_S_NUMBER, "default": 0.1, "description": "BH FDR level for significance"},
            "method": {**_S_STRING, "default": "default", "description": "gRNA assignment: default (mixture high MOI / maximum low MOI), mixture, thresholding, maximum"},
            "positive_control_pairs": {**_S_STRING, "description": "optional TSV of positive-control pairs (grna_target/response_id or intended_target_name/gene_id)"},
            "discovery_pairs": {**_S_STRING, "description": "optional TSV of discovery pairs (default: uns pairs_to_test)"},
            "pair_set": {**_S_STRING, "description": "construct discovery pairs instead: cis (targets within 500 kb of the gene TSS) or trans (all targets x all genes)"},
            "n_calibration_pairs": {**_S_INTEGER, "description": "number of NT pairs (default: number of discovery pairs)"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "calibration_check", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["sceptre-igvf", "calibration-check"],
        flag_map={"mudata": "--mudata", "label": "--label", "side": "--side", "control_group": "--control-group", "resampling_mechanism": "--resampling-mechanism", "formula": "--formula", "grna_integration_strategy": "--grna-integration-strategy", "alpha": "--alpha", "method": "--method", "positive_control_pairs": "--positive-control-pairs", "discovery_pairs": "--discovery-pairs", "pair_set": "--pair-set", "n_calibration_pairs": "--n-calibration-pairs", "no_plots": "--no-plots"},
        bool_flags={'no_plots'},
    ),

    _T(
        "sceptre_igvf_power_check",
        "Runs SCEPTRE's power check on positive-control pairs (from pairs_to_test pair_type == positive_control, uns positive_control_pairs, or a TSV). Input: .h5mu. Output: power_check_results.tsv with p-values and log2 fold changes, report and summary.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu) with mod/gene, mod/guide (intended_target_name, targeting, uns moi) and uns pairs_to_test"},
            "side": {**_S_STRING, "default": "both", "description": "test sidedness: both, left (knockdown) or right"},
            "control_group": {**_S_STRING, "default": "default", "description": "default, complement (all other cells) or nt_cells (non-targeting cells)"},
            "resampling_mechanism": {**_S_STRING, "default": "default", "description": "default, crt (conditional resampling on covariates) or permutations"},
            "formula": {**_S_STRING, "description": "covariate formula, e.g. '~ log(response_n_umis) + batch'; 'default' = upstream default"},
            "grna_integration_strategy": {**_S_STRING, "default": "union", "description": "union (per element), singleton (per guide) or bonferroni"},
            "alpha": {**_S_NUMBER, "default": 0.1, "description": "BH FDR level for significance"},
            "method": {**_S_STRING, "default": "default", "description": "gRNA assignment: default (mixture high MOI / maximum low MOI), mixture, thresholding, maximum"},
            "positive_control_pairs": {**_S_STRING, "description": "optional TSV of positive-control pairs (grna_target/response_id or intended_target_name/gene_id)"},
            "discovery_pairs": {**_S_STRING, "description": "optional TSV of discovery pairs (default: uns pairs_to_test)"},
            "pair_set": {**_S_STRING, "description": "construct discovery pairs instead: cis (targets within 500 kb of the gene TSS) or trans (all targets x all genes)"},
            "label": {**_S_STRING, "default": "power_check", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["sceptre-igvf", "power-check"],
        flag_map={"mudata": "--mudata", "label": "--label", "side": "--side", "control_group": "--control-group", "resampling_mechanism": "--resampling-mechanism", "formula": "--formula", "grna_integration_strategy": "--grna-integration-strategy", "alpha": "--alpha", "method": "--method", "positive_control_pairs": "--positive-control-pairs", "discovery_pairs": "--discovery-pairs", "pair_set": "--pair-set"},
    ),

    _T(
        "sceptre_igvf_run",
        "Runs the whole SCEPTRE workflow on an IGVF MuData: gRNA assignment, QC, calibration check, power check and discovery analysis, then writes a MuData in the sceptreIGVF sceptre_object_to_mudata schema (pairs_to_test with pair_type, test_results with p_value and log2_fc, guide_assignment layer). Output: .h5mu, result TSVs, figures, report and summary.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu) with mod/gene, mod/guide (intended_target_name, targeting, uns moi) and uns pairs_to_test"},
            "side": {**_S_STRING, "default": "both", "description": "test sidedness: both, left (knockdown) or right"},
            "control_group": {**_S_STRING, "default": "default", "description": "default, complement (all other cells) or nt_cells (non-targeting cells)"},
            "resampling_mechanism": {**_S_STRING, "default": "default", "description": "default, crt (conditional resampling on covariates) or permutations"},
            "formula": {**_S_STRING, "description": "covariate formula, e.g. '~ log(response_n_umis) + batch'; 'default' = upstream default"},
            "grna_integration_strategy": {**_S_STRING, "default": "union", "description": "union (per element), singleton (per guide) or bonferroni"},
            "alpha": {**_S_NUMBER, "default": 0.1, "description": "BH FDR level for significance"},
            "method": {**_S_STRING, "default": "default", "description": "gRNA assignment: default (mixture high MOI / maximum low MOI), mixture, thresholding, maximum"},
            "positive_control_pairs": {**_S_STRING, "description": "optional TSV of positive-control pairs (grna_target/response_id or intended_target_name/gene_id)"},
            "discovery_pairs": {**_S_STRING, "description": "optional TSV of discovery pairs (default: uns pairs_to_test)"},
            "pair_set": {**_S_STRING, "description": "construct discovery pairs instead: cis (targets within 500 kb of the gene TSS) or trans (all targets x all genes)"},
            "skip_calibration": {**_S_BOOLEAN, "default": False, "description": "skip the NT calibration check"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "sceptre_run", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["sceptre-igvf", "run"],
        flag_map={"mudata": "--mudata", "label": "--label", "side": "--side", "control_group": "--control-group", "resampling_mechanism": "--resampling-mechanism", "formula": "--formula", "grna_integration_strategy": "--grna-integration-strategy", "alpha": "--alpha", "method": "--method", "positive_control_pairs": "--positive-control-pairs", "discovery_pairs": "--discovery-pairs", "pair_set": "--pair-set", "skip_calibration": "--skip-calibration", "no_plots": "--no-plots"},
        bool_flags={'no_plots', 'skip_calibration'},
    ),

    _T(
        "sceptre_igvf_export_mudata",
        "Builds the eight benchmark MuData files of sceptreIGVF sceptre_object_to_mudata_inputs_outputs (inference/guide-assignment inputs and outputs, full and minimal): positive controls with results plus a sample of significant and non-significant discovery pairs, re-analysed with default settings. Output: guide_assignment/ and inference/ folders of .h5mu files, report and summary.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu) with mod/gene, mod/guide (intended_target_name, targeting, uns moi) and uns pairs_to_test"},
            "num_discovery_pairs": {**_S_INTEGER, "default": 100, "description": "discovery pairs to keep (half significant when possible)"},
            "prefix": {**_S_STRING, "default": "", "description": "file-name prefix, e.g. gasperini_"},
            "out_dir": {**_S_STRING, "description": "output directory (default: run dir)"},
            "guide_capture_method": {**_S_STRING, "description": "stored in guide uns capture_method"},
            "label": {**_S_STRING, "default": "export_mudata", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["sceptre-igvf", "export-mudata"],
        flag_map={"mudata": "--mudata", "label": "--label", "num_discovery_pairs": "--num-discovery-pairs", "prefix": "--prefix", "out_dir": "--out-dir", "guide_capture_method": "--guide-capture-method"},
    ),

    _T(
        "sceptre_igvf_r_runner",
        "Runs the real sceptreIGVF R package (assign_grnas_sceptre or inference_sceptre) through Rscript when R, MuData and sceptreIGVF are installed; otherwise writes the R script and prints the exact command plus install instructions. Output: the R script, R output .h5mu and test_results TSV when run, report and summary.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu) with mod/gene, mod/guide (intended_target_name, targeting, uns moi) and uns pairs_to_test"},
            "step": {**_S_STRING, "description": "assign or inference"},
            "side": {**_S_STRING, "description": "inference side (both/left/right)"},
            "formula": {**_S_STRING, "description": "R formula for inference, e.g. '~ log(response_n_umis)'"},
            "label": {**_S_STRING, "default": "r_runner", "description": "run label for the output directory"}},
         "required": ["mudata", "step"]},
        cli=["sceptre-igvf", "r-runner"],
        flag_map={"mudata": "--mudata", "label": "--label", "step": "--step", "side": "--side", "formula": "--formula"},
    ),

    _T(
        "sceptre_igvf_compare_results",
        "Compares two SCEPTRE result tables (e.g. this Python port vs the R package) on shared pair keys: Spearman correlation of -log10 p, log2 fold-change correlation and BH significance agreement. Input: two TSV/CSV tables. Output: merged table, report and summary.",
        {"type": "object", "properties": {
            "a": {**_S_STRING, "description": "first results table"},
            "b": {**_S_STRING, "description": "second results table"},
            "keys": {**_S_ARRAY_S, "default": ["intended_target_name", "gene_id"], "description": "join columns"},
            "label": {**_S_STRING, "default": "compare", "description": "run label for the output directory"}},
         "required": ["a", "b"]},
        cli=["sceptre-igvf", "compare-results"],
        flag_map={"a": "--a", "b": "--b", "keys": "--keys", "label": "--label"},
        flag_repeat={'keys'},
    ),

    _T(
        "flowfish_pipeline_run",
        "Runs the whole CRISPRi-FlowFISH Snakemake workflow (EngreitzLab/crispri-flowfish) from a sample sheet: sample-sheet validation and experiment keys, read mapping (bowtie or an exact Python mapper) and guide counts, count tables per PCR and experimental replicate, replicate correlations, per-guide MLE effect sizes from sorted bins, real-space normalisation to negative controls, 10-guide windows, TSS qPCR scaling, collapse to candidate elements and Mann-Whitney / t-test peak calling. Inputs are the sample sheet, guide design, sort-params directory ({Batch}_{Sample}.txt), FASTQ directory (or pre-computed count files), and optionally qPCR, gene list and enhancer BED. Outputs a results/ tree with upstream file names (raw_effects, real_space, windows, scaled, FullEnhancerScore, PeakCallingSummary, KnownEnhancers) plus report.md and summary.json.",
        {"type": "object", "properties": {
            "sample_sheet": {**_S_STRING, "description": "tab-delimited sample sheet: SampleID Bin PCRRep GeneSymbol qPCRGene + key columns (+ Batch, Sample for the sort-params file)"},
            "design": {**_S_STRING, "description": "guide design file (OligoID, MappingSequence, GuideSequence, OffTargetScore, target with 'negative_control', chr/start/end)"},
            "sortparams_dir": {**_S_STRING, "description": "directory of sort-params files named {Batch}_{Sample}.txt"},
            "experiment_keycols": {**_S_STRING, "default": "", "description": "comma-separated experiment key columns"},
            "replicate_keycols": {**_S_STRING, "default": "", "description": "comma-separated replicate key columns"},
            "fastq_dir": {**_S_STRING, "description": "directory with {SampleID}_*_R1_*.fastq.gz"},
            "counts_dir": {**_S_STRING, "description": "pre-computed {SampleID}.count.txt files (skips mapping)"},
            "bowtie_index": {**_S_STRING, "description": "bowtie index prefix (used when bowtie is on PATH)"},
            "python_mapper": {**_S_BOOLEAN, "default": False, "description": "map with the exact-match Python mapper even if bowtie is available"},
            "qpcr": {**_S_STRING, "description": "TSS qPCR table (qPCRGene, name, TSS_qPCR[, TSS_Override])"},
            "genelist": {**_S_STRING, "description": "gene list with name, chr, tss"},
            "enhancers": {**_S_STRING, "description": "candidate element BED (first 3 columns used)"},
            "cell_line": {**_S_STRING, "default": "", "description": "cell type label for KnownEnhancers"},
            "sort_params_format": {**_S_STRING, "default": "bigfoot_nototal", "description": "bigfoot_nototal (upstream default) | bigfoot | astrios | astrios_nototal"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "reproduce R 3.6 factor-code Count indexing and positional bin pairing"},
            "window": {**_S_INTEGER, "default": 10, "description": "guides per sliding window"},
            "max_span": {**_S_NUMBER, "default": 750, "description": "max window span (bp)"},
            "minsum": {**_S_NUMBER, "default": 0.05, "description": "min sorted cells per guide"},
            "minbins": {**_S_INTEGER, "default": 4, "description": "min bins a guide is observed in"},
            "min_guides": {**_S_INTEGER, "default": 3, "description": "min guides per element"},
            "fdr": {**_S_NUMBER, "default": 0.05, "description": "FDR threshold"},
            "power": {**_S_BOOLEAN, "default": False, "description": "also run the power analysis per replicate"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "flowfish", "description": "run label (output folder suffix)"}},
         "required": ["sample_sheet", "design", "sortparams_dir"]},
        cli=["flowfish-pipeline", "run"],
        flag_map={"sample_sheet": "--sample-sheet", "design": "--design", "sortparams_dir": "--sortparams-dir", "experiment_keycols": "--experiment-keycols", "replicate_keycols": "--replicate-keycols", "fastq_dir": "--fastq-dir", "counts_dir": "--counts-dir", "bowtie_index": "--bowtie-index", "python_mapper": "--python-mapper", "qpcr": "--qpcr", "genelist": "--genelist", "enhancers": "--enhancers", "cell_line": "--cell-line", "sort_params_format": "--sort-params-format", "upstream_compat": "--upstream-compat", "window": "--window", "max_span": "--max-span", "minsum": "--minsum", "minbins": "--minbins", "min_guides": "--min-guides", "fdr": "--fdr", "power": "--power", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots', 'upstream_compat', 'python_mapper', 'power'},
    ),

    _T(
        "flowfish_pipeline_estimate_effects",
        "Computes per-guide expression estimates from a FlowFISH bin_counts table exactly as estimate_effect_sizes.R: reads rescaled to sorted-cell counts, bin-midpoint weighted average, and a log-normal MLE with an EM-estimated seventh (ungated) bin, falling back to the weighted average when the fit fails. Inputs are a bin_counts.txt and a sort-params file. Outputs <name>.raw_effects.txt (WeightedAvg, logMean, logSD, method) with report.md and summary.json.",
        {"type": "object", "properties": {
            "counts": {**_S_STRING, "description": "bin_counts.txt (OligoID + one column per bin, optional All)"},
            "sort_params": {**_S_STRING, "description": "sort-params file"},
            "sort_params_format": {**_S_STRING, "default": "bigfoot_nototal", "description": "bigfoot_nototal | bigfoot | astrios | astrios_nototal"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "R 3.6 factor-code Count indexing"},
            "mu_seed": {**_S_NUMBER, "description": "override the log10 mean seed"},
            "sd_seed": {**_S_NUMBER, "description": "override the log10 sd seed"},
            "name": {**_S_STRING, "description": "output file prefix"},
            "label": {**_S_STRING, "default": "estimate_effects", "description": "run label (output folder suffix)"}},
         "required": ["counts", "sort_params"]},
        cli=["flowfish-pipeline", "estimate-effects"],
        flag_map={"counts": "--counts", "sort_params": "--sort-params", "sort_params_format": "--sort-params-format", "upstream_compat": "--upstream-compat", "mu_seed": "--mu-seed", "sd_seed": "--sd-seed", "name": "--name", "label": "--label"},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "flowfish_pipeline_real_space",
        "Converts MLE log10 guide estimates to real space and normalises them to negative controls as convert_to_real_space.py does: guide filters (length 18-21, off-target >= 50, <= 10 G, min cells, min bins), division by the control median, clamp to [0, 5], and rescale to control mean 1. Inputs are raw_effects.txt and the guide design. Outputs real_space.txt, bedgraphs and a convert log.",
        {"type": "object", "properties": {
            "mle": {**_S_STRING, "description": "raw_effects.txt"},
            "design": {**_S_STRING, "description": "guide design file"},
            "minsum": {**_S_NUMBER, "default": 0.0, "description": "min sumcells (the Snakefile uses 0.05)"},
            "minbins": {**_S_INTEGER, "default": 0, "description": "min observed bins (the Snakefile uses 4)"},
            "clamp": {**_S_NUMBER, "default": 5.0, "description": "max effect"},
            "minofftarget": {**_S_NUMBER, "default": 50, "description": "min OffTargetScore"},
            "name": {**_S_STRING, "description": "output prefix"},
            "label": {**_S_STRING, "default": "real_space", "description": "run label (output folder suffix)"}},
         "required": ["mle", "design"]},
        cli=["flowfish-pipeline", "real-space"],
        flag_map={"mle": "--mle", "design": "--design", "minsum": "--minsum", "minbins": "--minbins", "clamp": "--clamp", "minofftarget": "--minofftarget", "name": "--name", "label": "--label"},
    ),

    _T(
        "flowfish_pipeline_windows",
        "Computes sliding-window guide statistics (CalculateTilingStatistic.R): windows of consecutive non-control guides spanning at most max_span bp, mean effect relative to negative controls, and Mann-Whitney / Welch t tests with BH FDR. Input is real_space.txt. Outputs windows.txt and a mean bedgraph.",
        {"type": "object", "properties": {
            "input": {**_S_STRING, "description": "real_space.txt"},
            "window": {**_S_INTEGER, "default": 10, "description": "guides per window"},
            "max_span": {**_S_NUMBER, "default": 750, "description": "max span bp"},
            "score_column": {**_S_STRING, "default": "mleAvg", "description": "score column"},
            "name": {**_S_STRING, "description": "output prefix"},
            "label": {**_S_STRING, "default": "windows", "description": "run label (output folder suffix)"}},
         "required": ["input"]},
        cli=["flowfish-pipeline", "windows"],
        flag_map={"input": "--input", "window": "--window", "max_span": "--max-span", "score_column": "--score-column", "name": "--name", "label": "--label"},
    ),

    _T(
        "flowfish_pipeline_tss_kd",
        "Finds the strongest guide window within +/-500 bp of the target gene TSS and computes FlowFISH_at_TSS and the background term (FF - qPCR)/(1 - FF), as FlowFISHtssKD.py does. Inputs are windows.txt, the gene list and an optional qPCR table. Outputs ScreenInfo.txt.",
        {"type": "object", "properties": {
            "gene": {**_S_STRING, "description": "qPCRGene (or Screen id containing '-')"},
            "windows": {**_S_STRING, "description": "windows.txt"},
            "genelist": {**_S_STRING, "description": "gene list (name, chr, tss)"},
            "qpcr": {**_S_STRING, "description": "qPCR table; omit to skip the adjustment"},
            "slop": {**_S_INTEGER, "default": 500, "description": "bp around the TSS"},
            "name": {**_S_STRING, "description": "output prefix"},
            "label": {**_S_STRING, "default": "tss_kd", "description": "run label (output folder suffix)"}},
         "required": ["gene", "windows", "genelist"]},
        cli=["flowfish-pipeline", "tss-kd"],
        flag_map={"gene": "--gene", "windows": "--windows", "genelist": "--genelist", "qpcr": "--qpcr", "slop": "--slop", "name": "--name", "label": "--label"},
    ),

    _T(
        "flowfish_pipeline_normalize_qpcr",
        "Linearly rescales guide effects so the TSS knockdown matches qPCR (normalize_flowfish_to_qpcr.py): slope = (1 - qPCR)/(1 - FF), mleAvg -> slope*(x - FF) + qPCR clipped to [0, 5], mleSD scaled by the slope. Inputs are ScreenInfo.txt and real_space.txt. Outputs scaled.txt and scaled.bedgraph.",
        {"type": "object", "properties": {
            "info": {**_S_STRING, "description": "ScreenInfo.txt"},
            "input": {**_S_STRING, "description": "real_space.txt"},
            "clamp": {**_S_NUMBER, "default": 5.0, "description": "max effect"},
            "name": {**_S_STRING, "description": "output prefix"},
            "label": {**_S_STRING, "default": "normalize_qpcr", "description": "run label (output folder suffix)"}},
         "required": ["info", "input"]},
        cli=["flowfish-pipeline", "normalize-qpcr"],
        flag_map={"info": "--info", "input": "--input", "clamp": "--clamp", "name": "--name", "label": "--label"},
    ),

    _T(
        "flowfish_pipeline_collapse",
        "Overlaps non-control guides with a candidate-element BED and collapses their scaled effects per element (the upstream bedtools map step, done in pandas). Inputs are scaled.txt and the element BED. Outputs collapse.bed.",
        {"type": "object", "properties": {
            "scaled": {**_S_STRING, "description": "scaled.txt"},
            "enhancers": {**_S_STRING, "description": "element BED"},
            "name": {**_S_STRING, "description": "output prefix"},
            "label": {**_S_STRING, "default": "collapse", "description": "run label (output folder suffix)"}},
         "required": ["scaled", "enhancers"]},
        cli=["flowfish-pipeline", "collapse"],
        flag_map={"scaled": "--scaled", "enhancers": "--enhancers", "name": "--name", "label": "--label"},
    ),

    _T(
        "flowfish_pipeline_score_enhancers",
        "Scores candidate elements against negative-control guides as ScoreEnhancers.py does: elements with >= min_guides guides, mean effect, Mann-Whitney (scipy 1.5 default one-sided) and Student t-test with BH FDR, Significant (t-test FDR) and Regulated calls. Inputs are collapse.bed and scaled.txt. Outputs FullEnhancerScore.txt, the significant-peaks BED and PeakCallingSummary.txt.",
        {"type": "object", "properties": {
            "collapsed": {**_S_STRING, "description": "collapse.bed"},
            "scaled": {**_S_STRING, "description": "scaled.txt"},
            "expt_name": {**_S_STRING, "default": "", "description": "experiment name for the summary"},
            "min_guides": {**_S_INTEGER, "default": 3, "description": "min guides per element"},
            "fdr": {**_S_NUMBER, "default": 0.05, "description": "FDR threshold"},
            "min_effect": {**_S_NUMBER, "default": 0.0, "description": "min |mean - control mean|"},
            "utest_mode": {**_S_STRING, "default": "legacy", "description": "legacy | two-sided"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "score_enhancers", "description": "run label (output folder suffix)"}},
         "required": ["collapsed", "scaled"]},
        cli=["flowfish-pipeline", "score-enhancers"],
        flag_map={"collapsed": "--collapsed", "scaled": "--scaled", "expt_name": "--expt-name", "min_guides": "--min-guides", "fdr": "--fdr", "min_effect": "--min-effect", "utest_mode": "--utest-mode", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "flowfish_pipeline_power",
        "Estimates the power of the element t-test (the ScoreEnhancers.py power block, with its indentation bug fixed): negative-control scores shifted by background-adjusted effect sizes, sampled at the median guides per element, and BH-tested alongside the real elements. Inputs are collapse.bed, scaled.txt and optionally ScreenInfo.txt. Outputs a power table (EffectSize, Power).",
        {"type": "object", "properties": {
            "collapsed": {**_S_STRING, "description": "collapse.bed"},
            "scaled": {**_S_STRING, "description": "scaled.txt"},
            "info": {**_S_STRING, "description": "ScreenInfo.txt (Background)"},
            "reps": {**_S_INTEGER, "default": 400, "description": "simulations per effect size"},
            "effect_sizes": {"type": "array", "items": {"type": "number"}, "default": [0.05, 0.1, 0.15, 0.2, 0.25, 0.3], "description": "effect sizes to test"},
            "expt_name": {**_S_STRING, "default": "", "description": "name"},
            "fdr": {**_S_NUMBER, "default": 0.05, "description": "FDR"},
            "label": {**_S_STRING, "default": "power", "description": "run label (output folder suffix)"}},
         "required": ["collapsed", "scaled"]},
        cli=["flowfish-pipeline", "power"],
        flag_map={"collapsed": "--collapsed", "scaled": "--scaled", "info": "--info", "reps": "--reps", "effect_sizes": "--effect-sizes", "expt_name": "--expt-name", "fdr": "--fdr", "label": "--label"},
        flag_repeat={'effect_sizes'},
    ),

    _T(
        "flowfish_pipeline_format_screen",
        "Writes ScreenData.txt and a KnownEnhancers-style table from FullEnhancerScore.txt: gene TSS, distance, fraction change in expression, adjusted p, Significant/Regulated, cell type. Inputs are the score table, the gene symbol and the gene list. Outputs ScreenData.txt and KnownEnhancers.FlowFISH.txt.",
        {"type": "object", "properties": {
            "score": {**_S_STRING, "description": "FullEnhancerScore.txt"},
            "gene": {**_S_STRING, "description": "target gene symbol"},
            "genelist": {**_S_STRING, "description": "gene list"},
            "cell_line": {**_S_STRING, "default": "", "description": "cell type"},
            "name": {**_S_STRING, "description": "output prefix"},
            "label": {**_S_STRING, "default": "format_screen", "description": "run label (output folder suffix)"}},
         "required": ["score", "gene"]},
        cli=["flowfish-pipeline", "format-screen"],
        flag_map={"score": "--score", "gene": "--gene", "genelist": "--genelist", "cell_line": "--cell-line", "name": "--name", "label": "--label"},
    ),

    _T(
        "flowfish_pipeline_count_tables",
        "Builds guide count tables for a FlowFISH screen: per PCR replicate and per experimental replicate, counts summed per sorting bin with empty bins dropped, bin frequencies, the flat GuideCounts table, and PCR-replicate correlations. Inputs are the sample sheet, key columns and a directory of {SampleID}.count.txt. Outputs bin_counts / bin_freq tables and ReplicateCorrelation tsv.",
        {"type": "object", "properties": {
            "sample_sheet": {**_S_STRING, "description": "sample sheet"},
            "counts_dir": {**_S_STRING, "description": "directory of {SampleID}.count.txt"},
            "experiment_keycols": {**_S_STRING, "default": "", "description": "key columns"},
            "replicate_keycols": {**_S_STRING, "default": "", "description": "replicate columns"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "count_tables", "description": "run label (output folder suffix)"}},
         "required": ["sample_sheet", "counts_dir"]},
        cli=["flowfish-pipeline", "count-tables"],
        flag_map={"sample_sheet": "--sample-sheet", "counts_dir": "--counts-dir", "experiment_keycols": "--experiment-keycols", "replicate_keycols": "--replicate-keycols", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "flowfish_pipeline_map_reads",
        "Maps one FASTQ of guide amplicon reads with bowtie -v0 --all when bowtie and an index are available, or with an exact-match Python mapper given the design, then counts reads per guide. Inputs are a FASTQ, a sample id and a bowtie index or design. Outputs the alignment file, {SampleID}.count.txt and alignment_stats; with neither, prints the bowtie command.",
        {"type": "object", "properties": {
            "fastq": {**_S_STRING, "description": "FASTQ(.gz)"},
            "sample_id": {**_S_STRING, "description": "sample id"},
            "index": {**_S_STRING, "description": "bowtie index prefix"},
            "design": {**_S_STRING, "description": "guide design (Python mapper)"},
            "python_mapper": {**_S_BOOLEAN, "default": False, "description": "force the Python mapper"},
            "label": {**_S_STRING, "default": "map_reads", "description": "run label (output folder suffix)"}},
         "required": ["fastq", "sample_id"]},
        cli=["flowfish-pipeline", "map-reads"],
        flag_map={"fastq": "--fastq", "sample_id": "--sample-id", "index": "--index", "design": "--design", "python_mapper": "--python-mapper", "label": "--label"},
        bool_flags={'python_mapper'},
    ),

    _T(
        "flowfish_pipeline_samplesheet",
        "Validates a CRISPRi-FlowFISH sample sheet (required columns, unique SampleID, no 'Water' bins), discovers FASTQs, and adds the ExperimentIDPCRRep / ExperimentIDReplicates / ExperimentID keys. Inputs are the sheet, key columns and FASTQ directory. Outputs SampleList.snakemake.tsv.",
        {"type": "object", "properties": {
            "sample_sheet": {**_S_STRING, "description": "sample sheet"},
            "experiment_keycols": {**_S_STRING, "default": "", "description": "key columns"},
            "replicate_keycols": {**_S_STRING, "default": "", "description": "replicate columns"},
            "fastq_dir": {**_S_STRING, "description": "FASTQ directory"},
            "label": {**_S_STRING, "default": "samplesheet", "description": "run label (output folder suffix)"}},
         "required": ["sample_sheet"]},
        cli=["flowfish-pipeline", "samplesheet"],
        flag_map={"sample_sheet": "--sample-sheet", "experiment_keycols": "--experiment-keycols", "replicate_keycols": "--replicate-keycols", "fastq_dir": "--fastq-dir", "label": "--label"},
    ),

    _T(
        "flowfish_pipeline_guide_count_plots",
        "Runs guide-count QC (PlotGuideCounts.R adapted to guide counts): PCR-replicate correlations of %Reads per replicate and bin, per-guide coefficient of variation, and per-bin frequency relative to the average bin. Inputs are GuideCounts.flat.tsv.gz and SampleList.snakemake.tsv. Outputs three TSVs and figures.",
        {"type": "object", "properties": {
            "guide_counts": {**_S_STRING, "description": "GuideCounts.flat.tsv.gz"},
            "sample_sheet": {**_S_STRING, "description": "SampleList.snakemake.tsv"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "guide_counts", "description": "run label (output folder suffix)"}},
         "required": ["guide_counts", "sample_sheet"]},
        cli=["flowfish-pipeline", "guide-count-plots"],
        flag_map={"guide_counts": "--guide-counts", "sample_sheet": "--sample-sheet", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "perturb_tools_run",
        "Runs the perturb-tools bulk pooled CRISPR screen pipeline end to end (the upstream tutorials): log2(RPM+1) normalisation, per-replicate log fold changes of cond1 vs cond2 with mean/median/sd aggregates, optional sorting-screen delta-LFC (regex-selected condit_1 / condit_2 / control samples, t-test p-values), sample QC (count distribution, Gini, count correlation, replicate LFC correlation for all and positive-control guides, outlier jackpot guides) and a MAGeCK count table. Input is a guides x samples count table or a saved screen. Writes guides_lfc.tsv, sample_qc.tsv, correlation tables, mageck_count.txt, the screen (.h5ad + CSV dir), figures, report.md and summary.json.",
        {"type": "object", "properties": {
            "screen": {**_S_STRING, "description": ".h5ad or screen CSV directory written by perturb-tools (alternative to counts)"},
            "counts": {**_S_STRING, "description": "guides x samples count table (TSV/CSV; first column guide id; text columns become guide annotations, e.g. TKO readcount files)"},
            "guides": {**_S_STRING, "description": "optional guide annotation table (first column guide id)"},
            "samples": {**_S_STRING, "description": "optional sample annotation table (column 'sample' or first column)"},
            "parse_tko_names": {**_S_BOOLEAN, "default": False, "description": "derive replicate (last char) and time (int of name[1:-1]) from sample names like T18A"},
            "exclude": {**_S_ARRAY_S, "description": "drop samples by obs value, e.g. replicate=0"},
            "cond1": {**_S_STRING, "description": "condition 1 value of compare_col (numerator), e.g. 18"},
            "cond2": {**_S_STRING, "description": "condition 2 value (denominator), e.g. 8"},
            "rep_col": {**_S_STRING, "default": "replicate"},
            "compare_col": {**_S_STRING, "default": "time"},
            "condit_1": {**_S_STRING, "description": "sorting screen: regex for condition-1 sample columns (e.g. high)"},
            "condit_2": {**_S_STRING, "description": "regex for condition-2 columns (e.g. low)"},
            "control": {**_S_STRING, "description": "regex for control columns (e.g. presort)"},
            "targets": {**_S_ARRAY_S, "description": "sorting screen: keep guides whose target is listed"},
            "target_col": {**_S_STRING, "description": "guide feature holding the target gene (e.g. GENE); enables MAGeCK export and control groups"},
            "pos_ctrl": {**_S_STRING, "description": "positive-control genes: file (first column) or comma list"},
            "neg_ctrl": {**_S_STRING, "description": "negative-control genes: file or comma list"},
            "method": {**_S_STRING, "default": "spearman", "description": "pearson or spearman"},
            "excel": {**_S_BOOLEAN, "default": False},
            "no_plots": {**_S_BOOLEAN, "default": False},
            "label": {**_S_STRING}}},
        cli=["perturb-tools", "run"],
        flag_map={"screen": "--screen", "counts": "--counts", "guides": "--guides", "samples": "--samples", "parse_tko_names": "--parse-tko-names", "exclude": "--exclude", "cond1": "--cond1", "cond2": "--cond2", "rep_col": "--rep-col", "compare_col": "--compare-col", "condit_1": "--condit-1", "condit_2": "--condit-2", "control": "--control", "targets": "--targets", "target_col": "--target-col", "pos_ctrl": "--pos-ctrl", "neg_ctrl": "--neg-ctrl", "method": "--method", "excel": "--excel", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'exclude', 'targets'},
        bool_flags={'excel', 'no_plots', 'parse_tko_names'},
    ),

    _T(
        "perturb_tools_make_screen",
        "Builds a perturb-tools screen object (AnnData layout: samples x guides with sample and guide annotation) from a guides x samples count table plus optional guide/sample tables, optionally parsing TKO-style sample names into replicate and time. Writes screen.h5ad (when anndata is installed), a screen CSV directory reusable as screen= by the other perturb-tools tools, guides.tsv and samples.tsv.",
        {"type": "object", "properties": {
            "screen": {**_S_STRING, "description": ".h5ad or screen CSV directory written by perturb-tools (alternative to counts)"},
            "counts": {**_S_STRING, "description": "guides x samples count table (TSV/CSV; first column guide id; text columns become guide annotations, e.g. TKO readcount files)"},
            "guides": {**_S_STRING, "description": "optional guide annotation table (first column guide id)"},
            "samples": {**_S_STRING, "description": "optional sample annotation table (column 'sample' or first column)"},
            "parse_tko_names": {**_S_BOOLEAN, "default": False, "description": "derive replicate (last char) and time (int of name[1:-1]) from sample names like T18A"},
            "exclude": {**_S_ARRAY_S, "description": "drop samples by obs value, e.g. replicate=0"},
            "log_norm": {**_S_BOOLEAN, "default": False, "description": "also add the lognorm_counts layer"},
            "label": {**_S_STRING}}},
        cli=["perturb-tools", "make-screen"],
        flag_map={"screen": "--screen", "counts": "--counts", "guides": "--guides", "samples": "--samples", "parse_tko_names": "--parse-tko-names", "exclude": "--exclude", "log_norm": "--log-norm", "label": "--label"},
        flag_repeat={'exclude'},
        bool_flags={'log_norm', 'parse_tko_names'},
    ),

    _T(
        "perturb_tools_lfc_reps",
        "Computes per-replicate guide log fold changes (log2(RPM+1) of the cond1 sample minus the cond2 sample within each replicate, columns '{rep}.{cond1}_{cond2}.lfc') and their mean / median / sd aggregates ('{cond1}_{cond2}.lfc.{fn}'), plus the replicate LFC correlation and, with target_col and control gene lists, the median LFC per guide group. Writes guides_lfc.tsv, lfc_replicate_correlation.tsv, the updated screen, figures and a report.",
        {"type": "object", "properties": {
            "screen": {**_S_STRING, "description": ".h5ad or screen CSV directory written by perturb-tools (alternative to counts)"},
            "counts": {**_S_STRING, "description": "guides x samples count table (TSV/CSV; first column guide id; text columns become guide annotations, e.g. TKO readcount files)"},
            "guides": {**_S_STRING, "description": "optional guide annotation table (first column guide id)"},
            "samples": {**_S_STRING, "description": "optional sample annotation table (column 'sample' or first column)"},
            "parse_tko_names": {**_S_BOOLEAN, "default": False, "description": "derive replicate (last char) and time (int of name[1:-1]) from sample names like T18A"},
            "exclude": {**_S_ARRAY_S, "description": "drop samples by obs value, e.g. replicate=0"},
            "cond1": {**_S_STRING},
            "cond2": {**_S_STRING},
            "rep_col": {**_S_STRING, "default": "replicate"},
            "compare_col": {**_S_STRING, "default": "sort"},
            "aggregate": {**_S_ARRAY_S, "description": "mean, median and/or sd (default median)"},
            "target_col": {**_S_STRING},
            "pos_ctrl": {**_S_STRING},
            "neg_ctrl": {**_S_STRING},
            "method": {**_S_STRING, "default": "spearman"},
            "no_plots": {**_S_BOOLEAN, "default": False},
            "label": {**_S_STRING}},
         "required": ["cond1", "cond2"]},
        cli=["perturb-tools", "lfc-reps"],
        flag_map={"screen": "--screen", "counts": "--counts", "guides": "--guides", "samples": "--samples", "parse_tko_names": "--parse-tko-names", "exclude": "--exclude", "cond1": "--cond1", "cond2": "--cond2", "rep_col": "--rep-col", "compare_col": "--compare-col", "aggregate": "--aggregate", "target_col": "--target-col", "pos_ctrl": "--pos-ctrl", "neg_ctrl": "--neg-ctrl", "method": "--method", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'exclude', 'aggregate'},
        bool_flags={'no_plots', 'parse_tko_names'},
    ),

    _T(
        "perturb_tools_sort_lfc",
        "Computes the sorting-screen (e.g. FACS high vs low bin) guide statistic of perturb-tools: log-normalised counts minus the control (presort) per replicate, lfc.mean = mean(condit_1 - control) - mean(condit_2 - control), lfc.stdev, a two-sided t-test p-value per guide and -log10(pval). Sample columns are selected by regex; guides can be restricted to targets. Writes guide_lfc.tsv, a volcano plot and, when guides carry Start/End/center positions, the guide-enrichment-along-locus plot (optional GTF exon track).",
        {"type": "object", "properties": {
            "screen": {**_S_STRING, "description": ".h5ad or screen CSV directory written by perturb-tools (alternative to counts)"},
            "counts": {**_S_STRING, "description": "guides x samples count table (TSV/CSV; first column guide id; text columns become guide annotations, e.g. TKO readcount files)"},
            "guides": {**_S_STRING, "description": "optional guide annotation table (first column guide id)"},
            "samples": {**_S_STRING, "description": "optional sample annotation table (column 'sample' or first column)"},
            "parse_tko_names": {**_S_BOOLEAN, "default": False, "description": "derive replicate (last char) and time (int of name[1:-1]) from sample names like T18A"},
            "exclude": {**_S_ARRAY_S, "description": "drop samples by obs value, e.g. replicate=0"},
            "condit_1": {**_S_STRING, "description": "regex for condition-1 samples"},
            "condit_2": {**_S_STRING, "description": "regex for condition-2 samples"},
            "control": {**_S_STRING, "description": "regex for control samples"},
            "targets": {**_S_ARRAY_S},
            "target_col": {**_S_STRING, "default": "target"},
            "gtf": {**_S_STRING},
            "gene": {**_S_STRING},
            "no_plots": {**_S_BOOLEAN, "default": False},
            "label": {**_S_STRING}},
         "required": ["condit_1", "condit_2", "control"]},
        cli=["perturb-tools", "sort-lfc"],
        flag_map={"screen": "--screen", "counts": "--counts", "guides": "--guides", "samples": "--samples", "parse_tko_names": "--parse-tko-names", "exclude": "--exclude", "condit_1": "--condit-1", "condit_2": "--condit-2", "control": "--control", "targets": "--targets", "target_col": "--target-col", "gtf": "--gtf", "gene": "--gene", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'exclude', 'targets'},
        bool_flags={'no_plots', 'parse_tko_names'},
    ),

    _T(
        "perturb_tools_qc",
        "Produces the perturb-tools sample / guide quality report of a pooled screen: per-sample total reads, median count, zero-count guides and Gini coefficient of guide coverage, the pairwise sample count correlation (mean / median per sample), replicate LFC correlation (all guides and positive controls) when cond1/cond2 are given, and outlier (jackpot) guides whose RPM exceeds mad_z_thres x the condition median and abs_rpm_thres. Writes sample_qc.tsv, correlation tables, outlier_guides.tsv, figures and a report.",
        {"type": "object", "properties": {
            "screen": {**_S_STRING, "description": ".h5ad or screen CSV directory written by perturb-tools (alternative to counts)"},
            "counts": {**_S_STRING, "description": "guides x samples count table (TSV/CSV; first column guide id; text columns become guide annotations, e.g. TKO readcount files)"},
            "guides": {**_S_STRING, "description": "optional guide annotation table (first column guide id)"},
            "samples": {**_S_STRING, "description": "optional sample annotation table (column 'sample' or first column)"},
            "parse_tko_names": {**_S_BOOLEAN, "default": False, "description": "derive replicate (last char) and time (int of name[1:-1]) from sample names like T18A"},
            "exclude": {**_S_ARRAY_S, "description": "drop samples by obs value, e.g. replicate=0"},
            "cond1": {**_S_STRING},
            "cond2": {**_S_STRING},
            "rep_col": {**_S_STRING, "default": "replicate"},
            "compare_col": {**_S_STRING, "default": "time"},
            "target_col": {**_S_STRING, "default": "target"},
            "pos_ctrl": {**_S_STRING},
            "outlier_cond_col": {**_S_STRING, "description": "obs column defining conditions for outlier guides (e.g. time)"},
            "mad_z_thres": {**_S_NUMBER, "default": 5.0},
            "abs_rpm_thres": {**_S_NUMBER, "default": 10000.0},
            "method": {**_S_STRING, "default": "spearman"},
            "no_plots": {**_S_BOOLEAN, "default": False},
            "label": {**_S_STRING}}},
        cli=["perturb-tools", "qc"],
        flag_map={"screen": "--screen", "counts": "--counts", "guides": "--guides", "samples": "--samples", "parse_tko_names": "--parse-tko-names", "exclude": "--exclude", "cond1": "--cond1", "cond2": "--cond2", "rep_col": "--rep-col", "compare_col": "--compare-col", "target_col": "--target-col", "pos_ctrl": "--pos-ctrl", "outlier_cond_col": "--outlier-cond-col", "mad_z_thres": "--mad-z-thres", "abs_rpm_thres": "--abs-rpm-thres", "method": "--method", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'exclude'},
        bool_flags={'no_plots', 'parse_tko_names'},
    ),

    _T(
        "perturb_tools_to_mageck",
        "Exports a screen as a MAGeCK count table (sgRNA, gene, one integer column per sample; spaces in sgRNA names replaced, guides without a target dropped) and prints the mageck test command to run next. Writes mageck_count.txt and a report.",
        {"type": "object", "properties": {
            "screen": {**_S_STRING, "description": ".h5ad or screen CSV directory written by perturb-tools (alternative to counts)"},
            "counts": {**_S_STRING, "description": "guides x samples count table (TSV/CSV; first column guide id; text columns become guide annotations, e.g. TKO readcount files)"},
            "guides": {**_S_STRING, "description": "optional guide annotation table (first column guide id)"},
            "samples": {**_S_STRING, "description": "optional sample annotation table (column 'sample' or first column)"},
            "parse_tko_names": {**_S_BOOLEAN, "default": False, "description": "derive replicate (last char) and time (int of name[1:-1]) from sample names like T18A"},
            "exclude": {**_S_ARRAY_S, "description": "drop samples by obs value, e.g. replicate=0"},
            "target_col": {**_S_STRING, "default": "target_id", "description": "guide feature with the gene name"},
            "sgrna_col": {**_S_STRING},
            "sample_prefix": {**_S_STRING, "default": ""},
            "label": {**_S_STRING}},
         "required": ["target_col"]},
        cli=["perturb-tools", "to-mageck"],
        flag_map={"screen": "--screen", "counts": "--counts", "guides": "--guides", "samples": "--samples", "parse_tko_names": "--parse-tko-names", "exclude": "--exclude", "target_col": "--target-col", "sgrna_col": "--sgrna-col", "sample_prefix": "--sample-prefix", "label": "--label"},
        flag_repeat={'exclude'},
        bool_flags={'parse_tko_names'},
    ),

    _T(
        "perturb_tools_read_poolq",
        "Reads a Broad PoolQ output directory (counts.txt or expected-counts.txt, lognormalized-counts.txt, barcode-counts.txt, unexpected-sequences.txt, correlation.txt, quality.txt, runinfo.txt) into a perturb-tools screen, parses the quality report into metadata / per-sample-barcode / common-barcode tables, merges optional sample metadata on Condition, and checks PoolQ's log-normalisation against log2(RPM+1). Writes the screen (.h5ad + CSV dir), the PoolQ tables as TSV and a report.",
        {"type": "object", "properties": {
            "poolq_dir": {**_S_STRING},
            "sample_metadata": {**_S_STRING},
            "merge_on": {**_S_STRING, "default": "Condition"},
            "label": {**_S_STRING}},
         "required": ["poolq_dir"]},
        cli=["perturb-tools", "read-poolq"],
        flag_map={"poolq_dir": "--poolq-dir", "sample_metadata": "--sample-metadata", "merge_on": "--merge-on", "label": "--label"},
    ),

    _T(
        "perturb_tools_annotate_guides",
        "Annotates a guide table (barcode, barcode_id): protospacer from 20-nt barcodes, target region by direct barcode_id->target pairs or name matching, and, given target regions (target, Chromosome, Start, End) and a reference FASTA, the strand and genomic Start/End/center of each protospacer; optional GC content and homopolymer flags. Writes guides_annotated.tsv (usable for sort_lfc position plots) and a report.",
        {"type": "object", "properties": {
            "guides_table": {**_S_STRING},
            "regions": {**_S_STRING},
            "fasta": {**_S_STRING},
            "units": {**_S_STRING, "default": "bp", "description": "region units bp or mb"},
            "direct_pairs": {**_S_STRING},
            "annotations": {**_S_ARRAY_S},
            "biochem": {**_S_BOOLEAN, "default": False},
            "label": {**_S_STRING}},
         "required": ["guides_table"]},
        cli=["perturb-tools", "annotate-guides"],
        flag_map={"guides_table": "--guides-table", "regions": "--regions", "fasta": "--fasta", "units": "--units", "direct_pairs": "--direct-pairs", "annotations": "--annotations", "biochem": "--biochem", "label": "--label"},
        flag_repeat={'annotations'},
        bool_flags={'biochem'},
    ),

    _T(
        "perturb_tools_design_library",
        "Designs a tiling sgRNA library (perturb-tools experimental design): takes exons of a gene from a GTF (overlapping/engulfing fragments resolved) or a regions table, scans both strands of the reference FASTA for the PAM (default NGG), and reports each guide's protospacer, PAM, strand, PAM position, distance to the region centre and 10-nt-flanked context, flags BsmbI motifs (filtered by default) and adds GC content and homopolymer flags. Writes guide_library.tsv, bsmbi_flagged_guides.tsv, regions.tsv and a report.",
        {"type": "object", "properties": {
            "fasta": {**_S_STRING},
            "gtf": {**_S_STRING},
            "gene": {**_S_STRING},
            "chrom": {**_S_STRING},
            "regions": {**_S_STRING},
            "pam": {**_S_STRING, "default": "NGG"},
            "widen": {**_S_INTEGER, "default": 0},
            "flank": {**_S_INTEGER, "default": 10},
            "keep_bsmbi": {**_S_BOOLEAN, "default": False},
            "label": {**_S_STRING}},
         "required": ["fasta"]},
        cli=["perturb-tools", "design-library"],
        flag_map={"fasta": "--fasta", "gtf": "--gtf", "gene": "--gene", "chrom": "--chrom", "regions": "--regions", "pam": "--pam", "widen": "--widen", "flank": "--flank", "keep_bsmbi": "--keep-bsmbi", "label": "--label"},
        bool_flags={'keep_bsmbi'},
    ),

    _T(
        "igvf_sc_barcode_revcomp_detect",
        "Detects whether a barcode FASTQ carries onlist barcodes forward or reverse-complemented (IGVF single-cell pipeline barcode_revcomp_detect). Reads the first num_reads reads, compares sequence[offset:] with the onlist and its reverse complement, and requires >= threshold of num_reads to match. Writes the QC text and the chromap read-format string bc:<offset>:-1[:-].",
        {"type": "object", "properties": {
            "fastq": {**_S_STRING, "description": "Barcode FASTQ (gz or plain)"},
            "onlist": {**_S_STRING, "description": "Barcode onlist / whitelist"},
            "offset": {**_S_INTEGER, "default": 0, "description": "Barcode start in the read (10x 0, 10x multiome ATAC 8)"},
            "num_reads": {**_S_INTEGER, "default": 100000, "description": "Reads to inspect (denominator of the match proportion)"},
            "threshold": {**_S_NUMBER, "default": 0.45, "description": "Minimum match proportion"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["fastq", "onlist"]},
        cli=["igvf-sc-pipeline", "barcode-revcomp-detect"],
        flag_map={"fastq": "--fastq", "onlist": "--onlist", "offset": "--offset", "num_reads": "--num-reads", "threshold": "--threshold", "label": "--label"},
    ),

    _T(
        "igvf_sc_correct_fastq",
        "Corrects SHARE-seq round-1/2/3 cell barcodes from the last 99 bp of read 2 (exact, 1-mismatch, +/-1 bp shift lookups) and tallies poly-G reads, as in IGVF correct_fastq.py. Writes corrected read 1 / read 2 / barcode FASTQs with barcodes (and UMI for RNA) in the read names, plus <prefix>_barcode_qc.txt (match, mismatch, poly_G_barcode, poly_G_in_first_10bp).",
        {"type": "object", "properties": {
            "read1": {**_S_STRING, "description": "Read 1 FASTQ"},
            "read2": {**_S_STRING, "description": "Read 2 FASTQ (barcodes in its last 99 bp)"},
            "whitelist": {**_S_STRING, "description": "R1R2R3 24-mer combinations, one per line"},
            "sample_type": {**_S_STRING, "description": "ATAC or RNA"},
            "prefix": {**_S_STRING, "description": "Library prefix for outputs"},
            "pkr": {**_S_STRING, "description": "Optional PKR / subpool tag added to the barcode"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "Reproduce upstream's shared r1/r2/r3 barcode set"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["read1", "read2", "whitelist", "sample_type", "prefix"]},
        cli=["igvf-sc-pipeline", "correct-fastq"],
        flag_map={"read1": "--read1", "read2": "--read2", "whitelist": "--whitelist", "sample_type": "--sample-type", "prefix": "--prefix", "pkr": "--pkr", "upstream_compat": "--upstream-compat", "label": "--label"},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "igvf_sc_trim_fastq",
        "Trims dovetailed (overlapping) ATAC read pairs: finds reverse_complement(read2[0:20]) in read 1 exactly or with Levenshtein distance <= 1 and cuts both reads at that position + 20. Writes trimmed FASTQs and a trimming-stats table (total, untrimmed, trimmed, %trimmed).",
        {"type": "object", "properties": {
            "read1": {**_S_STRING, "description": "Read 1 FASTQ"},
            "read2": {**_S_STRING, "description": "Read 2 FASTQ"},
            "prefix": {**_S_STRING, "default": "sample", "description": "Output prefix"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["read1", "read2"]},
        cli=["igvf-sc-pipeline", "trim-fastq"],
        flag_map={"read1": "--read1", "read2": "--read2", "prefix": "--prefix", "label": "--label"},
    ),

    _T(
        "igvf_sc_tss_enrichment",
        "Computes the bulk TSS enrichment score and profile plus per-barcode ArchR-style TSS enrichment from an ATAC fragment file and a TSS bed (IGVF compute_tss_enrichment_bulk). Uses tabix via pysam when the file is indexed, otherwise reads gzip/plain fragments with pandas. Writes <prefix>.tss_score_bulk.txt, the bulk profile TSV/PNG and <prefix>.tss_enrichment_barcode_stats.tsv (barcode, reads_unique, reads_promoter, reads_tss, tss_enrichment).",
        {"type": "object", "properties": {
            "fragments": {**_S_STRING, "description": "Fragment file (chr, start, end, barcode)"},
            "regions": {**_S_STRING, "description": "TSS bed (chr, TSS, end, ..., strand)"},
            "flank": {**_S_INTEGER, "default": 2000, "description": "Bases each side of the TSS"},
            "window": {**_S_INTEGER, "default": 20, "description": "Smoothing window"},
            "strand_col": {**_S_INTEGER, "default": 4, "description": "1-based strand column"},
            "prefix": {**_S_STRING, "default": "sample", "description": "Output prefix"},
            "engine": {**_S_STRING, "default": "auto", "description": "auto, tabix or pandas"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip the PNG"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["fragments", "regions"]},
        cli=["igvf-sc-pipeline", "tss-enrichment"],
        flag_map={"fragments": "--fragments", "regions": "--regions", "flank": "--flank", "window": "--window", "strand_col": "--strand-col", "prefix": "--prefix", "engine": "--engine", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "igvf_sc_snapatac2_tsse",
        "Runs snapatac2 import, TSS enrichment and FRiP (tss_frac, promoter_frac) on a fragment file, as in upstream snapatac2-tss-enrichment.py; skipped with a message when snapatac2 is not installed. Writes per-barcode stats TSV, TSSe plots and an h5ad.",
        {"type": "object", "properties": {
            "fragments": {**_S_STRING, "description": "Fragment file"},
            "gtf": {**_S_STRING, "description": "Compressed GTF"},
            "chrom_sizes": {**_S_STRING, "description": "Chromosome sizes TSV"},
            "tss_bed": {**_S_STRING, "description": "TSS bed"},
            "promoter_bed": {**_S_STRING, "description": "Promoter bed"},
            "min_frag_cutoff": {**_S_INTEGER, "default": 100, "description": "Minimum fragments per barcode"},
            "prefix": {**_S_STRING, "default": "sample", "description": "Output prefix"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip plots"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["fragments", "gtf", "chrom_sizes", "tss_bed", "promoter_bed"]},
        cli=["igvf-sc-pipeline", "snapatac2-tsse"],
        flag_map={"fragments": "--fragments", "gtf": "--gtf", "chrom_sizes": "--chrom-sizes", "tss_bed": "--tss-bed", "promoter_bed": "--promoter-bed", "min_frag_cutoff": "--min-frag-cutoff", "prefix": "--prefix", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "igvf_sc_rna_qc_metrics",
        "Extracts per-barcode RNA QC metrics (total_counts, genes) from a kallisto|bustools h5ad; for the nac workflow the mature, ambiguous and nascent layers are summed, as in qc_rna_extract_metrics.py. Optionally suffixes barcodes with _<subpool>. Writes rna_barcode_metadata.tsv.",
        {"type": "object", "properties": {
            "h5ad": {**_S_STRING, "description": "kb h5ad"},
            "kb_workflow": {**_S_STRING, "default": "standard", "description": "standard or nac"},
            "subpool": {**_S_STRING, "default": "none", "description": "Subpool suffix or 'none'"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["h5ad"]},
        cli=["igvf-sc-pipeline", "rna-qc-metrics"],
        flag_map={"h5ad": "--h5ad", "kb_workflow": "--kb-workflow", "subpool": "--subpool", "label": "--label"},
    ),

    _T(
        "igvf_sc_mtx_to_h5ad",
        "Converts a Matrix Market count matrix plus barcode and gene lists into an h5ad (write_h5ad_from_mtx.py). Writes output.h5ad (or --out).",
        {"type": "object", "properties": {
            "mtx": {**_S_STRING, "description": "Matrix Market file (cells x genes)"},
            "barcodes": {**_S_STRING, "description": "Barcodes, one per line"},
            "genes": {**_S_STRING, "description": "Genes, one per line"},
            "out": {**_S_STRING, "description": "Output h5ad path"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["mtx", "barcodes", "genes"]},
        cli=["igvf-sc-pipeline", "mtx-to-h5ad"],
        flag_map={"mtx": "--mtx", "barcodes": "--barcodes", "genes": "--genes", "out": "--out", "label": "--label"},
    ),

    _T(
        "igvf_sc_modify_barcode_h5",
        "Appends _<suffix> to every barcode in /obs/barcode of an h5ad in place (modify_barcode_h5.py). The file stays readable by anndata.",
        {"type": "object", "properties": {
            "h5ad": {**_S_STRING, "description": "h5ad to modify in place"},
            "suffix": {**_S_STRING, "description": "Suffix (subpool) without the underscore"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["h5ad", "suffix"]},
        cli=["igvf-sc-pipeline", "modify-barcode-h5"],
        flag_map={"h5ad": "--h5ad", "suffix": "--suffix", "label": "--label"},
    ),

    _T(
        "igvf_sc_joint_qc",
        "Joint RNA x ATAC cell calling for multiome data (joint_cell_plotting.py): merges RNA metrics (UMIs, genes) and ATAC metrics (fragments = reads/2, TSS enrichment) per barcode and labels each barcode both / RNA only / ATAC only / neither by the four cutoffs. Writes joint_barcode_metadata.csv and scatter + density plots.",
        {"type": "object", "properties": {
            "rna_metrics": {**_S_STRING, "description": "RNA metrics TSV (barcode, total_counts, genes)"},
            "atac_metrics": {**_S_STRING, "description": "ATAC metrics TSV (barcode, reads, ..., TSSe in column 5)"},
            "remove_low_yielding_cells": {**_S_INTEGER, "default": 10, "description": "Drop barcodes below this many UMIs / fragments"},
            "min_umis": {**_S_INTEGER, "default": 100, "description": "RNA UMI cutoff"},
            "min_genes": {**_S_INTEGER, "default": 200, "description": "RNA gene cutoff"},
            "min_tss": {**_S_NUMBER, "default": 4, "description": "ATAC TSS enrichment cutoff"},
            "min_frags": {**_S_INTEGER, "default": 100, "description": "ATAC fragment cutoff"},
            "pkr": {**_S_STRING, "description": "Sample / PKR name for titles"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip plots"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["rna_metrics", "atac_metrics"]},
        cli=["igvf-sc-pipeline", "joint-qc"],
        flag_map={"rna_metrics": "--rna-metrics", "atac_metrics": "--atac-metrics", "remove_low_yielding_cells": "--remove-low-yielding-cells", "min_umis": "--min-umis", "min_genes": "--min-genes", "min_tss": "--min-tss", "min_frags": "--min-frags", "pkr": "--pkr", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "igvf_sc_barcode_rank",
        "Finds the elbow and knee of a barcode-rank curve (barcode_rank_functions.R): farthest point from the end-to-end line on (rank, log10 count), then a smooth.spline(spar=1) second-derivative slice to locate the knee among the top-ranked barcodes. Works on any column of a metrics table. Writes the ranked TSV, a plot and the points in summary.json.",
        {"type": "object", "properties": {
            "metrics": {**_S_STRING, "description": "Metrics table (TSV/CSV)"},
            "column": {**_S_STRING, "default": "1", "description": "Column name or 0-based index"},
            "cutoff": {**_S_NUMBER, "default": 10, "description": "Minimum value kept"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip plot"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["metrics"]},
        cli=["igvf-sc-pipeline", "barcode-rank"],
        flag_map={"metrics": "--metrics", "column": "--column", "cutoff": "--cutoff", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "igvf_sc_atac_qc_plots",
        "ATAC fragment barcode-rank plots with elbow/knee (atac_qc_plots.R). Reads the fragment column (unique / reads_unique / fragments) of a barcode metadata table and writes the two-panel plot and the elbow/knee ranks.",
        {"type": "object", "properties": {
            "metrics": {**_S_STRING, "description": "ATAC barcode metadata"},
            "fragment_cutoff": {**_S_NUMBER, "default": 10, "description": "Minimum fragments"},
            "column": {**_S_STRING, "description": "Override the fragment column"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip plot"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["metrics"]},
        cli=["igvf-sc-pipeline", "atac-qc-plots"],
        flag_map={"metrics": "--metrics", "fragment_cutoff": "--fragment-cutoff", "column": "--column", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "igvf_sc_rna_qc_plots",
        "RNA QC plots (rna_qc_plots.R): UMI and gene barcode-rank plots with elbow/knee and a genes-vs-UMIs scatter from rna_barcode_metadata.tsv. Writes PNGs and the elbow/knee ranks in summary.json.",
        {"type": "object", "properties": {
            "metrics": {**_S_STRING, "description": "RNA metrics TSV (total_counts, genes)"},
            "umi_cutoff": {**_S_NUMBER, "default": 10, "description": "Minimum UMIs"},
            "gene_cutoff": {**_S_NUMBER, "default": 10, "description": "Minimum genes"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip plots"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["metrics"]},
        cli=["igvf-sc-pipeline", "rna-qc-plots"],
        flag_map={"metrics": "--metrics", "umi_cutoff": "--umi-cutoff", "gene_cutoff": "--gene-cutoff", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "igvf_sc_insert_size_hist",
        "Parses a Picard CollectInsertSizeMetrics histogram and plots it (plot_insert_size_hist.py). Writes the histogram TSV and PNG with the modal insert size.",
        {"type": "object", "properties": {
            "histogram": {**_S_STRING, "description": "Picard histogram text file"},
            "pkr": {**_S_STRING, "description": "Sample name for the title"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip plot"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["histogram"]},
        cli=["igvf-sc-pipeline", "insert-size-hist"],
        flag_map={"histogram": "--histogram", "pkr": "--pkr", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "igvf_sc_tenx_barcode_map",
        "Builds the 10x multiome ATAC-to-RNA barcode conversion dictionary for chromap --barcode-translate (10x_create_barcode_mapping.wdl): rows <rna> <revcomp(atac)> then <rna> <atac>, only when both onlists have equal length. Writes barcode_conversion_dict.tsv.",
        {"type": "object", "properties": {
            "atac_onlist": {**_S_STRING, "description": "ATAC onlist (e.g. 737K-arc-v1 ATAC)"},
            "rna_onlist": {**_S_STRING, "description": "RNA onlist (same order)"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["atac_onlist", "rna_onlist"]},
        cli=["igvf-sc-pipeline", "tenx-barcode-map"],
        flag_map={"atac_onlist": "--atac-onlist", "rna_onlist": "--rna-onlist", "label": "--label"},
    ),

    _T(
        "igvf_sc_log_atac",
        "Turns a chromap alignment log plus barcode summary CSV into <prefix>_qc_metrics.json (task_log_atac.wdl): 'Number of ...' counts, '#' statistics and percentage_duplicates = 100*dups/(total - unmapped - lowmapq).",
        {"type": "object", "properties": {
            "alignment_log": {**_S_STRING, "description": "chromap log (<prefix>.log.txt)"},
            "barcode_summary": {**_S_STRING, "description": "chromap --summary CSV"},
            "prefix": {**_S_STRING, "default": "sample", "description": "Output prefix"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["alignment_log", "barcode_summary"]},
        cli=["igvf-sc-pipeline", "log-atac"],
        flag_map={"alignment_log": "--alignment-log", "barcode_summary": "--barcode-summary", "prefix": "--prefix", "label": "--label"},
    ),

    _T(
        "igvf_sc_kb_count",
        "Builds and (if kb is on PATH) runs the IGVF kallisto|bustools count command (atomic-workflows run_kallisto quantify nac): interleaved lanes, nac index files, --sum=total, optional replacement list; then appends the subpool to the h5ad barcodes and moves it to <out>.h5ad. Without kb the exact command is printed and the tool exits 0.",
        {"type": "object", "properties": {
            "read1": {**_S_ARRAY_S, "description": "Read 1 FASTQs, one per lane"},
            "read2": {**_S_ARRAY_S, "description": "Read 2 FASTQs, one per lane"},
            "read_barcode": {**_S_ARRAY_S, "description": "Optional barcode FASTQs per lane"},
            "index_dir": {**_S_STRING, "description": "Extracted kb index folder"},
            "read_format": {**_S_STRING, "description": "kb -x technology string"},
            "onlist": {**_S_STRING, "description": "Barcode onlist (.gz ok)"},
            "kb_mode": {**_S_STRING, "default": "nac", "description": "nac or standard"},
            "strand": {**_S_STRING, "default": "unstranded", "description": "forward, reverse or unstranded"},
            "replacement_list": {**_S_STRING, "description": "Optional barcode replacement list"},
            "subpool": {**_S_STRING, "default": "none", "description": "Subpool suffix or 'none'"},
            "threads": {**_S_INTEGER, "default": 4, "description": "Threads"},
            "output_dir": {**_S_STRING, "description": "kb output directory"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["read1", "read2", "index_dir", "read_format", "onlist"]},
        cli=["igvf-sc-pipeline", "kb-count"],
        flag_map={"read1": "--read1", "read2": "--read2", "read_barcode": "--read-barcode", "index_dir": "--index-dir", "read_format": "--read-format", "onlist": "--onlist", "kb_mode": "--kb-mode", "strand": "--strand", "replacement_list": "--replacement-list", "subpool": "--subpool", "threads": "--threads", "output_dir": "--output-dir", "label": "--label"},
        flag_repeat={'read1', 'read2', 'read_barcode'},
    ),

    _T(
        "igvf_sc_kb_index",
        "Builds and optionally runs the kb ref command for a nac (or standard) kallisto index from a genome FASTA and GTF (task_kb_index.wdl); prints the command when kb is missing.",
        {"type": "object", "properties": {
            "genome_fasta": {**_S_STRING, "description": "Genome FASTA"},
            "gtf": {**_S_STRING, "description": "Gene GTF"},
            "kb_mode": {**_S_STRING, "default": "nac", "description": "nac or standard"},
            "output_dir": {**_S_STRING, "description": "Index output directory"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["genome_fasta", "gtf"]},
        cli=["igvf-sc-pipeline", "kb-index"],
        flag_map={"genome_fasta": "--genome-fasta", "gtf": "--gtf", "kb_mode": "--kb-mode", "output_dir": "--output-dir", "label": "--label"},
    ),

    _T(
        "igvf_sc_chromap_align",
        "Builds and (if chromap is on PATH) runs the IGVF chromap alignment producing fragments or a BAM (task_chromap.wdl / task_chromap_bam.wdl) with the pipeline's fixed flags; appends the subpool to barcodes, then bgzips and tabix-indexes the fragments. Without chromap the exact command is printed and the tool exits 0.",
        {"type": "object", "properties": {
            "output": {**_S_STRING, "default": "fragments", "description": "fragments or bam"},
            "index_dir": {**_S_STRING, "description": "Extracted chromap index folder"},
            "read_format": {**_S_STRING, "description": "chromap --read-format (e.g. bc:0:15 or bc:8:-1:-)"},
            "reference_fasta": {**_S_STRING, "description": "Genome FASTA"},
            "onlist": {**_S_STRING, "description": "Barcode onlist"},
            "read1": {**_S_ARRAY_S, "description": "Read 1 FASTQs"},
            "read2": {**_S_ARRAY_S, "description": "Read 2 FASTQs"},
            "read_barcode": {**_S_ARRAY_S, "description": "Barcode FASTQs"},
            "barcode_translate": {**_S_STRING, "description": "10x conversion dictionary"},
            "subpool": {**_S_STRING, "default": "none", "description": "Subpool suffix or 'none'"},
            "threads": {**_S_INTEGER, "default": 8, "description": "Threads"},
            "prefix": {**_S_STRING, "default": "sample.atac", "description": "Output prefix"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["index_dir", "read_format", "reference_fasta", "onlist", "read1", "read2", "read_barcode"]},
        cli=["igvf-sc-pipeline", "chromap-align"],
        flag_map={"output": "--output", "index_dir": "--index-dir", "read_format": "--read-format", "reference_fasta": "--reference-fasta", "onlist": "--onlist", "read1": "--read1", "read2": "--read2", "read_barcode": "--read-barcode", "barcode_translate": "--barcode-translate", "subpool": "--subpool", "threads": "--threads", "prefix": "--prefix", "label": "--label"},
        flag_repeat={'read1', 'read2', 'read_barcode'},
    ),

    _T(
        "igvf_sc_chromap_index",
        "Builds and optionally runs the chromap index command for a genome FASTA (task_chromap_index.wdl); prints it when chromap is missing.",
        {"type": "object", "properties": {
            "genome_fasta": {**_S_STRING, "description": "Genome FASTA"},
            "output_dir": {**_S_STRING, "description": "Index output directory"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["genome_fasta"]},
        cli=["igvf-sc-pipeline", "chromap-index"],
        flag_map={"genome_fasta": "--genome-fasta", "output_dir": "--output-dir", "label": "--label"},
    ),

    _T(
        "igvf_sc_subpool_fragments",
        "Appends _<subpool> to the barcodes in an ATAC fragment file (column 4) and a chromap barcode summary (column 1), as run_chromap does, optionally re-bgzipping and tabix-indexing the fragments.",
        {"type": "object", "properties": {
            "subpool": {**_S_STRING, "description": "Subpool suffix"},
            "fragments": {**_S_STRING, "description": "Uncompressed fragment TSV"},
            "summary": {**_S_STRING, "description": "chromap barcode summary CSV"},
            "bgzip": {**_S_BOOLEAN, "default": False, "description": "bgzip + tabix afterwards"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["subpool"]},
        cli=["igvf-sc-pipeline", "subpool-fragments"],
        flag_map={"subpool": "--subpool", "fragments": "--fragments", "summary": "--summary", "bgzip": "--bgzip", "label": "--label"},
        bool_flags={'bgzip'},
    ),

    _T(
        "igvf_sc_genome_tsv",
        "Parses the pipeline genome TSV (fasta, kb_nac_idx_tar, chromap_idx_tar) and resolves the reference files, letting explicit paths override it.",
        {"type": "object", "properties": {
            "genome_tsv": {**_S_STRING, "description": "Genome TSV"},
            "fasta": {**_S_STRING, "description": "Override FASTA"},
            "kb_index": {**_S_STRING, "description": "Override kb index tarball"},
            "chromap_index": {**_S_STRING, "description": "Override chromap index tarball"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["genome_tsv"]},
        cli=["igvf-sc-pipeline", "genome-tsv"],
        flag_map={"genome_tsv": "--genome-tsv", "fasta": "--fasta", "kb_index": "--kb-index", "chromap_index": "--chromap-index", "label": "--label"},
    ),

    _T(
        "igvf_sc_check_inputs",
        "Classifies pipeline inputs as in task_check_inputs.wdl: gs:// kept, syn* fetched with synapse get, https fetched from the IGVF portal with IGVF_API_KEY/IGVF_SECRET_KEY from the environment, and local paths checked for existence. By default it only plans the downloads; download=true fetches them.",
        {"type": "object", "properties": {
            "paths": {**_S_ARRAY_S, "description": "Input paths / URIs"},
            "download": {**_S_BOOLEAN, "default": False, "description": "Actually download syn*/https inputs"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["paths"]},
        cli=["igvf-sc-pipeline", "check-inputs"],
        flag_map={"paths": "--paths", "download": "--download", "label": "--label"},
        flag_repeat={'paths'},
        bool_flags={'download'},
    ),

    _T(
        "igvf_sc_sample_fastqs",
        "Keeps the first N reads (default 10 million) of each FASTQ and re-gzips them under the same basename (task_sample_fastqs.wdl).",
        {"type": "object", "properties": {
            "fastqs": {**_S_ARRAY_S, "description": "FASTQ files"},
            "n_reads": {**_S_INTEGER, "default": 10000000, "description": "Reads to keep"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["fastqs"]},
        cli=["igvf-sc-pipeline", "sample-fastqs"],
        flag_map={"fastqs": "--fastqs", "n_reads": "--n-reads", "label": "--label"},
        flag_repeat={'fastqs'},
    ),

    _T(
        "igvf_sc_portal_download",
        "Downloads IGVF portal files over https using IGVF_API_KEY / IGVF_SECRET_KEY from the environment (download_with_credentials.py). Set dry_run to only print what would be fetched.",
        {"type": "object", "properties": {
            "urls": {**_S_ARRAY_S, "description": "Portal @@download URLs"},
            "dest": {**_S_STRING, "description": "Destination directory"},
            "dry_run": {**_S_BOOLEAN, "default": False, "description": "Print only"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["urls"]},
        cli=["igvf-sc-pipeline", "portal-download"],
        flag_map={"urls": "--urls", "dest": "--dest", "dry_run": "--dry-run", "label": "--label"},
        flag_repeat={'urls'},
        bool_flags={'dry_run'},
    ),

    _T(
        "igvf_sc_synapse_manifest",
        "Builds a Synapse upload manifest (path, parent, name, used, executed, activityName) from a pipeline results table, as upload-to-synapse-y3-results.py does, with provenance links to the WDL tasks. Dry run by default; execute=true resolves folders and entities on Synapse with SYNAPSE_AUTH_TOKEN.",
        {"type": "object", "properties": {
            "table": {**_S_STRING, "description": "Results TSV (Subpool, raw inputs and output columns)"},
            "project": {**_S_STRING, "description": "Synapse project / parent ID"},
            "local_root": {**_S_STRING, "description": "Local root replacing remote_root"},
            "remote_root": {**_S_STRING, "description": "Remote prefix, e.g. gs://"},
            "output": {**_S_STRING, "description": "Manifest path"},
            "execute": {**_S_BOOLEAN, "default": False, "description": "Resolve on Synapse"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["table", "project"]},
        cli=["igvf-sc-pipeline", "synapse-manifest"],
        flag_map={"table": "--table", "project": "--project", "local_root": "--local-root", "remote_root": "--remote-root", "output": "--output", "execute": "--execute", "label": "--label"},
        bool_flags={'execute'},
    ),

    _T(
        "igvf_sc_synapse_upload",
        "Uploads a Synapse manifest in batches grouped by parent (batch_upload_synapse.py), recording failed parents. Dry run unless execute=true and synapseclient + SYNAPSE_AUTH_TOKEN are available.",
        {"type": "object", "properties": {
            "manifest": {**_S_STRING, "description": "Manifest TSV"},
            "execute": {**_S_BOOLEAN, "default": False, "description": "Really upload"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}},
         "required": ["manifest"]},
        cli=["igvf-sc-pipeline", "synapse-upload"],
        flag_map={"manifest": "--manifest", "execute": "--execute", "label": "--label"},
        bool_flags={'execute'},
    ),

    _T(
        "igvf_sc_synapse_annotations",
        "Computes the Synapse file annotations (file_type, status, main_output, data_type, file_description) the IGVF pipeline assigns from file names (update-file-annotations.py); applies them under root_folder only with execute=true.",
        {"type": "object", "properties": {
            "names": {**_S_ARRAY_S, "description": "File names"},
            "names_file": {**_S_STRING, "description": "File with one name per line"},
            "root_folder": {**_S_STRING, "description": "Synapse folder to annotate"},
            "execute": {**_S_BOOLEAN, "default": False, "description": "Apply on Synapse"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}}},
        cli=["igvf-sc-pipeline", "synapse-annotations"],
        flag_map={"names": "--names", "names_file": "--names-file", "root_folder": "--root-folder", "execute": "--execute", "label": "--label"},
        flag_repeat={'names'},
        bool_flags={'execute'},
    ),

    _T(
        "igvf_sc_html_report",
        "Writes the pipeline's tabbed HTML summary (Joint / RNA / ATAC plots, statistics, logs) and matching CSV from PNGs, name,value metric files and log paths (write_html.py).",
        {"type": "object", "properties": {
            "images": {**_S_ARRAY_S, "description": "PNG files"},
            "logs": {**_S_ARRAY_S, "description": "Log paths / links"},
            "atac_metrics": {**_S_STRING, "description": "ATAC name,value CSV"},
            "rna_metrics": {**_S_STRING, "description": "RNA name,value CSV"},
            "prefix": {**_S_STRING, "default": "summary", "description": "Output prefix"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}}},
        cli=["igvf-sc-pipeline", "html-report"],
        flag_map={"images": "--images", "logs": "--logs", "atac_metrics": "--atac-metrics", "rna_metrics": "--rna-metrics", "prefix": "--prefix", "label": "--label"},
        flag_repeat={'logs', 'images'},
    ),

    _T(
        "igvf_sc_pipeline",
        "Plans (and runs where kb / chromap are installed) the whole IGVF single_cell_pipeline.wdl for a SHARE-seq or 10x multiome sample: genome TSV, input checks, 10x barcode mapping, kb count, chromap fragments + BAM, chromap log metrics, then RNA, ATAC (TSS enrichment) and joint QC plus an HTML summary when h5ad / fragment outputs exist. Accepts a Cromwell inputs JSON or flags; writes plan.json, report.md and all QC tables.",
        {"type": "object", "properties": {
            "inputs_json": {**_S_STRING, "description": "Cromwell inputs JSON with single_cell_pipeline.<key> entries"},
            "prefix": {**_S_STRING, "description": "Analysis-set prefix"},
            "subpool": {**_S_STRING, "description": "Subpool suffix"},
            "genome_tsv": {**_S_STRING, "description": "Genome TSV"},
            "atac_read1": {**_S_ARRAY_S, "description": "ATAC read 1 FASTQs"},
            "atac_read2": {**_S_ARRAY_S, "description": "ATAC read 2 FASTQs"},
            "fastq_barcode": {**_S_ARRAY_S, "description": "ATAC barcode FASTQs"},
            "atac_onlist": {**_S_STRING, "description": "ATAC onlist"},
            "atac_read_format": {**_S_STRING, "description": "chromap read format"},
            "rna_read1": {**_S_ARRAY_S, "description": "RNA read 1 FASTQs"},
            "rna_read2": {**_S_ARRAY_S, "description": "RNA read 2 FASTQs"},
            "rna_onlist": {**_S_STRING, "description": "RNA onlist"},
            "rna_read_format": {**_S_STRING, "description": "kb -x string"},
            "kb_mode": {**_S_STRING, "description": "nac or standard"},
            "create_onlist_mapping": {**_S_BOOLEAN, "default": False, "description": "Build the 10x ATAC<->RNA mapping"},
            "rna_h5ad": {**_S_STRING, "description": "Existing kb h5ad for QC"},
            "fragments": {**_S_STRING, "description": "Existing fragment file for QC"},
            "tss_bed": {**_S_STRING, "description": "TSS bed for ATAC QC"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip plots"},
            "label": {**_S_STRING, "description": "Run label for the output directory name"}}},
        cli=["igvf-sc-pipeline", "pipeline"],
        flag_map={"inputs_json": "--inputs-json", "prefix": "--prefix", "subpool": "--subpool", "genome_tsv": "--genome-tsv", "atac_read1": "--atac-read1", "atac_read2": "--atac-read2", "fastq_barcode": "--fastq-barcode", "atac_onlist": "--atac-onlist", "atac_read_format": "--atac-read-format", "rna_read1": "--rna-read1", "rna_read2": "--rna-read2", "rna_onlist": "--rna-onlist", "rna_read_format": "--rna-read-format", "kb_mode": "--kb-mode", "create_onlist_mapping": "--create-onlist-mapping", "rna_h5ad": "--rna-h5ad", "fragments": "--fragments", "tss_bed": "--tss-bed", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'atac_read2', 'fastq_barcode', 'rna_read1', 'atac_read1', 'rna_read2'},
        bool_flags={'create_onlist_mapping', 'no_plots'},
    ),

    _T(
        "gwas_e2g_setup",
        "Fetches the small EngreitzLab/GWAS_E2G_benchmarking resources from the pinned commit into Data/GWASE2G/resources: the UK Biobank variant key, the gene prioritisation (silver-standard + PoPS) table (~9.5 MB), the hg38 partition (~50 MB), the TSS reference and chromosome sizes. traits also fetches per-trait SuSiE variant lists ('all' = 94 traits, ~247 MB); synapse=true pulls the 1000G background SNPs (syn52264319, needs SYNAPSE_AUTH_TOKEN).",
        {"type": "object", "properties": {
            "traits": {**_S_ARRAY_S, "description": "UK Biobank traits whose variant lists to fetch, e.g. RBC MCV Lym, or 'all'"},
            "synapse": {**_S_BOOLEAN, "default": False, "description": "also download the background SNPs from Synapse"},
            "force": {**_S_BOOLEAN, "default": False, "description": "re-download files that exist"}}},
        cli=["gwas-e2g", "setup"],
        flag_map={"traits": "--traits"},
        flag_repeat={'traits'},
        bool_flags={'force', 'synapse'},
    ),

    _T(
        "gwas_e2g_validate_config",
        "Validates a GWAS E2G benchmark configuration (upstream config.yml and/or flags): methods present in the methods and predictions tables, biosample and trait groups defined and not clashing with biosample/trait names, every comparisons-table biosample and trait defined (ALL allowed), and input files present. Prints the resolved resource paths and, per method, its threshold (negated for inverse predictors), biosamples and usable biosample groups; exits non-zero on errors.",
        {"type": "object", "properties": {
            "config": {**_S_STRING, "description": "upstream-style config.yml (or .json): methodsTable, predictionsTable, methods, comparisonsTable, biosampleGroups, traitGroups, thresholds, resource paths; flags override it"},
            "variant_key": {**_S_STRING, "description": "TSV trait, variant_file (variant files: chr,start,end,rsid,pip,CredibleSet,trait); default Data/GWASE2G/resources/UKBB_variant_key.tsv"},
            "bg_variants": {**_S_STRING, "description": "1000G background SNP bed chr,start,end,rsid, no header (Synapse syn52264319)"},
            "partition": {**_S_STRING, "description": "PartitionCombined.bed (default: fetched resource)"},
            "chr_sizes": {**_S_STRING, "description": "chromosome sizes TSV (default: fetched resource)"},
            "tss": {**_S_STRING, "description": "TSS reference bed; 4th column = gene universe (default: fetched resource)"},
            "gene_prioritization_table": {**_S_STRING, "description": "silver-standard table: CredibleSet, Disease, TargetGene, truth, POPS.Score, PromoterDistanceToBestSNP (default: UKBiobank.ABCGene.anyabc.tsv)"},
            "threshold_pip": {**_S_NUMBER, "default": 0.1, "description": "variants kept when pip > this"},
            "trait_group": {**_S_ARRAY_S, "description": "trait groups as NAME=trait1,trait2 (merged and de-duplicated; ALL is added automatically)"},
            "threshold_pval": {**_S_NUMBER, "default": 0.05, "description": "alpha for CIs and significance"},
            "num_pops_genes": {**_S_INTEGER, "default": 2, "description": "PoPS rank cut-off"},
            "methods_table": {**_S_STRING, "description": "TSV method, boolean, inverse_predictor, pred_name_long, threshold, score_col, color"},
            "predictions_table": {**_S_STRING, "description": "TSV biosample + one column per method holding the prediction file (chr, start, end, TargetGene, score column); blank = no predictions"},
            "comparisons_table": {**_S_STRING, "description": "TSV name, biosample (biosample, biosample group or ALL), trait (trait, trait group or ALL)"},
            "methods": {**_S_ARRAY_S, "description": "methods to benchmark (default: config or all methods in both tables)"},
            "biosample_group": {**_S_ARRAY_S, "description": "biosample groups as NAME=b1,b2 (used for a method only when all members have predictions; ALL is automatic)"},
            "n_threshold_steps": {**_S_INTEGER, "default": 25, "description": "enrichment-recall curve steps"},
            "num_pred_genes": {**_S_INTEGER, "default": 2, "description": "predicted-gene rank cut-off for gene linking"},
            "plot_fixed_scale": {**_S_BOOLEAN, "default": False, "description": "same heatmap colour scale across methods"},
            "no_file_checks": {**_S_BOOLEAN, "default": False, "description": "skip existence checks of prediction and variant files"}}},
        cli=["gwas-e2g", "validate-config"],
        flag_map={"config": "--config", "variant_key": "--variant-key", "bg_variants": "--bg-variants", "partition": "--partition", "chr_sizes": "--chr-sizes", "tss": "--tss", "gene_prioritization_table": "--gene-prioritization-table", "threshold_pip": "--threshold-pip", "trait_group": "--trait-group", "threshold_pval": "--threshold-pval", "num_pops_genes": "--num-pops-genes", "methods_table": "--methods-table", "predictions_table": "--predictions-table", "comparisons_table": "--comparisons-table", "methods": "--methods", "biosample_group": "--biosample-group", "n_threshold_steps": "--n-threshold-steps", "num_pred_genes": "--num-pred-genes", "plot_fixed_scale": "--plot-fixed-scale"},
        flag_repeat={'methods', 'trait_group', 'biosample_group'},
        bool_flags={'plot_fixed_scale', 'no_file_checks'},
    ),

    _T(
        "gwas_e2g_variants",
        "Processes fine-mapped GWAS variants as the upstream process_variants rules do: per trait keeps pip > threshold_pip on the chrSizes chromosomes in distal noncoding sequence (partition ABC/AllPeaks/Other/OtherIntron), merges and de-duplicates trait groups plus ALL, and filters the 1000G background SNPs to the same partition. Writes a variants directory (filteredGWASVariants.merged.sorted.tsv.gz, bgVariants.distalNoncoding.bed.gz, counts) that gwas_e2g_run reuses via variants_dir, plus report.md and summary.json with nVariantsTotal per trait.",
        {"type": "object", "properties": {
            "config": {**_S_STRING, "description": "upstream-style config.yml (or .json): methodsTable, predictionsTable, methods, comparisonsTable, biosampleGroups, traitGroups, thresholds, resource paths; flags override it"},
            "variant_key": {**_S_STRING, "description": "TSV trait, variant_file (variant files: chr,start,end,rsid,pip,CredibleSet,trait); default Data/GWASE2G/resources/UKBB_variant_key.tsv"},
            "bg_variants": {**_S_STRING, "description": "1000G background SNP bed chr,start,end,rsid, no header (Synapse syn52264319)"},
            "partition": {**_S_STRING, "description": "PartitionCombined.bed (default: fetched resource)"},
            "chr_sizes": {**_S_STRING, "description": "chromosome sizes TSV (default: fetched resource)"},
            "tss": {**_S_STRING, "description": "TSS reference bed; 4th column = gene universe (default: fetched resource)"},
            "gene_prioritization_table": {**_S_STRING, "description": "silver-standard table: CredibleSet, Disease, TargetGene, truth, POPS.Score, PromoterDistanceToBestSNP (default: UKBiobank.ABCGene.anyabc.tsv)"},
            "threshold_pip": {**_S_NUMBER, "default": 0.1, "description": "variants kept when pip > this"},
            "trait_group": {**_S_ARRAY_S, "description": "trait groups as NAME=trait1,trait2 (merged and de-duplicated; ALL is added automatically)"},
            "threshold_pval": {**_S_NUMBER, "default": 0.05, "description": "alpha for CIs and significance"},
            "num_pops_genes": {**_S_INTEGER, "default": 2, "description": "PoPS rank cut-off"},
            "label": {**_S_STRING, "default": "gwas_variants", "description": "run label"}}},
        cli=["gwas-e2g", "variants"],
        flag_map={"config": "--config", "variant_key": "--variant-key", "bg_variants": "--bg-variants", "partition": "--partition", "chr_sizes": "--chr-sizes", "tss": "--tss", "gene_prioritization_table": "--gene-prioritization-table", "threshold_pip": "--threshold-pip", "trait_group": "--trait-group", "threshold_pval": "--threshold-pval", "num_pops_genes": "--num-pops-genes", "label": "--label"},
        flag_repeat={'trait_group'},
    ),

    _T(
        "gwas_e2g_baseline",
        "Computes the gene-linking baselines of evaluate_baseline_predictors.R from a gene prioritisation table: PoPS (top num_pops_genes by POPS.Score per credible set) and distanceToTSS (top num_pops_genes by PromoterDistanceToBestSNP) precision and recall against the silver-standard truth genes, per trait and over all traits, with Agresti-Coull CIs; the table is restricted to the TSS gene universe. Outputs baseline/gene_linking/precisionRecall.byTrait.tsv.gz, report.md and summary.json.",
        {"type": "object", "properties": {
            "config": {**_S_STRING, "description": "upstream-style config.yml (or .json): methodsTable, predictionsTable, methods, comparisonsTable, biosampleGroups, traitGroups, thresholds, resource paths; flags override it"},
            "gene_prioritization_table": {**_S_STRING, "description": "silver-standard table: CredibleSet, Disease, TargetGene, truth, POPS.Score, PromoterDistanceToBestSNP (default: UKBiobank.ABCGene.anyabc.tsv)"},
            "tss": {**_S_STRING, "description": "TSS reference bed; 4th column = gene universe (default: fetched resource)"},
            "num_pops_genes": {**_S_INTEGER, "default": 2, "description": "PoPS rank cut-off"},
            "threshold_pval": {**_S_NUMBER, "default": 0.05, "description": "alpha for CIs and significance"},
            "label": {**_S_STRING, "default": "gwas_baseline", "description": "run label"}}},
        cli=["gwas-e2g", "baseline"],
        flag_map={"config": "--config", "gene_prioritization_table": "--gene-prioritization-table", "tss": "--tss", "num_pops_genes": "--num-pops-genes", "threshold_pval": "--threshold-pval", "label": "--label"},
    ),

    _T(
        "gwas_e2g_run",
        "Runs the full GWAS enhancer-gene benchmark (port of EngreitzLab/GWAS_E2G_benchmarking, the scE2G / ENCODE-rE2G GWAS benchmark) for a methods table, predictions table and comparisons table. Per method and trait x biosample (biosample groups and ALL included) it computes enrichment of fine-mapped distal-noncoding GWAS variants in thresholded predicted enhancers vs 1000G SNPs (log-RR CI, hypergeometric p, Bonferroni), recall with Agresti-Coull CIs, enhancer set sizes, a quantile threshold span and enrichment-recall curves, and gene-linking precision/recall against silver-standard credible-set genes (top num_pred_genes, alone and intersected with PoPS top genes), plus PoPS and distance-to-TSS baselines. It then writes the metric ranges, clustered heatmap matrices, per-comparison curve values and thresholded performance comparison tables, figures, report.md and summary.json under Docs/GWASE2G/<run>. Upstream quirks (SE formula, span over distinct scores, unthresholded single-biosample linking, unfiltered gene table, ALL double counting) are corrected unless upstream_compat.",
        {"type": "object", "properties": {
            "config": {**_S_STRING, "description": "upstream-style config.yml (or .json): methodsTable, predictionsTable, methods, comparisonsTable, biosampleGroups, traitGroups, thresholds, resource paths; flags override it"},
            "variant_key": {**_S_STRING, "description": "TSV trait, variant_file (variant files: chr,start,end,rsid,pip,CredibleSet,trait); default Data/GWASE2G/resources/UKBB_variant_key.tsv"},
            "bg_variants": {**_S_STRING, "description": "1000G background SNP bed chr,start,end,rsid, no header (Synapse syn52264319)"},
            "partition": {**_S_STRING, "description": "PartitionCombined.bed (default: fetched resource)"},
            "chr_sizes": {**_S_STRING, "description": "chromosome sizes TSV (default: fetched resource)"},
            "tss": {**_S_STRING, "description": "TSS reference bed; 4th column = gene universe (default: fetched resource)"},
            "gene_prioritization_table": {**_S_STRING, "description": "silver-standard table: CredibleSet, Disease, TargetGene, truth, POPS.Score, PromoterDistanceToBestSNP (default: UKBiobank.ABCGene.anyabc.tsv)"},
            "threshold_pip": {**_S_NUMBER, "default": 0.1, "description": "variants kept when pip > this"},
            "trait_group": {**_S_ARRAY_S, "description": "trait groups as NAME=trait1,trait2 (merged and de-duplicated; ALL is added automatically)"},
            "threshold_pval": {**_S_NUMBER, "default": 0.05, "description": "alpha for CIs and significance"},
            "num_pops_genes": {**_S_INTEGER, "default": 2, "description": "PoPS rank cut-off"},
            "methods_table": {**_S_STRING, "description": "TSV method, boolean, inverse_predictor, pred_name_long, threshold, score_col, color"},
            "predictions_table": {**_S_STRING, "description": "TSV biosample + one column per method holding the prediction file (chr, start, end, TargetGene, score column); blank = no predictions"},
            "comparisons_table": {**_S_STRING, "description": "TSV name, biosample (biosample, biosample group or ALL), trait (trait, trait group or ALL)"},
            "methods": {**_S_ARRAY_S, "description": "methods to benchmark (default: config or all methods in both tables)"},
            "biosample_group": {**_S_ARRAY_S, "description": "biosample groups as NAME=b1,b2 (used for a method only when all members have predictions; ALL is automatic)"},
            "n_threshold_steps": {**_S_INTEGER, "default": 25, "description": "enrichment-recall curve steps"},
            "num_pred_genes": {**_S_INTEGER, "default": 2, "description": "predicted-gene rank cut-off for gene linking"},
            "plot_fixed_scale": {**_S_BOOLEAN, "default": False, "description": "same heatmap colour scale across methods"},
            "variants_dir": {**_S_STRING, "description": "output directory of gwas_e2g_variants (skips variant processing)"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "reproduce the upstream quirks for number-for-number comparison with a pipeline run"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "gwas_e2g", "description": "run label"}}},
        cli=["gwas-e2g", "run"],
        flag_map={"config": "--config", "variant_key": "--variant-key", "bg_variants": "--bg-variants", "partition": "--partition", "chr_sizes": "--chr-sizes", "tss": "--tss", "gene_prioritization_table": "--gene-prioritization-table", "threshold_pip": "--threshold-pip", "trait_group": "--trait-group", "threshold_pval": "--threshold-pval", "num_pops_genes": "--num-pops-genes", "methods_table": "--methods-table", "predictions_table": "--predictions-table", "comparisons_table": "--comparisons-table", "methods": "--methods", "biosample_group": "--biosample-group", "n_threshold_steps": "--n-threshold-steps", "num_pred_genes": "--num-pred-genes", "plot_fixed_scale": "--plot-fixed-scale", "variants_dir": "--variants-dir", "label": "--label"},
        flag_repeat={'methods', 'trait_group', 'biosample_group'},
        bool_flags={'plot_fixed_scale', 'upstream_compat', 'no_plots'},
    ),

    _T(
        "gwas_e2g_plot",
        "Re-runs the visualisation stage of the GWAS E2G benchmark from an existing gwas_e2g_run directory: colour palette, metric ranges, clustered enrichment/recall and precision/recall heatmaps (Ward on 1 - correlation), enrichment-recall curves and thresholded performance comparison scatter plots per comparison, optionally with a fixed colour scale across methods. Writes into <run>/plots with report.md and summary.json.",
        {"type": "object", "properties": {
            "run_dir": {**_S_STRING, "description": "Docs/GWASE2G/<run> produced by gwas_e2g_run"},
            "variants_dir": {**_S_STRING, "description": "variants dir, needed only when the run used variants_dir"},
            "plot_fixed_scale": {**_S_BOOLEAN, "default": False, "description": "same heatmap colour scale across methods"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "reproduce upstream quirks in the comparison tables"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "tables only"}},
         "required": ["run_dir"]},
        cli=["gwas-e2g", "plot"],
        flag_map={"run_dir": "--run-dir", "variants_dir": "--variants-dir"},
        bool_flags={'plot_fixed_scale', 'upstream_compat', 'no_plots'},
    ),

    _T(
        "encode_re2g_models",
        "Lists the nine pretrained ENCODE-rE2G logistic-regression models embedded from the upstream pickles (DNase/ATAC x H3K27ac x intact Hi-C/megamap, plus the 45-feature extended model): features, coefficients, intercept and CRISPR 70%-recall threshold. Optionally exports each as a model directory (feature_table.tsv, threshold_<t>, model.json) usable by apply. Writes report.md and summary.json.",
        {"type": "object", "properties": {
            "export": {**_S_STRING, "description": "directory to write one model dir per model (default: the run dir)"},
            "label": {**_S_STRING, "default": "models", "description": "run label (output dir suffix)"}}},
        cli=["encode-re2g", "models"],
        flag_map={"export": "--export", "label": "--label"},
    ),

    _T(
        "encode_re2g_select_model",
        "Chooses the ENCODE-rE2G model for each biosample of an ABC biosample config exactly as upstream: <dhs|atac>[_h3k27ac]_<intact_hic|megamap|avg_hic|powerlaw> from default_accessibility_feature, the H3K27ac column and HiC_file/HiC_type, or the model_dir column. Rejects untrained combinations (powerlaw, avg_hic). Writes config_biosamples_models.tsv with model dirs and thresholds.",
        {"type": "object", "properties": {
            "biosample_config": {**_S_STRING, "description": "ABC biosample config TSV (biosample, DHS, ATAC, H3K27ac, default_accessibility_feature, HiC_file, HiC_type[, model_dir])"},
            "model_root": {**_S_STRING, "description": "directory of model folders (default: embedded models)"},
            "label": {**_S_STRING, "default": "select_model", "description": "run label (output dir suffix)"}},
         "required": ["biosample_config"]},
        cli=["encode-re2g", "select-model"],
        flag_map={"biosample_config": "--biosample-config", "model_root": "--model-root", "label": "--label"},
    ),

    _T(
        "encode_re2g_features",
        "Builds the genome-wide ENCODE-rE2G feature table from ABC outputs (EnhancerPredictionsAllPutative + EnhancerList): numCandidateEnhGene, numTSSEnhGene, numNearbyEnhancers/sumNearbyEnhancers (5 kb), ABC numerator/denominator, gene classes, external features joined by overlap or TargetGene, interaction terms, renames and NA fills from the models' feature tables. Writes the upstream intermediate files and genomewide_features.tsv.gz.",
        {"type": "object", "properties": {
            "abc_dir": {**_S_STRING, "description": "ABC biosample dir with Predictions/ and Neighborhoods/"},
            "abc_predictions": {**_S_STRING, "description": "EnhancerPredictionsAllPutative.tsv.gz (instead of abc_dir)"},
            "enhancer_list": {**_S_STRING, "description": "EnhancerList.txt (instead of abc_dir)"},
            "model": {**_S_ARRAY_S, "description": "embedded model names or model dirs whose features to build, e.g. dhs_intact_hic"},
            "feature_table": {**_S_ARRAY_S, "description": "extra feature_table.tsv files"},
            "external_features_config": {**_S_STRING, "description": "external features config TSV (input_col, source_col, aggregate_function, join_by, source_file)"},
            "tss": {**_S_STRING, "description": "TSS500bp reference bed (default: setup resources)"},
            "chr_sizes": {**_S_STRING, "description": "chromosome sizes TSV"},
            "gene_classes": {**_S_STRING, "description": "gene_promoter_class TSV"},
            "biosample": {**_S_STRING, "description": "biosample name for the output subfolder"},
            "dedupe_nearby": {**_S_BOOLEAN, "default": False, "description": "count each neighbouring element once (upstream intent; pretrained models used the duplicated count)"},
            "label": {**_S_STRING, "default": "features", "description": "run label (output dir suffix)"}}},
        cli=["encode-re2g", "features"],
        flag_map={"abc_dir": "--abc-dir", "abc_predictions": "--abc-predictions", "enhancer_list": "--enhancer-list", "model": "--model", "feature_table": "--feature-table", "external_features_config": "--external-features-config", "tss": "--tss", "chr_sizes": "--chr-sizes", "gene_classes": "--gene-classes", "biosample": "--biosample", "dedupe_nearby": "--dedupe-nearby", "label": "--label"},
        flag_repeat={'model', 'feature_table'},
        bool_flags={'dedupe_nearby'},
    ),

    _T(
        "encode_re2g_apply",
        "Scores a genome-wide ENCODE-rE2G feature table with one or more models: X = log(|x| + 0.01), ENCODE-rE2G.Score = logistic regression probability (embedded coefficients reproduce the upstream pickles to 1e-15). Thresholds at the model threshold keeping self-promoters, writes encode_e2g_predictions.tsv.gz, the thresholded table, an IGV bedpe and the per-prediction stats.",
        {"type": "object", "properties": {
            "features": {**_S_STRING, "description": "genomewide_features.tsv.gz from features"},
            "model": {**_S_ARRAY_S, "description": "embedded model names (e.g. dhs_intact_hic, atac_megamap, extended) or model dirs"},
            "threshold": {**_S_NUMBER, "description": "override the model threshold"},
            "exclude_self_promoter": {**_S_BOOLEAN, "default": False, "description": "drop self-promoters (include_self_promoter = False)"},
            "accessibility": {**_S_ARRAY_S, "description": "accessibility BAM(s) for num_sequencing_reads"},
            "biosample": {**_S_STRING, "description": "biosample name"},
            "label": {**_S_STRING, "default": "apply", "description": "run label (output dir suffix)"}},
         "required": ["features", "model"]},
        cli=["encode-re2g", "apply"],
        flag_map={"features": "--features", "model": "--model", "threshold": "--threshold", "exclude_self_promoter": "--exclude-self-promoter", "accessibility": "--accessibility", "biosample": "--biosample", "label": "--label"},
        flag_repeat={'accessibility', 'model'},
        bool_flags={'exclude_self_promoter'},
    ),

    _T(
        "encode_re2g_run",
        "Runs the full ENCODE-rE2G apply workflow for every biosample of an ABC biosample config: model selection, feature generation from each biosample's ABC outputs, scoring, thresholding, bedpe, stats and QC plots across biosamples. Writes one folder per biosample and model plus report.md and summary.json.",
        {"type": "object", "properties": {
            "biosample_config": {**_S_STRING, "description": "ABC biosample config (optional columns model_dir, ABC_directory, external_features_config)"},
            "abc_results": {**_S_STRING, "description": "root dir with one ABC output dir per biosample"},
            "model_root": {**_S_STRING, "description": "directory of model folders (default: embedded)"},
            "tss": {**_S_STRING, "description": "TSS500bp reference bed"},
            "chr_sizes": {**_S_STRING, "description": "chromosome sizes"},
            "gene_classes": {**_S_STRING, "description": "gene classes TSV"},
            "exclude_self_promoter": {**_S_BOOLEAN, "default": False, "description": "drop self-promoters"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "run", "description": "run label (output dir suffix)"}},
         "required": ["biosample_config"]},
        cli=["encode-re2g", "run"],
        flag_map={"biosample_config": "--biosample-config", "abc_results": "--abc-results", "model_root": "--model-root", "tss": "--tss", "chr_sizes": "--chr-sizes", "gene_classes": "--gene-classes", "exclude_self_promoter": "--exclude-self-promoter", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'exclude_self_promoter', 'no_plots'},
    ),

    _T(
        "encode_re2g_stats",
        "Computes the ENCODE-rE2G per-prediction statistics of a thresholded prediction table: number of enhancers, genes and links, genes with >= 1 enhancer, mean genes per enhancer, mean enhancers per gene (with and without promoters), mean log10 distance to TSS, mean element size, and sequencing reads from BAMs. Writes the _stats.tsv table.",
        {"type": "object", "properties": {
            "predictions": {**_S_STRING, "description": "thresholded predictions TSV(.gz)"},
            "accessibility": {**_S_ARRAY_S, "description": "accessibility BAM(s)"},
            "label": {**_S_STRING, "default": "stats", "description": "run label (output dir suffix)"}},
         "required": ["predictions"]},
        cli=["encode-re2g", "stats"],
        flag_map={"predictions": "--predictions", "accessibility": "--accessibility", "label": "--label"},
        flag_repeat={'accessibility'},
    ),

    _T(
        "encode_re2g_qc_plots",
        "Summarises many ENCODE-rE2G _stats.tsv files (one per biosample/model): distributions of every metric, metric versus sequencing depth, and top/bottom-5 outlier datasets, optionally split into cells and tissues with ENCODE metadata. Writes figures, stats_table.tsv and outlier_stats.tsv.",
        {"type": "object", "properties": {
            "stats": {**_S_ARRAY_S, "description": "*_stats.tsv files"},
            "encode_metadata": {**_S_STRING, "description": "ENCODE metadata TSV (DNase Experiment accession, Biosample type)"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "qc", "description": "run label (output dir suffix)"}},
         "required": ["stats"]},
        cli=["encode-re2g", "qc-plots"],
        flag_map={"stats": "--stats", "encode_metadata": "--encode-metadata", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'stats'},
        bool_flags={'no_plots'},
    ),

    _T(
        "encode_re2g_crispr_features",
        "Builds ENCODE-rE2G training data: overlaps a genome-wide feature table with CRISPR E-G pairs (EPCrisprBenchmark format) by gene and closed-interval overlap, aggregating features with each feature's aggregate_function, recomputes missing distanceToTSS, fills NAs, drops genes outside the TSS universe and Regulated NA. Writes the merged, missing and for_training tables.",
        {"type": "object", "properties": {
            "features": {**_S_STRING, "description": "genomewide_features.tsv.gz"},
            "crispr": {**_S_STRING, "description": "CRISPR benchmark TSV (EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz after setup)"},
            "model": {**_S_ARRAY_S, "description": "models whose feature tables define the features"},
            "feature_table": {**_S_ARRAY_S, "description": "feature_table.tsv files"},
            "tss": {**_S_STRING, "description": "TSS500bp reference bed"},
            "no_na_fill": {**_S_BOOLEAN, "default": False, "description": "keep NAs (upstream NAnotfilled)"},
            "dataset": {**_S_STRING, "description": "dataset name for file names"},
            "label": {**_S_STRING, "default": "crispr_features", "description": "run label (output dir suffix)"}},
         "required": ["features", "crispr"]},
        cli=["encode-re2g", "crispr-features"],
        flag_map={"features": "--features", "crispr": "--crispr", "model": "--model", "feature_table": "--feature-table", "tss": "--tss", "no_na_fill": "--no-na-fill", "dataset": "--dataset", "label": "--label"},
        flag_repeat={'model', 'feature_table'},
        bool_flags={'no_na_fill'},
    ),

    _T(
        "encode_re2g_train",
        "Trains an ENCODE-rE2G logistic regression on CRISPR training data: full model plus leave-one-chromosome-out models on log(|x| + 0.01) features (optional degree-2 polynomial terms, overridable sklearn parameters). Reports per-chromosome and pooled log loss, AUROC and AUPRC, coefficients and the 70%-recall threshold, and writes a model_dir usable by apply.",
        {"type": "object", "properties": {
            "crispr_features": {**_S_STRING, "description": "for_training.*.tsv.gz from crispr-features"},
            "feature_table": {**_S_STRING, "description": "feature_table.tsv"},
            "model": {**_S_STRING, "description": "use this embedded model's feature table"},
            "polynomial": {**_S_BOOLEAN, "default": False, "description": "degree-2 polynomial features"},
            "override_params": {**_S_STRING, "description": "dict string of LogisticRegression overrides, e.g. {'penalty': 'l2', 'C': 1}"},
            "dataset": {**_S_STRING, "description": "dataset name"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "reproduce the upstream log_loss-in-AUROC bug of the pooled row"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "train", "description": "run label (output dir suffix)"}},
         "required": ["crispr_features"]},
        cli=["encode-re2g", "train"],
        flag_map={"crispr_features": "--crispr-features", "feature_table": "--feature-table", "model": "--model", "polynomial": "--polynomial", "override_params": "--override-params", "dataset": "--dataset", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots', 'upstream_compat', 'polynomial'},
    ),

    _T(
        "encode_re2g_feature_selection",
        "Runs ENCODE-rE2G forward or backward sequential feature selection on CRISPR training data with leave-one-chromosome-out CV AUPRC, then refits the selected order with paired BCa bootstraps of delta AUPRC and delta precision at 70% recall and bootstrap p-values. Writes forward/backward_feature_selection.tsv and bar plots.",
        {"type": "object", "properties": {
            "crispr_features": {**_S_STRING, "description": "for_training.*.tsv.gz"},
            "feature_table": {**_S_STRING, "description": "feature_table.tsv"},
            "model": {**_S_STRING, "description": "use this model's feature table"},
            "direction": {**_S_STRING, "default": "forward", "description": "forward or backward"},
            "n_boot": {**_S_INTEGER, "default": 1000, "description": "bootstrap resamples"},
            "seed": {**_S_INTEGER, "default": 0, "description": "random seed"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "feature_selection", "description": "run label (output dir suffix)"}},
         "required": ["crispr_features"]},
        cli=["encode-re2g", "feature-selection"],
        flag_map={"crispr_features": "--crispr-features", "feature_table": "--feature-table", "model": "--model", "direction": "--direction", "n_boot": "--n-boot", "seed": "--seed", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "encode_re2g_permutation_importance",
        "Computes ENCODE-rE2G permutation feature importance: each feature is shuffled n_repeats times, the model is refit with leave-one-chromosome-out CV, and delta AUPRC and delta precision at 70% recall versus the full model are recorded. Writes permutation_feature_importance.tsv and plots.",
        {"type": "object", "properties": {
            "crispr_features": {**_S_STRING, "description": "for_training.*.tsv.gz"},
            "feature_table": {**_S_STRING, "description": "feature_table.tsv"},
            "model": {**_S_STRING, "description": "use this model's feature table"},
            "n_repeats": {**_S_INTEGER, "default": 20, "description": "permutations per feature"},
            "seed": {**_S_INTEGER, "default": 0, "description": "random seed"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "permutation_importance", "description": "run label (output dir suffix)"}},
         "required": ["crispr_features"]},
        cli=["encode-re2g", "permutation-importance"],
        flag_map={"crispr_features": "--crispr-features", "feature_table": "--feature-table", "model": "--model", "n_repeats": "--n-repeats", "seed": "--seed", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "encode_re2g_all_feature_sets",
        "Evaluates every non-empty subset of an ENCODE-rE2G feature table (2^n - 1 models, n < 14) with leave-one-chromosome-out CV and bootstrapped AUPRC with 95% CI, sorted best first. Writes all_feature_sets.tsv.",
        {"type": "object", "properties": {
            "crispr_features": {**_S_STRING, "description": "for_training.*.tsv.gz"},
            "feature_table": {**_S_STRING, "description": "feature_table.tsv"},
            "model": {**_S_STRING, "description": "use this model's feature table"},
            "n_boot": {**_S_INTEGER, "default": 1000, "description": "bootstrap resamples"},
            "seed": {**_S_INTEGER, "default": 0, "description": "random seed"},
            "label": {**_S_STRING, "default": "all_feature_sets", "description": "run label (output dir suffix)"}},
         "required": ["crispr_features"]},
        cli=["encode-re2g", "all-feature-sets"],
        flag_map={"crispr_features": "--crispr-features", "feature_table": "--feature-table", "model": "--model", "n_boot": "--n-boot", "seed": "--seed", "label": "--label"},
    ),

    _T(
        "encode_re2g_compare_models",
        "Compares trained ENCODE-rE2G models on their CRISPR cross-validation predictions (missing CRISPR pairs scored 0) and a distance-to-TSS baseline: bootstrapped AUPRC and precision at 70% recall with 95% CIs, threshold and fraction missing. Writes performance_across_models.tsv and bar plots.",
        {"type": "object", "properties": {
            "train_dirs": {**_S_ARRAY_S, "description": "run dirs of encode_re2g_train"},
            "crispr": {**_S_STRING, "description": "raw CRISPR table for the distance baseline"},
            "n_boot": {**_S_INTEGER, "default": 1000, "description": "bootstrap resamples"},
            "seed": {**_S_INTEGER, "default": 0, "description": "random seed"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "compare_models", "description": "run label (output dir suffix)"}},
         "required": ["train_dirs"]},
        cli=["encode-re2g", "compare-models"],
        flag_map={"train_dirs": "--train-dirs", "crispr": "--crispr", "n_boot": "--n-boot", "seed": "--seed", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'train_dirs'},
        bool_flags={'no_plots'},
    ),

    _T(
        "encode_re2g_verify_upstream",
        "Re-scores the upstream ENCODE_rE2G CircleCI expected output (K562 chr22, dhs_intact_hic, 1.68M pairs) with the embedded coefficients and compares scores, the thresholded link set and the stats table; with an upstream clone it also checks every model pickle against the embedded coefficients. Needs `setup --test-data` or a clone.",
        {"type": "object", "properties": {
            "upstream_dir": {**_S_STRING, "description": "a clone of EngreitzLab/ENCODE_rE2G"},
            "expected_dir": {**_S_STRING, "description": "dir with the three expected files"},
            "label": {**_S_STRING, "default": "verify_upstream", "description": "run label (output dir suffix)"}}},
        cli=["encode-re2g", "verify-upstream"],
        flag_map={"upstream_dir": "--upstream-dir", "expected_dir": "--expected-dir", "label": "--label"},
    ),

    _T(
        "bulk_crispr_run",
        "Runs the whole Perturb-seq chain of IGVF-CRISPR/bulk_crispr_pipeline on kb-count output directories: per-lane cell QC (knee, mito, Scrublet doublets), guide binarisation, optional MULTI-seq demultiplexing, cis gene selection (+/- 1 Mb), SCEPTRE-style conditional-resampling tests per guide x gene, BH and Fisher aggregation, genome-browser tracks. Inputs: GTF, <sample>_L<lane>_ks_transcripts_out and _ks_guide_out dirs (or guide FASTQs). Outputs report.md, summary.json, result_guides.tsv, result_elements.tsv and a results modality bundle under Docs/BulkCRISPR.",
        {"type": "object", "properties": {
            "gtf": {**_S_STRING, "description": "GTF with gene records (gene_id, gene_name)"},
            "rna_dirs": {**_S_ARRAY_S, "description": "kb count cDNA output dirs, one per lane (<sample>_L<lane>_ks_transcripts_out)"},
            "guide_dirs": {**_S_ARRAY_S, "description": "kb count guide output dirs, one per lane"},
            "guide_fastqs": {**_S_ARRAY_S, "description": "per lane 'R1,R2' guide FASTQs (counted here instead of guide_dirs)"},
            "guide_names": {**_S_ARRAY_S, "description": "lane names for guide_fastqs, e.g. S1_L1"},
            "guide_table": {**_S_STRING, "description": "guide table (xlsx/tsv) when counting FASTQs"},
            "chemistry": {**_S_STRING, "default": "10XV3", "description": "assay (10XV2/10XV3/5PE) or kallisto bc:umi:seq string"},
            "whitelist": {**_S_STRING, "description": "barcode whitelist"},
            "expected_cell_number": {**_S_INTEGER, "default": 8000, "description": "EXPECTED_CELL_NUMBER: cells kept per lane via knee[expected]"},
            "mito_expected_percentage": {**_S_NUMBER, "default": 0.2, "description": "maximum mitochondrial fraction"},
            "transcripts_umi_threshold": {**_S_INTEGER, "default": 100, "description": "TRANSCRIPTS_UMI_TRHESHOLD (applied as minimum genes per cell, as upstream)"},
            "mito_prefix": {**_S_STRING, "default": "MT-", "description": "mitochondrial gene-name prefix"},
            "percentage_of_cells_to_include_transcript": {**_S_NUMBER, "default": 0.01, "description": "genes kept if detected in >= this fraction of cells"},
            "guide_umi_limit": {**_S_INTEGER, "default": 5, "description": "guide called when UMIs > this"},
            "merge": {**_S_BOOLEAN, "default": False, "description": "sum guides of the same target before testing"},
            "multiseq_r1": {**_S_STRING, "description": "MULTI-seq R1 FASTQ"},
            "multiseq_r2": {**_S_STRING, "description": "MULTI-seq R2 FASTQ"},
            "multiseq_barcodes": {**_S_STRING, "description": "MULTI-seq barcode list"},
            "distance_neighbors": {**_S_INTEGER, "default": 1000000, "description": "cis window (bp) around the element TSS"},
            "in_trans": {**_S_STRING, "default": "FALSE", "description": "TRUE tests every gene"},
            "add_gene_names": {**_S_STRING, "default": "", "description": "comma-separated genes always tested"},
            "engine": {**_S_STRING, "default": "sceptre", "description": "sceptre | mannwhitney | sceptre-r"},
            "direction": {**_S_STRING, "default": "both", "description": "left | right | both"},
            "B": {**_S_INTEGER, "default": 500, "description": "conditional resamples per guide"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "reproduce the upstream bugs (see module docstring)"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "perturbseq", "description": "run label"}},
         "required": ["gtf", "rna_dirs"]},
        cli=["bulk-crispr", "run"],
        flag_map={"gtf": "--gtf", "rna_dirs": "--rna-dirs", "guide_dirs": "--guide-dirs", "guide_fastqs": "--guide-fastqs", "guide_names": "--guide-names", "guide_table": "--guide-table", "chemistry": "--chemistry", "whitelist": "--whitelist", "expected_cell_number": "--expected-cell-number", "mito_expected_percentage": "--mito-expected-percentage", "transcripts_umi_threshold": "--transcripts-umi-threshold", "mito_prefix": "--mito-prefix", "percentage_of_cells_to_include_transcript": "--percentage-of-cells-to-include-transcript", "guide_umi_limit": "--guide-umi-limit", "merge": "--merge", "multiseq_r1": "--multiseq-r1", "multiseq_r2": "--multiseq-r2", "multiseq_barcodes": "--multiseq-barcodes", "distance_neighbors": "--distance-neighbors", "in_trans": "--in-trans", "add_gene_names": "--add-gene-names", "engine": "--engine", "direction": "--direction", "B": "--B", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'guide_names', 'guide_dirs', 'rna_dirs', 'guide_fastqs'},
        bool_flags={'upstream_compat', 'no_plots', 'merge'},
    ),

    _T(
        "bulk_crispr_count_guides",
        "Counts guide UMIs per cell barcode from guide-capture FASTQs the way kallisto|bustools kite does: exact and Hamming-1 guide matching (colliding variants dropped), barcode/UMI/sequence positions from the kallisto technology string, whitelist correction at Hamming 1, UMIs counted once per barcode x guide. Inputs: guide_features.txt or guide table, R1/R2. Outputs <name>_ks_guide_out/counts_unfiltered/adata.h5ad and kite_stats.json.",
        {"type": "object", "properties": {
            "guides": {**_S_STRING, "description": "guide_features.txt or guide table"},
            "r1": {**_S_STRING, "description": "R1 FASTQ"},
            "r2": {**_S_STRING, "description": "R2 FASTQ"},
            "chemistry": {**_S_STRING, "default": "10XV3", "description": "assay or kallisto technology string"},
            "whitelist": {**_S_STRING, "description": "barcode whitelist"},
            "name": {**_S_STRING, "default": "S1_L1", "description": "<sample>_L<lane>"},
            "max_reads": {**_S_INTEGER, "default": 0, "description": "stop after N reads (0 = all)"},
            "label": {**_S_STRING, "description": "run label"}},
         "required": ["guides", "r1", "r2"]},
        cli=["bulk-crispr", "count-guides"],
        flag_map={"guides": "--guides", "r1": "--r1", "r2": "--r2", "chemistry": "--chemistry", "whitelist": "--whitelist", "name": "--name", "max_reads": "--max-reads", "label": "--label"},
    ),

    _T(
        "bulk_crispr_composition",
        "Computes the per-position compositional bias (sd between the A/C/T/G counts) of the first 10,000 FASTQ lines of R1 and R2, as the pipeline's read-composition QC. Inputs: R1 (and R2) FASTQ. Outputs a per-position TSV and a line plot per read.",
        {"type": "object", "properties": {
            "r1": {**_S_STRING, "description": "R1 FASTQ"},
            "r2": {**_S_STRING, "description": "R2 FASTQ"},
            "n_lines": {**_S_INTEGER, "default": 10000, "description": "FASTQ lines read"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "upstream line filter"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "composition", "description": "run label"}},
         "required": ["r1"]},
        cli=["bulk-crispr", "composition"],
        flag_map={"r1": "--r1", "r2": "--r2", "n_lines": "--n-lines", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'upstream_compat', 'no_plots'},
    ),

    _T(
        "bulk_crispr_guide_table",
        "Converts a guide design sheet (sgRNA_ID, Target_name, sgRNA_sequences, chr, start, end; xlsx/tsv/csv) into the pipeline's guide ids <target>|<n>_sgrna_<chr>:<start>:<end> and guide_features.txt for counting. Outputs guide_features.txt and the processed table.",
        {"type": "object", "properties": {
            "guides": {**_S_STRING, "description": "guide sheet"},
            "out_dir": {**_S_STRING, "description": "output directory"},
            "label": {**_S_STRING, "default": "guide_table", "description": "run label"}},
         "required": ["guides"]},
        cli=["bulk-crispr", "guide-table"],
        flag_map={"guides": "--guides", "out_dir": "--out-dir", "label": "--label"},
    ),

    _T(
        "bulk_crispr_prefilter",
        "Per-lane single-cell QC: knee plot, minimum genes (the UMI threshold, as upstream), knee[EXPECTED_CELL_NUMBER] count floor, mitochondrial fraction filter, Scrublet doublet removal, guide/scRNA barcode intersection, lane concatenation and the low-expression gene filter. Inputs: a lane manifest or kb count dirs. Outputs filtered scRNA and guide h5ad files plus lane_qc.tsv.",
        {"type": "object", "properties": {
            "manifest": {**_S_STRING, "description": "initial_preprocessing_file_names.txt"},
            "rna_dirs": {**_S_ARRAY_S, "description": "kb count cDNA dirs"},
            "guide_dirs": {**_S_ARRAY_S, "description": "kb count guide dirs"},
            "expected_cell_number": {**_S_INTEGER, "default": 8000, "description": "EXPECTED_CELL_NUMBER: cells kept per lane via knee[expected]"},
            "mito_expected_percentage": {**_S_NUMBER, "default": 0.2, "description": "maximum mitochondrial fraction"},
            "transcripts_umi_threshold": {**_S_INTEGER, "default": 100, "description": "TRANSCRIPTS_UMI_TRHESHOLD (applied as minimum genes per cell, as upstream)"},
            "mito_prefix": {**_S_STRING, "default": "MT-", "description": "mitochondrial gene-name prefix"},
            "percentage_of_cells_to_include_transcript": {**_S_NUMBER, "default": 0.01, "description": "genes kept if detected in >= this fraction of cells"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "ignore the gene-filter fraction as upstream"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "prefilter", "description": "run label"}}},
        cli=["bulk-crispr", "prefilter"],
        flag_map={"manifest": "--manifest", "rna_dirs": "--rna-dirs", "guide_dirs": "--guide-dirs", "expected_cell_number": "--expected-cell-number", "mito_expected_percentage": "--mito-expected-percentage", "transcripts_umi_threshold": "--transcripts-umi-threshold", "mito_prefix": "--mito-prefix", "percentage_of_cells_to_include_transcript": "--percentage-of-cells-to-include-transcript", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'guide_dirs', 'rna_dirs'},
        bool_flags={'upstream_compat', 'no_plots'},
    ),

    _T(
        "bulk_crispr_multiseq",
        "Demultiplexes MULTI-seq sample hashes with a port of deMULTIplex: tag-to-barcode alignment at Hamming <= 1, UMI counting per cell x barcode, KDE-maxima thresholds with the quantile sweep and two rounds of negative removal. Inputs: a modality bundle, MULTI-seq R1/R2 and the barcode list. Outputs bar_table.csv, final_class.csv and a bundle without doublets and negatives.",
        {"type": "object", "properties": {
            "muon_data": {**_S_STRING, "description": "modality bundle dir or .h5mu"},
            "r1": {**_S_STRING, "description": "MULTI-seq R1"},
            "r2": {**_S_STRING, "description": "MULTI-seq R2"},
            "barcodes": {**_S_STRING, "description": "barcode list"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "classify nUMI_total and keep doublets as upstream"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "multiseq", "description": "run label"}},
         "required": ["muon_data", "r1", "r2", "barcodes"]},
        cli=["bulk-crispr", "multiseq"],
        flag_map={"muon_data": "--muon-data", "r1": "--r1", "r2": "--r2", "barcodes": "--barcodes", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'upstream_compat', 'no_plots'},
    ),

    _T(
        "bulk_crispr_de",
        "Tests every guide of each element against its cis genes with a SCEPTRE-style conditional-resampling NB score test (or Mann-Whitney, or the upstream R sceptre script). Guides need more than 30 cells. Inputs: the perturbdata directory from perturb-loader and the bundle. Outputs per-element results.txt files and de_results.tsv with p_value, z_value, log_fold_change.",
        {"type": "object", "properties": {
            "perturbdata": {**_S_STRING, "description": "perturbdata/ directory"},
            "muon_data": {**_S_STRING, "description": "bundle (default: the one perturb-loader used)"},
            "engine": {**_S_STRING, "default": "sceptre", "description": "sceptre | mannwhitney | sceptre-r"},
            "direction": {**_S_STRING, "default": "both", "description": "left | right | both"},
            "B": {**_S_INTEGER, "default": 500, "description": "resamples"},
            "min_cells_per_guide": {**_S_INTEGER, "default": 30, "description": "guides need > this many cells"},
            "label": {**_S_STRING, "default": "de", "description": "run label"}},
         "required": ["perturbdata"]},
        cli=["bulk-crispr", "de"],
        flag_map={"perturbdata": "--perturbdata", "muon_data": "--muon-data", "engine": "--engine", "direction": "--direction", "B": "--B", "min_cells_per_guide": "--min-cells-per-guide", "label": "--label"},
    ),

    _T(
        "bulk_crispr_results",
        "Aggregates differential-perturbation results: BH adj_pvalue over all tests, guide x gene layers (p, adj p, z, log fold change, significant at adj < 0.01) and element-level Fisher-combined p per gene (sig_not_adj at 0.05). Inputs: sceptre_out directory and bundle. Outputs result_guides.tsv, result_elements.tsv and a 4-modality results bundle.",
        {"type": "object", "properties": {
            "sceptre_dir": {**_S_STRING, "description": "directory with results.txt files"},
            "muon_data": {**_S_STRING, "description": "bundle"},
            "alpha_guide": {**_S_NUMBER, "default": 0.01, "description": "guide-level BH threshold"},
            "alpha_element": {**_S_NUMBER, "default": 0.05, "description": "element-level threshold"},
            "label": {**_S_STRING, "default": "results", "description": "run label"}},
         "required": ["sceptre_dir", "muon_data"]},
        cli=["bulk-crispr", "results"],
        flag_map={"sceptre_dir": "--sceptre-dir", "muon_data": "--muon-data", "alpha_guide": "--alpha-guide", "alpha_element": "--alpha-element", "label": "--label"},
    ),

    _T(
        "bulk_crispr_tracks",
        "Writes genome-browser tracks of the tested guides: a BED of guide intervals scored by -log10 p and a pyGenomeTracks links file joining guides to their gene TSS, plus tracks.ini and an arc plot. Inputs: result_guides.tsv and the bundle. Outputs tracks_dir/.",
        {"type": "object", "properties": {
            "results": {**_S_STRING, "description": "result_guides.tsv"},
            "muon_data": {**_S_STRING, "description": "bundle"},
            "alpha": {**_S_NUMBER, "default": 0.05, "description": "adj p cut-off"},
            "all_pairs": {**_S_BOOLEAN, "default": False, "description": "write every test"},
            "region": {**_S_STRING, "description": "chr:start-end to render with pyGenomeTracks"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "tracks", "description": "run label"}},
         "required": ["results", "muon_data"]},
        cli=["bulk-crispr", "tracks"],
        flag_map={"results": "--results", "muon_data": "--muon-data", "alpha": "--alpha", "all_pairs": "--all-pairs", "region": "--region", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'all_pairs', 'no_plots'},
    ),

    _T(
        "bulk_crispr_assign_guides",
        "Calls guides per cell with a depth-aware two-component Poisson mixture next to the fixed UMI threshold, storing both as layers. Inputs: a modality bundle with guide counts. Outputs guide_assignment.tsv (rates, agreement) and a bundle with a 'binarized' layer.",
        {"type": "object", "properties": {
            "muon_data": {**_S_STRING, "description": "bundle"},
            "guide_umi_limit": {**_S_INTEGER, "default": 5, "description": "threshold rule and mixture initialisation"},
            "label": {**_S_STRING, "default": "assign_guides", "description": "run label"}},
         "required": ["muon_data"]},
        cli=["bulk-crispr", "assign-guides"],
        flag_map={"muon_data": "--muon-data", "guide_umi_limit": "--guide-umi-limit", "label": "--label"},
    ),

    _T(
        "bulk_crispr_perturb_loader",
        "Builds the element x guide membership and element x tested-gene sets: genes within DISTANCE_NEIGHBORS of the element TSS (or of the guide for putative enhancers), all genes with in_trans TRUE, 10 random genes for controls, plus ADDGENENAMES; GUIDE_TYPE from the element name. Inputs: bundle and GTF. Outputs perturbdata/ TSVs.",
        {"type": "object", "properties": {
            "muon_data": {**_S_STRING, "description": "bundle"},
            "gtf": {**_S_STRING, "description": "GTF"},
            "distance_from_guide": {**_S_INTEGER, "default": 1000000, "description": "cis window bp"},
            "in_trans": {**_S_STRING, "default": "FALSE", "description": "TRUE or FALSE"},
            "add_gene_names": {**_S_STRING, "default": "", "description": "comma-separated genes"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "random genes for enhancers as upstream"},
            "label": {**_S_STRING, "default": "perturb_loader", "description": "run label"}},
         "required": ["muon_data", "gtf"]},
        cli=["bulk-crispr", "perturb-loader"],
        flag_map={"muon_data": "--muon-data", "gtf": "--gtf", "distance_from_guide": "--distance-from-guide", "in_trans": "--in-trans", "add_gene_names": "--add-gene-names", "upstream_compat": "--upstream-compat", "label": "--label"},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "bulk_crispr_config",
        "Parses and validates a Nextflow perturb.config (params.X = ...): required parameters, FASTQ/name list lengths, <sample>_L<lane> names, chemistry resolution, MULTI-seq parameters. Input: the config file. Output: JSON of params, resolved defaults and issues.",
        {"type": "object", "properties": {
            "config": {**_S_STRING, "description": "perturb.config"},
            "out": {**_S_STRING, "description": "write the JSON here"}},
         "required": ["config"]},
        cli=["bulk-crispr", "config"],
        flag_map={"config": "--config", "out": "--out"},
    ),

    _T(
        "bulk_crispr_map_rna",
        "Builds and runs the kallisto|bustools commands for the cDNA library (kb ref -d <reference>, kb count with the resolved chemistry and whitelist) when kb is on PATH; otherwise prints the exact commands. Inputs: FASTQs, chemistry, whitelist. Output: <name>_ks_transcripts_out.",
        {"type": "object", "properties": {
            "fastq": {**_S_ARRAY_S, "description": "R1 R2 FASTQs"},
            "chemistry": {**_S_STRING, "default": "10XV3", "description": "assay or technology string"},
            "whitelist": {**_S_STRING, "description": "barcode whitelist"},
            "name": {**_S_STRING, "default": "S1_L1", "description": "<sample>_L<lane>"},
            "reference": {**_S_STRING, "default": "human", "description": "kb ref -d reference"},
            "dry_run": {**_S_BOOLEAN, "default": False, "description": "print only"},
            "label": {**_S_STRING, "description": "run label"}},
         "required": ["fastq"]},
        cli=["bulk-crispr", "map-rna"],
        flag_map={"fastq": "--fastq", "chemistry": "--chemistry", "whitelist": "--whitelist", "name": "--name", "reference": "--reference", "dry_run": "--dry-run", "label": "--label"},
        flag_repeat={'fastq'},
        bool_flags={'dry_run'},
    ),

    _T(
        "bulk_crispr_assay_spec",
        "Returns the kallisto technology string and 10x whitelist URL for an assay name (10XV2, 10XV3, 5PE) or validates a custom bc:umi:seq string. Input: assay name. Output: 'chemistry,whitelist' and a JSON layout.",
        {"type": "object", "properties": {
            "assay": {**_S_STRING, "description": "10XV2, 10XV3, 5PE, custom, or a technology string"},
            "chemistry": {**_S_STRING, "description": "custom technology string"},
            "whitelist": {**_S_STRING, "description": "custom whitelist"}},
         "required": ["assay"]},
        cli=["bulk-crispr", "assay-spec"],
        flag_map={"assay": "--assay", "chemistry": "--chemistry", "whitelist": "--whitelist"},
    ),

    _T(
        "bulk_crispr_cellranger_inputs",
        "Writes Cell Ranger feature-barcode inputs from the guide sheet: feature_ref.csv (CRISPR Guide Capture) and library.csv, and prints (or with execute runs) cellranger count. Inputs: guide sheet, FASTQ dirs, sample names, transcriptome. Outputs the two CSVs.",
        {"type": "object", "properties": {
            "guides": {**_S_STRING, "description": "guide sheet"},
            "rna_fastq_dir": {**_S_STRING, "description": "cDNA FASTQ dir"},
            "guide_fastq_dir": {**_S_STRING, "description": "guide FASTQ dir"},
            "rna_sample": {**_S_STRING, "description": "cDNA sample name"},
            "guide_sample": {**_S_STRING, "description": "guide sample name"},
            "transcriptome": {**_S_STRING, "description": "Cell Ranger reference"},
            "pattern": {**_S_STRING, "default": "(BC)", "description": "feature pattern"},
            "execute": {**_S_BOOLEAN, "default": False, "description": "run cellranger when installed"},
            "label": {**_S_STRING, "default": "cellranger_inputs", "description": "run label"}},
         "required": ["guides", "rna_fastq_dir", "guide_fastq_dir", "rna_sample", "guide_sample", "transcriptome"]},
        cli=["bulk-crispr", "cellranger-inputs"],
        flag_map={"guides": "--guides", "rna_fastq_dir": "--rna-fastq-dir", "guide_fastq_dir": "--guide-fastq-dir", "rna_sample": "--rna-sample", "guide_sample": "--guide-sample", "transcriptome": "--transcriptome", "pattern": "--pattern", "execute": "--execute", "label": "--label"},
        bool_flags={'execute'},
    ),

    _T(
        "bulk_crispr_inspect",
        "Summarises a Perturb-seq modality bundle or .h5mu: modality shapes, obs/var columns, layers, guides per target element and guide types; with target, the cells carrying that element's guides and mean expression of chosen genes with versus without. Input: bundle. Output: JSON on stdout.",
        {"type": "object", "properties": {
            "muon_data": {**_S_STRING, "description": "bundle or .h5mu"},
            "target": {**_S_STRING, "description": "target element"},
            "genes": {**_S_ARRAY_S, "description": "genes to compare"}},
         "required": ["muon_data"]},
        cli=["bulk-crispr", "inspect"],
        flag_map={"muon_data": "--muon-data", "target": "--target", "genes": "--genes"},
        flag_repeat={'genes'},
    ),

    _T(
        "crisprdevtools_new_module",
        "Scaffolds a new Nextflow DSL2 module directory in the IGVF CRISPR pipeline layout, as crisprdevtools.create_new_module_nextflow does: bin/, conda_envs/, example_data/, processes/, test/, README.md, input.config, main.nf, bin/<name>.py, conda_envs/<name>.yaml and processes/<name>.nf. The 'full' template adds an include plus workflow, a params block, test/test.nf, example data and an argparse script. Outputs the module directory and a report with the layout check.",
        {"type": "object", "properties": {
            "name": {**_S_STRING, "description": "module name (Nextflow identifier)"},
            "dest": {**_S_STRING, "description": "parent directory for the module (default: the run directory)"},
            "template": {**_S_STRING, "default": "minimal", "description": "minimal (upstream contents) or full"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "reproduce upstream exactly: indented (invalid) YAML, process named seqSpecParser, no chmod"},
            "force": {**_S_BOOLEAN, "default": False, "description": "overwrite an existing module"},
            "label": {**_S_STRING, "description": "run label"}},
         "required": ["name"]},
        cli=["crisprdevtools", "new-module"],
        flag_map={"name": "--name", "dest": "--dest", "template": "--template", "upstream_compat": "--upstream-compat", "force": "--force", "label": "--label"},
        bool_flags={'upstream_compat', 'force'},
    ),

    _T(
        "crisprdevtools_check_module",
        "Validates an existing Nextflow module directory against the IGVF CRISPR module layout. It checks the required dirs and the DSL2 header, that the process declarations and includes resolve, that conda env files exist and parse, that bin/ scripts have a shebang and are executable, and that the params main.nf uses are set in input.config. Input: module path. Outputs a PASS/WARN/FAIL table in report.md and summary.json; exit 1 on any FAIL.",
        {"type": "object", "properties": {
            "path": {**_S_STRING, "description": "module directory"},
            "label": {**_S_STRING, "description": "run label"}},
         "required": ["path"]},
        cli=["crisprdevtools", "check-module"],
        flag_map={"path": "--path", "label": "--label"},
    ),

    _T(
        "crispr_fg_assay_spec",
        "Resolves a single-cell CRISPR assay name (10xv2, 10xv3, 10x5p, 10x5p-pe, dropseq, celseq2, or custom) to the kallisto bus -x read-format string and barcode whitelist, as Task 1 of the CRISPR-FG jamboree specifies. Prints the 'chemistry,whitelist' line the pipeline consumes and writes report.md + summary.json; unknown assays exit 2 (the pipeline stops). Optionally downloads the 10x whitelist.",
        {"type": "object", "properties": {
            "assay": {**_S_STRING, "description": "Assay name or 'custom'"},
            "chemistry": {**_S_STRING, "description": "kallisto -x string for custom"},
            "whitelist": {**_S_STRING, "description": "Whitelist path for custom"},
            "download": {**_S_BOOLEAN, "default": False, "description": "Download the whitelist"},
            "label": {**_S_STRING, "default": "assay_spec", "description": "Run label"}},
         "required": ["assay"]},
        cli=["crispr-fg-jamboree", "assay-spec"],
        flag_map={"assay": "--assay", "chemistry": "--chemistry", "whitelist": "--whitelist", "download": "--download", "label": "--label"},
        bool_flags={'download'},
    ),

    _T(
        "crispr_fg_cellranger_inputs",
        "Builds Cell Ranger feature-barcode inputs from a guide table (xlsx/csv/tsv with Target_name, chr, start, end, sgRNA_sequences): guides are renamed Target|n and given pipeline ids '<Target|n>_sgrna_<chr>:<start>:<end>'. Writes feature_ref.csv, library.csv (from the RNA and guide FASTQ paths), guide_features.txt and the cellranger count command (run only with run=true and cellranger on PATH).",
        {"type": "object", "properties": {
            "guide_table": {**_S_STRING, "description": "Guide table path"},
            "rna_fastqs": {**_S_ARRAY_S, "description": "RNA FASTQ paths"},
            "guide_fastqs": {**_S_ARRAY_S, "description": "Guide FASTQ paths"},
            "transcriptome": {**_S_STRING, "default": "refdata-gex-GRCh38-2020-A", "description": "Cell Ranger reference dir"},
            "read": {**_S_STRING, "default": "R2", "description": "feature_ref read"},
            "pattern": {**_S_STRING, "default": "(BC)", "description": "feature_ref pattern"},
            "run": {**_S_BOOLEAN, "default": False, "description": "Run cellranger if installed"},
            "label": {**_S_STRING, "default": "cellranger_inputs", "description": "Run label"}},
         "required": ["guide_table"]},
        cli=["crispr-fg-jamboree", "cellranger-inputs"],
        flag_map={"guide_table": "--guide-table", "rna_fastqs": "--rna-fastqs", "guide_fastqs": "--guide-fastqs", "transcriptome": "--transcriptome", "read": "--read", "pattern": "--pattern", "run": "--run", "label": "--label"},
        flag_repeat={'guide_fastqs', 'rna_fastqs'},
        bool_flags={'run'},
    ),

    _T(
        "crispr_fg_pipeline_config",
        "Writes the perturb.config of the LucasSilvaFerreira/pipeline_perturbseq_like Nextflow pipeline with the jamboree defaults (GUIDE_UMI_LIMIT 5, DISTANCE_NEIGHBORS 1 Mb, MITO_EXPECTED_PERCENTAGE 0.2, EXPECTED_CELL_NUMBER 10000, ...) plus KEY=VALUE overrides, and prints the nextflow launch command. Tower tokens are only read from TOWER_ACCESS_TOKEN, never written; launches only with run=true and nextflow on PATH.",
        {"type": "object", "properties": {
            "guide_features": {**_S_STRING, "description": "Guide table path for GUIDE_FEATURES"},
            "set": {**_S_ARRAY_S, "description": "KEY=VALUE overrides"},
            "run": {**_S_BOOLEAN, "default": False, "description": "Launch nextflow if installed"},
            "label": {**_S_STRING, "default": "pipeline_config", "description": "Run label"}}},
        cli=["crispr-fg-jamboree", "pipeline-config"],
        flag_map={"guide_features": "--guide-features", "set": "--set", "run": "--run", "label": "--label"},
        flag_repeat={'set'},
        bool_flags={'run'},
    ),

    _T(
        "crispr_fg_preprocess",
        "Runs the pipeline's cell QC on per-lane RNA and guide count matrices (.h5ad or 10x matrix directories): n_genes >= 100, knee cut at the expected cell number, percent_mito < 0.2, optional Scrublet doublets, shared barcodes, genes in >= 1% of cells. Adds guide coordinates parsed from the pipeline ids, gene TSS from a GTF/TSV and the sceptre covariates, and writes raw_mudata_guide_and_transcripts.h5mu (guides + scRNA) with QC tables and figures.",
        {"type": "object", "properties": {
            "rna": {**_S_ARRAY_S, "description": "Per-lane RNA counts"},
            "guides": {**_S_ARRAY_S, "description": "Per-lane guide counts"},
            "gene_table": {**_S_STRING, "description": "GTF or gene TSV for TSS coordinates"},
            "expected_cells": {**_S_INTEGER, "default": 10000, "description": "EXPECTED_CELL_NUMBER"},
            "mito_max": {**_S_NUMBER, "default": 0.2, "description": "Max mito fraction"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "Reproduce upstream covariate definitions"},
            "label": {**_S_STRING, "default": "preprocess", "description": "Run label"}},
         "required": ["rna", "guides"]},
        cli=["crispr-fg-jamboree", "preprocess"],
        flag_map={"rna": "--rna", "guides": "--guides", "gene_table": "--gene-table", "expected_cells": "--expected-cells", "mito_max": "--mito-max", "upstream_compat": "--upstream-compat", "label": "--label"},
        flag_repeat={'rna', 'guides'},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "crispr_fg_mudata_qc",
        "Computes the Task 2 QC of a guides/scRNA MuData: per-cell genes/UMIs/mito, highest-expressed genes, the n_genes < 6500 & mito < 0.1 filter, MOI histogram, guide and element coverage, positive/negative control classification and (with scanpy) highly variable genes. Writes TSVs, a QC figure, report.md and summary.json.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "Input .h5mu"},
            "max_genes": {**_S_INTEGER, "default": 6500, "description": "n_genes upper bound"},
            "max_mito": {**_S_NUMBER, "default": 0.1, "description": "mito fraction upper bound"},
            "label": {**_S_STRING, "default": "mudata_qc", "description": "Run label"}},
         "required": ["mudata"]},
        cli=["crispr-fg-jamboree", "mudata-qc"],
        flag_map={"mudata": "--mudata", "max_genes": "--max-genes", "max_mito": "--max-mito", "label": "--label"},
    ),

    _T(
        "crispr_fg_guide_gene_distance",
        "Computes guide x gene distances as sccrispr-tools does (|guide start - gene TSS| for pairs on the same chromosome) and flags cis pairs within the distance (default 1 Mb). Input is a guides/scRNA MuData with guide ids and transcript coordinates; writes guide_gene_distance.tsv.gz and cis_pairs.tsv.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "Input .h5mu"},
            "distance": {**_S_INTEGER, "default": 1000000, "description": "cis window in bp"},
            "how": {**_S_STRING, "default": "start", "description": "start or midpoint"},
            "label": {**_S_STRING, "default": "guide_gene_distance", "description": "Run label"}},
         "required": ["mudata"]},
        cli=["crispr-fg-jamboree", "guide-gene-distance"],
        flag_map={"mudata": "--mudata", "distance": "--distance", "how": "--how", "label": "--label"},
    ),

    _T(
        "crispr_fg_subset",
        "Subsets a guides/scRNA MuData to guides of chosen elements, control categories (POSITIVE_CONTROL, NEGATIVE_CONTROL, PUTATIVE_ENHANCER) or ids, and/or keeps only some modalities (e.g. drops result_* modalities). For a single guide it also exports its presence per cell and the counts of its cis genes. Writes subset.h5mu.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "Input .h5mu"},
            "element": {**_S_ARRAY_S, "description": "Target elements"},
            "category": {**_S_ARRAY_S, "description": "Guide categories"},
            "guide": {**_S_ARRAY_S, "description": "Guide ids"},
            "keep_modalities": {**_S_ARRAY_S, "description": "Modalities to keep"},
            "label": {**_S_STRING, "default": "subset", "description": "Run label"}},
         "required": ["mudata"]},
        cli=["crispr-fg-jamboree", "subset"],
        flag_map={"mudata": "--mudata", "element": "--element", "category": "--category", "guide": "--guide", "keep_modalities": "--keep-modalities", "label": "--label"},
        flag_repeat={'category', 'guide', 'element', 'keep_modalities'},
    ),

    _T(
        "crispr_fg_assign_guides",
        "Calls guides per cell and stores the binary matrix as the 'binarized' layer of the guides modality: method umi uses the pipeline rule UMI > GUIDE_UMI_LIMIT (5); poisson-mixture fits a per-guide Poisson GLM on the cell covariates and a two-component mixture (posterior >= 0.8). Optional merging of guides per element. Writes mu_with_binary.h5mu and a per-guide summary.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "Input .h5mu"},
            "method": {**_S_STRING, "default": "umi", "description": "umi or poisson-mixture"},
            "guide_umi_limit": {**_S_INTEGER, "default": 5, "description": "UMI limit"},
            "merge": {**_S_BOOLEAN, "default": False, "description": "Sum guides per element first"},
            "label": {**_S_STRING, "default": "assign_guides", "description": "Run label"}},
         "required": ["mudata"]},
        cli=["crispr-fg-jamboree", "assign-guides"],
        flag_map={"mudata": "--mudata", "method": "--method", "guide_umi_limit": "--guide-umi-limit", "merge": "--merge", "label": "--label"},
        bool_flags={'merge'},
    ),

    _T(
        "crispr_fg_differential",
        "Tests guide-gene pairs for differential expression (Task 3): per element, genes within DISTANCE_NEIGHBORS of the element gene's TSS (or all genes in trans), guides with > 30 assigned cells. The default sceptre-nb test is a SCEPTRE-style NB GLM score test with conditional resampling; mannwhitney and welch-t are alternative modules; engine r writes and runs the upstream run_sceptre_high_moi inputs. Writes differential_results.tsv.gz and mudata_results.h5mu (result_guides with adj_pvalue/log_fold_change/significant/z_value layers, result_elements Fisher-combined).",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "Input .h5mu with a binarized guide layer"},
            "distance": {**_S_INTEGER, "default": 1000000, "description": "cis window in bp"},
            "in_trans": {**_S_BOOLEAN, "default": False, "description": "Test all genes"},
            "test": {**_S_STRING, "default": "sceptre-nb", "description": "sceptre-nb, mannwhitney or welch-t"},
            "side": {**_S_STRING, "default": "both", "description": "both, left or right"},
            "engine": {**_S_STRING, "default": "python", "description": "python or r"},
            "min_cells_per_guide": {**_S_INTEGER, "default": 30, "description": "Minimum assigned cells (strictly greater)"},
            "gene_table": {**_S_STRING, "description": "GTF/TSV if scRNA.var lacks coordinates"},
            "label": {**_S_STRING, "default": "differential", "description": "Run label"}},
         "required": ["mudata"]},
        cli=["crispr-fg-jamboree", "differential"],
        flag_map={"mudata": "--mudata", "distance": "--distance", "in_trans": "--in-trans", "test": "--test", "side": "--side", "engine": "--engine", "min_cells_per_guide": "--min-cells-per-guide", "gene_table": "--gene-table", "label": "--label"},
        bool_flags={'in_trans'},
    ),

    _T(
        "crispr_fg_tracks",
        "Turns a mudata_results.h5mu into genome-browser tracks (Task 4): a BED of every guide-gene test (name guide|gene, score 100 x -log10 p), a bedGraph of the best p per guide, pyGenomeTracks links (BEDPE-like) from guides to gene TSSs for BH-significant pairs, a tracks.ini, per-element arc figures and the pyGenomeTracks commands (run when installed).",
        {"type": "object", "properties": {
            "results": {**_S_STRING, "description": "mudata_results.h5mu"},
            "all_links": {**_S_BOOLEAN, "default": False, "description": "Links for every test"},
            "fdr": {**_S_NUMBER, "default": 0.01, "description": "BH cut-off for links"},
            "label": {**_S_STRING, "default": "tracks", "description": "Run label"}},
         "required": ["results"]},
        cli=["crispr-fg-jamboree", "tracks"],
        flag_map={"results": "--results", "all_links": "--all-links", "fdr": "--fdr", "label": "--label"},
        bool_flags={'all_links'},
    ),

    _T(
        "crispr_fg_run",
        "Runs assign-guides -> differential -> tracks on a raw guides/scRNA MuData in one call, with the jamboree defaults (UMI > 5, 1 Mb cis window, SCEPTRE-style NB test). Produces the three run directories with their reports, the results MuData and the track files.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "Raw .h5mu"},
            "method": {**_S_STRING, "default": "umi", "description": "umi or poisson-mixture"},
            "distance": {**_S_INTEGER, "default": 1000000, "description": "cis window in bp"},
            "label": {**_S_STRING, "default": "jamboree", "description": "Run label"}},
         "required": ["mudata"]},
        cli=["crispr-fg-jamboree", "run"],
        flag_map={"mudata": "--mudata", "method": "--method", "distance": "--distance", "label": "--label"},
    ),

    _T(
        "fishash_build_inputs",
        "Builds the per-sample inputs of the Fishash Table 2 reproduction from raw GSE272457-style 10x barnyard triplets (<prefix>_features.tsv.gz, _barcodes.tsv.gz, _matrix.mtx.gz). Writes the guide x cell Matrix Market file (CLEANSER input), a guide .h5ad, <sample>_meta.csv (human/mouse/mito UMI sums, detected genes) and <sample>_guides.csv (nt_k <= 100 -> homo_guide).",
        {"type": "object", "properties": {
            "raw_dir": {**_S_STRING, "description": "Directory with the GSE272457 files"},
            "raw_prefix": {**_S_STRING, "description": "One raw triplet prefix"},
            "sample": {**_S_STRING, "description": "Sample name for raw_prefix"},
            "samples": {**_S_ARRAY_S, "description": "Samples to build with raw_dir"},
            "work_dir": {**_S_STRING, "description": "Output work directory"}},
         "required": ["work_dir"]},
        cli=["fishash-table2", "build-inputs"],
        flag_map={"raw_dir": "--raw-dir", "raw_prefix": "--raw-prefix", "sample": "--sample", "samples": "--samples", "work_dir": "--work-dir"},
        flag_repeat={'samples'},
    ),

    _T(
        "fishash_cleanser",
        "Computes CLEANSER 1.2.1 guide-assignment posteriors in Python: the zero-truncated Poisson/NB (CROP-seq, cs) or NB/NB (direct capture, dc) mixture per guide with CLEANSER's priors, bounds and low-pass library-size normalisation, sampled by adaptive Metropolis (or MAP) and summarised as the median posterior PZi. Input is the guide x cell .mtx; output is CLEANSER's coordinate posterior file plus a log carrying the 'Random seed:' marker. engine=cleanser runs the real binary when installed.",
        {"type": "object", "properties": {
            "input": {**_S_STRING, "description": "Guide x cell Matrix Market"},
            "mode": {**_S_STRING, "description": "cs or dc"},
            "out": {**_S_STRING, "description": "Posterior output path"},
            "method": {**_S_STRING, "default": "mcmc", "description": "mcmc or map"},
            "chains": {**_S_INTEGER, "default": 4, "description": "MCMC chains"},
            "num_samples": {**_S_INTEGER, "default": 1000, "description": "Draws per chain"},
            "num_warmup": {**_S_INTEGER, "default": 300, "description": "Warm-up iterations"},
            "seed": {**_S_INTEGER, "default": 20260810, "description": "Seed"},
            "engine": {**_S_STRING, "default": "python", "description": "python or cleanser"}},
         "required": ["input", "mode", "out"]},
        cli=["fishash-table2", "cleanser"],
        flag_map={"input": "--input", "mode": "--mode", "out": "--out", "method": "--method", "chains": "--chains", "num_samples": "--num-samples", "num_warmup": "--num-warmup", "seed": "--seed", "engine": "--engine"},
    ),

    _T(
        "fishash_sceptre_mixture",
        "Assigns guides with a Python approximation of SCEPTRE 0.10.3 assign_grnas(method='mixture'): per-gRNA Poisson GLM on the default high-MOI covariates, a reduced two-component EM with 5 restarts and posterior >= 0.8, with the count >= 5 backup for gRNAs with < 10 nonzero cells. engine=r writes (and with Rscript + sceptre runs) the upstream-equivalent R script instead. Writes <sample>_sceptre_mixture.mtx and per-guide fits.",
        {"type": "object", "properties": {
            "work_dir": {**_S_STRING, "default": ".", "description": "Work directory"},
            "sample": {**_S_STRING, "description": "Sample name"},
            "raw_prefix": {**_S_STRING, "description": "Raw 10x triplet prefix"},
            "raw_dir": {**_S_STRING, "description": "Directory with GSE272457 files"},
            "engine": {**_S_STRING, "default": "python", "description": "python or r"},
            "probability_threshold": {**_S_NUMBER, "default": 0.8, "description": "Posterior cut-off"}},
         "required": ["sample"]},
        cli=["fishash-table2", "sceptre-mixture"],
        flag_map={"work_dir": "--work-dir", "sample": "--sample", "raw_prefix": "--raw-prefix", "raw_dir": "--raw-dir", "engine": "--engine", "probability_threshold": "--probability-threshold"},
    ),

    _T(
        "fishash_score",
        "Scores species assignment exactly as score_table2.R: the Table 2 cohort (mito < 15%, 1500-6000 genes, 3500-20000 UMIs, species purity > 90%), prediction = human-only vs mouse-only guide calls with both/neither counted wrong, accuracy and its stderr for CLEANSER cs 0.80, dc 0.50, SCEPTRE mixture and the 0.95 public-code cut-offs. Compares each row with the published Table 2 value and writes extended_table2_reproduction.csv, a figure, report.md and summary.json.",
        {"type": "object", "properties": {
            "work_dir": {**_S_STRING, "description": "Work directory with inputs and caller outputs"},
            "samples": {**_S_ARRAY_S, "description": "Samples to score"},
            "strict": {**_S_BOOLEAN, "default": False, "description": "Stop on missing outputs"},
            "label": {**_S_STRING, "default": "table2", "description": "Run label"}},
         "required": ["work_dir"]},
        cli=["fishash-table2", "score"],
        flag_map={"work_dir": "--work-dir", "samples": "--samples", "strict": "--strict", "label": "--label"},
        flag_repeat={'samples'},
        bool_flags={'strict'},
    ),

    _T(
        "fishash_compare",
        "Compares a scored CSV (sample, method, accuracy) with the published Fishash Table 2 accuracies and with the upstream repository's own reproduced accuracies, flagging rows within a tolerance (default 0.001). Writes comparison.csv, report.md and summary.json.",
        {"type": "object", "properties": {
            "scored": {**_S_STRING, "description": "Scored CSV"},
            "tolerance": {**_S_NUMBER, "default": 0.001, "description": "Tolerance"},
            "label": {**_S_STRING, "default": "table2_compare", "description": "Run label"}},
         "required": ["scored"]},
        cli=["fishash-table2", "compare"],
        flag_map={"scored": "--scored", "tolerance": "--tolerance", "label": "--label"},
    ),

    _T(
        "fishash_run",
        "Runs the whole Fishash Table 2 reproduction on raw barnyard files: build-inputs, SCEPTRE-style mixture assignment, CLEANSER cs and dc posteriors, then scoring against Table 2. Outputs are cached in the work directory, like run_reproduction.sh, so reruns skip finished steps.",
        {"type": "object", "properties": {
            "raw_dir": {**_S_STRING, "description": "Directory with GSE272457 files"},
            "raw_prefix": {**_S_STRING, "description": "One raw triplet prefix"},
            "sample": {**_S_STRING, "description": "Sample name for raw_prefix"},
            "work_dir": {**_S_STRING, "description": "Work directory"},
            "cleanser_method": {**_S_STRING, "default": "mcmc", "description": "mcmc or map"},
            "chains": {**_S_INTEGER, "default": 4, "description": "CLEANSER chains"},
            "label": {**_S_STRING, "default": "table2", "description": "Run label"}},
         "required": ["work_dir"]},
        cli=["fishash-table2", "run"],
        flag_map={"raw_dir": "--raw-dir", "raw_prefix": "--raw-prefix", "sample": "--sample", "work_dir": "--work-dir", "cleanser_method": "--cleanser-method", "chains": "--chains", "label": "--label"},
    ),

    _T(
        "crispr_seqspec_catalog",
        "Catalogues every seqspec in IGVF-CRISPR/CRISPR-SeqSpec (CROP-seq + 10x v3 + MULTI-seq rna/guide/multiseq specs, TAP-seq Schraivogel 2020): assay, seqspec version, modalities, reads, regions, onlists, kb read-format string and seqspec-check error counts. Uses the embedded index of the pinned commit (or a fetched copy / a local directory). Writes catalog.tsv, reads.tsv, regions.tsv, onlists.tsv, catalog.json and report.md.",
        {"type": "object", "properties": {
            "live": {**_S_BOOLEAN, "default": False, "description": "Re-read the specs from the pinned upstream commit"},
            "from_dir": {**_S_STRING, "description": "Scan a local directory of seqspec YAMLs instead"},
            "label": {**_S_STRING, "default": "catalog", "description": "Run label for the output directory"}}},
        cli=["crispr-seqspec", "catalog"],
        flag_map={"live": "--live", "from_dir": "--from-dir", "label": "--label"},
        bool_flags={'live'},
    ),

    _T(
        "crispr_seqspec_check",
        "Validates a seqspec file the way seqspec check does: required keys and types, enums (modality, region_type, sequence_type, strand, onlist location), modality/library consistency, unique ids, primer ids, onlist file presence and md5, sequence length and alphabet by sequence_type, parent lengths vs children, read length vs available sequence, FASTQ presence. Handles seqspec 0.2.x and legacy 0.0.x. Writes issues.tsv, summary.json and report.md.",
        {"type": "object", "properties": {
            "spec": {**_S_STRING, "description": "seqspec YAML path, or a catalogue entry / alias: rna, guide, multiseq, tapseq (multiseq_10xv3_cropseq/*.yml, TAPseq_Schraivogel2020/spec.yaml)"},
            "spec_dir": {**_S_STRING, "description": "Directory holding onlists/FASTQs (default: the spec's directory)"},
            "no_md5": {**_S_BOOLEAN, "default": False, "description": "Skip md5 verification of local onlists"},
            "strict": {**_S_BOOLEAN, "default": False, "description": "Exit non-zero when errors are found"},
            "offline": {**_S_BOOLEAN, "default": False, "description": "Never fetch from GitHub; use the embedded index"},
            "label": {**_S_STRING, "default": "check", "description": "Run label for the output directory"}},
         "required": ["spec"]},
        cli=["crispr-seqspec", "check"],
        flag_map={"spec": "--spec", "spec_dir": "--spec-dir", "no_md5": "--no-md5", "strict": "--strict", "offline": "--offline", "label": "--label"},
        bool_flags={'strict', 'no_md5', 'offline'},
    ),

    _T(
        "crispr_seqspec_index",
        "Computes the read-format string of a seqspec modality for a tool, equivalent to seqspec index -t kb|chromap: projects each read onto the library (primer + strand, cut at read length) and reports barcode/UMI/feature coordinates as kb 'bc:umi:feature' triples, chromap --read-format, STARsolo CB/UMI arguments or a region table. Writes index_regions.tsv, summary.json and report.md and prints the string.",
        {"type": "object", "properties": {
            "spec": {**_S_STRING, "description": "seqspec YAML path, or a catalogue entry / alias: rna, guide, multiseq, tapseq (multiseq_10xv3_cropseq/*.yml, TAPseq_Schraivogel2020/spec.yaml)"},
            "modality": {**_S_STRING, "description": "Modality (default: all in the spec)"},
            "tool": {**_S_STRING, "default": "kb", "description": "kb, chromap, starsolo or tab"},
            "dedupe": {**_S_BOOLEAN, "default": False, "description": "Keep each region only in the first read covering it (reads longer than the insert)"},
            "feature_types": {**_S_ARRAY_S, "description": "Region types to report as the feature (default cdna/gdna, or sgRNA regions for crispr)"},
            "offline": {**_S_BOOLEAN, "default": False, "description": "Never fetch from GitHub; use the embedded index"},
            "label": {**_S_STRING, "default": "index", "description": "Run label for the output directory"}},
         "required": ["spec"]},
        cli=["crispr-seqspec", "index"],
        flag_map={"spec": "--spec", "modality": "--modality", "tool": "--tool", "dedupe": "--dedupe", "feature_types": "--feature-types", "offline": "--offline", "label": "--label"},
        flag_repeat={'feature_types'},
        bool_flags={'offline', 'dedupe'},
    ),

    _T(
        "crispr_seqspec_onlist",
        "Lists the onlist (barcode whitelist / guide list) filenames of a seqspec per modality and region, with location, md5 and the resolved local path when present (seqspec onlist equivalent). Writes onlists.tsv, summary.json and report.md.",
        {"type": "object", "properties": {
            "spec": {**_S_STRING, "description": "seqspec YAML path, or a catalogue entry / alias: rna, guide, multiseq, tapseq (multiseq_10xv3_cropseq/*.yml, TAPseq_Schraivogel2020/spec.yaml)"},
            "modality": {**_S_STRING, "description": "Restrict to one modality"},
            "offline": {**_S_BOOLEAN, "default": False, "description": "Never fetch from GitHub; use the embedded index"},
            "label": {**_S_STRING, "default": "onlist", "description": "Run label for the output directory"}},
         "required": ["spec"]},
        cli=["crispr-seqspec", "onlist"],
        flag_map={"spec": "--spec", "modality": "--modality", "offline": "--offline", "label": "--label"},
        bool_flags={'offline'},
    ),

    _T(
        "crispr_seqspec_info",
        "Prints the library structure of a seqspec: region tree with types, lengths, sequences and onlists, and each read's region layout in read coordinates. Stdout only.",
        {"type": "object", "properties": {
            "spec": {**_S_STRING, "description": "seqspec YAML path, or a catalogue entry / alias: rna, guide, multiseq, tapseq (multiseq_10xv3_cropseq/*.yml, TAPseq_Schraivogel2020/spec.yaml)"},
            "offline": {**_S_BOOLEAN, "default": False, "description": "Never fetch from GitHub; use the embedded index"}},
         "required": ["spec"]},
        cli=["crispr-seqspec", "info"],
        flag_map={"spec": "--spec", "offline": "--offline"},
        bool_flags={'offline'},
    ),

    _T(
        "crispr_seqspec_check_reads",
        "Checks FASTQs against a seqspec (the upstream submission rule of one read pair per modality): FASTQ presence, read lengths within [min_len, max_len], fixed regions matching at their projected coordinates (Hamming <= 1), and barcode/guide slices found in their onlists. Writes read_checks.tsv, summary.json and report.md.",
        {"type": "object", "properties": {
            "spec": {**_S_STRING, "description": "seqspec YAML path, or a catalogue entry / alias: rna, guide, multiseq, tapseq (multiseq_10xv3_cropseq/*.yml, TAPseq_Schraivogel2020/spec.yaml)"},
            "spec_dir": {**_S_STRING, "description": "Directory with the onlists"},
            "fastq_dir": {**_S_STRING, "description": "Directory with the FASTQs (default: spec dir)"},
            "n_reads": {**_S_INTEGER, "default": 10000, "description": "Reads sampled per FASTQ"},
            "strict": {**_S_BOOLEAN, "default": False, "description": "Exit non-zero on failed checks"},
            "offline": {**_S_BOOLEAN, "default": False, "description": "Never fetch from GitHub; use the embedded index"},
            "label": {**_S_STRING, "default": "check_reads", "description": "Run label for the output directory"}},
         "required": ["spec"]},
        cli=["crispr-seqspec", "check-reads"],
        flag_map={"spec": "--spec", "spec_dir": "--spec-dir", "fastq_dir": "--fastq-dir", "n_reads": "--n-reads", "strict": "--strict", "offline": "--offline", "label": "--label"},
        bool_flags={'strict', 'offline'},
    ),

    _T(
        "crispr_seqspec_fetch",
        "Downloads the upstream CRISPR-SeqSpec files (specs, 1k-read FASTQs, onlists; about 20 MB) at the pinned commit into Data/CRISPRSeqSpec, md5-verified, so check and check-reads can use local files.",
        {"type": "object", "properties": {
            "dest": {**_S_STRING, "description": "Destination directory"},
            "specs_only": {**_S_BOOLEAN, "default": False, "description": "Only the YAML specs and readmes"},
            "force": {**_S_BOOLEAN, "default": False, "description": "Re-download existing files"}}},
        cli=["crispr-seqspec", "fetch"],
        flag_map={"dest": "--dest", "specs_only": "--specs-only", "force": "--force"},
        bool_flags={'specs_only', 'force'},
    ),

    _T(
        "gasperini_pipeline_setup",
        "Downloads the upstream guide table (df_from_gasperini_tss.xlsx), the sample pilot MuData (28 MB) and the GEO GSE120861 pilot DE results (19 MB) into Data/GasperiniPipeline; optionally the 10x 737K whitelist. Prints, but never runs, the multi-GB SRA BAM + bamtofastq commands.",
        {"type": "object", "properties": {
            "dest": {**_S_STRING, "description": "Destination directory"},
            "no_geo": {**_S_BOOLEAN, "default": False, "description": "Skip the GEO files"},
            "whitelist": {**_S_BOOLEAN, "default": False, "description": "Also fetch 737K-august-2016.txt"},
            "print_raw": {**_S_BOOLEAN, "default": False, "description": "Print the raw-data download commands"},
            "force": {**_S_BOOLEAN, "default": False, "description": "Re-download"}}},
        cli=["gasperini-pipeline", "setup"],
        flag_map={"dest": "--dest", "no_geo": "--no-geo", "whitelist": "--whitelist", "print_raw": "--print-raw", "force": "--force"},
        bool_flags={'whitelist', 'no_geo', 'force', 'print_raw'},
    ),

    _T(
        "gasperini_pipeline_run",
        "Runs the pipeline's analysis stages on a MuData (default: the fetched Gasperini 2019 sample pilot, 7,314 cells x 98 guides x 2,127 genes): guide assignment at UMI > 3, covariates, cis pairs within 1 Mb, per-guide SCEPTRE-style conditional randomisation tests (NB score statistic, B resamples, skew-normal null), BH and Fisher aggregation to element-gene results. It can also run the original Gasperini NB likelihood-ratio test and compare with GSE120861. Writes processed and results MuData, results TSVs, a QQ/volcano figure, report.md and summary.json.",
        {"type": "object", "properties": {
            "h5mu": {**_S_STRING, "description": "MuData (.h5mu) with 'guides' (raw UMI counts) and 'scRNA' modalities"},
            "guide_umi_limit": {**_S_INTEGER, "default": 3, "description": "Guide assigned when UMI > this"},
            "distance": {**_S_INTEGER, "default": 1000000, "description": "Cis window (bp) around the element's TSS"},
            "in_trans": {**_S_BOOLEAN, "default": False, "description": "Test every gene"},
            "engine": {**_S_STRING, "default": "python", "description": "python (CRT port) or r (SCEPTRE via Rscript)"},
            "resamples": {**_S_INTEGER, "default": 500, "description": "CRT resamples"},
            "direction": {**_S_STRING, "default": "both", "description": "both, left or right"},
            "gasperini_test": {**_S_BOOLEAN, "default": False, "description": "Also run the original NB LRT"},
            "compare": {**_S_STRING, "description": "Original results table (GSE120861_all_deg_results.pilot.txt.gz)"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "Reproduce upstream quirks (random genes for non-gene elements, log(n_genes+1) covariate)"},
            "max_elements": {**_S_INTEGER, "description": "Test only the first N elements"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip figures"},
            "label": {**_S_STRING, "default": "gasperini_run", "description": "Run label for the output directory"}}},
        cli=["gasperini-pipeline", "run"],
        flag_map={"h5mu": "--h5mu", "guide_umi_limit": "--guide-umi-limit", "distance": "--distance", "in_trans": "--in-trans", "engine": "--engine", "resamples": "--resamples", "direction": "--direction", "gasperini_test": "--gasperini-test", "compare": "--compare", "upstream_compat": "--upstream-compat", "max_elements": "--max-elements", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots', 'upstream_compat', 'in_trans', 'gasperini_test'},
    ),

    _T(
        "gasperini_pipeline_qc_filter",
        "Per-lane cell QC and lane merge from kb count outputs (dirs or .h5ad): 200-UMI and knee[EXPECTED_CELL_NUMBER] cuts, MT- fraction < 0.2, scrublet-style doublet removal, barcode intersection with the guide library, batch_number, genes in >= 1% of cells. Writes full_raw_scrna/guide_ann_data.h5ad, lane_qc.tsv and report.md.",
        {"type": "object", "properties": {
            "rna": {**_S_ARRAY_S, "description": "kb transcript outputs, one per lane"},
            "guide": {**_S_ARRAY_S, "description": "kb guide outputs, one per lane"},
            "lanes": {**_S_ARRAY_S, "description": "Lane numbers"},
            "expected_cell_number": {**_S_INTEGER, "default": 10000, "description": "Knee rank"},
            "mito_expected_percentage": {**_S_NUMBER, "default": 0.2, "description": "Max mito fraction"},
            "transcripts_umi_threshold": {**_S_INTEGER, "default": 200, "description": "Min UMIs per cell"},
            "no_doublets": {**_S_BOOLEAN, "default": False, "description": "Skip doublet detection"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "Apply the UMI cutoff as min_genes like upstream"},
            "label": {**_S_STRING, "default": "qc_filter", "description": "Run label for the output directory"}},
         "required": ["rna", "guide"]},
        cli=["gasperini-pipeline", "qc-filter"],
        flag_map={"rna": "--rna", "guide": "--guide", "lanes": "--lanes", "expected_cell_number": "--expected-cell-number", "mito_expected_percentage": "--mito-expected-percentage", "transcripts_umi_threshold": "--transcripts-umi-threshold", "no_doublets": "--no-doublets", "upstream_compat": "--upstream-compat", "label": "--label"},
        flag_repeat={'rna', 'guide', 'lanes'},
        bool_flags={'no_doublets', 'upstream_compat'},
    ),

    _T(
        "gasperini_pipeline_mudata",
        "Builds the raw MuData from QC'd guide and scRNA AnnData: guide chr/start/end/number/target parsed from '<element>|<n>_sgrna_<chr>:<start>:<end>' names, guide sequences, uns elements and varm guide_by_element, and scRNA TSS coordinates from a GTF. Writes raw_mudata_guide_and_transcripts.h5mu.",
        {"type": "object", "properties": {
            "rna_h5ad": {**_S_STRING, "description": "scRNA AnnData"},
            "guide_h5ad": {**_S_STRING, "description": "Guide AnnData"},
            "gtf": {**_S_STRING, "description": "GTF(.gz) for gene coordinates"},
            "guide_features": {**_S_STRING, "description": "guide_features.txt for sequences"},
            "label": {**_S_STRING, "default": "mudata", "description": "Run label for the output directory"}},
         "required": ["rna_h5ad", "guide_h5ad"]},
        cli=["gasperini-pipeline", "mudata"],
        flag_map={"rna_h5ad": "--rna-h5ad", "guide_h5ad": "--guide-h5ad", "gtf": "--gtf", "guide_features": "--guide-features", "label": "--label"},
    ),

    _T(
        "gasperini_pipeline_guide_table",
        "Turns a guide table (xlsx/tsv with Target_name, sgRNA_sequences, chr, start, end; e.g. df_from_gasperini_tss.xlsx) into the pipeline's guide ids '<Target>|<n>_sgrna_<chr>:<start>:<end>' and guide_features.txt for the kite index.",
        {"type": "object", "properties": {
            "table": {**_S_STRING, "description": "Guide table path"},
            "label": {**_S_STRING, "default": "guide_table", "description": "Run label for the output directory"}},
         "required": ["table"]},
        cli=["gasperini-pipeline", "guide-table"],
        flag_map={"table": "--table", "label": "--label"},
    ),

    _T(
        "gasperini_pipeline_composition",
        "Computes per-position nucleotide composition of FASTQs (A/C/T/G counts and their SD, the pipeline's compositional-bias QC) from the first 10,000 lines; writes a TSV and figure per FASTQ.",
        {"type": "object", "properties": {
            "fastq": {**_S_ARRAY_S, "description": "FASTQ(.gz) files"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip figures"},
            "label": {**_S_STRING, "default": "composition", "description": "Run label for the output directory"}},
         "required": ["fastq"]},
        cli=["gasperini-pipeline", "composition"],
        flag_map={"fastq": "--fastq", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'fastq'},
        bool_flags={'no_plots'},
    ),

    _T(
        "gasperini_pipeline_multiseq",
        "Demultiplexes MULTI-seq sample barcodes (deMULTIplex port): extracts cell/UMI/tag from R1/R2, matches tags with Hamming <= 1, counts unique UMIs per cell x tag, classifies cells by KDE thresholds with a quantile sweep over two rounds (singlet tag, Doublet, Negative) and optionally filters a MuData. Writes bar_table.tsv, final_class.tsv and report.md.",
        {"type": "object", "properties": {
            "r1": {**_S_STRING, "description": "MULTI-seq R1 FASTQ"},
            "r2": {**_S_STRING, "description": "MULTI-seq R2 FASTQ"},
            "cell_barcodes": {**_S_STRING, "description": "Cell barcodes, one per line"},
            "bar_ref": {**_S_STRING, "description": "CSV of MULTI-seq barcodes"},
            "h5mu": {**_S_STRING, "description": "MuData to filter"},
            "label": {**_S_STRING, "default": "multiseq", "description": "Run label for the output directory"}},
         "required": ["r1", "r2", "cell_barcodes", "bar_ref"]},
        cli=["gasperini-pipeline", "multiseq"],
        flag_map={"r1": "--r1", "r2": "--r2", "cell_barcodes": "--cell-barcodes", "bar_ref": "--bar-ref", "h5mu": "--h5mu", "label": "--label"},
    ),

    _T(
        "gasperini_pipeline_compare",
        "Compares element-gene (or NB) results with the original Gasperini 2019 table (GSE120861 all_deg_results): joins on gene symbol and element + '_TSS', reports Spearman of -log10 p, BH hit overlap with Fisher exact test, effect sign agreement and self-TSS recovery. Writes matched_pairs.tsv, summary.json and report.md.",
        {"type": "object", "properties": {
            "results": {**_S_STRING, "description": "element_gene_results.tsv, results.tsv or gasperini_nb_results.tsv"},
            "original": {**_S_STRING, "description": "Original results table"},
            "fdr": {**_S_NUMBER, "default": 0.1, "description": "BH threshold"},
            "label": {**_S_STRING, "default": "compare", "description": "Run label for the output directory"}},
         "required": ["results", "original"]},
        cli=["gasperini-pipeline", "compare"],
        flag_map={"results": "--results", "original": "--original", "fdr": "--fdr", "label": "--label"},
    ),

    _T(
        "gasperini_pipeline_inspect",
        "Checks a MuData against the IGVF sample-pilot schema (guides obs/var/uns/varm fields, scRNA obs/var fields, shared barcodes, integer guide counts, guide_by_element consistency, guide naming). Writes schema_checks.tsv and report.md.",
        {"type": "object", "properties": {
            "h5mu": {**_S_STRING, "description": "MuData path (default: fetched sample pilot)"},
            "label": {**_S_STRING, "default": "inspect", "description": "Run label for the output directory"}}},
        cli=["gasperini-pipeline", "inspect"],
        flag_map={"h5mu": "--h5mu", "label": "--label"},
    ),

    _T(
        "abc_pipeline_setup",
        "Downloads the upstream ABC reference files (abc_thresholds.tsv, the K562 quantile-normalisation reference, UbiquitouslyExpressedGenes, hg38/hg19 gene bounds, TSS500bp include-lists, chromosome sizes, blocklists) from the pinned broadinstitute/ABC-Enhancer-Gene-Prediction commit into Data/ABCPipeline/reference, where the other abc-pipeline tools pick them up as defaults. Use dry_run to list the URLs only.",
        {"type": "object", "properties": {
            "dest": {**_S_STRING, "description": "destination directory (default Data/ABCPipeline/reference)"},
            "dry_run": {**_S_BOOLEAN, "default": False, "description": "list files without downloading"},
            "force": {**_S_BOOLEAN, "default": False, "description": "re-download present files"}}},
        cli=["abc-pipeline", "setup"],
        flag_map={"dest": "--dest", "dry_run": "--dry-run", "force": "--force"},
        bool_flags={'dry_run', 'force'},
    ),

    _T(
        "abc_pipeline_call_peaks",
        "Calls accessibility peaks for ABC exactly as rules/macs2.smk (macs2 callpeak -p 0.1 --shift -75 --extsize 150 --nomodel --keep-dup all --call-summits), then keeps peaks on the chrom-sizes chromosomes and sorts them (macs2_peaks.narrowPeak.sorted). Without macs2 on PATH it prints the exact command and exits 0, or with python_fallback runs a MACS-like Python caller (approximation). Inputs: DNase BAM or ATAC tagAlign.gz; output: narrowPeak + report.",
        {"type": "object", "properties": {
            "accessibility": {**_S_ARRAY_S, "description": "DNase BAM / ATAC tagAlign.gz / fragments files"},
            "chrom_sizes": {**_S_STRING, "description": "chromosome sizes TSV"},
            "pval": {**_S_NUMBER, "default": 0.1, "description": "macs2 -p"},
            "genome_size": {**_S_STRING, "default": "hs", "description": "macs2 -g (hs, mm, or a number)"},
            "python_fallback": {**_S_BOOLEAN, "default": False, "description": "use the Python MACS-like caller when macs2 is absent"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["accessibility", "chrom_sizes"]},
        cli=["abc-pipeline", "call-peaks"],
        flag_map={"accessibility": "--accessibility", "chrom_sizes": "--chrom-sizes", "pval": "--pval", "genome_size": "--genome-size", "python_fallback": "--python-fallback", "label": "--label", "outdir": "--outdir"},
        flag_repeat={'accessibility'},
        bool_flags={'python_fallback'},
    ),

    _T(
        "abc_pipeline_fragments_to_tagalign",
        "Converts a 10x scATAC fragments file into an ABC tagAlign.gz: each fragment becomes two tags (start..mid, +) and (mid+1..end, -), as in the upstream scATAC walkthrough. Input: fragments.tsv.gz; output: tagAlign.gz usable as ATAC input for every abc-pipeline step.",
        {"type": "object", "properties": {
            "fragments": {**_S_STRING, "description": "fragments.tsv(.gz)"},
            "out": {**_S_STRING, "description": "output tagAlign.gz path"}},
         "required": ["fragments", "out"]},
        cli=["abc-pipeline", "fragments-to-tagalign"],
        flag_map={"fragments": "--fragments", "out": "--out"},
    ),

    _T(
        "abc_pipeline_candidate_regions",
        "Defines ABC candidate elements exactly as makeCandidateRegions.py: counts reads in every narrowPeak peak (averaging several files; chrX/Y doubled), keeps the n_strongest_peaks (150000) with most reads, extends summits +/-250 bp, merges, removes blocklisted regions and forces in the TSS include-list. Inputs: narrowPeak with summits, accessibility reads, chrom sizes; outputs: <peaks>.candidateRegions.bed and the Counts.bed.",
        {"type": "object", "properties": {
            "narrowpeak": {**_S_STRING, "description": "macs2 narrowPeak(.sorted) with summit column"},
            "accessibility": {**_S_ARRAY_S, "description": "DNase/ATAC BAM, tagAlign or fragments"},
            "chrom_sizes": {**_S_STRING, "description": "chromosome sizes"},
            "includelist": {**_S_STRING, "description": "TSS include-list BED (e.g. CollapsedGeneBounds.hg38.TSS500bp.bed)"},
            "blocklist": {**_S_STRING, "description": "blocklist BED"},
            "n_strongest_peaks": {**_S_INTEGER, "default": 150000, "description": "peaks kept by read count"},
            "peak_extend_from_summit": {**_S_INTEGER, "default": 250, "description": "bp each side of the summit"},
            "ignore_summits": {**_S_BOOLEAN, "default": False, "description": "extend whole peaks instead of summits"},
            "min_peak_width": {**_S_INTEGER, "default": 500, "description": "minimum width with ignore_summits"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["narrowpeak", "accessibility", "chrom_sizes"]},
        cli=["abc-pipeline", "candidate-regions"],
        flag_map={"narrowpeak": "--narrowpeak", "accessibility": "--accessibility", "chrom_sizes": "--chrom-sizes", "includelist": "--includelist", "blocklist": "--blocklist", "n_strongest_peaks": "--n-strongest-peaks", "peak_extend_from_summit": "--peak-extend-from-summit", "ignore_summits": "--ignore-summits", "min_peak_width": "--min-peak-width", "label": "--label", "outdir": "--outdir"},
        flag_repeat={'accessibility'},
        bool_flags={'ignore_summits'},
    ),

    _T(
        "abc_pipeline_neighborhoods",
        "Runs run.neighborhoods.py: counts DHS/ATAC/H3K27ac reads (BAM, tagAlign, fragments, bigWig or bedGraph) in candidate elements, gene bodies and TSS +/-500 bp; computes RPM/RPKM/quantiles, promoter/genic/intergenic classes, quantile-normalises to the K562 reference (rank method, promoters and non-promoters separately) and activity_base = sqrt(normalized H3K27ac x normalized accessibility). Outputs EnhancerList.txt, GeneList.txt (PromoterActivityQuantile, is_ue, Expression) and the CountReads bedgraphs.",
        {"type": "object", "properties": {
            "candidate_regions": {**_S_STRING, "description": "candidateRegions.bed"},
            "genes": {**_S_STRING, "description": "gene BED6 + Ensembl_ID + gene_type"},
            "chrom_sizes": {**_S_STRING, "description": "chromosome sizes"},
            "dhs": {**_S_ARRAY_S, "description": "DNase files"},
            "atac": {**_S_ARRAY_S, "description": "ATAC files"},
            "h3k27ac": {**_S_ARRAY_S, "description": "H3K27ac files"},
            "default_accessibility_feature": {**_S_STRING, "description": "DHS or ATAC"},
            "ubiquitously_expressed_genes": {**_S_STRING, "description": "ubiquitous gene list"},
            "expression_table": {**_S_STRING, "description": "comma-separated gene<TAB>value tables"},
            "qnorm": {**_S_STRING, "description": "quantile-normalisation reference (default: set-up K562 reference)"},
            "no_qnorm": {**_S_BOOLEAN, "default": False, "description": "skip quantile normalisation"},
            "qnorm_method": {**_S_STRING, "default": "rank", "description": "rank or quantile"},
            "cell_type": {**_S_STRING, "description": "cell type label"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "read the ubiquitous list with a header as upstream"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["candidate_regions", "genes", "chrom_sizes"]},
        cli=["abc-pipeline", "neighborhoods"],
        flag_map={"candidate_regions": "--candidate-regions", "genes": "--genes", "chrom_sizes": "--chrom-sizes", "dhs": "--dhs", "atac": "--atac", "h3k27ac": "--h3k27ac", "default_accessibility_feature": "--default-accessibility-feature", "ubiquitously_expressed_genes": "--ubiquitously-expressed-genes", "expression_table": "--expression-table", "qnorm": "--qnorm", "no_qnorm": "--no-qnorm", "qnorm_method": "--qnorm-method", "cell_type": "--cell-type", "upstream_compat": "--upstream-compat", "label": "--label", "outdir": "--outdir"},
        flag_repeat={'h3k27ac', 'dhs', 'atac'},
        bool_flags={'no_qnorm', 'upstream_compat'},
    ),

    _T(
        "abc_pipeline_predict",
        "Computes ABC scores exactly as predict.py: every element within 5 Mb of each TSS, power-law contact exp(scale - gamma ln(d+1)) with the config gamma/scale, or Hi-C contact (hic via hicstraw, juicebox, bedpe or avg directories) made doubly stochastic, diagonal-corrected, rescaled to the reference power law and given the power-law pseudocount; ABC.Score = activity x contact / per-gene sum, self-promoters = 1. Inputs: EnhancerList.txt, GeneList.txt; outputs: EnhancerPredictionsAllPutative(.NonExpressedGenes).tsv.gz with upstream columns and the variant-overlap file.",
        {"type": "object", "properties": {
            "enhancers": {**_S_STRING, "description": "EnhancerList.txt"},
            "genes": {**_S_STRING, "description": "GeneList.txt"},
            "chrom_sizes": {**_S_STRING, "description": "chromosome sizes"},
            "accessibility_feature": {**_S_STRING, "description": "DHS or ATAC"},
            "cell_type": {**_S_STRING, "description": "CellType column"},
            "hic_file": {**_S_STRING, "description": ".hic file/URL or Hi-C directory"},
            "hic_type": {**_S_STRING, "default": "hic", "description": "hic, juicebox, bedpe or avg"},
            "hic_resolution": {**_S_INTEGER, "description": "Hi-C resolution (5000)"},
            "hic_gamma": {**_S_NUMBER, "default": 1.024238616787792, "description": "power-law gamma"},
            "hic_scale": {**_S_NUMBER, "default": 5.9594510043736655, "description": "power-law scale"},
            "score_column": {**_S_STRING, "default": "ABC.Score", "description": "score for the variant-overlap file"},
            "window": {**_S_INTEGER, "default": 5000000, "description": "max TSS distance"},
            "no_scale_hic_using_powerlaw": {**_S_BOOLEAN, "default": False, "description": "do not rescale Hi-C to the reference power law"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip the figure"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["enhancers", "genes", "chrom_sizes", "accessibility_feature"]},
        cli=["abc-pipeline", "predict"],
        flag_map={"enhancers": "--enhancers", "genes": "--genes", "chrom_sizes": "--chrom-sizes", "accessibility_feature": "--accessibility-feature", "cell_type": "--cell-type", "hic_file": "--hic-file", "hic_type": "--hic-type", "hic_resolution": "--hic-resolution", "hic_gamma": "--hic-gamma", "hic_scale": "--hic-scale", "score_column": "--score-column", "window": "--window", "no_scale_hic_using_powerlaw": "--no-scale-hic-using-powerlaw", "no_plots": "--no-plots", "label": "--label", "outdir": "--outdir"},
        bool_flags={'no_plots', 'no_scale_hic_using_powerlaw'},
    ),

    _T(
        "abc_pipeline_filter",
        "Thresholds ABC predictions exactly as filter_predictions.py: picks the threshold from abc_thresholds.tsv by accessibility, H3K27ac and Hi-C type (e.g. DHS+H3K27ac intact Hi-C 0.027, power law 0.017; 0.02 otherwise) unless given, keeps score > threshold, drops promoters except self-promoters. Outputs EnhancerPredictionsFull_threshold<t>_self_promoter.tsv/.bedpe.gz, the slim EnhancerPredictions file and GenePredictionStats.",
        {"type": "object", "properties": {
            "pred_file": {**_S_STRING, "description": "EnhancerPredictionsAllPutative.tsv.gz"},
            "pred_nonexpressed_file": {**_S_STRING, "description": "NonExpressedGenes file (default: alongside)"},
            "score_column": {**_S_STRING, "default": "ABC.Score", "description": "score column"},
            "accessibility_feature": {**_S_STRING, "default": "DHS", "description": "DHS or ATAC"},
            "has_h3k27ac": {**_S_BOOLEAN, "default": False, "description": "activity used H3K27ac"},
            "hic_type": {**_S_STRING, "description": "hic, avg, juicebox, bedpe; omit for power law"},
            "threshold": {**_S_NUMBER, "description": "explicit threshold"},
            "exclude_self_promoter": {**_S_BOOLEAN, "default": False, "description": "drop self-promoters"},
            "only_expressed_genes": {**_S_BOOLEAN, "default": False, "description": "expressed genes only"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "reproduce upstream's avg_hic threshold lookup"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["pred_file"]},
        cli=["abc-pipeline", "filter"],
        flag_map={"pred_file": "--pred-file", "pred_nonexpressed_file": "--pred-nonexpressed-file", "score_column": "--score-column", "accessibility_feature": "--accessibility-feature", "has_h3k27ac": "--has-h3k27ac", "hic_type": "--hic-type", "threshold": "--threshold", "exclude_self_promoter": "--exclude-self-promoter", "only_expressed_genes": "--only-expressed-genes", "upstream_compat": "--upstream-compat", "label": "--label", "outdir": "--outdir"},
        bool_flags={'has_h3k27ac', 'upstream_compat', 'only_expressed_genes', 'exclude_self_promoter'},
    ),

    _T(
        "abc_pipeline_variant_overlap",
        "Writes the upstream EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz from an AllPutative table: non-promoters with score > 0.015 or promoters with score > 0.1, within 2 Mb, each element trimmed by 150 bp on both ends for variant intersection.",
        {"type": "object", "properties": {
            "all_putative": {**_S_STRING, "description": "EnhancerPredictionsAllPutative.tsv.gz"},
            "score_column": {**_S_STRING, "default": "ABC.Score", "description": "score column"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["all_putative"]},
        cli=["abc-pipeline", "variant-overlap"],
        flag_map={"all_putative": "--all-putative", "score_column": "--score-column", "label": "--label", "outdir": "--outdir"},
    ),

    _T(
        "abc_pipeline_qc",
        "Computes the upstream ABC QC (grabMetrics.py/metrics.py): enhancers per gene, genes per enhancer, pairs per chromosome, E-G distance quantiles, peak and candidate-region widths, read counts per feature, and distribution / Hi-C power-law plots. Inputs: thresholded predictions, candidate regions, Neighborhoods dir; outputs QCSummary_<name>.tsv, per-gene tables, PNGs and QCPlots PDF.",
        {"type": "object", "properties": {
            "preds_file": {**_S_STRING, "description": "EnhancerPredictionsFull_*.tsv"},
            "candidate_regions": {**_S_STRING, "description": "candidateRegions.bed"},
            "neighborhood_dir": {**_S_STRING, "description": "Neighborhoods directory"},
            "chrom_sizes": {**_S_STRING, "description": "chromosome sizes"},
            "narrowpeak": {**_S_STRING, "description": "optional real peaks for the peak-width metrics"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip plots"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["preds_file", "candidate_regions", "chrom_sizes"]},
        cli=["abc-pipeline", "qc"],
        flag_map={"preds_file": "--preds-file", "candidate_regions": "--candidate-regions", "neighborhood_dir": "--neighborhood-dir", "chrom_sizes": "--chrom-sizes", "narrowpeak": "--narrowpeak", "no_plots": "--no-plots", "label": "--label", "outdir": "--outdir"},
        bool_flags={'no_plots'},
    ),

    _T(
        "abc_pipeline_powerlaw_fit",
        "Fits the Hi-C contact-distance power law as compute_powerlaw_fit_from_hic.py: sums contacts per distance between 5 kb and 1 Mb, normalises by the number of entries at the resolution, and regresses log contact on log distance; hic_gamma = -slope, hic_scale = intercept (feed them to predict). Inputs: juicebox, bedpe or avg Hi-C directory; outputs hic.powerlaw.tsv and hic.mean_var.tsv.",
        {"type": "object", "properties": {
            "hic_dir": {**_S_STRING, "description": "Hi-C directory (<chr>/<chr>.KRobserved.gz etc.)"},
            "hic_type": {**_S_STRING, "default": "juicebox", "description": "juicebox, bedpe or avg"},
            "hic_resolution": {**_S_INTEGER, "default": 5000, "description": "resolution"},
            "min_window": {**_S_INTEGER, "default": 5000, "description": "min distance"},
            "max_window": {**_S_INTEGER, "default": 1000000, "description": "max distance"},
            "chr": {**_S_STRING, "default": "all", "description": "comma-separated chromosomes or all"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "avg: upstream's distance-in-bins"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip plot"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["hic_dir"]},
        cli=["abc-pipeline", "powerlaw-fit"],
        flag_map={"hic_dir": "--hic-dir", "hic_type": "--hic-type", "hic_resolution": "--hic-resolution", "min_window": "--min-window", "max_window": "--max-window", "chr": "--chr", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label", "outdir": "--outdir"},
        bool_flags={'no_plots', 'upstream_compat'},
    ),

    _T(
        "abc_pipeline_average_hic",
        "Builds a cell-type-averaged Hi-C matrix for one chromosome as makeAverageHiC.py: each cell type's KR matrix made doubly stochastic, rescaled to the reference power law with its own fitted gamma/scale, outer-joined and averaged over non-missing entries, NaN when fewer than min_cell_types_required. Output <chr>/<chr>.avg.gz.",
        {"type": "object", "properties": {
            "celltypes": {**_S_STRING, "description": "comma-separated cell types under basedir/<ct>/5kb_resolution_intra"},
            "basedir": {**_S_STRING, "description": "base directory"},
            "celltype_dirs": {**_S_ARRAY_S, "description": "explicit juicebox dirs with powerlaw/hic.powerlaw.txt"},
            "chromosome": {**_S_STRING, "description": "chromosome"},
            "min_cell_types_required": {**_S_INTEGER, "default": 3, "description": "minimum cell types per entry"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["chromosome"]},
        cli=["abc-pipeline", "average-hic"],
        flag_map={"celltypes": "--celltypes", "basedir": "--basedir", "celltype_dirs": "--celltype-dirs", "chromosome": "--chromosome", "min_cell_types_required": "--min-cell-types-required", "label": "--label", "outdir": "--outdir"},
        flag_repeat={'celltype_dirs'},
    ),

    _T(
        "abc_pipeline_split_avg_hic",
        "Splits a genome-wide average Hi-C bed.gz (e.g. ENCODE ENCFF134PUN) into AvgHiC/<chr>/<chr>.bed.gz per chromosome, the directory layout predict uses with hic_type avg (extract_avg_hic.py).",
        {"type": "object", "properties": {
            "avg_hic_bed_file": {**_S_STRING, "description": "genome-wide avg Hi-C bed.gz"},
            "output_dir": {**_S_STRING, "default": ".", "description": "output directory"}},
         "required": ["avg_hic_bed_file"]},
        cli=["abc-pipeline", "split-avg-hic"],
        flag_map={"avg_hic_bed_file": "--avg-hic-bed-file", "output_dir": "--output-dir"},
    ),

    _T(
        "abc_pipeline_juicebox_dump",
        "Dumps a .hic file to the juicebox directory layout (KR observed and KR norm per chromosome at 5 kb) with juicer_tools, as juicebox_dump.py; without juicer_tools it prints the exact commands and exits 0.",
        {"type": "object", "properties": {
            "hic_file": {**_S_STRING, "description": ".hic path or URL"},
            "juicebox": {**_S_STRING, "default": "java -jar juicer_tools.jar", "description": "juicer_tools command"},
            "resolution": {**_S_INTEGER, "default": 5000, "description": "resolution"},
            "outdir": {**_S_STRING, "default": ".", "description": "output directory"},
            "chromosomes": {**_S_STRING, "default": "all", "description": "comma-separated or all"},
            "include_raw": {**_S_BOOLEAN, "default": False, "description": "also dump raw observed"}},
         "required": ["hic_file"]},
        cli=["abc-pipeline", "juicebox-dump"],
        flag_map={"hic_file": "--hic-file", "juicebox": "--juicebox", "resolution": "--resolution", "outdir": "--outdir", "chromosomes": "--chromosomes", "include_raw": "--include-raw"},
        bool_flags={'include_raw'},
    ),

    _T(
        "abc_pipeline_hic_bedgraph",
        "Writes one bedGraph per gene of the Hi-C row at its TSS bin within 5 Mb (make_bedgraph_from_HiC.py), with bins below the KR-norm cutoff interpolated; for browsing contacts around a gene.",
        {"type": "object", "properties": {
            "genes": {**_S_STRING, "description": "gene BED"},
            "hic_dir": {**_S_STRING, "description": "juicebox Hi-C directory"},
            "outdir": {**_S_STRING, "description": "output directory"},
            "resolution": {**_S_INTEGER, "default": 5000, "description": "resolution"},
            "kr_cutoff": {**_S_NUMBER, "default": 0.1, "description": "KR norm cutoff"}},
         "required": ["genes", "hic_dir", "outdir"]},
        cli=["abc-pipeline", "hic-bedgraph"],
        flag_map={"genes": "--genes", "hic_dir": "--hic-dir", "outdir": "--outdir", "resolution": "--resolution", "kr_cutoff": "--kr-cutoff"},
    ),

    _T(
        "abc_pipeline_compare",
        "Compares two ABC prediction tables the way the upstream end-to-end test does: rows matched on chr/start/end/TargetGene, the tested columns (name, class, ABC.Score.Numerator, ABC.Score, powerlaw.Score) compared within atol; reports max differences, correlations and a scatter plot.",
        {"type": "object", "properties": {
            "test": {**_S_STRING, "description": "prediction file under test"},
            "expected": {**_S_STRING, "description": "reference prediction file"},
            "atol": {**_S_NUMBER, "default": 1e-06, "description": "absolute tolerance"},
            "fail_on_diff": {**_S_BOOLEAN, "default": False, "description": "exit 1 when different"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip plot"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}},
         "required": ["test", "expected"]},
        cli=["abc-pipeline", "compare"],
        flag_map={"test": "--test", "expected": "--expected", "atol": "--atol", "fail_on_diff": "--fail-on-diff", "no_plots": "--no-plots", "label": "--label", "outdir": "--outdir"},
        bool_flags={'no_plots', 'fail_on_diff'},
    ),

    _T(
        "abc_pipeline_run",
        "Runs the whole upstream ABC Snakefile for one biosample or an upstream-format biosamples table: peaks (macs2 or the Python fallback, or a given narrowPeak), candidate regions, neighborhoods, predictions (power law or Hi-C), automatic threshold, filtering and QC, laid out as <biosample>/Peaks, Neighborhoods, Predictions, Metrics with report.md and summary.json. References default to the files fetched by abc_pipeline_setup.",
        {"type": "object", "properties": {
            "biosamples_table": {**_S_STRING, "description": "TSV with biosample, DHS, ATAC, H3K27ac, default_accessibility_feature, HiC_file, HiC_type, HiC_resolution, alt_TSS, alt_genes"},
            "biosample": {**_S_STRING, "description": "biosample name (single-sample mode)"},
            "dhs": {**_S_ARRAY_S, "description": "DNase files"},
            "atac": {**_S_ARRAY_S, "description": "ATAC files"},
            "h3k27ac": {**_S_ARRAY_S, "description": "H3K27ac files"},
            "hic_file": {**_S_STRING, "description": "Hi-C file or directory"},
            "hic_type": {**_S_STRING, "description": "hic, juicebox, bedpe or avg"},
            "hic_resolution": {**_S_INTEGER, "description": "Hi-C resolution (5000)"},
            "genes": {**_S_STRING, "description": "gene BED"},
            "tss": {**_S_STRING, "description": "TSS include-list"},
            "chrom_sizes": {**_S_STRING, "description": "chromosome sizes"},
            "blocklist": {**_S_STRING, "description": "blocklist BED"},
            "qnorm": {**_S_STRING, "description": "qnorm reference"},
            "narrowpeak": {**_S_STRING, "description": "use this narrowPeak instead of calling peaks"},
            "python_fallback": {**_S_BOOLEAN, "default": False, "description": "Python peak caller when macs2 is absent"},
            "threshold": {**_S_NUMBER, "description": "explicit threshold"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip plots"},
            "label": {**_S_STRING, "description": "run label; output goes to Docs/ABCPipeline/<timestamp>_<label>"},
            "outdir": {**_S_STRING, "description": "explicit output directory"}}},
        cli=["abc-pipeline", "run"],
        flag_map={"biosamples_table": "--biosamples-table", "biosample": "--biosample", "dhs": "--dhs", "atac": "--atac", "h3k27ac": "--h3k27ac", "hic_file": "--hic-file", "hic_type": "--hic-type", "hic_resolution": "--hic-resolution", "genes": "--genes", "tss": "--tss", "chrom_sizes": "--chrom-sizes", "blocklist": "--blocklist", "qnorm": "--qnorm", "narrowpeak": "--narrowpeak", "python_fallback": "--python-fallback", "threshold": "--threshold", "no_plots": "--no-plots", "label": "--label", "outdir": "--outdir"},
        flag_repeat={'h3k27ac', 'dhs', 'atac'},
        bool_flags={'no_plots', 'python_fallback'},
    ),

    _T(
        "abc_pipeline_selftest",
        "Runs the abc-pipeline self-test: a synthetic chromosome with planted enhancers, reads, a power-law Hi-C map with a planted loop, then every subcommand with assertions on the upstream definitions (candidate regions, counting, qnorm, ABC formula, Hi-C loaders, thresholds, QC). Returns pass/fail per check.",
        {"type": "object", "properties": {
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"}}},
        cli=["abc-pipeline", "selftest"],
        flag_map={"no_plots": "--no-plots"},
        bool_flags={'no_plots'},
    ),

    _T(
        "sce2g_pipeline_setup",
        "Downloads the data the scE2G port needs from the pinned upstream commits: the four v3 models' qnorm references into Data/scE2G/models/<model>/, plus TSS500 / gene bounds / chrom sizes / gene classes / Sheth-Qiu QC reference into Data/scE2G/pipeline_resources. Optional flags add the CRISPR benchmark, the GENCODE v43 GTF (53 MB, needed to map RNA gene names) and the upstream chr22 test fixture.",
        {"type": "object", "properties": {
            "crispr": {**_S_BOOLEAN, "default": False, "description": "Also fetch the EPCrisprBenchmark CRISPR table."},
            "gtf": {**_S_BOOLEAN, "default": False, "description": "Also fetch the GENCODE v43 GTF (53 MB)."},
            "example": {**_S_BOOLEAN, "default": False, "description": "Also fetch the chr22 test fixture for validate-example."},
            "force": {**_S_BOOLEAN, "default": False, "description": "Re-download files already present."},
            "dry_run": {**_S_BOOLEAN, "default": False, "description": "Only list what would be fetched."}}},
        cli=["sce2g-pipeline", "setup"],
        flag_map={"crispr": "--crispr", "gtf": "--gtf", "example": "--example", "force": "--force", "dry_run": "--dry-run"},
        bool_flags={'dry_run', 'gtf', 'crispr', 'force', 'example'},
    ),

    _T(
        "sce2g_pipeline_models",
        "Lists the four embedded scE2G v3 models (feature lists, logistic-regression weights of the full and 23 held-out-chromosome models, score and TPM thresholds, qnorm reference status). Writes model_weights.tsv and a report; with export, writes upstream-shaped model directories usable by predict.",
        {"type": "object", "properties": {
            "export": {**_S_STRING, "description": "Directory to write model directories (feature_table.tsv, threshold files, model_coefficients.json)."},
            "label": {**_S_STRING, "default": "models", "description": "Run label (output directory suffix)."}}},
        cli=["sce2g-pipeline", "models"],
        flag_map={"export": "--export", "label": "--label"},
    ),

    _T(
        "sce2g_pipeline_run",
        "Runs the whole scE2G pipeline for one cluster (or a config_cell_clusters.tsv): fragments -> tagAlign and counts, Kendall peak-gene pairs, peak x cell ATAC matrix, Kendall tau-b against log-normalised RNA, ARC-E2G, ENCODE_rE2G activity features, genome-wide features, model application with qnorm and TPM filter, thresholded links, gene/element lists, per-cluster stats and QC vs the Sheth-Qiu reference. Needs the cluster's ABC outputs (--abc-dir with Peaks/, Predictions/, Neighborhoods/); without them it writes the tagAlign and prints the ABC command. With crispr it also benchmarks against CRISPR pairs.",
        {"type": "object", "properties": {
            "cluster": {**_S_STRING, "default": "cluster", "description": "Cluster / biosample name."},
            "fragments": {**_S_STRING, "description": "10x ATAC fragments file (.tsv.gz)."},
            "rna_matrix": {**_S_STRING, "description": "RNA counts: .h5ad, .csv.gz (genes x cells) or 10x directory."},
            "abc_dir": {**_S_STRING, "description": "ABC results directory for the cluster."},
            "cluster_config": {**_S_STRING, "description": "Upstream config_cell_clusters.tsv with an abc_dir column (multi-cluster run)."},
            "models": {**_S_ARRAY_S, "description": "Embedded models (multiome_powerlaw_v3, multiome_megamap_v3, scATAC_powerlaw_v3, scATAC_megamap_v3) or model directories; default multiome_powerlaw_v3 scATAC_powerlaw_v3."},
            "hic": {**_S_BOOLEAN, "default": False, "description": "Cluster has Hi-C: ARC uses ABC.Score instead of powerlaw.Score."},
            "gtf": {**_S_STRING, "description": "GENCODE GTF for RNA gene-name mapping (default: setup --gtf)."},
            "no_gene_mapping": {**_S_BOOLEAN, "default": False, "description": "RNA gene names already match the TSS500 reference."},
            "rna_unfiltered": {**_S_BOOLEAN, "default": False, "description": "RNA matrix contains extra cells; keep those present in the fragments."},
            "crispr": {**_S_STRING, "description": "CRISPR benchmark table; adds cv scores, CRISPR features and the benchmark."},
            "igv_tracks": {**_S_BOOLEAN, "default": False, "description": "Also write ATAC_norm coverage and bedpe links."},
            "n_boot": {**_S_INTEGER, "default": 1000, "description": "Bootstrap resamples for the benchmark."},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip figures."},
            "label": {**_S_STRING, "default": "run", "description": "Run label (output directory suffix)."}}},
        cli=["sce2g-pipeline", "run"],
        flag_map={"cluster": "--cluster", "fragments": "--fragments", "rna_matrix": "--rna-matrix", "abc_dir": "--abc-dir", "cluster_config": "--cluster-config", "models": "--models", "hic": "--hic", "gtf": "--gtf", "no_gene_mapping": "--no-gene-mapping", "rna_unfiltered": "--rna-unfiltered", "crispr": "--crispr", "igv_tracks": "--igv-tracks", "n_boot": "--n-boot", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'models'},
        bool_flags={'hic', 'no_gene_mapping', 'igv_tracks', 'rna_unfiltered', 'no_plots'},
    ),

    _T(
        "sce2g_pipeline_frag_to_tagalign",
        "Converts an ATAC fragment file into the scE2G/ABC inputs: fragments restricted to chrom-sizes chromosomes (bgzip + tabix), fragment_count.txt, cell_barcodes.txt and a sorted tagAlign with two tags per fragment split at the fragment midpoint.",
        {"type": "object", "properties": {
            "fragments": {**_S_STRING, "description": "Fragments file."},
            "chrom_sizes": {**_S_STRING, "description": "Chromosome sizes (default: setup)."},
            "cluster": {**_S_STRING, "default": "cluster", "description": "Cluster name."},
            "preprocessed": {**_S_BOOLEAN, "default": False, "description": "Fragments already sorted and filtered."},
            "label": {**_S_STRING, "default": "tagalign", "description": "Run label (output directory suffix)."}},
         "required": ["fragments"]},
        cli=["sce2g-pipeline", "frag-to-tagalign"],
        flag_map={"fragments": "--fragments", "chrom_sizes": "--chrom-sizes", "cluster": "--cluster", "preprocessed": "--preprocessed", "label": "--label"},
        bool_flags={'preprocessed'},
    ),

    _T(
        "sce2g_pipeline_kendall",
        "Computes scE2G's Kendall tau-b between each candidate peak's binarised accessibility and its gene's log-normalised expression across cells (exact closed form of the Sheth-Qiu algorithm), plus RNA pseudobulk TPM, percent cells detected and mean log-normalised expression per gene. Inputs: Kendall pairs (from kendall-pairs), an ATAC matrix or fragments, and the RNA matrix; outputs Pairs.Kendall.tsv.gz, gene_expression_metrics.tsv.gz, umi_count and cell_count.",
        {"type": "object", "properties": {
            "pairs": {**_S_STRING, "description": "Kendall/Pairs.tsv.gz."},
            "rna_matrix": {**_S_STRING, "description": "RNA counts matrix."},
            "fragments": {**_S_STRING, "description": "Fragments file (if no atac_matrix)."},
            "atac_matrix": {**_S_STRING, "description": "Prefix written by atac-matrix."},
            "gtf": {**_S_STRING, "description": "GENCODE GTF for gene-name mapping."},
            "no_gene_mapping": {**_S_BOOLEAN, "default": False, "description": "RNA names already are TSS500 names."},
            "cluster": {**_S_STRING, "default": "cluster", "description": "Cluster name."},
            "label": {**_S_STRING, "default": "kendall", "description": "Run label (output directory suffix)."}},
         "required": ["pairs", "rna_matrix"]},
        cli=["sce2g-pipeline", "kendall"],
        flag_map={"pairs": "--pairs", "rna_matrix": "--rna-matrix", "fragments": "--fragments", "atac_matrix": "--atac-matrix", "gtf": "--gtf", "no_gene_mapping": "--no-gene-mapping", "cluster": "--cluster", "label": "--label"},
        bool_flags={'no_gene_mapping'},
    ),

    _T(
        "sce2g_pipeline_kendall_pairs",
        "Builds scE2G's candidate peak-gene pairs: MACS2 narrowPeak peaks crossed with the ABC EnhancerPredictionsAllPutative pairs they overlap, with PeakName and PairName; writes Kendall/Pairs.tsv.gz.",
        {"type": "object", "properties": {
            "narrowpeak": {**_S_STRING, "description": "ABC Peaks/macs2_peaks.narrowPeak.sorted."},
            "abc_predictions": {**_S_STRING, "description": "ABC Predictions/EnhancerPredictionsAllPutative.tsv.gz."},
            "cluster": {**_S_STRING, "default": "cluster", "description": "Cluster name."},
            "label": {**_S_STRING, "default": "kendall_pairs", "description": "Run label (output directory suffix)."}},
         "required": ["narrowpeak", "abc_predictions"]},
        cli=["sce2g-pipeline", "kendall-pairs"],
        flag_map={"narrowpeak": "--narrowpeak", "abc_predictions": "--abc-predictions", "cluster": "--cluster", "label": "--label"},
    ),

    _T(
        "sce2g_pipeline_arc",
        "Computes ARC-E2G: joins the maximum overlapping Kendall to every ABC enhancer-gene pair and combines it with the ABC (powerlaw) score as exp(log ABC + Kendall x sd(log ABC)/sd(Kendall) x r). Writes EnhancerPredictionsAllPutative_ARC.tsv.gz with the RNA columns attached.",
        {"type": "object", "properties": {
            "abc_predictions": {**_S_STRING, "description": "ABC EnhancerPredictionsAllPutative.tsv.gz."},
            "kendall": {**_S_STRING, "description": "Pairs.Kendall.tsv.gz."},
            "hic": {**_S_BOOLEAN, "default": False, "description": "Use ABC.Score (Hi-C) instead of powerlaw.Score."},
            "cluster": {**_S_STRING, "default": "cluster", "description": "Cluster name."},
            "label": {**_S_STRING, "default": "arc", "description": "Run label (output directory suffix)."}},
         "required": ["abc_predictions", "kendall"]},
        cli=["sce2g-pipeline", "arc"],
        flag_map={"abc_predictions": "--abc-predictions", "kendall": "--kendall", "hic": "--hic", "cluster": "--cluster", "label": "--label"},
        bool_flags={'hic'},
    ),

    _T(
        "sce2g_pipeline_activity_features",
        "Computes the ENCODE_rE2G features the scE2G models read from ABC outputs: numCandidateEnhGene, numTSSEnhGene, numNearbyEnhancers / sumNearbyEnhancers and the gene-class columns, and writes ActivityOnly_features.tsv.gz for the chosen models.",
        {"type": "object", "properties": {
            "abc_predictions": {**_S_STRING, "description": "ABC EnhancerPredictionsAllPutative.tsv.gz."},
            "enhancer_list": {**_S_STRING, "description": "ABC Neighborhoods/EnhancerList.txt."},
            "gene_classes": {**_S_STRING, "description": "Gene promoter-class table (default: setup)."},
            "models": {**_S_ARRAY_S, "description": "Embedded models (multiome_powerlaw_v3, multiome_megamap_v3, scATAC_powerlaw_v3, scATAC_megamap_v3) or model directories; default multiome_powerlaw_v3 scATAC_powerlaw_v3."},
            "cluster": {**_S_STRING, "default": "cluster", "description": "Cluster name."},
            "label": {**_S_STRING, "default": "activity_features", "description": "Run label (output directory suffix)."}},
         "required": ["abc_predictions"]},
        cli=["sce2g-pipeline", "activity-features"],
        flag_map={"abc_predictions": "--abc-predictions", "enhancer_list": "--enhancer-list", "gene_classes": "--gene-classes", "models": "--models", "cluster": "--cluster", "label": "--label"},
        flag_repeat={'models'},
    ),

    _T(
        "sce2g_pipeline_features",
        "Assembles the genome-wide scE2G feature table: the union of the models' feature tables plus scE2G's ARC rows, ARC/Kendall/RNA columns merged by (chr, gene) overlap, renamed to model feature names and NA-filled. Writes feature_table.tsv, external_features_config.tsv and genomewide_features.tsv.gz.",
        {"type": "object", "properties": {
            "activity_features": {**_S_STRING, "description": "ActivityOnly_features.tsv.gz."},
            "arc": {**_S_STRING, "description": "EnhancerPredictionsAllPutative_ARC.tsv.gz."},
            "models": {**_S_ARRAY_S, "description": "Embedded models (multiome_powerlaw_v3, multiome_megamap_v3, scATAC_powerlaw_v3, scATAC_megamap_v3) or model directories; default multiome_powerlaw_v3 scATAC_powerlaw_v3."},
            "cluster": {**_S_STRING, "default": "cluster", "description": "Cluster name."},
            "label": {**_S_STRING, "default": "features", "description": "Run label (output directory suffix)."}},
         "required": ["activity_features"]},
        cli=["sce2g-pipeline", "features"],
        flag_map={"activity_features": "--activity-features", "arc": "--arc", "models": "--models", "cluster": "--cluster", "label": "--label"},
        flag_repeat={'models'},
    ),

    _T(
        "sce2g_pipeline_predict",
        "Applies scE2G models to a genome-wide feature table: E2G.Score from the embedded logistic-regression weights (optionally E2G.Score.cv from held-out-chromosome models), quantile normalisation against the model's reference (E2G.Score.qnorm), the TPM filter (E2G.Score.qnorm.ignoreTPM), thresholded links and optional bedpe per model.",
        {"type": "object", "properties": {
            "features": {**_S_STRING, "description": "genomewide_features.tsv.gz."},
            "models": {**_S_ARRAY_S, "description": "Embedded models (multiome_powerlaw_v3, multiome_megamap_v3, scATAC_powerlaw_v3, scATAC_megamap_v3) or model directories; default multiome_powerlaw_v3 scATAC_powerlaw_v3."},
            "cluster": {**_S_STRING, "default": "cluster", "description": "Cluster name."},
            "crispr_benchmarking": {**_S_BOOLEAN, "default": False, "description": "Also held-out-chromosome cv scores."},
            "no_qnorm": {**_S_BOOLEAN, "default": False, "description": "Skip quantile normalisation."},
            "bedpe": {**_S_BOOLEAN, "default": False, "description": "Also write bedpe links."},
            "label": {**_S_STRING, "default": "predict", "description": "Run label (output directory suffix)."}},
         "required": ["features"]},
        cli=["sce2g-pipeline", "predict"],
        flag_map={"features": "--features", "models": "--models", "cluster": "--cluster", "crispr_benchmarking": "--crispr-benchmarking", "no_qnorm": "--no-qnorm", "bedpe": "--bedpe", "label": "--label"},
        flag_repeat={'models'},
        bool_flags={'no_qnorm', 'bedpe', 'crispr_benchmarking'},
    ),

    _T(
        "sce2g_pipeline_qc",
        "Combines scE2G per-cluster statistics files into all_qc_stats.tsv, flags clusters below 2e6 fragments / 100 cells / 1e6 UMIs, compares each metric with the Sheth-Qiu 2024 reference clusters (percentiles) and draws the upstream QC panels.",
        {"type": "object", "properties": {
            "stats": {**_S_ARRAY_S, "description": "scE2G_predictions_threshold<t>_stats.tsv files."},
            "reference": {**_S_STRING, "description": "Reference QC table (default: setup)."},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip figures."},
            "label": {**_S_STRING, "default": "qc", "description": "Run label (output directory suffix)."}},
         "required": ["stats"]},
        cli=["sce2g-pipeline", "qc"],
        flag_map={"stats": "--stats", "reference": "--reference", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'stats'},
        bool_flags={'no_plots'},
    ),

    _T(
        "sce2g_pipeline_crispr_features",
        "Overlaps CRISPR element-gene pairs with scE2G predictions (apply mode) or genome-wide features (training mode) per gene, aggregates overlapping predictions, recomputes distance to TSS and fills missing values; writes the EPCrisprBenchmark ..._features_NAfilled table used by benchmark and train.",
        {"type": "object", "properties": {
            "predictions": {**_S_STRING, "description": "scE2G_predictions.tsv.gz or genomewide_features.tsv.gz."},
            "crispr": {**_S_STRING, "description": "CRISPR benchmark table (default: setup --crispr)."},
            "models": {**_S_ARRAY_S, "description": "Embedded models (multiome_powerlaw_v3, multiome_megamap_v3, scATAC_powerlaw_v3, scATAC_megamap_v3) or model directories; default multiome_powerlaw_v3 scATAC_powerlaw_v3."},
            "mode": {**_S_STRING, "default": "apply", "description": "apply or training."},
            "model_name": {**_S_STRING, "description": "Model name (output subdirectory)."},
            "cluster": {**_S_STRING, "default": "cluster", "description": "Cluster name."},
            "label": {**_S_STRING, "default": "crispr_features", "description": "Run label (output directory suffix)."}},
         "required": ["predictions"]},
        cli=["sce2g-pipeline", "crispr-features"],
        flag_map={"predictions": "--predictions", "crispr": "--crispr", "models": "--models", "mode": "--mode", "model_name": "--model-name", "cluster": "--cluster", "label": "--label"},
        flag_repeat={'models'},
    ),

    _T(
        "sce2g_pipeline_benchmark",
        "Benchmarks scE2G CRISPR feature tables as upstream benchmark_performance.py does: AUPRC, precision at 50/70% recall and precision/recall at the model threshold with BCa bootstrap 95% CIs, plus the scATAC_ABC, noTPMfilter and distanceToTSS comparison rows. Writes crispr_benchmarking_performance_summary.tsv and an AUPRC figure.",
        {"type": "object", "properties": {
            "crispr_features": {**_S_ARRAY_S, "description": "Paths .../<cluster>/<model>/EPCrispr..._NAfilled.tsv.gz."},
            "model_names": {**_S_ARRAY_S, "description": "Model names (default: parent directory)."},
            "model_thresholds": {**_S_ARRAY_S, "description": "Thresholds (default: embedded model thresholds)."},
            "n_boot": {**_S_INTEGER, "default": 1000, "description": "Bootstrap resamples."},
            "ci_method": {**_S_STRING, "default": "BCa", "description": "BCa, percentile or basic."},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip the figure."},
            "label": {**_S_STRING, "default": "benchmark", "description": "Run label (output directory suffix)."}},
         "required": ["crispr_features"]},
        cli=["sce2g-pipeline", "benchmark"],
        flag_map={"crispr_features": "--crispr-features", "model_names": "--model-names", "model_thresholds": "--model-thresholds", "n_boot": "--n-boot", "ci_method": "--ci-method", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'model_thresholds', 'crispr_features', 'model_names'},
        bool_flags={'no_plots'},
    ),

    _T(
        "sce2g_pipeline_train",
        "Trains a new scE2G-style model: unpenalised logistic regression on log(|x| + 0.01) of CRISPR-labelled features, a full model plus one model per held-out chromosome, with coefficients, per-chromosome log-loss/AUROC/AUPRC and threshold files. Output is a model directory that predict accepts.",
        {"type": "object", "properties": {
            "crispr_features": {**_S_STRING, "description": "NAfilled CRISPR x features table."},
            "features": {**_S_STRING, "description": "genomewide_features.tsv.gz (merged with crispr here)."},
            "crispr": {**_S_STRING, "description": "CRISPR benchmark table."},
            "model": {**_S_STRING, "default": "multiome_powerlaw_v3", "description": "Embedded model whose feature table to use."},
            "feature_table": {**_S_STRING, "description": "Custom feature_table.tsv."},
            "name": {**_S_STRING, "default": "custom_model", "description": "Model name."},
            "polynomial": {**_S_BOOLEAN, "default": False, "description": "Add degree-2 terms."},
            "score_threshold": {**_S_NUMBER, "description": "Score threshold (default: 70% recall of the cv score)."},
            "tpm_threshold": {**_S_NUMBER, "description": "TPM threshold file to write."},
            "label": {**_S_STRING, "default": "train", "description": "Run label (output directory suffix)."}}},
        cli=["sce2g-pipeline", "train"],
        flag_map={"crispr_features": "--crispr-features", "features": "--features", "crispr": "--crispr", "model": "--model", "feature_table": "--feature-table", "name": "--name", "polynomial": "--polynomial", "score_threshold": "--score-threshold", "tpm_threshold": "--tpm-threshold", "label": "--label"},
        bool_flags={'polynomial'},
    ),

    _T(
        "sce2g_pipeline_validate_example",
        "Recomputes the upstream chr22 test fixture's Pairs.Kendall.tsv.gz (Kendall and RNA columns) from its fragments and RNA matrix and reports the agreement (NA pattern, max absolute difference). Needs setup --example --gtf.",
        {"type": "object", "properties": {
            "count_mode": {**_S_STRING, "default": "fragment", "description": "fragment (Signac, default) or insertion."},
            "label": {**_S_STRING, "default": "validate_example", "description": "Run label (output directory suffix)."}}},
        cli=["sce2g-pipeline", "validate-example"],
        flag_map={"count_mode": "--count-mode", "label": "--label"},
    ),

    _T(
        "sce2g_pipeline_validate_release",
        "Recomputes a released IGVF scE2G prediction set from its pairs, elements and genes files (numCandidateEnhGene, numTSSEnhGene, numNearbyEnhancers rebuilt; model applied with the embedded weights, qnorm and TPM filter) and reports how closely Score.ignoreTPM and Score are reproduced.",
        {"type": "object", "properties": {
            "pairs": {**_S_STRING, "description": "All-pairs file (e.g. IGVFFI1706PNVV.tsv.gz)."},
            "elements": {**_S_STRING, "description": "Element file (e.g. IGVFFI3094BAXH.tsv.gz)."},
            "genes": {**_S_STRING, "description": "Gene file (e.g. IGVFFI9905RPTO.tsv.gz)."},
            "model": {**_S_STRING, "default": "multiome_powerlaw_v3", "description": "Model."},
            "chrom": {**_S_STRING, "description": "Restrict to one chromosome."},
            "label": {**_S_STRING, "default": "validate_release", "description": "Run label (output directory suffix)."}},
         "required": ["pairs", "elements", "genes"]},
        cli=["sce2g-pipeline", "validate-release"],
        flag_map={"pairs": "--pairs", "elements": "--elements", "genes": "--genes", "model": "--model", "chrom": "--chrom", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_run",
        "Runs the whole IGVF Perturb-seq CRISPR_Pipeline from count matrices: knee barcode filter + QC, MuData creation (guide metadata, non-targeting buckets, target keys), guide assignment, cis pairs within 1 Mb, SCEPTRE/PerTurbo-style cis + trans inference, additional QC, evaluate_controls, IGV tracks and an HTML dashboard. Inputs are kb-style <batch>_ks_transcripts_out / <batch>_ks_guide_out directories (or .h5ad), the guide design TSV and a GTF. Writes inference_mudata.h5mu (mod/gene, mod/guide with guide_assignment layer, uns cis_/trans_ per_guide/per_element results), TSVs, QC report and dashboard.html.",
        {"type": "object", "properties": {
            "rna": {**_S_ARRAY_S, "description": "RNA count inputs: <batch>_ks_transcripts_out dirs or one .h5ad."},
            "guide": {**_S_ARRAY_S, "description": "Guide count inputs: <batch>_ks_guide_out dirs or one .h5ad."},
            "guide_metadata": {**_S_STRING, "description": "Guide design TSV (guide_id, spacer, targeting, type, guide_chr/start/end, intended_target_name/chr/start/end)."},
            "gtf": {**_S_STRING, "description": "GENCODE GTF(.gz) for gene coordinates and distance-based pairs."},
            "hashing": {**_S_ARRAY_S, "description": "Optional <batch>_ks_hashing_out dirs (cell hashing)."},
            "barcode_filter": {**_S_STRING, "default": "knee2", "description": "none | knee | knee2 (default knee2)."},
            "assignment_method": {**_S_STRING, "default": "sceptre", "description": "sceptre | cleanser | threshold."},
            "pairing_strategy": {**_S_STRING, "default": "default", "description": "default | by_distance | predefined_pairs | all_by_all."},
            "pairs": {**_S_STRING, "description": "pairs_to_test CSV for predefined_pairs."},
            "inference_method": {**_S_STRING, "default": "default", "description": "sceptre | perturbo | sceptre,perturbo | default (cis both + trans PerTurbo)."},
            "engine": {**_S_STRING, "default": "python", "description": "python (method ports) | external (upstream SCEPTRE R / PerTurbo when installed)."},
            "n_resamples": {**_S_INTEGER, "default": 500, "description": "SCEPTRE-style permutation resamples per test."},
            "encode_bed_dir": {**_S_STRING, "description": "Directory with the 5 ENCODE ChIP BED files for the TF benchmark."},
            "moi": {**_S_STRING, "default": "high", "description": "high | low | auto (mean guides per cell > 1.5)."},
            "scrublet": {**_S_BOOLEAN, "description": "Remove doublets (scrublet or fallback) before assignment."},
            "dual_guide": {**_S_BOOLEAN, "description": "Collapse cells with two guides on one element (DUAL_GUIDE)."},
            "no_plots": {**_S_BOOLEAN, "description": "Skip figures."},
            "label": {**_S_STRING, "default": "crispr_run", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["rna", "guide", "guide_metadata"]},
        cli=["crispr-pipeline", "run"],
        flag_map={"rna": "--rna", "guide": "--guide", "guide_metadata": "--guide-metadata", "gtf": "--gtf", "hashing": "--hashing", "barcode_filter": "--barcode-filter", "assignment_method": "--assignment-method", "pairing_strategy": "--pairing-strategy", "pairs": "--pairs", "inference_method": "--inference-method", "engine": "--engine", "n_resamples": "--n-resamples", "encode_bed_dir": "--encode-bed-dir", "moi": "--moi", "scrublet": "--scrublet", "dual_guide": "--dual-guide", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'hashing', 'rna', 'guide'},
        bool_flags={'no_plots', 'dual_guide', 'scrublet'},
    ),

    _T(
        "crispr_pipeline_seqspec",
        "Parses a seqspec YAML and computes the kb count technology string (seqspec index -t kb) per modality, as the pipeline's seqSpecParser does. Writes <modality>_parsed_seqSpec.txt (modality, representation, barcode_whitelist, seqspec_file), the region table and the spacer-tag variant.",
        {"type": "object", "properties": {
            "yaml": {**_S_STRING, "description": "seqspec YAML path."},
            "modalities": {**_S_ARRAY_S, "default": ["rna"], "description": "Modalities to index (rna, guide/crispr, hashing/tag)."},
            "whitelist": {**_S_STRING, "description": "Barcode onlist path recorded in the output."},
            "label": {**_S_STRING, "default": "seqspec", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["yaml"]},
        cli=["crispr-pipeline", "seqspec"],
        flag_map={"yaml": "--yaml", "modalities": "--modalities", "whitelist": "--whitelist", "label": "--label"},
        flag_repeat={'modalities'},
    ),

    _T(
        "crispr_pipeline_seqspec_check",
        "Scans guide or hashtag FASTQs to find where the feature sits in R1/R2 and in which orientation (seqSpecCheck.py score = 3*HitRatio + 2*PosPurity + FlankPurity + 1-Gini). Inputs: R1/R2 FASTQs and the guide metadata; outputs position_table.csv, best_config.csv and a reverse-complement / spacer-tag recommendation.",
        {"type": "object", "properties": {
            "read1": {**_S_ARRAY_S, "description": "R1 FASTQ(s)."},
            "read2": {**_S_ARRAY_S, "description": "R2 FASTQ(s)."},
            "metadata": {**_S_STRING, "description": "Guide metadata (spacer column) or hash table (sequence)."},
            "max_reads": {**_S_INTEGER, "default": 100000, "description": "Reads examined per file."},
            "label": {**_S_STRING, "default": "seqspec_check", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["read1", "read2", "metadata"]},
        cli=["crispr-pipeline", "seqspec-check"],
        flag_map={"read1": "--read1", "read2": "--read2", "metadata": "--metadata", "max_reads": "--max-reads", "label": "--label"},
        flag_repeat={'read2', 'read1'},
    ),

    _T(
        "crispr_pipeline_feature_ref",
        "Builds the guide or hashing feature reference: guide_features.txt / hashing_table.txt (optionally reverse-complemented or spacer-tag prefixed), the kite 1-mismatch FASTA and t2g, and runs kb ref --workflow kite when kb is installed (otherwise prints the command).",
        {"type": "object", "properties": {
            "kind": {**_S_STRING, "default": "guide", "description": "guide | hash."},
            "table": {**_S_STRING, "description": "Guide metadata TSV or hash metadata TSV."},
            "rev_comp": {**_S_BOOLEAN, "description": "Reverse-complement spacers."},
            "spacer": {**_S_STRING, "description": "Spacer tag prepended to guides (switches kb to k=31)."},
            "label": {**_S_STRING, "default": "feature_ref", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["table"]},
        cli=["crispr-pipeline", "feature-ref"],
        flag_map={"kind": "--kind", "table": "--table", "rev_comp": "--rev-comp", "spacer": "--spacer", "label": "--label"},
        bool_flags={'rev_comp'},
    ),

    _T(
        "crispr_pipeline_map",
        "Maps a guide, hashing or RNA library: runs kb count with the seqspec technology (kite / kite:10xFB for features) when kb is installed, or counts guide/hashtag UMIs per cell barcode in pure Python (onlist, 1-mismatch barcode and feature correction, spacer-tag whole-read search). Writes a kb-style <batch>_ks_<modality>_out directory (counts_unfiltered/adata.h5ad, inspect.json, run_info.json).",
        {"type": "object", "properties": {
            "modality": {**_S_STRING, "description": "rna | guide | hash."},
            "fastqs": {**_S_ARRAY_S, "description": "FASTQs in kb order (R1 R2 [R1 R2 ...])."},
            "technology": {**_S_STRING, "description": "kb -x string, e.g. 0,0,16:0,16,28:1,0,0."},
            "parsed_seqspec": {**_S_STRING, "description": "<modality>_parsed_seqSpec.txt from the seqspec tool."},
            "batch": {**_S_STRING, "description": "Batch / measurement set name (output dir prefix)."},
            "features": {**_S_STRING, "description": "Guide metadata or hash table (Python engine)."},
            "barcodes": {**_S_STRING, "description": "Cell barcode onlist."},
            "spacer": {**_S_STRING, "description": "Spacer tag before the guide."},
            "is_10x3v3": {**_S_BOOLEAN, "description": "Use kb's 10XV3 + kite:10xFB instead of the seqspec string."},
            "index": {**_S_STRING, "description": "kb index (kb engine)."},
            "t2g": {**_S_STRING, "description": "kb t2g (kb engine)."},
            "engine": {**_S_STRING, "default": "auto", "description": "auto | kb | python."},
            "label": {**_S_STRING, "default": "map", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["modality", "fastqs", "batch"]},
        cli=["crispr-pipeline", "map"],
        flag_map={"modality": "--modality", "fastqs": "--fastqs", "technology": "--technology", "parsed_seqspec": "--parsed-seqspec", "batch": "--batch", "features": "--features", "barcodes": "--barcodes", "spacer": "--spacer", "is_10x3v3": "--is-10x3v3", "index": "--index", "t2g": "--t2g", "engine": "--engine", "label": "--label"},
        flag_repeat={'fastqs'},
        bool_flags={'is_10x3v3'},
    ),

    _T(
        "crispr_pipeline_concat",
        "Concatenates per-batch kb outputs (<batch>_ks_*_out) into one AnnData, adding the batch column from the directory name and joining covariates (anndata_concat.py). Writes concatenated_adata.h5ad.",
        {"type": "object", "properties": {
            "inputs": {**_S_ARRAY_S, "description": "<batch>_ks_*_out directories."},
            "covariates": {**_S_STRING, "description": "parse_covariate.csv (optional)."},
            "label": {**_S_STRING, "default": "concat", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["inputs"]},
        cli=["crispr-pipeline", "concat"],
        flag_map={"inputs": "--inputs", "covariates": "--covariates", "label": "--label"},
        flag_repeat={'inputs'},
    ),

    _T(
        "crispr_pipeline_preprocess",
        "Preprocesses concatenated scRNA counts as the pipeline does: knee / knee2 barcode filter on the barcode-rank curve (or min genes), scanpy QC metrics, genes in >= 10 cells, pct_counts_mt < QC_pct_mito. Writes filtered_anndata.h5ad and the knee / violin / scatter figures.",
        {"type": "object", "properties": {
            "adata": {**_S_STRING, "description": "Concatenated RNA .h5ad."},
            "gene_names": {**_S_STRING, "description": "cells_x_genes.genes.names.txt or a ks_transcripts_out directory."},
            "min_genes": {**_S_INTEGER, "default": 500, "description": "Min genes per cell when barcode_filter is none."},
            "pct_mito": {**_S_NUMBER, "default": 15.0, "description": "Max mitochondrial percent."},
            "reference": {**_S_STRING, "default": "human", "description": "human | mouse."},
            "barcode_filter": {**_S_STRING, "default": "knee2", "description": "none | knee | knee2."},
            "label": {**_S_STRING, "default": "preprocess", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["adata"]},
        cli=["crispr-pipeline", "preprocess"],
        flag_map={"adata": "--adata", "gene_names": "--gene-names", "min_genes": "--min-genes", "pct_mito": "--pct-mito", "reference": "--reference", "barcode_filter": "--barcode-filter", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_create_mudata",
        "Creates mudata.h5mu from filtered RNA and guide counts (create_mdata.py): merges the guide design, fills type, buckets controls into non-targeting|k, builds intended_target_key, sets moi / capture_method, adds GTF gene coordinates and intersects barcodes (optionally with hashing).",
        {"type": "object", "properties": {
            "rna": {**_S_STRING, "description": "filtered_anndata.h5ad."},
            "guide": {**_S_STRING, "description": "Concatenated guide .h5ad."},
            "guide_metadata": {**_S_STRING, "description": "Guide design TSV."},
            "gtf": {**_S_STRING, "description": "GTF(.gz)."},
            "hashing": {**_S_STRING, "description": "Demultiplexed hashing .h5ad (optional)."},
            "moi": {**_S_STRING, "default": "high", "description": "high | low | auto."},
            "capture_method": {**_S_STRING, "default": "CROP-seq", "description": "CROP-seq | direct capture."},
            "label": {**_S_STRING, "default": "create_mudata", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["rna", "guide", "guide_metadata"]},
        cli=["crispr-pipeline", "create-mudata"],
        flag_map={"rna": "--rna", "guide": "--guide", "guide_metadata": "--guide-metadata", "gtf": "--gtf", "hashing": "--hashing", "moi": "--moi", "capture_method": "--capture-method", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_hashing",
        "Demultiplexes cell-hashing counts: keeps barcodes in the filtered RNA, runs GMM-demux per batch when installed or a per-HTO Gaussian-mixture fallback, labels hto_type / hto_type_split and drops negatives and multiplets. Writes concatenated_hashing_demux.h5ad and the unfiltered version.",
        {"type": "object", "properties": {
            "hashing": {**_S_STRING, "description": "Hashing counts .h5ad (obs batch)."},
            "rna": {**_S_STRING, "description": "Filtered RNA .h5ad."},
            "label": {**_S_STRING, "default": "hashing", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["hashing", "rna"]},
        cli=["crispr-pipeline", "hashing"],
        flag_map={"hashing": "--hashing", "rna": "--rna", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_doublets",
        "Removes doublets from a MuData with scrublet (or a scrublet-method fallback: simulated doublets, PCA, kNN doublet score, automatic threshold) and re-intersects modalities. Writes mdata_doublets.h5mu and the score histogram.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "mudata.h5mu."},
            "label": {**_S_STRING, "default": "doublets", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mudata"]},
        cli=["crispr-pipeline", "doublets"],
        flag_map={"mudata": "--mudata", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_assign_guides",
        "Assigns guides to cells per batch with the SCEPTRE Poisson-mixture method (EM restarts, posterior >= 0.8), a CLEANSER-style ambient/NB mixture or a UMI threshold, then keeps genes expressed in more than 5% of cells and optionally collapses dual guides. engine external runs the cleanser binary or the upstream SCEPTRE R script when installed. Writes concat_mudata.h5mu with the guide_assignment layer.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "mudata.h5mu."},
            "method": {**_S_STRING, "default": "sceptre", "description": "sceptre | cleanser | threshold."},
            "threshold": {**_S_NUMBER, "description": "Posterior threshold (sceptre 0.8, cleanser 1.0)."},
            "n_em_rep": {**_S_INTEGER, "default": 5, "description": "EM restarts (SCEPTRE mixture)."},
            "umi_threshold": {**_S_INTEGER, "default": 5, "description": "UMI threshold for the threshold method."},
            "min_cells_fraction": {**_S_NUMBER, "default": 0.05, "description": "Gene filter fraction (QC_min_cells_per_gene)."},
            "dual_guide": {**_S_BOOLEAN, "description": "Collapse dual-guide cells."},
            "label": {**_S_STRING, "default": "assign_guides", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."},
            "engine": {**_S_STRING, "default": "python", "description": "python (method ports) | external (cleanser binary / upstream assign_grnas_sceptre.R; falls back to python when absent)."},
            "upstream_bin": {**_S_STRING, "description": "Upstream bin/ directory with assign_grnas_sceptre.R (engine external)."}},
         "required": ["mudata"]},
        cli=["crispr-pipeline", "assign-guides"],
        flag_map={"mudata": "--mudata", "method": "--method", "threshold": "--threshold", "n_em_rep": "--n-em-rep", "umi_threshold": "--umi-threshold", "min_cells_fraction": "--min-cells-fraction", "dual_guide": "--dual-guide", "label": "--label", "engine": "--engine", "upstream_bin": "--upstream-bin"},
        bool_flags={'dual_guide'},
    ),

    _T(
        "crispr_pipeline_pairs",
        "Builds pairs_to_test (each guide x GTF genes within 1 Mb of the guide, user pairs, or all genes) and prepares the inference MuData (target metadata, gene_id keys, cis subset of tested genes, tested + control guides and targeted cells for the default strategy). Writes pairs_to_test.csv and mudata_inference_input.h5mu.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "concat_mudata.h5mu."},
            "strategy": {**_S_STRING, "default": "default", "description": "default | by_distance | predefined_pairs | all."},
            "gtf": {**_S_STRING, "description": "GTF(.gz)."},
            "pairs": {**_S_STRING, "description": "User pairs CSV (guide_id, gene_name)."},
            "limit": {**_S_INTEGER, "default": 1000000, "description": "Max guide-to-gene-start distance (bp); -1 for all genes."},
            "label": {**_S_STRING, "default": "pairs", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mudata"]},
        cli=["crispr-pipeline", "pairs"],
        flag_map={"mudata": "--mudata", "strategy": "--strategy", "gtf": "--gtf", "pairs": "--pairs", "limit": "--limit", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_inference",
        "Tests perturbation effects per guide and per element: the SCEPTRE slot (NB score statistic calibrated by label permutations with a skew-normal fit, complement control) and the PerTurbo slot (NB maximum-likelihood fold change with Wald p-value and log2_fc_std), or the upstream R/torch packages with engine external. Method default = cis SCEPTRE + PerTurbo on pairs and trans PerTurbo on all pairs. Writes inference_mudata.h5mu, per-table TSVs and per_guide_output / per_element_output.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "mudata_inference_input.h5mu (pairs_to_test in uns)."},
            "trans_mudata": {**_S_STRING, "description": "Full concat_mudata.h5mu for trans (default method)."},
            "method": {**_S_STRING, "default": "default", "description": "sceptre | perturbo | sceptre,perturbo | default."},
            "engine": {**_S_STRING, "default": "python", "description": "python | external."},
            "upstream_bin": {**_S_STRING, "description": "Upstream bin/ directory for the external engine."},
            "n_resamples": {**_S_INTEGER, "default": 500, "description": "Permutation resamples."},
            "side": {**_S_STRING, "default": "both", "description": "both | left | right."},
            "label": {**_S_STRING, "default": "inference", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mudata"]},
        cli=["crispr-pipeline", "inference"],
        flag_map={"mudata": "--mudata", "trans_mudata": "--trans-mudata", "method": "--method", "engine": "--engine", "upstream_bin": "--upstream-bin", "n_resamples": "--n-resamples", "side": "--side", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_merge_results",
        "Merges pipeline result tables into a MuData: SCEPTRE + PerTurbo (outer merges, coordinate-aware element keys), cis + trans uns keys, per-chunk SCEPTRE outputs, legacy test_results, or exports per_guide_output / per_element_output with average target expression.",
        {"type": "object", "properties": {
            "mode": {**_S_STRING, "description": "methods | cis-trans | chunks | legacy | add-test-results | export."},
            "base_mudata": {**_S_STRING, "description": "Base MuData."},
            "sceptre_per_guide": {**_S_STRING, "description": "SCEPTRE per-guide TSV."},
            "sceptre_per_element": {**_S_STRING, "description": "SCEPTRE per-element TSV."},
            "perturbo_per_guide": {**_S_STRING, "description": "PerTurbo per-guide TSV."},
            "perturbo_per_element": {**_S_STRING, "description": "PerTurbo per-element TSV."},
            "cis_per_guide": {**_S_STRING, "description": "Cis per-guide TSV."},
            "cis_per_element": {**_S_STRING, "description": "Cis per-element TSV."},
            "trans_per_guide": {**_S_STRING, "description": "Trans per-guide TSV."},
            "trans_per_element": {**_S_STRING, "description": "Trans per-element TSV."},
            "per_guide_files": {**_S_ARRAY_S, "description": "Chunk per-guide TSVs."},
            "per_element_files": {**_S_ARRAY_S, "description": "Chunk per-element TSVs."},
            "label": {**_S_STRING, "default": "merge_results", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mode"]},
        cli=["crispr-pipeline", "merge-results"],
        flag_map={"mode": "--mode", "base_mudata": "--base-mudata", "sceptre_per_guide": "--sceptre-per-guide", "sceptre_per_element": "--sceptre-per-element", "perturbo_per_guide": "--perturbo-per-guide", "perturbo_per_element": "--perturbo-per-element", "cis_per_guide": "--cis-per-guide", "cis_per_element": "--cis-per-element", "trans_per_guide": "--trans-per-guide", "trans_per_element": "--trans-per-element", "per_guide_files": "--per-guide-files", "per_element_files": "--per-element-files", "label": "--label"},
        flag_repeat={'per_guide_files', 'per_element_files'},
    ),

    _T(
        "crispr_pipeline_qc",
        "Computes the pipeline's additional QC on an inference MuData: gene and guide metric tables overall and per batch, intended-target knockdown (strong = log2FC <= log2(0.4), p < 0.05) with AUROC/AUPRC against non-targeting guides, and trans QC after BH correction (per-guide significant counts, optional validated links). Writes TSVs, plots and qc_report.html.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "inference_mudata.h5mu."},
            "results_key": {**_S_STRING, "default": "auto", "description": "uns key or auto (trans_per_guide_results first)."},
            "fc_threshold": {**_S_NUMBER, "default": 0.4, "description": "Fold-change threshold for strong knockdown."},
            "pval_threshold": {**_S_NUMBER, "default": 0.05, "description": "Significance threshold."},
            "fdr_method": {**_S_STRING, "default": "fdr_bh", "description": "fdr_bh | bonferroni | none."},
            "validated_links": {**_S_STRING, "description": "TSV of validated trans links (guide_target, gene)."},
            "label": {**_S_STRING, "default": "qc", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mudata"]},
        cli=["crispr-pipeline", "qc"],
        flag_map={"mudata": "--mudata", "results_key": "--results-key", "fc_threshold": "--fc-threshold", "pval_threshold": "--pval-threshold", "fdr_method": "--fdr-method", "validated_links": "--validated-links", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_evaluate",
        "Evaluates an inference MuData as evaluate_controls.py does (direct-target tests vs non-targeting guides on the same genes: volcano, PR/ROC, AUROC/AUPRC) and writes IGV bedpe/bedgraph tracks, volcano plots and central-node network plots for the cis and trans per-element results.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "inference_mudata.h5mu."},
            "gtf": {**_S_STRING, "description": "GTF(.gz) for gene names."},
            "num_nodes": {**_S_INTEGER, "default": 1, "description": "Central nodes per network plot."},
            "central_nodes": {**_S_ARRAY_S, "description": "Explicit central nodes (intended_target_name)."},
            "label": {**_S_STRING, "default": "evaluate", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mudata"]},
        cli=["crispr-pipeline", "evaluate"],
        flag_map={"mudata": "--mudata", "gtf": "--gtf", "num_nodes": "--num-nodes", "central_nodes": "--central-nodes", "label": "--label"},
        flag_repeat={'central_nodes'},
    ),

    _T(
        "crispr_pipeline_expression",
        "Returns per-cell expression of one gene in cells carrying a guide versus cells with only non-targeting guides (utlis.get_perturbation_expression). Writes the table and a violin plot.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "MuData with guide_assignment."},
            "guide_id": {**_S_STRING, "description": "Guide id."},
            "gene_id": {**_S_STRING, "description": "Gene id (gene.var_names)."},
            "label": {**_S_STRING, "default": "expression", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mudata", "guide_id", "gene_id"]},
        cli=["crispr-pipeline", "expression"],
        flag_map={"mudata": "--mudata", "guide_id": "--guide-id", "gene_id": "--gene-id", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_tf_benchmark",
        "Benchmarks trans results against ENCODE ChIP-seq: for each TF element, one-sided Fisher test of significant trans genes (p <= 0.001) vs genes whose TSS +/- w bp overlaps a peak, w in 500..10000. Inputs: inference MuData, GTF and the 5 bundled ENCODE BED files or TF=peaks pairs. Writes enrichment tables and tf_benchmark.png.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "inference_mudata.h5mu with trans_per_element_results."},
            "gtf": {**_S_STRING, "description": "GTF(.gz)."},
            "encode_bed_dir": {**_S_STRING, "description": "Directory with the ENCODE BED files (setup)."},
            "peaks": {**_S_ARRAY_S, "description": "TF_ID=peaks.bed(.gz) overrides."},
            "windows": {"type": "array", "items": {"type": "integer"}, "default": [500, 1000, 2500, 5000, 10000], "description": "Promoter half-widths (bp)."},
            "cutoff": {**_S_NUMBER, "default": 0.001, "description": "Significance cutoff on p_value."},
            "upstream_compat": {**_S_BOOLEAN, "description": "Reproduce upstream END-as-TSS behaviour."},
            "label": {**_S_STRING, "default": "tf_benchmark", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mudata", "gtf"]},
        cli=["crispr-pipeline", "tf-benchmark"],
        flag_map={"mudata": "--mudata", "gtf": "--gtf", "encode_bed_dir": "--encode-bed-dir", "peaks": "--peaks", "windows": "--windows", "cutoff": "--cutoff", "upstream_compat": "--upstream-compat", "label": "--label"},
        flag_repeat={'windows', 'peaks'},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "crispr_pipeline_dashboard",
        "Builds the pipeline dashboard: UMI / detected-gene / guide-UMI threshold bar plots, guides-per-cell and cells-per-guide histograms, HTO plots, kb mapping summaries, the barcode filtering flow and top inference tables, as one self-contained dashboard.html.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "inference_mudata.h5mu."},
            "ks_dirs": {**_S_ARRAY_S, "description": "<batch>_ks_*_out directories for mapping JSON."},
            "gene_ann": {**_S_STRING, "description": "Unfiltered concatenated RNA .h5ad."},
            "gene_ann_filtered": {**_S_STRING, "description": "filtered_anndata.h5ad."},
            "guide_ann": {**_S_STRING, "description": "Concatenated guide .h5ad."},
            "hashing_unfiltered_demux": {**_S_STRING, "description": "Unfiltered hashing demux .h5ad."},
            "image_dirs": {**_S_ARRAY_S, "description": "Extra directories whose PNGs are embedded (QC, evaluation, benchmark)."},
            "label": {**_S_STRING, "default": "dashboard", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mudata"]},
        cli=["crispr-pipeline", "dashboard"],
        flag_map={"mudata": "--mudata", "ks_dirs": "--ks-dirs", "gene_ann": "--gene-ann", "gene_ann_filtered": "--gene-ann-filtered", "guide_ann": "--guide-ann", "hashing_unfiltered_demux": "--hashing-unfiltered-demux", "image_dirs": "--image-dirs", "label": "--label"},
        flag_repeat={'ks_dirs', 'image_dirs'},
    ),

    _T(
        "crispr_pipeline_samplesheet",
        "Validates a CRISPR_Pipeline samplesheet (R1/R2, file_modality, measurement_sets, seqspec, onlist, guide design, hashtag map): groups runs per modality and measurement set, writes parse_covariate.csv, the covariate formula and FASTQ batches; with from_portal maps a portal per-sample TSV first.",
        {"type": "object", "properties": {
            "samplesheet": {**_S_STRING, "description": "Samplesheet CSV/TSV."},
            "from_portal": {**_S_BOOLEAN, "description": "Input is a portal per-sample TSV."},
            "label": {**_S_STRING, "default": "samplesheet", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["samplesheet"]},
        cli=["crispr-pipeline", "samplesheet"],
        flag_map={"samplesheet": "--samplesheet", "from_portal": "--from-portal", "label": "--label"},
        bool_flags={'from_portal'},
    ),

    _T(
        "crispr_pipeline_portal_samplesheet",
        "Builds the per-sample TSV for an IGVF analysis set from the portal (measurement and auxiliary sets, validated seqspecs, onlists, guide library, hashtag map) and the pipeline samplesheet.",
        {"type": "object", "properties": {
            "accession": {**_S_STRING, "description": "Analysis set accession (IGVFDS...)."},
            "rna_seqspec": {**_S_STRING, "description": "Fallback RNA seqspec."},
            "sgrna_seqspec": {**_S_STRING, "description": "Fallback sgRNA seqspec."},
            "hash_seqspec": {**_S_STRING, "description": "Fallback hashing seqspec."},
            "label": {**_S_STRING, "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["accession"]},
        cli=["crispr-pipeline", "portal-samplesheet"],
        flag_map={"accession": "--accession", "rna_seqspec": "--rna-seqspec", "sgrna_seqspec": "--sgrna-seqspec", "hash_seqspec": "--hash-seqspec", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_portal_download",
        "Plans (and with execute, downloads and md5-verifies) the FASTQ, seqspec, onlist, guide design and hashtag map files of a portal per-sample TSV. Writes per_sample_local.tsv.",
        {"type": "object", "properties": {
            "sample": {**_S_STRING, "description": "Per-sample TSV."},
            "file_types": {**_S_STRING, "default": "all", "description": "fastq | other | all."},
            "output_dir": {**_S_STRING, "description": "Download directory."},
            "execute": {**_S_BOOLEAN, "description": "Actually download (large FASTQs)."},
            "label": {**_S_STRING, "default": "portal_download", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["sample"]},
        cli=["crispr-pipeline", "portal-download"],
        flag_map={"sample": "--sample", "file_types": "--file-types", "output_dir": "--output-dir", "execute": "--execute", "label": "--label"},
        bool_flags={'execute'},
    ),

    _T(
        "crispr_pipeline_chunk",
        "Splits an inference MuData by genes for SCEPTRE (auto/off/force by n_cells x n_genes) or PerTurbo (balanced ranges) and writes the chunk files and chunk_manifest.tsv.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "Inference input MuData."},
            "mode": {**_S_STRING, "default": "sceptre", "description": "sceptre | perturbo."},
            "chunk_size": {**_S_INTEGER, "default": 1000, "description": "Genes per chunk."},
            "chunk_mode": {**_S_STRING, "default": "auto", "description": "auto | off | force."},
            "label": {**_S_STRING, "default": "chunk", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["mudata"]},
        cli=["crispr-pipeline", "chunk"],
        flag_map={"mudata": "--mudata", "mode": "--mode", "chunk_size": "--chunk-size", "chunk_mode": "--chunk-mode", "label": "--label"},
    ),

    _T(
        "crispr_pipeline_nextflow",
        "Writes (and with execute, runs) the exact upstream Nextflow command for pinellolab/CRISPR_Pipeline at the pinned commit with the given samplesheet, profile and params.",
        {"type": "object", "properties": {
            "input": {**_S_STRING, "description": "Samplesheet."},
            "outdir": {**_S_STRING, "description": "Output directory."},
            "profile": {**_S_STRING, "default": "local", "description": "local | slurm | google."},
            "param": {**_S_ARRAY_S, "description": "NAME=VALUE pipeline params."},
            "execute": {**_S_BOOLEAN, "description": "Run nextflow if installed."},
            "label": {**_S_STRING, "default": "nextflow", "description": "Run label used in the output directory name (Docs/CRISPRPipeline/<timestamp>_<label>)."}},
         "required": ["input", "outdir"]},
        cli=["crispr-pipeline", "nextflow"],
        flag_map={"input": "--input", "outdir": "--outdir", "profile": "--profile", "param": "--param", "execute": "--execute", "label": "--label"},
        flag_repeat={'param'},
        bool_flags={'execute'},
    ),

    _T(
        "crispr_pipeline_setup",
        "Downloads the 5 ENCODE ChIP-seq BED files used by the TF benchmark and the example seqspec / guide design files from the pinned upstream commit into Data/CRISPRPipeline (under 1 MB).",
        {"type": "object", "properties": {
            "dest": {**_S_STRING, "description": "Destination directory."},
            "force": {**_S_BOOLEAN, "description": "Re-download."}}},
        cli=["crispr-pipeline", "setup"],
        flag_map={"dest": "--dest", "force": "--force"},
        bool_flags={'force'},
    ),

    _T(
        "crispr_jamboree2_inspect",
        "Summarises a jamboree-format Perturb-seq MuData (modalities, shapes, obs/var columns, layers, MOI, capture method, pairs_to_test and test_results with pair_type counts) and checks it against the jamboree input/output specification (eval_setup.ipynb). Input: .h5mu. Output: report.md + summary.json with the schema check.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu): mod/gene counts, mod/guide counts + layers['guide_assignment'], guide.var targeting / intended_target_name, uns pairs_to_test"},
            "strict": {**_S_BOOLEAN, "default": False, "description": "exit non-zero when a required field is missing"},
            "label": {**_S_STRING, "default": "inspect", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree2", "inspect"],
        flag_map={"mudata": "--mudata", "strict": "--strict", "label": "--label"},
        bool_flags={'strict'},
    ),

    _T(
        "crispr_jamboree2_guide_reference",
        "Builds the STARsolo guide reference of CountGuides.ipynb / parsing_guide_medata.py: one pseudo-chromosome per protospacer (chr_<sgRNA_ID>_peseudo_chr_<seq>), BED12, isoforms.txt, a bed2gtf-style GTF, pseudo_genome.fa and a protospacer FASTA. Runs STAR genomeGenerate when STAR is on PATH, otherwise prints the command. Input: guide metadata (xlsx/tsv/csv with sgRNA_ID + sgRNA_sequences or guide_id + protospacer).",
        {"type": "object", "properties": {
            "guides": {**_S_STRING, "description": "guide metadata table"},
            "nnn_revcomp": {**_S_BOOLEAN, "default": False, "description": "use 'NNN' + reverse complement of each protospacer"},
            "sa_index_nbases": {**_S_INTEGER, "default": 5, "description": "STAR --genomeSAindexNbases"},
            "no_star": {**_S_BOOLEAN, "default": False, "description": "do not run STAR even when installed"},
            "label": {**_S_STRING, "default": "guide_reference", "description": "run label for the output directory"}},
         "required": ["guides"]},
        cli=["crispr-jamboree2", "guide-reference"],
        flag_map={"guides": "--guides", "nnn_revcomp": "--nnn-revcomp", "sa_index_nbases": "--sa-index-nbases", "no_star": "--no-star", "label": "--label"},
        bool_flags={'no_star', 'nnn_revcomp'},
    ),

    _T(
        "crispr_jamboree2_count_guides",
        "Counts gRNA UMIs per cell from paired FASTQs in pure Python with the starSoloGuide.nf layout (CB = R1 1-16, UMI = R1 17-28): whitelist 1-mismatch barcode correction, protospacer (or its reverse complement) found in R2, UMI de-duplication with 1-mismatch collapsing. Alternatively parses an existing STARsolo Solo.out directory (soloOutGuideParsing.py). Output: guides_star_solo.h5ad (cells x guides), MatrixMarket files and per-guide UMI totals.",
        {"type": "object", "properties": {
            "guides": {**_S_STRING, "description": "guide metadata table"},
            "r1": {**_S_STRING, "description": "R1 FASTQ(.gz)"},
            "r2": {**_S_STRING, "description": "R2 FASTQ(.gz)"},
            "whitelist": {**_S_STRING, "description": "cell barcode whitelist"},
            "solo_out": {**_S_STRING, "description": "existing STARsolo Solo.out directory to parse instead"},
            "cb_start": {**_S_INTEGER, "default": 1, "description": "1-based CB start"},
            "cb_len": {**_S_INTEGER, "default": 16, "description": "CB length"},
            "umi_start": {**_S_INTEGER, "default": 17, "description": "1-based UMI start"},
            "umi_len": {**_S_INTEGER, "default": 12, "description": "UMI length"},
            "max_reads": {**_S_INTEGER, "description": "stop after this many read pairs"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "write guides x cells, as soloOutGuideParsing.py"},
            "label": {**_S_STRING, "default": "count_guides", "description": "run label for the output directory"}},
         "required": ["guides"]},
        cli=["crispr-jamboree2", "count-guides"],
        flag_map={"guides": "--guides", "r1": "--r1", "r2": "--r2", "whitelist": "--whitelist", "solo_out": "--solo-out", "cb_start": "--cb-start", "cb_len": "--cb-len", "umi_start": "--umi-start", "umi_len": "--umi-len", "max_reads": "--max-reads", "upstream_compat": "--upstream-compat", "label": "--label"},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "crispr_jamboree2_seqspec_index",
        "Computes `seqspec index` for a seqspec YAML (seqSpecExtract.ipynb): the positions of barcode, UMI and cDNA/guide regions on each read, as the kb `-x` technology string and STARsolo --soloCB*/--soloUMI* flags, plus the barcode whitelist filename. Negative-strand reads are walked from their primer backwards. Output: regions table + report.",
        {"type": "object", "properties": {
            "spec": {**_S_STRING, "description": "seqspec YAML"},
            "modality": {**_S_STRING, "description": "modality, e.g. crispr / guide / rna"},
            "reads": {**_S_STRING, "description": "comma-separated read ids, e.g. R1,R2"},
            "tool": {**_S_STRING, "default": "kb", "description": "kb or starsolo"},
            "label": {**_S_STRING, "default": "seqspec", "description": "run label for the output directory"}},
         "required": ["spec", "modality", "reads"]},
        cli=["crispr-jamboree2", "seqspec-index"],
        flag_map={"spec": "--spec", "modality": "--modality", "reads": "--reads", "tool": "--tool", "label": "--label"},
    ),

    _T(
        "crispr_jamboree2_assign_guides",
        "Writes the guide_assignment layer of a MuData (igvf_guide_assignment.py): UMI threshold (count >= t, default 5), CLEANSER (zero-truncated Poisson+NB for CROP-seq or NB+NB for direct capture, Stan priors, MAP fit per gRNA; posterior PZi, or 1 where PZi >= t), or the sceptre two-component Poisson-GLM mixture (posterior > 0.8). Output: .h5mu with the layer, cells-per-guide table, CLEANSER parameters.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu): mod/gene counts, mod/guide counts + layers['guide_assignment'], guide.var targeting / intended_target_name, uns pairs_to_test"},
            "method": {**_S_STRING, "default": "umi-threshold", "description": "umi-threshold, cleanser or sceptre"},
            "threshold": {**_S_NUMBER, "description": "UMI threshold, or CLEANSER posterior threshold (omit to store posteriors)"},
            "capture_method": {**_S_STRING, "description": "override guide.uns capture_method: CROP-seq or direct capture"},
            "out": {**_S_STRING, "description": "output .h5mu path"},
            "label": {**_S_STRING, "default": "assign_guides", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree2", "assign-guides"],
        flag_map={"mudata": "--mudata", "method": "--method", "threshold": "--threshold", "capture_method": "--capture-method", "out": "--out", "label": "--label"},
    ),

    _T(
        "crispr_jamboree2_infer",
        "Tests every pairs_to_test (element, gene) pair with the jamboree-2 inference modules: R wilcox.test (low MOI: NT-cell controls; high MOI: complement), MASS glm.nb with a log(total UMIs) offset, scanpy rank_genes_groups (wilcoxon, t-test, t-test_overestim_var on log1p counts), a sceptre re-implementation (NB score test + conditional randomisation / permutations, skew-normal p) and PerTurbo when installed. Output: test_results.<method>.tsv per method (upstream columns p_value / log2_fc) and optionally an .h5mu with uns test_results.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu): mod/gene counts, mod/guide counts + layers['guide_assignment'], guide.var targeting / intended_target_name, uns pairs_to_test"},
            "methods": {**_S_ARRAY_S, "default": ["wilcoxon", "negbinom", "scanpy-wilcoxon", "sceptre"], "description": "wilcoxon, negbinom, scanpy-wilcoxon, scanpy-t-test, scanpy-t-test-overestim-var, sceptre, perturbo"},
            "side": {**_S_STRING, "default": "both", "description": "left, right or both"},
            "moi": {**_S_STRING, "description": "override MOI: low or high"},
            "resamples": {**_S_INTEGER, "default": 499, "description": "sceptre resamples"},
            "seed": {**_S_INTEGER, "default": 4, "description": "random seed"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "negbinom writes `l2fc` = exp(b1) as upstream"},
            "out": {**_S_STRING, "description": "output .h5mu (first method)"},
            "label": {**_S_STRING, "default": "infer", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree2", "infer"],
        flag_map={"mudata": "--mudata", "methods": "--methods", "side": "--side", "moi": "--moi", "resamples": "--resamples", "seed": "--seed", "upstream_compat": "--upstream-compat", "out": "--out", "label": "--label"},
        flag_repeat={'methods'},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "crispr_jamboree2_simulate",
        "Simulates a Perturb-seq screen from a real MuData (simulation_engine R scripts): DESeq2-style poscounts size factors and parametric-trend MAP dispersions per gene, then NB counts with planted effect sizes for chosen (gene, element) pairs and guide-to-guide variability. Optional null pairs (effect_size 1) make the result usable by evaluate. Output: simulation .h5mu (pairs_to_test with effect_size), gene_dispersions.tsv.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu): mod/gene counts, mod/guide counts + layers['guide_assignment'], guide.var targeting / intended_target_name, uns pairs_to_test"},
            "pairs": {**_S_STRING, "description": "TSV gene_id, intended_target_name[, effect_size]"},
            "perturb": {**_S_ARRAY_S, "description": "gene_id:intended_target_name:effect_size entries"},
            "effect_size": {**_S_NUMBER, "default": 0.5, "description": "default effect size (fold change)"},
            "n_null_pairs": {**_S_INTEGER, "default": 0, "description": "random null pairs to add"},
            "guide_var": {**_S_NUMBER, "default": 0.1, "description": "SD of guide-level effect sizes"},
            "seed": {**_S_INTEGER, "default": 1, "description": "random seed"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "write the unsimulated counts, reproducing the upstream writeH5MU bug"},
            "out": {**_S_STRING, "description": "output .h5mu"},
            "label": {**_S_STRING, "default": "simulate", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree2", "simulate"],
        flag_map={"mudata": "--mudata", "pairs": "--pairs", "perturb": "--perturb", "effect_size": "--effect-size", "n_null_pairs": "--n-null-pairs", "guide_var": "--guide-var", "seed": "--seed", "upstream_compat": "--upstream-compat", "out": "--out", "label": "--label"},
        flag_repeat={'perturb'},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "crispr_jamboree2_evaluate",
        "Scores inference results as eval_and_clustergram.ipynb: truth = effect_size != 1 (simulations) or positive vs negative control pair_type; AUPRC (trapezoid auc of the PR curve) and AUROC on -log10(p), Pearson r and MSE of log2 effect sizes, plus the gRNA x gene log2 mean-ratio clustergram (average linkage). Inputs: an .h5mu with uns test_results and/or result tables. Output: metrics, curves, clustergram TSV/PNG.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "MuData with uns test_results (also used for the clustergram)"},
            "results": {**_S_ARRAY_S, "description": "test_results TSV/CSV files"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "score by the raw p-value as the notebook"},
            "no_clustergram": {**_S_BOOLEAN, "default": False, "description": "skip the clustergram"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "evaluate", "description": "run label for the output directory"}}},
        cli=["crispr-jamboree2", "evaluate"],
        flag_map={"mudata": "--mudata", "results": "--results", "upstream_compat": "--upstream-compat", "no_clustergram": "--no-clustergram", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'results'},
        bool_flags={'upstream_compat', 'no_clustergram', 'no_plots'},
    ),

    _T(
        "crispr_jamboree2_volcano",
        "Volcano and IGV export of volcano_and_igv.ipynb: down / up-regulated pairs at |log2_fc| >= threshold and p (or posterior probability) <= threshold, a volcano figure, an IGV bedgraph for promoter pairs (element == gene) and a bedpe of element -> gene links with coordinates from gene.var / guide.var. Output: volcano_down/up.tsv, igv_links_bedpe.tsv, igv_promoter_bedgraph.tsv.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu): mod/gene counts, mod/guide counts + layers['guide_assignment'], guide.var targeting / intended_target_name, uns pairs_to_test"},
            "results": {**_S_STRING, "description": "test_results table (default: mudata uns test_results)"},
            "log2fc_threshold": {**_S_NUMBER, "default": 2.0, "description": "|log2 FC| threshold"},
            "sig_threshold": {**_S_NUMBER, "default": 0.01, "description": "p / posterior threshold"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip the figure"},
            "label": {**_S_STRING, "default": "volcano", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree2", "volcano"],
        flag_map={"mudata": "--mudata", "results": "--results", "log2fc_threshold": "--log2fc-threshold", "sig_threshold": "--sig-threshold", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "crispr_jamboree2_network",
        "Element -> gene regulatory network of network.ipynb / evaluation.py: edges with |log2_fc| >= min_weight around each central element (default: the most significant ones), node size -log10 p, circular-layout figure, plus mean_nontargeting_expression per gene (cells with only non-targeting gRNAs). Output: network_edges.tsv, mean_nontargeting_expression.tsv, network_plot.png.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu): mod/gene counts, mod/guide counts + layers['guide_assignment'], guide.var targeting / intended_target_name, uns pairs_to_test"},
            "results": {**_S_STRING, "description": "test_results table"},
            "central_nodes": {**_S_ARRAY_S, "description": "central elements"},
            "n_central": {**_S_INTEGER, "default": 2, "description": "central nodes when none given"},
            "min_weight": {**_S_NUMBER, "default": 0.01, "description": "minimum |log2_fc|"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip the figure"},
            "label": {**_S_STRING, "default": "network", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree2", "network"],
        flag_map={"mudata": "--mudata", "results": "--results", "central_nodes": "--central-nodes", "n_central": "--n-central", "min_weight": "--min-weight", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'central_nodes'},
        bool_flags={'no_plots'},
    ),

    _T(
        "crispr_jamboree2_run",
        "Runs the jamboree-2 benchmark on one MuData: guide assignment if the layer is missing, the chosen inference methods, evaluation against planted effects or control pairs (AUPRC / AUROC per method), volcano / IGV tracks and a network for the top elements. Output: one run directory with every table, figures, report.md and summary.json.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "IGVF CRISPR MuData (.h5mu): mod/gene counts, mod/guide counts + layers['guide_assignment'], guide.var targeting / intended_target_name, uns pairs_to_test"},
            "methods": {**_S_ARRAY_S, "default": ["wilcoxon", "negbinom", "scanpy-wilcoxon", "sceptre"], "description": "inference methods"},
            "assign": {**_S_STRING, "description": "(re)assign guides first: umi-threshold, cleanser or sceptre"},
            "threshold": {**_S_NUMBER, "default": 5.0, "description": "assignment threshold"},
            "side": {**_S_STRING, "default": "both", "description": "left, right or both"},
            "resamples": {**_S_INTEGER, "default": 499, "description": "sceptre resamples"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "reproduce upstream quirks"},
            "log2fc_threshold": {**_S_NUMBER, "default": 2.0, "description": "volcano |log2 FC| threshold"},
            "sig_threshold": {**_S_NUMBER, "default": 0.01, "description": "volcano significance threshold"},
            "out": {**_S_STRING, "description": "output .h5mu with test_results"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "jamboree2_run", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree2", "run"],
        flag_map={"mudata": "--mudata", "methods": "--methods", "assign": "--assign", "threshold": "--threshold", "side": "--side", "resamples": "--resamples", "upstream_compat": "--upstream-compat", "log2fc_threshold": "--log2fc-threshold", "sig_threshold": "--sig-threshold", "out": "--out", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'methods'},
        bool_flags={'upstream_compat', 'no_plots'},
    ),

    _T(
        "crispr_jamboree3_configure",
        "Turns the pipeline's configuration table (configuration.csv from the web form: tab_name, variable, variable_value, batch_name, read1, read2) into the Nextflow pipeline_input.config exactly as parse_interface_configuration.py (byte-identical on the upstream example), plus the covariate list (parse_covariate.csv) and the sceptre covariate formula of prepare_formula.py (cov_string.txt).",
        {"type": "object", "properties": {
            "config_table": {**_S_STRING, "description": "configuration.csv"},
            "mapping_tabs": {**_S_ARRAY_S, "default": ["scRNA", "Guides", "Hash"], "description": "mapping tabs"},
            "label": {**_S_STRING, "default": "configure", "description": "run label for the output directory"}},
         "required": ["config_table"]},
        cli=["crispr-jamboree3", "configure"],
        flag_map={"config_table": "--config-table", "mapping_tabs": "--mapping-tabs", "label": "--label"},
        flag_repeat={'mapping_tabs'},
    ),

    _T(
        "crispr_jamboree3_seqspec_check",
        "Checks FASTQ layouts as seqSpecCheck.py: per-position A/C/G/T frequencies of R1 and R2 (first max_reads reads) and a position table of where the metadata sequences (2nd column: guide or HTO sequences) start in R2. Inputs: FASTQ pairs + metadata. Output: position_table.csv, frequency tables, seqSpec_check_plots.png.",
        {"type": "object", "properties": {
            "read1": {**_S_ARRAY_S, "description": "R1 FASTQ(.gz) files"},
            "read2": {**_S_ARRAY_S, "description": "R2 FASTQ(.gz) files"},
            "metadata": {**_S_STRING, "description": "guide / hashing metadata (sequences in column 2)"},
            "max_reads": {**_S_INTEGER, "default": 100000, "description": "reads per file"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip the figure"},
            "label": {**_S_STRING, "default": "seqspec_check", "description": "run label for the output directory"}},
         "required": ["read1", "read2", "metadata"]},
        cli=["crispr-jamboree3", "seqspec-check"],
        flag_map={"read1": "--read1", "read2": "--read2", "metadata": "--metadata", "max_reads": "--max-reads", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'read1', 'read2'},
        bool_flags={'no_plots'},
    ),

    _T(
        "crispr_jamboree3_seqspec_parse",
        "Parses a seqspec YAML as parsing_guide_metadata.py: the reads of each modality, the kb `-x` representation (barcode:umi:cdna positions) and the barcode whitelist filename, written as parsed_seqSpec.txt. Input: seqspec YAML. Output: TSV + report.",
        {"type": "object", "properties": {
            "yaml": {**_S_STRING, "description": "seqspec YAML"},
            "modality": {**_S_ARRAY_S, "description": "rna, guide or hashing"},
            "directory": {**_S_STRING, "description": "seqspec directory recorded in the table"},
            "label": {**_S_STRING, "default": "seqspec_parse", "description": "run label for the output directory"}},
         "required": ["yaml", "modality"]},
        cli=["crispr-jamboree3", "seqspec-parse"],
        flag_map={"yaml": "--yaml", "modality": "--modality", "directory": "--directory", "label": "--label"},
        flag_repeat={'modality'},
    ),

    _T(
        "crispr_jamboree3_map",
        "Writes the kallisto | bustools mapping commands of the pipeline (kb ref -d human for RNA, kb ref --workflow kite with guide_features.txt / hashing_table.txt, kb count --h5ad -x <chemistry> -o <batch>_ks_<modality>_out) and runs them with execute when kb is on PATH. Chemistry and whitelist come from the seqspec YAML when given.",
        {"type": "object", "properties": {
            "modality": {**_S_STRING, "description": "rna, guide or hashing"},
            "seqspec_yaml": {**_S_STRING, "description": "seqspec YAML"},
            "chemistry": {**_S_STRING, "description": "kb -x string"},
            "whitelist": {**_S_STRING, "description": "barcode whitelist"},
            "metadata": {**_S_STRING, "description": "guide / hashing spreadsheet"},
            "fastqs": {**_S_ARRAY_S, "description": "'batch:R1 R2' strings per batch"},
            "genome": {**_S_STRING, "default": "genome.fa.gz", "description": "genome FASTA"},
            "species": {**_S_STRING, "default": "human", "description": "kb ref -d species"},
            "threads": {**_S_INTEGER, "default": 4, "description": "threads"},
            "execute": {**_S_BOOLEAN, "default": False, "description": "run the commands (needs kb)"},
            "label": {**_S_STRING, "default": "map", "description": "run label for the output directory"}},
         "required": ["modality"]},
        cli=["crispr-jamboree3", "map"],
        flag_map={"modality": "--modality", "seqspec_yaml": "--seqspec-yaml", "chemistry": "--chemistry", "whitelist": "--whitelist", "metadata": "--metadata", "fastqs": "--fastqs", "genome": "--genome", "species": "--species", "threads": "--threads", "execute": "--execute", "label": "--label"},
        flag_repeat={'fastqs'},
        bool_flags={'execute'},
    ),

    _T(
        "crispr_jamboree3_concat",
        "Concatenates per-batch kb outputs as anndata_concat.py: <batch>_ks_*/counts_unfiltered/adata.h5ad sorted by name, batch from the directory prefix joined to the covariate table, outer join with index_unique '_'. Output: concatenated_adata.h5ad.",
        {"type": "object", "properties": {
            "inputs": {**_S_ARRAY_S, "description": "<batch>_ks_*_out directories or .h5ad files"},
            "covariates": {**_S_STRING, "description": "parse_covariate.csv"},
            "out": {**_S_STRING, "description": "output .h5ad"},
            "label": {**_S_STRING, "default": "concat", "description": "run label for the output directory"}},
         "required": ["inputs"]},
        cli=["crispr-jamboree3", "concat"],
        flag_map={"inputs": "--inputs", "covariates": "--covariates", "out": "--out", "label": "--label"},
        flag_repeat={'inputs'},
    ),

    _T(
        "crispr_jamboree3_preprocess",
        "scRNA QC of preprocess_adata.py: gene symbols, version-stripped Ensembl ids, batch_number, mt / ribo flags, scanpy calculate_qc_metrics columns, filters min_genes (500), min_cells (3) and pct_counts_mt < pct_mito (20), knee plot. Input: RNA counts .h5ad. Output: filtered_anndata.h5ad + unfiltered QC AnnData.",
        {"type": "object", "properties": {
            "adata": {**_S_STRING, "description": "RNA counts .h5ad"},
            "gene_names": {**_S_STRING, "description": "cells_x_genes.genes.names.txt"},
            "min_genes": {**_S_INTEGER, "default": 500, "description": "minimum genes per cell"},
            "min_cells": {**_S_INTEGER, "default": 3, "description": "minimum cells per gene"},
            "pct_mito": {**_S_NUMBER, "default": 20.0, "description": "maximum percent mitochondrial counts"},
            "reference": {**_S_STRING, "default": "human", "description": "human or mouse"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "always use the MT- prefix as upstream"},
            "out": {**_S_STRING, "description": "output .h5ad"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip the knee plot"},
            "label": {**_S_STRING, "default": "preprocess", "description": "run label for the output directory"}},
         "required": ["adata"]},
        cli=["crispr-jamboree3", "preprocess"],
        flag_map={"adata": "--adata", "gene_names": "--gene-names", "min_genes": "--min-genes", "min_cells": "--min-cells", "pct_mito": "--pct-mito", "reference": "--reference", "upstream_compat": "--upstream-compat", "out": "--out", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'upstream_compat', 'no_plots'},
    ),

    _T(
        "crispr_jamboree3_create_mudata",
        "Builds the pipeline MuData as create_mdata(_HASHING).py: guide var from the guide metadata (intended target, coordinates, sequence, targeting), var_names sgRNA_ID|sequence, uns moi / capture_method, guide obs totals, gene coordinates from the GTF, renamed RNA QC columns, barcodes intersected across RNA / guide [/ hashing]. Output: mudata.h5mu.",
        {"type": "object", "properties": {
            "rna": {**_S_STRING, "description": "filtered RNA .h5ad"},
            "guide": {**_S_STRING, "description": "guide counts .h5ad"},
            "hashing": {**_S_STRING, "description": "demultiplexed hashing .h5ad"},
            "guide_metadata": {**_S_STRING, "description": "guide metadata (sgRNA_ID, sgRNA_sequences, Target_name, chr, start, end)"},
            "gtf": {**_S_STRING, "description": "GENCODE GTF(.gz)"},
            "moi": {**_S_STRING, "default": "high", "description": "low or high"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "positional metadata assignment and targeting TRUE for all, as upstream"},
            "out": {**_S_STRING, "description": "output .h5mu"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip the guide knee plot"},
            "label": {**_S_STRING, "default": "create_mudata", "description": "run label for the output directory"}},
         "required": ["rna", "guide", "guide_metadata"]},
        cli=["crispr-jamboree3", "create-mudata"],
        flag_map={"rna": "--rna", "guide": "--guide", "hashing": "--hashing", "guide_metadata": "--guide-metadata", "gtf": "--gtf", "moi": "--moi", "upstream_compat": "--upstream-compat", "out": "--out", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'upstream_compat', 'no_plots'},
    ),

    _T(
        "crispr_jamboree3_doublets",
        "Scrublet doublet detection on the gene modality (doublets.py defaults: 2x simulated doublets, expected rate 0.1, kNN doublet score, threshold_minimum on the simulated scores), removes predicted doublets from every modality. Output: mdata_doublets.h5mu + histogram.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "MuData (.h5mu) from create-mudata / earlier stages"},
            "seed": {**_S_INTEGER, "default": 0, "description": "random seed"},
            "out": {**_S_STRING, "description": "output .h5mu"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip the histogram"},
            "label": {**_S_STRING, "default": "doublets", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree3", "doublets"],
        flag_map={"mudata": "--mudata", "seed": "--seed", "out": "--out", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "crispr_jamboree3_demultiplex",
        "Demultiplexes cell hashing like GMM-Demux + demultiplex_filter.py + filter_hashing.py: per-HTO two-component Gaussian mixture on log1p counts, hto_type (positive HTO set), hto_type_split ('multiplets' when > 1), negatives dropped, barcodes intersected with the filtered RNA, per batch. Output: concatenated_hashing_demux.h5ad and GMM_full.csv / .config per batch.",
        {"type": "object", "properties": {
            "hashing": {**_S_STRING, "description": "hashing counts .h5ad"},
            "rna": {**_S_STRING, "description": "filtered RNA .h5ad to intersect with"},
            "out": {**_S_STRING, "description": "output .h5ad"},
            "label": {**_S_STRING, "default": "demultiplex", "description": "run label for the output directory"}},
         "required": ["hashing"]},
        cli=["crispr-jamboree3", "demultiplex"],
        flag_map={"hashing": "--hashing", "rna": "--rna", "out": "--out", "label": "--label"},
    ),

    _T(
        "crispr_jamboree3_assign_guides",
        "Assigns gRNAs to cells: CLEANSER Stan mixtures (per gRNA; upstream's single pooled fit with upstream_compat), the sceptre Poisson mixture, or a UMI threshold. Output: .h5mu with layers['guide_assignment'] and guide_assignment.csv.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "MuData (.h5mu) from create-mudata / earlier stages"},
            "method": {**_S_STRING, "default": "cleanser", "description": "cleanser, sceptre or umi-threshold"},
            "threshold": {**_S_NUMBER, "description": "posterior or UMI threshold"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "CLEANSER pooled fit as cleanser.py"},
            "out": {**_S_STRING, "description": "output .h5mu"},
            "label": {**_S_STRING, "default": "assign_guides", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree3", "assign-guides"],
        flag_map={"mudata": "--mudata", "method": "--method", "threshold": "--threshold", "upstream_compat": "--upstream-compat", "out": "--out", "label": "--label"},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "crispr_jamboree3_prepare_inference",
        "Builds uns['pairs_to_test'] as prepare_inference.py (user pairs restricted to genes and guides in the MuData, gene_name -> gene_id) or create_pairs_to_test.py (every GTF gene within distance_from_center of each guide). Output: mudata_inference_input.h5mu.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "MuData (.h5mu) from create-mudata / earlier stages"},
            "pairs": {**_S_STRING, "description": "user pairs CSV: guide_id, gene_name, intended_target_name, pair_type"},
            "gtf": {**_S_STRING, "description": "GTF for distance-based pairs"},
            "distance_from_center": {**_S_INTEGER, "default": 1000000, "description": "maximum guide-gene distance (bp)"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "GTF pairs carry gene symbols as upstream"},
            "out": {**_S_STRING, "description": "output .h5mu"},
            "label": {**_S_STRING, "default": "prepare_inference", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree3", "prepare-inference"],
        flag_map={"mudata": "--mudata", "pairs": "--pairs", "gtf": "--gtf", "distance_from_center": "--distance-from-center", "upstream_compat": "--upstream-compat", "out": "--out", "label": "--label"},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "crispr_jamboree3_infer",
        "Perturbation inference as inference_sceptre.R: a sceptre re-implementation (formula log(response_n_nonzero) + log(response_n_umis), optional non-redundant covariates, CRT / permutation null with skew-normal p) with element results joined onto every guide pair, or PerTurbo when installed. Output: test_results.csv and inference_mudata.h5mu.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "MuData with uns pairs_to_test"},
            "method": {**_S_STRING, "default": "sceptre", "description": "sceptre or perturbo"},
            "side": {**_S_STRING, "default": "both", "description": "left, right or both"},
            "formula": {**_S_STRING, "default": "default", "description": "'default' or another value to include the covariates"},
            "cov_string": {**_S_STRING, "description": "covariates, e.g. 'batch + cov1'"},
            "cov_string_file": {**_S_STRING, "description": "cov_string.txt"},
            "resamples": {**_S_INTEGER, "default": 499, "description": "resamples"},
            "seed": {**_S_INTEGER, "default": 4, "description": "seed"},
            "moi": {**_S_STRING, "description": "override MOI"},
            "out": {**_S_STRING, "description": "output .h5mu"},
            "label": {**_S_STRING, "default": "infer", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree3", "infer"],
        flag_map={"mudata": "--mudata", "method": "--method", "side": "--side", "formula": "--formula", "cov_string": "--cov-string", "cov_string_file": "--cov-string-file", "resamples": "--resamples", "seed": "--seed", "moi": "--moi", "out": "--out", "label": "--label"},
    ),

    _T(
        "crispr_jamboree3_evaluate",
        "Evaluation sub-workflow: volcano_plot.py (|log2_fc| >= 1, p <= 0.05, top-10 annotations), select_nodes.py + network_plot.py (edges |log2_fc| >= 0.1 around the most significant targets) and igv.py (promoter bedgraph when the target is the gene's GTF name, else bedpe). Output: evaluation_output/ with tables, output.bedpe/.bedgraph and figures.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "inference_mudata.h5mu"},
            "gtf": {**_S_STRING, "description": "GTF for gene names"},
            "central_nodes": {**_S_ARRAY_S, "description": "central targets (default: top by p)"},
            "num_nodes": {**_S_INTEGER, "default": 2, "description": "central nodes to select"},
            "log2_fc": {**_S_NUMBER, "default": 1.0, "description": "volcano |log2 FC|"},
            "p_value": {**_S_NUMBER, "default": 0.05, "description": "volcano p threshold"},
            "min_weight": {**_S_NUMBER, "default": 0.1, "description": "network minimum |log2_fc|"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "evaluate", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree3", "evaluate"],
        flag_map={"mudata": "--mudata", "gtf": "--gtf", "central_nodes": "--central-nodes", "num_nodes": "--num-nodes", "log2_fc": "--log2-fc", "p_value": "--p-value", "min_weight": "--min-weight", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'central_nodes'},
        bool_flags={'no_plots'},
    ),

    _T(
        "crispr_jamboree3_dashboard",
        "Builds the pipeline dashboard (create_dashboard_plots.py / create_dashboard_df.py / process_json.py): UMI, detected-gene and guide-UMI threshold counts, sgRNA frequencies, guides-per-cell / cells-per-guide histograms, highlight blocks (barcode intersections, cells after filtering, significant and direct-targeting pairs), kb mapping tables. Output: self-contained dashboard.html, blocks TSV, figures.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "inference_mudata.h5mu"},
            "gene_ann": {**_S_STRING, "description": "unfiltered RNA .h5ad"},
            "gene_ann_filtered": {**_S_STRING, "description": "filtered RNA .h5ad"},
            "guide_ann": {**_S_STRING, "description": "guide counts .h5ad"},
            "guide_fq_tbl": {**_S_STRING, "description": "seqspec-check position_table.csv"},
            "kb_dirs": {**_S_ARRAY_S, "description": "<batch>_ks_*_out directories with inspect.json/run_info.json"},
            "extra_figures": {**_S_ARRAY_S, "description": "extra PNGs to embed"},
            "title": {**_S_STRING, "default": "CRISPR Perturb-seq pipeline dashboard", "description": "dashboard title"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "dashboard", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree3", "dashboard"],
        flag_map={"mudata": "--mudata", "gene_ann": "--gene-ann", "gene_ann_filtered": "--gene-ann-filtered", "guide_ann": "--guide-ann", "guide_fq_tbl": "--guide-fq-tbl", "kb_dirs": "--kb-dirs", "extra_figures": "--extra-figures", "title": "--title", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'kb_dirs', 'extra_figures'},
        bool_flags={'no_plots'},
    ),

    _T(
        "crispr_jamboree3_run",
        "Runs the jamboree-3 single-cell pipeline from count matrices: preprocess -> [hashing demultiplex] -> create-mudata -> Scrublet doublets -> guide assignment -> pairs to test -> sceptre inference -> evaluation (volcano / network / IGV) -> dashboard. Inputs: RNA and guide (optionally hashing) count AnnDatas, guide metadata, GTF, pairs CSV. Output: one run directory with inference_mudata.h5mu, test_results.csv, evaluation_output/, dashboard.html, report.md, summary.json.",
        {"type": "object", "properties": {
            "rna": {**_S_STRING, "description": "RNA counts .h5ad (obs batch)"},
            "guide": {**_S_STRING, "description": "guide counts .h5ad"},
            "hashing": {**_S_STRING, "description": "hashing counts .h5ad"},
            "guide_metadata": {**_S_STRING, "description": "guide metadata table"},
            "gene_names": {**_S_STRING, "description": "gene symbols file (default: from the GTF)"},
            "gtf": {**_S_STRING, "description": "GTF(.gz)"},
            "pairs": {**_S_STRING, "description": "user pairs CSV (default: GTF distance pairs)"},
            "moi": {**_S_STRING, "default": "high", "description": "low or high"},
            "min_genes": {**_S_INTEGER, "default": 500, "description": "minimum genes per cell"},
            "min_cells": {**_S_INTEGER, "default": 3, "description": "minimum cells per gene"},
            "pct_mito": {**_S_NUMBER, "default": 20.0, "description": "maximum percent mito"},
            "reference": {**_S_STRING, "default": "human", "description": "human or mouse"},
            "assignment_method": {**_S_STRING, "default": "sceptre", "description": "cleanser, sceptre or umi-threshold"},
            "threshold": {**_S_NUMBER, "default": 1.0, "description": "assignment threshold"},
            "side": {**_S_STRING, "default": "both", "description": "left, right or both"},
            "formula": {**_S_STRING, "default": "default", "description": "sceptre formula ('default')"},
            "cov_string": {**_S_STRING, "description": "covariates"},
            "resamples": {**_S_INTEGER, "default": 499, "description": "sceptre resamples"},
            "distance_from_center": {**_S_INTEGER, "default": 1000000, "description": "GTF pair distance"},
            "central_nodes": {**_S_ARRAY_S, "description": "network centres"},
            "num_nodes": {**_S_INTEGER, "default": 2, "description": "central nodes"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "reproduce upstream quirks"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip figures"},
            "label": {**_S_STRING, "default": "jamboree3_run", "description": "run label for the output directory"}},
         "required": ["rna", "guide", "guide_metadata"]},
        cli=["crispr-jamboree3", "run"],
        flag_map={"rna": "--rna", "guide": "--guide", "hashing": "--hashing", "guide_metadata": "--guide-metadata", "gene_names": "--gene-names", "gtf": "--gtf", "pairs": "--pairs", "moi": "--moi", "min_genes": "--min-genes", "min_cells": "--min-cells", "pct_mito": "--pct-mito", "reference": "--reference", "assignment_method": "--assignment-method", "threshold": "--threshold", "side": "--side", "formula": "--formula", "cov_string": "--cov-string", "resamples": "--resamples", "distance_from_center": "--distance-from-center", "central_nodes": "--central-nodes", "num_nodes": "--num-nodes", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'central_nodes'},
        bool_flags={'upstream_compat', 'no_plots'},
    ),

    _T(
        "crispr_jamboree4_prepare",
        "Prepares the inference checkpoint of the 2025 jamboree task (prepare_user_guide_inference / prepare_inference.py): the pairs-to-test sheet (guide_id, gene_name, intended_target_name, pair_type) restricted to genes and guides present in the guide-assigned MuData, gene_name -> gene_id, stored as uns['pairs_to_test']. Output: mudata_inference_input.h5mu.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "guide-assigned MuData (.h5mu)"},
            "pairs": {**_S_STRING, "description": "pairs CSV: guide_id, gene_name, intended_target_name, pair_type"},
            "out": {**_S_STRING, "description": "output .h5mu"},
            "label": {**_S_STRING, "default": "prepare", "description": "run label for the output directory"}},
         "required": ["mudata", "pairs"]},
        cli=["crispr-jamboree4", "prepare"],
        flag_map={"mudata": "--mudata", "pairs": "--pairs", "out": "--out", "label": "--label"},
    ),

    _T(
        "crispr_jamboree4_infer",
        "Runs the benchmarked inference methods on uns['pairs_to_test']: a sceptre re-implementation (inference_sceptre.R defaults), a glm.nb comparator, and PerTurbo when the package is installed. Output: test_results.<method>.tsv with p_value / log2_fc per pair.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "MuData with uns pairs_to_test"},
            "methods": {**_S_ARRAY_S, "default": ["sceptre", "perturbo"], "description": "sceptre, negbinom, perturbo"},
            "side": {**_S_STRING, "default": "both", "description": "left, right or both"},
            "resamples": {**_S_INTEGER, "default": 499, "description": "sceptre resamples"},
            "seed": {**_S_INTEGER, "default": 4, "description": "seed"},
            "formula": {**_S_STRING, "default": "default", "description": "'default' or another value to include the covariates"},
            "cov_string": {**_S_STRING, "description": "covariates, e.g. 'batch'"},
            "label": {**_S_STRING, "default": "infer", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree4", "infer"],
        flag_map={"mudata": "--mudata", "methods": "--methods", "side": "--side", "resamples": "--resamples", "seed": "--seed", "formula": "--formula", "cov_string": "--cov-string", "label": "--label"},
        flag_repeat={'methods'},
    ),

    _T(
        "crispr_jamboree4_merge",
        "Merges method results as mergedResults / export_output_multiple.py: sceptre_* and perturbo_* columns inner-joined on (guide_id, gene_id), inference_mudata.h5mu with uns test_results, per_guide_output.tsv (upstream column set, cell_number, avg_gene_expression) and per_element_output.tsv (guide ids joined per element).",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "mudata_inference_input.h5mu"},
            "sceptre_results": {**_S_STRING, "description": "sceptre test_results table or .h5mu"},
            "perturbo_results": {**_S_STRING, "description": "PerTurbo test_results table or .h5mu"},
            "extra_results": {**_S_ARRAY_S, "description": "name=path for other methods"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "group the per-element table by gene_id only, as upstream"},
            "label": {**_S_STRING, "default": "merge", "description": "run label for the output directory"}},
         "required": ["mudata"]},
        cli=["crispr-jamboree4", "merge"],
        flag_map={"mudata": "--mudata", "sceptre_results": "--sceptre-results", "perturbo_results": "--perturbo-results", "extra_results": "--extra-results", "upstream_compat": "--upstream-compat", "label": "--label"},
        flag_repeat={'extra_results'},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "crispr_jamboree4_evaluate",
        "Reconstructs evaluation_curve.py: true labels from pair_type (positive / direct-targeting controls vs negative / non-targeting controls), AUPRC (trapezoid auc of the PR curve), AUROC and average precision for every <method>_p_value column scored by -log10 p, for one or more datasets (e.g. H9 and WTC11). Output: evaluation_metrics.tsv and the PR / ROC grid (evaluation_curves.png).",
        {"type": "object", "properties": {
            "results": {**_S_ARRAY_S, "description": "[name=]path to inference_mudata.h5mu / test_results / per_guide_output.tsv"},
            "pairs": {**_S_STRING, "description": "pairs sheet joined when a table lacks pair_type"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip the figure"},
            "label": {**_S_STRING, "default": "evaluate", "description": "run label for the output directory"}},
         "required": ["results"]},
        cli=["crispr-jamboree4", "evaluate"],
        flag_map={"results": "--results", "pairs": "--pairs", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'results'},
        bool_flags={'no_plots'},
    ),

    _T(
        "crispr_jamboree4_run",
        "Runs the 2025 jamboree benchmark task end to end on one dataset: prepare the pairs, run the chosen methods, merge the results into the pipeline outputs and report AUPRC / AUROC of each method on the control pairs. Inputs: guide-assigned MuData + pairs sheet. Output: run directory with inference_mudata.h5mu, per_guide/per_element outputs, curves, report.md, summary.json.",
        {"type": "object", "properties": {
            "mudata": {**_S_STRING, "description": "guide-assigned MuData"},
            "pairs": {**_S_STRING, "description": "pairs CSV"},
            "methods": {**_S_ARRAY_S, "default": ["sceptre", "perturbo"], "description": "sceptre, negbinom, perturbo"},
            "side": {**_S_STRING, "default": "both", "description": "left, right or both"},
            "resamples": {**_S_INTEGER, "default": 499, "description": "sceptre resamples"},
            "seed": {**_S_INTEGER, "default": 4, "description": "seed"},
            "formula": {**_S_STRING, "default": "default", "description": "sceptre formula"},
            "cov_string": {**_S_STRING, "description": "covariates"},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "per-element table by gene_id"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "skip the figure"},
            "label": {**_S_STRING, "default": "jamboree4_run", "description": "run label for the output directory"}},
         "required": ["mudata", "pairs"]},
        cli=["crispr-jamboree4", "run"],
        flag_map={"mudata": "--mudata", "pairs": "--pairs", "methods": "--methods", "side": "--side", "resamples": "--resamples", "seed": "--seed", "formula": "--formula", "cov_string": "--cov-string", "upstream_compat": "--upstream-compat", "no_plots": "--no-plots", "label": "--label"},
        flag_repeat={'methods'},
        bool_flags={'upstream_compat', 'no_plots'},
    ),

    _T(
        "principal_pseudobulks_fetch",
        "Finds, processed-first, the IGVF Portal inputs for principal pseudobulks from any accession: walks the lineage (portal_lineage) to the principal analysis set's cell-annotation table and, per uniform-pipeline intermediate analysis set (one 10x lane), the best ATAC fragments file and RNA h5ad. Writes portal_manifest.tsv (lane, analysis set, local paths) and downloads within a GB budget via portal_download; dry_run only lists what would be fetched.",
        {"type": "object", "properties": {
            "accession": {**_S_STRING, "description": "Any IGVF accession: principal/intermediate analysis set, measurement set or sample."},
            "max_gb": {**_S_NUMBER, "default": 20.0, "description": "Total download budget in GB."},
            "dry_run": {**_S_BOOLEAN, "default": False, "description": "Plan and list only, download nothing."},
            "include_lab": {**_S_BOOLEAN, "default": False, "description": "Also use lab (non-uniform) analysis sets."},
            "dest": {**_S_STRING, "description": "Download directory."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["accession"]},
        cli=["principal-pseudobulks", "fetch"],
        flag_map={"accession": "--accession", "max_gb": "--max-gb", "dry_run": "--dry-run", "include_lab": "--include-lab", "dest": "--dest", "label": "--label"},
        bool_flags={'dry_run', 'include_lab'},
    ),

    _T(
        "principal_pseudobulks_build",
        "Builds the primary pseudobulk layout (annotation-{cell_type}-{subsample}/ with fragments.tsv.gz, rna_counts_mtx.h5ad and per_cell_qc.tsv) from a cell-annotation table plus per-lane fragments and RNA h5ad files. Barcodes become {16-bp}_{IGVF lane accession} (GEM -1 stripped, optional ATAC->RNA multiome barcode translation); per-cell QC computes UMIs, genes, % mito/ribo, fragments, % duplicates, nucleosomal signal, TSS enrichment (with a TSS bed) and FRiP (with peaks).",
        {"type": "object", "properties": {
            "annotations": {**_S_STRING, "description": "Cell annotation table (.tsv/.csv/.h5ad obs) with barcode and cell type columns."},
            "manifest": {**_S_STRING, "description": "Lane manifest TSV from fetch (lane, fragments, rna_matrix)."},
            "dataset": {**_S_STRING, "description": "Dataset name used in the output layout."},
            "annotation_column": {**_S_STRING, "description": "Cell-type column (auto-detected)."},
            "subsample_column": {**_S_STRING, "description": "Subsample column (auto-detected)."},
            "lane_column": {**_S_STRING, "description": "Lane column when barcodes carry no lane suffix."},
            "atac_rna_map": {**_S_STRING, "description": "Two-column ATAC barcode -> RNA barcode table."},
            "tss": {**_S_STRING, "description": "TSS bed for tss_enrichment."},
            "peaks": {**_S_STRING, "description": "Peak bed for FRiP."},
            "min_cells": {**_S_INTEGER, "default": 1, "description": "Skip groups with fewer cells."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["annotations", "dataset"]},
        cli=["principal-pseudobulks", "build-pseudobulks"],
        flag_map={"annotations": "--annotations", "manifest": "--manifest", "dataset": "--dataset", "annotation_column": "--annotation-column", "subsample_column": "--subsample-column", "lane_column": "--lane-column", "atac_rna_map": "--atac-rna-map", "tss": "--tss", "peaks": "--peaks", "min_cells": "--min-cells", "label": "--label"},
    ),

    _T(
        "principal_pseudobulks_qc_datatable",
        "Concatenates the per_cell_qc.tsv of every annotation-{cell_type}-IGVF* subsample directory into one cluster-level per-cell QC datatable (upstream Step 0), skipping empty files and failing on a header mismatch. Output is the input to explore and qc-filter.",
        {"type": "object", "properties": {
            "pseudobulks": {**_S_STRING, "description": "The dataset's pseudobulks/ directory."},
            "cell_type": {**_S_STRING, "description": "Cell type string after 'annotation-'."},
            "out": {**_S_STRING, "description": "Output TSV path."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["pseudobulks", "cell_type"]},
        cli=["principal-pseudobulks", "qc-datatable"],
        flag_map={"pseudobulks": "--pseudobulks", "cell_type": "--cell-type", "out": "--out", "label": "--label"},
    ),

    _T(
        "principal_pseudobulks_explore",
        "Reports, for one or more candidate QC threshold sets, how many cells each threshold removes in total and alone (cells saved if only that threshold were dropped), the most stringent threshold and subsamples with surviving cells (explore_qc_thresholds.R; >= / <= semantics). Writes explore_report.txt and TSV tables; optional per-subsample breakdown.",
        {"type": "object", "properties": {
            "meta": {**_S_STRING, "description": "Per-cell QC datatable."},
            "rna_min": {**_S_NUMBER, "description": "Minimum RNA reads (UMIs) per cell (default 1000)."},
            "rna_max": {**_S_NUMBER, "description": "Maximum RNA reads per cell (default Inf)."},
            "gene_min": {**_S_NUMBER, "description": "Minimum genes per cell (explore default 1000, qc-filter default 0)."},
            "gene_max": {**_S_NUMBER, "description": "Maximum genes per cell."},
            "pct_mt_max": {**_S_NUMBER, "description": "Maximum % mitochondrial UMIs (default 30)."},
            "pct_ribo_max": {**_S_NUMBER, "description": "Maximum % ribosomal UMIs (default 100)."},
            "atac_min": {**_S_NUMBER, "description": "Minimum ATAC fragments per cell (default 1000)."},
            "atac_max": {**_S_NUMBER, "description": "Maximum ATAC fragments per cell."},
            "tss_enr_min": {**_S_NUMBER, "description": "Minimum TSS enrichment (default 3)."},
            "nuc_signal_max": {**_S_NUMBER, "description": "Maximum nucleosomal signal (default 1.5)."},
            "pct_dup_max": {**_S_NUMBER, "description": "Maximum % duplicated ATAC reads (default 100)."},
            "frip_min": {**_S_NUMBER, "description": "Minimum FRiP (default 0)."},
            "sets": {**_S_ARRAY_S, "description": "Threshold sets, each a flag string, e.g. \"--tss-min 5 --pct-mt-max 20\"."},
            "show_subsamples": {**_S_BOOLEAN, "default": False, "description": "Add the per-subsample breakdown."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["meta"]},
        cli=["principal-pseudobulks", "explore"],
        flag_map={"meta": "--meta", "rna_min": "--rna-min", "rna_max": "--rna-max", "gene_min": "--gene-min", "gene_max": "--gene-max", "pct_mt_max": "--pct-mt-max", "pct_ribo_max": "--pct-ribo-max", "atac_min": "--atac-min", "atac_max": "--atac-max", "tss_enr_min": "--tss-min", "nuc_signal_max": "--nuc-max", "pct_dup_max": "--pct-dup-max", "frip_min": "--frip-min", "sets": "--sets", "show_subsamples": "--show-subsamples", "label": "--label"},
        flag_repeat={'sets'},
        bool_flags={'show_subsamples'},
    ),

    _T(
        "principal_pseudobulks_qc_filter",
        "Applies the final per-cell QC thresholds with plot_per_cell_qc.R's strict inequalities and writes the QC guide filtered_barcodes_with_subsamples.tsv.gz (barcode, subsample, analysis_accession), qc_thresholds.tsv, filtered_cell_subsample_metrics.tsv and PNG RNA/ATAC QC figures plus cells-per-subsample bars after RNA, ATAC and all QC.",
        {"type": "object", "properties": {
            "meta": {**_S_STRING, "description": "Per-cell QC datatable."},
            "rna_min": {**_S_NUMBER, "description": "Minimum RNA reads (UMIs) per cell (default 1000)."},
            "rna_max": {**_S_NUMBER, "description": "Maximum RNA reads per cell (default Inf)."},
            "gene_min": {**_S_NUMBER, "description": "Minimum genes per cell (explore default 1000, qc-filter default 0)."},
            "gene_max": {**_S_NUMBER, "description": "Maximum genes per cell."},
            "pct_mt_max": {**_S_NUMBER, "description": "Maximum % mitochondrial UMIs (default 30)."},
            "pct_ribo_max": {**_S_NUMBER, "description": "Maximum % ribosomal UMIs (default 100)."},
            "atac_min": {**_S_NUMBER, "description": "Minimum ATAC fragments per cell (default 1000)."},
            "atac_max": {**_S_NUMBER, "description": "Maximum ATAC fragments per cell."},
            "tss_enr_min": {**_S_NUMBER, "description": "Minimum TSS enrichment (default 3)."},
            "nuc_signal_max": {**_S_NUMBER, "description": "Maximum nucleosomal signal (default 1.5)."},
            "pct_dup_max": {**_S_NUMBER, "description": "Maximum % duplicated ATAC reads (default 100)."},
            "frip_min": {**_S_NUMBER, "description": "Minimum FRiP (default 0)."},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip figures."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["meta"]},
        cli=["principal-pseudobulks", "qc-filter"],
        flag_map={"meta": "--meta", "rna_min": "--rna-min", "rna_max": "--rna-max", "gene_min": "--gene-min", "gene_max": "--gene-max", "pct_mt_max": "--pct-mt-max", "pct_ribo_max": "--pct-ribo-max", "atac_min": "--atac-min", "atac_max": "--atac-max", "tss_enr_min": "--tss-min", "nuc_signal_max": "--nuc-max", "pct_dup_max": "--pct-dup-max", "frip_min": "--frip-min", "no_plots": "--no-plots", "label": "--label"},
        bool_flags={'no_plots'},
    ),

    _T(
        "principal_pseudobulks_filter_atac",
        "Filters every fragments.tsv.gz of a cell type's pseudobulk directories to the QC-guide barcodes (full barcode match), drops chromosomes absent from the chrom sizes, sorts in chrom-sizes order and writes a bgzipped, tabix-indexed fragment file. Fails if any guide barcode has no fragments or a record lacks 5 fields.",
        {"type": "object", "properties": {
            "qc_guide": {**_S_STRING, "description": "QC guide TSV.gz."},
            "pseudobulks": {**_S_STRING, "description": "pseudobulks/ directory."},
            "cell_type": {**_S_STRING, "description": "Cell type (comma-separated to merge)."},
            "chrom_sizes": {**_S_STRING, "description": "Chrom sizes (default: embedded IGVF DACC GRCh38)."},
            "out": {**_S_STRING, "description": "Output .tsv.gz."},
            "clean": {**_S_BOOLEAN, "default": False, "description": "Write the 16-bp barcode only."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["qc_guide", "pseudobulks", "cell_type"]},
        cli=["principal-pseudobulks", "filter-atac"],
        flag_map={"qc_guide": "--qc-guide", "pseudobulks": "--pseudobulks", "cell_type": "--cell-type", "chrom_sizes": "--chrom-sizes", "out": "--out", "clean": "--clean", "label": "--label"},
        bool_flags={'clean'},
    ),

    _T(
        "principal_pseudobulks_filter_rna",
        "Filters and concatenates the cell type's rna_counts_mtx.h5ad matrices to the QC-guide barcodes, maps versioned Ensembl IDs to gene symbols through the GTF (hard failure on any unmatched ID; optional standard chromosomes only), sums IDs sharing a symbol, checks no duplicate barcodes/genes and cell count == guide, and writes a genes x cells Matrix Market directory (or .csv.gz / .h5ad).",
        {"type": "object", "properties": {
            "qc_guide": {**_S_STRING, "description": "QC guide TSV.gz."},
            "pseudobulks": {**_S_STRING, "description": "pseudobulks/ directory."},
            "cell_type": {**_S_STRING, "description": "Cell type."},
            "out": {**_S_STRING, "description": "Output .mtx (directory), .csv.gz or .h5ad."},
            "gtf": {**_S_STRING, "description": "GENCODE 43 GTF (IGVFFI9573KOZR)."},
            "ensembl_ids_as_genes": {**_S_BOOLEAN, "default": False, "description": "Keep Ensembl IDs."},
            "standard_chromosomes_only": {**_S_BOOLEAN, "default": False, "description": "chr1-22, X, Y, M genes only."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["qc_guide", "pseudobulks", "cell_type"]},
        cli=["principal-pseudobulks", "filter-rna"],
        flag_map={"qc_guide": "--qc-guide", "pseudobulks": "--pseudobulks", "cell_type": "--cell-type", "out": "--out", "gtf": "--gtf", "ensembl_ids_as_genes": "--ensembl-ids-as-genes", "standard_chromosomes_only": "--standard-chromosomes-only", "label": "--label"},
        bool_flags={'standard_chromosomes_only', 'ensembl_ids_as_genes'},
    ),

    _T(
        "principal_pseudobulks_package_rna",
        "Packages an RNA count matrix directory into the distributed .tar.gz with decompressed matrix.mtx, barcodes.tsv and features.tsv members, as FILE_SPEC_RNA_COUNT_MATRIX.md describes.",
        {"type": "object", "properties": {
            "rna_matrix": {**_S_STRING, "description": "RNA matrix directory."},
            "out": {**_S_STRING, "description": "Output .tar.gz."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["rna_matrix"]},
        cli=["principal-pseudobulks", "package-rna"],
        flag_map={"rna_matrix": "--rna-matrix", "out": "--out", "label": "--label"},
    ),

    _T(
        "principal_pseudobulks_config_table",
        "Writes the per-dataset cluster config table (cluster, rna_matrix_file, atac_frag_file, HiC columns, alt_TSS, alt_genes, model_dir) that scE2G / sce2g-pipeline --cluster-config reads, for cell types already filtered under out_root.",
        {"type": "object", "properties": {
            "out_root": {**_S_STRING, "description": "The run's out_dir."},
            "dataset": {**_S_STRING, "description": "Dataset."},
            "cell_types": {**_S_ARRAY_S, "description": "Cluster names."},
            "out": {**_S_STRING, "description": "Output TSV."},
            "upstream_compat": {**_S_BOOLEAN, "default": False, "description": "model_dir as upstream's models/... paths."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["out_root", "dataset", "cell_types"]},
        cli=["principal-pseudobulks", "config-table"],
        flag_map={"out_root": "--out-root", "dataset": "--dataset", "cell_types": "--cell-types", "out": "--out", "upstream_compat": "--upstream-compat", "label": "--label"},
        flag_repeat={'cell_types'},
        bool_flags={'upstream_compat'},
    ),

    _T(
        "principal_pseudobulks_run",
        "Runs the whole upstream Snakefile: for every (dataset, cell_type) of a config_QC_pseudobulks.yaml/json, filter-atac and filter-rna into {out}/{dataset}/{cell_type}/ and the per-dataset config table; auto_qc builds the QC datatable and QC guide where missing; dry_run lists jobs with memory estimates.",
        {"type": "object", "properties": {
            "config": {**_S_STRING, "description": "config_QC_pseudobulks.yaml or .json."},
            "dry_run": {**_S_BOOLEAN, "default": False, "description": "List jobs only."},
            "auto_qc": {**_S_BOOLEAN, "default": False, "description": "Run qc-filter where the QC guide is missing."},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip figures."},
            "label": {**_S_STRING, "description": "Run-directory label."},
            "rna_min": {**_S_NUMBER, "description": "Minimum RNA reads (UMIs) per cell (default 1000)."},
            "rna_max": {**_S_NUMBER, "description": "Maximum RNA reads per cell (default Inf)."},
            "gene_min": {**_S_NUMBER, "description": "Minimum genes per cell (explore default 1000, qc-filter default 0)."},
            "gene_max": {**_S_NUMBER, "description": "Maximum genes per cell."},
            "pct_mt_max": {**_S_NUMBER, "description": "Maximum % mitochondrial UMIs (default 30)."},
            "pct_ribo_max": {**_S_NUMBER, "description": "Maximum % ribosomal UMIs (default 100)."},
            "atac_min": {**_S_NUMBER, "description": "Minimum ATAC fragments per cell (default 1000)."},
            "atac_max": {**_S_NUMBER, "description": "Maximum ATAC fragments per cell."},
            "tss_enr_min": {**_S_NUMBER, "description": "Minimum TSS enrichment (default 3)."},
            "nuc_signal_max": {**_S_NUMBER, "description": "Maximum nucleosomal signal (default 1.5)."},
            "pct_dup_max": {**_S_NUMBER, "description": "Maximum % duplicated ATAC reads (default 100)."},
            "frip_min": {**_S_NUMBER, "description": "Minimum FRiP (default 0)."}},
         "required": ["config"]},
        cli=["principal-pseudobulks", "run"],
        flag_map={"config": "--config", "dry_run": "--dry-run", "auto_qc": "--auto-qc", "no_plots": "--no-plots", "label": "--label", "rna_min": "--rna-min", "rna_max": "--rna-max", "gene_min": "--gene-min", "gene_max": "--gene-max", "pct_mt_max": "--pct-mt-max", "pct_ribo_max": "--pct-ribo-max", "atac_min": "--atac-min", "atac_max": "--atac-max", "tss_enr_min": "--tss-min", "nuc_signal_max": "--nuc-max", "pct_dup_max": "--pct-dup-max", "frip_min": "--frip-min"},
        bool_flags={'dry_run', 'no_plots', 'auto_qc'},
    ),

    _T(
        "principal_pseudobulks_validate",
        "Checks a QC guide, RNA count matrix directory and fragment file against every guarantee of the two upstream file specs: header and barcode format {16-bp}_{IGVF lane}, unique barcodes, integer genes x cells Matrix Market, gene symbols without duplicates, barcode sets identical, fragments sorted and tabix-indexed. Writes spec_checks.tsv; exits 1 on any failure.",
        {"type": "object", "properties": {
            "qc_guide": {**_S_STRING, "description": "QC guide TSV.gz."},
            "rna_matrix": {**_S_STRING, "description": "RNA matrix directory."},
            "fragments": {**_S_STRING, "description": "Filtered fragments .tsv.gz."},
            "label": {**_S_STRING, "description": "Run-directory label."}}},
        cli=["principal-pseudobulks", "validate"],
        flag_map={"qc_guide": "--qc-guide", "rna_matrix": "--rna-matrix", "fragments": "--fragments", "label": "--label"},
    ),

    _T(
        "principal_pseudobulks_sce2g_prep",
        "Turns a principal-pseudobulk config table into sce2g-pipeline inputs: per cluster tagAlign, fragment_count and cell barcodes (via sce2g_pipeline_skill.frag_to_tagalign) and RNA pseudobulk TPM / detection / mean log-norm (via rna_features), plus config_cell_clusters.tsv for sce2g-pipeline run --cluster-config.",
        {"type": "object", "properties": {
            "config_table": {**_S_STRING, "description": "{dataset}_config.tsv from run or config-table."},
            "chrom_sizes": {**_S_STRING, "description": "Chrom sizes (default embedded)."},
            "label": {**_S_STRING, "description": "Run-directory label."}},
         "required": ["config_table"]},
        cli=["principal-pseudobulks", "sce2g-prep"],
        flag_map={"config_table": "--config-table", "chrom_sizes": "--chrom-sizes", "label": "--label"},
    ),

    _T(
        "principal_pseudobulks_selftest",
        "Runs the offline self-test: synthetic 10x lanes, subsamples, cell types and genes with planted low-quality cells, a cross-lane barcode collision, ATAC->RNA barcode translation and GTF edge cases; every subcommand is asserted, including the fixture-portal fetch and the sce2g-pipeline hand-off.",
        {"type": "object", "properties": {
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "Skip figures."}}},
        cli=["principal-pseudobulks", "selftest"],
        flag_map={"no_plots": "--no-plots"},
        bool_flags={'no_plots'},
    ),

    _T(
        "e2g_qc_resolve_exclusions",
        "Computes the per-cluster quality gate of the E2G QC-and-Predictions pipeline from a pipeline config: cell / fragment / UMI totals from filtered_cell_subsample_metrics.tsv (or the prefiltered per-cell QC join, or CATlas prefiltered metrics) against min_cell_count 100, min_fragments_total 2e6, min_umi_count 1e6. Reports included, upload-eligible and manifest-eligible clusters with exclusion reasons. Writes <dataset>_cluster_stats.tsv, <dataset>_pipeline_plan.tsv, report.md and summary.json.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "exclusions", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "config": {**_S_STRING, "description": "*_pipeline_config.yaml (or .json)"},
            "output_dir": {**_S_STRING, "description": "override the config's output_dir"},
            "legacy_implicit_source": {**_S_BOOLEAN, "default": False, "description": "synapse-submission / CATlas stats-source rule"},
            "write_to_output_dir": {**_S_BOOLEAN, "default": False, "description": "write cluster_stats/ and igvf_metadata/ under output_dir"}},
         "required": ["config"]},
        cli=["e2g-qc-predictions", "resolve-exclusions"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "config": "--config", "output_dir": "--output-dir", "legacy_implicit_source": "--legacy-implicit-source", "write_to_output_dir": "--write-to-output-dir"},
        bool_flags={'write_to_output_dir', 'legacy_implicit_source'},
    ),

    _T(
        "e2g_qc_merge_metrics",
        "Merges the filtered_cell_subsample_metrics.tsv files of a merged cluster's components by subsample: counts summed, per-cell ratios recomputed, FRiP / TSS enrichment n_cells-weighted. Writes the merged metrics TSV the quality gate reads.",
        {"type": "object", "properties": {
            "metrics": {**_S_ARRAY_S, "description": "metrics"},
            "out": {**_S_STRING, "description": "out"},
            "force": {**_S_BOOLEAN, "default": False, "description": "force"}},
         "required": ["metrics", "out"]},
        cli=["e2g-qc-predictions", "merge-metrics"],
        flag_map={"metrics": "--metrics", "out": "--out", "force": "--force"},
        flag_repeat={'metrics'},
        bool_flags={'force'},
    ),

    _T(
        "e2g_qc_prefiltered_metrics",
        "Derives a filtered_cell_subsample_metrics.tsv row (distinct barcodes, fragment lines) from an already-filtered ATAC fragments file, as used for CATlas clusters. Writes the metrics TSV.",
        {"type": "object", "properties": {
            "frag_file": {**_S_STRING, "description": "frag file"},
            "subsample_name": {**_S_STRING, "description": "subsample name"},
            "out": {**_S_STRING, "description": "out"}},
         "required": ["frag_file", "subsample_name", "out"]},
        cli=["e2g-qc-predictions", "prefiltered-metrics"],
        flag_map={"frag_file": "--frag-file", "subsample_name": "--subsample-name", "out": "--out"},
    ),

    _T(
        "e2g_qc_build_qc_datatables",
        "Builds <dataset>_data/<cluster>_per_cell_qc.tsv tables by concatenating each pseudobulk's per_cell_qc file (sorted subsample directories, header asserted identical), merging comma-listed annotations per the pipeline config. Writes one TSV per cluster under QC_datatables/.",
        {"type": "object", "properties": {
            "pseudobulks_root": {**_S_STRING, "description": "pseudobulks root"},
            "dest": {**_S_STRING, "description": "QC_datatables/ is created under here"},
            "datasets": {**_S_ARRAY_S, "description": "datasets"},
            "config": {**_S_STRING, "description": "pipeline config supplying merged-cluster annotations"},
            "force": {**_S_BOOLEAN, "default": False, "description": "force"}},
         "required": ["pseudobulks_root", "dest"]},
        cli=["e2g-qc-predictions", "build-qc-datatables"],
        flag_map={"pseudobulks_root": "--pseudobulks-root", "dest": "--dest", "datasets": "--datasets", "config": "--config", "force": "--force"},
        flag_repeat={'datasets'},
        bool_flags={'force'},
    ),

    _T(
        "e2g_qc_filter_atac",
        "Filters ATAC fragments of one cluster's pseudobulk directories to the QC guide's full barcodes, drops chromosomes absent from the chrom-sizes file, sorts in chrom-sizes order and writes a bgzipped, tabix-indexed fragments file. Fails if any guide barcode is never seen (legacy branches: --allow-missing-barcodes).",
        {"type": "object", "properties": {
            "qc_guide": {**_S_STRING, "description": "qc guide"},
            "pseudobulks": {**_S_STRING, "description": "pseudobulks"},
            "cell_type": {**_S_STRING, "description": "cell type"},
            "chrom_sizes": {**_S_STRING, "description": "chrom sizes"},
            "out": {**_S_STRING, "description": "out"},
            "allow_missing_barcodes": {**_S_BOOLEAN, "default": False, "description": "legacy branches: warn instead of failing"},
            "lenient_fields": {**_S_BOOLEAN, "default": False, "description": "legacy branches: skip records with < 4 fields"},
            "clean": {**_S_BOOLEAN, "default": False, "description": "write the 16-bp barcode only (legacy option)"}},
         "required": ["qc_guide", "pseudobulks", "cell_type", "chrom_sizes", "out"]},
        cli=["e2g-qc-predictions", "filter-atac"],
        flag_map={"qc_guide": "--qc-guide", "pseudobulks": "--pseudobulks", "cell_type": "--cell-type", "chrom_sizes": "--chrom-sizes", "out": "--out", "allow_missing_barcodes": "--allow-missing-barcodes", "lenient_fields": "--lenient-fields", "clean": "--clean"},
        bool_flags={'clean', 'allow_missing_barcodes', 'lenient_fields'},
    ),

    _T(
        "e2g_qc_filter_rna",
        "Filters and concatenates one cluster's rna_counts_mtx.h5ad files to the QC guide barcodes, collapses Ensembl IDs to gene symbols by summing (embedded gene_symbol or an exact versioned GTF match, optional standard chromosomes only), runs the four sanity checks and writes an .mtx directory, .h5ad or .csv.gz; optionally the Portal tar package.",
        {"type": "object", "properties": {
            "qc_guide": {**_S_STRING, "description": "qc guide"},
            "pseudobulks": {**_S_STRING, "description": "pseudobulks"},
            "cell_type": {**_S_STRING, "description": "cell type"},
            "out": {**_S_STRING, "description": "out"},
            "ensembl_ids_as_genes": {**_S_BOOLEAN, "default": False, "description": "ensembl ids as genes"},
            "gtf": {**_S_STRING, "description": "gtf"},
            "standard_chromosomes_only": {**_S_BOOLEAN, "default": False, "description": "standard chromosomes only"},
            "log": {**_S_STRING, "description": "log"},
            "package_tar": {**_S_STRING, "description": "also write the Filtered Matrix File tar.gz here (mtx output only)"}},
         "required": ["qc_guide", "pseudobulks", "cell_type", "out"]},
        cli=["e2g-qc-predictions", "filter-rna"],
        flag_map={"qc_guide": "--qc-guide", "pseudobulks": "--pseudobulks", "cell_type": "--cell-type", "out": "--out", "ensembl_ids_as_genes": "--ensembl-ids-as-genes", "gtf": "--gtf", "standard_chromosomes_only": "--standard-chromosomes-only", "log": "--log", "package_tar": "--package-tar"},
        bool_flags={'ensembl_ids_as_genes', 'standard_chromosomes_only'},
    ),

    _T(
        "e2g_qc_package_rna",
        "Packages an RNA count matrix directory (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz) as the flat tar.gz of decompressed files the IGVF Filtered Matrix File spec requires. Writes the tarball.",
        {"type": "object", "properties": {
            "matrix_dir": {**_S_STRING, "description": "matrix dir"},
            "out": {**_S_STRING, "description": "out"}},
         "required": ["matrix_dir", "out"]},
        cli=["e2g-qc-predictions", "package-rna"],
        flag_map={"matrix_dir": "--matrix-dir", "out": "--out"},
    ),

    _T(
        "e2g_qc_sce2g_config",
        "Writes the per-dataset scE2G cell_clusters table and cluster metadata table (cell type / ontology id joined from a lab annotation table), the scE2G config overlay JSON and every cluster-model score threshold; prints the snakemake command upstream would run. Writes TSVs, JSON, report.md and summary.json.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "sce2g_config", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "config": {**_S_STRING, "description": "*_pipeline_config.yaml (or .json)"},
            "output_dir": {**_S_STRING, "description": "override the config's output_dir"},
            "lab_annotations": {**_S_STRING, "description": "lab_annotations_with_cl.tsv (dataset, lab_celltype, CL term, qualifier, CL_ID)"}},
         "required": ["config"]},
        cli=["e2g-qc-predictions", "sce2g-config"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "config": "--config", "output_dir": "--output-dir", "lab_annotations": "--lab-annotations"},
    ),

    _T(
        "e2g_qc_reformat",
        "Reformats scE2G outputs into the IGVF Portal sharing format: full or thresholded predictions (consortium column order, portal metadata header with SampleTermName / SampleTermID / CellAnnotation / ScoreThreshold / Metadata link), gene and element lists, the EnhancerList element BED (bgzip + tabix) or the thresholded bedpe (sorted, bgzip + tabix). --format synapse-legacy and --catlas reproduce the other branches.",
        {"type": "object", "properties": {
            "kind": {**_S_STRING, "default": "auto", "description": "auto: full / thresholded / gene list / element list from the input name and --threshold (choices: auto, element-bed, bedpe)"},
            "input": {**_S_STRING, "description": "input"},
            "out": {**_S_STRING, "description": "out"},
            "model": {**_S_STRING, "default": "multiome_powerlaw_v3", "description": "model"},
            "method": {**_S_STRING, "description": "default scE2G_<model>"},
            "version": {**_S_STRING, "default": "1.2", "description": "version"},
            "cell_type": {**_S_STRING, "description": "SampleTermName (portal cell_type.term_name)"},
            "term_id": {**_S_STRING, "description": "SampleTermID (CURIE)"},
            "summary": {**_S_STRING, "description": "CellAnnotation (SampleSummaryShort in --format synapse-legacy)"},
            "threshold": {**_S_STRING, "description": "score threshold; header 'ScoreThreshold: Score >= t'"},
            "portal_link": {**_S_STRING, "description": "the file's IGVF alias"},
            "all_columns": {**_S_BOOLEAN, "default": False, "description": "all columns"},
            "format": {**_S_STRING, "default": "portal", "description": "format (choices: portal, synapse-legacy)"},
            "catlas": {**_S_BOOLEAN, "default": False, "description": "CATlas branch: scATAC full files without normalizedATAC_enh"}},
         "required": ["input", "out"]},
        cli=["e2g-qc-predictions", "reformat"],
        flag_map={"kind": "--kind", "input": "--input", "out": "--out", "model": "--model", "method": "--method", "version": "--version", "cell_type": "--cell-type", "term_id": "--term-id", "summary": "--summary", "threshold": "--threshold", "portal_link": "--portal-link", "all_columns": "--all-columns", "format": "--format", "catlas": "--catlas"},
        bool_flags={'catlas', 'all_columns'},
    ),

    _T(
        "e2g_qc_candidates",
        "Extracts the candidate enhancer-gene pairs table (element coordinates, gene, SampleSummaryShort, E2G_Distance, isSelfPromoter) from an scE2G predictions table. Writes a gzipped TSV.",
        {"type": "object", "properties": {
            "input": {**_S_STRING, "description": "input"},
            "out": {**_S_STRING, "description": "out"},
            "summary": {**_S_STRING, "description": "summary"}},
         "required": ["input", "out", "summary"]},
        cli=["e2g-qc-predictions", "candidates"],
        flag_map={"input": "--input", "out": "--out", "summary": "--summary"},
    ),

    _T(
        "e2g_qc_features",
        "Writes the scE2G feature table for sharing: '# Source: scE2G <model>' header block, fixed leading element / gene columns and every remaining feature column. Writes a gzipped TSV.",
        {"type": "object", "properties": {
            "input": {**_S_STRING, "description": "input"},
            "out": {**_S_STRING, "description": "out"},
            "model": {**_S_STRING, "description": "model"},
            "cell_type": {**_S_STRING, "description": "cell type"},
            "term_id": {**_S_STRING, "description": "term id"},
            "summary": {**_S_STRING, "description": "summary"}},
         "required": ["input", "out", "model"]},
        cli=["e2g-qc-predictions", "features"],
        flag_map={"input": "--input", "out": "--out", "model": "--model", "cell_type": "--cell-type", "term_id": "--term-id", "summary": "--summary"},
    ),

    _T(
        "e2g_qc_aggregate_qc",
        "Rebuilds a dataset's all_qc_stats.tsv from the newest scE2G_predictions_threshold*_stats.tsv of every cluster directory on disk (cell_count patched from metrics) and draws the QC plots (enhancer, E-G, gene, distance metrics vs depth with 2e6 / 100 / 1e6 lines). Writes the TSV, PNG figures, report.md and summary.json.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "aggregate_qc", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "dataset_dir": {**_S_STRING, "description": "dataset dir"},
            "plots_dir": {**_S_STRING, "description": "plots/<dataset> (metrics files for cell_count patching)"},
            "out": {**_S_STRING, "description": "out"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "no plots"}},
         "required": ["dataset_dir", "plots_dir"]},
        cli=["e2g-qc-predictions", "aggregate-qc"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "dataset_dir": "--dataset-dir", "plots_dir": "--plots-dir", "out": "--out", "no_plots": "--no-plots"},
        bool_flags={'no_plots'},
    ),

    _T(
        "e2g_qc_stale_reformats",
        "Finds reformatted files whose embedded CellAnnotation / SampleTermName / SampleTermID header no longer matches the Cell Annotation cache and optionally deletes them (with their .tbi) so they are regenerated.",
        {"type": "object", "properties": {
            "results_dir": {**_S_STRING, "description": "<output_dir>/uniformly_processed"},
            "datasets": {**_S_ARRAY_S, "description": "datasets"},
            "state_db": {**_S_STRING, "default": "/Users/hzhou/Research/Projects/Compute/IGVFagent/Docs/E2GQCPredictions/igvf_metadata_state.db", "description": "state db"},
            "annotations_tsv": {**_S_STRING, "description": "annotations tsv"},
            "delete": {**_S_BOOLEAN, "default": False, "description": "delete"}},
         "required": ["results_dir", "datasets"]},
        cli=["e2g-qc-predictions", "stale-reformats"],
        flag_map={"results_dir": "--results-dir", "datasets": "--datasets", "state_db": "--state-db", "annotations_tsv": "--annotations-tsv", "delete": "--delete"},
        flag_repeat={'datasets'},
        bool_flags={'delete'},
    ),

    _T(
        "e2g_qc_cell_metadata",
        "Builds the Cell Annotation cache from one IGVF Portal PseudobulkSet multireport (live authenticated GET or a saved JSON): classifies primary vs principal pseudobulks, caches every single-sample primary, derives each cluster's annotation, sample term and qualifier from its QC-guide subsamples, and writes the digest-checked snapshot and a per-cluster status TSV. --catlas-seed seeds the WashU CATlas clusters.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "cell_metadata", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "config": {**_S_STRING, "description": "*_pipeline_config.yaml (or .json)"},
            "output_dir": {**_S_STRING, "description": "override the config's output_dir"},
            "multireport_json": {**_S_STRING, "description": "saved PseudobulkSet multireport (offline)"},
            "state_db": {**_S_STRING, "default": "/Users/hzhou/Research/Projects/Compute/IGVFagent/Docs/E2GQCPredictions/igvf_metadata_state.db", "description": "state db"},
            "offline": {**_S_BOOLEAN, "default": False, "description": "never GET; derive from the cache"},
            "ttl_hours": {**_S_NUMBER, "default": 24, "description": "ttl hours"},
            "write_snapshot": {**_S_BOOLEAN, "default": False, "description": "write the snapshot under output_dir/igvf_metadata"},
            "catlas_seed": {**_S_BOOLEAN, "default": False, "description": "seed CATlas (WashU) clusters from the lab-filtered multireport"},
            "dataset": {**_S_STRING, "description": "dataset label for --catlas-seed (default catlas)"}}},
        cli=["e2g-qc-predictions", "cell-metadata"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "config": "--config", "output_dir": "--output-dir", "multireport_json": "--multireport-json", "state_db": "--state-db", "offline": "--offline", "ttl_hours": "--ttl-hours", "write_snapshot": "--write-snapshot", "catlas_seed": "--catlas-seed", "dataset": "--dataset"},
        bool_flags={'write_snapshot', 'offline', 'catlas_seed'},
    ),

    _T(
        "e2g_qc_cell_annotation_report",
        "Writes the IGVF cell annotation report: one row per dataset-cluster of primary pseudobulks with CellAnnotation, SampleTermID, SampleTermName, CellQualifier, subsamples and status; with --qc-guide-dir resolves a single value per cluster from the most-contributing subsample. Input: live multireport or saved JSON.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "cell_annotation_report", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "multireport_json": {**_S_STRING, "description": "multireport json"},
            "qc_guide_dir": {**_S_STRING, "description": "qc guide dir"},
            "out": {**_S_STRING, "description": "out"}}},
        cli=["e2g-qc-predictions", "cell-annotation-report"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "multireport_json": "--multireport-json", "qc_guide_dir": "--qc-guide-dir", "out": "--out"},
    ),

    _T(
        "e2g_qc_dataset_accessions",
        "Maps each informal dataset label (from primary pseudobulk aliases) to its principal analysis set accession(s) on the IGVF Portal and warns when a dataset maps to more than one. Writes the mapping JSON.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "dataset_accessions", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "multireport_json": {**_S_STRING, "description": "multireport json"},
            "json_out": {**_S_STRING, "description": "json out"}}},
        cli=["e2g-qc-predictions", "dataset-accessions"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "multireport_json": "--multireport-json", "json_out": "--json-out"},
    ),

    _T(
        "e2g_qc_washu_report",
        "Flattens the WashU (Yang Li lab, CATlas) PseudobulkSets into a TSV with annotations, sample terms and the fragments file alias / href / md5. Input: live lab-filtered multireport or saved JSON.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "washu_report", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "multireport_json": {**_S_STRING, "description": "multireport json"},
            "lab": {**_S_STRING, "default": "/labs/yang-li/", "description": "lab"},
            "out": {**_S_STRING, "description": "out"}}},
        cli=["e2g-qc-predictions", "washu-report"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "multireport_json": "--multireport-json", "lab": "--lab", "out": "--out"},
    ),

    _T(
        "e2g_qc_verify_fragments",
        "Verifies local fragments files listed in a WashU report against their Portal md5sum and, with --download, fetches missing or mismatched ones through the authenticated Portal download. Writes a per-file verification TSV.",
        {"type": "object", "properties": {
            "report": {**_S_STRING, "description": "report"},
            "outdir": {**_S_STRING, "description": "outdir"},
            "download": {**_S_BOOLEAN, "default": False, "description": "download"},
            "force": {**_S_BOOLEAN, "default": False, "description": "force"},
            "summary": {**_S_STRING, "description": "summary"}},
         "required": ["report", "outdir"]},
        cli=["e2g-qc-predictions", "verify-fragments"],
        flag_map={"report": "--report", "outdir": "--outdir", "download": "--download", "force": "--force", "summary": "--summary"},
        bool_flags={'force', 'download'},
    ),

    _T(
        "e2g_qc_portal_files",
        "Discovers the primary-pseudobulk fragments / RNA matrix / per-cell QC files on the IGVF Portal, resolves where each belongs on disk from submitted_file_name (recovering irregular annotation_0 sets from the alias), flags review reasons and plans (or with --download performs, md5-verified) the download. Writes portal_files.tsv, report.md and summary.json.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "portal_files", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "multireport_json": {**_S_STRING, "description": "multireport json"},
            "lab": {**_S_STRING, "default": "/labs/anshul-kundaje/", "description": "lab"},
            "content_types": {**_S_ARRAY_S, "description": "content types"},
            "datasets": {**_S_ARRAY_S, "description": "datasets"},
            "download_root": {**_S_STRING, "description": "download root"},
            "download": {**_S_BOOLEAN, "default": False, "description": "fetch files marked todo (md5-verified)"}}},
        cli=["e2g-qc-predictions", "portal-files"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "multireport_json": "--multireport-json", "lab": "--lab", "content_types": "--content-types", "datasets": "--datasets", "download_root": "--download-root", "download": "--download"},
        flag_repeat={'content_types', 'datasets'},
        bool_flags={'download'},
    ),

    _T(
        "e2g_qc_compare_archive",
        "Compares Portal-downloaded pseudobulk files with an existing archive by raw md5, decompressed-content md5, TSV and h5ad structure, and rolls up a per-cluster release decision (shareable, changed, novel). Writes per-file and per-cluster TSVs.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "compare_archive", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "discovery": {**_S_STRING, "description": "portal_files.tsv from portal-files"},
            "archive_root": {**_S_STRING, "description": "archive root"},
            "no_characterize": {**_S_BOOLEAN, "default": False, "description": "no characterize"}},
         "required": ["discovery", "archive_root"]},
        cli=["e2g-qc-predictions", "compare-archive"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "discovery": "--discovery", "archive_root": "--archive-root", "no_characterize": "--no-characterize"},
        bool_flags={'no_characterize'},
    ),

    _T(
        "e2g_qc_manifest",
        "Builds the IGVF Portal submission for the chosen clusters: eleven metadata tables / fifteen variants (principal pseudobulk set, prediction set, tabular / matrix / index / signal files, QC document) with aliases, derived_from / file_set / input_file_sets links, controlled vocabulary and required-field validation; dependency rounds, iu_register TSVs, REST payloads, upload_plan.json and an optional lineage check of every external reference. DRY-RUN unless execute=true, which POSTs/PATCHes in round order to the IGVF sandbox (production only with production=true) using IGVF_ACCESS_KEY / IGVF_SECRET_ACCESS_KEY.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "manifest", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "config": {**_S_STRING, "description": "*_pipeline_config.yaml (or .json)"},
            "output_dir": {**_S_STRING, "description": "override the config's output_dir"},
            "cluster_keys": {**_S_STRING, "description": "'dataset/cluster,...' or a bare dataset (manifest-eligible expansion)"},
            "excluded_cluster_keys": {**_S_STRING, "default": "", "description": "excluded cluster keys"},
            "state_db": {**_S_STRING, "description": "state db"},
            "manifest_dir": {**_S_STRING, "description": "manifest dir"},
            "tables": {**_S_ARRAY_S, "description": "tables"},
            "offline": {**_S_BOOLEAN, "default": False, "description": "no alias lookups (every row not in the ledger is a POST)"},
            "live_aliases": {**_S_STRING, "description": "JSON {alias: record} treated as already live (offline PATCH planning)"},
            "lineage": {**_S_BOOLEAN, "default": False, "description": "resolve every external reference with an authenticated GET"},
            "lineage_walk": {**_S_BOOLEAN, "default": False, "description": "also walk the Portal graph from igvf.principal_analysis_sets"},
            "lineage_fixture": {**_S_STRING, "description": "JSON {path: object} answering lineage GETs offline"},
            "execute": {**_S_BOOLEAN, "default": False, "description": "POST/PATCH for real (credentials from IGVF_ACCESS_KEY / IGVF_SECRET_ACCESS_KEY)"},
            "production": {**_S_BOOLEAN, "default": False, "description": "target api.data.igvf.org instead of the sandbox"}},
         "required": ["config", "cluster_keys"]},
        cli=["e2g-qc-predictions", "manifest"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "config": "--config", "output_dir": "--output-dir", "cluster_keys": "--cluster-keys", "excluded_cluster_keys": "--excluded-cluster-keys", "state_db": "--state-db", "manifest_dir": "--manifest-dir", "tables": "--tables", "offline": "--offline", "live_aliases": "--live-aliases", "lineage": "--lineage", "lineage_walk": "--lineage-walk", "lineage_fixture": "--lineage-fixture", "execute": "--execute", "production": "--production"},
        flag_repeat={'tables'},
        bool_flags={'execute', 'lineage_walk', 'lineage', 'offline', 'production'},
    ),

    _T(
        "e2g_qc_patch_submitter_comment",
        "Plans the 'Version 1' submitter_comment backfill for live Prediction Sets (set when empty, prefix once when free text exists, never stack) and writes the PATCH TSV; patches only with execute=true.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "submitter_comment", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "multireport_json": {**_S_STRING, "description": "multireport json"},
            "dataset": {**_S_ARRAY_S, "description": "dataset"},
            "execute": {**_S_BOOLEAN, "default": False, "description": "execute"},
            "production": {**_S_BOOLEAN, "default": False, "description": "production"}}},
        cli=["e2g-qc-predictions", "patch-submitter-comment"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "multireport_json": "--multireport-json", "dataset": "--dataset", "execute": "--execute", "production": "--production"},
        flag_repeat={'dataset'},
        bool_flags={'execute', 'production'},
    ),

    _T(
        "e2g_qc_report",
        "Writes the per-cluster coverage report.tsv: quality gate result, Cell Annotation availability, predictions / reformat status and the manifest roll-up (ready, blocked, excluded, not-generated) from each dataset's manifest_coverage.tsv.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "coverage_report", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "config": {**_S_STRING, "description": "*_pipeline_config.yaml (or .json)"},
            "output_dir": {**_S_STRING, "description": "override the config's output_dir"},
            "state_db": {**_S_STRING, "default": "/Users/hzhou/Research/Projects/Compute/IGVFagent/Docs/E2GQCPredictions/igvf_metadata_state.db", "description": "state db"},
            "annotations_tsv": {**_S_STRING, "description": "annotations tsv"},
            "manifest_dir": {**_S_STRING, "description": "manifest dir"},
            "out": {**_S_STRING, "description": "out"}},
         "required": ["config"]},
        cli=["e2g-qc-predictions", "report"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "config": "--config", "output_dir": "--output-dir", "state_db": "--state-db", "annotations_tsv": "--annotations-tsv", "manifest_dir": "--manifest-dir", "out": "--out"},
    ),

    _T(
        "e2g_qc_synapse_manifest",
        "Diffs the files a Synapse product (filtered_data, predictions, candidates, features) should hold against an inventory of owned entities: new, overwrite and stale rows (other models' files invisible, longest cluster-name match), preserving other runs' manifest rows. Writes the path/parent manifest; syncs only with execute=true.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "synapse", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "product": {**_S_STRING, "description": "product (choices: candidates, features, filtered_data, predictions)"},
            "cluster_keys": {**_S_STRING, "description": "cluster keys"},
            "files": {**_S_ARRAY_S, "description": "files"},
            "inventory": {**_S_STRING, "description": "TSV of owned Synapse entities (dataset, cluster, name, id, modifiedBy, modifiedOn)"},
            "manifest_out": {**_S_STRING, "description": "manifest out"},
            "parent_id": {**_S_STRING, "description": "parent id"},
            "legacy_first_match": {**_S_BOOLEAN, "default": False, "description": "legacy first match"},
            "execute": {**_S_BOOLEAN, "default": False, "description": "execute"},
            "confirm_delete": {**_S_BOOLEAN, "default": False, "description": "confirm delete"},
            "confirm_overwrite": {**_S_BOOLEAN, "default": False, "description": "confirm overwrite"}},
         "required": ["product", "cluster_keys"]},
        cli=["e2g-qc-predictions", "synapse-manifest"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "product": "--product", "cluster_keys": "--cluster-keys", "files": "--files", "inventory": "--inventory", "manifest_out": "--manifest-out", "parent_id": "--parent-id", "legacy_first_match": "--legacy-first-match", "execute": "--execute", "confirm_delete": "--confirm-delete", "confirm_overwrite": "--confirm-overwrite"},
        flag_repeat={'files'},
        bool_flags={'legacy_first_match', 'execute', 'confirm_delete', 'confirm_overwrite'},
    ),

    _T(
        "e2g_qc_synapse_orphans",
        "Read-only listing of Synapse folders (cluster mode, against report.tsv) or files (file mode, against a predictions manifest) that have no in-scope replacement. Writes synapse_orphans.tsv.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "synapse_orphans", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "mode": {**_S_STRING, "default": "cluster", "description": "mode (choices: cluster, file)"},
            "inventory": {**_S_STRING, "description": "inventory"},
            "in_scope_tsv": {**_S_STRING, "description": "report.tsv (cluster mode)"},
            "manifest": {**_S_STRING, "description": "predictions manifest TSV (file mode)"},
            "exclude_top": {**_S_ARRAY_S, "description": "exclude top"},
            "only": {**_S_ARRAY_S, "description": "only"},
            "user_id": {**_S_STRING, "description": "user id"}},
         "required": ["inventory"]},
        cli=["e2g-qc-predictions", "synapse-orphans"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "mode": "--mode", "inventory": "--inventory", "in_scope_tsv": "--in-scope-tsv", "manifest": "--manifest", "exclude_top": "--exclude-top", "only": "--only", "user_id": "--user-id"},
        flag_repeat={'only', 'exclude_top'},
    ),

    _T(
        "e2g_qc_distance_depth",
        "Relates distance to TSS of predicted enhancer-gene links to sequencing depth (CATlas analysis): per cluster mean / median / IQR distance vs total fragments and fragments per cell, a log10 trend fitted over clusters above 2e6 fragments, neuron vs non-neuron groups. Writes a TSV, PNG figures, report.md and summary.json.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "distance_depth", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "stats": {**_S_ARRAY_S, "description": "all_qc_stats.tsv file(s)"},
            "predictions_map": {**_S_STRING, "description": "TSV cluster, path (thresholded predictions) for median / IQR"},
            "group_map": {**_S_STRING, "description": "TSV cluster, group"},
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "no plots"}},
         "required": ["stats"]},
        cli=["e2g-qc-predictions", "distance-depth"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "stats": "--stats", "predictions_map": "--predictions-map", "group_map": "--group-map", "no_plots": "--no-plots"},
        flag_repeat={'stats'},
        bool_flags={'no_plots'},
    ),

    _T(
        "e2g_qc_run",
        "Runs the QC-and-Predictions driver for a pipeline config: preflight quality gate, Cell Annotation warm-up (snapshot), local packaging of existing scE2G outputs (portal reformatting, element BED, bedpe index, RNA tarball, candidates, features, QC stats), a dry-run IGVF manifest and the audit report.tsv; exit 0 only when every manifest-eligible cluster is complete. Never submits.",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "default": "run", "description": "run label (output dir suffix)"},
            "out_dir": {**_S_STRING, "description": "write here instead of Docs/E2GQCPredictions/<ts>_<label>"},
            "config": {**_S_STRING, "description": "*_pipeline_config.yaml (or .json)"},
            "output_dir": {**_S_STRING, "description": "override the config's output_dir"},
            "mode": {**_S_STRING, "description": "mode (choices: default, local_only)"},
            "multireport_json": {**_S_STRING, "description": "multireport json"},
            "state_db": {**_S_STRING, "description": "state db"},
            "offline": {**_S_BOOLEAN, "default": False, "description": "offline"},
            "dry_run": {**_S_BOOLEAN, "default": False, "description": "dry run"},
            "legacy_implicit_source": {**_S_BOOLEAN, "default": False, "description": "legacy implicit source"},
            "lab_annotations": {**_S_STRING, "description": "lab annotations"}},
         "required": ["config"]},
        cli=["e2g-qc-predictions", "run"],
        flag_map={"label": "--label", "out_dir": "--out-dir", "config": "--config", "output_dir": "--output-dir", "mode": "--mode", "multireport_json": "--multireport-json", "state_db": "--state-db", "offline": "--offline", "dry_run": "--dry-run", "legacy_implicit_source": "--legacy-implicit-source", "lab_annotations": "--lab-annotations"},
        bool_flags={'offline', 'dry_run', 'legacy_implicit_source'},
    ),

    _T(
        "e2g_qc_selftest",
        "Runs the synthetic end-to-end self-test of every e2g-qc-predictions subcommand (planted QC failures, annotation disagreements, merged clusters, manifest rounds, in-memory execute) and reports pass / fail.",
        {"type": "object", "properties": {
            "no_plots": {**_S_BOOLEAN, "default": False, "description": "no plots"}}},
        cli=["e2g-qc-predictions", "selftest"],
        flag_map={"no_plots": "--no-plots"},
        bool_flags={'no_plots'},
    ),

    _T(
        "alphagenome_setup",
        "Checks that AlphaGenome can be used here: whether the SDK is importable (it needs Python 3.10+), whether a 3.10+ interpreter with the SDK exists, and whether an API key is configured (never printed). install=true creates that environment; ping=true makes one API call. Run this first when an AlphaGenome tool reports the SDK or key missing.",
        {"type": "object", "properties": {
            "install": {**_S_BOOLEAN, "description": "Create a Python 3.10+ environment with the SDK if missing."},
            "ping": {**_S_BOOLEAN, "description": "Call the API once to confirm the key works."}}},
        cli=["alphagenome", "setup"],
        bool_flags={'install', 'ping'},
    ),

    _T(
        "alphagenome_metadata",
        "Lists every track AlphaGenome predicts (output type, assay, ontology CURIE, biosample) so you can pick ontology terms for the other AlphaGenome tools. search filters by tissue/cell-type text (e.g. 'heart', 'hepatocyte'). Writes tracks.tsv, biosamples.tsv and a report.",
        {"type": "object", "properties": {
            "search": {**_S_STRING, "description": "Substring of biosample name, CURIE or track name."},
            "ontology": {**_S_ARRAY_S, "description": "Ontology CURIEs to restrict tracks, e.g. UBERON:0000948 (heart), UBERON:0002107 (liver), CL:0000746 (cardiac muscle cell). Find them with alphagenome_metadata."},
            "output_type": {**_S_ARRAY_S, "description": "Restrict to output types (RNA_SEQ, DNASE, ATAC, CAGE, CHIP_HISTONE, CHIP_TF, SPLICE_*, CONTACT_MAPS, PROCAP)."},
            "organism": {**_S_STRING, "default": "HOMO_SAPIENS", "description": "HOMO_SAPIENS (hg38, default) or MUS_MUSCULUS (mm10)."},
            "label": {**_S_STRING, "description": "Short run label."}}},
        cli=["alphagenome", "metadata"],
        flag_map={"output_type": "--output-type"},
        flag_repeat={'output_type', 'ontology'},
    ),

    _T(
        "alphagenome_predict_interval",
        "Predicts AlphaGenome tracks (expression, accessibility, histone and TF ChIP, splicing, contact maps) over a region or gene, in chosen tissues/cell types. The region is centred in a supported window. Writes a per-track summary TSV (mean/max/sum over the region), a track figure and a report. ontology is required unless all_tracks.",
        {"type": "object", "properties": {
            "interval": {**_S_STRING, "description": "Region chr:start-end, 0-based half-open (BED convention)."},
            "gene": {**_S_STRING, "description": "Gene symbol instead of an interval; coordinates come from the IGVF Catalog."},
            "outputs": {**_S_ARRAY_S, "description": "Output types to predict, e.g. RNA_SEQ, DNASE."},
            "ontology": {**_S_ARRAY_S, "description": "Ontology CURIEs to restrict tracks, e.g. UBERON:0000948 (heart), UBERON:0002107 (liver), CL:0000746 (cardiac muscle cell). Find them with alphagenome_metadata."},
            "length": {**_S_STRING, "default": "1MB", "description": "Prediction window: 16KB, 100KB, 500KB or 1MB (default 1MB, the model's full context)."},
            "organism": {**_S_STRING, "default": "HOMO_SAPIENS", "description": "HOMO_SAPIENS (hg38, default) or MUS_MUSCULUS (mm10)."},
            "all_tracks": {**_S_BOOLEAN, "description": "Return every tissue (large)."},
            "label": {**_S_STRING, "description": "Short run label."}},
         "required": ["outputs"]},
        cli=["alphagenome", "predict-interval"],
        flag_map={"all_tracks": "--all-tracks"},
        flag_repeat={'outputs', 'ontology'},
        bool_flags={'all_tracks'},
    ),

    _T(
        "alphagenome_predict_variant",
        "Predicts REF and ALT tracks for one variant (rsID, SPDI, chr-pos-ref-alt or chr:pos:ref>alt; rsIDs resolved through the IGVF Catalog) and ranks tracks by the ALT vs REF change near the variant. Writes ref_alt_summary.tsv, a REF/ALT overlay figure and a report.",
        {"type": "object", "properties": {
            "variant": {**_S_STRING, "description": "The variant, any notation."},
            "outputs": {**_S_ARRAY_S, "description": "Output types, e.g. RNA_SEQ, DNASE, SPLICE_SITES."},
            "ontology": {**_S_ARRAY_S, "description": "Ontology CURIEs to restrict tracks, e.g. UBERON:0000948 (heart), UBERON:0002107 (liver), CL:0000746 (cardiac muscle cell). Find them with alphagenome_metadata."},
            "radius": {**_S_INTEGER, "default": 1000, "description": "Summarise within +/- this many bp of the variant (default 1000)."},
            "length": {**_S_STRING, "default": "1MB", "description": "Prediction window: 16KB, 100KB, 500KB or 1MB (default 1MB, the model's full context)."},
            "organism": {**_S_STRING, "default": "HOMO_SAPIENS", "description": "HOMO_SAPIENS (hg38, default) or MUS_MUSCULUS (mm10)."},
            "label": {**_S_STRING, "description": "Short run label."}},
         "required": ["variant", "outputs"]},
        cli=["alphagenome", "predict-variant"],
        flag_repeat={'outputs', 'ontology'},
    ),

    _T(
        "alphagenome_score_variants",
        "Scores variant effects with AlphaGenome's recommended variant scorers (or chosen ones: DNASE, ATAC, CHIP_TF, CHIP_HISTONE, CAGE, PROCAP, RNA_SEQ, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, POLYADENYLATION, CONTACT_MAPS and *_ACTIVE). Variants in any notation, or a file/VCF; rsIDs resolved through the IGVF Catalog. Writes the tidy score table (raw and quantile scores per track and gene), a per-variant summary and a report of the strongest effects. Use ontology or search to focus on a tissue.",
        {"type": "object", "properties": {
            "variants": {**_S_ARRAY_S, "description": "Variants, any notation."},
            "input": {**_S_STRING, "description": "Path to a file of variants or a VCF."},
            "scorers": {**_S_ARRAY_S, "description": "Scorer names; default the recommended set."},
            "ontology": {**_S_ARRAY_S, "description": "Ontology CURIEs to restrict tracks, e.g. UBERON:0000948 (heart), UBERON:0002107 (liver), CL:0000746 (cardiac muscle cell). Find them with alphagenome_metadata."},
            "search": {**_S_STRING, "description": "Keep rows whose biosample/track/tissue contains this text."},
            "length": {**_S_STRING, "default": "1MB", "description": "Prediction window: 16KB, 100KB, 500KB or 1MB (default 1MB, the model's full context)."},
            "organism": {**_S_STRING, "default": "HOMO_SAPIENS", "description": "HOMO_SAPIENS (hg38, default) or MUS_MUSCULUS (mm10)."},
            "max_variants": {**_S_INTEGER, "default": 1000, "description": "Safety cap (default 1000)."},
            "label": {**_S_STRING, "description": "Short run label."}}},
        cli=["alphagenome", "score-variants"],
        flag_map={"max_variants": "--max-variants"},
        flag_repeat={'ontology', 'scorers', 'variants'},
    ),

    _T(
        "alphagenome_score_interval",
        "Scores a region or gene with AlphaGenome's recommended interval scorers (gene-level expression and other summaries), optionally focused on tissues. Writes interval_scores.tsv and a report.",
        {"type": "object", "properties": {
            "interval": {**_S_STRING, "description": "Region chr:start-end, 0-based half-open (BED convention)."},
            "gene": {**_S_STRING, "description": "Gene symbol instead of an interval; coordinates come from the IGVF Catalog."},
            "ontology": {**_S_ARRAY_S, "description": "Ontology CURIEs to restrict tracks, e.g. UBERON:0000948 (heart), UBERON:0002107 (liver), CL:0000746 (cardiac muscle cell). Find them with alphagenome_metadata."},
            "search": {**_S_STRING, "description": "Tissue/biosample text filter."},
            "length": {**_S_STRING, "default": "1MB", "description": "Prediction window: 16KB, 100KB, 500KB or 1MB (default 1MB, the model's full context)."},
            "organism": {**_S_STRING, "default": "HOMO_SAPIENS", "description": "HOMO_SAPIENS (hg38, default) or MUS_MUSCULUS (mm10)."},
            "label": {**_S_STRING, "description": "Short run label."}}},
        cli=["alphagenome", "score-interval"],
        flag_repeat={'ontology'},
    ),

    _T(
        "alphagenome_ism",
        "In silico mutagenesis: scores every alternative base at every position of a short interval (3 variants per bp, capped at max_bp) with AlphaGenome and reports the most sensitive positions. Writes ism_scores.tsv, a position x base matrix, a heatmap and a report.",
        {"type": "object", "properties": {
            "ism_interval": {**_S_STRING, "description": "chr:start-end to mutate, 0-based half-open (keep it short)."},
            "scorer": {**_S_ARRAY_S, "description": "Scorer name(s), e.g. DNASE; default recommended set."},
            "ontology": {**_S_ARRAY_S, "description": "Ontology CURIEs to restrict tracks, e.g. UBERON:0000948 (heart), UBERON:0002107 (liver), CL:0000746 (cardiac muscle cell). Find them with alphagenome_metadata."},
            "search": {**_S_STRING, "description": "Tissue text filter."},
            "length": {**_S_STRING, "default": "1MB", "description": "Prediction window: 16KB, 100KB, 500KB or 1MB (default 1MB, the model's full context)."},
            "max_bp": {**_S_INTEGER, "default": 100, "description": "Width cap (default 100 bp)."},
            "organism": {**_S_STRING, "default": "HOMO_SAPIENS", "description": "HOMO_SAPIENS (hg38, default) or MUS_MUSCULUS (mm10)."},
            "label": {**_S_STRING, "description": "Short run label."}},
         "required": ["ism_interval"]},
        cli=["alphagenome", "ism"],
        flag_map={"ism_interval": "--ism-interval", "max_bp": "--max-bp"},
        flag_repeat={'scorer', 'ontology'},
    ),

    _T(
        "alphagenome_atlas_scorers",
        "Lists the scorers available in the AlphaGenome Atlas, the pre-computed genome-wide variant-effect resource (including the AlphaGenome Variant Impact score).",
        {"type": "object", "properties": {
            "label": {**_S_STRING, "description": "Short run label."}}},
        cli=["alphagenome", "atlas-scorers"],
    ),

    _T(
        "alphagenome_atlas_variants",
        "Looks up pre-computed AlphaGenome Atlas scores (including AVI) for variants in any notation, without running the model, so it is faster and suits more variants than score_variants. Optional scorer, ontology and gene filters. Writes atlas_scores.tsv and a report.",
        {"type": "object", "properties": {
            "variants": {**_S_ARRAY_S, "description": "Variants, any notation."},
            "input": {**_S_STRING, "description": "File of variants or VCF."},
            "scorers": {**_S_ARRAY_S, "description": "Atlas scorer names; default all (see alphagenome_atlas_scorers)."},
            "ontology": {**_S_ARRAY_S, "description": "Ontology CURIEs to restrict tracks, e.g. UBERON:0000948 (heart), UBERON:0002107 (liver), CL:0000746 (cardiac muscle cell). Find them with alphagenome_metadata."},
            "genes": {**_S_ARRAY_S, "description": "Restrict gene-level scores to these gene symbols."},
            "label": {**_S_STRING, "description": "Short run label."}}},
        cli=["alphagenome", "atlas-variants"],
        flag_repeat={'ontology', 'genes', 'scorers', 'variants'},
    ),

    _T(
        "alphagenome_atlas_interval",
        "Retrieves pre-computed AlphaGenome Atlas scores for every possible variant in a region or gene (3 per bp, capped at max_bp). Writes atlas_scores.tsv and a report of the strongest effects.",
        {"type": "object", "properties": {
            "interval": {**_S_STRING, "description": "Region chr:start-end, 0-based half-open (BED convention)."},
            "gene": {**_S_STRING, "description": "Gene symbol instead of an interval; coordinates come from the IGVF Catalog."},
            "scorers": {**_S_ARRAY_S, "description": "Atlas scorer names; default all."},
            "ontology": {**_S_ARRAY_S, "description": "Ontology CURIEs to restrict tracks, e.g. UBERON:0000948 (heart), UBERON:0002107 (liver), CL:0000746 (cardiac muscle cell). Find them with alphagenome_metadata."},
            "genes": {**_S_ARRAY_S, "description": "Gene symbols for gene-level scores."},
            "max_bp": {**_S_INTEGER, "default": 20000, "description": "Width cap (default 20000 bp)."},
            "label": {**_S_STRING, "description": "Short run label."}}},
        cli=["alphagenome", "atlas-interval"],
        flag_map={"max_bp": "--max-bp"},
        flag_repeat={'genes', 'scorers', 'ontology'},
    ),

    _T(
        "alphagenome_selftest",
        "Runs the AlphaGenome skill's offline self-test (notation parsing, coordinates, credentials handling and every subcommand against a fake backend). Needs no key and no network.",
        {"type": "object", "properties": {
}},
        cli=["alphagenome", "selftest"],
    ),

    _T(
        "sce2g_selftest",
        "Self-test of the scE2G workbench on a fake checkout and synthetic "
        "feature / prediction / CRISPR files: patches, config rows, checks and "
        "benchmark all asserted. Run before touching a real checkout.",
        {"type": "object", "properties": {"no_plots": {**_S_BOOLEAN, "default": False}}},
        cli=["sce2g", "selftest"], bool_flags={"no_plots"},
    ),

    _T(
        "sce2g_predict",
        "★ ENHANCER→GENE LINKS FROM PAIRED SINGLE-CELL ATAC + RNA ★ "
        "(scE2G-style). Computes, for every peak-gene pair in a window: the "
        "KENDALL rank correlation between the peak's accessibility and the "
        "gene's expression across metacells, the ABC share, and distance. The "
        "correlation is the thing ABC cannot see — abc_score uses bulk "
        "activity and distance, so it cannot tell a peak that CO-VARIES with "
        "the gene from one that merely sits nearby and is busy. Requires the "
        "SAME cells in both matrices (multiome or matched). NOTE: these are "
        "FEATURES plus a transparent combined score, NOT the output of "
        "scE2G's trained model — that model's asset is its fitted weights, "
        "which are not reproduced here. To retrieve published scE2G links "
        "instead, use sce2g_kg_pull.",
        {
            "type": "object",
            "properties": {
                "rna":       {**_S_STRING, "description": "Gene x cell .h5ad."},
                "atac":      {**_S_STRING, "description":
                               "Peak x cell .h5ad, SAME cell barcodes."},
                "peaks":     {**_S_STRING, "description":
                               "Peak BED; names match the ATAC var_names."},
                "genes":     {**_S_STRING, "description":
                               "Gene TSS BED; names match the RNA var_names."},
                "metacells": {**_S_INTEGER, "default": 50,
                               "description": "Cells are pooled before "
                               "correlating — scATAC is near-binary and a "
                               "per-cell correlation mostly measures dropout."},
                "window":    {**_S_INTEGER, "default": 1000000},
                "label":     {**_S_STRING},
            },
            "required": ["rna", "atac", "peaks", "genes"],
        },
        cli=["sce2g-predict", "predict"],
        flag_map={"rna": "--rna", "atac": "--atac", "peaks": "--peaks",
                   "genes": "--genes", "metacells": "--metacells",
                   "window": "--window", "label": "--label"},
    ),

    _T(
        "sctransform_run",
        "★ SCTRANSFORM NORMALISATION ★ — variance stabilisation by regularised "
        "negative-binomial regression (Hafemeister & Satija 2019, Seurat's "
        "SCTransform). Use INSTEAD OF log-normalisation when clustering looks "
        "driven by sequencing depth: log-normalising assumes every gene "
        "scales with depth the same way and low-expressed genes do not, so "
        "depth leaks into the PCA. Returns Pearson residuals (NOT counts — do "
        "not feed them to a count model such as sc_crispr_de_test or "
        "mpra_activity; use them for PCA / clustering / HVG selection). "
        "Residual variance also ranks genes, replacing a log-normalised HVG "
        "list.",
        {
            "type": "object",
            "properties": {
                "input":     {**_S_STRING, "description": "Raw-count .h5ad."},
                "min_cells": {**_S_INTEGER, "default": 5},
                "bandwidth": {**_S_NUMBER, "default": 0.3,
                               "description": "Kernel width over log10(gene "
                               "mean) for the regularisation step."},
                "label":     {**_S_STRING},
            },
            "required": ["input"],
        },
        cli=["sctransform", "run"],
        flag_map={"input": "--input", "min_cells": "--min-cells",
                   "bandwidth": "--bandwidth", "label": "--label"},
    ),

    _T(
        "guide_map",
        "★ MAP FASTQ READS TO A GUIDE LIBRARY ★ with imperfect matching. Point "
        "it at any FASTQ and any guide-library CSV/TSV and get per-guide "
        "counts. Handles TWO different kinds of imperfect match, which are "
        "not interchangeable: `mismatches` is a Hamming budget for SEQUENCING "
        "ERROR, while `edit` (ABE or CBE) masks the editor's product so a "
        "base-EDITED read still matches its own guide without spending that "
        "budget. Use for a screen whose counts you want independently of the "
        "full pipelines — crispr_screen_analyze and base_editing_screen_"
        "analyze do this internally for a whole IGVF screen.",
        {
            "type": "object",
            "properties": {
                "fastq":      {**_S_STRING},
                "library":    {**_S_STRING, "description":
                                "Guide library CSV/TSV; sequence and id "
                                "columns auto-detected."},
                "seq_col":    {**_S_STRING},
                "id_col":     {**_S_STRING},
                "mismatches": {**_S_INTEGER, "default": 1,
                                "description": "Hamming budget for sequencing "
                                "error; 0 = exact only."},
                "edit":       {**_S_STRING, "description":
                                "ABE or CBE — mask the editor's product."},
                "max_reads":  {**_S_INTEGER},
                "label":      {**_S_STRING},
            },
            "required": ["fastq", "library"],
        },
        cli=["guide-map", "map"],
        flag_map={"fastq": "--fastq", "library": "--library",
                   "seq_col": "--seq-col", "id_col": "--id-col",
                   "mismatches": "--mismatches", "edit": "--edit",
                   "max_reads": "--max-reads", "label": "--label"},
    ),

    _T(
        "abc_score",
        "★ PREDICT WHICH ENHANCERS REGULATE WHICH GENES ★ — the Activity-by-"
        "Contact model (Fulco 2019 / Nasser 2021). GENERATES predictions from "
        "your own data, where enhancer_gene_overview and "
        "catalog_variant_enhancers only RETRIEVE predictions others computed. "
        "Needs candidate elements (BED), gene TSSs (BED) and an ATAC/DNase "
        "bigWig; H3K27ac is optional but makes activity much better. Hi-C is "
        "NOT required — without it, contact falls back to the genome-wide "
        "power law, which Nasser 2021 showed performs close to the Hi-C "
        "version. NOTE the score is a SHARE: every gene's predictions sum to "
        "1 across its neighbourhood, so a strong element among stronger "
        "neighbours scores low. It is not 'how active is this enhancer'.",
        {
            "type": "object",
            "properties": {
                "elements":  {**_S_STRING, "description":
                               "Candidate element BED (peaks/cCREs). This tool "
                               "does NOT call peaks."},
                "genes":     {**_S_STRING, "description":
                               "Gene/TSS BED, name in column 4, strand in 6."},
                "atac":      {**_S_STRING, "description": "ATAC/DNase bigWig."},
                "h3k27ac":   {**_S_STRING, "description": "H3K27ac bigWig."},
                "hic":       {**_S_STRING, "description":
                               "Optional contact file; power law if omitted."},
                "gamma":     {**_S_NUMBER, "default": -0.87,
                               "description": "Power-law exponent (Fulco 2019)."},
                "threshold": {**_S_NUMBER, "default": 0.02},
                "label":     {**_S_STRING},
            },
            "required": ["elements", "genes", "atac"],
        },
        cli=["abc", "score"],
        flag_map={"elements": "--elements", "genes": "--genes",
                   "atac": "--atac", "h3k27ac": "--h3k27ac", "hic": "--hic",
                   "gamma": "--gamma", "threshold": "--threshold",
                   "label": "--label"},
    ),

    _T(
        "crispr_surf_deconvolve",
        "★ TILING SCREEN → REGULATORY REGIONS ★ (CRISPR-SURF method). In a "
        "tiling screen each guide perturbs a WINDOW, not a point — Cas9 cuts "
        "locally, dCas9-KRAB spreads hundreds of bases — so one element makes "
        "every nearby guide look active and the signal is SMEARED. This "
        "deconvolves that smear back into discrete elements by solving "
        "y = A·beta with an L1 penalty, and calls significance against an "
        "empirical null built from the screen's own negative-control guides. "
        "USE THIS when per-guide scores give a broad blur and you need the "
        "actual element boundaries, or to separate two nearby elements. "
        "Emits a bedgraph, a per-bin table and significant_regions.csv.",
        {
            "type": "object",
            "properties": {
                "guides":     {**_S_STRING, "description":
                                "CSV/TSV: chrom,start[,stop],score column(s)"
                                "[,class]. `class` marks negative_control "
                                "guides — without them significance is "
                                "uncalibrated and the run says so."},
                "score_cols": {**_S_STRING, "description":
                                "Comma list of replicate score columns; "
                                "auto-detected when omitted."},
                "nuclease":   {**_S_STRING, "description":
                                "cas9 | cpf1 | crispri | crispra. Sets the "
                                "perturbation range: 20 bp for a nuclease, "
                                "250 bp for CRISPRi/a spreading."},
                "range":      {**_S_INTEGER, "description":
                                "Perturbation range in bp; overrides nuclease."},
                "bin_size":   {**_S_INTEGER, "default": 10},
                "lam":        {**_S_NUMBER, "default": 0.1,
                                "description": "L1 strength; higher = sparser."},
                "fdr":        {**_S_NUMBER, "default": 0.05},
                "label":      {**_S_STRING},
            },
            "required": ["guides"],
        },
        cli=["crispr-surf", "deconvolve"],
        flag_map={"guides": "--guides", "score_cols": "--score-cols",
                   "nuclease": "--nuclease", "range": "--range",
                   "bin_size": "--bin-size", "lam": "--lam", "fdr": "--fdr",
                   "label": "--label"},
    ),

    _T(
        "crispresso_analyze",
        "★ QUANTIFY GENOME-EDITING OUTCOMES FROM AMPLICON READS ★ — CRISPResso2 "
        "(Clement et al., Nat Biotechnol 2019). Aligns amplicon sequencing "
        "reads against a reference and quantifies indels, substitutions and "
        "HDR, including base-editing outcomes and the per-position "
        "substitution table a bystander analysis reads. USE THIS for "
        "\"what edits did this amplicon get\", \"quantify indels\", "
        "\"base-editing efficiency from FASTQ\", or bystander/reporter allele "
        "questions. It is a separate program run as a subprocess; pass "
        "`extra_args` for any flag not surfaced here.",
        {
            "type": "object",
            "properties": {
                "fastq_r1":   {**_S_STRING, "description": "Reads (fastq or fastq.gz)."},
                "fastq_r2":   {**_S_STRING, "description": "Optional mate."},
                "amplicon":   {**_S_STRING, "description":
                                "Reference amplicon SEQUENCE (not a path)."},
                "guide":      {**_S_STRING, "description":
                                "sgRNA spacer sequence, no PAM."},
                "base_editor": {**_S_BOOLEAN, "description":
                                 "Base-editor output: per-position "
                                 "substitution quantification."},
                "conversion": {**_S_STRING, "description":
                                "For base editors, e.g. 'A,G' for ABE or "
                                "'C,T' for CBE."},
                "name":       {**_S_STRING, "description": "Run name."},
                "output_dir": {**_S_STRING},
                "extra_args": {**_S_STRING, "description":
                                "Any further CRISPResso flags, verbatim."},
            },
            "required": ["fastq_r1", "amplicon"],
        },
        cli=["crispresso", "analyze"],
        flag_map={"fastq_r1": "--fastq-r1", "fastq_r2": "--fastq-r2",
                   "amplicon": "--amplicon", "guide": "--guide",
                   "conversion": "--conversion", "name": "--name",
                   "output_dir": "--output-dir", "extra_args": "--extra-args"},
        bool_flags={"base_editor"},
    ),

    _T(
        "bean_paper_benchmark",
        "Benchmark against ONE SPECIFIC PAPER: Ryu et al., CRISPR-BEAN, "
        "Nat Genet 56:925-937 (2024), base-editing screens of LDLR / LDL-C "
        "variants. DO NOT USE IT FOR ANY OTHER PAPER — for a general "
        "\"reproduce this paper\" request use `paper_benchmark`, which "
        "routes to the correct chain. Runs on the AUTHORS' OWN deposited "
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
                "max_download_gb": {**_S_NUMBER, "default": 200,
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
                "max_download_gb": {**_S_NUMBER, "default": 200,
                    "description": "Transfer ceiling in GB. A guard, not a "
                        "sampling limit: over it the run refuses rather than "
                        "keeping a partial file. Raise it for large runs."},
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
        bool_flags=("verbose",),
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
        bool_flags=("verbose",),
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
        bool_flags=("urls_only",),
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
        bool_flags=("fetch",),
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
        bool_flags=("list",),
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


# Names that came from user manifests, so a refresh can tell "re-authored
# extension" (replace) from "manifest shadowing a built-in" (skip).
_USER_TOOL_NAMES: "set[str]" = set()


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
        if spec["name"] in _BY_NAME and spec["name"] not in _USER_TOOL_NAMES:
            # Recorded for `igvfagent extensions`, logged at INFO. As a
            # WARNING it reached stderr of EVERY tool subprocess (no handler
            # is configured at import time, so Python's last-resort handler
            # printed it), and the agent's failure banner then showed this
            # line instead of the tool's own error.
            msg = (f"{spec.get('source')}: user tool `{spec['name']}` shadows "
                   f"a built-in tool of the same name; skipped (remove the "
                   f"manifest -- the built-in already does this)")
            logger.info(msg)
            try:
                _userext._note(msg)
            except Exception:
                pass
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
        if tool.name in _USER_TOOL_NAMES:
            # Re-authored extension: replace the stale manifest in place so a
            # long-lived process (the Streamlit UI) sees the new parameters.
            # Built-ins are never replaced -- the check above skips them.
            _TOOLS[:] = [t if t.name != tool.name else tool for t in _TOOLS]
        else:
            _TOOLS.append(tool)
        _BY_NAME[tool.name] = tool
        _USER_TOOL_NAMES.add(tool.name)


_merge_user_tools()


def refresh_user_tools() -> int:
    """Re-scan the extension directories and absorb any new user tools.

    For long-lived processes (the Streamlit UI) where the import-time
    merge already happened: call after a manifest is added on disk.
    Built-in names are never redefined; a user tool whose manifest was
    re-authored is replaced in place. Returns the number of tools added.
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


def _flag_value(flag: str, value: str) -> "list[str]":
    """``--flag=value`` when the value starts with '-' (a Markdown list, a
    negative number, "-1"): argparse would read ``--flag -x`` as two options."""
    return [f"{flag}={value}"] if value.startswith("-") and flag.startswith("--") else [flag, value]


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
                argv.extend(_flag_value(flag, _coerce_value(v)))
            continue
        flag = tool.flag_map.get(name, "--" + name.replace("_", "-"))
        # Convention: ``flag_map={name: ""}`` means the argument is
        # positional. Skip emitting an empty-string flag token before
        # the value — that would make argparse choke with exit_code=2.
        if flag == "":
            argv.append(_coerce_value(value))
        else:
            argv.extend(_flag_value(flag, _coerce_value(value)))
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
    # Carry the signed-in user into the subprocess. A thread-local cannot
    # cross a fork, and without this every history_recall / history_search the
    # agent makes would run unattributed -- which the history store reads as
    # "no accounts here, show everything", leaking one user's private sessions
    # into another user's answer.
    try:
        from . import _history as _h

        acting = _h.actor()
        if acting:
            env["IGVF_ACTING_USER"] = acting
        env["IGVF_ACTING_ADMIN"] = os.environ.get("IGVF_ACTING_ADMIN", "0")
    except Exception:
        pass
    if extra_env:
        env.update(extra_env)
    logger.info("tool=%s argv=%s", name, shlex.join(argv))
    try:
        # stdin=DEVNULL: a tool must never inherit ours. Under the MCP server
        # our stdin IS the JSON-RPC stream, and a child that reads it (an
        # agent-authored extension did, when an argument was missing) eats
        # the requests and hangs every later call.
        proc = subprocess.run(
            argv, cwd=cwd, env=env, capture_output=True, stdin=subprocess.DEVNULL,
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
