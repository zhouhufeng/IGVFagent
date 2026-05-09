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
        "phenotypes, biosamples, predictions. Accepts rsID/SPDI/HGVS.",
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
        "Discover IGVF Parse SPLiT-seq AnalysisSets by lab/tissue/taxa. "
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
        "Pull MPRA/STARR/BlueSTARR metadata from Catalog or Portal.",
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
        "crispri_pull",
        "Pull CRISPRi/CRISPR-FACS/Perturb-seq evidence from the Catalog.",
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
        "Recommend a study workflow + cognate published studies + matching "
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
        "optional SCREEN cCRE overlay.",
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
            },
            "required": ["region"],
        },
        cli=["encode", "browser"],
        flag_map={"region": "--region", "tracks": "--track",
                   "ccre_bed": "--ccre-bed", "width": "--width",
                   "label": "--label"},
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
        "End-to-end super-enhancer → target-gene pipeline. For a chosen "
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
        "End-to-end bulk RNA-seq: QC + PCA + differential expression "
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

]


_BY_NAME = {t.name: t for t in _TOOLS}


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
    argv tokens."""
    argv = ["igvfagent", *tool.cli]
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
        argv.extend([flag, _coerce_value(value)])
    return argv


# --------------------------- Execution --------------------------------------

_REPORT_RE   = re.compile(r"^(?:Report|Lit manifest|IGVF manifest|Manifest|"
                           r"Evidence pack|Local KG|Plot|Wrote|Wrote report|"
                           r"Wrote: |Playbook):\s*(.+?)\s*$", re.M)


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
    argv = _build_argv(tool, arguments or {})
    # Replace the leading "igvfagent" placeholder with whatever runner we
    # have available (binary or `python -m igvfagent.cli`).
    head = _resolve_igvfagent()
    argv = head + argv[1:]

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
    for m in _REPORT_RE.finditer(stdout or ""):
        line = m.group(0)
        key, _, value = line.partition(":")
        key = key.strip().lower().replace(" ", "_")
        artefacts.setdefault(key, []).append(value.strip())
    return artefacts


# --------------------------- Pretty render ---------------------------------

def render_tool_summary(tool: Tool) -> str:
    return (f"{tool.name}\n  cli: igvfagent {' '.join(tool.cli)}\n"
            f"  desc: {tool.description}\n"
            f"  params: {json.dumps(tool.parameters.get('properties', {}), default=str)[:200]}")


__all__ = [
    "Tool", "list_tools", "get_tool",
    "to_anthropic_schema", "to_openai_schema",
    "execute",
]
