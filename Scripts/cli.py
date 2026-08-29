"""Top-level CLI dispatcher for the ``igvfagent`` console script.

Exposes every existing skill behind a short, memorable subcommand:

  igvfagent kg gene APOE --depth 2
  igvfagent splitseq retrieve --limit 50
  igvfagent ref learn --topic '10x multiome human putamen'
  igvfagent portal-kg pull --tissue macrophage --limit 100

The dispatcher does not duplicate any argparse logic; it imports the
target module (which still has its own ``main()`` and argparse parser)
and delegates with a synthesized ``sys.argv``. This keeps a single
source of truth for each skill's CLI surface.
"""

from __future__ import annotations

import importlib
import os
import sys
from typing import Optional

from . import __version__

# Skill name -> (module path, one-line description)
SKILLS: "dict[str, tuple[str, str]]" = {
    "client":           ("igvfagent.igvf_client",
                          "IGVF Portal / Catalog / Knowledge Graph / ENCODE client"),
    "data":             ("igvfagent.igvf_data_skills",
                          "Catalog / Portal / ENCODE overview + smoke summaries"),
    "frontpage":        ("igvfagent.igvf_frontpage_summary",
                          "Refresh front-page Portal + Knowledge Graph stats"),
    "specialized":      ("igvfagent.igvf_specialized_data_skills",
                          "Specialized IGVF assay catalog (Parse SPLiT-seq, etc.)"),
    "explain":          ("igvfagent.data_illustration_interpretation",
                          "Explain an IGVF / ENCODE accession or search URL"),
    "variant":          ("igvfagent.annotate_variant_list",
                          "Annotate variants against IGVF Catalog evidence"),
    "advanced-variant": ("igvfagent.advanced_variant_analysis",
                          "Integrated variant scoring + logistic + report"),
    "enhancer":         ("igvfagent.enhancer_gene_linkage_skills",
                          "Enhancer-gene linkage retrieval and comparison"),
    "ccre":             ("igvfagent.ccre_linkage_annotation_skills",
                          "cCRE / FAVOR / linkage annotations"),
    "mpra":             ("igvfagent.mpra_data_skills",
                          "MPRA / STARR / BlueSTARR retrieval and analysis"),
    "starrseq":         ("igvfagent.starr_seq_skill",
                          "STARR-seq allelic test (mpralm clean-room rewrite)"),
    "share":            ("igvfagent.share_seq_skill",
                          "SHARE-seq joint scATAC+scRNA QC (Ma 2020 / Broad pipeline)"),
    "flowfish":         ("igvfagent.flowfish_crispr_skill",
                          "CRISPRi Flow-FISH screen analysis (Fulco/Nasser/Engreitz)"),
    "crispri":          ("igvfagent.crispri_data_skills",
                          "CRISPRi / CRISPR-FACS / Perturb-seq evidence"),
    "singlecell":       ("igvfagent.single_cell_data_skills",
                          "Single-cell discovery and example analysis"),
    "multiome":         ("igvfagent.multiome_10x_pipeline",
                          "10x Multiome retrieval pipeline"),
    "splitseq":         ("igvfagent.splitseq_pipeline",
                          "Parse SPLiT-seq end-to-end pipeline"),
    "grn":              ("igvfagent.catalog_grn_skill",
                          "Gene regulatory network (dEx) + protein-variant "
                          "effects from the Catalog"),
    "kg":               ("igvfagent.kg_traversal_skill",
                          "IGVF Knowledge Graph multi-hop traversal"),
    "portal-kg":        ("igvfagent.portal_to_kg_skill",
                          "Portal → local Knowledge Graph ETL"),
    "sce2g-kg":         ("igvfagent.sce2g_kg_skill",
                          "scE2G element→gene linkages → local KG (bulk, "
                          "adaptive region tiling, resumable, runs to "
                          "completion regardless of size)"),
    "kg-mirror":        ("igvfagent.kg_mirror_skill",
                          "Local IGVF KG mirror (Arango -> Parquet + DuckDB)"),
    "ref":              ("igvfagent.reference_skill",
                          "Literature retrieval / validation / study design"),
    "bench":            ("igvfagent.benchmark_skill",
                          "Paper → reproduction benchmark. Resolve a publication "
                          "from a title / URL / DOI / PMID / author+journal+year, "
                          "harvest its Data Availability statement + accessions "
                          "from the full text, route them onto an IGVFagent "
                          "analysis chain, and scaffold a runnable "
                          "Benchmarks/<paper-id>/ that concordance.py scores. "
                          "Subcommands: resolve, harvest, route, scaffold, run, "
                          "score, report, pipeline, selftest, list-routes."),
    "mavedb":           ("igvfagent.mavedb_mapping_skill",
                          "MaveDB scoreset → genomic coords (chr/pos/ref/alt)"),
    "encode":           ("igvfagent.encode_pipeline",
                          "ENCODE ChIP/ATAC/DNase/Hi-C/ChIA-PET pipeline"),
    "se-targets":       ("igvfagent.se_target_pipeline",
                          "Super-enhancer → target-gene pipeline (ENCODE)"),
    "geo":              ("igvfagent.geo_retrieval_skill",
                          "NCBI GEO retrieval (search / metadata / download)"),
    "rnaseq":           ("igvfagent.rnaseq_analysis_skill",
                          "Bulk RNA-seq QC / PCA / DEG / DEG→cCRE linkage"),
    "proteomics":       ("igvfagent.proteomics_skill",
                          "Proteomics & PPI: BioGRID/IntAct/HuRI/Reactome/KEGG/"
                          "IGVF integration, KG, viz, literature survey"),
    "sc-analyze":       ("igvfagent.singlecell_analysis",
                          "Single-cell analysis: QC, PCA, UMAP, t-SNE, Leiden, "
                          "markers, publication figures (Scanpy-driven)"),
    "perturb-catalog":  ("igvfagent.perturbation_catalog_skill",
                          "Perturbation Catalogue retrieval: MAVE / CRISPR-screen "
                          "/ Perturb-seq datasets, per-row effects, GSEA, downloads"),
    "multiseq":         ("igvfagent.multiseq_analysis_skill",
                          "MULTI-seq / Cell Hashing demultiplexing "
                          "(Python port of deMULTIplex2)"),
    "skillcard":        ("igvfagent.skillcard_skill",
                          "Machine-readable skill specifications (skill cards) "
                          "with validation status derived from the benchmarks"),
    "pathway-viz":      ("igvfagent.pathway_viz_skill",
                          "Pathway / PPI network figure for a gene list "
                          "(STRING + Reactome + local KG)"),
    "mcp":              ("igvfagent.mcp_server_skill",
                          "Serve the skill registry over the Model Context "
                          "Protocol so other agents can call IGVF skills"),
    "playbook-freeze":  ("igvfagent.playbook_freeze_skill",
                          "Freeze a recorded session into a deterministic "
                          "playbook with pinned args and artefact hashes"),
    "document":         ("igvfagent.document_ingest_skill",
                          "Read an uploaded manuscript (PDF/DOCX/text) and "
                          "derive accessions, assays and a reproduction plan"),
    "eqtl-mpra":        ("igvfagent.eqtl_mpra_concordance_skill",
                          "Tissue-filtered fine-mapped eQTLs and MPRA-vs-eQTL "
                          "correlation / sign-concordance"),
    "eval-tiers":       ("igvfagent.eval_tiers_skill",
                          "Tier 2 (planning / tool selection) and Tier 3 "
                          "(conclusion validity) evaluation"),
    "variant-verify":   ("igvfagent.variant_verify_skill",
                          "Cross-source verification of variant annotations "
                          "(FAVOR snapshot vs live ClinVar)"),
    "variant-list":     ("igvfagent.variant_list_skill",
                          "Annotate a pasted/file variant list in any notation "
                          "(chr-pos-ref-alt, rsID, SPDI, HGVS, VCF) via FAVOR + "
                          "IGVF Catalog, and add the results to the local KG"),
    "artifact":         ("igvfagent.artifact_read_skill",
                          "Read back reports/manifests the agent produced "
                          "(workspace-contained: read / grep / ls)"),
    "extauthor":        ("igvfagent.ext_author_skill",
                          "Author new IGVFagent tools/skills from the agent "
                          "itself (opt-in: IGVF_ALLOW_AGENT_AUTHORING=1)"),
    "warehouse":        ("igvfagent.warehouse_skill",
                          "Central DuckDB warehouse — Silver tier for the "
                          "IGVF integrated data layer"),
    "network":          ("igvfagent.network_integration_skill",
                          "Network integration — clean-room MILP for "
                          "context-specific subnetworks (CARNIVAL / Steiner)"),
    "enrich":           ("igvfagent.enrichment_skill",
                          "GO + Pathway enrichment validation (ORA via Enrichr; "
                          "GSEA preranked) over GO_BP/MF/CC + Reactome + KEGG "
                          "+ WikiPathways + MSigDB Hallmark"),
    "portal":           ("igvfagent.portal_query_skill",
                          "IGVF Portal canonical-query layer — faceted search, "
                          "report.tsv, batch-download, endpoint-param introspection "
                          "(clean-room reimpl of IGVF-DACC/igvf-portal-mcp)"),
    "catalog":          ("igvfagent.catalog_query_skill",
                          "IGVF Catalog (KG) canonical-query layer — universal "
                          "get-entity, search-region, find-associations (semantic), "
                          "find-ld, resolve-id, list-sources "
                          "(clean-room reimpl of IGVF-DACC/igvf-catalog-mcp)"),
    "chipatlas":        ("igvfagent.chipatlas_skill",
                          "ChIP-Atlas (Ohta/Oki) reprocessed ChIP/ATAC/DNase/Bisulfite "
                          "peak archive — browse + per-experiment files + assembled "
                          "BEDs + Target-Genes + WABI Enrichment "
                          "(clean-room reimpl of inutano/chip-atlas)"),
    "synapse":          ("igvfagent.synapse_skill",
                          "Sage Bionetworks Synapse retrieval — entity / "
                          "children / walk / search / download. "
                          "Anonymous for public deposits; PAT-authenticated "
                          "for PsychENCODE, AMP-AD, AMP-PD, and IGVF-controlled "
                          "Synapse cohorts (set SYNAPSE_AUTH_TOKEN). "
                          "Pure urllib + json, no synapseclient dep."),
    "sceps":            ("igvfagent.sceps_skill",
                          "scEPS single-cell disease-neighborhood statistics — "
                          "clean-room reimpl of Genentech/sceps. Integrates GWAS "
                          "(MAGMA Z) + single-cell atlas; per-neighborhood "
                          "variance-component d-statistic (GWAS vs matched "
                          "control genes). Subcommands: estimate. Runs under the "
                          "scEPS env (scanpy/anndata/statsmodels)."),
    "calibrate":        ("igvfagent.excalibr_skill",
                          "Functional-assay calibration to ACMG/AMP evidence — "
                          "clean-room reimpl of rosstewart/exCALIBR (Zeiberg "
                          "et al. bioRxiv 2025). Bootstrap constrained "
                          "skew-normal mixture EM + Bayesian (Tavtigian) "
                          "calibration turns MAVE / VAMP-seq / SGE scores into "
                          "PS3 / BS3 evidence strengths. Subcommands: "
                          "thresholds, prepare, run, assign, selftest."),
    "open4gene":        ("igvfagent.open4gene_skill",
                          "Open4Gene peak-to-gene linkage — clean-room reimpl of "
                          "hbliu/Open4Gene (Liu et al. Science 2025). Hurdle model "
                          "(logit zero + zero-truncated NB count) linking snATAC "
                          "peaks to snRNA genes across multiome cells, per cell "
                          "type, with covariates. Subcommands: link. statsmodels-"
                          "based; validated vs the R pscl::hurdle reference."),
    "figshare":         ("igvfagent.figshare_skill",
                          "figshare retrieval — article / files / download / "
                          "search. Resolves numeric id, DOI, article URL, or "
                          "private /s/<token> share link (as printed in paper "
                          "Data Availability statements). Downloads with md5 "
                          "verify. Pure urllib + json; the general-purpose "
                          "counterpart to Zenodo for research-data deposits."),
}

# All reserved namespaces are now wired. Kept as a (currently empty)
# anchor so future placeholders have a documented home.
RESERVED: "dict[str, str]" = {}

# Introspection commands wired in step 2: lists the LLM provider router's
# registered backends, the tool registry the agent runtime will see, and
# (step 5) the locally installed Ollama models when the daemon is running.
# `extensions` inspects user-supplied skills/tools (see _userext.py).
INTROSPECTION = ("backends", "tools", "models", "extensions")

# Top-level commands wired in step 3 (`ask`) and step 4 (`ui`).
TOP_LEVEL = ("ask", "ui", "playbook", "eval", "localstore", "consistency")


def _user_skills() -> "dict[str, dict]":
    """User-supplied skills discovered from the extension directories
    (see ``_userext``), minus any name shadowed by a built-in command.

    Never raises: extension discovery is best-effort so a broken user
    file can't take the whole CLI down.
    """
    try:
        from . import _userext
        found = _userext.discover_skills()
    except Exception:
        return {}
    reserved = (set(SKILLS) | set(RESERVED) | set(INTROSPECTION)
                | set(TOP_LEVEL))
    return {name: entry for name, entry in found.items()
            if name not in reserved}


def _print_counts(*, json_out: bool = False) -> int:
    """Authoritative inventory: skills, tools, benchmarks, checks.

    Derived at call time from the CLI registry, the tool registry and
    Benchmarks/taxonomy.py — the same sources the software itself uses — so
    the paper and the repository cannot disagree.
    """
    import json as _json
    import subprocess as _sp
    from pathlib import Path as _P

    root = _P(os.environ.get("IGVF_PROJECT_ROOT")
              or _P(__file__).resolve().parents[1]).resolve()

    counts = {"version": __version__,
              "skills": len(SKILLS),
              "top_level_commands": len(TOP_LEVEL)}
    try:
        from igvfagent import _tools, _llm
    except Exception:
        import _tools, _llm  # type: ignore
    tools = _tools.list_tools()
    dicts = [{"name": t.name, "description": t.description or ""} for t in tools]
    exposed = _llm.canonical_tools(dicts) or []
    counts["tools_registered"] = len(tools)
    counts["tools_exposed_default"] = len(exposed)
    counts["tool_cap"] = int(os.environ.get("IGVF_LLM_MAX_TOOLS",
                                             _llm._DEFAULT_MAX_TOOLS))
    counts["tools_starred"] = sum(
        1 for t in tools if (t.description or "").lstrip().startswith("★"))

    tax_path = root / "Benchmarks" / "taxonomy.py"
    if tax_path.is_file():
        try:
            out = _sp.run([sys.executable, str(tax_path), "--json"],
                          capture_output=True, text=True, timeout=120).stdout
            tax = _json.loads(out)
            counts["benchmarks"] = tax["n_benchmarks"]
            counts["benchmark_checks"] = sum(b["n_checks"] for b in tax["benchmarks"])
            counts["benchmarks_by_class"] = tax["summary"]
        except Exception:
            pass

    try:
        counts["git_commit"] = _sp.run(
            ["git", "-C", str(root), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=20).stdout.strip() or None
    except Exception:
        counts["git_commit"] = None

    if json_out:
        print(_json.dumps(counts, indent=2))
        return 0
    print(f"igvfagent {counts['version']}"
          + (f"  (commit {counts['git_commit']})" if counts.get("git_commit") else ""))
    print(f"  skills                 {counts['skills']}"
          f"   (+{counts['top_level_commands']} top-level commands)")
    print(f"  tools registered       {counts['tools_registered']}")
    print(f"  tools exposed default  {counts['tools_exposed_default']}"
          f"   (cap IGVF_LLM_MAX_TOOLS={counts['tool_cap']}, "
          f"{counts['tools_starred']} starred kept first)")
    if "benchmarks" in counts:
        c = counts["benchmarks_by_class"]
        print(f"  benchmarks             {counts['benchmarks']}"
              f"   ({counts['benchmark_checks']} machine-checked criteria)")
        print(f"    quantitative reproduction {c.get('A', 0)}"
              f"  ·  method validation {c.get('B', 0)}"
              f"  ·  retrieval {c.get('C', 0)}"
              f"  ·  artefact {c.get('D', 0)}")
    return 0


def _print_help() -> None:
    print(f"igvfagent {__version__}")
    print()
    print("Local, auditable AI agent for IGVF / ENCODE data analysis with")
    print("a built-in Plan → Action → Results → Evaluation orchestrator.")
    print()
    print("Usage:")
    print("  igvfagent <skill> [subcommand] [args ...]")
    print("  igvfagent --version | --help | --list")
    print()
    user = _user_skills()
    print("Available skills:")
    width = max(len(name) for name in [*SKILLS, *user])
    for name, (_, doc) in SKILLS.items():
        print(f"  {name:{width}}  {doc}")
    if user:
        print()
        print("User skills (from ~/.igvfagent + UserExtensions/):")
        for name, entry in user.items():
            print(f"  {name:{width}}  {entry['description']}")
    print()
    print("Top-level commands:")
    print(f"  {'ask':{width}}  Natural-language ReAct agent (LLM-driven)")
    print(f"  {'ui':{width}}  Launch the Streamlit browser UI")
    print(f"  {'playbook':{width}}  Run a deterministic YAML playbook (Docs/Playbooks/)")
    print(f"  {'eval':{width}}  Backend-comparison eval harness")
    print()
    print("Introspection:")
    print(f"  {'backends':{width}}  List configured LLM provider backends")
    print(f"  {'tools':{width}}  List the tool registry the agent runtime sees")
    print(f"  {'models':{width}}  List Ollama models installed on the local daemon")
    print(f"  {'extensions':{width}}  List user-supplied skills / tools and where "
          "they were found")
    if RESERVED:
        print()
        print("Reserved (not yet wired):")
        for name, doc in RESERVED.items():
            print(f"  {name:{width}}  {doc}")
    print()
    print("Use `igvfagent <skill> --help` for a skill's subcommands.")
    print("Documentation:  https://github.com/zhouhufeng/IGVFagent")


def _print_list() -> None:
    """Machine-readable skill list (one name per line)."""
    for name in SKILLS:
        print(name)
    for name in _user_skills():
        print(name)


def main(argv: Optional["list[str]"] = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("--help", "-h", "help"):
        _print_help()
        return 0
    if args[0] in ("--version", "-V", "version"):
        # `--counts` prints authoritative, machine-derived inventory numbers so
        # a manuscript can cite one frozen artefact instead of hand-copied
        # figures that drift. Every number here is read from the running
        # system, never typed in.
        if "--counts" in args[1:]:
            return _print_counts(json_out="--json" in args[1:])
        print(__version__)
        return 0
    if args[0] in ("--list", "list"):
        _print_list()
        return 0

    skill = args[0]
    if skill in RESERVED:
        sys.stderr.write(
            f"`{skill}` is reserved but not yet wired. {RESERVED[skill]}\n"
        )
        return 64
    if skill in INTROSPECTION:
        return _run_introspection(skill, args[1:])
    if skill in TOP_LEVEL:
        return _run_top_level(skill, args[1:])
    if skill in SKILLS:
        module_name, _ = SKILLS[skill]
        try:
            mod = importlib.import_module(module_name)
        except ImportError as exc:
            sys.stderr.write(f"failed to load skill `{skill}`: {exc}\n")
            return 1
    else:
        # Not a built-in — fall back to user-supplied skills discovered
        # under ~/.igvfagent/skills/ and UserExtensions/skills/.
        entry = _user_skills().get(skill)
        if entry is None:
            sys.stderr.write(f"unknown skill: {skill}\n")
            sys.stderr.write("Run `igvfagent --help` for the skill list.\n")
            return 2
        try:
            from . import _userext
            mod = _userext.load_skill(skill, entry)
        except Exception as exc:
            sys.stderr.write(f"failed to load user skill `{skill}` "
                             f"({entry['path']}): {exc}\n")
            return 1
    if not hasattr(mod, "main"):
        sys.stderr.write(f"skill `{skill}` has no main() entrypoint\n")
        return 1

    # Hand off to the underlying skill's argparse parser. We synthesize
    # argv[0] so the skill's own --help text shows the friendly skill
    # name rather than the dispatcher path.
    sys.argv = [f"igvfagent {skill}"] + args[1:]
    rc_code = 0
    try:
        rc = mod.main()
        rc_code = int(rc) if isinstance(rc, int) else 0
    except SystemExit as exc:  # argparse + sys.exit propagation
        code = exc.code
        if code is None:
            rc_code = 0
        elif isinstance(code, int):
            rc_code = code
        else:
            # Non-int code (e.g. sys.exit("some error message")) — print to
            # stderr and use exit 1, matching Python's own default behavior.
            sys.stderr.write(f"{code}\n")
            rc_code = 1
    finally:
        # Core default: every skill run grows the local KG/DB from its new
        # on-disk outputs. Best-effort; never affects the command's exit code.
        _post_run_harvest(skill)
    return rc_code


def _post_run_harvest(skill: str) -> None:
    if os.environ.get("IGVF_LOCALSTORE", "1") == "0":
        return
    try:
        try:
            from igvfagent import _localstore as ls
        except Exception:
            import _localstore as ls  # type: ignore
        ls.harvest()
    except Exception:
        pass  # the skill already did its job; growth is a bonus


def _run_top_level(skill: str, args: "list[str]") -> int:
    if skill == "ask":
        return _run_ask(args)
    if skill == "ui":
        return _run_ui(args)
    if skill == "playbook":
        return _run_playbook(args)
    if skill == "eval":
        return _run_eval(args)
    if skill == "localstore":
        return _run_localstore(args)
    if skill == "consistency":
        try:
            from igvfagent import consistency_check as cc
        except Exception:
            import consistency_check as cc  # type: ignore
        return cc.main(args)
    return 2


def _run_localstore(args: "list[str]") -> int:
    """Inspect / grow the local Knowledge Graph + database.

    Subcommands:
      stats     show current KG/DB size (nodes, edges, downloads, analyses)
      harvest   scan Docs/ + Benchmarks/_data and ingest anything new
      backfill  re-ingest table CONTENT from runs the ledger already saw
                (harvest skips those forever, so artefacts produced before
                table ingestion existed need this once)
    """
    import json
    try:
        from igvfagent import _localstore as ls
    except Exception:
        sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
        import _localstore as ls  # type: ignore
    sub = args[0] if args else "stats"
    if sub == "harvest":
        print("Harvesting on-disk downloads + analyses into the local KG/DB …")
        print(json.dumps(ls.harvest(), indent=2, default=str))
        print(json.dumps(ls.stats(), indent=2, default=str))
        return 0
    if sub == "backfill":
        print("Back-filling entity edges from existing artefact tables …")
        print(json.dumps(ls.backfill(), indent=2, default=str))
        print(json.dumps(ls.stats(), indent=2, default=str))
        return 0
    if sub in ("stats", "", "--help", "-h"):
        print(json.dumps(ls.stats(), indent=2, default=str))
        return 0
    sys.stderr.write(f"unknown localstore subcommand: {sub}\n"
                     "  use: igvfagent localstore [stats|harvest]\n")
    return 2


def _run_ask(args: "list[str]") -> int:
    import argparse
    parser = argparse.ArgumentParser(
        prog="igvfagent ask",
        description="Natural-language ReAct agent. The LLM plans, calls "
                    "IGVFagent skills as tools, and writes a final answer "
                    "with file paths to the artefacts.",
    )
    parser.add_argument("query", nargs="+",
                         help="Your question (any natural-language text).")
    parser.add_argument("--backend", default=None,
                         help="LLM backend (anthropic / openai / codex / "
                              "ollama / vllm / tgi / groq / together / "
                              "deepinfra / huggingface / custom). Default: "
                              "auto-detect.")
    parser.add_argument("--model", default=None,
                         help="Model name. Default: backend-specific (e.g. "
                              "claude-sonnet-4-5, gpt-4o-mini, qwen3:8b).")
    parser.add_argument("--max-iterations", type=int, default=8,
                         help="Cap on Plan→Action loop iterations.")
    parser.add_argument("--max-tokens", type=int, default=4096)
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--tool", action="append", default=None,
                         help="Restrict the agent to a subset of tools. "
                              "Repeat the flag for multiple tools.")
    parser.add_argument("--system-prompt-file", default=None,
                         help="Path to a custom system prompt.")
    parser.add_argument("--quiet", action="store_true",
                         help="Suppress the per-step progress trace.")
    parser.add_argument("--no-persist", action="store_true",
                         help="Don't write the transcript / report files.")
    parsed = parser.parse_args(args)

    # Lazy-import the agent module so `igvfagent --help` doesn't pay the
    # SDK-import cost.
    try:
        from igvfagent import _agent
    except ImportError:
        import _agent  # type: ignore[no-redef]

    system_prompt = None
    if parsed.system_prompt_file:
        from pathlib import Path as _P
        system_prompt = _P(parsed.system_prompt_file).read_text()

    cb = None if parsed.quiet else _agent._print_callback
    try:
        result = _agent.run(
            " ".join(parsed.query),
            backend=parsed.backend, model=parsed.model,
            max_iterations=parsed.max_iterations,
            max_tokens=parsed.max_tokens,
            temperature=parsed.temperature,
            tools_subset=parsed.tool,
            system_prompt=system_prompt,
            callback=cb,
            persist=not parsed.no_persist,
        )
    except RuntimeError as e:
        sys.stderr.write(f"agent failed to start: {e}\n")
        sys.stderr.write(
            "Hint: install an LLM SDK and set credentials. For local "
            "free use:\n"
            "  pip install 'igvfagent[llm]'\n"
            "  ollama serve   # in another terminal\n"
            "  ollama pull qwen3:8b\n"
            "Then re-run `igvfagent ask ...`.\n"
        )
        return 1

    return 0 if result.stop_reason == "complete" else 2


def _run_playbook(args: "list[str]") -> int:
    """Top-level: run a deterministic YAML playbook of tool calls."""
    import argparse, json, sys, time
    parser = argparse.ArgumentParser(
        prog="igvfagent playbook",
        description="Execute a pinned YAML playbook of tool calls "
                    "deterministically. Same artefacts across LLM backends; "
                    "only the final synthesis prose may vary.",
    )
    parser.add_argument("playbook",
                         help="Path to YAML playbook (or short name under "
                              "Docs/Playbooks/, e.g. `apoe_evidence_pack`).")
    parser.add_argument("--param", action="append", default=[],
                         help="Override a playbook parameter: key=value "
                              "(repeat for multiple).")
    parser.add_argument("--backend", default=None)
    parser.add_argument("--model", default=None)
    parser.add_argument("--max-tokens", type=int, default=4096)
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--no-synthesize", action="store_true",
                         help="Skip the final LLM synthesis step.")
    parser.add_argument("--no-persist", action="store_true",
                         help="Don't write the run summary to Docs/.")
    parser.add_argument("--quiet", action="store_true")
    parsed = parser.parse_args(args)

    try:
        from igvfagent import _playbook
    except ImportError:
        import _playbook  # type: ignore[no-redef]
    from pathlib import Path as _P

    # Resolve playbook path: explicit file path OR short name under Docs/Playbooks/
    pb_arg = parsed.playbook
    pb_path = _P(pb_arg)
    if not pb_path.is_file():
        for cand in [
            _P("Docs/Playbooks") / pb_arg,
            _P("Docs/Playbooks") / f"{pb_arg}.yaml",
            _P("Docs/Playbooks") / f"{pb_arg}.yml",
        ]:
            if cand.is_file():
                pb_path = cand
                break
    if not pb_path.is_file():
        print(f"Playbook not found: {pb_arg}", file=sys.stderr)
        print("Available playbooks under Docs/Playbooks/:")
        d = _P("Docs/Playbooks")
        if d.is_dir():
            for f in sorted(d.glob("*.y*ml")):
                print(f"  {f.stem}")
        return 2

    params = _playbook.parse_param_args(parsed.param)

    def _on_step(rec):
        if parsed.quiet:
            return
        flag = "✅" if rec["exit_code"] == 0 else "❌"
        print(f"  [{rec['idx']+1}] {flag} {rec['tool']}  exit={rec['exit_code']}  "
              f"artefacts={rec['n_artefacts']}  {rec['elapsed_sec']:.1f}s")

    if not parsed.quiet:
        print(f"▶ Running playbook: {pb_path}")
        if params:
            print(f"  params: {params}")
    result = _playbook.run_playbook(
        pb_path, params=params,
        backend=parsed.backend, model=parsed.model,
        synthesize=not parsed.no_synthesize,
        max_tokens=parsed.max_tokens, temperature=parsed.temperature,
        on_step=_on_step,
    )
    if not parsed.quiet:
        print(f"\n✓ {len(result['steps'])} steps, "
              f"{len(result['artefacts'])} artefacts, "
              f"{result['elapsed_sec']:.1f}s")

    # Persist a run summary alongside the agent transcripts
    if not parsed.no_persist:
        from datetime import datetime
        out_dir = _P("Docs/Playbooks/Runs")
        out_dir.mkdir(parents=True, exist_ok=True)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        stem = pb_path.stem
        json_path = out_dir / f"{ts}_{stem}_run.json"
        md_path = out_dir / f"{ts}_{stem}_run.md"
        json_path.write_text(json.dumps(result, indent=2, default=str))
        lines = [f"# Playbook run: {result.get('study')}", "",
                 f"Playbook: `{pb_path}`",
                 f"Params: `{json.dumps(result.get('params') or {})}`",
                 f"Elapsed: {result['elapsed_sec']:.1f}s", "",
                 "## Steps", "",
                 "| # | tool | exit | artefacts | sec |",
                 "|---|---|---:|---:|---:|"]
        for s in result["steps"]:
            lines.append(f"| {s['idx']+1} | `{s['tool']}` | {s['exit_code']} | "
                          f"{s['n_artefacts']} | {s['elapsed_sec']} |")
        lines += ["", "## Artefacts", ""]
        for a in result["artefacts"]:
            lines.append(f"- `{a}`")
        if result.get("final_answer"):
            lines += ["", "## Synthesis", "", result["final_answer"]]
        md_path.write_text("\n".join(lines))
        if not parsed.quiet:
            print(f"  summary: {md_path}")
            print(f"  json:    {json_path}")
    return 0


def _run_eval(args: "list[str]") -> int:
    """Top-level: backend-comparison eval harness."""
    import argparse, json, sys, time
    parser = argparse.ArgumentParser(
        prog="igvfagent eval",
        description="Run the same query through multiple (backend, model) "
                    "pairs and report pairwise divergence. Useful for "
                    "deciding which questions are any-backend-safe.",
    )
    parser.add_argument("query", nargs="+",
                         help="The question to compare across backends.")
    parser.add_argument("--backends", required=True,
                         help="Comma-separated backends (e.g. anthropic,ollama).")
    parser.add_argument("--models", default=None,
                         help="Comma-separated models, one per backend. "
                              "If omitted, each backend uses its default.")
    parser.add_argument("--runs", type=int, default=1,
                         help="Replicates per (backend, model) pair.")
    parser.add_argument("--max-iterations", type=int, default=8)
    parser.add_argument("--max-tokens", type=int, default=4096)
    parser.add_argument("--label", default=None)
    parser.add_argument("--quiet", action="store_true")
    parsed = parser.parse_args(args)

    try:
        from igvfagent import _eval
    except ImportError:
        import _eval  # type: ignore[no-redef]
    from pathlib import Path as _P
    from datetime import datetime

    backends = [b.strip() for b in parsed.backends.split(",") if b.strip()]
    models_arg = parsed.models or ""
    models = [m.strip() or None for m in models_arg.split(",")] if models_arg else [None] * len(backends)
    if len(models) == 1 and len(backends) > 1:
        models = models * len(backends)
    if len(backends) != len(models):
        print(f"--backends ({len(backends)}) and --models ({len(models)}) must align",
              file=sys.stderr)
        return 2
    pairs = list(zip(backends, models))
    question = " ".join(parsed.query)

    print(f"▶ Eval: '{question}'")
    print(f"  pairs: {pairs}")
    print(f"  runs/pair: {parsed.runs}")

    records = _eval.run_eval(
        question, pairs=pairs, runs=parsed.runs,
        max_iterations=parsed.max_iterations,
        max_tokens=parsed.max_tokens,
        quiet=parsed.quiet,
    )
    diffs = _eval.diff_matrix(records)

    out_dir = _P("Docs/Eval")
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    label = parsed.label or "eval"
    md_path = out_dir / f"{ts}_{label}_eval.md"
    json_path = out_dir / f"{ts}_{label}_eval.json"
    _eval.write_eval_report(question, records, diffs, md_path)
    json_path.write_text(json.dumps(
        {"question": question, "records": records, "diffs": diffs}, indent=2,
        default=str))

    print(f"\n✓ Report: {md_path}")
    print(f"  JSON:   {json_path}")
    print()
    print(f"  {len(records)} runs, {len(diffs)} pairwise comparisons")
    if diffs:
        avg_jac = sum(d["jaccard_artefacts"] for d in diffs) / len(diffs)
        avg_cos = sum(d["cos_final_answer"] for d in diffs) / len(diffs)
        print(f"  avg Jaccard(artefacts) = {avg_jac:.3f}")
        print(f"  avg cos(final_answer)   = {avg_cos:.3f}")
    return 0




def _run_ui(args: "list[str]") -> int:
    import argparse, shutil, subprocess
    parser = argparse.ArgumentParser(
        prog="igvfagent ui",
        description="Launch the Streamlit browser UI for IGVFagent.",
    )
    parser.add_argument("--port", type=int, default=8501)
    # Default to 127.0.0.1 (explicit IPv4) instead of "localhost" because
    # most modern browsers resolve "localhost" to IPv6 (::1) first while
    # Streamlit binds IPv4 only — leading to "site can't be reached" on
    # http://localhost:8501. The IPv4 form is unambiguous everywhere.
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--no-browser", action="store_true",
                         help="Don't auto-open the browser tab.")
    parser.add_argument("--", dest="extra", nargs=argparse.REMAINDER,
                         help="Forward arbitrary `streamlit run` flags.")
    parsed = parser.parse_args(args)

    # Prefer ``python -m streamlit run`` so we always invoke the streamlit
    # installed in the SAME venv as igvfagent, regardless of whether the
    # venv's bin/ is on PATH.
    try:
        import streamlit  # noqa: F401  -- presence check only
    except ImportError:
        sys.stderr.write(
            "streamlit is not installed in this environment.\n"
            "Install with: pip install 'igvfagent[ui]'\n"
            "Then re-run: igvfagent ui\n"
        )
        return 1

    # Prefer the streamlit_app.py shipped inside the installed package;
    # fall back to the file next to this dispatcher when running out of
    # a checkout without `pip install`. We resolve the path *without*
    # importing the module so Streamlit's import-time hooks
    # (set_page_config et al.) don't warn outside a real run context.
    from importlib.util import find_spec
    from pathlib import Path as _P
    spec = find_spec("igvfagent.streamlit_app")
    if spec and spec.origin:
        app_path = spec.origin
    else:
        app_path = str(_P(__file__).resolve().parent / "streamlit_app.py")

    # On first launch Streamlit prints an interactive "Email:" telemetry
    # consent prompt that blocks startup until stdin is touched — which
    # silently hangs `igvfagent ui` for first-time users. Pre-creating
    # an empty credentials file declines the consent and skips the
    # prompt entirely. Idempotent on re-runs.
    cred_dir = _P.home() / ".streamlit"
    cred_file = cred_dir / "credentials.toml"
    if not cred_file.exists():
        try:
            cred_dir.mkdir(parents=True, exist_ok=True)
            cred_file.write_text('[general]\nemail = ""\n')
        except Exception:
            pass  # not fatal; worst case the user sees the prompt

    cmd = [
        sys.executable, "-m", "streamlit", "run", app_path,
        "--server.port", str(parsed.port),
        "--server.address", parsed.host,
        "--browser.gatherUsageStats", "false",
        "--server.headless", "true" if parsed.no_browser else "false",
    ]
    if parsed.no_browser:
        cmd.extend(["--server.headless", "true"])
    if parsed.extra:
        cmd.extend(parsed.extra)

    sys.stderr.write(
        f"\n  → Streamlit UI starting at http://{parsed.host}:{parsed.port}\n"
        f"  → Source: {app_path}\n"
        f"  → Backends: `igvfagent backends`  ·  Tools: `igvfagent tools`\n\n"
    )
    return subprocess.call(cmd)


def _run_introspection(skill: str, args: "list[str]") -> int:
    if skill == "backends":
        from . import _llm
        for name in _llm.list_backends():
            d = _llm.describe_backend(name)
            tag = d.get("sdk", "")
            extras = ""
            if "base_url" in d:
                extras += f"  base_url={d['base_url']}"
            if "api_key_env" in d:
                extras += f"  api_key_env={d['api_key_env']}"
            print(f"{name:14} {tag:10}{extras}")
        return 0
    if skill == "tools":
        from . import _tools
        if "--json" in args:
            import json as _json
            print(_json.dumps([t.to_dict() for t in _tools.list_tools()],
                                indent=2))
            return 0
        for t in _tools.list_tools():
            print(_tools.render_tool_summary(t))
            print()
        return 0
    if skill == "extensions":
        from . import _userext
        user_tools = _userext.discover_tools()
        user_skills = _user_skills()
        problems = _userext.problems()
        if "--json" in args:
            import json as _json
            print(_json.dumps({
                "directories": [str(d) for d in _userext.extension_dirs()],
                "tools":       user_tools,
                "skills":      user_skills,
                "problems":    problems,
            }, indent=2))
            return 0
        print("Extension directories (searched in order):")
        for d in _userext.extension_dirs():
            mark = "found" if d.is_dir() else "absent"
            print(f"  {d}  [{mark}]")
        print()
        print(f"User tools ({len(user_tools)}):")
        for t in user_tools:
            runner = " ".join(t["command"]) if t["command"] \
                else "igvfagent " + " ".join(t["cli"])
            print(f"  {t['name']:28} {runner}")
            print(f"  {'':28} from {t['source']}")
        print()
        print(f"User skills ({len(user_skills)}):")
        for name, entry in user_skills.items():
            print(f"  {name:28} {entry['description']}")
            print(f"  {'':28} from {entry['path']}")
        if problems:
            print()
            print("Problems (skipped definitions):")
            for p in problems:
                print(f"  ! {p}")
        return 0
    if skill == "models":
        from . import _llm
        if "--json" in args:
            import json as _json
            print(_json.dumps(_llm.list_ollama_models(), indent=2))
            return 0
        models = _llm.list_ollama_models()
        if not models:
            sys.stderr.write(
                "No Ollama models found.\n"
                "  - Is `ollama serve` running?\n"
                "  - Or set OLLAMA_HOST_BASE to your daemon "
                "(e.g. http://localhost:11434/v1).\n"
            )
            return 1
        print(f"{'name':52} {'size':>8}  family")
        print("-" * 78)
        for m in sorted(models, key=lambda r: (r.get("name") or "")):
            size = (f"{m['size_gb']:>5.2f}GB"
                    if m.get("size_gb") is not None else "      ?")
            print(f"{(m.get('name') or ''):52} {size:>8}  "
                  f"{m.get('family','') or ''}")
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
