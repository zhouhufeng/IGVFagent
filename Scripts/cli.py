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
    "crispri":          ("igvfagent.crispri_data_skills",
                          "CRISPRi / CRISPR-FACS / Perturb-seq evidence"),
    "singlecell":       ("igvfagent.single_cell_data_skills",
                          "Single-cell discovery and example analysis"),
    "multiome":         ("igvfagent.multiome_10x_pipeline",
                          "10x Multiome retrieval pipeline"),
    "splitseq":         ("igvfagent.splitseq_pipeline",
                          "Parse SPLiT-seq end-to-end pipeline"),
    "kg":               ("igvfagent.kg_traversal_skill",
                          "IGVF Knowledge Graph multi-hop traversal"),
    "portal-kg":        ("igvfagent.portal_to_kg_skill",
                          "Portal → local Knowledge Graph ETL"),
    "ref":              ("igvfagent.reference_skill",
                          "Literature retrieval / validation / study design"),
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
    "warehouse":        ("igvfagent.warehouse_skill",
                          "Central DuckDB warehouse — Silver tier for the "
                          "IGVF integrated data layer"),
}

# All reserved namespaces are now wired. Kept as a (currently empty)
# anchor so future placeholders have a documented home.
RESERVED: "dict[str, str]" = {}

# Introspection commands wired in step 2: lists the LLM provider router's
# registered backends, the tool registry the agent runtime will see, and
# (step 5) the locally installed Ollama models when the daemon is running.
INTROSPECTION = ("backends", "tools", "models")

# Top-level commands wired in step 3 (`ask`) and step 4 (`ui`).
TOP_LEVEL = ("ask", "ui")


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
    print("Available skills:")
    width = max(len(name) for name in SKILLS)
    for name, (_, doc) in SKILLS.items():
        print(f"  {name:{width}}  {doc}")
    print()
    print("Top-level commands:")
    print(f"  {'ask':{width}}  Natural-language ReAct agent (LLM-driven)")
    print(f"  {'ui':{width}}  Launch the Streamlit browser UI")
    print()
    print("Introspection:")
    print(f"  {'backends':{width}}  List configured LLM provider backends")
    print(f"  {'tools':{width}}  List the tool registry the agent runtime sees")
    print(f"  {'models':{width}}  List Ollama models installed on the local daemon")
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


def main(argv: Optional["list[str]"] = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("--help", "-h", "help"):
        _print_help()
        return 0
    if args[0] in ("--version", "-V", "version"):
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
    if skill not in SKILLS:
        sys.stderr.write(f"unknown skill: {skill}\n")
        sys.stderr.write("Run `igvfagent --help` for the skill list.\n")
        return 2

    module_name, _ = SKILLS[skill]
    try:
        mod = importlib.import_module(module_name)
    except ImportError as exc:
        sys.stderr.write(f"failed to load skill `{skill}`: {exc}\n")
        return 1
    if not hasattr(mod, "main"):
        sys.stderr.write(f"skill `{skill}` has no main() entrypoint\n")
        return 1

    # Hand off to the underlying skill's argparse parser. We synthesize
    # argv[0] so the skill's own --help text shows the friendly skill
    # name rather than the dispatcher path.
    sys.argv = [f"igvfagent {skill}"] + args[1:]
    try:
        rc = mod.main()
    except SystemExit as exc:  # argparse + sys.exit propagation
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        # Non-int code (e.g. sys.exit("some error message")) — print to
        # stderr and use exit 1, matching Python's own default behavior.
        sys.stderr.write(f"{code}\n")
        return 1
    return int(rc) if isinstance(rc, int) else 0


def _run_top_level(skill: str, args: "list[str]") -> int:
    if skill == "ask":
        return _run_ask(args)
    if skill == "ui":
        return _run_ui(args)
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
