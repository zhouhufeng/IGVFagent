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
}

# Reserved namespaces for the upcoming agent + UI layers (steps 3-4 in
# the shipping plan). They print a "not yet wired" message until those
# steps land — keeping the CLI surface stable for users now.
RESERVED = {
    "ask": "Natural-language agent (LLM-driven). Coming in step 3.",
    "ui":  "Browser UI (Streamlit). Coming in step 4.",
}

# Introspection commands wired in step 2: lists the LLM provider router's
# registered backends and the tool registry the agent runtime will see.
INTROSPECTION = ("backends", "tools")


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
    print("Introspection:")
    print(f"  {'backends':{width}}  List configured LLM provider backends")
    print(f"  {'tools':{width}}  List the tool registry the agent runtime sees")
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
        return int(exc.code) if exc.code is not None else 0
    return int(rc) if isinstance(rc, int) else 0


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
    return 2


if __name__ == "__main__":
    sys.exit(main())
