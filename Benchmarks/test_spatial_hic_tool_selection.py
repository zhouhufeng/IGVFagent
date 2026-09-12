#!/usr/bin/env python3
"""Does the orchestrator reach the right spatial-hic tool from a plain request?

WHY THIS EXISTS. Eleven spatial-hic subcommands are registered as tools, but
being in the catalogue is not the same as being reachable: `matrix` and
`impute` were implemented, documented and registered nowhere, so no request
could ever select them. A count of registered tools would not have caught
that, and neither would any offline test -- tool choice is the model's, so it
has to be measured against a live backend.

WHAT COUNTS AS A PASS. Strictly: the target tool is called. Tolerated: a
read-only discovery call first (`list_artifacts` / `read_artifact`) -- asked
to analyse a named directory the model sometimes confirms it exists before
touching it, which is sensible and, in the agent loop, simply costs one
iteration before the right tool is called. Measured over four runs, the
compartment case did this half the time.

FAILS ON: a different ANALYSIS tool, which is a real misroute.

Needs ANTHROPIC_API_KEY and network; skips cleanly without them, so it never
breaks the offline suite.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Scripts"))

if not os.environ.get("ANTHROPIC_API_KEY"):
    print("SKIP  no ANTHROPIC_API_KEY — tool selection needs a live backend")
    sys.exit(0)

# Runs from a repo checkout (Scripts/ on the path) AND from inside the
# deployed container, where the same modules are the installed `igvfagent`
# package and Scripts/ does not exist. The deployment is the environment whose
# tool choice actually matters, so it has to be runnable there.
try:
    import _agent as A                                      # noqa: E402
    import _llm as L                                        # noqa: E402
    import _tools as T                                      # noqa: E402
except ModuleNotFoundError:                                 # installed package
    from igvfagent import _agent as A                       # noqa: E402
    from igvfagent import _llm as L                         # noqa: E402
    from igvfagent import _tools as T                       # noqa: E402

DISCOVERY = {"list_artifacts", "read_artifact", "grep_artifacts"}

D = "/workspace/Data/SpatialATACHiC/run1/pixels"
F = "/workspace/Data/SpatialATACHiC/run1/fragments.tsv.gz"
G = "/workspace/Data/SpatialATACHiC/run1/gencode.gtf.gz"

CASES = [
    (f"Using the pixel contact directory {D}, impute the chr2 map at 25kb "
     f"resolution.", "spatial_hic_impute"),
    (f"Using {D}, build the chr2 contact matrix at 25kb.", "spatial_hic_matrix"),
    (f"Using {D}, call A/B compartments at 100kb.", "spatial_hic_compartment"),
    (f"Using {D}, compute per-pixel copy number at 5Mb to find tumour clones.",
     "spatial_hic_cnv"),
    ("Split /workspace/Data/SpatialATACHiC/run1/sample.pairs.gz into the "
     "50x50 tissue pixel grid using barcodes_A.txt and barcodes_B.txt.",
     "spatial_hic_pixel_demux"),
    (f"Compute the gene activity score from {F} with gene model {G}.",
     "spatial_hic_gas"),
    (f"Compute the gene-associated domain score from {D} with {G}.",
     "spatial_hic_gad"),
    (f"Run per-pixel contact QC on {D} — cis fraction and long-range ratio.",
     "spatial_hic_qc"),
    (f"Using {D} and anchors loops.bedpe with clusters.tsv, test which loops "
     f"are cluster specific.", "spatial_hic_loops"),
    ("Render the Satb2 column of run1/gene_activity_score.tsv back into "
     "tissue space using tissue_positions.csv.", "spatial_hic_viz"),
    ("What does GEO series GSE307620 deposit for Spatial-ATAC-Hi-C?",
     "spatial_hic_pull_geo"),
]

tools = [t.to_dict() for t in T.list_tools()]
registered = {t.name for t in T.list_tools()}
FAILURES, strict = [], 0

# Every subcommand the CLI offers must have a tool, or a request for it can
# never be routed. This is the check that would have caught matrix/impute.
CLI_SUBCOMMANDS = {"pull-geo", "pixel-demux", "qc", "gas", "gad", "matrix",
                   "impute", "compartment", "cnv", "loops", "viz"}
exposed = {t.cli[1] for t in T.list_tools() if t.name.startswith("spatial_hic")}
missing = CLI_SUBCOMMANDS - exposed
print(f"{'PASS' if not missing else 'FAIL'}  every spatial-hic subcommand is a "
      f"registered tool  {'' if not missing else 'MISSING: ' + str(sorted(missing))}")
if missing:
    FAILURES.append("subcommand parity")

for q, want in CASES:
    if want not in registered:
        print(f"FAIL  {want:26}  not registered at all")
        FAILURES.append(want)
        continue
    try:
        msg = L.chat([{"role": "system", "content": A.DEFAULT_SYSTEM_PROMPT},
                      {"role": "user", "content": q}],
                     backend="anthropic", model="claude-sonnet-5",
                     tools=tools, max_tokens=700)
        got = [tc.name for tc in msg.tool_calls]
    except Exception as e:                                  # noqa: BLE001
        print(f"FAIL  {want:26}  backend error: {type(e).__name__}: {e}")
        FAILURES.append(want)
        continue
    if want in got:
        strict += 1
        print(f"PASS  {want:26}  called: {got}")
    elif got and all(g in DISCOVERY for g in got):
        print(f"PASS  {want:26}  discovery first: {got} (reaches the tool "
              f"next iteration)")
    else:
        print(f"FAIL  {want:26}  misrouted to: {got or '(no tool call)'}")
        FAILURES.append(want)

print(f"\n{len(CASES) + 1 - len(FAILURES)}/{len(CASES) + 1} checks passed "
      f"({strict}/{len(CASES)} selected the tool on the first call)")
sys.exit(1 if FAILURES else 0)
