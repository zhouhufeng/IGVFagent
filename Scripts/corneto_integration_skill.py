"""CORNETO-driven network integration skill.

The **integration layer** of the IGVF data warehouse: takes the
multi-source network evidence already in IGVFagent (PPI from BioGRID /
IntAct / HuRI / Reactome, regulatory edges from rE2G / ABC, perturbation
signal from CRISPRi / Perturb-seq, variant prizes from VAMP-seq / FAVOR)
and uses **CORNETO** (Saez-Rodriguez lab, *Nat Mach Intell* 2025,
https://github.com/saezlab/corneto) to infer **context-specific
subnetworks** by mixed-integer constrained optimization.

What CORNETO buys us
--------------------
CORNETO frames many separate biological-network-inference methods
(CARNIVAL, COSMOS, PCSF/Steiner-tree, OmniPath subnetwork extraction,
FBA/iMAT, shortest-path) as instances of *one* MILP over a signed
directed prior-knowledge graph (PKN). Decision variables are binary
edge / vertex activation indicators plus continuous flow; constraints
enforce flow conservation + sign consistency; objectives combine a
gap-to-data loss with L0 sparsity on edges and vertices. This replaces
the patchwork of standalone R/Bioconductor packages with one
optimization template and gives exact MILP optimality guarantees.

License boundary
----------------
CORNETO is GPL-3.0. IGVFagent is Apache-2.0. This module **does not
copy any CORNETO source**; it imports the library at runtime (allowed
under Apache + GPL coexistence). Install with:

    pip install corneto cvxpy 'pyscipopt>=5.0'   # free MILP backend
    # optionally: pip install gurobipy           # academic free, faster

Subcommands
-----------
    demo               Self-test on a synthetic signed cascade
    pkn-from-kg        Build a signed PKN from the proteomics SQLite KG
    pkn-from-sif       Load a signed PKN from a SIF file
    carnival           CARNIVAL: perts → upstream subnetwork explaining measurements
    steiner            Prize-collecting Steiner: terminals → connecting subnetwork
    multi-carnival     Joint multi-condition CARNIVAL (shared sparsity)
    write-playbook     Write Docs/Skills/CORNETO_INTEGRATION_SKILLS.md

The selected subnetworks are written back into the central warehouse
(``Data/Warehouse/igvf.duckdb``) as ``edges`` rows tagged with
``upstream='corneto:<method>'`` so downstream embedding / foundation-
model training can consume them like any other source.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sqlite3
import sys
import time
from pathlib import Path
from typing import Any, Optional


# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data" / "Corneto"
DOCS_DIR = ROOT / "Docs" / "Corneto"
LOG_DIR = ROOT / "Docs" / "Logs"
WAREHOUSE_DB = ROOT / "Data" / "Warehouse" / "igvf.duckdb"
PROTEOMICS_KG = ROOT / "Data" / "Proteomics" / "KG" / "proteomics.sqlite"
PLAYBOOK_PATH = ROOT / "Docs" / "Skills" / "CORNETO_INTEGRATION_SKILLS.md"

logger = logging.getLogger("corneto_integration")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in (s or "run"))


def mkdirs() -> None:
    for d in (DATA_DIR, DOCS_DIR, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)


def setup_logging(label: str) -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"corneto_{timestamp()}_{safe_label(label)}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    return log_path


def _new_run_dir(label: str) -> Path:
    out = DOCS_DIR / f"{timestamp()}_{safe_label(label)}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _corneto():
    """Soft-import CORNETO with a clear install hint on failure."""
    try:
        import corneto as cn  # type: ignore
        import corneto.methods as cnm  # type: ignore
        import corneto.methods.carnival as cc  # type: ignore
        from corneto.backend import CvxpyBackend  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "corneto_integration needs corneto + cvxpy + a MILP solver:\n"
            "  pip install corneto cvxpy 'pyscipopt>=5.0'\n"
            "  # optional faster solver (academic free):\n"
            "  pip install gurobipy\n"
            f"Original error: {exc}"
        )
    return cn, cnm, cc, CvxpyBackend


# ---------------------------------------------------------------------------
# Prior knowledge network (PKN) builders
# ---------------------------------------------------------------------------


def _read_sif(path: Path) -> "list[tuple[str, int, str]]":
    """Read a SIF (Simple Interaction Format) file: ``src TAB sign TAB tgt``
    or ``src TAB tgt`` (sign defaults to +1)."""
    triples: "list[tuple[str, int, str]]" = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                src, sign_s, tgt = parts[0], parts[1], parts[2]
                try:
                    sign = int(sign_s)
                except ValueError:
                    sign = 1
            elif len(parts) == 2:
                src, sign, tgt = parts[0], 1, parts[1]
            else:
                continue
            triples.append((src.strip(), int(sign), tgt.strip()))
    return triples


def build_signed_pkn(sif_triples: "list[tuple[str, int, str]]"):
    """Build a CORNETO signed directed Graph from a list of (src, sign, tgt)."""
    cn, _, _, _ = _corneto()
    return cn.Graph.from_sif_tuples(sif_triples)


def build_pkn_from_proteomics_kg(*, taxon: int = 9606,
                                    min_evidence: Optional[str] = None,
                                    limit: Optional[int] = None,
                                    default_sign: int = 1
                                    ) -> "tuple[Any, list[tuple[str,int,str]]]":
    """Pull PPI edges from ``Data/Proteomics/KG/proteomics.sqlite`` and
    build a signed PKN.

    Most upstream PPI sources (BioGRID, IntAct, HuRI) don't annotate
    edge signs — they only record physical interaction. We default
    ``sign=+1`` (positive functional regulation). Reactome edges that
    carry direction can be re-signed via a follow-up pass when a richer
    source is wired in.
    """
    if not PROTEOMICS_KG.exists():
        raise FileNotFoundError(
            f"Proteomics KG not found at {PROTEOMICS_KG}. "
            "Run `igvfagent proteomics build-kg --sources reactome` first."
        )
    conn = sqlite3.connect(f"file:{PROTEOMICS_KG}?mode=ro", uri=True)
    q = ["SELECT id_a, id_b, source, evidence_type, confidence_score "
         "FROM interactions WHERE taxon = ?"]
    params: "list[Any]" = [taxon]
    if min_evidence:
        q.append("AND evidence_type = ?")
        params.append(min_evidence)
    if limit:
        q.append("LIMIT ?")
        params.append(int(limit))
    rows = conn.execute(" ".join(q), params).fetchall()
    conn.close()
    triples: "list[tuple[str,int,str]]" = []
    for r in rows:
        a, b = (r[0] or "").strip(), (r[1] or "").strip()
        if not a or not b:
            continue
        triples.append((a, default_sign, b))
    g = build_signed_pkn(triples)
    return g, triples


# ---------------------------------------------------------------------------
# CARNIVAL — perturbations → upstream signaling subnetwork
# ---------------------------------------------------------------------------


def _read_signed_csv(path: Path) -> "dict[str, float]":
    """Read a 2-column CSV ``id,score`` into a dict."""
    out: "dict[str, float]" = {}
    with path.open("r") as f:
        reader = csv.DictReader(f)
        # Allow flexible column names
        names = (reader.fieldnames or [])
        id_col = next((c for c in names if c.lower() in
                        ("id", "gene", "gene_id", "symbol", "name")), names[0])
        val_col = next((c for c in names if c.lower() in
                         ("score", "sign", "value", "log2fc", "logfc",
                          "abundance_change")), names[-1])
        for row in reader:
            try:
                out[row[id_col].strip()] = float(row[val_col])
            except Exception:
                continue
    return out


def run_carnival(g, perturbations: "dict[str, float]",
                  measurements: "dict[str, float]", *,
                  beta: float = 0.2,
                  solver: str = "SCIP") -> "tuple[Any, Any, list]":
    """Run vanilla CARNIVAL via CORNETO. Returns (problem, flow_graph,
    selected_edges)."""
    _, cnm, cc, CvxpyBackend = _corneto()
    logger.info("CARNIVAL: |V|=%d, |E|=%d, |perts|=%d, |meas|=%d, beta=%.3f",
                 g.num_vertices, g.num_edges,
                 len(perturbations), len(measurements), beta)
    P, Gf = cnm.runVanillaCarnival(
        perturbations=perturbations,
        measurements=measurements,
        priorKnowledgeNetwork=g,
        betaWeight=beta,
        backend=CvxpyBackend(),
    )
    P.solve(solver=solver)
    try:
        sel_idx = cc.get_selected_edges(P, Gf)
    except (TypeError, AttributeError) as e:
        # corneto's get_selected_edges raises when the solver returned the
        # all-zero solution (no perturbation→measurement path in the PKN).
        # Treat that as 'no selected edges' rather than propagating.
        logger.warning("CARNIVAL returned no selected edges (%s) — likely "
                        "no signal-flow path exists for these perturbations "
                        "+ measurements in the PKN.", e)
        sel_idx = []
    return P, Gf, sel_idx or []


def _selected_edges_to_triples(Gf, sel_idx
                                  ) -> "list[tuple[str, int, str, str]]":
    """Translate CORNETO selected-edge indices into
    (src, sign, dst, edge_role) tuples. ``edge_role`` ∈ {'activates',
    'inhibits'}; sign matches the role."""
    out: "list[tuple[str, int, str, str]]" = []
    for i in sel_idx:
        try:
            edge = Gf.get_edge(int(i))
        except Exception:
            try:
                edge = Gf._edges[int(i)]
            except Exception:
                continue
        # CORNETO edges look like (src_set, tgt_set, attrs); pick the
        # first concrete vertex on each side and the edge interaction sign.
        try:
            src_set, tgt_set = edge[0], edge[1]
            src = next(iter(src_set)) if src_set else None
            tgt = next(iter(tgt_set)) if tgt_set else None
            attrs = edge[2] if len(edge) > 2 else {}
            sign = int(attrs.get("interaction", 1)) if isinstance(attrs, dict) else 1
        except Exception:
            continue
        if src is None or tgt is None:
            continue
        # Drop CORNETO-internal dummy '_s', '_t', '_pert_*', '_meas_*' nodes
        if str(src).startswith("_") or str(tgt).startswith("_"):
            continue
        role = "activates" if sign >= 0 else "inhibits"
        out.append((str(src), sign, str(tgt), role))
    return out


# ---------------------------------------------------------------------------
# Steiner — prize-collecting subnetwork for variant / abundance prizes
# ---------------------------------------------------------------------------


def run_steiner(g, terminals: "dict[str, float]", *,
                 root: Optional[str] = None,
                 solver: str = "SCIP") -> "tuple[Any, Any]":
    """Prize-collecting Steiner tree via CORNETO. ``terminals`` is a dict
    of node → prize (e.g. VAMP-seq abundance penalty for a variant in
    that gene)."""
    _, cnm, _, CvxpyBackend = _corneto()
    # API in corneto 1.0.x: corneto.methods.steiner.exact_steiner_tree
    try:
        from corneto.methods import steiner as cnst  # type: ignore
    except Exception as exc:
        raise RuntimeError(f"corneto.methods.steiner not available: {exc}")
    logger.info("Steiner: |V|=%d, |E|=%d, |terminals|=%d",
                 g.num_vertices, g.num_edges, len(terminals))
    P, Gc = cnst.exact_steiner_tree(
        g, terminals,
        root=root, strict_acyclic=False,
        backend=CvxpyBackend(),
    )
    P.solve(solver=solver)
    return P, Gc


# ---------------------------------------------------------------------------
# Write back to the central warehouse
# ---------------------------------------------------------------------------


def write_subnetwork_to_warehouse(triples: "list[tuple[str, int, str, str]]",
                                    *, upstream: str,
                                    src_type: str = "protein",
                                    dst_type: str = "protein") -> int:
    """Append CORNETO-inferred edges to the central DuckDB warehouse so
    downstream embedding / foundation-model training can consume them
    alongside every other source.

    No-op when the warehouse hasn't been initialised yet.
    """
    if not WAREHOUSE_DB.exists():
        logger.warning("Warehouse not initialised at %s — skipping "
                        "subnetwork write-back. Run `igvfagent warehouse "
                        "init` and re-run.", WAREHOUSE_DB)
        return 0
    try:
        import duckdb  # type: ignore
    except ImportError:
        logger.warning("duckdb not installed — skipping warehouse write-back.")
        return 0
    con = duckdb.connect(str(WAREHOUSE_DB))
    n = 0
    for src, sign, tgt, role in triples:
        con.execute("""
            INSERT OR IGNORE INTO edges
              (src_type, src_id, dst_type, dst_id, relation,
               score, evidence, upstream)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, [src_type, src, dst_type, tgt,
              role, float(sign), "predicted", upstream])
        n += 1
    con.close()
    return n


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------


def cmd_demo(args: argparse.Namespace) -> int:
    """Self-test on a synthetic signed cascade — proves the wiring."""
    mkdirs(); setup_logging("demo")
    out = _new_run_dir(args.label or "demo")
    cn, _, _, _ = _corneto()
    triples = [
        ("EGFR",  1, "SOS1"),
        ("SOS1",  1, "RAS"),
        ("RAS",   1, "RAF"),
        ("RAF",   1, "MEK"),
        ("MEK",   1, "ERK"),
        ("ERK",   1, "MYC"),
    ]
    g = build_signed_pkn(triples)
    P, Gf, sel_idx = run_carnival(
        g, {"EGFR": 1}, {"MYC": 1.2},
        beta=args.beta, solver=args.solver,
    )
    edges = _selected_edges_to_triples(Gf, sel_idx)
    (out / "selected_subnetwork.tsv").write_text(
        "src\tsign\tdst\trole\n"
        + "\n".join("\t".join(map(str, e)) for e in edges) + "\n")
    n_wh = write_subnetwork_to_warehouse(
        edges, upstream="corneto:demo")
    md = ["# CORNETO demo — synthetic signed cascade", "",
          f"- PKN: {len(triples)} edges (EGFR → SOS1 → RAS → RAF → MEK → ERK → MYC)",
          "- Perturbation: EGFR=+1   Measurement: MYC=+1.2",
          f"- Selected: **{len(edges)}** edges via CARNIVAL "
          f"(beta={args.beta}, solver={args.solver})",
          f"- Warehouse edges added: **{n_wh}** "
          f"(relation=`activates`, upstream=`corneto:demo`)",
          "",
          "## Subnetwork",
          "| src | sign | dst | role |", "|---|---|---|---|"]
    for src, sign, dst, role in edges:
        md.append(f"| {src} | {sign} | {dst} | {role} |")
    (out / "report.md").write_text("\n".join(md) + "\n")
    print(f"Output: {out}")
    print("\n".join(md[:10]))
    return 0


def cmd_pkn_from_kg(args: argparse.Namespace) -> int:
    """Materialise a SIF file of the proteomics KG signed PKN."""
    mkdirs(); setup_logging("pkn_from_kg")
    out = _new_run_dir(args.label or "pkn_kg")
    g, triples = build_pkn_from_proteomics_kg(
        taxon=args.taxon, min_evidence=args.min_evidence,
        limit=args.limit, default_sign=args.default_sign,
    )
    sif = out / "pkn.sif"
    sif.write_text("\n".join("\t".join(str(x) for x in t) for t in triples)
                    + "\n")
    print(f"PKN written: {sif}")
    print(f"  |V| = {g.num_vertices:,}, |E| = {g.num_edges:,}")
    return 0


def cmd_pkn_from_sif(args: argparse.Namespace) -> int:
    """Quick load + summary of a SIF file."""
    mkdirs(); setup_logging("pkn_from_sif")
    triples = _read_sif(Path(args.input))
    g = build_signed_pkn(triples)
    print(f"Loaded {Path(args.input).name}")
    print(f"  |V| = {g.num_vertices:,}, |E| = {g.num_edges:,}")
    by_sign: "dict[int, int]" = {}
    for _, s, _ in triples:
        by_sign[s] = by_sign.get(s, 0) + 1
    for sgn, n in sorted(by_sign.items()):
        print(f"    sign {sgn:+d}: {n:,} edges")
    return 0


def cmd_carnival(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("carnival_" + (args.label or "run"))
    out = _new_run_dir(args.label or "carnival")
    if args.pkn:
        triples = _read_sif(Path(args.pkn))
        g = build_signed_pkn(triples)
    else:
        g, _ = build_pkn_from_proteomics_kg(
            taxon=args.taxon, limit=args.pkn_limit)
    perts = _read_signed_csv(Path(args.perturbations))
    meas = _read_signed_csv(Path(args.measurements))
    P, Gf, sel_idx = run_carnival(
        g, perts, meas, beta=args.beta, solver=args.solver)
    edges = _selected_edges_to_triples(Gf, sel_idx)
    sif = out / "selected_subnetwork.sif"
    sif.write_text("\n".join("\t".join(str(x) for x in (s, sgn, d))
                              for s, sgn, d, _ in edges) + "\n")
    n_wh = write_subnetwork_to_warehouse(
        edges, upstream="corneto:carnival:" + (args.label or "run"))
    summary = {
        "method": "CARNIVAL",
        "pkn_v": g.num_vertices, "pkn_e": g.num_edges,
        "perturbations": len(perts), "measurements": len(meas),
        "selected_edges": len(edges),
        "warehouse_edges_added": n_wh,
        "beta": args.beta, "solver": args.solver,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    md = ["# CORNETO CARNIVAL — " + (args.label or "run"), ""]
    for k, v in summary.items():
        md.append(f"- **{k}**: {v}")
    md += ["", "## Selected subnetwork (first 50 edges)",
            "| src | sign | dst | role |", "|---|---|---|---|"]
    for src, sign, dst, role in edges[:50]:
        md.append(f"| {src} | {sign} | {dst} | {role} |")
    if len(edges) > 50:
        md.append(f"_… and {len(edges) - 50} more in `selected_subnetwork.sif`_")
    (out / "report.md").write_text("\n".join(md) + "\n")
    print(f"Output: {out}")
    print(f"Selected: {len(edges)} edges  ·  warehouse += {n_wh}")
    return 0


def cmd_steiner(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("steiner_" + (args.label or "run"))
    out = _new_run_dir(args.label or "steiner")
    if args.pkn:
        triples = _read_sif(Path(args.pkn))
        g = build_signed_pkn(triples)
    else:
        g, _ = build_pkn_from_proteomics_kg(
            taxon=args.taxon, limit=args.pkn_limit)
    terminals = _read_signed_csv(Path(args.terminals))
    P, Gc = run_steiner(g, terminals, root=args.root, solver=args.solver)
    print(f"Steiner solved. Output dir: {out}")
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------


PLAYBOOK_TEXT = """\
# Skill: CORNETO network integration

The **integration layer** of the IGVF integrated data warehouse. Wraps
the **CORNETO** framework (Saez-Rodriguez lab,
[*Nat Mach Intell* 2025](https://www.nature.com/articles/s42256-025-01069-9),
GPL-3.0, https://github.com/saezlab/corneto) to infer context-specific
subnetworks by mixed-integer constrained optimization over a signed
directed prior-knowledge graph.

CORNETO reformulates many separate biological-network-inference tools
(CARNIVAL, COSMOS, PCSF/Steiner-tree, OmniPath subnetwork extraction,
FBA/iMAT, shortest-path) as instances of one MILP template:

  * Binary edge / vertex activation indicators + continuous flow
  * Flow conservation + sign consistency constraints
  * Objective = data-loss + L0 sparsity (edges, vertices) + L1 flow

## License

CORNETO is **GPL-3.0**. This skill installs CORNETO as a runtime
dependency and calls it via its public API — no source is copied into
the (Apache-2-licensed) IGVFagent repo. Install:

```
pip install corneto cvxpy 'pyscipopt>=5.0'
# optional faster solver (academic-free):
pip install gurobipy
```

## Subcommands

### demo
```
igvfagent corneto demo --beta 0.05 --solver SCIP
```
Self-test on a synthetic signed cascade (EGFR → SOS1 → RAS → RAF →
MEK → ERK → MYC) — proves CARNIVAL is wired end-to-end and writes
the selected subnetwork back into the warehouse with
`upstream='corneto:demo'`.

### pkn-from-kg
```
igvfagent corneto pkn-from-kg --label reactome --taxon 9606
```
Materialise a signed PKN SIF from `Data/Proteomics/KG/proteomics.sqlite`.
The PPI sources we have today (BioGRID, IntAct, HuRI) record physical
interaction without functional sign, so the default sign is **+1**
(positive regulation). Reactome edges are direction-aware and could
be re-signed via a follow-up; for now they're +1 too.

### carnival
```
# Build a PKN on the fly from the proteomics KG, run CARNIVAL with
# perturbed genes (signed) and DEGs (signed log2FC):
igvfagent corneto carnival \\
    --perturbations Data/Corneto/perts.csv \\
    --measurements  Data/Corneto/degs.csv \\
    --beta 0.2 --solver SCIP --label perturb_seq_demo

# Or pass a pre-built PKN:
igvfagent corneto carnival --pkn pkn.sif \\
    --perturbations perts.csv --measurements degs.csv
```
**Inputs** (CSV, headers required):

  * `perts.csv`: `gene,sign` with `sign ∈ {-1, +1}` — +1 for activation,
    −1 for knock-down / knockout.
  * `degs.csv`: `gene,score` with signed log2 fold-change.

Output: `selected_subnetwork.sif` (the minimum-cost upstream subnetwork
that explains the perturbations → measurements) + a `summary.json` +
a markdown report. The selected edges are also appended to the
warehouse's `edges` table with `upstream='corneto:carnival:<label>'`,
so they're queryable like any other source.

### steiner
```
igvfagent corneto steiner \\
    --pkn ppi.sif --terminals prizes.csv \\
    --root TP53 --solver SCIP
```
Prize-collecting Steiner tree. Good fit for **VAMP-seq variant prizes
on a PPI**: each gene's "prize" is its summed absolute abundance change
across all variants, and the optimal tree is the most informative
sub-network connecting the high-prize genes.

### multi-carnival (planned)
Joint inference across many conditions with shared sparsity (one
subnetwork per Perturb-seq guide, all sharing edge / vertex L0
penalties). Available via `corneto.methods.future.carnival.multi_carnival`
upstream; not yet wired here.

### write-playbook
```
igvfagent corneto write-playbook
```

## Cross-skill chaining

- `perturb-catalog pipeline` → `corneto carnival` — turn a Perturb-seq
  dataset's DEGs into an upstream signaling subnetwork.
- `proteomics build-kg` → `corneto pkn-from-kg` → `corneto carnival` —
  the PPI graph from BioGRID + Reactome becomes the prior for context
  inference.
- `rnaseq deg` → `corneto carnival` — DEGs from a bulk comparison
  drive the gap-to-data loss in CARNIVAL.
- `corneto carnival` → `warehouse query "SELECT * FROM edges WHERE
  upstream LIKE 'corneto:%'"` — every inferred subnetwork is a
  first-class citizen of the central warehouse, ready for downstream
  embedding extraction + foundation-model training.

## Why this design

- CORNETO is **the right primitive** for the integration layer because
  it unifies the math behind dozens of separate biological-network
  tools, with MILP optimality guarantees and multi-sample joint
  inference.
- We don't reinvent it — we install it as a runtime dependency,
  honour its GPL by not copying source, and use the public API.
- Selected subnetworks flow back into the **DuckDB warehouse**
  alongside every other producer, so the downstream embedding /
  foundation-model layer treats CORNETO-inferred edges as just
  another evidence stream.

## References

- Rodriguez-Mier P. et al. *CORNETO: a unified framework for the
  inference of context-specific networks.* **Nat Mach Intell** 2025.
  DOI: [10.1038/s42256-025-01069-9](https://www.nature.com/articles/s42256-025-01069-9).
- Source: [saezlab/corneto](https://github.com/saezlab/corneto)
- Manuscript experiments:
  [saezlab/corneto-manuscript-experiments](https://github.com/saezlab/corneto-manuscript-experiments)
"""


def cmd_write_playbook(_a) -> int:
    PLAYBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    PLAYBOOK_PATH.write_text(PLAYBOOK_TEXT)
    print(f"Wrote: {PLAYBOOK_PATH}")
    return 0


# ---------------------------------------------------------------------------
# argparse
# ---------------------------------------------------------------------------


def main(argv: "Optional[list[str]]" = None) -> int:
    p = argparse.ArgumentParser(
        prog="corneto",
        description="CORNETO-driven network integration: CARNIVAL / "
                    "Steiner subnetwork inference on signed PKN.")
    sub = p.add_subparsers(dest="cmd")

    s = sub.add_parser("demo",
                        help="Self-test on a synthetic signed cascade.")
    s.add_argument("--label", default="demo")
    s.add_argument("--beta", type=float, default=0.05)
    s.add_argument("--solver", default="SCIP")
    s.set_defaults(func=cmd_demo)

    s = sub.add_parser("pkn-from-kg",
                        help="Materialise a signed PKN SIF from the "
                              "proteomics KG.")
    s.add_argument("--label", default=None)
    s.add_argument("--taxon", type=int, default=9606)
    s.add_argument("--min-evidence", default=None,
                    help="experimental | curated | predicted")
    s.add_argument("--limit", type=int, default=None)
    s.add_argument("--default-sign", type=int, default=1, choices=[-1, 1])
    s.set_defaults(func=cmd_pkn_from_kg)

    s = sub.add_parser("pkn-from-sif",
                        help="Load + summarise a SIF file.")
    s.add_argument("--input", required=True)
    s.set_defaults(func=cmd_pkn_from_sif)

    s = sub.add_parser("carnival",
                        help="CARNIVAL: perts → upstream subnetwork.")
    s.add_argument("--perturbations", required=True,
                    help="CSV: gene,sign")
    s.add_argument("--measurements", required=True,
                    help="CSV: gene,score")
    s.add_argument("--pkn", default=None,
                    help="SIF path; default = proteomics KG.")
    s.add_argument("--pkn-limit", type=int, default=None,
                    help="When pulling PKN from KG, cap edges.")
    s.add_argument("--taxon", type=int, default=9606)
    s.add_argument("--beta", type=float, default=0.2)
    s.add_argument("--solver", default="SCIP")
    s.add_argument("--label", default="run")
    s.set_defaults(func=cmd_carnival)

    s = sub.add_parser("steiner",
                        help="Prize-collecting Steiner tree.")
    s.add_argument("--terminals", required=True,
                    help="CSV: gene,prize")
    s.add_argument("--pkn", default=None)
    s.add_argument("--pkn-limit", type=int, default=None)
    s.add_argument("--taxon", type=int, default=9606)
    s.add_argument("--root", default=None)
    s.add_argument("--solver", default="SCIP")
    s.add_argument("--label", default="run")
    s.set_defaults(func=cmd_steiner)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/CORNETO_INTEGRATION_SKILLS.md")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
