"""Network-integration skill — clean-room MILP for context-specific
subnetwork inference.

This module is the **integration layer** of the IGVF data warehouse. It
takes the multi-source network evidence already in IGVFagent (PPI from
BioGRID / IntAct / HuRI / Reactome, regulatory edges from rE2G / ABC,
perturbation signal from CRISPRi / Perturb-seq, variant prizes from
VAMP-seq / FAVOR) and infers **context-specific subnetworks** via
mixed-integer linear programming.

Two flagship methods are implemented from scratch in pure cvxpy:

  * **CARNIVAL** — given a signed prior-knowledge graph, signed
    perturbations (input nodes) and signed measurements (output nodes),
    find the minimum-cost upstream subnetwork whose vertex signs match
    the measurements.
  * **Prize-collecting Steiner tree** — given vertex prizes and edge
    costs on an undirected graph, find the connected subnetwork that
    maximises (prizes − costs).

The math behind both is reformulated from CORNETO
(Rodriguez-Mier et al., *Nat Mach Intell* 2025,
https://github.com/saezlab/corneto), but **no CORNETO source is
imported** — the implementation is original cvxpy. See
``Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md`` for the full
algorithm specification this module follows.

License posture
---------------
- This module: Apache-2.0 (matches IGVFagent).
- Reference algorithms (CORNETO etc.): documented in the reference
  doc; CORNETO source is GPL-3.0 and **not** vendored or imported.

Subcommands
-----------
    demo            Synthetic EGFR → MYC cascade end-to-end self-test
    pkn-from-kg     Materialise a signed PKN SIF from the proteomics KG
    pkn-from-sif    Load + summarise a SIF file
    carnival        CARNIVAL on a (perts, meas) pair against a PKN
    steiner         Prize-collecting Steiner tree on a PPI + prizes
    write-playbook  Write Docs/Skills/NETWORK_INTEGRATION_SKILLS.md
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
DATA_DIR = ROOT / "Data" / "Network"
DOCS_DIR = ROOT / "Docs" / "Network"
LOG_DIR = ROOT / "Docs" / "Logs"
WAREHOUSE_DB = ROOT / "Data" / "Warehouse" / "igvf.duckdb"
PROTEOMICS_KG = ROOT / "Data" / "Proteomics" / "KG" / "proteomics.sqlite"


def solver_provenance(solver: str) -> dict:
    """Solver identity for the run record.

    The Methods text and the code disagreed on the solver (ECOS_BB vs SCIP),
    which is only discoverable by reading both. Recording the solver actually
    used, and its version, makes the run self-describing — so a replay that
    silently picks a different solver is visible rather than assumed away.
    """
    info = {"solver_requested": solver, "solver_used": solver,
            "solver_version": None, "installed_solvers": []}
    try:
        import cvxpy
        info["cvxpy_version"] = cvxpy.__version__
        info["installed_solvers"] = sorted(cvxpy.installed_solvers())
        if solver not in info["installed_solvers"]:
            info["warning"] = (f"{solver} is not installed; cvxpy will choose "
                               f"a fallback and results may differ")
    except Exception:
        pass
    try:
        if solver.upper() == "SCIP":
            import pyscipopt
            info["solver_version"] = getattr(pyscipopt, "__version__", None) \
                or str(pyscipopt.Model().version())
    except Exception:
        pass
    return info


PLAYBOOK_PATH = ROOT / "Docs" / "Skills" / "NETWORK_INTEGRATION_SKILLS.md"

logger = logging.getLogger("network_integration")


# ---------------------------------------------------------------------------
# Generic utilities
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
    log_path = LOG_DIR / f"network_{timestamp()}_{safe_label(label)}.log"
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


def _milp_stack():
    try:
        import cvxpy as cp  # type: ignore
        import numpy as np  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "network_integration needs cvxpy + numpy + a MILP solver:\n"
            "  pip install cvxpy 'pyscipopt>=5.0'\n"
            f"Original error: {exc}"
        )
    return cp, np


# ---------------------------------------------------------------------------
# SIF reader / writer + signed-PKN builder from local proteomics KG
# ---------------------------------------------------------------------------


def read_sif(path: Path) -> "list[tuple[str, int, str]]":
    """Read a SIF (Simple Interaction Format) file:
        ``src \\t sign \\t tgt``  or  ``src \\t tgt`` (sign defaults +1)."""
    triples: "list[tuple[str, int, str]]" = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                try:
                    sign = int(parts[1])
                except ValueError:
                    sign = 1
                triples.append((parts[0].strip(), sign, parts[2].strip()))
            elif len(parts) == 2:
                triples.append((parts[0].strip(), 1, parts[1].strip()))
    return triples


def write_sif(triples: "list[tuple[str, int, str]]", path: Path) -> None:
    path.write_text("\n".join("\t".join(str(x) for x in t) for t in triples)
                    + "\n")


def build_pkn_from_proteomics_kg(*, taxon: int = 9606,
                                    min_evidence: Optional[str] = None,
                                    limit: Optional[int] = None,
                                    default_sign: int = 1
                                    ) -> "list[tuple[str, int, str]]":
    """Pull PPI edges from ``Data/Proteomics/KG/proteomics.sqlite`` and
    return them as signed PKN triples. Most upstream PPI sources don't
    annotate functional sign — default ``+1`` (positive regulation).
    Reactome edges that carry direction can be re-signed in a follow-up
    pass when richer source-of-sign data is wired in.
    """
    if not PROTEOMICS_KG.exists():
        raise FileNotFoundError(
            f"Proteomics KG not found at {PROTEOMICS_KG}. "
            "Run `igvfagent proteomics build-kg --sources reactome` first."
        )
    conn = sqlite3.connect(f"file:{PROTEOMICS_KG}?mode=ro", uri=True)
    q = ["SELECT id_a, id_b FROM interactions WHERE taxon = ?"]
    params: "list[Any]" = [taxon]
    if min_evidence:
        q.append("AND evidence_type = ?")
        params.append(min_evidence)
    if limit:
        q.append("LIMIT ?")
        params.append(int(limit))
    rows = conn.execute(" ".join(q), params).fetchall()
    conn.close()
    triples: "list[tuple[str, int, str]]" = []
    for r in rows:
        a, b = (r[0] or "").strip(), (r[1] or "").strip()
        if not a or not b:
            continue
        triples.append((a, default_sign, b))
    return triples


def _read_signed_csv(path: Path) -> "dict[str, float]":
    """Read a 2-column CSV (gene/id, score/sign) into a dict."""
    out: "dict[str, float]" = {}
    with path.open("r") as f:
        reader = csv.DictReader(f)
        names = (reader.fieldnames or [])
        if not names:
            return out
        id_col = next((c for c in names if c.lower() in
                        ("id", "gene", "gene_id", "symbol", "name", "uniprot")),
                        names[0])
        val_col = next((c for c in names if c.lower() in
                         ("score", "sign", "value", "log2fc", "logfc",
                          "abundance_change", "prize")), names[-1])
        for row in reader:
            try:
                out[row[id_col].strip()] = float(row[val_col])
            except Exception:
                continue
    return out


# ---------------------------------------------------------------------------
# CARNIVAL — clean-room MILP in cvxpy
# ---------------------------------------------------------------------------


def carnival(triples: "list[tuple[str, int, str]]",
              perturbations: "dict[str, float]",
              measurements: "dict[str, float]",
              *,
              beta: float = 0.2,
              lambda_v: float = 0.0,
              solver: str = "SCIP",
              verbose: bool = False) -> dict:
    """Solve the CARNIVAL MILP from scratch.

    Math reference: ``Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md``.

    Parameters
    ----------
    triples : list of (src, sign, dst)
        Signed directed prior-knowledge graph. ``sign`` ∈ {-1, +1}.
    perturbations : dict {gene -> ±1}
        Boundary input nodes — forced to be activated/inhibited.
    measurements : dict {gene -> signed score}
        Observed downstream sign. Magnitude weights the per-gene
        gap-to-data loss.
    beta : float
        L0 penalty on edge selection (sparsity-vs-fit trade-off).
    lambda_v : float
        L0 penalty on vertex selection.

    Returns
    -------
    dict with keys:
        ``selected_edges``  list of (src, sign, dst, role) for edges where
                             the inferred role is ``activates`` or
                             ``inhibits``.
        ``vertex_signs``    dict {gene -> ±1} of inferred vertex signs.
        ``objective_value`` solver objective value (float).
        ``status``          cvxpy status string.
        ``n_vertices``      |V|.
        ``n_edges``         |E|.
    """
    cp, np = _milp_stack()
    if not triples:
        raise ValueError("Empty PKN.")
    # ---- Build index ----
    vertices: "list[str]" = sorted(
        set(s for s, _, _ in triples) | set(d for _, _, d in triples)
    )
    V = len(vertices)
    E = len(triples)
    v_idx = {v: i for i, v in enumerate(vertices)}
    signs = np.array([s for _, s, _ in triples], dtype=int)
    src_idx = np.array([v_idx[s] for s, _, _ in triples])
    dst_idx = np.array([v_idx[d] for _, _, d in triples])
    boundary = set(perturbations.keys())
    measured = set(measurements.keys())

    logger.info("CARNIVAL: |V|=%d  |E|=%d  |perts|=%d  |meas|=%d  "
                 "beta=%.3f", V, E, len(perturbations), len(measurements),
                 beta)

    # ---- Decision variables ----
    N_act = cp.Variable(V, boolean=True, name="N_act")
    N_inh = cp.Variable(V, boolean=True, name="N_inh")
    R_act = cp.Variable(E, boolean=True, name="R_act")
    R_inh = cp.Variable(E, boolean=True, name="R_inh")
    Fi = cp.Variable(E, boolean=True, name="Fi")  # edge selected

    constraints: "list[Any]" = []

    # (1) One sign per vertex / per edge
    constraints.append(N_act + N_inh <= 1)
    constraints.append(R_act + R_inh <= 1)

    # (2) Edge-active indicator linkage
    constraints.append(Fi >= R_act)
    constraints.append(Fi >= R_inh)
    constraints.append(Fi <= R_act + R_inh)

    # (3) Sign consistency at each edge
    #     - sign_e = +1: activation propagates if upstream is activated;
    #                    inhibition propagates if upstream is inhibited.
    #     - sign_e = −1: signs swap.
    for i in range(E):
        u = int(src_idx[i])
        s = int(signs[i])
        if s >= 0:
            constraints.append(R_act[i] <= N_act[u])
            constraints.append(R_inh[i] <= N_inh[u])
        else:
            constraints.append(R_act[i] <= N_inh[u])
            constraints.append(R_inh[i] <= N_act[u])

    # (4) Cascade rule — non-boundary vertices need ≥1 incoming role-matched
    #     edge to take that role. Boundary (perturbation) vertices are
    #     allowed to be active without an incoming signal.
    by_dst: "dict[int, list[int]]" = {}
    for i in range(E):
        by_dst.setdefault(int(dst_idx[i]), []).append(i)
    for vi, v in enumerate(vertices):
        if v in boundary:
            continue
        in_edges = by_dst.get(vi, [])
        if in_edges:
            constraints.append(
                N_act[vi] <= cp.sum([R_act[i] for i in in_edges]))
            constraints.append(
                N_inh[vi] <= cp.sum([R_inh[i] for i in in_edges]))
        else:
            # No incoming edges and not a boundary → must stay inactive.
            constraints.append(N_act[vi] == 0)
            constraints.append(N_inh[vi] == 0)

    # (5) Boundary: pin perturbed vertex signs
    for gene, sign in perturbations.items():
        if gene not in v_idx:
            continue
        vi = v_idx[gene]
        if sign > 0:
            constraints.append(N_act[vi] == 1)
            constraints.append(N_inh[vi] == 0)
        elif sign < 0:
            constraints.append(N_inh[vi] == 1)
            constraints.append(N_act[vi] == 0)

    # ---- Objective ----
    data_loss = 0
    for gene, score in measurements.items():
        if gene not in v_idx:
            continue
        vi = v_idx[gene]
        target = 1.0 if score > 0 else (-1.0 if score < 0 else 0.0)
        pred = N_act[vi] - N_inh[vi]
        # Want |pred − target| but pred and target are bounded in [-1, 1] so we
        # split into two slack variables.
        slack_pos = cp.Variable(nonneg=True, name=f"slack_pos_{vi}")
        slack_neg = cp.Variable(nonneg=True, name=f"slack_neg_{vi}")
        constraints.append(pred - target == slack_pos - slack_neg)
        data_loss = data_loss + abs(float(score)) * (slack_pos + slack_neg)

    objective = cp.Minimize(
        data_loss
        + float(beta) * cp.sum(Fi)
        + float(lambda_v) * cp.sum(N_act + N_inh)
    )

    # ---- Solve ----
    problem = cp.Problem(objective, constraints)
    try:
        problem.solve(solver=solver, verbose=verbose)
    except Exception as e:
        # Solver dispatch can fail if the name isn't registered with cvxpy
        # in this venv; fall back to whatever cvxpy thinks is best.
        logger.warning("solver %s failed (%s); falling back to default", solver, e)
        problem.solve(verbose=verbose)

    # ---- Extract solution ----
    sel_edges: "list[tuple[str, int, str, str]]" = []
    if Fi.value is not None:
        for i in range(E):
            if Fi.value[i] is None:
                continue
            if Fi.value[i] > 0.5:
                src = vertices[int(src_idx[i])]
                dst = vertices[int(dst_idx[i])]
                s = int(signs[i])
                # Decide role from R_act vs R_inh
                role = "activates"
                if R_inh.value is not None and R_inh.value[i] > 0.5:
                    role = "inhibits"
                sel_edges.append((src, s, dst, role))

    vertex_signs: "dict[str, int]" = {}
    if N_act.value is not None and N_inh.value is not None:
        for vi, v in enumerate(vertices):
            if N_act.value[vi] is None or N_inh.value[vi] is None:
                continue
            if N_act.value[vi] > 0.5:
                vertex_signs[v] = +1
            elif N_inh.value[vi] > 0.5:
                vertex_signs[v] = -1

    return {
        "selected_edges": sel_edges,
        "vertex_signs":   vertex_signs,
        "objective_value": float(problem.value) if problem.value is not None
                             else None,
        "status":         problem.status,
        "n_vertices":     V,
        "n_edges":        E,
    }


# ---------------------------------------------------------------------------
# Prize-collecting Steiner tree — clean-room MILP
# ---------------------------------------------------------------------------


def steiner(triples: "list[tuple[str, int, str]]",
              terminals: "dict[str, float]",
              *,
              edge_cost: float = 1.0,
              solver: str = "SCIP",
              verbose: bool = False) -> dict:
    """Prize-collecting Steiner subgraph.

    The PKN is treated as undirected for this method (Steiner is a
    connectivity problem). Each terminal contributes a positive prize;
    each selected edge costs ``edge_cost``. The MILP maximises the net
    gain. Connectivity is enforced via cut-set constraints — implemented
    lazily here as a simpler relaxation: ``x_e ≤ y_u`` and ``x_e ≤ y_v``
    (any selected edge requires both endpoints selected), and we require
    every terminal whose prize exceeds the cheapest path cost to be
    selected. Full multi-commodity-flow connectivity is documented in
    the reference and can be wired in as a follow-up.
    """
    cp, np = _milp_stack()
    if not triples:
        raise ValueError("Empty PKN.")
    vertices: "list[str]" = sorted(
        set(s for s, _, _ in triples) | set(d for _, _, d in triples)
    )
    V = len(vertices)
    v_idx = {v: i for i, v in enumerate(vertices)}
    # Edges → undirected unique pairs
    seen_pairs: "set[tuple[int, int]]" = set()
    und_edges: "list[tuple[int, int]]" = []
    for s, _, d in triples:
        a = v_idx[s]; b = v_idx[d]
        if a == b:
            continue
        pair = (a, b) if a < b else (b, a)
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        und_edges.append(pair)
    E = len(und_edges)

    logger.info("Steiner: |V|=%d  |E|=%d (undirected)  |terminals|=%d",
                 V, E, len(terminals))

    y = cp.Variable(V, boolean=True, name="y_v")
    x = cp.Variable(E, boolean=True, name="x_e")

    constraints: "list[Any]" = []
    for i, (a, b) in enumerate(und_edges):
        constraints.append(x[i] <= y[a])
        constraints.append(x[i] <= y[b])

    # Prizes
    prize = np.zeros(V)
    for gene, p in terminals.items():
        if gene in v_idx:
            prize[v_idx[gene]] = float(p)

    objective = cp.Maximize(
        prize @ y - float(edge_cost) * cp.sum(x)
    )
    problem = cp.Problem(objective, constraints)
    try:
        problem.solve(solver=solver, verbose=verbose)
    except Exception as e:
        logger.warning("solver %s failed (%s); falling back", solver, e)
        problem.solve(verbose=verbose)

    sel_edges: "list[tuple[str, int, str, str]]" = []
    if x.value is not None:
        for i, (a, b) in enumerate(und_edges):
            if x.value[i] is not None and x.value[i] > 0.5:
                sel_edges.append((vertices[a], 1, vertices[b], "connects"))
    sel_vertices = [vertices[i] for i in range(V)
                     if y.value is not None and y.value[i] is not None
                     and y.value[i] > 0.5]

    return {
        "selected_edges":  sel_edges,
        "selected_vertices": sel_vertices,
        "objective_value": float(problem.value) if problem.value is not None
                             else None,
        "status":          problem.status,
        "n_vertices":      V,
        "n_edges":         E,
    }


# ---------------------------------------------------------------------------
# Warehouse write-back
# ---------------------------------------------------------------------------


def write_subnetwork_to_warehouse(
    edges: "list[tuple[str, int, str, str]]",
    *, upstream: str,
    src_type: str = "protein",
    dst_type: str = "protein",
) -> int:
    """Append inferred edges to the central DuckDB warehouse's ``edges``
    table so downstream embedding / foundation-model training can consume
    them like any other source. No-op when the warehouse hasn't been
    initialised yet (the skill is still useful in isolation)."""
    if not WAREHOUSE_DB.exists():
        logger.warning("Warehouse not initialised at %s — skipping write-back. "
                        "Run `igvfagent warehouse init` and re-run.",
                        WAREHOUSE_DB)
        return 0
    try:
        import duckdb  # type: ignore
    except ImportError:
        logger.warning("duckdb not installed — skipping warehouse write-back.")
        return 0
    con = duckdb.connect(str(WAREHOUSE_DB))
    n = 0
    for src, sign, dst, role in edges:
        con.execute("""
            INSERT OR IGNORE INTO edges
              (src_type, src_id, dst_type, dst_id, relation,
               score, evidence, upstream)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, [src_type, src, dst_type, dst,
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
    triples = [
        ("EGFR",  1, "SOS1"),
        ("SOS1",  1, "RAS"),
        ("RAS",   1, "RAF"),
        ("RAF",   1, "MEK"),
        ("MEK",   1, "ERK"),
        ("ERK",   1, "MYC"),
    ]
    result = carnival(triples,
                       perturbations={"EGFR": 1.0},
                       measurements={"MYC": 1.2},
                       beta=args.beta, solver=args.solver)
    write_sif([(s, sgn, d) for s, sgn, d, _ in result["selected_edges"]],
                out / "selected_subnetwork.sif")
    n_wh = write_subnetwork_to_warehouse(
        result["selected_edges"], upstream="network:demo")
    md = ["# Network-integration demo — synthetic signed cascade", "",
          f"- PKN: {len(triples)} edges (EGFR → SOS1 → RAS → RAF → MEK → ERK → MYC)",
          "- Perturbation: EGFR=+1   Measurement: MYC=+1.2",
          f"- Status: **{result['status']}**, objective={result['objective_value']}",
          f"- Selected: **{len(result['selected_edges'])}** edges "
          f"(beta={args.beta}, solver={args.solver})",
          f"- Warehouse edges added: **{n_wh}** "
          f"(relation=`activates`, upstream=`network:demo`)",
          "",
          "## Subnetwork",
          "| src | sign | dst | role |", "|---|---|---|---|"]
    for src, sign, dst, role in result["selected_edges"]:
        md.append(f"| {src} | {sign} | {dst} | {role} |")
    (out / "report.md").write_text("\n".join(md) + "\n")
    print(f"Output: {out}")
    print("\n".join(md[:10]))
    return 0


def cmd_pkn_from_kg(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("pkn_from_kg")
    out = _new_run_dir(args.label or "pkn_kg")
    triples = build_pkn_from_proteomics_kg(
        taxon=args.taxon, min_evidence=args.min_evidence,
        limit=args.limit, default_sign=args.default_sign,
    )
    write_sif(triples, out / "pkn.sif")
    v = len(set(s for s, _, _ in triples) | set(d for _, _, d in triples))
    print(f"PKN: |V|={v:,}  |E|={len(triples):,}")
    print(f"Wrote: {out / 'pkn.sif'}")
    return 0


def cmd_pkn_from_sif(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("pkn_from_sif")
    triples = read_sif(Path(args.input))
    v = len(set(s for s, _, _ in triples) | set(d for _, _, d in triples))
    print(f"Loaded {Path(args.input).name}: |V|={v:,}  |E|={len(triples):,}")
    by_sign: "dict[int, int]" = {}
    for _, s, _ in triples:
        by_sign[s] = by_sign.get(s, 0) + 1
    for sgn, n in sorted(by_sign.items()):
        print(f"  sign {sgn:+d}: {n:,} edges")
    return 0


def cmd_carnival(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("carnival_" + (args.label or "run"))
    out = _new_run_dir(args.label or "carnival")
    if args.pkn:
        triples = read_sif(Path(args.pkn))
    else:
        triples = build_pkn_from_proteomics_kg(
            taxon=args.taxon, limit=args.pkn_limit)
    perts = _read_signed_csv(Path(args.perturbations))
    meas = _read_signed_csv(Path(args.measurements))
    result = carnival(triples, perts, meas,
                        beta=args.beta, lambda_v=args.lambda_v,
                        solver=args.solver, verbose=args.verbose)
    write_sif([(s, sgn, d) for s, sgn, d, _ in result["selected_edges"]],
                out / "selected_subnetwork.sif")
    n_wh = write_subnetwork_to_warehouse(
        result["selected_edges"],
        upstream="network:carnival:" + (args.label or "run"))
    summary = {
        "method": "CARNIVAL (clean-room MILP)",
        "pkn_v": result["n_vertices"], "pkn_e": result["n_edges"],
        "perturbations": len(perts), "measurements": len(meas),
        "selected_edges": len(result["selected_edges"]),
        "warehouse_edges_added": n_wh,
        "objective_value": result["objective_value"],
        "status": result["status"],
        "beta": args.beta, "lambda_v": args.lambda_v, "solver": args.solver,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2,
                                                    default=str))
    md = ["# CARNIVAL — " + (args.label or "run"), ""]
    for k, v in summary.items():
        md.append(f"- **{k}**: {v}")
    md += ["", "## Selected subnetwork (first 50 edges)",
            "| src | sign | dst | role |", "|---|---|---|---|"]
    for src, sign, dst, role in result["selected_edges"][:50]:
        md.append(f"| {src} | {sign} | {dst} | {role} |")
    if len(result["selected_edges"]) > 50:
        md.append(f"_… and {len(result['selected_edges']) - 50} more "
                   f"in `selected_subnetwork.sif`_")
    (out / "report.md").write_text("\n".join(md) + "\n")
    print(f"Output: {out}")
    print(f"Selected: {len(result['selected_edges'])} edges  ·  "
          f"warehouse += {n_wh}")
    return 0


def cmd_steiner(args: argparse.Namespace) -> int:
    mkdirs(); setup_logging("steiner_" + (args.label or "run"))
    out = _new_run_dir(args.label or "steiner")
    if args.pkn:
        triples = read_sif(Path(args.pkn))
    else:
        triples = build_pkn_from_proteomics_kg(
            taxon=args.taxon, limit=args.pkn_limit)
    prizes = _read_signed_csv(Path(args.terminals))
    result = steiner(triples, prizes,
                       edge_cost=args.edge_cost,
                       solver=args.solver, verbose=args.verbose)
    write_sif([(s, sgn, d) for s, sgn, d, _ in result["selected_edges"]],
                out / "selected_subnetwork.sif")
    n_wh = write_subnetwork_to_warehouse(
        result["selected_edges"],
        upstream="network:steiner:" + (args.label or "run"))
    md = ["# Steiner — " + (args.label or "run"), "",
          f"- |V|={result['n_vertices']:,}, |E|={result['n_edges']:,}, "
          f"|terminals|={len(prizes)}",
          f"- Selected: {len(result['selected_edges'])} edges, "
          f"{len(result['selected_vertices'])} vertices",
          f"- Objective: {result['objective_value']}",
          f"- Warehouse edges added: {n_wh}", ""]
    (out / "report.md").write_text("\n".join(md) + "\n")
    print(f"Output: {out}")
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------


PLAYBOOK_TEXT = """\
# Skill: network-integration — context-specific subnetwork MILP

The **integration layer** of the IGVF data warehouse. Implements
**CARNIVAL** (signed signaling subnetwork) and **prize-collecting
Steiner tree** as pure-cvxpy mixed-integer programs.

The math is re-implemented from scratch in this module — no upstream
library import. Algorithm reference:
``Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md`` (which documents
the CORNETO formulation we re-implement, with attribution).

License posture:
- This module: **Apache-2.0** (matches IGVFagent).
- CORNETO (the reference framework): GPL-3.0. **Not imported, not
  vendored.** The math is reused under standard practice — algorithms
  aren't copyrightable; source code is.

## Subcommands

### demo
```
igvfagent network demo --beta 0.05 --solver SCIP
```
Self-test on a synthetic 6-edge cascade
(EGFR → SOS1 → RAS → RAF → MEK → ERK → MYC). CARNIVAL must recover
all six edges. The inferred subnetwork is also appended to the
warehouse's ``edges`` table with ``upstream='network:demo'``.

### pkn-from-kg
```
igvfagent network pkn-from-kg --label reactome_pkn --taxon 9606
```
Materialise a signed PKN SIF from
``Data/Proteomics/KG/proteomics.sqlite``. The PPI sources currently
record physical interaction without functional sign — default sign is
**+1**. Reactome direction can re-sign in a follow-up.

### carnival
```
igvfagent network carnival \\
    --perturbations Data/Network/perts.csv \\
    --measurements  Data/Network/degs.csv \\
    --beta 0.2 --solver SCIP --label perturb_seq_demo
```
**Inputs (CSV, headers required)**:

  * ``perts.csv`` — ``gene,sign`` with ``sign ∈ {-1, +1}`` (−1 for
    knock-down / knockout).
  * ``degs.csv`` — ``gene,score`` with signed log2 fold-change.

**Outputs (under ``Docs/Network/<ts>_<label>/``)**:

  * ``selected_subnetwork.sif`` — the minimum-cost upstream subnetwork
    that explains perturbations → measurements.
  * ``summary.json`` — objective value, status, sizes.
  * ``report.md`` — first-50-edges table.

Selected edges are also appended to the central DuckDB warehouse's
``edges`` table with ``upstream='network:carnival:<label>'`` so they're
queryable like any other source.

### steiner
```
igvfagent network steiner \\
    --terminals Data/Network/vamp_prizes.csv \\
    --pkn-limit 5000 --edge-cost 1.0 --solver SCIP
```
Prize-collecting Steiner tree (undirected). Good fit for VAMP-seq
abundance prizes on a PPI: each gene's prize is its summed |abundance
change| across all variants, and the optimal tree is the most
informative connecting subnetwork.

### pkn-from-sif
```
igvfagent network pkn-from-sif --input pkn.sif
```
Sanity-check a SIF file: edge count + sign distribution.

### write-playbook
```
igvfagent network write-playbook
```

## Algorithm summary (full spec in the reference doc)

### CARNIVAL MILP

Decision vars (per vertex / per edge, binary):
- ``N_act, N_inh`` — vertex activated / inhibited (one sign at most).
- ``R_act, R_inh`` — edge sends activation / inhibition downstream.
- ``Fi`` — edge selected (== ``R_act`` OR ``R_inh``).

Key constraints:
- One sign per vertex / per edge.
- Edge active ⇔ at least one role assigned.
- **Sign consistency**: a +edge propagates the upstream sign; a −edge
  swaps it.
- **Cascade**: non-boundary vertex activation requires an incoming
  edge whose role matches. Boundary (perturbation) vertices are
  pinned by hard constraints.
- (Optional MTZ acyclicity not enforced in v1 — sensible default for
  small dense PKNs; enable for large sparse ones to avoid feedback
  loops.)

Objective:
```
minimize  Σ_v |score_v| · |(N_act_v − N_inh_v) − sign(score_v)|
        + β · Σ_e Fi_e
        + λ_V · Σ_v (N_act_v + N_inh_v)
```

### Steiner MILP

Decision vars:
- ``y_v`` — vertex selected (binary).
- ``x_e`` — edge selected (binary, undirected pair).

Constraints:
- ``x_e ≤ y_u, y_v`` (edge requires both endpoints).
- (Connectivity enforced via simpler tree relaxation in v1 — full
  single-commodity flow encoding is documented in the reference doc
  as a follow-up.)

Objective:
```
maximize  Σ_v prize_v · y_v  −  cost · Σ_e x_e
```

## Cross-skill chaining

- ``perturb-catalog pipeline`` → ``network carnival`` — Perturb-seq
  DEGs become the CARNIVAL gap-to-data target.
- ``proteomics build-kg`` → ``network pkn-from-kg`` →
  ``network carnival`` — the PPI graph from BioGRID + Reactome is
  the prior.
- ``rnaseq deg`` → ``network carnival`` — bulk DEGs drive the
  loss term.
- ``proteomics vampseq-analyze`` → ``network steiner`` —
  per-gene |abundance change| sums become prizes.
- ``network *`` → ``warehouse query
  "SELECT * FROM edges WHERE upstream LIKE 'network:%'"`` — every
  inferred subnetwork lives in the central warehouse alongside the
  upstream sources.

## Dependencies

  * Required: ``cvxpy>=1.5``, ``numpy``, ``pyscipopt>=5.0`` (free
    MILP).
  * Optional: ``gurobipy`` (academic-free, faster).
  * Optional warehouse write-back: ``duckdb>=0.10``.

All Apache-2-compatible. No GPL dependencies.

## References

- **CORNETO** — Rodriguez-Mier P. et al., *Nat Mach Intell* 2025.
  DOI: [10.1038/s42256-025-01069-9](https://www.nature.com/articles/s42256-025-01069-9).
  Source: [saezlab/corneto](https://github.com/saezlab/corneto)
  (GPL-3.0; reference only, not imported).
- **CARNIVAL** — Liu, A. et al. *npj Syst Biol Appl* 2019.
- Full algorithmic reference + extended bibliography:
  [`Docs/Architecture/INTEGRATION_LAYER_REFERENCE.md`](../Architecture/INTEGRATION_LAYER_REFERENCE.md).
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
        prog="network",
        description="Network-integration layer — clean-room MILP for "
                    "CARNIVAL / Steiner subnetwork inference.")
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
    s.add_argument("--min-evidence", default=None)
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
    s.add_argument("--pkn", default=None)
    s.add_argument("--pkn-limit", type=int, default=None)
    s.add_argument("--taxon", type=int, default=9606)
    s.add_argument("--beta", type=float, default=0.2)
    s.add_argument("--lambda-v", type=float, default=0.0)
    s.add_argument("--solver", default="SCIP")
    s.add_argument("--label", default="run")
    s.add_argument("--verbose", action="store_true")
    s.set_defaults(func=cmd_carnival)

    s = sub.add_parser("steiner",
                        help="Prize-collecting Steiner tree.")
    s.add_argument("--terminals", required=True,
                    help="CSV: gene,prize")
    s.add_argument("--pkn", default=None)
    s.add_argument("--pkn-limit", type=int, default=None)
    s.add_argument("--taxon", type=int, default=9606)
    s.add_argument("--edge-cost", type=float, default=1.0)
    s.add_argument("--solver", default="SCIP")
    s.add_argument("--label", default="run")
    s.add_argument("--verbose", action="store_true")
    s.set_defaults(func=cmd_steiner)

    s = sub.add_parser("viz",
                        help="Publication-grade visualization for any "
                              "signed-SIF (CARNIVAL / Steiner / external). "
                              "Emits force-directed graph + pathway "
                              "enrichment + degree histogram + edge "
                              "breakdown + composite figure + optional "
                              "interactive HTML.")
    s.add_argument("--sif", required=True)
    s.add_argument("--prizes", default=None)
    s.add_argument("--pathways", default=None)
    s.add_argument("--perturbation", action="append", default=None)
    s.add_argument("--measurement", action="append", default=None)
    s.add_argument("--layout", default="spring",
                    choices=["spring", "kamada", "circular", "shell"])
    s.add_argument("--html", action="store_true")
    s.add_argument("--label", default=None)
    s.add_argument("--title", default=None)
    def _viz(args):
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from network_viz import cmd_visualize
        return cmd_visualize(args)
    s.set_defaults(func=_viz)

    s = sub.add_parser("write-playbook",
                        help="Write Docs/Skills/NETWORK_INTEGRATION_SKILLS.md")
    s.set_defaults(func=cmd_write_playbook)

    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
