"""Local Knowledge Graph visualization module.

Renders the SQLite-backed local Knowledge Graphs that IGVFagent builds:

  1. ``Data/Proteomics/KG/proteomics.sqlite``  — proteins / interactions /
     pathways / pathway_membership / igvf_evidence (built by
     ``proteomics build-kg``).
  2. ``Data/KG/portal_kg.sqlite``              — generic nodes / edges
     (built by ``portal-kg pull``).

The module is **driver-agnostic**: every backend exposes the same
``KGAccessor`` interface (``stats``, ``search``, ``neighborhood``,
``hubs``, ``sample_rows``). A Streamlit panel (``render_streamlit_panel``)
ties everything together — schema browser, search box, network plot,
stats. The same accessors can also be used from a CLI or notebook.

Visualization uses matplotlib + networkx (already in the analysis
extras) so this module does not introduce a new dependency.
"""

from __future__ import annotations

import json
import logging
import sqlite3
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
import os
from pathlib import Path
from typing import Any, Optional

# Streamlit width-argument compat (see Scripts/_stcompat.py). Kept as a
# dual-mode import so the module works installed or from a checkout.
try:
    from igvfagent._stcompat import fit
except Exception:  # pragma: no cover - checkout / direct-run fallback
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent))
    from _stcompat import fit

logger = logging.getLogger("kg_visualizer")

# ---------------------------------------------------------------------------
# KG discovery
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()


@dataclass
class KGDescriptor:
    """One entry in the local-KG registry."""
    key: str               # short stable id (used in UI selectors)
    label: str             # human-readable name
    path: Path             # SQLite file
    kind: str              # "proteomics" | "portal" | "generic"
    description: str       # one-line description
    enabled: bool          # exists & non-empty


def _sqlite_table_names(p: Path) -> "list[str]":
    if not p.exists():
        return []
    try:
        conn = sqlite3.connect(f"file:{p}?mode=ro", uri=True)
        try:
            rows = conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' "
                "AND name NOT LIKE 'sqlite_%'"
            ).fetchall()
            return [r[0] for r in rows]
        finally:
            conn.close()
    except Exception as e:
        logger.warning("Could not read %s: %s", p, e)
        return []


def discover_local_kgs(root: Optional[Path] = None) -> "list[KGDescriptor]":
    root = (root or PROJECT_ROOT).resolve()
    out: "list[KGDescriptor]" = []

    prot = root / "Data" / "Proteomics" / "KG" / "proteomics.sqlite"
    tables = _sqlite_table_names(prot)
    out.append(KGDescriptor(
        key="proteomics",
        label="Proteomics PPI-KG",
        path=prot,
        kind="proteomics",
        description=(
            "BioGRID / IntAct / HuRI / Reactome / KEGG + IGVF protein "
            "evidence integrated into one SQLite KG by "
            "`proteomics build-kg`."
        ),
        enabled=bool(tables) and "interactions" in tables,
    ))

    # The growing integrated graph. NOTE the filename: every writer
    # (_localstore, portal_to_kg_skill, sce2g_kg_skill, variant-list) uses
    # local_kg.sqlite. This discovery previously looked for portal_kg.sqlite,
    # a name nothing in the codebase writes — so the main knowledge graph was
    # invisible in the UI while quietly accumulating on disk.
    local = root / "Data" / "KG" / "local_kg.sqlite"
    tables = _sqlite_table_names(local)
    out.append(KGDescriptor(
        key="local",
        label="IGVF integrated KG",
        path=local,
        kind="portal",
        description=(
            "The growing local graph: annotated variants, genes, regulatory "
            "elements, diseases and pathways, plus analysis provenance. Fed "
            "by `variant-list annotate`, `portal-kg pull`, `sce2g`, and every "
            "agent tool call."
        ),
        enabled=bool(tables) and "nodes" in tables and "edges" in tables,
    ))

    # Legacy path, kept so an older checkout's data is still reachable.
    portal = root / "Data" / "KG" / "portal_kg.sqlite"
    tables = _sqlite_table_names(portal)
    if tables:
        out.append(KGDescriptor(
            key="portal",
            label="IGVF Portal → KG (legacy path)",
            path=portal,
            kind="portal",
            description="Older portal_kg.sqlite mirror, if present.",
            enabled="nodes" in tables and "edges" in tables,
        ))

    return out


# ---------------------------------------------------------------------------
# Accessors
# ---------------------------------------------------------------------------


class _Accessor:
    """Common interface every KG accessor implements."""

    kind = "generic"
    label = "Generic KG"

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.conn = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True)
        self.conn.row_factory = sqlite3.Row

    def close(self) -> None:
        try:
            self.conn.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # Each subclass implements these:
    def stats(self) -> dict:                       # noqa: D401
        raise NotImplementedError

    def search(self, term: str, *, limit: int = 50) -> "list[dict]":
        raise NotImplementedError

    def neighborhood(self, center: str, *, depth: int = 1,
                      max_neighbors: int = 60) -> dict:
        """Return {nodes: [{id,kind,...}], edges: [{a,b,kind,...}]}."""
        raise NotImplementedError

    def hubs(self, *, k: int = 30) -> "list[tuple[str, int]]":
        raise NotImplementedError

    def sample_rows(self, table: str, *, limit: int = 20) -> "list[dict]":
        try:
            rows = self.conn.execute(
                f"SELECT * FROM {table} LIMIT ?", (limit,)
            ).fetchall()
            return [dict(r) for r in rows]
        except Exception as e:
            logger.warning("sample_rows(%s): %s", table, e)
            return []

    def table_names(self) -> "list[str]":
        rows = self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name NOT LIKE 'sqlite_%' ORDER BY name"
        ).fetchall()
        return [r[0] for r in rows]

    def table_count(self, table: str) -> int:
        try:
            return self.conn.execute(
                f"SELECT COUNT(*) FROM {table}"
            ).fetchone()[0]
        except Exception:
            return 0


# --- Proteomics accessor ----------------------------------------------------


class ProteomicsAccessor(_Accessor):
    kind = "proteomics"
    label = "Proteomics PPI-KG"

    def stats(self) -> dict:
        c = self.conn.cursor()
        out: "dict[str, Any]" = {}
        out["interactions_total"] = c.execute(
            "SELECT COUNT(*) FROM interactions").fetchone()[0]
        out["proteins_distinct"] = c.execute(
            "SELECT COUNT(DISTINCT id) FROM ("
            "SELECT id_a AS id FROM interactions UNION "
            "SELECT id_b AS id FROM interactions)"
        ).fetchone()[0]
        out["pathways_total"] = c.execute(
            "SELECT COUNT(*) FROM pathways").fetchone()[0]
        out["pathway_membership"] = c.execute(
            "SELECT COUNT(*) FROM pathway_membership").fetchone()[0]
        out["igvf_evidence"] = c.execute(
            "SELECT COUNT(*) FROM igvf_evidence").fetchone()[0]
        out["per_source"] = dict(c.execute(
            "SELECT source, COUNT(*) FROM interactions "
            "GROUP BY source ORDER BY 2 DESC"
        ).fetchall())
        out["per_id_type"] = dict(c.execute(
            "SELECT id_type, COUNT(*) FROM interactions GROUP BY id_type"
        ).fetchall())
        out["per_evidence_type"] = dict(c.execute(
            "SELECT evidence_type, COUNT(*) FROM interactions "
            "GROUP BY evidence_type"
        ).fetchall())
        out["per_detection_method"] = dict(c.execute(
            "SELECT detection_method, COUNT(*) FROM interactions "
            "GROUP BY detection_method "
            "ORDER BY 2 DESC LIMIT 12"
        ).fetchall())
        return out

    def search(self, term: str, *, limit: int = 50) -> "list[dict]":
        t = (term or "").strip()
        if not t:
            return []
        like = f"%{t}%"
        # Find proteins (any node touched by an interaction) that match.
        rows = self.conn.execute(
            """
            SELECT id, MAX(deg) AS degree, MAX(src) AS sources FROM (
              SELECT id_a AS id, COUNT(*) AS deg,
                       GROUP_CONCAT(DISTINCT source) AS src
                FROM interactions
                WHERE id_a LIKE ? OR id_a = ?
                GROUP BY id_a
              UNION ALL
              SELECT id_b AS id, COUNT(*) AS deg,
                       GROUP_CONCAT(DISTINCT source) AS src
                FROM interactions
                WHERE id_b LIKE ? OR id_b = ?
                GROUP BY id_b
            )
            GROUP BY id
            ORDER BY degree DESC
            LIMIT ?
            """,
            (like, t, like, t, limit),
        ).fetchall()
        return [{"id": r[0], "degree": r[1], "sources": r[2]} for r in rows]

    def neighborhood(self, center: str, *, depth: int = 1,
                      max_neighbors: int = 60) -> dict:
        if not center:
            return {"nodes": [], "edges": []}
        nodes: "dict[str, dict]" = {center: {"id": center, "kind": "center"}}
        edges: "list[dict]" = []
        frontier = {center}
        for _ in range(max(1, depth)):
            new_frontier: "set[str]" = set()
            placeholders = ",".join(["?"] * len(frontier))
            rows = self.conn.execute(
                f"""
                SELECT id_a, id_b, source, detection_method, evidence_type
                FROM interactions
                WHERE id_a IN ({placeholders}) OR id_b IN ({placeholders})
                LIMIT ?
                """,
                list(frontier) + list(frontier) + [max_neighbors * 4],
            ).fetchall()
            for r in rows:
                a, b = r[0], r[1]
                src, det, ev = r[2], r[3], r[4]
                edges.append({"a": a, "b": b, "kind": src or "?",
                                "detection": det or "", "evidence": ev or ""})
                for n in (a, b):
                    if n not in nodes:
                        nodes[n] = {"id": n, "kind": "neighbor"}
                        new_frontier.add(n)
                if len(nodes) >= max_neighbors:
                    break
            frontier = new_frontier
            if not frontier or len(nodes) >= max_neighbors:
                break
        return {"nodes": list(nodes.values()), "edges": edges}

    def hubs(self, *, k: int = 30) -> "list[tuple[str, int]]":
        return self.conn.execute(
            """
            SELECT id, COUNT(*) AS deg FROM (
              SELECT id_a AS id FROM interactions UNION ALL
              SELECT id_b AS id FROM interactions
            ) GROUP BY id ORDER BY deg DESC LIMIT ?
            """,
            (k,),
        ).fetchall()

    def pathway_membership(self, node: str) -> "list[dict]":
        rows = self.conn.execute(
            """
            SELECT pm.pathway_id, p.name, p.source
            FROM pathway_membership pm
            LEFT JOIN pathways p ON p.pathway_id = pm.pathway_id
            WHERE pm.member_id = ?
            LIMIT 50
            """,
            (node,),
        ).fetchall()
        return [{"pathway_id": r[0], "name": r[1], "source": r[2]}
                 for r in rows]


# --- Portal nodes/edges accessor -------------------------------------------


class PortalKGAccessor(_Accessor):
    kind = "portal"
    label = "IGVF Portal → KG"

    def stats(self) -> dict:
        c = self.conn.cursor()
        out: "dict[str, Any]" = {}
        out["nodes_total"] = c.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
        out["edges_total"] = c.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
        out["per_node_type"] = dict(c.execute(
            "SELECT node_type, COUNT(*) FROM nodes GROUP BY node_type "
            "ORDER BY 2 DESC LIMIT 30"
        ).fetchall())
        out["per_edge_type"] = dict(c.execute(
            "SELECT edge_type, COUNT(*) FROM edges GROUP BY edge_type "
            "ORDER BY 2 DESC LIMIT 30"
        ).fetchall())
        out["per_source"] = dict(c.execute(
            "SELECT source, COUNT(*) FROM nodes GROUP BY source"
        ).fetchall())
        return out

    def search(self, term: str, *, limit: int = 50) -> "list[dict]":
        t = (term or "").strip()
        if not t:
            return []
        like = f"%{t}%"
        rows = self.conn.execute(
            """
            SELECT id, node_type, label, source FROM nodes
            WHERE id = ? OR label LIKE ? OR id LIKE ?
            ORDER BY length(label)
            LIMIT ?
            """,
            (t, like, like, limit),
        ).fetchall()
        return [{"id": r[0], "node_type": r[1],
                  "label": r[2] or r[0], "source": r[3]} for r in rows]

    def neighborhood(self, center: str, *, depth: int = 1,
                      max_neighbors: int = 60) -> dict:
        if not center:
            return {"nodes": [], "edges": []}
        # Resolve center to id (label search supported)
        row = self.conn.execute(
            "SELECT id, node_type, label FROM nodes WHERE id=? OR label=?",
            (center, center),
        ).fetchone()
        if row is None:
            return {"nodes": [], "edges": []}
        center_id = row[0]
        nodes: "dict[str, dict]" = {
            center_id: {"id": center_id, "kind": row[1] or "node",
                          "label": row[2] or center_id, "center": True}}
        edges: "list[dict]" = []
        frontier = {center_id}
        for _ in range(max(1, depth)):
            new_frontier: "set[str]" = set()
            placeholders = ",".join(["?"] * len(frontier))
            rows = self.conn.execute(
                f"""
                SELECT from_node, to_node, edge_type, source
                FROM edges
                WHERE from_node IN ({placeholders})
                   OR to_node   IN ({placeholders})
                LIMIT ?
                """,
                list(frontier) + list(frontier) + [max_neighbors * 4],
            ).fetchall()
            for r in rows:
                a, b, etype, src = r
                edges.append({"a": a, "b": b, "kind": etype, "source": src})
                for n in (a, b):
                    if n not in nodes:
                        nrow = self.conn.execute(
                            "SELECT node_type, label FROM nodes WHERE id=?",
                            (n,),
                        ).fetchone()
                        nodes[n] = {
                            "id": n,
                            "kind": (nrow[0] if nrow else "node"),
                            "label": (nrow[1] if nrow and nrow[1] else n),
                        }
                        new_frontier.add(n)
                if len(nodes) >= max_neighbors:
                    break
            frontier = new_frontier
            if not frontier or len(nodes) >= max_neighbors:
                break
        return {"nodes": list(nodes.values()), "edges": edges}

    def hubs(self, *, k: int = 30) -> "list[tuple[str, int]]":
        return self.conn.execute(
            """
            SELECT id, deg FROM (
              SELECT from_node AS id, COUNT(*) AS deg FROM edges
                GROUP BY from_node
              UNION ALL
              SELECT to_node AS id, COUNT(*) AS deg FROM edges
                GROUP BY to_node
            ) GROUP BY id ORDER BY SUM(deg) DESC LIMIT ?
            """,
            (k,),
        ).fetchall()


def accessor_for(desc: KGDescriptor) -> _Accessor:
    if desc.kind == "proteomics":
        return ProteomicsAccessor(desc.path)
    if desc.kind == "portal":
        return PortalKGAccessor(desc.path)
    raise ValueError(f"Unknown KG kind: {desc.kind}")


# ---------------------------------------------------------------------------
# Subgraph rendering (matplotlib + networkx)
# ---------------------------------------------------------------------------

# A reproducible categorical palette, colorblind-friendly.
_KIND_PALETTE = [
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#666666", "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
    "#9467BD", "#8C564B", "#E377C2", "#BCBD22", "#17BECF",
]


def _color_for_kind(kind: str, kind_index: "dict[str, int]") -> str:
    if kind not in kind_index:
        kind_index[kind] = len(kind_index)
    return _KIND_PALETTE[kind_index[kind] % len(_KIND_PALETTE)]


def render_subgraph_figure(subgraph: dict, *,
                            center: Optional[str] = None,
                            title: Optional[str] = None,
                            figsize: tuple = (8, 8)) -> Any:
    """Render a subgraph dict (nodes + edges) as a matplotlib Figure."""
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        import networkx as nx  # type: ignore
    except Exception as e:
        raise RuntimeError(
            f"render_subgraph_figure requires matplotlib + networkx: {e}"
        )

    nodes = subgraph.get("nodes", [])
    edges = subgraph.get("edges", [])

    if not nodes:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No nodes returned for this query.\n"
                          "Try a different gene/term.",
                 ha="center", va="center", fontsize=11, alpha=0.7)
        ax.axis("off")
        return fig

    G = nx.MultiGraph()
    by_id = {n["id"]: n for n in nodes}
    for n in nodes:
        G.add_node(n["id"], **{k: v for k, v in n.items() if k != "id"})
    edge_kinds = []
    for e in edges:
        if e["a"] in by_id and e["b"] in by_id and e["a"] != e["b"]:
            G.add_edge(e["a"], e["b"], kind=e.get("kind", "?"))
            edge_kinds.append(e.get("kind", "?"))

    # Layout — spring for small graphs, kamada_kawai for medium.
    n = G.number_of_nodes()
    if n <= 25:
        pos = nx.spring_layout(G, k=1.2, seed=7)
    elif n <= 80:
        pos = nx.spring_layout(G, k=0.9, seed=7, iterations=80)
    else:
        pos = nx.spring_layout(G, k=0.6, seed=7, iterations=50)

    fig, ax = plt.subplots(figsize=figsize)
    kind_index: "dict[str, int]" = {}
    node_colors = []
    node_sizes = []
    for nid in G.nodes():
        kind = by_id[nid].get("kind", "node")
        is_center = (nid == center) or by_id[nid].get("center", False)
        node_colors.append("#E6550D" if is_center
                            else _color_for_kind(kind, kind_index))
        node_sizes.append(900 if is_center else 220)

    # Edges first (so nodes paint over them)
    edge_kind_index: "dict[str, int]" = {}
    edge_color_map = []
    for u, v, k in G.edges(keys=True):
        d = G.get_edge_data(u, v, k)
        edge_color_map.append(_color_for_kind(d.get("kind", "?"),
                                                edge_kind_index))
    nx.draw_networkx_edges(G, pos, alpha=0.5, edge_color=edge_color_map,
                            width=0.9, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                            node_size=node_sizes, alpha=0.95,
                            edgecolors="white", linewidths=0.7, ax=ax)
    labels = {nid: (by_id[nid].get("label") or nid)[:18]
               for nid in G.nodes()
               if G.degree(nid) >= max(2, n // 18)
               or nid == center
               or by_id[nid].get("center", False)}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=7, ax=ax)

    # Legend — node kinds and edge kinds, separate columns
    from matplotlib.patches import Patch  # type: ignore
    legend_handles = []
    for kind, idx in kind_index.items():
        legend_handles.append(Patch(
            facecolor=_KIND_PALETTE[idx % len(_KIND_PALETTE)],
            label=f"node: {kind}"))
    for kind, idx in edge_kind_index.items():
        legend_handles.append(Patch(
            facecolor=_KIND_PALETTE[idx % len(_KIND_PALETTE)],
            label=f"edge: {kind}"))
    if legend_handles:
        ax.legend(handles=legend_handles[:14], loc="upper right",
                   fontsize=7, framealpha=0.85)

    ax.axis("off")
    ax.set_title(title or (f"Ego graph: {center} (n={n}, e={G.number_of_edges()})"
                            if center else
                            f"Subgraph (n={n}, e={G.number_of_edges()})"),
                  fontsize=11)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Streamlit panel
# ---------------------------------------------------------------------------


def _ensure_in_session(st, key: str, default):
    if key not in st.session_state:
        st.session_state[key] = default


def render_streamlit_panel(st) -> None:
    """Drop-in Streamlit panel.

    Call from streamlit_app.py under a tab. Discovers local KGs, lets the
    user pick one, browse schema, search for a node, and renders an
    interactive subgraph figure plus a summary-stats column.
    """
    st.markdown(
        "### 🕸 Knowledge Graph Explorer\n"
        "Browse the **local** SQLite KGs that IGVFagent builds. Pick a KG, "
        "search for a node, and explore its neighborhood. The same KGs are "
        "queryable from the CLI (`igvfagent proteomics kg-stats`, "
        "`igvfagent kg gene <SYM>`)."
    )

    descs = discover_local_kgs()
    enabled = [d for d in descs if d.enabled]

    if not enabled:
        st.warning(
            "No local KGs found yet. Build one with:\n\n"
            "```\n"
            "igvfagent proteomics build-kg --sources all\n"
            "# or\n"
            "igvfagent portal-kg pull --tissue macrophage --limit 100\n"
            "```\n"
            "After the build completes, refresh this tab."
        )
        st.subheader("Expected locations")
        for d in descs:
            present = "✅" if d.path.exists() else "—"
            st.write(f"{present} **{d.label}**  ·  `{d.path}`")
            st.caption(d.description)
        return

    # KG selector
    options = {d.key: d for d in enabled}
    keys = list(options.keys())
    sel = st.selectbox(
        "Knowledge graph",
        keys,
        format_func=lambda k: options[k].label,
        key="kg_viz_kg_select",
    )
    desc = options[sel]
    st.caption(desc.description + f"  ·  `{desc.path}`")

    try:
        acc = accessor_for(desc)
    except Exception as e:
        st.error(f"Could not open {desc.label}: {e}")
        return

    try:
        # --- Top stats strip ---------------------------------------------------
        s = acc.stats()
        cols = st.columns(4)
        if desc.kind == "proteomics":
            cols[0].metric("Interactions", f"{s['interactions_total']:,}")
            cols[1].metric("Distinct proteins", f"{s['proteins_distinct']:,}")
            cols[2].metric("Pathways", f"{s['pathways_total']:,}")
            cols[3].metric("IGVF evidence files",
                            f"{s['igvf_evidence']:,}")
        else:
            cols[0].metric("Nodes", f"{s['nodes_total']:,}")
            cols[1].metric("Edges", f"{s['edges_total']:,}")
            cols[2].metric("Node types", len(s.get("per_node_type", {})))
            cols[3].metric("Edge types", len(s.get("per_edge_type", {})))

        # --- Source / type breakdowns -----------------------------------------
        with st.expander("Source / type breakdowns", expanded=False):
            if desc.kind == "proteomics":
                bc1, bc2 = st.columns(2)
                with bc1:
                    st.markdown("**Interactions per source**")
                    if s.get("per_source"):
                        st.bar_chart(s["per_source"])
                with bc2:
                    st.markdown("**Per evidence type**")
                    if s.get("per_evidence_type"):
                        st.bar_chart(s["per_evidence_type"])
                if s.get("per_detection_method"):
                    st.markdown("**Top detection methods**")
                    st.bar_chart(s["per_detection_method"])
            else:
                bc1, bc2 = st.columns(2)
                with bc1:
                    st.markdown("**Node types**")
                    if s.get("per_node_type"):
                        st.bar_chart(s["per_node_type"])
                with bc2:
                    st.markdown("**Edge types**")
                    if s.get("per_edge_type"):
                        st.bar_chart(s["per_edge_type"])

        # --- Search + neighborhood -------------------------------------------
        st.divider()
        st.subheader("Explore a node")
        _ensure_in_session(st, "kg_viz_center", "")
        col_q, col_d, col_n = st.columns([3, 1, 1])
        with col_q:
            term = st.text_input(
                "Search by gene symbol / UniProt / node-id / label",
                value=st.session_state.kg_viz_center,
                placeholder=("e.g. TP53, P04637, /genes/IGVFGE…  ↩  "
                              "or pick a hub below"),
                key="kg_viz_search",
            )
        with col_d:
            depth = st.selectbox("Depth", [1, 2, 3],
                                   index=0, key="kg_viz_depth")
        with col_n:
            cap = st.selectbox(
                "Max neighbors",
                [30, 60, 120, 250],
                index=1, key="kg_viz_max_neighbors")

        # Suggestions row — top hubs
        hubs = acc.hubs(k=12)
        if hubs:
            st.caption("Top hubs (click to use as the center node):")
            hub_cols = st.columns(min(6, len(hubs)))
            for i, (hid, deg) in enumerate(hubs[:6]):
                with hub_cols[i]:
                    if st.button(f"{hid[:14]}\n(d={deg})",
                                  key=f"kg_viz_hub_{i}"):
                        st.session_state.kg_viz_center = hid
                        st.rerun()

        if term and term != st.session_state.kg_viz_center:
            # If the user typed something but it doesn't match a top hub,
            # offer a suggestion list before drawing.
            hits = acc.search(term, limit=10)
            if hits:
                center_candidate = hits[0].get("id")
                if len(hits) > 1:
                    st.write(
                        "Search hits — first match used as center; click "
                        "another to switch:")
                    sug_cols = st.columns(min(5, len(hits)))
                    for i, h in enumerate(hits[:5]):
                        with sug_cols[i]:
                            lab = h.get("id") or h.get("label") or "?"
                            extra = (h.get("node_type")
                                      or h.get("sources") or "")
                            if st.button(
                                f"{lab[:14]}\n{extra[:18]}",
                                key=f"kg_viz_hit_{i}",
                            ):
                                st.session_state.kg_viz_center = h["id"]
                                st.rerun()
                st.session_state.kg_viz_center = center_candidate

        center = st.session_state.kg_viz_center or (
            term if term else (hubs[0][0] if hubs else None))
        if not center:
            st.info("Type a node id or click a hub above to render a "
                     "neighborhood.")
        else:
            sub = acc.neighborhood(
                center,
                depth=int(depth),
                max_neighbors=int(cap),
            )
            n_n = len(sub.get("nodes", []))
            n_e = len(sub.get("edges", []))
            st.caption(
                f"**Center:** `{center}`  ·  Depth {depth}  ·  "
                f"{n_n} nodes, {n_e} edges (cap {cap})"
            )
            fig = render_subgraph_figure(sub, center=center,
                                            title=f"{center} — depth {depth}")
            st.pyplot(fig, **fit(st.pyplot))

            # Show edge table
            with st.expander("Edges in this subgraph", expanded=False):
                if sub.get("edges"):
                    st.dataframe(sub["edges"], **fit(st.dataframe),
                                  height=240)
                else:
                    st.write("(no edges in this subgraph)")

            # Pathway annotation for the proteomics KG
            if desc.kind == "proteomics":
                pw = ProteomicsAccessor.pathway_membership(acc, center)
                if pw:
                    with st.expander(
                        f"Pathway membership for `{center}` (Reactome/KEGG)",
                        expanded=False,
                    ):
                        st.dataframe(pw, **fit(st.dataframe),
                                      height=240)

        # --- Schema browser ---------------------------------------------------
        with st.expander("Schema & sample rows", expanded=False):
            tables = acc.table_names()
            tab_sel = st.selectbox(
                "Table",
                tables,
                key="kg_viz_table_sel",
            )
            n = acc.table_count(tab_sel)
            st.caption(f"`{tab_sel}` — {n:,} rows total")
            rows = acc.sample_rows(tab_sel, limit=25)
            if rows:
                st.dataframe(rows, **fit(st.dataframe), height=320)
            else:
                st.write("(empty)")

        # --- Top hubs detail table --------------------------------------------
        with st.expander("Top hubs (degree)", expanded=False):
            big = acc.hubs(k=30)
            if big:
                st.dataframe(
                    [{"node_id": h, "degree": d} for h, d in big],
                    **fit(st.dataframe), height=320,
                )
            else:
                st.write("(no hub data)")

    finally:
        acc.close()


# ---------------------------------------------------------------------------
# CLI smoke entrypoint (debug only)
# ---------------------------------------------------------------------------


def _cli_smoke() -> int:
    """`python kg_visualizer.py` -> print stats for any local KG."""
    for d in discover_local_kgs():
        line = (f"  {'✅' if d.enabled else '—'}  {d.label}: {d.path}"
                f"  [{d.kind}]")
        print(line)
        if not d.enabled:
            continue
        try:
            with accessor_for(d) as acc:
                s = acc.stats()
                print(f"      stats: {json.dumps(s, default=str)[:300]}...")
                hubs = acc.hubs(k=5)
                print(f"      top-5 hubs: {hubs}")
        except Exception as e:
            print(f"      ERROR: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(_cli_smoke())
