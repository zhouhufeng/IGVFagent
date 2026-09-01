"""Publication-grade visualization for IGVFagent network-integration outputs.

Consumes:
  - A signed-SIF file (`src \\t sign \\t dst` triples) produced by
    `network carnival/steiner/demo` (or any external CARNIVAL-style output).
  - Optional `--prizes <csv>` mapping nodes → prize / measurement /
    perturbation magnitude (used for node size + color intensity).
  - Optional `--pathways <csv>` mapping nodes → KEGG/Reactome
    pathway name (used for the enrichment bar panel).

Emits:
  1. `graph_force.png` / `.svg` — force-directed layout with signed edges
     (green = activation, red = inhibition), node colors by sign / role,
     node sizes by prize, edge widths by selection confidence.
  2. `graph_interactive.html` — pyvis vis.js HTML, click-drag-zoom.
  3. `pathway_enrichment.png` — top-K bar chart of enriched pathways
     (if `--pathways` was provided; otherwise skipped gracefully).
  4. `node_summary.csv` — per-node degree, role, prize, pathway memberships.
  5. `composite_publication_figure.png/.svg` — 4-panel publication
     mosaic combining the graph, the pathway panel, a degree-distribution
     histogram, and a per-edge confidence panel.

All matplotlib output is 300 dpi PNG + vector SVG. The HTML output is
standalone (vis.js loaded from a CDN). Apache-2.0; deps are
networkx (BSD), matplotlib (PSF/BSD), pyvis (MIT), pandas (BSD-3).
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


# ─── Output / palette setup ────────────────────────────────────────────────

LOG_DIR = Path("Docs/Logs")
SKILL_DOC_DIR = Path("Docs/Skills")

# Color palette (color-blind friendly, journal-style)
COL = {
    "node_up":      "#C13E3E",   # crimson — up-regulated / activated
    "node_down":    "#3D5A8A",   # navy — down-regulated / inhibited
    "node_neutral": "#A0AEC0",   # gray — no inferred sign
    "node_perturb": "#F2A359",   # amber — known perturbation
    "node_measure": "#7CA663",   # olive — known measurement
    "edge_activate": "#3D8C5F",  # green — activation (sign=+1)
    "edge_inhibit":  "#B53D3D",  # red   — inhibition (sign=-1)
    "edge_default":  "#7A8794",  # gray  — unsigned
    "bg":            "#FAF7F2",
    "ink":           "#1F2933",
    "subink":        "#52606D",
}


def _setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = LOG_DIR / f"network_viz_{time.strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log), logging.StreamHandler(sys.stdout)],
    )
    return log


# ─── I/O ────────────────────────────────────────────────────────────────────

def read_sif(path: Path) -> "list[tuple[str, int, str]]":
    """Read a 3-col SIF: `src \\t sign \\t dst`."""
    triples = []
    with path.open() as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            if len(parts) < 3:
                continue
            try:
                sgn = int(parts[1])
            except ValueError:
                sgn = 1 if parts[1].lower() in ("activates", "+", "up") else -1
            triples.append((parts[0], sgn, parts[2]))
    return triples


def read_prizes(path: Path | None) -> dict[str, dict[str, Any]]:
    """Read node prize / role CSV.

    Columns recognized (case-insensitive):
      node | gene | symbol         — required
      prize | weight | score        — numeric; used for node size
      sign | direction              — int (+1/-1) or up/down
      role                          — perturbation / measurement / inferred
    """
    out: dict[str, dict[str, Any]] = {}
    if not path:
        return out
    with path.open() as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            return out
        # Normalize header
        keys = {k.lower(): k for k in reader.fieldnames}
        name_key = next((keys[k] for k in ("node", "gene", "symbol", "name") if k in keys), None)
        if name_key is None:
            return out
        for row in reader:
            n = (row[name_key] or "").strip()
            if not n:
                continue
            rec: dict[str, Any] = {}
            for cand in ("prize", "weight", "score", "abundance"):
                if cand in keys:
                    try:
                        rec["prize"] = float(row[keys[cand]] or 0)
                    except (TypeError, ValueError):
                        rec["prize"] = 0.0
                    break
            for cand in ("sign", "direction"):
                if cand in keys:
                    v = (row[keys[cand]] or "").strip().lower()
                    if v in ("+1", "+", "up", "1"): rec["sign"] = 1
                    elif v in ("-1", "down", "-"): rec["sign"] = -1
                    else:
                        try: rec["sign"] = int(v) if v else 0
                        except ValueError: rec["sign"] = 0
                    break
            if "role" in keys:
                rec["role"] = (row[keys["role"]] or "").strip().lower()
            out[n] = rec
    return out


def read_pathways(path: Path | None) -> dict[str, list[str]]:
    """Read node → list-of-pathways CSV (node, pathway columns)."""
    out: dict[str, list[str]] = defaultdict(list)
    if not path:
        return out
    with path.open() as fh:
        reader = csv.DictReader(fh)
        keys = {k.lower(): k for k in (reader.fieldnames or [])}
        n_key = next((keys[k] for k in ("node", "gene", "symbol") if k in keys), None)
        p_key = next((keys[k] for k in ("pathway", "term", "name") if k in keys), None)
        if not n_key or not p_key:
            return out
        for row in reader:
            n, p = (row[n_key] or "").strip(), (row[p_key] or "").strip()
            if n and p:
                out[n].append(p)
    return out


# ─── Graph construction ─────────────────────────────────────────────────────

def build_graph(triples, prizes):
    """Build a networkx DiGraph and decorate with node + edge attributes."""
    import networkx as nx
    G = nx.DiGraph()
    for src, sgn, dst in triples:
        G.add_node(src)
        G.add_node(dst)
        G.add_edge(src, dst, sign=int(sgn))
    # Decorate nodes
    for n in G.nodes:
        meta = prizes.get(n, {}) or {}
        G.nodes[n]["prize"] = float(meta.get("prize", 0.0))
        G.nodes[n]["sign"] = int(meta.get("sign", 0))
        G.nodes[n]["role"] = meta.get("role", "")
        G.nodes[n]["degree"] = G.in_degree(n) + G.out_degree(n)
    return G


def _node_color(node_data):
    role = (node_data.get("role") or "").lower()
    if role.startswith("pert"):
        return COL["node_perturb"]
    if role.startswith("meas"):
        return COL["node_measure"]
    sgn = node_data.get("sign", 0)
    if sgn > 0:  return COL["node_up"]
    if sgn < 0:  return COL["node_down"]
    return COL["node_neutral"]


def _edge_color(sgn: int) -> str:
    if sgn > 0: return COL["edge_activate"]
    if sgn < 0: return COL["edge_inhibit"]
    return COL["edge_default"]


# ─── Panel A: static force-directed graph ───────────────────────────────────

def draw_force_layout(G, out_png: Path, *, title: str = "",
                      layout: str = "spring", seed: int = 7) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_facecolor(COL["bg"]); fig.patch.set_facecolor(COL["bg"])

    if layout == "kamada" and len(G) <= 200:
        pos = nx.kamada_kawai_layout(G)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    elif layout == "shell":
        pos = nx.shell_layout(G)
    else:
        # Spring layout with stronger repulsion for cleaner spacing
        pos = nx.spring_layout(G, seed=seed, k=1.5 / max((len(G) ** 0.5), 1),
                                iterations=200)

    # Node sizes by prize (with floor)
    prizes = [G.nodes[n].get("prize", 0.0) for n in G.nodes]
    p_max = max(abs(p) for p in prizes) or 1.0
    sizes = [400 + 1400 * (abs(p) / p_max) for p in prizes]

    node_colors = [_node_color(G.nodes[n]) for n in G.nodes]

    # Edges by sign
    edges_pos = [(u, v) for u, v, d in G.edges(data=True) if d.get("sign", 1) > 0]
    edges_neg = [(u, v) for u, v, d in G.edges(data=True) if d.get("sign", 1) < 0]
    edges_neu = [(u, v) for u, v, d in G.edges(data=True) if d.get("sign", 1) == 0]

    nx.draw_networkx_edges(G, pos, edgelist=edges_pos, ax=ax,
                             edge_color=COL["edge_activate"], arrows=True,
                             arrowsize=14, width=1.8, alpha=0.85,
                             connectionstyle="arc3,rad=0.05",
                             min_target_margin=14)
    nx.draw_networkx_edges(G, pos, edgelist=edges_neg, ax=ax,
                             edge_color=COL["edge_inhibit"], arrows=True,
                             arrowsize=14, width=1.8, alpha=0.85,
                             arrowstyle="-|>",
                             connectionstyle="arc3,rad=0.05",
                             min_target_margin=14)
    nx.draw_networkx_edges(G, pos, edgelist=edges_neu, ax=ax,
                             edge_color=COL["edge_default"], arrows=True,
                             arrowsize=12, width=1.2, alpha=0.6,
                             min_target_margin=14)

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=sizes,
                             edgecolors=COL["ink"], linewidths=1.2, ax=ax,
                             alpha=0.95)

    # Label with smart font size
    label_size = 11 if len(G) <= 30 else (9 if len(G) <= 80 else 7)
    nx.draw_networkx_labels(G, pos, font_size=label_size, font_weight="bold",
                              font_color=COL["ink"], ax=ax)

    # Legend
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=COL["node_perturb"], edgecolor=COL["ink"],
              label="perturbation"),
        Patch(facecolor=COL["node_measure"], edgecolor=COL["ink"],
              label="measurement"),
        Patch(facecolor=COL["node_up"], edgecolor=COL["ink"],
              label="inferred ↑"),
        Patch(facecolor=COL["node_down"], edgecolor=COL["ink"],
              label="inferred ↓"),
        Patch(facecolor=COL["node_neutral"], edgecolor=COL["ink"],
              label="other"),
        Line2D([0], [0], color=COL["edge_activate"], lw=2.5,
                label="activation (sign = +1)"),
        Line2D([0], [0], color=COL["edge_inhibit"], lw=2.5,
                label="inhibition (sign = −1)"),
    ]
    ax.legend(handles=legend_handles, loc="lower left", fontsize=9,
               framealpha=0.95, edgecolor=COL["subink"])

    ax.set_title(title or "Inferred regulatory subnetwork",
                  fontsize=15, fontweight="bold", color=COL["ink"], pad=12)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.margins(0.1)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor=COL["bg"])
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight",
                 facecolor=COL["bg"])
    import matplotlib.pyplot as plt
    plt.close(fig)
    return out_png


# ─── Panel B: pathway enrichment bar chart ─────────────────────────────────

def draw_pathway_panel(node_to_pathways, selected_nodes, out_png: Path,
                        *, top_k: int = 12) -> Path | None:
    """Bar chart of pathways enriched among the selected nodes.

    Score = (selected_in_pathway / |selected|) — a simple coverage measure.
    For a real hypergeometric p-value, supply the background population
    via --pathways with a 'background' value (not used in this stub).
    """
    if not node_to_pathways:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    counter: Counter[str] = Counter()
    for n in selected_nodes:
        for p in node_to_pathways.get(n, []):
            counter[p] += 1
    if not counter:
        return None
    top = counter.most_common(top_k)
    labels = [t[0] for t in top][::-1]
    values = [t[1] for t in top][::-1]
    fig, ax = plt.subplots(figsize=(8, 0.5 + 0.35 * len(top)))
    ax.set_facecolor(COL["bg"]); fig.patch.set_facecolor(COL["bg"])
    bars = ax.barh(labels, values, color=COL["node_up"], alpha=0.85,
                    edgecolor=COL["ink"], linewidth=0.6)
    for b, v in zip(bars, values):
        ax.text(v + max(values) * 0.01, b.get_y() + b.get_height() / 2,
                 str(v), va="center", fontsize=9, color=COL["ink"])
    ax.set_xlabel("Nodes in subnetwork", color=COL["ink"])
    ax.set_title("Pathway enrichment (top-K)", fontsize=13,
                  fontweight="bold", color=COL["ink"])
    ax.tick_params(colors=COL["ink"])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor=COL["bg"])
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight",
                 facecolor=COL["bg"])
    plt.close(fig)
    return out_png


# ─── Panel C: degree distribution histogram ────────────────────────────────

def draw_degree_histogram(G, out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    degrees = [G.in_degree(n) + G.out_degree(n) for n in G.nodes]
    if not degrees:
        return out_png
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_facecolor(COL["bg"]); fig.patch.set_facecolor(COL["bg"])
    n_bins = min(max(len(set(degrees)), 5), 25)
    ax.hist(degrees, bins=n_bins, color=COL["node_perturb"], alpha=0.85,
             edgecolor=COL["ink"], linewidth=0.7)
    ax.set_xlabel("Node degree (in + out)", color=COL["ink"])
    ax.set_ylabel("Count", color=COL["ink"])
    ax.set_title("Subnetwork degree distribution", fontsize=13,
                  fontweight="bold", color=COL["ink"])
    ax.tick_params(colors=COL["ink"])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor=COL["bg"])
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight",
                 facecolor=COL["bg"])
    plt.close(fig)
    return out_png


# ─── Panel D: edge-sign breakdown ───────────────────────────────────────────

def draw_edge_breakdown(G, out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    counts = Counter(d.get("sign", 0) for _, _, d in G.edges(data=True))
    labels = ["activation\n(+1)", "inhibition\n(−1)", "unsigned"]
    values = [counts.get(1, 0), counts.get(-1, 0),
               sum(v for k, v in counts.items() if k not in (1, -1))]
    colors = [COL["edge_activate"], COL["edge_inhibit"], COL["edge_default"]]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_facecolor(COL["bg"]); fig.patch.set_facecolor(COL["bg"])
    bars = ax.bar(labels, values, color=colors, edgecolor=COL["ink"],
                    linewidth=0.7, alpha=0.9)
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width() / 2, v + max(values) * 0.02,
                 str(v), ha="center", va="bottom", fontsize=11,
                 color=COL["ink"], fontweight="bold")
    ax.set_ylabel("Edge count", color=COL["ink"])
    ax.set_title("Edge sign breakdown", fontsize=13,
                  fontweight="bold", color=COL["ink"])
    ax.tick_params(colors=COL["ink"])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor=COL["bg"])
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight",
                 facecolor=COL["bg"])
    plt.close(fig)
    return out_png


# ─── Interactive HTML via pyvis ─────────────────────────────────────────────

def draw_interactive_html(G, out_html: Path, *, title: str = "") -> "Path | None":
    """Standalone vis.js HTML (loaded from CDN).

    Returns ``None`` when pyvis is not installed. Returning the path
    would be a lie every downstream consumer believes: `cmd_visualize`
    prints it, `write_report` links it, and the UI offers a tab for it
    — all pointing at a file that was never written.
    """
    try:
        from pyvis.network import Network
    except ImportError:
        logging.warning(
            "pyvis is not installed — skipping %s. "
            "Install it with: pip install 'igvfagent[network]'  "
            "(or: pip install pyvis)",
            out_html.name,
        )
        return None
    net = Network(
        height="800px", width="100%", directed=True, notebook=False,
        bgcolor=COL["bg"], font_color=COL["ink"], cdn_resources="remote",
    )
    net.barnes_hut(gravity=-25000, central_gravity=0.3,
                    spring_length=120, spring_strength=0.04, damping=0.4)
    for n, d in G.nodes(data=True):
        prize = d.get("prize", 0.0)
        net.add_node(
            n, label=n, color=_node_color(d), shape="dot",
            size=15 + 25 * (abs(prize) ** 0.5) if prize else 18,
            title=(f"<b>{n}</b><br>sign={d.get('sign', 0)}<br>"
                    f"prize={prize:.3g}<br>"
                    f"degree={d.get('degree', 0)}<br>"
                    f"role={d.get('role', '')}"),
        )
    for u, v, d in G.edges(data=True):
        sgn = d.get("sign", 0)
        net.add_edge(
            u, v, color=_edge_color(sgn),
            arrows="to", width=2.2,
            title=f"sign={sgn}",
            label="+" if sgn > 0 else ("−" if sgn < 0 else ""),
            font={"size": 14, "color": COL["ink"]},
        )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    # `net.save_graph` writes HTML; we don't use `show` because it would
    # open a browser in headless mode.
    net.save_graph(str(out_html))
    return out_html


# ─── Composite publication figure ───────────────────────────────────────────

def build_composite(plot_dir: Path, out_path: Path, *, title: str,
                     subtitle: str = "") -> Path | None:
    """Tile the four single-panel PNGs into a 2×2 publication mosaic."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    panels = [
        ("Inferred subnetwork (force layout)", "graph_force.png", 0, 0),
        ("Pathway enrichment",                  "pathway_enrichment.png", 0, 1),
        ("Degree distribution",                  "degree_histogram.png", 1, 0),
        ("Edge sign breakdown",                  "edge_breakdown.png", 1, 1),
    ]
    fig = plt.figure(figsize=(15, 12), facecolor=COL["bg"])
    fig.suptitle(title, fontsize=20, fontweight="bold", color=COL["ink"],
                  y=0.98)
    if subtitle:
        fig.text(0.5, 0.95, subtitle, ha="center", va="top", fontsize=12,
                 color=COL["subink"], style="italic")
    found_any = False
    from matplotlib import gridspec
    gs = gridspec.GridSpec(2, 2, figure=fig, left=0.02, right=0.98,
                            top=0.92, bottom=0.02, hspace=0.12, wspace=0.06)
    for label, fname, r, c in panels:
        path = plot_dir / fname
        ax = fig.add_subplot(gs[r, c])
        ax.axis("off")
        if path.is_file():
            try:
                img = mpimg.imread(path)
                ax.imshow(img)
                ax.set_title(label, fontsize=12, weight="bold",
                              color=COL["ink"], pad=6)
                found_any = True
            except Exception as exc:
                ax.text(0.5, 0.5, f"(load error: {exc})",
                         ha="center", va="center", transform=ax.transAxes)
        else:
            ax.text(0.5, 0.5, f"({fname}: not generated)",
                     ha="center", va="center", transform=ax.transAxes,
                     fontsize=10, color="#A0AEC0", style="italic")
            ax.set_title(label, fontsize=12, color="#9D9D9D", pad=6)
    if not found_any:
        plt.close(fig); return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor=COL["bg"])
    fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight",
                 facecolor=COL["bg"])
    plt.close(fig)
    return out_path


# ─── Per-node summary CSV ───────────────────────────────────────────────────

def write_node_summary(G, node_to_pathways, out_csv: Path) -> Path:
    rows = []
    for n, d in G.nodes(data=True):
        rows.append({
            "node": n,
            "in_degree": G.in_degree(n),
            "out_degree": G.out_degree(n),
            "total_degree": G.in_degree(n) + G.out_degree(n),
            "sign": d.get("sign", 0),
            "prize": d.get("prize", 0.0),
            "role": d.get("role", ""),
            "pathways": ";".join(node_to_pathways.get(n, [])),
        })
    rows.sort(key=lambda r: -r["total_degree"])
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows
                            else ["node"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return out_csv


# ─── Narrative report ───────────────────────────────────────────────────────

def write_report(out_dir: Path, *, sif_path: Path, G, composite,
                  pathway_panel, interactive_html, node_csv,
                  title: str) -> Path:
    counts = Counter(d.get("sign", 0) for _, _, d in G.edges(data=True))
    perturbs = [n for n, d in G.nodes(data=True)
                 if (d.get("role") or "").startswith("pert")]
    measures = [n for n, d in G.nodes(data=True)
                 if (d.get("role") or "").startswith("meas")]
    top_hubs = sorted(G.nodes, key=lambda n: -(G.in_degree(n) + G.out_degree(n)))[:5]
    lines = [
        f"# Network visualization — {title}", "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}",
        f"Source SIF: `{sif_path}`",
        f"|V| = **{G.number_of_nodes()}**  ·  |E| = **{G.number_of_edges()}** "
        f"(+1: {counts.get(1, 0)}, −1: {counts.get(-1, 0)}, "
        f"unsigned: {sum(v for k, v in counts.items() if k not in (1, -1))})",
        "",
    ]
    if composite:
        lines += ["## Publication composite figure", "",
                   f"![composite](Plots/{composite.name})", ""]
    lines += ["## Individual panels", "",
               "- `Plots/graph_force.png` — force-directed signed graph",
               "- `Plots/degree_histogram.png` — degree distribution",
               "- `Plots/edge_breakdown.png` — edge-sign counts"]
    if pathway_panel:
        lines.append("- `Plots/pathway_enrichment.png` — top enriched pathways")
    if interactive_html:
        lines.append(f"- `{interactive_html.name}` — interactive vis.js HTML "
                      f"(double-click to open in a browser)")
    lines.append(f"- `{node_csv.name}` — per-node degree / prize / role table")
    lines += [
        "",
        "## Top hubs (by total degree)",
        "",
        "| node | in | out | total |",
        "|---|---:|---:|---:|",
    ]
    for n in top_hubs:
        lines.append(f"| `{n}` | {G.in_degree(n)} | {G.out_degree(n)} | "
                      f"{G.in_degree(n) + G.out_degree(n)} |")
    if perturbs:
        lines += ["", f"**Perturbation nodes:** {', '.join(perturbs)}"]
    if measures:
        lines.append(f"**Measurement nodes:** {', '.join(measures)}")
    out = out_dir / "viz_report.md"
    out.write_text("\n".join(lines))
    return out


# ─── CLI ────────────────────────────────────────────────────────────────────

def cmd_visualize(args: argparse.Namespace) -> int:
    _setup_logging()
    sif_path = Path(args.sif)
    if not sif_path.is_file():
        raise SystemExit(f"SIF not found: {sif_path}")
    triples = read_sif(sif_path)
    if not triples:
        raise SystemExit(f"No edges parsed from {sif_path}")
    prizes = read_prizes(Path(args.prizes) if args.prizes else None)
    node_pw = read_pathways(Path(args.pathways) if args.pathways else None)

    # Tag perturbation / measurement nodes from CLI flags (if provided)
    for n in (args.perturbation or []):
        prizes.setdefault(n, {})["role"] = "perturbation"
    for n in (args.measurement or []):
        prizes.setdefault(n, {})["role"] = "measurement"

    G = build_graph(triples, prizes)
    logging.info("loaded subnetwork: |V|=%d  |E|=%d", G.number_of_nodes(),
                  G.number_of_edges())

    label = args.label or sif_path.stem
    out_dir = Path("Docs/Network") / f"{time.strftime('%Y%m%d_%H%M%S')}_{label}_viz"
    plots_dir = out_dir / "Plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    title = args.title or f"Inferred subnetwork — {label}"

    graph_png = draw_force_layout(G, plots_dir / "graph_force.png",
                                    title=title, layout=args.layout)
    selected_nodes = list(G.nodes)
    pw_panel = draw_pathway_panel(node_pw, selected_nodes,
                                    plots_dir / "pathway_enrichment.png")
    deg_panel = draw_degree_histogram(G, plots_dir / "degree_histogram.png")
    edge_panel = draw_edge_breakdown(G, plots_dir / "edge_breakdown.png")
    html_path = draw_interactive_html(G, out_dir / "graph_interactive.html",
                                        title=title) if args.html else None
    node_csv = write_node_summary(G, node_pw, out_dir / "node_summary.csv")

    composite = build_composite(
        plots_dir,
        out_path=plots_dir / "composite_publication_figure.png",
        title=title,
        subtitle=f"Source: {sif_path}  ·  |V|={G.number_of_nodes()}  |E|={G.number_of_edges()}",
    )

    report = write_report(out_dir, sif_path=sif_path, G=G,
                           composite=composite, pathway_panel=pw_panel,
                           interactive_html=html_path, node_csv=node_csv,
                           title=label)

    print(f"Report: {report}")
    print(f"Output dir: {out_dir}")
    if composite:
        print(f"Composite figure: {composite}")
        print(f"Composite (SVG): {composite.with_suffix('.svg')}")
    print(f"Graph: {graph_png}")
    if pw_panel:
        print(f"Pathway panel: {pw_panel}")
    print(f"Degree histogram: {deg_panel}")
    print(f"Edge breakdown: {edge_panel}")
    if html_path:
        print(f"Interactive HTML: {html_path}")
    print(f"Node summary CSV: {node_csv}")
    return 0


def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "NETWORK_VISUALIZATION_SKILL.md"
    path.write_text("""# Skill: Network Integration Visualization

Publication-grade visualization for `network` skill outputs (CARNIVAL /
Steiner / external SIFs).

## Inputs

- `--sif <path>` — signed-SIF triples `src \\t sign \\t dst`. Required.
- `--prizes <csv>` (optional) — per-node `prize / sign / role` table.
  Columns recognized (case-insensitive): `node|gene|symbol`, `prize|weight|score`,
  `sign|direction`, `role` (perturbation / measurement / inferred).
- `--pathways <csv>` (optional) — `node, pathway` long-form table used
  for the pathway-enrichment bar panel.
- `--perturbation NAME` / `--measurement NAME` — tag specific nodes
  (repeatable) without needing a CSV.

## Outputs

All under `Docs/Network/<ts>_<label>_viz/`:

- `Plots/graph_force.png/.svg` — force-directed signed graph
- `Plots/pathway_enrichment.png/.svg` — top-K enriched pathways
- `Plots/degree_histogram.png/.svg` — degree distribution
- `Plots/edge_breakdown.png/.svg` — edge sign counts
- `Plots/composite_publication_figure.png/.svg` — 2×2 publication mosaic
- `graph_interactive.html` — interactive vis.js HTML (with `--html`)
- `node_summary.csv` — per-node degree, sign, prize, role, pathways
- `viz_report.md` — narrative report linking to all artefacts

## Example

```bash
# Visualize the EGFR -> MYC demo cascade
igvfagent network viz \\
    --sif Docs/Network/<ts>_demo/selected_subnetwork.sif \\
    --perturbation EGFR --measurement MYC \\
    --layout spring --html --label demo
```

License: Apache-2.0. Deps: networkx (BSD), matplotlib (PSF), pyvis (MIT), pandas (BSD-3).
""")
    return path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="network_viz",
        description="Publication-grade visualization for network-"
                     "integration outputs.")
    sub = parser.add_subparsers(dest="cmd")

    p = sub.add_parser("visualize",
        help="Build force-directed graph + composite figure from a SIF.")
    p.add_argument("--sif", required=True,
                    help="Signed-SIF input (output of `network carnival/"
                          "steiner/demo`).")
    p.add_argument("--prizes",   default=None,
                    help="Optional per-node prize / role CSV.")
    p.add_argument("--pathways", default=None,
                    help="Optional node→pathway CSV for the enrichment panel.")
    p.add_argument("--perturbation", action="append", default=None,
                    help="Tag a node as perturbation (repeatable).")
    p.add_argument("--measurement", action="append", default=None,
                    help="Tag a node as measurement (repeatable).")
    p.add_argument("--layout", default="spring",
                    choices=["spring", "kamada", "circular", "shell"])
    p.add_argument("--html", action="store_true",
                    help="Also emit interactive vis.js HTML.")
    p.add_argument("--label", default=None)
    p.add_argument("--title", default=None)
    p.set_defaults(func=cmd_visualize)

    p = sub.add_parser("write-playbook",
        help="Write Docs/Skills/NETWORK_VISUALIZATION_SKILL.md")
    p.set_defaults(func=lambda a: (print(f"Wrote: {write_playbook()}"), 0)[1])

    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help(); return 2
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
