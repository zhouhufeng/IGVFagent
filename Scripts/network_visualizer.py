"""Local network-integration visualizer for the Streamlit UI.

Renders the publication-grade outputs produced by `igvfagent network viz`
(under ``Docs/Network/<ts>_<label>_viz/``):

  - Composite publication figure (4-panel mosaic)
  - Individual panel images
  - Interactive vis.js HTML (embedded with proper height)
  - Per-node summary CSV (sortable DataFrame)
  - Top hubs / perturbation / measurement nodes panel
  - SIF browser for re-generating viz on demand

Module is driver-agnostic: ``render_streamlit_panel(st)`` is the
single entry point called from ``streamlit_app.py`` under the
🔗 Network tab. The visualizer reads directly from disk — no LLM
involvement — so it's deterministic and survives backend changes.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

# Streamlit width-argument compat (see Scripts/_stcompat.py). Kept as a
# dual-mode import so the module works installed or from a checkout.
try:
    from igvfagent._stcompat import fit
except Exception:  # pragma: no cover - checkout / direct-run fallback
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent))
    from _stcompat import fit


PROJECT_ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
NETWORK_DIR = PROJECT_ROOT / "Docs" / "Network"


def discover_viz_runs() -> "list[dict[str, Any]]":
    """List every ``Docs/Network/<ts>_<label>_viz/`` directory.

    Returns a list of dicts sorted newest-first:
        {path, label, ts, has_composite, has_html, has_csv, n_nodes, n_edges}
    """
    runs: "list[dict]" = []
    if not NETWORK_DIR.is_dir():
        return runs
    for sub in NETWORK_DIR.iterdir():
        if not sub.is_dir() or not sub.name.endswith("_viz"):
            continue
        plots = sub / "Plots"
        composite = plots / "composite_publication_figure.png"
        html = sub / "graph_interactive.html"
        csv = sub / "node_summary.csv"
        report = sub / "viz_report.md"
        n_nodes = n_edges = None
        if csv.is_file():
            try:
                with csv.open() as fh:
                    n_nodes = sum(1 for _ in fh) - 1
            except Exception:
                pass
        if report.is_file():
            try:
                txt = report.read_text(errors="ignore")
                import re
                m = re.search(r"\|V\|\s*=\s*\*\*(\d+)\*\*\s*·\s*\|E\|\s*=\s*\*\*(\d+)", txt)
                if m:
                    n_nodes = int(m.group(1))
                    n_edges = int(m.group(2))
            except Exception:
                pass
        ts, _, lbl = sub.name[:15], None, sub.name[16:-4] if sub.name.endswith("_viz") else sub.name
        runs.append({
            "path": sub,
            "name": sub.name,
            "label": lbl,
            "ts": sub.name[:15],
            "composite": composite if composite.is_file() else None,
            "html": html if html.is_file() else None,
            "csv": csv if csv.is_file() else None,
            "report": report if report.is_file() else None,
            "plots_dir": plots if plots.is_dir() else None,
            "n_nodes": n_nodes,
            "n_edges": n_edges,
        })
    runs.sort(key=lambda r: r["ts"], reverse=True)
    return runs


def discover_sif_files() -> "list[Path]":
    """List every ``selected_subnetwork.sif`` (and other .sif files)
    available for visualization under ``Docs/Network/``.

    Returns newest-first list.
    """
    sifs: "list[Path]" = []
    if not NETWORK_DIR.is_dir():
        return sifs
    for f in NETWORK_DIR.rglob("*.sif"):
        sifs.append(f)
    sifs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return sifs


# ─── Streamlit panel ────────────────────────────────────────────────────────

def render_streamlit_panel(st) -> None:
    """Drop-in panel for the 🔗 Network tab."""
    st.markdown(
        "### 🔗 Network Integration Explorer\n"
        "Browse and render publication-grade visualizations of "
        "CARNIVAL / Steiner / external subnetworks. Each viz run produces "
        "a 4-panel composite figure, an interactive vis.js HTML, and a "
        "per-node summary CSV. The same outputs are produced from the "
        "CLI (`igvfagent network viz --sif <path> --html`)."
    )

    runs = discover_viz_runs()
    sifs = discover_sif_files()

    # ── Two sub-views: browse existing runs, or generate a new one ────
    if runs:
        st.markdown(f"##### Existing viz runs ({len(runs)})")
        labels = [
            f"{r['label']}  ·  {r['ts']}  "
            f"({'|V|=' + str(r['n_nodes']) if r['n_nodes'] is not None else ''}"
            f"{', |E|=' + str(r['n_edges']) if r['n_edges'] is not None else ''})"
            for r in runs
        ]
        idx = st.selectbox(
            "Pick a viz run to inspect",
            range(len(runs)),
            format_func=lambda i: labels[i],
            key="network_viz_run_pick",
        )
        run = runs[idx]
        _render_run(st, run)
    else:
        st.info(
            "No `Docs/Network/*_viz/` runs found yet. "
            "Generate one from a SIF file below, "
            "or run from the CLI:\n\n"
            "```\n"
            "igvfagent network viz --sif <path>/selected_subnetwork.sif "
            "--html --label demo\n"
            "```"
        )

    st.divider()

    # ── Generate a new viz from any SIF available on disk ─────────────
    with st.expander("➕ Generate a new visualization from a SIF file",
                     expanded=False):
        if not sifs:
            st.caption(
                "No .sif files found under `Docs/Network/`. Run "
                "`igvfagent network demo`, `network carnival`, or "
                "`network steiner` first."
            )
        else:
            sif_labels = [
                f"{s.relative_to(PROJECT_ROOT)}  "
                f"({s.stat().st_size:,} B)"
                for s in sifs
            ]
            sif_idx = st.selectbox(
                "Pick a SIF",
                range(len(sifs)),
                format_func=lambda i: sif_labels[i],
                key="network_viz_sif_pick",
            )
            label = st.text_input("Label",
                                    value=sifs[sif_idx].parent.name + "_viz_ui",
                                    key="network_viz_label")
            colp, colm = st.columns(2)
            with colp:
                perts = st.text_input(
                    "Perturbation nodes (comma-separated)", "",
                    key="network_viz_perts",
                    help="Tag specific nodes as perturbations (orange)."
                )
            with colm:
                meas = st.text_input(
                    "Measurement nodes (comma-separated)", "",
                    key="network_viz_meas",
                    help="Tag specific nodes as measurements (olive)."
                )
            layout = st.radio(
                "Layout",
                ["spring", "kamada", "circular", "shell"],
                index=0, horizontal=True, key="network_viz_layout",
            )
            do_html = st.checkbox("Also generate interactive HTML",
                                    value=True, key="network_viz_html")
            if st.button("Generate", type="primary",
                          key="network_viz_generate"):
                _run_viz(st, sifs[sif_idx], label, perts, meas, layout,
                           do_html)


def _render_run(st, run: "dict") -> None:
    """Render a single viz run's artefacts inside the panel."""
    # Header metrics
    cols = st.columns(4)
    cols[0].metric("Label", run["label"])
    cols[1].metric("|V| nodes", run["n_nodes"] if run["n_nodes"] is not None
                    else "?")
    cols[2].metric("|E| edges", run["n_edges"] if run["n_edges"] is not None
                    else "?")
    cols[3].metric("Run dir", run["ts"])

    tabs = st.tabs([
        "🖼  Composite figure", "🎯 Interactive HTML", "📊 Node summary",
        "📑 Individual panels", "📄 Report",
    ])

    # Composite figure
    with tabs[0]:
        if run["composite"]:
            try:
                st.image(str(run["composite"]),
                          caption=run["composite"].name,
                          **fit(st.image))
            except Exception as exc:
                st.error(f"Could not render composite: {exc}")
            with st.expander("Download links", expanded=False):
                _download_button(st, run["composite"], key="dl_comp_png")
                svg = run["composite"].with_suffix(".svg")
                if svg.is_file():
                    _download_button(st, svg, key="dl_comp_svg")
        else:
            st.info("No composite figure in this run.")

    # Interactive HTML (pyvis / vis.js)
    with tabs[1]:
        if run["html"]:
            try:
                html_text = run["html"].read_text(encoding="utf-8")
                st.components.v1.html(html_text, height=820, scrolling=True)
                st.caption(f"Standalone file: `{run['html']}`  ·  "
                            f"Open this `.html` in a browser for a full-window view.")
                _download_button(st, run["html"], key="dl_html",
                                  mime="text/html")
            except Exception as exc:
                st.error(f"Could not render HTML: {exc}")
        else:
            st.info(
                "No `graph_interactive.html` in this run. Generate one by "
                "re-running with `--html`."
            )

    # Node summary CSV
    with tabs[2]:
        if run["csv"]:
            try:
                import pandas as pd
                df = pd.read_csv(run["csv"])
                st.dataframe(df, **fit(st.dataframe), height=420)
                if not df.empty and "total_degree" in df.columns:
                    st.markdown("**Top 5 hubs by degree:**")
                    top = df.nlargest(5, "total_degree")
                    st.dataframe(top, **fit(st.dataframe), height=200)
                _download_button(st, run["csv"], key="dl_csv")
            except Exception as exc:
                st.error(f"Could not load CSV: {exc}")
        else:
            st.info("No `node_summary.csv` in this run.")

    # Individual panels
    with tabs[3]:
        if run["plots_dir"]:
            pngs = sorted(run["plots_dir"].glob("*.png"))
            if not pngs:
                st.info("No PNG plots in this run.")
            else:
                ncol = 2
                for i in range(0, len(pngs), ncol):
                    row = st.columns(ncol)
                    for j, png in enumerate(pngs[i:i + ncol]):
                        with row[j]:
                            try:
                                st.image(str(png), caption=png.name,
                                          **fit(st.image))
                            except Exception as exc:
                                st.warning(f"{png.name}: {exc}")
        else:
            st.info("No Plots/ directory in this run.")

    # Narrative report
    with tabs[4]:
        if run["report"]:
            try:
                st.markdown(run["report"].read_text())
            except Exception as exc:
                st.error(f"Could not read report: {exc}")
        else:
            st.info("No `viz_report.md` in this run.")


def _download_button(st, path: Path, *, key: str, mime: str | None = None):
    try:
        data = path.read_bytes()
        st.download_button(
            label=f"⬇ Download {path.name}",
            data=data, file_name=path.name,
            mime=mime,
            key=key,
        )
    except Exception as exc:
        st.caption(f"(download unavailable: {exc})")


def _run_viz(st, sif_path: Path, label: str, perts: str, meas: str,
              layout: str, do_html: bool) -> None:
    """Invoke `network_viz.cmd_visualize` in-process and refresh."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        import network_viz as nv
    except ImportError as exc:
        st.error(f"network_viz module not importable: {exc}")
        return

    class _NS:
        pass

    args = _NS()
    args.sif = str(sif_path)
    args.prizes = None
    args.pathways = None
    args.perturbation = [s.strip() for s in (perts or "").split(",")
                         if s.strip()] or None
    args.measurement = [s.strip() for s in (meas or "").split(",")
                        if s.strip()] or None
    args.layout = layout
    args.html = do_html
    args.label = label or "ui_run"
    args.title = None
    try:
        with st.spinner("Running network_viz …"):
            nv.cmd_visualize(args)
    except Exception as exc:
        st.error(f"Visualization failed: {exc}")
        return
    st.success("Done. Refreshing the run list…")
    try:
        st.rerun()
    except Exception:
        try:
            st.experimental_rerun()
        except Exception:
            pass
