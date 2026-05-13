"""Single-cell embedding visualizer — Streamlit panel.

Surfaces the `.h5ad` outputs from ``sc-analyze`` inside the IGVFagent
UI: pick a dataset, choose UMAP or t-SNE, colour by cluster / sample /
QC metric / gene expression. Matches the pattern used by
``kg_visualizer.py``.

The module is intentionally Streamlit-aware on its public surface
(``render_streamlit_panel``) but the underlying load + render
functions are reusable from a notebook or the CLI.
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger("sc_visualizer")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCAN_DIRS = [
    PROJECT_ROOT / "Docs" / "SingleCell",
    PROJECT_ROOT / "Data" / "SingleCell",
    PROJECT_ROOT / "Data",
]


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


@dataclass
class H5adDescriptor:
    path: Path
    label: str
    n_obs: int = 0
    n_vars: int = 0
    obs_cols: "list[str]" = field(default_factory=list)
    embeddings: "list[str]" = field(default_factory=list)
    has_leiden: bool = False
    has_markers: bool = False
    size_mb: float = 0.0
    error: Optional[str] = None


def _inspect_h5ad(path: Path) -> H5adDescriptor:
    """Lightweight peek at an h5ad without loading the whole matrix."""
    desc = H5adDescriptor(path=path, label=path.parent.name)
    try:
        desc.size_mb = path.stat().st_size / (1 << 20)
        # anndata supports a backed-mode read which avoids loading X.
        import anndata as ad  # type: ignore
        a = ad.read_h5ad(path, backed="r")
        desc.n_obs = a.n_obs
        desc.n_vars = a.n_vars
        desc.obs_cols = list(a.obs.columns)
        desc.embeddings = [k for k in a.obsm.keys()
                            if k.startswith("X_")]
        desc.has_leiden = "leiden" in a.obs.columns
        desc.has_markers = (a.uns is not None
                              and "rank_genes_groups" in a.uns)
        try:
            a.file.close()
        except Exception:
            pass
    except Exception as e:
        desc.error = str(e)
    return desc


def discover_h5ad_files() -> "list[H5adDescriptor]":
    """Walk the conventional output dirs + a few extras for .h5ad files."""
    found: "list[Path]" = []
    seen: "set[Path]" = set()
    for d in SCAN_DIRS:
        if not d.exists():
            continue
        for p in d.rglob("*.h5ad"):
            if p in seen:
                continue
            seen.add(p)
            found.append(p)
    # Sort newest first by mtime
    found.sort(key=lambda p: p.stat().st_mtime if p.exists() else 0,
                reverse=True)
    return [_inspect_h5ad(p) for p in found]


# ---------------------------------------------------------------------------
# Plot rendering — matplotlib + scanpy
# ---------------------------------------------------------------------------


def _scanpy():
    try:
        import matplotlib  # type: ignore
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        import scanpy as sc  # type: ignore
        import anndata as ad  # type: ignore
        import numpy as np  # type: ignore
        import pandas as pd  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "single-cell visualization requires scanpy + anndata + "
            "matplotlib (+ optional networkx + leidenalg). Install:\n"
            "  pip install scanpy 'anndata>=0.10' umap-learn leidenalg "
            "python-igraph matplotlib\n"
            f"Original error: {exc}"
        )
    return sc, ad, np, pd, plt


def load_h5ad(path: Path):
    sc, ad, np, pd, plt = _scanpy()
    return ad.read_h5ad(str(path))


def render_embedding_figure(adata, *, embedding: str = "umap",
                              color_keys: "list[str]",
                              ncols: int = 2,
                              point_size: int = 12) -> Any:
    """Render UMAP / t-SNE / PCA coloured by the chosen obs cols / genes."""
    sc, ad, np, pd, plt = _scanpy()
    key = f"X_{embedding.lower()}"
    if key not in adata.obsm:
        # Build the figure with a friendly message
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5,
                 f"No `{key}` embedding in this h5ad.\n"
                 "Run `sc-analyze umap` / `tsne` first.",
                 ha="center", va="center", fontsize=11, alpha=0.75)
        ax.axis("off")
        return fig

    valid = []
    for c in color_keys:
        if c in adata.obs.columns:
            valid.append(c)
        elif c in adata.var_names:
            valid.append(c)
        elif adata.raw is not None and c in adata.raw.var_names:
            valid.append(c)
    if not valid:
        valid = ["leiden"] if "leiden" in adata.obs.columns else [
            adata.obs.columns[0]]

    # scanpy plotting API: `sc.pl.umap(adata, color=...)` etc.
    plotter = getattr(sc.pl, embedding.lower())
    try:
        plotter(
            adata, color=valid, ncols=min(ncols, len(valid)),
            show=False, use_raw=True, frameon=False, size=point_size,
        )
    except TypeError:
        # Older scanpy without `size` kw
        plotter(adata, color=valid, ncols=min(ncols, len(valid)),
                 show=False, use_raw=True, frameon=False)
    fig = plt.gcf()
    fig.tight_layout()
    return fig


def render_qc_figure(adata) -> Any:
    sc, ad, np, pd, plt = _scanpy()
    cols = [c for c in ("n_genes_by_counts", "total_counts",
                          "pct_counts_mt")
             if c in adata.obs.columns]
    if not cols:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No QC metrics in this h5ad.\n"
                          "Run `sc-analyze qc` first.",
                 ha="center", va="center", alpha=0.7)
        ax.axis("off")
        return fig
    fig, axes = plt.subplots(1, len(cols), figsize=(3.5 * len(cols), 3.2))
    if len(cols) == 1:
        axes = [axes]
    colors = ["#1F77B4", "#2CA02C", "#D62728"]
    for ax, col, color in zip(axes, cols, colors):
        vals = adata.obs[col].values
        ax.violinplot([vals], showmedians=True, widths=0.7)
        ax.scatter(np.ones_like(vals)
                    + np.random.uniform(-0.15, 0.15, len(vals)),
                    vals, s=2, alpha=0.25, color=color)
        ax.set_title(col); ax.set_xticks([])
        ax.grid(alpha=0.25, linestyle=":")
    fig.suptitle(f"QC violins  ·  n_cells = {adata.n_obs:,}",
                  fontsize=11)
    fig.tight_layout()
    return fig


def cluster_composition(adata, group_key: str = "leiden",
                          breakdown_by: Optional[str] = None) -> dict:
    sc, ad, np, pd, plt = _scanpy()
    if group_key not in adata.obs.columns:
        return {}
    if breakdown_by and breakdown_by in adata.obs.columns:
        ct = pd.crosstab(adata.obs[group_key], adata.obs[breakdown_by])
        return {"crosstab": ct.to_dict(),
                 "cluster_sizes": adata.obs[group_key].value_counts().to_dict()}
    return {"cluster_sizes": adata.obs[group_key].value_counts().to_dict()}


def top_markers_df(adata, n_top: int = 20):
    """Extract top-N markers per cluster from anndata.uns."""
    sc, ad, np, pd, plt = _scanpy()
    if "rank_genes_groups" not in adata.uns:
        return None
    res = adata.uns["rank_genes_groups"]
    groups = res["names"].dtype.names
    rows = []
    for g in groups:
        for i in range(min(n_top, len(res["names"][g]))):
            rows.append({
                "cluster": g,
                "rank": i + 1,
                "gene": str(res["names"][g][i]),
                "logfc": float(res["logfoldchanges"][g][i]),
                "padj": float(res["pvals_adj"][g][i]),
                "score": float(res["scores"][g][i]),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Streamlit panel
# ---------------------------------------------------------------------------


def _ensure_in_session(st, key: str, default):
    if key not in st.session_state:
        st.session_state[key] = default


def render_streamlit_panel(st) -> None:
    """Drop-in Streamlit panel for single-cell visualization."""
    st.markdown(
        "### 🔬 Single-cell embedding viewer\n"
        "Browse the `.h5ad` outputs that `sc-analyze` writes under "
        "`Docs/SingleCell/<run>/processed.h5ad`. Pick a dataset, choose "
        "an embedding (UMAP or t-SNE), and colour by Leiden cluster, "
        "any obs column, or gene expression."
    )

    descs = discover_h5ad_files()
    if not descs:
        st.warning(
            "No `.h5ad` files found. Build one with:\n\n"
            "```\n"
            "igvfagent sc-analyze pipeline --input <counts> --label demo\n"
            "```\n"
            "Or for an instant demo, generate Scanpy's PBMC3k tutorial:\n"
            "```\n"
            "python -c 'import scanpy as sc; "
            "sc.datasets.pbmc3k().write_h5ad(\"/tmp/pbmc3k.h5ad\")'\n"
            "igvfagent sc-analyze pipeline --input /tmp/pbmc3k.h5ad "
            "--label pbmc3k_demo\n"
            "```\n"
            "Then refresh this tab."
        )
        return

    # File selector
    labels = [f"{d.label}  ·  {d.n_obs:,} cells × {d.n_vars:,} genes  "
              f"·  {d.size_mb:.1f} MB" for d in descs]
    idx = st.selectbox(
        "Dataset",
        list(range(len(descs))),
        format_func=lambda i: labels[i],
        key="sc_viz_dataset",
    )
    desc = descs[idx]
    st.caption(f"`{desc.path}`")

    if desc.error:
        st.error(f"Could not open: {desc.error}")
        return

    # Feature snapshot strip
    cols = st.columns(4)
    cols[0].metric("Cells", f"{desc.n_obs:,}")
    cols[1].metric("Genes", f"{desc.n_vars:,}")
    embs = [e.replace("X_", "") for e in desc.embeddings]
    cols[2].metric("Embeddings", ", ".join(embs) or "—")
    cols[3].metric(
        "Has clusters / markers",
        ("✅" if desc.has_leiden else "—")
        + " / "
        + ("✅" if desc.has_markers else "—"),
    )

    if not desc.embeddings:
        st.info(
            "This dataset has no UMAP/t-SNE embedding yet. Run:\n\n"
            f"```\nigvfagent sc-analyze umap --input {desc.path}\n```"
        )
        return

    # Load the full anndata once and cache in session for fast re-render
    cache_key = f"_sc_adata_{desc.path}"
    if cache_key not in st.session_state:
        with st.spinner(f"Loading {desc.path.name}…"):
            try:
                st.session_state[cache_key] = load_h5ad(desc.path)
            except Exception as e:
                st.error(f"Load failed: {e}")
                return
    adata = st.session_state[cache_key]

    # Embedding + coloring picker
    st.divider()
    pick_col, pick_size, pick_ncols = st.columns([1, 1, 1])
    available_embs = [e.replace("X_", "") for e in desc.embeddings
                       if e in ("X_umap", "X_tsne", "X_pca")]
    with pick_col:
        embedding = st.radio(
            "Embedding",
            options=available_embs,
            horizontal=True,
            key="sc_viz_embedding",
        )
    with pick_size:
        point_size = st.slider("Point size", 4, 30, 12,
                                 key="sc_viz_point_size")
    with pick_ncols:
        ncols = st.slider("Columns per panel", 1, 4, 2,
                            key="sc_viz_ncols")

    # Color picker: obs columns + gene search
    obs_choices = list(adata.obs.columns)
    default_color = (["leiden"] if "leiden" in obs_choices
                      else [obs_choices[0]] if obs_choices else [])
    color_keys = st.multiselect(
        "Colour by (obs columns)",
        options=obs_choices,
        default=default_color,
        key="sc_viz_color_obs",
    )
    gene_input = st.text_input(
        "Add genes to overlay (comma-separated symbols)",
        placeholder="e.g. CD3D, CD8A, MS4A1, LYZ, APOE",
        key="sc_viz_color_genes",
    )
    gene_keys: "list[str]" = []
    if gene_input.strip():
        candidates = [g.strip() for g in gene_input.split(",") if g.strip()]
        # Validate against var_names (use .raw if HVGs were filtered)
        var_set = set(adata.var_names)
        raw_set = set(adata.raw.var_names) if adata.raw is not None else set()
        missing = []
        for g in candidates:
            if g in var_set or g in raw_set:
                gene_keys.append(g)
            else:
                missing.append(g)
        if missing:
            st.caption(f"⚠️ Not found in dataset: {', '.join(missing)}")

    final_color = (color_keys or []) + gene_keys
    if not final_color:
        final_color = ["leiden"] if "leiden" in obs_choices else (
            [obs_choices[0]] if obs_choices else [])

    # Render
    with st.spinner(f"Rendering {embedding.upper()} ({len(final_color)} panels)…"):
        try:
            fig = render_embedding_figure(
                adata, embedding=embedding,
                color_keys=final_color, ncols=ncols,
                point_size=point_size,
            )
            st.pyplot(fig, use_container_width=True)
        except Exception as exc:
            import traceback
            st.error(f"Render error: {exc}")
            with st.expander("Traceback"):
                st.code(traceback.format_exc())

    # QC + markers + cluster composition (expandable)
    with st.expander("QC violins", expanded=False):
        try:
            st.pyplot(render_qc_figure(adata), use_container_width=True)
        except Exception as exc:
            st.warning(f"QC plot unavailable: {exc}")

    if desc.has_markers:
        with st.expander("Top markers per cluster", expanded=False):
            n_top = st.slider("Top markers per cluster", 5, 50, 15,
                                key="sc_viz_n_markers")
            df = top_markers_df(adata, n_top=n_top)
            if df is not None and len(df):
                # Pivoted view — one column per cluster, top-N rows
                pivot = (df.pivot(index="rank", columns="cluster",
                                    values="gene"))
                st.dataframe(pivot, use_container_width=True, height=360)
                st.caption("Full ranked table:")
                st.dataframe(df, use_container_width=True, height=240)

    if desc.has_leiden:
        with st.expander("Cluster composition", expanded=False):
            obs_cols = list(adata.obs.columns)
            split_by = st.selectbox(
                "Optional breakdown",
                ["(none)"] + [c for c in obs_cols
                                if c != "leiden" and adata.obs[c].nunique() <= 50],
                key="sc_viz_split_by",
            )
            sb = split_by if split_by != "(none)" else None
            comp = cluster_composition(adata, group_key="leiden",
                                          breakdown_by=sb)
            if comp.get("cluster_sizes"):
                sizes = {str(k): v for k, v in comp["cluster_sizes"].items()}
                st.bar_chart(sizes)
            if sb and comp.get("crosstab"):
                st.dataframe(comp["crosstab"], use_container_width=True,
                              height=240)

    with st.expander("Obs columns table (cell metadata sample)",
                      expanded=False):
        st.dataframe(adata.obs.head(50), use_container_width=True,
                      height=320)


# ---------------------------------------------------------------------------
# CLI smoke entrypoint (debug only)
# ---------------------------------------------------------------------------


def _cli_smoke() -> int:
    for d in discover_h5ad_files():
        line = (f"  {d.path}  ({d.n_obs:,}×{d.n_vars:,}  "
                f"embs={d.embeddings}  leiden={d.has_leiden}  "
                f"markers={d.has_markers})")
        print(line)
        if d.error:
            print(f"    ERROR: {d.error}")
    return 0


if __name__ == "__main__":
    sys.exit(_cli_smoke())
