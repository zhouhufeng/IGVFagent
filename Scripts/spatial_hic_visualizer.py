"""Spatial-ATAC-Hi-C browser panel for the Streamlit UI.

Renders whatever ``igvfagent spatial-hic`` has written under
``Docs/SpatialATACHiC/<ts>_<label>/``:

  - Per-pixel contact QC, against the bands Wang 2026 reports
  - The 50 x 50 tissue grid for any per-pixel value or marker gene
  - A/B compartment tracks
  - Copy-number profiles and per-pixel clone structure
  - Cluster-specific loops and the APA pileup
  - Any figure the run produced

Module is driver-agnostic: ``render_streamlit_panel(st)`` is the single
entry point called from ``streamlit_app.py`` under the 🧬 Spatial tab.
Everything is read straight off disk — no LLM involvement — so the panel
is deterministic and survives backend changes.

The grid rendering is done here rather than by shelling out to
``spatial-hic viz`` so the tab stays responsive: picking a different
gene redraws from an already-loaded matrix instead of re-running a CLI
command.
"""

from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path
from typing import Any, Optional

# Streamlit width-argument compat (see Scripts/_stcompat.py). Dual-mode
# import so the module works installed or from a checkout.
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
SPATIAL_DIR = PROJECT_ROOT / "Docs" / "SpatialATACHiC"

GRID_N = 50

# Bands Wang 2026 reports across its samples, used to colour the QC
# metrics green or amber rather than leaving the reader to look them up.
PAPER_BANDS = {
    "median_cis_fraction": (0.881, 0.903),
    "median_long_range_ratio": (0.240, 0.333),
    "median_total_contacts": (25343, 58403),
}

# Which stage produced a run, inferred from the artefacts it left.
_STAGE_MARKERS = [
    ("demux", "demux_summary.json"),
    ("qc", "qc_summary.json"),
    ("gas", "gas_summary.json"),
    ("gad", "gad_summary.json"),
    ("compartment", "compartment_summary.json"),
    ("cnv", "cnv_summary.json"),
    ("loops", "loops_summary.json"),
    ("impute", "impute_summary.json"),
    ("matrix", "matrix_summary.json"),
    ("viz", "viz_summary.json"),
    ("geo", "geo_summary.json"),
]

_PIXEL_ID_RE = re.compile(r"^(\d{1,3})x(\d{1,3})$")


# ─── Discovery ──────────────────────────────────────────────────────────────

def discover_runs() -> "list[dict[str, Any]]":
    """Every ``Docs/SpatialATACHiC/<ts>_<label>/`` directory, newest first."""
    runs: "list[dict]" = []
    if not SPATIAL_DIR.is_dir():
        return runs
    for sub in sorted(SPATIAL_DIR.iterdir(), reverse=True):
        if not sub.is_dir() or not sub.name[:2].isdigit():
            continue
        stages = [name for name, marker in _STAGE_MARKERS
                  if (sub / marker).is_file()]
        if not stages:
            continue
        # Run dirs are `YYYYMMDD_HHMMSS_<label>` — the timestamp itself
        # contains an underscore, so split on the *second* one.
        parts = sub.name.split("_", 2)
        if len(parts) == 3 and len(parts[0]) == 8 and len(parts[1]) == 6:
            ts, label = f"{parts[0]}_{parts[1]}", parts[2]
        else:
            ts, _, label = sub.name.partition("_")
        runs.append({
            "path": sub,
            "name": sub.name,
            "ts": ts,
            "label": label or sub.name,
            "stages": stages,
        })
    return runs


def _read_json(path: Path) -> "Optional[dict]":
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _read_tsv(path: Path, limit: int = 0) -> "list[dict]":
    try:
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if limit:
                return [r for _i, r in zip(range(limit), reader)]
            return list(reader)
    except Exception:
        return []


def _f(v: Any) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


# ─── Grid rendering ─────────────────────────────────────────────────────────

def _grid_from_pairs(pairs: "list[tuple[str, float]]",
                     positions: "Optional[dict[str, tuple[int, int]]]" = None):
    """Place (pixel_id, value) onto the 50 x 50 tissue grid."""
    import numpy as np
    grid = np.full((GRID_N, GRID_N), np.nan)
    placed = 0
    for name, val in pairs:
        rc = None
        if positions and name in positions:
            rc = positions[name]
        else:
            m = _PIXEL_ID_RE.match(name)
            if m:
                rc = (int(m.group(2)), int(m.group(1)))
        if rc is None:
            continue
        r, c = rc
        if 1 <= r <= GRID_N and 1 <= c <= GRID_N:
            r, c = r - 1, c - 1
        if 0 <= r < GRID_N and 0 <= c < GRID_N:
            grid[r, c] = val
            placed += 1
    return grid, placed


def _render_grid(st, grid, *, title: str, cmap: str, label: str):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        st.info("No finite values to plot.")
        return
    lo, hi = float(np.percentile(finite, 2)), float(np.percentile(finite, 98))
    if lo == hi:
        lo, hi = float(finite.min()), float(finite.max()) or lo + 1

    fig, ax = plt.subplots(figsize=(5.4, 4.8), dpi=170)
    im = ax.imshow(grid, cmap=cmap, vmin=lo, vmax=hi,
                   interpolation="nearest", origin="upper")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("barcode A channel")
    ax.set_ylabel("barcode B channel")
    ax.set_xticks([0, GRID_N // 2, GRID_N - 1])
    ax.set_xticklabels([1, GRID_N // 2 + 1, GRID_N], fontsize=8)
    ax.set_yticks([0, GRID_N // 2, GRID_N - 1])
    ax.set_yticklabels([1, GRID_N // 2 + 1, GRID_N], fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.8, label=label)
    fig.tight_layout()
    st.pyplot(fig, clear_figure=True)
    plt.close(fig)


# ─── Streamlit panel ────────────────────────────────────────────────────────

def render_streamlit_panel(st) -> None:
    """Drop-in panel for the 🧬 Spatial tab."""
    st.markdown(
        "### 🧬 Spatial-ATAC-Hi-C Explorer\n"
        "Browse spatially resolved 3D genome + chromatin accessibility runs "
        "— per-pixel contact QC, gene activity and gene-associated domain "
        "scores in tissue space, A/B compartments, copy-number clones, and "
        "cluster-specific loops. Produced by "
        "`igvfagent spatial-hic` (Wang et al., *Nat Methods* 2026; "
        "GEO GSE307620)."
    )

    runs = discover_runs()
    if not runs:
        st.info(
            "No `Docs/SpatialATACHiC/` runs found yet. Build one from the "
            "CLI — the benchmark below runs offline in about 20 seconds and "
            "produces every artefact this tab renders:\n\n"
            "```\n"
            "bash Benchmarks/wang2026_spatial_atac_hic/run.sh\n"
            "```\n\n"
            "Or start from real data:\n\n"
            "```\n"
            "igvfagent spatial-hic pull-geo --gse GSE307620 \\\n"
            "    --download 'fragments|positions'\n"
            "igvfagent spatial-hic qc --pairs-dir <dir> --label mysample\n"
            "```"
        )
        return

    labels = [f"{r['label']}  ·  {r['ts']}  ({', '.join(r['stages'])})"
              for r in runs]
    idx = st.selectbox("Pick a run", range(len(runs)),
                       format_func=lambda i: labels[i],
                       key="spatial_hic_run_pick")
    run = runs[idx]
    path = run["path"]

    tabs = st.tabs([
        "📋 QC", "🗺  Tissue map", "🧭 Compartments",
        "🧬 Copy number", "🔁 Loops", "🖼  Figures", "📁 Files",
    ])

    with tabs[0]:
        _render_qc(st, path)
    with tabs[1]:
        _render_tissue_map(st, path, runs)
    with tabs[2]:
        _render_compartments(st, path)
    with tabs[3]:
        _render_cnv(st, path)
    with tabs[4]:
        _render_loops(st, path)
    with tabs[5]:
        _render_figures(st, path)
    with tabs[6]:
        _render_files(st, path)


def _render_qc(st, path: Path) -> None:
    demux = _read_json(path / "demux_summary.json")
    if demux:
        st.markdown("##### Pixel demultiplex")
        cols = st.columns(4)
        frac = _f(demux.get("assigned_fraction"))
        cols[0].metric("Reads assigned", f"{frac:.1%}")
        cols[1].metric("Pixels with data", demux.get("pixels_with_data", "?"))
        cols[2].metric("Grid capacity", demux.get("grid_capacity", "?"))
        cols[3].metric("Median contacts",
                       f"{_f(demux.get('median_contacts_per_pixel')):,.0f}")
        if frac == frac and frac < 0.5:
            st.warning(
                f"Only {frac:.1%} of reads resolved to a pixel. The usual "
                "cause is the barcode concatenation order — re-run "
                "`pixel-demux` with `--layout AB`."
            )

    qc = _read_json(path / "qc_summary.json")
    if not qc:
        if not demux:
            st.info("No QC summary in this run.")
        return

    st.markdown("##### Per-pixel contact QC")
    st.caption("Green = inside the band Wang 2026 reports across its samples.")
    rows = []
    for key, nice in (("median_total_contacts", "Median total contacts"),
                      ("median_cis_fraction", "Median cis fraction"),
                      ("median_long_range_ratio", "Long-range ratio (≥10 kb / cis)")):
        if key not in qc:
            continue
        v = _f(qc[key])
        lo, hi = PAPER_BANDS.get(key, (None, None))
        inside = lo is not None and lo <= v <= hi
        rows.append({
            "Metric": nice,
            "Value": f"{v:,.4f}" if v < 1 else f"{v:,.0f}",
            "Paper band": f"{lo:,.3f} – {hi:,.3f}" if lo is not None and lo < 1
                          else (f"{lo:,.0f} – {hi:,.0f}" if lo is not None else "—"),
            "In band": "✅" if inside else ("⚠️" if lo is not None else "—"),
        })
    if "tss_enrichment" in qc:
        rows.append({"Metric": "TSS enrichment (ArchR)",
                     "Value": f"{_f(qc['tss_enrichment']):,.2f}",
                     "Paper band": "—", "In band": "—"})
    if rows:
        st.dataframe(rows, **fit(st.dataframe), hide_index=True)

    c = st.columns(3)
    c[0].metric("Pixels", qc.get("n_pixels", "?"))
    c[1].metric("Passing threshold", qc.get("n_pixels_passing", "?"))
    c[2].metric("Total contacts", f"{_f(qc.get('total_contacts')):,.0f}")

    tsv = path / "pixel_qc.tsv"
    if tsv.is_file():
        with st.expander("Per-pixel table", expanded=False):
            rows = _read_tsv(tsv)
            st.dataframe(rows, **fit(st.dataframe), height=340)
            _download(st, tsv, key="spatial_dl_qc")


def _render_tissue_map(st, path: Path, runs: "list[dict]") -> None:
    """Render any per-pixel value, or a gene row, onto the tissue grid."""
    # Score matrices are gene x pixel; per-pixel tables are pixel x metric.
    sources: "dict[str, Path]" = {}
    for r in runs:
        for fname, kind in (("gene_activity_score.tsv", "GAS (ATAC)"),
                            ("gene_associated_domain_score.tsv", "GAD (Hi-C)"),
                            ("pixel_qc.tsv", "per-pixel QC")):
            f = r["path"] / fname
            if f.is_file():
                sources[f"{kind} — {r['label']} · {r['ts']}"] = f
    if not sources:
        st.info(
            "No per-pixel table found in any run. Produce one with "
            "`spatial-hic gas`, `gad`, or `qc`."
        )
        return

    keys = sorted(sources)
    pick = st.selectbox("Source table", keys, key="spatial_map_source")
    table = sources[pick]

    try:
        with table.open() as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            body = list(reader)
    except Exception as exc:
        st.error(f"Could not read {table.name}: {exc}")
        return
    if not body:
        st.info(f"{table.name} has no data rows.")
        return

    is_matrix = table.name.endswith("_score.tsv")
    if is_matrix:
        genes = [r[0] for r in body if r]
        gene = st.selectbox("Gene", genes, key="spatial_map_gene")
        row = next((r for r in body if r and r[0] == gene), None)
        if row is None:
            st.info("Gene not found.")
            return
        pairs = [(header[i], _f(row[i]))
                 for i in range(1, min(len(header), len(row)))]
        title, unit = f"{gene} — {pick.split(' — ')[0]}", "score"
    else:
        numeric = [h for h in header[1:]
                   if any(_f(r[header.index(h)]) == _f(r[header.index(h)])
                          for r in body[:20] if len(r) > header.index(h))]
        if not numeric:
            st.info("No numeric column to plot.")
            return
        col = st.selectbox("Column", numeric, key="spatial_map_col")
        ci = header.index(col)
        pairs = [(r[0], _f(r[ci])) for r in body if len(r) > ci]
        title, unit = col, col

    cmap = st.selectbox("Colour map",
                        ["magma", "viridis", "inferno", "cividis", "RdBu_r"],
                        key="spatial_map_cmap")

    try:
        grid, placed = _grid_from_pairs(pairs)
    except Exception as exc:
        st.error(f"Could not build the grid: {exc}")
        return
    if placed == 0:
        st.warning(
            "No pixel could be placed on the grid. This table's keys are not "
            "`AAxBB` pixel ids — for a GAS matrix, re-run `spatial-hic gas` "
            "with `--barcode-a` / `--barcode-b` so its columns are pixel ids."
        )
        return

    _render_grid(st, grid, title=title, cmap=cmap, label=unit)
    st.caption(f"Placed {placed}/{len(pairs)} pixels on the "
               f"{GRID_N}×{GRID_N} grid.")


def _render_compartments(st, path: Path) -> None:
    summary = _read_json(path / "compartment_summary.json")
    tsv = path / "compartments.tsv"
    if not summary and not tsv.is_file():
        st.info("No compartment call in this run. Produce one with "
                "`spatial-hic compartment`.")
        return
    if summary:
        c = st.columns(4)
        c[0].metric("Resolution", f"{summary.get('resolution', 0):,} bp")
        c[1].metric("Chromosomes", summary.get("n_chroms", "?"))
        c[2].metric("Bins called",
                    f"{summary.get('n_bins_called', 0):,} / "
                    f"{summary.get('n_bins', 0):,}")
        c[3].metric("A compartment", f"{_f(summary.get('frac_A')):.1%}")
        oriented = summary.get("oriented_by", "")
        if "unoriented" in str(oriented):
            st.warning(
                "This track is **unoriented** — the eigenvector sign is "
                "arbitrary and can flip between chromosomes. Re-run with "
                "`--gene-model` to orient A toward gene-dense bins."
            )
        else:
            st.caption(f"Oriented by {oriented}.")

    if not tsv.is_file():
        return
    rows = _read_tsv(tsv)
    chroms = sorted({r["chrom"] for r in rows})
    chrom = st.selectbox("Chromosome", chroms, key="spatial_comp_chrom")
    sel = [r for r in rows if r["chrom"] == chrom and r.get("pc1") not in ("", None)]
    if not sel:
        st.info(f"No called bins on {chrom}.")
        return

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        xs = [int(r["start"]) / 1e6 for r in sel]
        ys = [_f(r["pc1"]) for r in sel]
        fig, ax = plt.subplots(figsize=(9, 2.4), dpi=170)
        ax.fill_between(xs, ys, 0, where=[y > 0 for y in ys],
                        color="#c0392b", step="mid", label="A")
        ax.fill_between(xs, ys, 0, where=[y <= 0 for y in ys],
                        color="#2874a6", step="mid", label="B")
        ax.axhline(0, color="#444", lw=0.6)
        ax.set_xlabel(f"{chrom} (Mb)")
        ax.set_ylabel("PC1")
        ax.legend(frameon=False, fontsize=8, ncol=2)
        fig.tight_layout()
        st.pyplot(fig, clear_figure=True)
        plt.close(fig)
    except Exception as exc:
        st.warning(f"Could not draw the track: {exc}")

    _download(st, tsv, key="spatial_dl_comp")


def _render_cnv(st, path: Path) -> None:
    summary = _read_json(path / "cnv_summary.json")
    pseudo = path / "cnv_pseudobulk.tsv"
    if not summary and not pseudo.is_file():
        st.info("No copy-number call in this run. Produce one with "
                "`spatial-hic cnv`.")
        return
    if summary:
        c = st.columns(4)
        c[0].metric("Resolution", f"{summary.get('resolution', 0):,} bp")
        c[1].metric("Ploidy baseline", summary.get("ploidy", "?"))
        c[2].metric("Bins called",
                    f"{summary.get('n_bins_called', 0):,} / "
                    f"{summary.get('n_bins', 0):,}")
        c[3].metric("Pixels", summary.get("n_pixels", "—"))
        flags = [k for k in ("gc_corrected", "mappability_corrected",
                             "segmented", "smoothed") if summary.get(k)]
        st.caption("Applied: " + (", ".join(flags) if flags else "no corrections") +
                   ". Bias correction is a linear fit, not NeoLoopFinder's "
                   "Poisson GLM — concordant, not bit-identical.")

    if pseudo.is_file():
        rows = _read_tsv(pseudo)
        chroms = sorted({r["chrom"] for r in rows})
        chrom = st.selectbox("Chromosome", chroms, key="spatial_cnv_chrom")
        sel = [r for r in rows if r["chrom"] == chrom
               and r.get("cn_ratio") not in ("", None)]
        if sel:
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                xs = [int(r["start"]) / 1e6 for r in sel]
                ys = [_f(r["cn_ratio"]) for r in sel]
                fig, ax = plt.subplots(figsize=(9, 2.6), dpi=170)
                ax.scatter(xs, ys, s=9, color="#333", zorder=3)
                if sel[0].get("cn_segment") not in ("", None):
                    ax.step(xs, [_f(r["cn_segment"]) for r in sel],
                            where="mid", color="#e67e22", lw=1.6,
                            label="HMM segment", zorder=2)
                    ax.legend(frameon=False, fontsize=8)
                ax.axhline(_f(summary.get("ploidy", 2)) if summary else 2,
                           color="#888", ls="--", lw=0.8)
                ax.set_xlabel(f"{chrom} (Mb)")
                ax.set_ylabel("CN ratio")
                fig.tight_layout()
                st.pyplot(fig, clear_figure=True)
                plt.close(fig)
            except Exception as exc:
                st.warning(f"Could not draw the profile: {exc}")
        _download(st, pseudo, key="spatial_dl_cnv")

    per_pixel = path / "cnv_per_pixel.tsv"
    if per_pixel.is_file():
        st.markdown("##### Per-pixel copy number (spatial clones)")
        try:
            with per_pixel.open() as fh:
                reader = csv.reader(fh, delimiter="\t")
                header = next(reader)
                body = list(reader)
            bins = header[1:]
            pick = st.selectbox("Genomic bin", bins, key="spatial_cnv_bin")
            bi = header.index(pick)
            pairs = [(r[0], _f(r[bi])) for r in body if len(r) > bi]
            grid, placed = _grid_from_pairs(pairs)
            if placed:
                _render_grid(st, grid, title=f"CN ratio — {pick}",
                             cmap="RdBu_r", label="CN ratio")
                st.caption(f"Placed {placed}/{len(pairs)} pixels. "
                           "Spatially coherent blocks are candidate clones.")
            else:
                st.info("Pixel ids in this matrix are not grid coordinates.")
        except Exception as exc:
            st.warning(f"Could not render the per-pixel matrix: {exc}")
        _download(st, per_pixel, key="spatial_dl_cnvpp")


def _render_loops(st, path: Path) -> None:
    summary = _read_json(path / "loops_summary.json")
    diff = path / "cluster_specific_loops.tsv"
    if not summary and not diff.is_file():
        st.info("No loop analysis in this run. Produce one with "
                "`spatial-hic loops --bedpe <anchors>`.")
        return
    if summary:
        c = st.columns(4)
        c[0].metric("Loops", summary.get("n_loops", "?"))
        c[1].metric("Pixels", summary.get("n_pixels", "?"))
        c[2].metric("Clusters", summary.get("n_clusters", "—"))
        c[3].metric("Cluster-specific",
                    summary.get("n_cluster_specific_loops", "—"))
        by = summary.get("specific_by_cluster") or {}
        if by:
            st.caption("Specific loops by cluster: " +
                       ", ".join(f"{k} = {v}" for k, v in sorted(by.items())))

    if diff.is_file():
        rows = _read_tsv(diff)
        only = st.checkbox("Show only cluster-specific loops", value=False,
                           key="spatial_loops_only")
        shown = [r for r in rows if r.get("specific") == "yes"] if only else rows
        st.dataframe(shown, **fit(st.dataframe), height=320, hide_index=True)
        _download(st, diff, key="spatial_dl_loops")

    apa = path / "apa_matrix.tsv"
    if apa.is_file():
        st.markdown("##### Aggregate peak analysis")
        enr = summary.get("apa_enrichment") if summary else None
        if enr is not None:
            st.caption(f"Centre / corner enrichment: **{_f(enr):.2f}×**. "
                       "A real loop set produces a hot centre pixel.")
        else:
            st.caption("Corner background is zero, so no enrichment ratio is "
                       "defined — the pileup is too sparse for a background "
                       "estimate.")
        try:
            import numpy as np
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            mat = np.loadtxt(apa, delimiter="\t")
            fig, ax = plt.subplots(figsize=(3.4, 3.0), dpi=170)
            im = ax.imshow(mat, cmap="Reds", interpolation="nearest")
            ax.set_xticks([]); ax.set_yticks([])
            fig.colorbar(im, ax=ax, shrink=0.8)
            fig.tight_layout()
            st.pyplot(fig, clear_figure=True)
            plt.close(fig)
        except Exception as exc:
            st.warning(f"Could not render the APA pileup: {exc}")


def _render_figures(st, path: Path) -> None:
    plots = path / "Plots"
    pngs = sorted(plots.glob("*.png")) if plots.is_dir() else []
    if not pngs:
        st.info("No figures in this run.")
        return
    ncol = 2
    for i in range(0, len(pngs), ncol):
        row = st.columns(ncol)
        for j, png in enumerate(pngs[i:i + ncol]):
            with row[j]:
                try:
                    st.image(str(png), caption=png.name, **fit(st.image))
                except Exception as exc:
                    st.warning(f"{png.name}: {exc}")


def _render_files(st, path: Path) -> None:
    files = sorted((p for p in path.rglob("*") if p.is_file()),
                   key=lambda p: p.name)
    if not files:
        st.info("This run directory is empty.")
        return
    rows = [{"file": str(p.relative_to(path)),
             "size": f"{p.stat().st_size:,} B"} for p in files[:400]]
    st.dataframe(rows, **fit(st.dataframe), height=360, hide_index=True)
    st.caption(f"`{path}`" + ("" if len(files) <= 400
                              else f"  ·  showing 400 of {len(files)} files"))


def _download(st, path: Path, *, key: str, mime: "Optional[str]" = None) -> None:
    try:
        st.download_button(label=f"⬇ Download {path.name}",
                           data=path.read_bytes(), file_name=path.name,
                           mime=mime, key=key)
    except Exception as exc:
        st.caption(f"(download unavailable: {exc})")
