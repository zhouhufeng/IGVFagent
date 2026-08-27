#!/usr/bin/env python3
"""Benchmark figure gallery for the IGVFagent Streamlit UI.

Renders every figure the reproducibility benchmark suite has produced —
the suite-level dashboard, each paper's plots, and the latest concordance
scoring run — directly in the browser, so the results are inspectable
without digging through ``Benchmarks/<paper-id>/figures/`` on disk.

Exposes ``render_streamlit_panel(st)``, matching the contract of
``kg_visualizer`` / ``sc_visualizer`` / ``network_visualizer``.

Design notes
------------
* **PNG by default, SVG on request.** Both are committed for every figure.
  PNG renders through ``st.image``; SVG needs the iframe embed that
  ``streamlit_app._render_svg`` provides (Streamlit >=1.50 dropped SVG
  support from ``st.image``), so we reuse that helper when available and
  fall back to a download button when it isn't.
* **Benchmarks with no figures are reported, not hidden.** A paper that
  produced no plots is a real gap in the suite; silently omitting it would
  make the gallery look complete when it isn't.
* **Concordance results are optional.** ``Benchmarks/results/`` is
  gitignored and regenerated per scoring run, so a fresh checkout has
  none. We say so and name the command instead of erroring.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# Streamlit width-argument compat (see Scripts/_stcompat.py) — same shim the
# other visualizers use, so this panel renders on old and new Streamlit alike.
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from igvfagent._stcompat import fit  # type: ignore
except Exception:
    try:
        from _stcompat import fit  # type: ignore
    except Exception:                       # last resort: no width kwarg
        def fit(fn, stretch: bool = True) -> "Dict[str, Any]":  # type: ignore
            return {}

__all__ = ["render_streamlit_panel", "find_benchmarks", "benchmarks_root"]

# Figures are committed as matched .png/.svg pairs; these are the raster
# formats st.image can load directly.
_RASTER_SUFFIXES = (".png", ".jpg", ".jpeg", ".gif", ".webp")

# Status icons mirror Benchmarks/concordance.py so the gallery and the
# terminal report read the same way.
_STATUS_ICON = {
    "ok": "✅", "partial": "🟡", "unreviewed": "⊘",
    "no_run_found": "⚪", "failed": "❌", "error": "❌",
}


def benchmarks_root() -> Optional[Path]:
    """Locate the repo's ``Benchmarks/`` directory, or None."""
    import os
    env = os.environ.get("IGVF_PROJECT_ROOT")
    candidates = []
    if env:
        candidates.append(Path(env) / "Benchmarks")
    # Scripts/benchmark_visualizer.py -> repo root is one level up.
    candidates.append(Path(__file__).resolve().parent.parent / "Benchmarks")
    candidates.append(Path.cwd() / "Benchmarks")
    for c in candidates:
        if c.is_dir():
            return c
    return None


def _headline(readme: Path) -> str:
    """First bolded claim or first prose paragraph of a benchmark README."""
    try:
        text = readme.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return ""
    # Prefer an explicit **Question.** / **Result.** lead-in.
    m = re.search(r"^\*\*(?:Question|Result|Headline)\.?\*\*\s*(.+?)(?:\n\n|\Z)",
                  text, re.M | re.S)
    if m:
        return " ".join(m.group(1).split())[:400]
    for para in text.split("\n\n"):
        p = para.strip()
        if p and not p.startswith(("#", "![", "|", "```")):
            return " ".join(p.split())[:400]
    return ""


def find_benchmarks(root: Path) -> List[Dict[str, Any]]:
    """Inventory every benchmark directory and the figures it produced."""
    out: List[Dict[str, Any]] = []
    for d in sorted(root.iterdir()):
        if not d.is_dir() or d.name in ("figures", "results", "__pycache__"):
            continue
        figdir = d / "figures"
        rasters = sorted(
            p for p in figdir.glob("*")
            if p.suffix.lower() in _RASTER_SUFFIXES
        ) if figdir.is_dir() else []
        svgs = sorted(figdir.glob("*.svg")) if figdir.is_dir() else []
        readme = d / "README.md"
        out.append({
            "id": d.name,
            "dir": d,
            "rasters": rasters,
            "svgs": svgs,
            "n_fig": len(rasters) or len(svgs),
            "readme": readme if readme.is_file() else None,
            "headline": _headline(readme) if readme.is_file() else "",
            "expected": (d / "expected.json") if (d / "expected.json").is_file() else None,
        })
    return out


def _latest_concordance(root: Path) -> Optional[Path]:
    """Newest ``results/<ts>_concordance.json``, if any scoring run exists."""
    rd = root / "results"
    if not rd.is_dir():
        return None
    runs = sorted(rd.glob("*_concordance.json"))
    return runs[-1] if runs else None


def _render_figure(st, path: Path, *, svg_mode: bool) -> None:
    """Render one figure, reusing the app's SVG embed helper when present."""
    if svg_mode and path.suffix.lower() == ".svg":
        render_svg = None
        try:                                    # same-package import
            from igvfagent.streamlit_app import _render_svg as render_svg  # type: ignore
        except Exception:
            try:                                # bare-checkout import
                from streamlit_app import _render_svg as render_svg  # type: ignore
            except Exception:
                render_svg = None
        if render_svg is not None:
            render_svg(str(path))
            return
        # No helper available: offer the file rather than a broken <img>.
        st.caption(f"`{path.name}` — SVG inline rendering unavailable")
        try:
            st.download_button("⬇ " + path.name, path.read_bytes(),
                                file_name=path.name, mime="image/svg+xml",
                                key=f"dl_{path}")
        except Exception:
            pass
        return
    try:
        st.image(str(path), caption=path.name, **fit(st.image))
    except Exception as e:
        st.caption(f"(could not render `{path.name}`: {e})")


def _render_checks(st, expected: Path) -> None:
    """Show a benchmark's declared ground-truth checks."""
    try:
        spec = json.loads(expected.read_text(encoding="utf-8"))
    except Exception as e:
        st.caption(f"(expected.json unreadable: {e})")
        return
    checks = spec.get("checks") or []
    if not checks:
        st.caption("_No declared checks._")
        return
    n_unconf = sum(1 for c in checks if c.get("confirmed") is False)
    rows = []
    for c in checks:
        # An explicit "confirmed": false marks a claim scraped from paper
        # prose that concordance.py reports but never scores.
        unconf = c.get("confirmed") is False
        rows.append({
            "": "⊘" if unconf else "•",
            "check": c.get("name", ""),
            "type": c.get("type", ""),
            "scored": "no" if unconf else "yes",
        })
    st.dataframe(rows, **fit(st.dataframe), hide_index=True)
    if n_unconf:
        st.caption(
            f"⊘ {n_unconf} of {len(checks)} check(s) are unconfirmed "
            "(extracted from paper prose) and are reported but **never "
            "scored** — see Benchmarks/OPERATIONS_GUIDE.md."
        )


def render_streamlit_panel(st) -> None:  # noqa: C901
    """Render the benchmark figure gallery."""
    root = benchmarks_root()
    if root is None:
        st.warning(
            "Could not locate the `Benchmarks/` directory. Set "
            "`IGVF_PROJECT_ROOT` to the repo root and restart the UI."
        )
        return

    benches = find_benchmarks(root)
    with_figs = [b for b in benches if b["n_fig"]]
    without = [b for b in benches if not b["n_fig"]]
    total_figs = sum(b["n_fig"] for b in with_figs)

    st.subheader("📊 Reproducibility benchmark figures")
    c1, c2, c3 = st.columns(3)
    c1.metric("Benchmarks", len(benches))
    c2.metric("With figures", len(with_figs))
    c3.metric("Figures", total_figs)

    fmt = st.radio(
        "Format", ["PNG (raster)", "SVG (vector)"], horizontal=True,
        key="_bench_fmt",
        help="Both are committed for every figure. SVG is vector "
             "(publication-ready) but renders in an iframe, so it is slower.",
    )
    svg_mode = fmt.startswith("SVG")

    # ---------------------------- suite dashboard ----------------------
    dash = root / "figures" / ("dashboard.svg" if svg_mode else "dashboard.png")
    if dash.is_file():
        st.markdown("### Suite dashboard")
        _render_figure(st, dash, svg_mode=svg_mode)

    # -------------------------- concordance report ---------------------
    st.markdown("### Latest concordance run")
    latest = _latest_concordance(root)
    if latest is None:
        st.info(
            "No scoring run found. `Benchmarks/results/` is gitignored and "
            "regenerated on each run — produce one with:\n\n"
            "```\npython Benchmarks/concordance.py --all\n```"
        )
    else:
        try:
            rep = json.loads(latest.read_text(encoding="utf-8"))
            per = rep.get("benchmarks") or rep.get("results") or []
            rows = []
            for r in per:
                status = str(r.get("status", ""))
                rows.append({
                    "": _STATUS_ICON.get(status, "•"),
                    "benchmark": r.get("benchmark") or r.get("paper_id", ""),
                    "status": status,
                    "passed": r.get("n_passed", ""),
                    "total": r.get("n_total", ""),
                })
            st.caption(f"`{latest.name}`")
            if rows:
                st.dataframe(rows, **fit(st.dataframe), hide_index=True)
            else:
                st.caption("_(report contained no per-benchmark rows)_")
        except Exception as e:
            st.caption(f"(could not read {latest.name}: {e})")

    # ---------------------------- per-paper figures --------------------
    st.markdown("### Per-paper figures")
    ids = [b["id"] for b in with_figs]
    picked = st.multiselect(
        "Papers to display", ids, default=ids, key="_bench_papers",
        help="All papers are shown by default. Narrow the list to speed up "
             "rendering or to focus on one reproduction.",
    )
    expand = st.checkbox("Expand all", value=True, key="_bench_expand")

    for b in with_figs:
        if b["id"] not in picked:
            continue
        figs = b["svgs"] if (svg_mode and b["svgs"]) else b["rasters"]
        if svg_mode and not b["svgs"]:
            figs = b["rasters"]        # no vector version committed
        with st.expander(f"**{b['id']}** — {len(figs)} figure(s)",
                          expanded=expand):
            if b["headline"]:
                st.caption(b["headline"])
            for p in figs:
                _render_figure(st, p, svg_mode=svg_mode)
            if b["expected"] is not None:
                st.markdown("**Declared checks**")
                _render_checks(st, b["expected"])
            if b["readme"] is not None:
                with st.expander("README", expanded=False):
                    try:
                        st.markdown(b["readme"].read_text(
                            encoding="utf-8", errors="replace"))
                    except Exception as e:
                        st.caption(f"(README unreadable: {e})")

    # --------------------------- honest coverage gap -------------------
    if without:
        st.markdown("### Benchmarks with no committed figures")
        st.caption(
            "These reproductions are part of the suite but produced no plots, "
            "so nothing is shown above for them — a real coverage gap, not a "
            "rendering failure."
        )
        for b in without:
            st.markdown(f"- `{b['id']}`"
                        + (f" — {b['headline'][:160]}" if b["headline"] else ""))


if __name__ == "__main__":
    r = benchmarks_root()
    if r is None:
        raise SystemExit("Benchmarks/ not found")
    bs = find_benchmarks(r)
    print(f"{len(bs)} benchmarks, {sum(b['n_fig'] for b in bs)} figures")
    for b in bs:
        print(f"  {b['id']:34} {b['n_fig']:2} fig")
