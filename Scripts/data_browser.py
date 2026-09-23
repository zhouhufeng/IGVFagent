"""One file browser over every run the agent has produced.

Before this, each kind of output had its own tab with its own run picker
(Single-cell, Spatial, Network), and everything else -- eQTL enrichment,
scE2G benchmarks, Portal fetches under ``Data/Processed/<accession>`` -- had
no viewer at all outside the chat transcript that produced it. This panel
lists runs from all skills in one place and routes a file to the viewer that
understands it: an ``.h5ad`` to the embedding viewer, a Spatial-ATAC-Hi-C run
to its browser, a network viz run to the network view, and everything else
to the chat's own artefact renderer.

What counts as a run
    ``Docs/<Skill>/<YYYYMMDD_HHMMSS>_<label>/`` -- the layout every skill's
    ``run_dir()`` writes -- and ``Data/Processed/<accession>/``, where
    ``processed-first fetch`` puts files downloaded from the IGVF Portal.
    Static documentation under ``Docs/`` (Architecture, Skills, ...) is not a
    run and is not listed.

Who sees what
    Without sign-in (a laptop) and for admins: every run on disk. For a
    signed-in user on a shared deployment: only runs referenced by chat
    sessions the history store says that user may see (their own, shared
    projects). A run directory is not per-user on disk, so listing the
    whole tree to everyone would show one person's results to another.
    Every file also passes ``_pathguard.is_safe_artifact`` -- secrets under
    the workspace never render, whoever is asking.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable, Optional

try:
    from igvfagent import _pathguard  # type: ignore
except Exception:  # running from a checkout
    import _pathguard  # type: ignore

ROOT = _pathguard.project_root()
DOCS = ROOT / "Docs"
PROCESSED = ROOT / "Data" / "Processed"

_RUN_NAME = re.compile(r"^\d{8}_\d{6}")
# Never listed as a skill directory even if a timestamped folder appears.
_SKIP_SKILLS = {"Logs", "Secret", "Secrets", "Figures", "Architecture",
                "Skills", "SkillCards", "Playbooks", "References"}
MAX_FILES = 400             # per run, in the file picker
RENDER_LIMIT = 50 * 2**20   # larger files get metadata + a CLI hint only


def _is_run_dir(p: Path) -> bool:
    return p.is_dir() and bool(_RUN_NAME.match(p.name))


def _run_of(path: Path) -> Optional[Path]:
    """The run directory a file belongs to, or None."""
    try:
        rel = path.resolve().relative_to(ROOT)
    except Exception:
        return None
    parts = rel.parts
    if len(parts) >= 3 and parts[0] == "Data" and parts[1] == "Processed":
        return ROOT.joinpath(*parts[:3])
    if parts and parts[0] == "Docs":
        for j in range(2, len(parts)):
            if _RUN_NAME.match(parts[j]):
                return ROOT.joinpath(*parts[:j + 1])
    return None


def all_runs() -> "list[dict[str, Any]]":
    """Every run on disk, newest first."""
    runs: "list[dict[str, Any]]" = []
    if DOCS.is_dir():
        for skill in DOCS.iterdir():
            if not skill.is_dir() or skill.name in _SKIP_SKILLS:
                continue
            try:
                subs = list(skill.iterdir())
            except OSError:
                continue
            for sub in subs:
                if _is_run_dir(sub) and _pathguard.why_blocked(
                        sub, require_file=False) is None:
                    runs.append(_describe(sub))
    if PROCESSED.is_dir():
        for sub in PROCESSED.iterdir():
            if sub.is_dir() and _pathguard.why_blocked(
                    sub, require_file=False) is None:
                runs.append(_describe(sub))
    runs.sort(key=lambda r: r["mtime"], reverse=True)
    return runs


def runs_for_viewer(viewer: str, history) -> "list[dict[str, Any]]":
    """Runs referenced by the sessions ``viewer`` may see."""
    import json
    dirs: "dict[Path, None]" = {}
    try:
        rows = history.recent_sessions(limit=300, viewer=viewer)
    except Exception:
        return []
    for r in rows:
        try:
            full = history.session(r["run_dir"], viewer=viewer) or {}
            arts = json.loads(full.get("artefacts") or "[]")
        except Exception:
            arts = []
        for a in arts:
            p = Path(a)
            if not p.is_absolute():
                p = ROOT / p
            run = _run_of(p)
            if (run is not None and run.is_dir()
                    and _pathguard.why_blocked(run, require_file=False) is None):
                dirs[run] = None
    runs = [_describe(d) for d in dirs]
    runs.sort(key=lambda r: r["mtime"], reverse=True)
    return runs


def _describe(d: Path) -> "dict[str, Any]":
    try:
        mtime = d.stat().st_mtime
    except OSError:
        mtime = 0.0
    rel = d.relative_to(ROOT)
    group = rel.parts[1] if len(rel.parts) > 1 else ""
    if rel.parts[:2] == ("Data", "Processed"):
        group = "Portal fetch"
    return {"path": d, "rel": str(rel), "group": group, "name": d.name,
            "mtime": mtime}


def run_files(run: Path) -> "list[Path]":
    """Viewable files in a run, reports first, capped at MAX_FILES."""
    out: "list[Path]" = []
    # Unsorted walk with an early stop: a run that holds a pipeline's full
    # output tree must not be enumerated in full on every page render.
    for p in run.rglob("*"):
        if len(out) >= MAX_FILES:
            break
        if p.is_file() and _pathguard.is_safe_artifact(p):
            out.append(p)
    order = {"report.md": 0, "summary.json": 1}
    out.sort(key=lambda p: (order.get(p.name, 2),
                            0 if p.suffix in (".md", ".png", ".svg") else 1,
                            str(p)))
    return out


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def _index_of(items: "list[Any]", match: Callable[[Any], bool]) -> Optional[int]:
    for i, it in enumerate(items):
        if match(it):
            return i
    return None


def render_streamlit_panel(st, *, render_file: Callable[[str], None],
                           viewer: Optional[str] = None,
                           is_admin: bool = False,
                           history=None, scviz=None, sphic=None,
                           nwviz=None) -> None:
    """The run browser. ``render_file`` is the chat's artefact renderer."""
    st.markdown(
        "### 📂 Runs and files\n"
        "Every analysis run and every file fetched from the IGVF Portal, in "
        "one place. Pick a run, then a file; files a dedicated viewer "
        "understands open there.")
    scoped = bool(viewer) and not is_admin and history is not None
    runs = runs_for_viewer(viewer, history) if scoped else all_runs()
    if not runs:
        st.info("No runs yet. Ask something in Chat, or run any "
                "`igvfagent <skill>` command; its output directory appears "
                "here." if not scoped else
                "No runs from your sessions yet. Results of the questions "
                "you ask in Chat appear here.")
        return
    if scoped:
        st.caption(f"Showing runs from your sessions ({len(runs)}).")

    groups = sorted({r["group"] for r in runs})
    c1, c2 = st.columns([1, 2])
    with c1:
        group = st.selectbox("Skill", ["All"] + groups, key="db_group")
        text = st.text_input("Filter", key="db_filter",
                             placeholder="label or accession")
    shown = [r for r in runs
             if (group == "All" or r["group"] == group)
             and (not text or text.lower() in r["rel"].lower())]
    with c2:
        if not shown:
            st.info("No run matches.")
            return
        ri = st.selectbox(
            f"Run ({len(shown)})", range(len(shown)),
            format_func=lambda i: f"{shown[i]['group']} · {shown[i]['name']}",
            key="db_run")
    run = shown[min(ri or 0, len(shown) - 1)]
    run_path: Path = run["path"]
    st.caption(f"`{run['rel']}`")

    _route_run(st, run_path, sphic=sphic, nwviz=nwviz)

    files = run_files(run_path)
    if not files:
        st.info("This run has no viewable files.")
        return
    labels = [f"{p.relative_to(run_path)}  ·  {_human(p.stat().st_size)}"
              for p in files]
    fi = st.selectbox(f"File ({len(files)}"
                      f"{'+' if len(files) >= MAX_FILES else ''})",
                      range(len(files)), format_func=lambda i: labels[i],
                      key=f"db_file_{run['rel']}")
    path = files[min(fi or 0, len(files) - 1)]
    size = path.stat().st_size

    if path.suffix == ".h5ad":
        _route_h5ad(st, path, scviz)
        return
    if size > RENDER_LIMIT:
        st.info(f"`{path.name}` is {_human(size)}, too large to show here. "
                f"On the server it is at `{path.relative_to(ROOT)}`.")
        return
    render_file(str(path))


def _route_h5ad(st, path: Path, scviz) -> None:
    st.markdown(f"`{path.name}` is an AnnData file.")
    if scviz is None:
        st.caption("The single-cell viewer is unavailable in this build.")
        return
    try:
        paths = scviz.discover_h5ad_paths()
    except Exception:
        paths = []
    idx = _index_of(paths, lambda p: Path(p).resolve() == path.resolve())
    if idx is None:
        st.caption("The single-cell viewer does not scan this directory, so "
                   "it cannot open this file.")
        return

    def _open() -> None:
        st.session_state["sc_viz_dataset"] = idx

    st.button("Open in the Single-cell viewer", on_click=_open,
              key=f"db_open_sc_{path}")
    st.caption("Then switch to the **Single-cell** sub-tab.")


def _route_run(st, run_path: Path, *, sphic, nwviz) -> None:
    """Offer the dedicated viewer for runs that have one."""
    rel = run_path.relative_to(ROOT).parts
    if len(rel) < 2:
        return
    if rel[1] == "SpatialATACHiC" and sphic is not None:
        try:
            runs = sphic.discover_runs()
        except Exception:
            runs = []
        idx = _index_of(runs, lambda r: Path(r["path"]).resolve()
                        == run_path.resolve())
        if idx is not None:
            def _open() -> None:
                st.session_state["spatial_hic_run_pick"] = idx
            st.button("Open in the Spatial-ATAC-Hi-C viewer", on_click=_open,
                      key=f"db_open_sp_{run_path}")
            st.caption("Then switch to the **Spatial-ATAC-Hi-C** sub-tab.")
    if rel[1] == "Network" and nwviz is not None:
        try:
            runs = nwviz.discover_viz_runs()
        except Exception:
            runs = []
        idx = _index_of(runs, lambda r: Path(r["path"]).resolve()
                        == run_path.resolve())
        if idx is not None:
            def _open() -> None:
                st.session_state["network_viz_run_pick"] = idx
            st.button("Open in Your networks", on_click=_open,
                      key=f"db_open_nw_{run_path}")
            st.caption("Then open **🕸 Knowledge & networks → Your networks**.")
