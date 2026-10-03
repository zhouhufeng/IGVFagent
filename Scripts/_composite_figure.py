"""Helper for assembling publication-grade composite figures from a
directory of per-result PNGs. Used by the *-showcase composer skills."""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional


def build_composite_figure(
    plots_dir: Path,
    *,
    out_path: Path,
    title: str,
    subtitle: Optional[str] = None,
    layout: Iterable[tuple],  # each tuple: (panel_label, png_filename, (row, col), optional rowspan, optional colspan)
    n_rows: int,
    n_cols: int,
    panel_w: float = 5.0,
    panel_h: float = 4.2,
    title_h: float = 0.8,
    facecolor: str = "white",
) -> Optional[Path]:
    """Tile existing PNGs into a single publication-grade mosaic.

    Returns the output path on success, or None if no panels were found.
    Saves a 300-dpi PNG and a vector SVG side-by-side. Missing PNGs
    leave their panels blank with a "(not generated)" caption rather
    than failing the whole composite.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        from matplotlib import gridspec
    except ImportError:
        return None

    fig_w = n_cols * panel_w
    fig_h = n_rows * panel_h + title_h
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=facecolor)
    gs = gridspec.GridSpec(
        n_rows, n_cols, figure=fig,
        left=0.02, right=0.98, top=1 - title_h / fig_h - 0.01,
        bottom=0.02, hspace=0.18, wspace=0.06,
    )

    # Title strip at top
    fig.suptitle(title, fontsize=18, fontweight="bold", y=0.985, color="#1F2933")
    if subtitle:
        fig.text(0.5, 1 - title_h / fig_h * 0.45, subtitle,
                 ha="center", va="top", fontsize=11, color="#52606D",
                 style="italic")

    found_any = False
    for entry in layout:
        if len(entry) == 3:
            panel_label, png_name, pos = entry
            rowspan, colspan = 1, 1
        elif len(entry) == 5:
            panel_label, png_name, pos, rowspan, colspan = entry
        else:
            continue
        r, c = pos
        ax = fig.add_subplot(gs[r:r + rowspan, c:c + colspan])
        ax.axis("off")
        path = plots_dir / png_name
        if not path.is_file():
            ax.text(0.5, 0.5, f"({png_name}: not generated)",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=9, color="#A0AEC0", style="italic")
            ax.set_title(panel_label, fontsize=10.5, weight="bold",
                         color="#9D9D9D", pad=6)
            continue
        try:
            img = mpimg.imread(path)
            ax.imshow(img)
            ax.set_title(panel_label, fontsize=10.5, weight="bold",
                         color="#1F2933", pad=6)
            found_any = True
        except Exception as exc:
            ax.text(0.5, 0.5, f"(load error: {exc.__class__.__name__})",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=9, color="#C13E3E")
            ax.set_title(panel_label, fontsize=10.5, weight="bold", pad=6)
    if not found_any:
        plt.close(fig)
        return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor=facecolor)
    # Also save SVG for vector use
    svg_path = out_path.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight", facecolor=facecolor)
    plt.close(fig)
    return out_path
