"""Generate a publication-quality flow chart of the IGVFagent architecture (v2)."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import matplotlib.patheffects as path_effects

# ─── Palette ────────────────────────────────────────────────────────────────
BG = "#FAF7F2"
INK = "#1F2933"; SUBINK = "#52606D"
EDGE = "#3D405B"; EDGE_LIGHT = "#A0AEC0"

L5 = "#9D6B9D"   # plum   — user entry
L4 = "#E89B3C"   # amber  — agent runtime
L2A = "#7A9070"  # sage   — filesystem
L2B = "#345267"  # slate  — duckdb
L1 = "#4A5A6E"   # gray-slate — upstream

SKILL_GROUPS = [
    ("Discovery & Retrieval", "#81B29A",
     "client · data · frontpage · specialized\nexplain · perturb-catalog · geo · ref"),
    ("Variant Analysis", "#E07A5F",
     "variant · advanced-variant · ccre"),
    ("Regulatory Genomics", "#F2A359",
     "enhancer · mpra · starrseq\ncrispri · flowfish · se-targets"),
    ("Single-cell & Multiomics", "#5C8DAA",
     "singlecell · sc-analyze · multiome\nsplitseq · share · multiseq"),
    ("Bulk & Proteomics", "#7C73A5",
     "rnaseq · encode · proteomics"),
    ("Knowledge Integration", "#C13E3E",
     "kg · portal-kg · kg-mirror\nwarehouse · ★ network ★"),
]

# Figure
fig, ax = plt.subplots(figsize=(16, 13.5))
ax.set_xlim(0, 100); ax.set_ylim(0, 100)
ax.set_aspect("equal"); ax.axis("off")
fig.patch.set_facecolor(BG); ax.set_facecolor(BG)


def text(x, y, s, *, size=10, weight="normal", color=INK, ha="center", va="center",
         style="normal"):
    return ax.text(x, y, s, fontsize=size, fontweight=weight, color=color,
                   ha=ha, va=va, family="DejaVu Sans", style=style, zorder=10)


def rounded_box(x, y, w, h, *, fc, ec=EDGE, lw=1.2, alpha=1.0, rounding=1.0,
                shadow=True):
    if shadow:
        ax.add_patch(FancyBboxPatch(
            (x + 0.25, y - 0.35), w, h,
            boxstyle=f"round,pad=0.0,rounding_size={rounding}",
            linewidth=0, fc="#000000", alpha=0.07, zorder=1))
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0.0,rounding_size={rounding}",
        linewidth=lw, edgecolor=ec, facecolor=fc, alpha=alpha, zorder=2))


def arrow(x1, y1, x2, y2, *, color=EDGE, lw=2.0, head=15, alpha=1.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                  mutation_scale=head, color=color,
                                  linewidth=lw, alpha=alpha, zorder=5,
                                  connectionstyle="arc3,rad=0"))


def layer_header(y_top, label_main, label_sub, fill, text_color):
    """A thin colored ribbon ABOVE a layer's content."""
    rounded_box(1.0, y_top - 2.5, 98.0, 2.4, fc=fill, ec=fill, lw=0,
                rounding=0.7, shadow=False)
    text(2.5, y_top - 1.3, label_main, size=10.5, weight="bold",
         color=text_color, ha="left")
    text(97.5, y_top - 1.3, label_sub, size=8.8, color=text_color,
         ha="right", style="italic")


# ─── Title ──────────────────────────────────────────────────────────────────
text(50, 97.5, "IGVFagent  ·  Architecture and Skill Topology",
     size=24, weight="bold", color=INK)
text(50, 95.0,
     "30 skills  ·  98 registered tools  ·  Apache-2.0  ·  clean-room "
     "reimplementations only",
     size=12, color=SUBINK, style="italic")
ax.add_patch(Rectangle((10, 93.6), 80, 0.18, color=EDGE_LIGHT, zorder=1))


# ─── LAYER 5 — User entry points ─────────────────────────────────────────────
layer_header(92.5, "Layer 5 · User entry points",
             "three faces of the same engine",
             L5, "#FFFFFF")
y5_top = 90; y5_bot = 84
for i, (label, sub) in enumerate([
    ("Terminal CLI", "igvfagent <skill> <subcmd>\ndeterministic · scriptable"),
    ("Natural-language agent", "igvfagent ask \"<query>\"\nLLM plans · chains tools"),
    ("Browser UI", "igvfagent ui\nStreamlit chat + live events"),
]):
    x = 5 + i * 31
    w = 28
    rounded_box(x, y5_bot, w, y5_top - y5_bot, fc="#FFFFFF", ec=L5, lw=1.8)
    text(x + w / 2, y5_top - 1.3, label, size=12, weight="bold", color=L5)
    text(x + w / 2, y5_top - 3.7, sub, size=8.8, color=SUBINK)

# Convergence arrows ↓
for i in range(3):
    x = 5 + i * 31 + 14
    arrow(x, y5_bot - 0.3, 50, 79.5, color=EDGE, lw=1.4, head=11, alpha=0.6)


# ─── LAYER 4 — Agent runtime ─────────────────────────────────────────────────
layer_header(82.0, "Layer 4 · Agent runtime & tool dispatch",
             "every skill is a Tool here", L4, "#1F2933")
y4_top = 79.5; y4_bot = 73
for i, (label, sub) in enumerate([
    ("Scripts/cli.py",
     "SKILLS dict (30 entries)\nshell command → skill module"),
    ("Scripts/_tools.py",
     "98 tool defs · JSON schemas\nexecute(name, args) → subprocess"),
    ("Scripts/_agent.py",
     "ReAct loop\nPlan → Action → Result → Evaluation\n10+ LLM backends"),
]):
    x = 5 + i * 31
    w = 28
    rounded_box(x, y4_bot, w, y4_top - y4_bot, fc="#FFFFFF", ec=L4, lw=1.8)
    text(x + w / 2, y4_top - 1.3, label, size=11.5, weight="bold", color=INK)
    text(x + w / 2, y4_top - 4.0, sub, size=8.5, color=SUBINK)

# Down arrow to Layer 3
arrow(50, y4_bot - 0.2, 50, 70.5, color=EDGE, lw=2.4, head=16)
text(52, 71.7, "execute(name, args) → parses stdout for artefacts",
     size=8.5, color=SUBINK, ha="left", style="italic")


# ─── LAYER 3 — 30 skills grouped by domain ───────────────────────────────────
layer_header(70.0, "Layer 3 · 30 skills, grouped by domain",
             "the producers", "#E8E5DD", INK)
y3_top = 67.5; y3_bot = 44

card_w = 28; card_h = 10.5
cols_x = [5, 36, 67]
rows_y = [y3_top - card_h, y3_bot]  # row 1 then row 2

for idx, (group_label, color, members) in enumerate(SKILL_GROUPS):
    cx = cols_x[idx % 3]
    cy = rows_y[idx // 3]
    is_integration = "★ network ★" in members
    border_lw = 3.2 if is_integration else 1.2
    rounded_box(cx, cy, card_w, card_h, fc=color, ec=("#7A1212" if is_integration else color),
                 lw=border_lw, alpha=0.92, rounding=0.9)
    text(cx + card_w / 2, cy + card_h - 2.0, group_label,
         size=12.5, weight="bold", color="#FFFFFF")
    text(cx + card_w / 2, cy + card_h / 2 - 0.8, members,
         size=9.2, color="#FFFFFF")
    if is_integration:
        # Star icon
        text(cx + 2.0, cy + card_h - 2.0, "★", size=18, color="#FFFF8A",
             weight="bold", ha="left")
        # Italic tag at the bottom of the card
        text(cx + card_w / 2, cy + 1.2,
             "apex of the skill DAG\n(CORNETO clean-room MILP)",
             size=8.0, color="#FFE6E6", style="italic")

# Arrow down to Layer 2
arrow(50, y3_bot - 0.2, 50, 40.5, color=EDGE, lw=2.4, head=16)
text(52, 41.7, "skills read from / write to local stores",
     size=8.5, color=SUBINK, ha="left", style="italic")


# ─── LAYER 2 — Local persistence ─────────────────────────────────────────────
layer_header(39.5, "Layer 2 · Local persistence",
             "two complementary stores", "#E2E8E0", INK)
y2_top = 37; y2_bot = 23

# Filesystem store (left)
rounded_box(3, y2_bot, 45, y2_top - y2_bot, fc=L2A, ec=L2A, lw=0, rounding=0.9)
text(25.5, y2_top - 1.7, "Filesystem (per-run, timestamped)",
     size=12, weight="bold", color="#FFFFFF")
text(25.5, y2_top - 4.2,
     "Data/Cache/   Data/Manifests/<Skill>/",
     size=9.2, color="#FFFFFF")
text(25.5, y2_top - 6.0,
     "Data/Input/   Data/Warehouse/KG/<collection>/",
     size=9.2, color="#FFFFFF")
text(25.5, y2_top - 9.4,
     "Docs/<Skill>/<ts>_*.md       ← human reports\n"
     "Docs/<Skill>/<ts>_*.svg/png  ← plots\n"
     "Docs/<Skill>/<ts>_*.csv/tsv  ← machine tables",
     size=8.6, color="#FFFFFF", ha="center")

# DuckDB warehouses (right)
rounded_box(52, y2_bot, 45, y2_top - y2_bot, fc=L2B, ec=L2B, lw=0, rounding=0.9)
text(74.5, y2_top - 1.4, "DuckDB warehouses",
     size=12, weight="bold", color="#FFFFFF")

# Silver tier inset
rounded_box(54.5, y2_bot + 7.0, 40, 5.0, fc="#FFFFFF", ec="#FFFFFF", lw=0,
            alpha=0.95, rounding=0.6)
text(56.5, y2_bot + 11.0, "SILVER", size=9.5, weight="bold",
     color="#7B5E2A", ha="left")
text(74.5, y2_bot + 11.0, "igvfagent_warehouse.duckdb",
     size=10.3, weight="bold", color=INK, ha="center")
text(74.5, y2_bot + 8.4,
     "proteomics_kg · perturb · multiseq\nmavedb · portal_kg · edges (inferred by network)",
     size=8.4, color=SUBINK, ha="center")

# Bronze tier inset
rounded_box(54.5, y2_bot + 0.8, 40, 5.0, fc="#FFFFFF", ec="#FFFFFF", lw=0,
            alpha=0.95, rounding=0.6)
text(56.5, y2_bot + 4.8, "BRONZE", size=9.5, weight="bold",
     color="#7A4A1F", ha="left")
text(74.5, y2_bot + 4.8, "igvf_kg_mirror.duckdb",
     size=10.3, weight="bold", color=INK, ha="center")
text(74.5, y2_bot + 2.2,
     "kg_<collection> views over zstd-Parquet shards\n"
     "55 Arango collections  (48/55 done · 1.0B rows)",
     size=8.4, color=SUBINK, ha="center")

# Bi-directional arrow to Layer 1
ax.add_patch(FancyArrowPatch((50, y2_bot - 0.2), (50, 20.5),
                              arrowstyle="<|-|>", mutation_scale=16,
                              color=EDGE, linewidth=2.4, zorder=5))
text(52, 21.4, "HTTP / HTTPS · cookie-auth Portal · guest Arango",
     size=8.5, color=SUBINK, ha="left", style="italic")


# ─── LAYER 1 — Upstream services ─────────────────────────────────────────────
layer_header(20.0, "Layer 1 · Upstream services",
             "no local copy — queried on demand", L1, "#FFFFFF")
y1_top = 17.5; y1_bot = 4.5

services = [
    ("data.igvf.org",          "IGVF Portal"),
    ("api.catalogkg.igvf.org", "IGVF Catalog REST"),
    ("db.catalog.igvf.org",    "IGVF KG (Arango)"),
    ("encodeproject.org",      "ENCODE"),
    ("api.genohub.org",        "FAVOR"),
    ("ncbi/geo · sra",         "NCBI GEO / SRA"),
    ("mavedb.org",              "MaveDB"),
    ("biogrid · intact · huri","PPI sources"),
    ("reactome · kegg",        "pathways"),
    ("chembl · pubmed",        "drugs · literature"),
]
# 2 rows × 5 columns
cell_w = 18; cell_h = 5.6; gap = 0.6
total_w = 5 * cell_w + 4 * gap
x0 = (100 - total_w) / 2
for i, (host, label) in enumerate(services):
    r = i // 5; c = i % 5
    x = x0 + c * (cell_w + gap)
    y = (y1_top - cell_h) - r * (cell_h + 1.0)
    rounded_box(x, y, cell_w, cell_h, fc="#FFFFFF", ec=L1, lw=1.4, rounding=0.5)
    text(x + cell_w / 2, y + cell_h - 1.7, host, size=9.0, weight="bold",
         color=INK)
    text(x + cell_w / 2, y + 1.7, label, size=8.0, color=SUBINK)


# ─── Footer ──────────────────────────────────────────────────────────────────
ax.add_patch(Rectangle((0, 0), 100, 3.2, color="#1F2933", zorder=0))
text(2.0, 1.8, "★  Knowledge integration",
     size=9.5, weight="bold", color="#F2CC8F", ha="left")
text(2.0, 0.7,
     "the `network` skill is the apex of the DAG — reads Silver + Bronze "
     "warehouses, runs CARNIVAL / Steiner MILP (pure cvxpy), writes "
     "inferred edges back",
     size=8.0, color="#D0D5DD", ha="left")
text(98.0, 1.6, "Apache-2.0",
     size=11, weight="bold", color="#F2CC8F", ha="right")
text(98.0, 0.6, "no GPL runtime deps",
     size=8.0, color="#D0D5DD", ha="right")


# ─── Save ────────────────────────────────────────────────────────────────────
out_dir = Path("Docs/Figures")
out_dir.mkdir(parents=True, exist_ok=True)
plt.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
for ext, kwargs in [("png", {"dpi": 300}), ("svg", {}), ("pdf", {})]:
    out = out_dir / f"IGVFagent_architecture_flow.{ext}"
    fig.savefig(out, bbox_inches="tight", facecolor=BG, **kwargs)
    print(f"Wrote {out}  ({out.stat().st_size:,} bytes)")
