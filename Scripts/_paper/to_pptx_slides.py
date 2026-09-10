#!/usr/bin/env python3
"""slides.html -> slides.pptx  (39 native PowerPoint slides)

The deck's slides are HTML fragments in a JS array. Each fragment is
parsed and re-laid out as real PowerPoint shapes — text frames, tables,
cards, stat blocks — rather than screenshotted, so every word stays
editable in PowerPoint.

Fonts are deliberately substituted for ones present on any machine
(Georgia / Calibri / Consolas) instead of the web faces, so the deck does
not silently fall back mid-talk.
"""
import sys, re, html as ihtml
from bs4 import BeautifulSoup, NavigableString
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

INK=RGBColor(0x14,0x16,0x2A); MUTED=RGBColor(0x5E,0x61,0x78)
INFER=RGBColor(0x3B,0x44,0xA6); MEAS=RGBColor(0xAF,0x5F,0x2B)
OK=RGBColor(0x2C,0x6A,0x4C); GROUND=RGBColor(0xF3,0xF4,0xF8)
SURF=RGBColor(0xFF,0xFF,0xFF); SURF2=RGBColor(0xF7,0xF8,0xFB)
RULE=RGBColor(0xDC,0xDD,0xE7); PAPER=RGBColor(0xFF,0xFF,0xFF)
SERIF="Georgia"; SANS="Calibri"; MONO="Consolas"

W,H = 13.333, 7.5
M = 0.62                      # side margin
CW = W - 2*M                  # content width


def txt_of(node):
    return re.sub(r"\s+"," ",node.get_text(" ",strip=True)).strip()


def emit_runs(par, node, bold=False, italic=False, mono=False, color=None, size=None):
    for ch in node.children:
        if isinstance(ch, NavigableString):
            t = re.sub(r"\s+"," ",str(ch))
            if not t.strip(): continue
            r = par.add_run(); r.text = t
            f = r.font
            f.bold, f.italic = bold, italic
            f.name = MONO if mono else SANS
            if size: f.size = Pt(size)
            f.color.rgb = color if color is not None else INK
            continue
        cls = ch.get("class") or []
        emit_runs(par, ch,
                  bold=bold or ch.name in ("b","strong"),
                  italic=italic or ch.name in ("em","i"),
                  mono=mono or ch.name=="code" or "mono" in cls,
                  color=(MEAS if "m" in cls and ch.name=="span" else color),
                  size=size)


def box(slide, x,y,w,h, fill=SURF, line=RULE, accent=None):
    sh = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x),Inches(y),Inches(w),Inches(h))
    sh.adjustments[0] = 0.045
    sh.fill.solid(); sh.fill.fore_color.rgb = fill
    if line is None: sh.line.fill.background()
    else: sh.line.color.rgb = line; sh.line.width = Pt(0.75)
    sh.shadow.inherit = False
    if accent is not None:
        bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x),Inches(y),Inches(0.045),Inches(h))
        bar.fill.solid(); bar.fill.fore_color.rgb = accent
        bar.line.fill.background(); bar.shadow.inherit=False
    return sh


def tb(slide, x,y,w,h):
    t = slide.shapes.add_textbox(Inches(x),Inches(y),Inches(w),Inches(h))
    tf = t.text_frame; tf.word_wrap = True
    tf.margin_left=tf.margin_right=tf.margin_top=tf.margin_bottom=0
    return tf


def line_para(tf, first, text, size, color, font=SANS, bold=False, space=4, align=None):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.space_after = Pt(space)
    r = p.add_run(); r.text = text
    r.font.size = Pt(size); r.font.color.rgb = color; r.font.name = font; r.font.bold = bold
    if align: p.alignment = align
    return p


def render_cards(slide, cards, x, y, w, ncol):
    gap = 0.22
    cw = (w - gap*(ncol-1))/ncol
    rows = (len(cards)+ncol-1)//ncol
    tall = any(c.find('div', class_='stat') for c in cards)
    ch = min(2.75 if tall else 2.05, (H - y - 0.70)/max(1,rows) - gap*0.6)
    for i, c in enumerate(cards):
        r, col = divmod(i, ncol)
        cx = x + col*(cw+gap); cy = y + r*(ch+gap)
        accent = INFER if "i" in (c.get("class") or []) else None
        box(slide, cx, cy, cw, ch, fill=SURF2, accent=accent)
        tf = tb(slide, cx+0.20, cy+0.15, cw-0.36, ch-0.28)
        first = True
        lbl = c.find("span", class_="lbl")
        if lbl:
            line_para(tf, first, txt_of(lbl).upper(), 8.5,
                      INFER if "i" in (lbl.get("class") or []) else MEAS, MONO, True, 3)
            first = False
        h3 = c.find("h3")
        if h3:
            line_para(tf, first, txt_of(h3), 13, INK, SANS, True, 3); first=False
        # A card may itself contain stat blocks (a headline number inside a
        # callout). Emit them, or the slide's most important figure vanishes.
        for sv in c.find_all("div", class_="stat"):
            sl = sv.find_next("div", class_="statlbl")
            line_para(tf, first, txt_of(sv), 20,
                      MEAS if "m" in (sv.get("class") or []) else INFER, SERIF, False, 1)
            first = False
            if sl and sl.find_parent("div", class_="card") is c:
                line_para(tf, False, txt_of(sl).upper(), 8, MUTED, MONO, False, 4)
        for p_ in c.find_all("p"):
            t = txt_of(p_)
            if t:
                line_para(tf, first, t, 10.5, MUTED, SANS, False, 3); first=False
    return y + rows*(ch+gap)


def render_stats(slide, groups, x, y, w):
    n = max(1,len(groups)); gap=0.3
    cw = (w-gap*(n-1))/n
    for i,(val,lbl,is_m) in enumerate(groups):
        cx = x + i*(cw+gap)
        tf = tb(slide, cx, y, cw, 1.6)
        line_para(tf, True, val, 34, MEAS if is_m else INFER, SERIF, False, 5)
        line_para(tf, False, lbl.upper(), 8.5, MUTED, MONO, False, 0)
    return y + 1.5


def render_table(slide, t, x, y, w):
    heads=[txt_of(th) for th in t.select("thead th")]
    body=t.select("tbody tr")
    rows=len(body)+1; cols=max(1,len(heads))
    th_ = min(0.34, (H-y-0.7)/rows)
    tbl = slide.shapes.add_table(rows, cols, Inches(x),Inches(y),Inches(w),Inches(th_*rows)).table
    for i,h in enumerate(heads):
        c=tbl.cell(0,i); c.text=""
        p=c.text_frame.paragraphs[0]; r=p.add_run(); r.text=h.upper()
        r.font.size=Pt(8); r.font.bold=True; r.font.name=MONO; r.font.color.rgb=MUTED
        c.fill.solid(); c.fill.fore_color.rgb=SURF2
    for ri,tr in enumerate(body,1):
        tds=tr.find_all(["td","th"]); istot="tot" in (tr.get("class") or [])
        for ci,td in enumerate(tds[:cols]):
            c=tbl.cell(ri,ci); c.text=""
            p=c.text_frame.paragraphs[0]
            emit_runs(p, td, size=9.5)
            for r in p.runs:
                r.font.size=Pt(9.5)
                if istot: r.font.bold=True
                if "yes" in (td.get("class") or []): r.font.color.rgb=OK
            c.fill.solid(); c.fill.fore_color.rgb = SURF2 if istot else PAPER
    for r in tbl.rows: r.height=Inches(th_)
    return y + th_*rows + 0.1


def render_term(slide, node, x, y, w):
    lines = node.get_text("\n").split("\n")
    lines = [l for l in lines if l.strip()!=""] or [""]
    h = min(3.0, 0.235*len(lines)+0.35)
    box(slide, x,y,w,h, fill=RGBColor(0x14,0x16,0x2A), line=None)
    tf = tb(slide, x+0.24, y+0.16, w-0.48, h-0.3)
    for i,l in enumerate(lines):
        col = RGBColor(0x8C,0x95,0xEC) if l.strip().startswith(">") else (
              RGBColor(0x63,0xB9,0x8C) if l.strip().startswith("$") else (
              RGBColor(0x9A,0x9C,0xB4) if l.strip().startswith("#") else RGBColor(0xDE,0x9A,0x62)))
        line_para(tf, i==0, l, 10.5, col, MONO, False, 1)
    return y+h+0.15


def render_slide(prs, title, frag, idx, total):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    s.background.fill.solid(); s.background.fill.fore_color.rgb = PAPER
    soup = BeautifulSoup(frag, "lxml")
    root = soup.body or soup
    node = root.find("div", class_="grow") or root
    y = 0.5

    kicker = node.find("p", class_="kicker")
    if kicker:
        tf = tb(s, M, y, CW, 0.3)
        line_para(tf, True, txt_of(kicker).upper(), 9.5,
                  INFER if "i" in (kicker.get("class") or []) else MEAS, MONO, True, 0)
        y += 0.36

    head = node.find(["h1","h2"])
    if head:
        t = txt_of(head)
        size = 40 if head.name=="h1" else (26 if len(t)>62 else 30)
        lines = max(1, len(t)//48 + 1)
        tf = tb(s, M, y, CW*0.92, 0.55*lines+0.2)
        line_para(tf, True, t, size, INK, SERIF, False, 0)
        y += 0.56*lines + 0.22

    for el in node.find_all(True, recursive=False):
        if el in (kicker, head) or y > H-0.95: continue
        cls = el.get("class") or []
        if el.name == "p":
            if "kicker" in cls: continue
            t = txt_of(el)
            if not t: continue
            big = "big" in cls; sub = "sub" in cls
            size = 17 if big else (12 if sub else 13)
            lines = max(1, len(t)//(78 if not big else 52) + 1)
            tf = tb(s, M, y, CW*0.88, 0.3*lines)
            p = tf.paragraphs[0]; p.space_after=Pt(0)
            emit_runs(p, el, size=size, color=MUTED if sub else INK)
            for r in p.runs:
                r.font.size = Pt(size); r.font.name = SERIF if big else SANS
                if not r.font.color or r.font.color.rgb is None: r.font.color.rgb = INK
            y += 0.27*lines + 0.14
        elif el.name == "ul":
            items = el.find_all("li", recursive=False)
            tf = tb(s, M, y, CW*0.9, 0.34*len(items))
            for i,li in enumerate(items):
                p = tf.paragraphs[0] if i==0 else tf.add_paragraph()
                p.space_after = Pt(6)
                b = p.add_run(); b.text = "— "
                b.font.color.rgb = INFER if "i" in cls else MEAS
                b.font.size = Pt(12.5); b.font.name = SANS
                emit_runs(p, li, size=12.5)
                for r in p.runs[1:]: r.font.size = Pt(12.5)
            y += 0.33*len(items) + 0.12
        elif el.name == "div" and ("cols" in cls or "row" in cls):
            cards = el.find_all("div", class_="card", recursive=False)
            stats = el.find_all("div", recursive=False)
            has_stat = any(d.find("div", class_="stat") for d in stats)
            if cards:
                ncol = 4 if "c4" in cls else (3 if "c3" in cls else (2 if "c2" in cls else 1))
                y = render_cards(s, cards, M, y, CW, ncol) + 0.1
            elif has_stat:
                # A single column may stack several stat pairs (e.g. a
                # 12.5 TB / 1.04 TB contrast separated by a rule), so collect
                # EVERY stat, not just the first in each child.
                groups=[]
                for d in stats:
                    for sv in d.find_all("div", class_="stat"):
                        sl = sv.find_next_sibling("div", class_="statlbl")
                        groups.append((txt_of(sv), txt_of(sl) if sl else "",
                                       "m" in (sv.get("class") or [])))
                if groups: y = render_stats(s, groups, M, y, CW) + 0.1
            else:
                for d in stats:
                    for sub in d.find_all(["h3","p","ul"], recursive=False):
                        t = txt_of(sub)
                        if not t: continue
                        tf = tb(s, M, y, CW*0.9, 0.3)
                        line_para(tf, True, t, 12, MUTED if sub.name=="p" else INK, SANS,
                                  sub.name=="h3", 0)
                        y += 0.3
        elif el.name == "div" and "card" in cls:
            y = render_cards(s, [el], M, y, CW, 1) + 0.1
        elif el.name == "div" and "term" in cls:
            y = render_term(s, el, M, y, CW)
        elif el.name == "table":
            y = render_table(s, el, M, y, CW)
        elif el.name == "div" and "rule" in cls:
            ln = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(M),Inches(y),Inches(CW*0.45),Inches(0.012))
            ln.fill.solid(); ln.fill.fore_color.rgb=RULE; ln.line.fill.background(); ln.shadow.inherit=False
            y += 0.22

    # footer
    tf = tb(s, M, H-0.42, CW, 0.25)
    p = tf.paragraphs[0]
    r = p.add_run(); r.text = f"IGVF Agent"
    r.font.size=Pt(8.5); r.font.color.rgb=MUTED; r.font.name=MONO
    tf2 = tb(s, W-M-1.2, H-0.42, 1.2, 0.25)
    p2 = tf2.paragraphs[0]; p2.alignment = PP_ALIGN.RIGHT
    r2 = p2.add_run(); r2.text = f"{idx:02d} / {total}"
    r2.font.size=Pt(8.5); r2.font.color.rgb=MUTED; r2.font.name=MONO
    return s


def main(src, out):
    raw = open(src, encoding="utf-8").read()
    blocks = re.findall(r'^S\("((?:[^"\\]|\\.)*)",\s*`(.*?)`\);\s*$', raw, re.S | re.M)
    prs = Presentation(); prs.slide_width=Inches(W); prs.slide_height=Inches(H)
    for i,(title, frag) in enumerate(blocks, 1):
        render_slide(prs, title, frag, i, len(blocks))
    prs.save(out)
    print(f"wrote {out} — {len(blocks)} slides")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
