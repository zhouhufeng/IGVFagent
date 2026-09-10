#!/usr/bin/env python3
"""Do the repository's markdown files actually render?

Added after unescaped pipes were found breaking tables in nine rows across
five files -- `|z| > 2.5`, `|abundance change|`, a regex alternation
`(NM|NP|ENST|ENSP)`, and shell pipes inside backticked commands. Markdown
reads every unescaped `|` as a column break, so those rows rendered with
more columns than their header and the table came apart. Two of them sat in
the root README's attribution table, which is the first thing a reader sees.

None of it is visible in a diff or caught by a test suite; it only shows up
in a rendered view that nobody opens for a benchmark's OPERATIONS.md.

    python3 Benchmarks/check_markdown.py $(git ls-files '*.md')

Checks: table rows whose unescaped-pipe count differs from their header,
unclosed code fences, in-page anchors matching no heading in the same file,
and relative links to files that do not exist. Exits non-zero if anything
fails, so it can gate a commit.
"""

import re
import sys
import pathlib

def cols(line):
    """Columns as Markdown sees them: unescaped pipes only."""
    return len(re.findall(r'(?<!\\)\|', line)) - 1

def check(path):
    lines = pathlib.Path(path).read_text().splitlines()
    bad, in_tbl, hdr, hdr_ln, in_fence = [], False, 0, 0, False
    for i, l in enumerate(lines, 1):
        if l.lstrip().startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        is_row = l.lstrip().startswith("|")
        if is_row and not in_tbl:
            in_tbl, hdr, hdr_ln = True, cols(l), i
            continue
        if is_row and in_tbl:
            if set(l) <= set("|-: "):        # separator row
                continue
            c = cols(l)
            if c != hdr:
                bad.append((i, c, hdr, hdr_ln, l.strip()[:70]))
        elif not is_row:
            in_tbl = False
    return bad

total = 0
for p in sys.argv[1:]:
    b = check(p)
    total += len(b)
    if b:
        print(f"\n{p}: {len(b)} malformed row(s)")
        for i, c, h, hl, t in b:
            print(f"   line {i}: {c} cols (header at {hl} has {h}) -> {t}")
print(f"\n{total} malformed table row(s) across {len(sys.argv)-1} file(s)")

# ---- fences, anchors, relative links -------------------------------------
unclosed, bad_anchor, bad_rel = [], [], []
for f in sys.argv[1:]:
    p = pathlib.Path(f)
    try:
        s = p.read_text()
    except OSError:
        continue
    if s.count("\n```") % 2 != 0:
        unclosed.append(f)
    heads = {re.sub(r"[^a-z0-9 -]", "", h.lower()).replace(" ", "-")
             for h in re.findall(r"^#{1,6} (.+)$", s, re.M)}
    bad_anchor += [(f, a) for a in set(re.findall(r"\]\(#([A-Za-z0-9_-]+)\)", s))
                   if a.lower() not in heads]
    for link in set(re.findall(r"\]\((?!https?:|#|mailto:)([^)\s]+)\)", s)):
        t = link.split("#")[0]
        if t and not (p.parent / t).exists() and not pathlib.Path(t).exists():
            bad_rel.append((f, t))

for label, items in (("unclosed code fences", unclosed),
                      ("broken in-page anchors", bad_anchor),
                      ("broken relative links", sorted(set(bad_rel)))):
    print(f"{label}: {len(items)}")
    for it in items[:10]:
        print(f"   {it}")

failures = total + len(unclosed) + len(bad_anchor) + len(set(bad_rel))
print(f"\nTOTAL problems: {failures}")
sys.exit(1 if failures else 0)
