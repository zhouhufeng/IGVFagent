#!/usr/bin/env python3
"""An artefact must be drawn once, not once per renderer.

Reported across three retests and "fixed" twice before this: figures appeared
twice in a completed turn. The cause was not the dedup logic, which was
correct -- it was that TWO renderers run over the same turn.
_render_markdown_with_images draws figures the answer cites with ![](...)
syntax; _render_artefacts draws the run's declared artefacts. A run that both
writes a figure and mentions it -- the normal shape -- drew it twice.

The mechanism to prevent that already existed: the inline renderer RETURNS
the paths it drew, documented as being "so callers can dedup against their
own linked artefacts list". Every chat call site discarded the return value.
"""
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Scripts"))

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:62} {detail}")
    if not ok:
        FAILURES.append(name)


_tmp = tempfile.TemporaryDirectory()
root = Path(_tmp.name).resolve()
os.environ["IGVF_PROJECT_ROOT"] = str(root)
(root / "Docs" / "Run1").mkdir(parents=True)
fig = root / "Docs" / "Run1" / "figure.png"
fig.write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * 64)
tbl = root / "Docs" / "Run1" / "table.csv"
tbl.write_text("a,b\n1,2\n")

import streamlit_app as app  # noqa: E402
app._PROJECT_ROOT = root

# The figure the answer cited inline must not reappear in the panel.
inline = {str(fig)}
shown = app._artefacts_to_show([str(fig), str(tbl)], inline)
check("a figure drawn inline is not drawn again in the panel",
      str(fig) not in shown, str(shown))
check("but other artefacts still appear", str(tbl) in shown)

# The spelling differences that defeated the earlier string-based fix.
rel = "Docs/Run1/figure.png"
shown = app._artefacts_to_show([rel], {str(fig)})
check("a RELATIVE path is matched against an absolute inline render",
      shown == [], str(shown))
shown = app._artefacts_to_show([str(fig)], {rel})
check("and the reverse spelling too", shown == [], str(shown))

# The whole-directory case: announcing the run dir AND a file in it.
shown = app._artefacts_to_show([str(root / "Docs" / "Run1"), str(fig)], set())
check("a directory and a file inside it collapse to one entry per file",
      len(shown) == len(set(shown)), str(shown))
check("and the directory itself expands to its files",
      any(s.endswith("figure.png") for s in shown))

# Nothing left to show must mean an empty list, so the caller can skip the
# panel entirely rather than render an empty one under a non-zero count.
shown = app._artefacts_to_show([str(fig)], {str(fig)})
check("everything-already-inline yields an empty list", shown == [])

# The header count and the contents come from the same call, so they cannot
# disagree -- the bug this would otherwise introduce.
src = Path(app.__file__).read_text()
check("the panel header counts the FILTERED list",
      "shown = _artefacts_to_show(artefacts, inline_rendered)" in src
      and "Artefacts ({len(shown)})" in src)
check("the replay path passes its inline renders too",
      "already_rendered=_inline" in src)
check("the live turn captures what it rendered inline",
      "inline_rendered = _render_markdown_with_images(" in src)
check("no chat call site discards the inline return any more",
      src.count("            _render_markdown_with_images(result.final_answer,") == 0)

print(f"\n{len(FAILURES)} failure(s)")
_tmp.cleanup()
sys.exit(1 if FAILURES else 0)
