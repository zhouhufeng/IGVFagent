#!/usr/bin/env python3
"""The bug-report tab must render, and must carry what makes a report usable.

Two failure modes, both silent on the site:

  * the tab raises and Streamlit shows a blank panel, which looks like the
    feature was never added;
  * it renders but omits the build id, so reports arrive that cannot be tied
    to a version of the code. The site is hot-copied several times a day, so
    "it gave a wrong answer" is unactionable without knowing which build.

The forum sets frame-ancestors 'self', so it cannot be embedded -- this is
asserted too, to stop anyone "improving" the tab into a blank iframe.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Scripts"))

FAILURES = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name:60} {detail}")
    if not ok:
        FAILURES.append(name)


class Stub:
    """Records every Streamlit call instead of rendering."""

    def __init__(self, sink):
        self._sink = sink
        self.session_state = {}

    def __getattr__(self, name):
        def f(*a, **k):
            self._sink.append((name, a[0] if a else None, k))
            return Stub(self._sink)
        return f


import streamlit_app as app  # noqa: E402

calls = []
try:
    app._render_feedback_tab(Stub(calls))
    rendered = True
except Exception as exc:                     # noqa: BLE001
    rendered = False
    print("   raised:", exc)
check("the tab renders without raising", rendered)

text = " ".join(str(c[1]) for c in calls if isinstance(c[1], str))
made = [c[0] for c in calls]

check("it links to the forum", "https://discussion.genohub.org/" in text)
check("the link is a button, not only prose", "link_button" in made)
check("it includes the build id, so a report names its code version",
      "build:" in text)
check("and the backend and model", "backend:" in text and "model:" in text)
check("and a session handle for log correlation", "session:" in text)
check("and a timestamp", "when:" in text)
check("it gives a structure to fill in",
      "What I asked" in text and "What happened instead" in text)
check("it says what the session id is and is not",
      "not an account" in text)
check("it does not promise to collect uploaded data",
      "No data you uploaded is included" in text)
check("it invites reports of false tool-unavailability claims",
      "unavailable" in text)

# The forum cannot be iframed; an embed would render blank. Scope this to
# the tab body -- the module legitimately iframes a data: URI elsewhere to
# preview a PDF, and a file-wide assertion flagged that instead.
import inspect  # noqa: E402
src = Path(app.__file__).read_text()
tab_src = inspect.getsource(app._render_feedback_tab)
# Ban the CALL, not the word: the docstring says "an iframe here would
# render a blank panel", which is the explanation worth keeping. An
# assertion on the word banned its own rationale.
code_only = "\n".join(
    l for l in tab_src.splitlines()
    if not l.lstrip().startswith("#") and '"""' not in l)
check("the tab does not try to iframe the forum",
      "components.iframe" not in code_only and "<iframe" not in code_only)
check("and says why, so nobody 'fixes' it into an embed",
      "frame-ancestors" in tab_src)
check("the URL is a module constant, not repeated inline",
      src.count('"https://discussion.genohub.org/"') == 1)
check("the sidebar also points at the forum, not only this tab",
      "Report it on the forum" in src)

print(f"\n{len(FAILURES)} failure(s)")
sys.exit(1 if FAILURES else 0)
