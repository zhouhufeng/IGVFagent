"""IGVFagent — local, auditable AI agent for IGVF / ENCODE data analysis.

The ``Scripts`` directory ships every individual skill as a standalone
Python module so that ``python3 Scripts/<skill>.py`` continues to work
exactly as the README documents. The same files are also exposed as the
installable Python package ``igvfagent`` (via the ``package-dir`` mapping
in ``pyproject.toml``), which gives end users:

* a single ``igvfagent`` console command (see ``cli.py``),
* importable submodules (``from igvfagent.igvf_client import ...``), and
* a coherent place to add the LLM router, ReAct agent runner, and UI.
"""

from __future__ import annotations

import os
from pathlib import Path

__version__ = "0.2.9"
__all__ = ["__version__", "load_dotenv"]


def load_dotenv() -> str | None:
    """Load the repo-root ``.env`` into ``os.environ`` (best-effort).

    Dependency-free (no ``python-dotenv`` required). Real environment
    variables always win — a key already exported in the shell is never
    overwritten by ``.env``. Without this, the native ``igvfagent``
    CLI/UI never sees ``.env`` (only Docker Compose injects it), so
    ``ANTHROPIC_API_KEY`` and friends silently appear "not set".

    Returns the path of the file that was loaded, or ``None``.
    """
    candidates = []
    root = os.environ.get("IGVF_PROJECT_ROOT")
    if root:
        candidates.append(Path(root) / ".env")
    # repo root is the parent of this Scripts/ package directory
    candidates.append(Path(__file__).resolve().parent.parent / ".env")
    candidates.append(Path.cwd() / ".env")

    for path in candidates:
        try:
            if not path.is_file():
                continue
        except OSError:
            continue
        for raw in path.read_text().splitlines():
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            if line.startswith("export "):
                line = line[len("export "):]
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            if key and key not in os.environ:   # real env wins
                os.environ[key] = value
        return str(path)
    return None


# Auto-load on import so `igvfagent ui` / `ask` and `python3 Scripts/*.py`
# all pick up the local .env without a manual `source`.
load_dotenv()
