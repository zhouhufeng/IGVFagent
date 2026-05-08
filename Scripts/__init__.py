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

__version__ = "0.1.0"
__all__ = ["__version__"]
