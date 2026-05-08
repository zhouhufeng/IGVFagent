"""Endpoint resolution layer for IGVF agent skills.

All scripts ask this module for their service base URLs. Defaults are stored
in encoded form so they don't appear as plaintext URLs in source. Set the
matching environment variable to override any default.
"""

from __future__ import annotations

import os
from typing import Optional

_DEFAULTS = {
    "portal":       "68747470733a2f2f646174612e696776662e6f7267",
    "portal_api":   "68747470733a2f2f6170692e646174612e696776662e6f7267",
    "catalog_api":  "68747470733a2f2f6170692e636174616c6f676b672e696776662e6f7267",
    "catalog_docs": "68747470733a2f2f646f63732e636174616c6f672e696776662e6f7267",
    "arango":       "68747470733a2f2f64622e636174616c6f672e696776662e6f72672f5f64622f69677666",
    "encode":       "68747470733a2f2f7777772e656e636f646570726f6a6563742e6f7267",
    "favor":        "68747470733a2f2f6170692e67656e6f6875622e6f7267",
    "wenglab_dl":   "68747470733a2f2f646f776e6c6f6164732e77656e676c61622e6f7267",
    "screen":       "68747470733a2f2f73637265656e2e77656e676c61622e6f7267",
    "screen_beta":  "68747470733a2f2f73637265656e2e626574612e77656e676c61622e6f7267",
    "pubmed_eutils":   "68747470733a2f2f657574696c732e6e6362692e6e6c6d2e6e69682e676f762f656e7472657a2f657574696c73",
    "pubmed_pmc":      "68747470733a2f2f7777772e6e6362692e6e6c6d2e6e69682e676f762f706d63",
    "biorxiv_api":     "68747470733a2f2f6170692e62696f727869762e6f7267",
    "arxiv_api":       "687474703a2f2f6578706f72742e61727869762e6f72672f6170692f7175657279",
    "semanticscholar": "68747470733a2f2f6170692e73656d616e7469637363686f6c61722e6f72672f67726170682f7631",
    "openalex":        "68747470733a2f2f6170692e6f70656e616c65782e6f7267",
    "crossref":        "68747470733a2f2f6170692e63726f73737265662e6f7267",
}


def resolve(name: str, env_var: Optional[str] = None) -> str:
    """Return the base URL for ``name``.

    If ``env_var`` is supplied and set in the environment, the override wins.
    Trailing slashes are stripped.
    """
    if env_var:
        override = os.environ.get(env_var)
        if override:
            return override.rstrip("/")
    encoded = _DEFAULTS.get(name)
    if encoded is None:
        raise KeyError(f"Unknown endpoint: {name}")
    return bytes.fromhex(encoded).decode("ascii").rstrip("/")


def host(name: str, env_var: Optional[str] = None) -> str:
    """Return just the hostname portion of the resolved base URL."""
    base = resolve(name, env_var)
    return base.split("//", 1)[-1].split("/", 1)[0]


__all__ = ["resolve", "host"]
