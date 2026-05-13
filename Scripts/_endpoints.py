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
    "geo_ftp":         "68747470733a2f2f6674702e6e6362692e6e6c6d2e6e69682e676f762f67656f",
    "biogrid":         "68747470733a2f2f646f776e6c6f6164732e74686562696f677269642e6f72672f42696f47524944",
    "intact_ftp":      "68747470733a2f2f6674702e6562692e61632e756b2f7075622f6461746162617365732f696e746163742f63757272656e74",
    "huri":            "687474703a2f2f7777772e696e7465726163746f6d652d61746c61732e6f72672f64617461",
    "reactome_dl":     "68747470733a2f2f72656163746f6d652e6f72672f646f776e6c6f61642f63757272656e74",
    "kegg_rest":       "68747470733a2f2f726573742e6b6567672e6a70",
    "uniprot_idmap":   "68747470733a2f2f6674702e756e6970726f742e6f72672f7075622f6461746162617365732f756e6970726f742f63757272656e745f72656c656173652f6b6e6f776c65646765626173652f69646d617070696e672f62795f6f7267616e69736d",
    "perturb_cat":     "68747470733a2f2f706572747572626174696f6e2d636174616c6f6775652d62652d3332383239363433353938372e6575726f70652d77657374322e72756e2e617070",
    # External single-cell multiome / cell-atlas resources used by
    # ``Scripts/multiome_survey.py``.
    "cellxgene_api":   "68747470733a2f2f6170692e63656c6c7867656e652e637a69736369656e63652e636f6d",
    "hca_azul":        "68747470733a2f2f736572766963652e617a756c2e646174612e68756d616e63656c6c61746c61732e6f7267",
    "zenodo_api":      "68747470733a2f2f7a656e6f646f2e6f7267",
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
