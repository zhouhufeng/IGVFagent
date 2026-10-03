#!/usr/bin/env python3
"""Resolve IGVF Portal credentials from the environment or a secret file.

The Portal accepts an access-key pair as HTTP Basic credentials. Until now
only ``IGVF_ACCESS_KEY`` / ``IGVF_SECRET_ACCESS_KEY`` were consulted, which
means an operator who had downloaded a key pair from the Portal UI still had
to hand-transcribe it into two exported variables before anything
authenticated worked -- and an un-exported shell left every skill silently
anonymous, visible only as missing unreleased records.

So a file is read as well. The default location is
``Docs/Secret/IGVFportalAPI.txt``, which ``.gitignore`` already excludes
via ``**/Secret/``.

Precedence is environment first, file second: a deployment sets real
environment variables, and those must not be shadowed by a stale file left
in a developer checkout.

Nothing here logs or prints a secret. ``describe()`` exists so callers can
report *where* credentials came from, and reveals only the key id -- the
public half of the pair.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Optional, Tuple

__all__ = ["portal_credentials", "describe", "credential_file"]

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()

DEFAULT_FILE = ROOT / "Docs" / "Secret" / "IGVFportalAPI.txt"

# Labels the Portal UI uses when it shows a freshly minted key pair. Matched
# on normalised words rather than exact strings, because the wording is not
# stable: the Portal writes "Access Key ID" / "Access Key Secret", while the
# ENCODE-derived docs and some DACC tooling say "Secret Access Key" and
# "secret_access_key". Word-order-independent matching covers all of them.
_KEY_WORDS = {"access", "key", "id", "user", "igvf", "keyid", "accesskeyid"}
_SECRET_WORD = "secret"


def credential_file() -> Path:
    """Path of the credential file, overridable for tests and deployments."""
    override = os.environ.get("IGVF_PORTAL_API_FILE")
    return Path(override).expanduser() if override else DEFAULT_FILE


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9 ]+", " ", s.strip().lower()).strip()


def _classify(label: str) -> Optional[str]:
    """Map a label to "key", "secret", or None.

    Checked for "secret" anywhere in the label first, since every secret
    label also carries the key words -- "Access Key Secret" and "Secret
    Access Key" are the same field with the words swapped, and an
    order-sensitive test silently classifies one of them as the key id.
    """
    n = _norm(label)
    if not n:
        return None
    words = set(n.split())
    if _SECRET_WORD in words:
        return "secret"
    # A key label must be made only of key words, so a stray line of prose
    # is not mistaken for a label whose next line is a credential.
    if words and words <= _KEY_WORDS:
        return "key"
    return None


def _parse(text: str) -> Tuple[Optional[str], Optional[str]]:
    """Pull a key id and secret out of any of the plausible layouts.

    Handles the Portal's own copy-paste block, where a label sits on one line
    and its value on the next:

        Access Key ID
        ABC12345

        Secret Access Key
        xxxxxxxxxxxxxxxx

    and also ``key: value`` / ``key = value`` on one line, and a JSON object
    as written by some DACC tooling.
    """
    stripped = text.strip()
    if stripped.startswith("{"):
        try:
            obj = json.loads(stripped)
        except ValueError:
            obj = None
        if isinstance(obj, dict):
            lowered = {_norm(k): v for k, v in obj.items()}
            key = secret = None
            for k, v in lowered.items():
                slot = _classify(k)
                if slot == "key" and not key:
                    key = str(v).strip()
                elif slot == "secret" and not secret:
                    secret = str(v).strip()
            if key and secret:
                return key, secret

    lines = [ln.strip() for ln in text.splitlines()]
    key = secret = None
    pending: Optional[str] = None
    for ln in lines:
        if not ln or ln.startswith("#"):
            continue
        inline = re.match(r"^([A-Za-z][A-Za-z0-9 _-]*)\s*[:=]\s*(\S.*)$", ln)
        if inline:
            slot = _classify(inline.group(1))
            if slot == "key" and not key:
                key = inline.group(2).strip()
                pending = None
                continue
            if slot == "secret" and not secret:
                secret = inline.group(2).strip()
                pending = None
                continue
        slot = _classify(ln)
        if slot:
            # A bare label: its value is the next non-empty line.
            pending = slot
            continue
        if pending == "key" and not key:
            key = ln
        elif pending == "secret" and not secret:
            secret = ln
        pending = None
    return key or None, secret or None


def _from_file() -> Tuple[Optional[str], Optional[str], Optional[str]]:
    path = credential_file()
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return None, None, None
    key, secret = _parse(text)
    if key and secret:
        return key, secret, str(path)
    return None, None, str(path) if text.strip() else None


def portal_credentials() -> Optional[Tuple[str, str]]:
    """Return ``(key_id, secret)`` for HTTP Basic, or None if unavailable."""
    key = os.environ.get("IGVF_ACCESS_KEY")
    secret = os.environ.get("IGVF_SECRET_ACCESS_KEY")
    if key and secret:
        return (key.strip(), secret.strip())
    fkey, fsecret, _ = _from_file()
    if fkey and fsecret:
        return (fkey, fsecret)
    return None


def describe() -> dict:
    """Where credentials came from, for status output. No secret is included."""
    env_key = os.environ.get("IGVF_ACCESS_KEY")
    env_secret = os.environ.get("IGVF_SECRET_ACCESS_KEY")
    if env_key and env_secret:
        return {"source": "env", "key_id": env_key.strip(),
                "detail": "IGVF_ACCESS_KEY + IGVF_SECRET_ACCESS_KEY"}
    fkey, fsecret, path = _from_file()
    if fkey and fsecret:
        return {"source": "file", "key_id": fkey, "detail": path}
    partial = bool(env_key or env_secret)
    return {"source": "none", "key_id": None,
            "detail": ("only one of IGVF_ACCESS_KEY / "
                        "IGVF_SECRET_ACCESS_KEY is set — both are required"
                        if partial else
                        f"no credentials in env or {credential_file()}")}
