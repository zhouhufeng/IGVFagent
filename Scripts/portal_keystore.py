"""Each web user's own IGVF Portal keys: caught in chat, verified, kept encrypted.

People submitting through igvfagent.genohub.org paste their Portal access key
pair into the chat. The pair must never reach the model provider, the saved
transcript, or the shared history that recalls past answers for other users.
So the web app runs `intercept()` on every message BEFORE anything else sees
it: the pair is cut out, checked against the Portal hosts (keys are per host:
a staging key gets HTTP 401 on prod), stored under that user, and replaced in
the message by a placeholder. Submission tools then receive the user's own
keys in their subprocess environment (`env_for`), never the operator's.

Storage: Data/Keystore/<sha256(user)>.json, the secret encrypted with Fernet
(the `cryptography` package) under IGVF_KEYSTORE_SECRET, file mode 0600. With
no `cryptography` or no IGVF_KEYSTORE_SECRET, keys are held in memory for the
life of the process only, and `persistent()` says so.
"""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

HOSTS = {"staging": "https://api.staging.igvf.org", "prod": "https://api.data.igvf.org"}

# Portal key shapes: Access Key ID = 8 of [A-Z0-9], secret = 16 of [a-z0-9].
_ID = re.compile(r"(?<![A-Za-z0-9])(?=[A-Z0-9]*[A-Z])(?=[A-Z0-9]*[0-9])[A-Z0-9]{8}(?![A-Za-z0-9])")
_SECRET = re.compile(r"(?<![A-Za-z0-9])(?=[a-z0-9]*[a-z])(?=[a-z0-9]*[0-9])[a-z0-9]{16}(?![A-Za-z0-9])")
_HINT = re.compile(r"\b(key|keys|secret|access|credential|credentials|igvf|portal)\b", re.I)

_MEM: "dict[str, dict]" = {}
_LOCK = threading.Lock()


def find_pairs(text: str) -> "list[tuple[str, str]]":
    """Key ID / secret pairs in a message, in the order they appear."""
    if not text or not _HINT.search(text):
        return []
    ids = [m.group(0) for m in _ID.finditer(text)]
    secrets = [m.group(0) for m in _SECRET.finditer(text)]
    return list(zip(ids, secrets))


def redact(text: str, pairs: "list[tuple[str, str]]") -> str:
    for k, s in pairs:
        text = text.replace(k, "[IGVF key ID]").replace(s, "[IGVF secret]")
    return text


def verify(key: str, secret: str, mode: str, timeout: float = 20) -> "dict | None":
    """The account the pair signs in as on `mode`, or None if the host refuses it."""
    token = base64.b64encode(f"{key}:{secret}".encode()).decode()
    req = urllib.request.Request(f"{HOSTS[mode]}/session-properties", headers={
        "Authorization": f"Basic {token}", "Accept": "application/json",
        "User-Agent": "igvfagent-keystore"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = json.loads(r.read().decode() or "{}")
    except (urllib.error.HTTPError, urllib.error.URLError, ValueError):
        return None
    user = data.get("user") or {}
    if not user:
        return None
    labs = [str(x).rstrip("/").rsplit("/", 1)[-1] for x in (user.get("submits_for") or [])]
    return {"account": user.get("title") or user.get("email") or "?", "submits_for": labs}


def identify(key: str, secret: str) -> "list[tuple[str, dict]]":
    """Every host the pair works on (normally exactly one)."""
    out = []
    for mode in HOSTS:
        info = verify(key, secret, mode)
        if info:
            out.append((mode, info))
    return out


# ------------------------------ storage -------------------------------------


def _root() -> Path:
    return Path(os.environ.get("IGVF_PROJECT_ROOT") or Path(__file__).resolve().parents[1])


def _fernet():
    master = os.environ.get("IGVF_KEYSTORE_SECRET", "")
    if not master:
        return None
    try:
        from cryptography.fernet import Fernet
    except ImportError:
        return None
    return Fernet(base64.urlsafe_b64encode(hashlib.sha256(master.encode()).digest()))


def persistent() -> bool:
    """True when keys are written to disk encrypted; False = memory only."""
    return _fernet() is not None


def _user_id(user: str) -> str:
    return hashlib.sha256((user or "local").encode()).hexdigest()[:32]


def _path(user: str) -> Path:
    return _root() / "Data" / "Keystore" / f"{_user_id(user)}.json"


def _read(user: str) -> dict:
    f = _fernet()
    if f is None:
        with _LOCK:
            return dict(_MEM.get(_user_id(user), {}))
    p = _path(user)
    if not p.is_file():
        return {}
    try:
        return json.loads(p.read_text())
    except ValueError:
        return {}


def _write(user: str, data: dict) -> None:
    f = _fernet()
    if f is None:
        with _LOCK:
            _MEM[_user_id(user)] = data
        return
    p = _path(user)
    p.parent.mkdir(parents=True, exist_ok=True)
    os.chmod(p.parent, 0o700)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(data))
    os.chmod(tmp, 0o600)
    tmp.replace(p)


def save(user: str, mode: str, key: str, secret: str, info: dict) -> None:
    data = _read(user)
    f = _fernet()
    blob = f.encrypt(secret.encode()).decode() if f else secret
    data[mode] = {"key_id": key, "secret": blob, "encrypted": bool(f),
                  "account": info.get("account"), "submits_for": info.get("submits_for", []),
                  "saved_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
    _write(user, data)


def register(user: str, mode: str, key: str, secret: str) -> "tuple[bool, str]":
    """Save a pair typed into the keys panel, once the chosen host accepts it.

    Each host has its own slot, so a staging pair and a prod pair are kept
    side by side; saving one never touches the other.
    """
    key, secret = (key or "").strip(), (secret or "").strip()
    if mode not in HOSTS:
        return False, f"Unknown Portal {mode!r}."
    if not (_ID.fullmatch(key) and _SECRET.fullmatch(secret)):
        return False, ("That does not look like a Portal key pair: the Access Key ID is 8 capital "
                       "letters and digits, and the Secret 16 lower-case letters and digits.")
    info = verify(key, secret, mode)
    if not info:
        other = next((m for m in HOSTS if m != mode and verify(key, secret, m)), None)
        return False, (f"{HOSTS[mode]} did not accept this pair." + (
            f" It works on {other} instead: keys belong to one host, so save it as the {other} key."
            if other else " Check it, or create a new pair in your profile on that host."))
    save(user, mode, key, secret, info)
    labs = ", ".join(info.get("submits_for") or []) or "no lab"
    where = "saved encrypted on the server" if persistent() else "kept for this server session only"
    return True, f"{mode} keys saved for {info['account']} (can submit for {labs}); {where}."


def load(user: str, mode: str) -> "tuple[str, str] | None":
    rec = _read(user).get(mode)
    if not rec:
        return None
    if rec.get("encrypted"):
        f = _fernet()
        if f is None:
            return None
        try:
            secret = f.decrypt(rec["secret"].encode()).decode()
        except Exception:
            return None
    else:
        secret = rec["secret"]
    return rec["key_id"], secret


def listing(user: str) -> "list[dict]":
    """What is on file, safe to show: no secret, the key ID masked."""
    return [{"mode": m, "key_id": r["key_id"][:2] + "******", "account": r.get("account"),
             "submits_for": r.get("submits_for", []), "saved_at": r.get("saved_at")}
            for m, r in sorted(_read(user).items())]


def delete(user: str, mode: "str | None" = None) -> None:
    data = _read(user)
    if mode:
        data.pop(mode, None)
    else:
        data = {}
    if data or _fernet() is None:
        _write(user, data)
    else:
        try:
            _path(user).unlink()
        except FileNotFoundError:
            pass


# Every name a Portal key pair travels under. Submission tools in a signed-in
# web session get the user's pair under the per-host names and none of these
# from the operator's own environment.
PORTAL_ENV = ("IGVF_API_KEY", "IGVF_SECRET_KEY", "IGVF_ACCESS_KEY", "IGVF_SECRET_ACCESS_KEY",
              "IGVF_PROD_API_KEY", "IGVF_PROD_SECRET_KEY", "IGVF_STAGING_API_KEY",
              "IGVF_STAGING_SECRET_KEY", "IGVF_PORTAL_API_FILE")


def env_for(user: str) -> dict:
    """Subprocess environment entries carrying this user's keys, per host."""
    out = {}
    for mode in HOSTS:
        pair = load(user, mode)
        if pair:
            out[f"IGVF_{mode.upper()}_API_KEY"], out[f"IGVF_{mode.upper()}_SECRET_KEY"] = pair
    return out


def intercept(text: str, user: str) -> "tuple[str, list[str]]":
    """Cut Portal keys out of a chat message, store the ones that work.

    Returns the message with the keys replaced, and one line per pair saying
    what happened (shown to the user; safe for the model and transcripts).
    """
    pairs = find_pairs(text)
    if not pairs:
        return text, []
    notes = []
    for key, secret in pairs:
        hosts = identify(key, secret)
        if not hosts:
            notes.append("A Portal key pair was found and removed from your message, but "
                         "neither staging nor prod accepted it, so it was not saved. Check "
                         "the pair, or issue a new one from your profile on that host.")
            continue
        for mode, info in hosts:
            save(user, mode, key, secret, info)
            labs = ", ".join(info.get("submits_for") or []) or "no lab"
            where = "saved encrypted on the server" if persistent() else "kept for this server session only"
            notes.append(f"IGVF {mode} keys received for {info['account']} (can submit for: "
                         f"{labs}); {where}. They were removed from your message and are not "
                         f"shown to the model.")
    return redact(text, pairs), notes
