"""DiscourseConnect authentication gate for IGVFagent.

Replaces the single shared password with real accounts. Identity comes from
the Genohub Discourse forum, which already has a signup flow, email
verification and staff moderation -- there is no reason to build a second
credential store, and a homegrown one on a public host is the thing that gets
breached.

The flow, once per browser session:

  1. nginx asks this service (``/_auth/verify``, via auth_request) whether the
     request carries a valid session cookie. No cookie -> 401.
  2. nginx turns that 401 into a redirect to ``/_auth/login``.
  3. ``/_auth/login`` signs a nonce with the shared secret and bounces the
     browser to Discourse's ``/session/sso_provider``.
  4. Discourse authenticates the person (logging in or signing up) and
     redirects back to ``/_auth/callback`` with a signed payload carrying
     their email, username and GROUP MEMBERSHIPS.
  5. We verify the signature, spend the nonce, and require membership of the
     approval group. Approved -> a signed session cookie. Not approved -> a
     page telling them exactly what to ask for.
  6. Every later request is settled by the cookie alone: one HMAC, no I/O, no
     network call to Discourse.

**Group membership is the approval step.** A person can sign up at the forum
freely; they cannot use IGVFagent until an admin adds them to the
``igvfagent`` group. That keeps approval where the admins already work instead
of in a bespoke queue nobody checks.

Authentication deliberately lives HERE and not inside the Streamlit app. The
app runs analysis CLIs as subprocesses on prompts from the public internet, so
it is the thing being protected; it must never also be the thing deciding who
gets in.

Stdlib only, on purpose: an auth boundary with no dependency tree is one that
cannot be compromised through its dependency tree.

Configuration (environment):

  IGVF_SSO_SECRET        shared with Discourse's `discourse connect provider
                         secrets` setting, as `igvfagent.genohub.org|<secret>`
  IGVF_SESSION_SECRET    signs our own session cookie; unrelated to the above
  IGVF_DISCOURSE_URL     https://discussion.genohub.org
  IGVF_PUBLIC_URL        https://igvfagent.genohub.org
  IGVF_APPROVAL_GROUP    Discourse group granting access (default: igvfagent)
  IGVF_SESSION_HOURS     cookie lifetime (default: 168 = one week)
  IGVF_DISCOURSE_API_KEY optional; only consulted if a payload arrives with no
                         groups field at all (see _groups_for)
  IGVF_DISCOURSE_API_USER  optional; username the API key acts as
  IGVF_BOOTSTRAP_TOKEN   optional break-glass; see BOOTSTRAP_TOKEN below.
                         Unset it once Discourse SSO is confirmed working.
"""
from __future__ import annotations

import base64
import binascii
import hashlib
import hmac
import json
import logging
import os
import secrets
import sys
import threading
import time
import urllib.parse
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

log = logging.getLogger("igvf-auth")

SSO_SECRET = os.environ.get("IGVF_SSO_SECRET", "").encode()
SESSION_SECRET = os.environ.get("IGVF_SESSION_SECRET", "").encode()
DISCOURSE_URL = os.environ.get(
    "IGVF_DISCOURSE_URL", "https://discussion.genohub.org").rstrip("/")
PUBLIC_URL = os.environ.get(
    "IGVF_PUBLIC_URL", "https://igvfagent.genohub.org").rstrip("/")
APPROVAL_GROUP = os.environ.get("IGVF_APPROVAL_GROUP", "igvfagent").strip()
SESSION_HOURS = int(os.environ.get("IGVF_SESSION_HOURS", "168"))
API_KEY = os.environ.get("IGVF_DISCOURSE_API_KEY", "")
API_USER = os.environ.get("IGVF_DISCOURSE_API_USER", "system")
# Break-glass. Optional and unset by default. Turning a live site over to a
# new identity provider in one step means that if anything about the Discourse
# side is not yet right -- the provider setting not enabled, the secret
# mistyped, the group not created -- NOBODY can get in, including the person
# who would fix it. This gives one URL that issues a session without Discourse,
# for verifying the stack before the cutover and for recovering after one.
# It is exactly as strong as the shared password it replaces, which is why it
# should be removed from the environment once SSO is confirmed working.
BOOTSTRAP_TOKEN = os.environ.get("IGVF_BOOTSTRAP_TOKEN", "")
COOKIE = "igvf_session"
PORT = int(os.environ.get("IGVF_AUTH_PORT", "9000"))

# Nonces we have issued but not yet seen come back. Single-use and
# short-lived, so a callback URL cannot be replayed out of a browser history
# or a referer log. One process, one dict -- this is a single-replica
# deployment and a shared store would be a dependency bought for nothing.
_NONCES: "dict[str, float]" = {}
_NONCE_TTL = 600.0
_NONCE_LOCK = threading.Lock()


def _now() -> float:
    return time.time()


def _issue_nonce() -> str:
    n = secrets.token_urlsafe(24)
    with _NONCE_LOCK:
        cutoff = _now() - _NONCE_TTL
        for k, born in list(_NONCES.items()):
            if born < cutoff:
                del _NONCES[k]
        _NONCES[n] = _now()
    return n


def _spend_nonce(n: str) -> bool:
    with _NONCE_LOCK:
        born = _NONCES.pop(n, None)
    return born is not None and (_now() - born) <= _NONCE_TTL


# ------------------------------ signing -------------------------------------


def _truthy(value) -> bool:
    """Discourse sends booleans as the STRINGS "true"/"false".

    ``bool("false")`` is True, so a plain bool() here would have marked every
    single user an administrator. Parse the wire format, do not coerce it.
    """
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("true", "1", "yes")


def _sig(secret: bytes, payload: bytes) -> str:
    return hmac.new(secret, payload, hashlib.sha256).hexdigest()


def _b64u(raw: bytes) -> str:
    return base64.urlsafe_b64encode(raw).decode().rstrip("=")


def _unb64u(text: str) -> bytes:
    return base64.urlsafe_b64decode(text + "=" * (-len(text) % 4))


def make_session(user: dict) -> str:
    """A self-contained signed cookie: identity plus an expiry, nothing else.

    No server-side session table, so a restart does not sign everyone out and
    there is no store to grow unbounded. The cost is that revocation waits for
    expiry -- acceptable for a lab tool, and rotating IGVF_SESSION_SECRET
    invalidates every session at once if it is ever needed.
    """
    body = json.dumps({
        "u": user.get("username", ""),
        "e": user.get("email", ""),
        "n": user.get("name", ""),
        "i": user.get("external_id", ""),
        "a": _truthy(user.get("admin")),
        "x": int(_now() + SESSION_HOURS * 3600),
    }, separators=(",", ":"), sort_keys=True).encode()
    return f"{_b64u(body)}.{_sig(SESSION_SECRET, body)}"


def read_session(raw: str) -> "dict | None":
    if not raw or "." not in raw:
        return None
    encoded, _, given = raw.rpartition(".")
    try:
        body = _unb64u(encoded)
    except (binascii.Error, ValueError):
        return None
    if not hmac.compare_digest(_sig(SESSION_SECRET, body), given):
        return None
    try:
        data = json.loads(body)
    except ValueError:
        return None
    if not isinstance(data, dict) or int(data.get("x", 0)) < _now():
        return None
    return data


# ------------------------------ Discourse -----------------------------------


def _groups_for(payload: "dict[str, str]") -> "list[str]":
    """Group names for the person who just authenticated.

    Discourse puts ``groups`` in the DiscourseConnect payload, and that is the
    path taken in practice: it is already signed by the same secret, so it
    costs nothing and cannot be forged. The API fallback exists only for the
    case where the field is absent entirely -- and if neither is available we
    return None-ish (an empty list) so the caller DENIES. Failing closed
    matters more here than convenience: a bug that silently grants access to
    everyone who can sign up at a public forum is the worst outcome this file
    can have.
    """
    raw = payload.get("groups")
    if raw is not None:
        return [g.strip() for g in raw.split(",") if g.strip()]
    if not API_KEY:
        log.error("Discourse sent no `groups` field and no API key is "
                  "configured -- denying. Set IGVF_DISCOURSE_API_KEY, or "
                  "check that this Discourse version includes groups in the "
                  "SSO provider payload.")
        return []
    username = payload.get("username", "")
    if not username:
        return []
    url = f"{DISCOURSE_URL}/u/{urllib.parse.quote(username)}.json"
    req = urllib.request.Request(url, headers={
        "Api-Key": API_KEY, "Api-Username": API_USER,
        "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            blob = json.loads(resp.read().decode())
    except Exception as exc:
        log.error("Discourse group lookup failed for %s: %s", username, exc)
        return []
    return [g.get("name", "") for g in (blob.get("user", {}).get("groups") or [])]


def approved(payload: "dict[str, str]") -> bool:
    if _truthy(payload.get("admin")):
        return True
    wanted = APPROVAL_GROUP.lower()
    return any(g.lower() == wanted for g in _groups_for(payload))


# ------------------------------ pages ---------------------------------------

_PAGE = """<!doctype html><meta charset=utf-8>
<meta name=viewport content="width=device-width,initial-scale=1">
<title>{title} — IGVF Agent</title>
<style>
 :root {{ color-scheme: light dark; }}
 body {{ margin:0; min-height:100vh; display:grid; place-items:center;
        font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
        background:#f6f7f8; color:#1b1b1b; padding:24px; }}
 @media (prefers-color-scheme:dark) {{ body {{ background:#14181b; color:#e8eaed; }} }}
 .card {{ max-width:30rem; text-align:center; }}
 h1 {{ font-size:1.35rem; margin:0 0 .75rem; }}
 p {{ margin:0 0 1rem; }}
 code {{ background:rgba(127,127,127,.18); padding:.1em .4em; border-radius:4px; }}
 a.btn {{ display:inline-block; background:#38707f; color:#fff; text-decoration:none;
         padding:.6em 1.3em; border-radius:6px; font-weight:600; }}
 .muted {{ opacity:.7; font-size:.9rem; }}
</style>
<div class=card>{body}</div>
"""


def page(title: str, body: str) -> bytes:
    return _PAGE.format(title=title, body=body).encode()


def _denied_page(payload: "dict[str, str]") -> bytes:
    who = payload.get("username") or payload.get("email") or "your account"
    return page("Access pending", f"""
<h1>Not yet approved for IGVF Agent</h1>
<p>You are signed in to the Genohub community as <code>{_esc(who)}</code>,
   but that account is not in the <code>{_esc(APPROVAL_GROUP)}</code> group
   yet.</p>
<p>Ask an IGVF Agent administrator to add you. Once they do, come back to this
   page — nothing else is needed on your side.</p>
<p><a class=btn href="/_auth/login">Try again</a></p>
<p class=muted>Signing up to the forum and being approved for the agent are
   two separate steps, on purpose: the agent runs analyses on a shared
   machine.</p>""")


def _esc(text: str) -> str:
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


# ------------------------------ handler -------------------------------------


class Handler(BaseHTTPRequestHandler):
    server_version = "igvf-auth"
    protocol_version = "HTTP/1.1"

    # The default logs a line per request to stderr, which at auth_request
    # volume means a line per subresource. nginx already has the access log.
    def log_message(self, fmt, *args):  # noqa: A003
        pass

    # -- helpers --
    def _cookies(self) -> "dict[str, str]":
        out: "dict[str, str]" = {}
        for part in (self.headers.get("Cookie") or "").split(";"):
            name, _, value = part.strip().partition("=")
            if name:
                out[name] = value
        return out

    def _send(self, code: int, body: bytes = b"", *,
              ctype: str = "text/html; charset=utf-8",
              headers: "list[tuple[str, str]] | None" = None) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for key, value in (headers or []):
            self.send_header(key, value)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _redirect(self, where: str,
                  headers: "list[tuple[str, str]] | None" = None) -> None:
        self._send(302, b"", headers=[("Location", where)] + list(headers or []))

    def _set_cookie(self, value: str, *, max_age: int) -> "tuple[str, str]":
        # Secure is unconditional: the only way in is through Cloudflare over
        # TLS, and a cookie that would also travel over plain HTTP is a
        # cookie that leaks on the one misconfiguration that matters.
        return ("Set-Cookie",
                f"{COOKIE}={value}; Path=/; HttpOnly; Secure; "
                f"SameSite=Lax; Max-Age={max_age}")

    # -- routes --
    def do_GET(self) -> None:  # noqa: N802
        path, _, query = self.path.partition("?")
        args = urllib.parse.parse_qs(query)
        if path == "/_auth/verify":
            return self._verify()
        if path == "/_auth/login":
            return self._login(args)
        if path == "/_auth/callback":
            return self._callback(args)
        if path == "/_auth/bootstrap":
            return self._bootstrap(args)
        if path == "/_auth/logout":
            return self._logout()
        if path == "/_auth/health":
            return self._send(200, b"ok", ctype="text/plain")
        self._send(404, page("Not found", "<h1>Not found</h1>"))

    do_HEAD = do_GET

    def _verify(self) -> None:
        """nginx auth_request target. Cookie only -- no network, no disk."""
        data = read_session(self._cookies().get(COOKIE, ""))
        if not data:
            return self._send(401, b"", ctype="text/plain")
        # Identity for the app, read back out of the subrequest by
        # auth_request_set and forwarded upstream.
        self._send(200, b"", ctype="text/plain", headers=[
            ("X-IGVF-User", data.get("u", "")),
            ("X-IGVF-Email", data.get("e", "")),
            ("X-IGVF-Name", data.get("n", "")),
            ("X-IGVF-Id", str(data.get("i", ""))),
            ("X-IGVF-Admin", "1" if data.get("a") else "0"),
        ])

    def _login(self, args) -> None:
        if not SSO_SECRET or not SESSION_SECRET:
            return self._send(500, page("Misconfigured", """
<h1>Authentication is not configured</h1>
<p>This deployment has no SSO secret set, so nobody can sign in. The operator
   needs to set <code>IGVF_SSO_SECRET</code> and
   <code>IGVF_SESSION_SECRET</code>.</p>"""))
        nonce = _issue_nonce()
        payload = urllib.parse.urlencode({
            "nonce": nonce,
            "return_sso_url": f"{PUBLIC_URL}/_auth/callback",
        }).encode()
        sso = base64.b64encode(payload).decode()
        url = (f"{DISCOURSE_URL}/session/sso_provider?"
               + urllib.parse.urlencode({"sso": sso,
                                          "sig": _sig(SSO_SECRET, sso.encode())}))
        self._redirect(url)

    def _callback(self, args) -> None:
        sso = (args.get("sso") or [""])[0]
        sig = (args.get("sig") or [""])[0]
        if not sso or not sig:
            return self._send(400, page("Bad request",
                                        "<h1>Missing SSO response</h1>"))
        if not hmac.compare_digest(_sig(SSO_SECRET, sso.encode()), sig):
            log.warning("rejected callback with a bad signature")
            return self._send(403, page("Rejected", """
<h1>That sign-in could not be verified</h1>
<p>The response did not carry a valid signature. If this keeps happening the
   shared secret here and in Discourse do not match.</p>"""))
        try:
            fields = urllib.parse.parse_qs(base64.b64decode(sso).decode())
        except (binascii.Error, ValueError, UnicodeDecodeError):
            return self._send(400, page("Bad request",
                                        "<h1>Malformed SSO response</h1>"))
        payload = {k: v[0] for k, v in fields.items() if v}

        if not _spend_nonce(payload.get("nonce", "")):
            # Either a replay or a login that sat on the Discourse page for
            # more than ten minutes. Both are fixed by starting again.
            return self._redirect("/_auth/login")

        if not approved(payload):
            log.info("denied %s -- not in %s",
                     payload.get("username", "?"), APPROVAL_GROUP)
            return self._send(403, _denied_page(payload))

        log.info("signed in %s <%s>", payload.get("username", "?"),
                 payload.get("email", "?"))
        self._redirect("/", headers=[
            self._set_cookie(make_session(payload),
                             max_age=SESSION_HOURS * 3600)])

    def _bootstrap(self, args) -> None:
        if not BOOTSTRAP_TOKEN:
            return self._send(404, page("Not found", "<h1>Not found</h1>"))
        given = (args.get("token") or [""])[0]
        if not hmac.compare_digest(BOOTSTRAP_TOKEN, given):
            log.warning("bootstrap attempt with a bad token from %s",
                        self.headers.get("X-Real-IP", "?"))
            return self._send(403, page("Rejected", "<h1>Rejected</h1>"))
        log.warning("BOOTSTRAP session issued -- this bypasses Discourse. "
                    "Unset IGVF_BOOTSTRAP_TOKEN once SSO is working.")
        self._redirect("/", headers=[
            self._set_cookie(make_session({"username": "bootstrap",
                                            "name": "Break-glass session",
                                            "admin": False}),
                             max_age=min(SESSION_HOURS, 12) * 3600)])

    def _logout(self) -> None:
        self._send(200, page("Signed out", f"""
<h1>Signed out</h1>
<p>Your IGVF Agent session on this browser has ended.</p>
<p><a class=btn href="/_auth/login">Sign in again</a></p>
<p class=muted>You are still signed in to the Genohub community itself — sign
   out there at <code>{_esc(DISCOURSE_URL)}</code> if you meant to do
   that too.</p>"""), headers=[self._set_cookie("", max_age=0)])


def main() -> int:
    logging.basicConfig(
        level=logging.INFO, stream=sys.stdout,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    missing = [n for n, v in (("IGVF_SSO_SECRET", SSO_SECRET),
                              ("IGVF_SESSION_SECRET", SESSION_SECRET)) if not v]
    if missing:
        # Refuse rather than start an open door. compose restarts it, the
        # operator sees the reason on every attempt.
        log.error("refusing to start: %s not set", ", ".join(missing))
        return 2
    if SSO_SECRET == SESSION_SECRET:
        log.error("refusing to start: IGVF_SSO_SECRET and "
                  "IGVF_SESSION_SECRET must differ -- reusing one secret lets "
                  "anyone who can mint an SSO payload mint a session cookie")
        return 2
    log.info("gate up on :%d  discourse=%s  group=%s  session=%dh",
             PORT, DISCOURSE_URL, APPROVAL_GROUP, SESSION_HOURS)
    ThreadingHTTPServer(("0.0.0.0", PORT), Handler).serve_forever()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
