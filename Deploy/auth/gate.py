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
import urllib.error
import urllib.parse
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

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
BRANDING = Path(os.environ.get("IGVF_BRANDING_DIR", "/app/branding"))
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


_PROVIDER_STATE: "dict[str, float | bool]" = {"at": 0.0, "ok": False}


def _provider_ready() -> bool:
    """Is Discourse actually configured to log people in?

    The provider endpoint 404s until `enable discourse connect provider` is
    switched on. Without this check the landing page would offer a Sign in
    button that dead-ends on somebody else's 404 with no explanation. Cached
    for five minutes so a page view is not a round trip, and any failure is
    treated as "not ready" -- an unnecessary warning is a far smaller problem
    than a button that silently does nothing.
    """
    now = _now()
    if now - float(_PROVIDER_STATE["at"]) < 300:
        return bool(_PROVIDER_STATE["ok"])
    # Only a definite 404 counts as "off". Anything else -- a timeout, a
    # network error, a bot challenge -- is inconclusive, and an inconclusive
    # probe must not put a scary banner on the front page of a working site.
    # (Discourse sits behind a CDN that answers 403 to the default
    # Python-urllib agent, which an earlier version read as "route exists" and
    # got exactly backwards. Hence the explicit agent and the 404-only rule.)
    ok = True
    try:
        req = urllib.request.Request(
            f"{DISCOURSE_URL}/session/sso_provider",
            headers={"User-Agent": "Mozilla/5.0 (compatible; IGVFagent-gate/1.0; "
                                    f"+{PUBLIC_URL})"})
        with urllib.request.urlopen(req, timeout=5) as resp:
            ok = resp.status != 404
    except urllib.error.HTTPError as exc:
        ok = exc.code != 404
    except Exception:
        ok = True
    _PROVIDER_STATE.update({"at": now, "ok": ok})
    return ok


def approved(payload: "dict[str, str]") -> bool:
    if _truthy(payload.get("admin")):
        return True
    wanted = APPROVAL_GROUP.lower()
    return any(g.lower() == wanted for g in _groups_for(payload))


# ------------------------------ pages ---------------------------------------
#
# The gate owns every page an unauthenticated visitor sees, so these are the
# first impression of the project for anyone arriving at the bare domain.
# Previously that first impression was an instant redirect into Discourse --
# no explanation of what the site is, and a raw 404 if the forum's SSO
# provider happened to be switched off. A landing page costs one HTTP response
# and replaces both problems.
#
# Styling is inline and dependency-free for the same reason the rest of this
# file is: the auth boundary does not get a build step or a CDN.

BRAND = "#38707f"

_SHELL = """<!doctype html><meta charset=utf-8>
<meta name=viewport content="width=device-width,initial-scale=1">
<title>{title}</title>
<link rel=icon href="/_auth/logo-mark.png">
<style>
 :root {{
   color-scheme: light dark;
   --brand: {brand};
   --bg: #f4f6f7; --card: #ffffff; --ink: #16202a; --muted: #5b6b78;
   --line: rgba(22,32,42,.10);
 }}
 @media (prefers-color-scheme: dark) {{
   :root {{ --bg: #0f1418; --card: #171e24; --ink: #e9edf0; --muted: #9bacb8;
            --line: rgba(255,255,255,.10); }}
 }}
 * {{ box-sizing: border-box; }}
 body {{
   margin: 0; min-height: 100vh; padding: 32px 16px;
   display: flex; align-items: center; justify-content: center;
   background:
     radial-gradient(1100px 520px at 50% -10%, color-mix(in srgb, var(--brand) 20%, transparent), transparent 70%),
     var(--bg);
   color: var(--ink);
   font: 16px/1.65 -apple-system, BlinkMacSystemFont, "Segoe UI", Inter,
         Roboto, "Helvetica Neue", Arial, sans-serif;
   -webkit-font-smoothing: antialiased;
 }}
 .card {{
   width: 100%; max-width: 40rem; background: var(--card);
   border: 1px solid var(--line); border-radius: 16px;
   padding: 40px 36px; text-align: center;
   box-shadow: 0 1px 2px rgba(0,0,0,.04), 0 12px 40px rgba(0,0,0,.07);
 }}
 .logo {{ width: 100%; max-width: 310px; height: auto; margin: 0 auto 6px; display: block; }}
 h1 {{ font-size: 1.3rem; line-height: 1.3; margin: 18px 0 10px; letter-spacing: -.01em; }}
 p {{ margin: 0 0 14px; color: var(--muted); }}
 p.lead {{ color: var(--ink); font-size: 1.02rem; }}
 .btn {{
   display: block; width: 100%; margin: 22px 0 10px; padding: .85em 1.2em;
   background: var(--brand); color: #fff; border: 0; border-radius: 9px;
   font: inherit; font-weight: 650; text-decoration: none; cursor: pointer;
 }}
 .btn:hover {{ filter: brightness(1.07); }}
 .btn.ghost {{
   background: transparent; color: var(--ink);
   border: 1px solid var(--line); font-weight: 550; margin-top: 0;
 }}
 .steps {{
   text-align: left; margin: 26px 0 0; padding: 18px 20px;
   border: 1px solid var(--line); border-radius: 11px;
   background: color-mix(in srgb, var(--brand) 5%, transparent);
 }}
 .steps ol {{ margin: 0; padding-left: 1.2em; }}
 .steps li {{ margin: 0 0 7px; color: var(--muted); }}
 .steps li:last-child {{ margin-bottom: 0; }}
 .steps b {{ color: var(--ink); font-weight: 600; }}
 .facts {{
   display: flex; flex-wrap: wrap; gap: 10px 26px; justify-content: center;
   margin: 22px 0 4px; padding: 0; list-style: none;
 }}
 .facts li {{ color: var(--muted); font-size: .88rem; }}
 .facts b {{ color: var(--ink); display: block; font-size: 1.15rem; font-weight: 650; }}
 code {{
   background: color-mix(in srgb, var(--ink) 9%, transparent);
   padding: .12em .42em; border-radius: 5px; font-size: .9em;
 }}
 .muted {{ font-size: .87rem; color: var(--muted); }}
 hr {{ border: 0; border-top: 1px solid var(--line); margin: 26px 0 20px; }}
 .warn {{
   text-align: left; margin: 20px 0 0; padding: 14px 16px; border-radius: 10px;
   border: 1px solid rgba(190,120,20,.35);
   background: rgba(220,150,40,.10); color: var(--ink); font-size: .92rem;
 }}
 @media (max-width: 460px) {{ .card {{ padding: 30px 22px; }} }}
</style>
<div class=card>{body}</div>
"""


def page(title: str, body: str) -> bytes:
    return _SHELL.format(title=title, body=body, brand=BRAND).encode()


def group_url() -> str:
    """The group page, where Discourse's own Request Membership button lives.

    Sending people here rather than telling them to "ask an administrator" is
    the difference between a request that gets filed and notifies the group's
    owners, and an email somebody has to remember to act on.
    """
    return f"{DISCOURSE_URL}/g/{APPROVAL_GROUP}"


def _logo(alt: str = "IGVF Agent") -> str:
    return f'<img class=logo src="/_auth/logo.png" alt="{alt}">'


def _welcome_page(provider_ready: bool) -> bytes:
    """The first thing anyone sees at the bare domain."""
    # Placed ABOVE the buttons, not below: a warning that a button will not
    # work is only useful before it is clicked.
    trouble = "" if provider_ready else f"""
<div class=warn><b>Sign-in is not available yet.</b> The Genohub Community has
  not finished being set up as this site's login provider, so the sign-in
  button below will not work. An administrator needs to enable
  <code>discourse connect provider</code> on
  {_esc(DISCOURSE_URL)}.</div>"""
    return page("IGVF Agent", f"""
{_logo()}
<p class=lead>An auditable AI agent for discovering, retrieving and analysing
   data across the IGVF ecosystem — Portal, Catalog and Knowledge Graph —
   alongside ENCODE and related public resources.</p>
<ul class=facts>
  <li><b>85</b> skills</li>
  <li><b>258</b> typed tools</li>
  <li><b>Local</b> execution</li>
</ul>
{trouble}
<a class=btn href="/_auth/login">Sign in with the Genohub Community</a>
<a class="btn ghost" href="{_esc(group_url())}">Request access</a>
<div class=steps>
  <ol>
    <li><b>Sign up</b> at the
        <a href="{_esc(DISCOURSE_URL)}/signup">Genohub Community</a>, if you
        have not already.</li>
    <li><b>Request access</b> on the
        <a href="{_esc(group_url())}">IGVF Agent Users</a> group page and say
        briefly who you are — having a community account does not grant access
        on its own.</li>
    <li><b>Sign in here</b> once an administrator has approved you. Nothing
        else is needed on your side.</li>
  </ol>
</div>
<hr>
<p class=muted>Analyses run on shared hardware, which is why access is
   approved rather than open. Prefer to run it yourself? IGVF Agent is open
   source and installs locally with your own API keys, or entirely free with
   local open-weight models.</p>""")


def _denied_page(payload: "dict[str, str]") -> bytes:
    who = payload.get("username") or payload.get("email") or "your account"
    return page("Access pending — IGVF Agent", f"""
{_logo()}
<h1>Not yet approved</h1>
<p>You are signed in to the Genohub Community as
   <code>{_esc(who)}</code>, but that account is not in the
   <code>{_esc(APPROVAL_GROUP)}</code> group yet.</p>
<p>Request access on the group page below and say briefly who you are. An
   administrator reviews it; once you are approved, reload this page.</p>
<a class=btn href="{_esc(group_url())}">Request access</a>
<a class="btn ghost" href="/_auth/login">I have been approved — try again</a>
<hr>
<p class=muted>Signing up to the community and being approved for the agent
   are two separate steps, on purpose: analyses run on shared hardware.</p>""")


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
        if path in ("/_auth", "/_auth/", "/_auth/welcome"):
            return self._welcome()
        if path in ("/_auth/logo.png", "/_auth/logo-mark.png"):
            return self._asset(path.rsplit("/", 1)[-1])
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
        if path == "/_auth/echo":
            return self._echo()
        self._send(404, page("Not found", "<h1>Not found</h1>"))

    do_HEAD = do_GET

    def _welcome(self) -> None:
        # If a session is already valid there is nothing to welcome anyone to;
        # send them straight into the app rather than showing a sign-in page
        # to somebody who is signed in.
        if read_session(self._cookies().get(COOKIE, "")):
            return self._redirect("/")
        self._send(200, _welcome_page(_provider_ready()))

    def _asset(self, name: str) -> None:
        path = BRANDING / name
        try:
            blob = path.read_bytes()
        except OSError:
            return self._send(404, b"", ctype="text/plain")
        self._send(200, blob, ctype="image/png",
                   headers=[("Cache-Control", "public, max-age=86400")])

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

    def _echo(self) -> None:
        """Reflect the X-IGVF-* headers this request arrived with.

        A diagnostic for one specific silent failure: if auth_request_set or
        proxy_set_header were wrong, the app would receive no identity, read
        that as "this deployment has no accounts", and show every user
        everyone's private work. Nothing here is privileged -- it echoes what
        the caller sent -- so it is safe to leave in place, and pointing a
        location at it through the same proxy_set_header lines the app gets is
        the only way to test that chain end to end.
        """
        body = "".join(
            f"{h}: {self.headers.get(h, '')}\n"
            for h in ("X-IGVF-User", "X-IGVF-Email", "X-IGVF-Name",
                      "X-IGVF-Admin"))
        self._send(200, body.encode(), ctype="text/plain; charset=utf-8")

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
        self._send(200, page("Signed out — IGVF Agent", f"""
{_logo()}
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
