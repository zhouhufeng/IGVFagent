#!/usr/bin/env python3
"""
igvf-sub -- preflight, audit and submission helper for the IGVF data portal.

Zero dependencies (Python 3.9+ stdlib only), so it runs on a laptop, a cluster
login node or inside a container without an environment to build first.

Why this exists
---------------
The IGVF submission loop is slow: you submit, the DACC runs audits, and a month
later a meeting tells you what was missing. The same handful of things are
missing every time. This tool encodes those recurring findings as checks you can
run in seconds, before the meeting, against the live portal.

Credentials
-----------
    export IGVF_API_KEY=...        # Access Key ID from the portal
    export IGVF_SECRET_KEY=...     # Access Key Secret
    export IGVF_MODE=prod          # or: sandbox, staging
Optionally IGVF_LAB and IGVF_AWARD, matching igvf_utils conventions.

Read access to your own unreleased objects REQUIRES the keys; released objects
are public.
"""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import urllib.error
import urllib.parse
import urllib.request

__version__ = "1.0.0"

# Portal hostnames, matching igvf_utils' IGVF_{PROD,SANDBOX,STAGING}_MODE.
MODES = {
    "prod": "https://api.data.igvf.org",
    "sandbox": "https://api.sandbox.igvf.org",
    "staging": "https://api.staging.igvf.org",
}

# api.sandbox.igvf.org now answers HTTP 410 with "is deprecated. Please update
# your client to use api.staging.igvf.org." Older submission docs and meeting
# notes still say sandbox, so redirect the intent rather than letting the call
# fail with a confusing 410.
DEPRECATED_MODES = {"sandbox": "staging"}
UI_HOSTS = {
    "prod": "https://data.igvf.org",
    "sandbox": "https://sandbox.igvf.org",
    "staging": "https://staging.igvf.org",
}

# Statuses that make an object unusable as an input. Discovered the hard way:
# a revoked input silently blocks release of everything downstream.
DEAD_STATUSES = {"revoked", "deleted", "replaced", "archived"}

# The generic GRCh38 reference the DACC repeatedly told submitters to use.
# NOTE: as of 2026-09 this is status=archived on prod -- see RULE ref-archived.
LEGACY_GRCH38_REF = "IGVFFI8743BGYR"


# ── submission dependency order (spec §2.1) ────────────────────────────────
# The DACC meetings kept surfacing the same failure: an object posted before
# the thing it points at, which the Portal rejects with a reference error that
# names the missing target rather than the ordering mistake. This encodes the
# order once so a plan can be checked before anything is sent.
SUBMISSION_ORDER = [
    ("software", "the tool itself"),
    ("software_version", "MUST point at a TAGGED GitHub release, linked by "
                          "alias (lab-name:tool)"),
    ("workflow", "ties the software versions together"),
    ("analysis_step", "a step within the workflow"),
    ("analysis_step_version", "what a FILE links to, via its "
                               "analysis_step_version property"),
    ("prediction_set", "or curated_set / analysis_set — the file set"),
    ("tabular_file", "the data itself; gzipped, with derived_from and "
                      "reference_files"),
    ("document", "the file format specification, linked from the file's "
                  "file_format_specifications (NOT documents)"),
    ("curated_set_external", "external-source tabular file under a curated "
                              "set: Name | URL | description"),
]

# Properties the Portal will not accept an object without, beyond whatever the
# live schema says. These come from the meeting record rather than the schema,
# because several are schema-optional and audit-required -- which is exactly
# the combination that gets a submission bounced a month later.
REQUIRED_BY_TYPE = {
    "prediction_set": ["lab", "award", "file_set_type", "description",
                        "samples", "input_file_sets"],
    "curated_set": ["lab", "award", "file_set_type", "description"],
    "analysis_set": ["lab", "award", "file_set_type", "description"],
    "tabular_file": ["lab", "award", "content_type", "file_format",
                      "derived_from", "reference_files",
                      "file_format_specifications", "analysis_step_version"],
    "document": ["lab", "award", "document_type"],
    "software_version": ["lab", "award", "version", "downloaded_url"],
    "analysis_step_version": ["lab", "award", "analysis_step"],
}


def plan_submission(objects: "list[dict]") -> dict:
    """Order a set of draft objects for submission, and say what is missing.

    `objects` is a list of dicts each carrying at least `type`, plus whatever
    properties are drafted so far. Nothing is sent; this is the step that is
    supposed to happen before anything is.
    """
    order = {t: i for i, (t, _) in enumerate(SUBMISSION_ORDER)}
    known, unknown = [], []
    for o in objects:
        t = str(o.get("type") or o.get("@type") or "").strip()
        (known if t in order else unknown).append(o)
    known.sort(key=lambda o: order[str(o.get("type"))])

    steps = []
    for i, o in enumerate(known, 1):
        t = str(o.get("type"))
        need = REQUIRED_BY_TYPE.get(t, ["lab", "award"])
        missing = [k for k in need if not o.get(k)]
        steps.append({
            "position": i, "type": t,
            "alias": o.get("aliases") or o.get("alias"),
            "missing_required": missing,
            "ready": not missing,
            "why_here": dict(SUBMISSION_ORDER)[t],
        })
    return {
        "steps": steps,
        "n_ready": sum(1 for s in steps if s["ready"]),
        "n_blocked": sum(1 for s in steps if not s["ready"]),
        "unrecognised_types": [o.get("type") for o in unknown],
        "order_reference": [t for t, _ in SUBMISSION_ORDER],
        "note": ("Nothing was sent. Post in the order above: the Portal "
                  "rejects an object whose reference does not exist yet, and "
                  "reports it as a missing target rather than as an ordering "
                  "mistake."),
    }


def _portal_key_pair():
    """The portal key pair, from either naming convention.

    This tool was written standalone against IGVF_API_KEY / IGVF_SECRET_KEY,
    the names igvf_utils uses. IGVFagent stores the same pair as
    IGVF_ACCESS_KEY / IGVF_SECRET_ACCESS_KEY -- that is what
    Deploy/make-live.sh upserts into .env.prod and what _credentials.py
    parses from Docs/Secret/IGVFportalAPI.txt. Reading only the first pair
    would report "no credentials" on a deployment that has them, which is
    the same failure as probing `bean --version` on a working install.

    Order: explicit igvf_utils names first, so a submitter who followed the
    IGVF-Submit README keeps working; then IGVFagent's own names; then
    _credentials.py, which also reads the on-disk key file.
    """
    key = os.environ.get("IGVF_API_KEY")
    secret = os.environ.get("IGVF_SECRET_KEY")
    if key and secret:
        return key, secret
    key = os.environ.get("IGVF_ACCESS_KEY")
    secret = os.environ.get("IGVF_SECRET_ACCESS_KEY")
    if key and secret:
        return key, secret
    try:                                    # optional: absent when standalone
        from . import _credentials          # type: ignore
    except ImportError:
        try:
            import _credentials             # type: ignore
        except ImportError:
            return None, None
    try:
        pair = _credentials.portal_credentials()
    except Exception:
        return None, None
    return pair if pair else (None, None)


class _StripAuthRedirect(urllib.request.HTTPRedirectHandler):
    """Drop Authorization when a redirect crosses hosts.

    File downloads redirect from the portal to a presigned S3 URL. S3 rejects a
    request carrying both a presigned signature and an Authorization header with
    HTTP 400, so the credential must not follow the redirect.
    """

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is not None:
            old_host = urllib.parse.urlsplit(req.full_url).netloc
            if old_host != urllib.parse.urlsplit(newurl).netloc:
                for store in (new.headers, new.unredirected_hdrs):
                    for key in [k for k in store if k.lower() == "authorization"]:
                        del store[key]
        return new


class PortalError(RuntimeError):
    pass


class Portal:
    """Thin authenticated client over the IGVF REST API."""

    def __init__(self, mode: str = None, api_key: str = None, secret_key: str = None):
        self.mode = mode or os.environ.get("IGVF_MODE", "prod")
        if self.mode in DEPRECATED_MODES:
            replacement = DEPRECATED_MODES[self.mode]
            print(f"note: '{self.mode}' is deprecated by the portal; using "
                  f"'{replacement}' instead. Rehearsal keys are per-host, so a "
                  f"{self.mode} key will not work on {replacement}.", file=sys.stderr)
            self.mode = replacement
        if self.mode in MODES:
            self.base = MODES[self.mode]
        elif self.mode.startswith("http"):
            self.base = self.mode.rstrip("/")
        else:
            self.base = f"https://api.{self.mode.lstrip('.')}".rstrip("/")
        # Fall back to the ambient key pair ONLY when the caller passed
        # neither. Substituting into a half-supplied pair would silently
        # replace an explicitly chosen credential with whatever is in the
        # environment -- and portal keys are per-host, so that is how a prod
        # key ends up being sent to staging.
        if api_key is None and secret_key is None:
            self.api_key, self.secret_key = _portal_key_pair()
        else:
            self.api_key, self.secret_key = api_key, secret_key

    @property
    def authenticated(self) -> bool:
        return bool(self.api_key and self.secret_key)

    @property
    def ui_base(self) -> str:
        return UI_HOSTS.get(self.mode, self.base.replace("api.", ""))

    def _request(self, method: str, path: str, params: dict = None, payload=None,
                 anon: bool = False):
        url = self.base.rstrip("/") + "/" + path.lstrip("/")
        if params:
            url += ("&" if "?" in url else "?") + urllib.parse.urlencode(params)
        body = None
        req_headers = {"Accept": "application/json", "User-Agent": f"igvf-sub/{__version__}"}
        if payload is not None:
            body = json.dumps(payload).encode()
            req_headers["Content-Type"] = "application/json"
        if self.authenticated and not anon:
            token = base64.b64encode(f"{self.api_key}:{self.secret_key}".encode()).decode()
            req_headers["Authorization"] = f"Basic {token}"
        req = urllib.request.Request(url, data=body, headers=req_headers, method=method)
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                raw = resp.read()
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", "replace")[:500]
            hint = ""
            if exc.code == 403:
                hint = ("\n  Hint: unreleased objects need IGVF_API_KEY and IGVF_SECRET_KEY, "
                        "and the keys must belong to a lab with access.")
            elif exc.code == 422:
                hint = "\n  Hint: 422 is a schema validation failure -- the payload violates the profile."
            raise PortalError(f"HTTP {exc.code} on {method} {url}\n  {detail}{hint}") from None
        except urllib.error.URLError as exc:
            raise PortalError(f"cannot reach {url}: {exc.reason}") from None
        return json.loads(raw.decode()) if raw else {}

    def get(self, ident: str, frame: str = "object") -> dict:
        params = {"frame": frame} if frame else None
        return self._request("GET", f"/{ident.strip('/')}/", params=params)

    def patch(self, ident: str, payload: dict) -> dict:
        return self._request("PATCH", f"/{ident.strip('/')}/", payload=payload)

    def post(self, profile: str, payload: dict) -> dict:
        return self._request("POST", f"/{profile}/", payload=payload)

    def profile(self, profile: str) -> dict:
        return self._request("GET", f"/profiles/{profile}")

    def search(self, item_type: str, limit="all", **filters) -> list:
        params = {"type": item_type, "limit": limit, "format": "json"}
        params.update(filters)
        res = self._request("GET", "/search/", params=params)
        return res.get("@graph", [])

    def download(self, obj: dict, dest: str) -> int:
        """Download a file object's payload via its href."""
        href = obj.get("href")
        if not href:
            raise PortalError(f"{obj.get('accession')} has no href (nothing uploaded?)")
        url = self.base.rstrip("/") + href
        req = urllib.request.Request(url, headers={"User-Agent": f"igvf-sub/{__version__}"})
        if self.authenticated:
            token = base64.b64encode(f"{self.api_key}:{self.secret_key}".encode()).decode()
            req.add_header("Authorization", f"Basic {token}")
        opener = urllib.request.build_opener(_StripAuthRedirect())
        with opener.open(req, timeout=600) as resp, open(dest, "wb") as out:
            total = 0
            while True:
                chunk = resp.read(1 << 16)
                if not chunk:
                    break
                out.write(chunk)
                total += len(chunk)
        return total


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

USE_COLOR = sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


def _c(text: str, code: str) -> str:
    return f"\033[{code}m{text}\033[0m" if USE_COLOR else text


def fail(msg: str) -> str:
    return _c("FAIL", "31;1") + " " + msg


def warn(msg: str) -> str:
    return _c("WARN", "33;1") + " " + msg


def ok(msg: str) -> str:
    return _c("OK  ", "32;1") + " " + msg


def info(msg: str) -> str:
    return _c("--  ", "2") + " " + msg


def die(msg: str, code: int = 1):
    print(_c("error:", "31;1") + " " + msg, file=sys.stderr)
    sys.exit(code)


def acc_of(ref) -> str:
    """'/tabular-files/IGVFFI0856UJHP/' -> 'IGVFFI0856UJHP'."""
    if isinstance(ref, dict):
        ref = ref.get("@id") or ref.get("accession") or ""
    return str(ref).strip("/").split("/")[-1]


def type_of(obj: dict) -> str:
    types = obj.get("@type") or []
    return types[0] if types else "?"


# ---------------------------------------------------------------------------
# Audit rules
#
# Each rule is a recurring DACC finding from the monthly IGVF submission
# meetings, turned into something machine-checkable. `since` records where the
# requirement came from, so a failing check can be traced to the conversation
# that produced it rather than looking like an arbitrary opinion.
# ---------------------------------------------------------------------------

FILE_TYPES = {"TabularFile", "ReferenceFile", "SequenceFile", "AlignmentFile",
              "SignalFile", "MatrixFile", "ImageFile", "ModelFile", "ConfigurationFile"}
FILESET_TYPES = {"PredictionSet", "CuratedSet", "AnalysisSet", "MeasurementSet",
                 "AuxiliarySet", "ConstructLibrarySet", "ModelSet"}


class Finding:
    def __init__(self, rule, level, message, fix=None, since=None):
        self.rule, self.level, self.message = rule, level, message
        self.fix, self.since = fix, since

    def render(self) -> str:
        head = {"fail": fail, "warn": warn, "ok": ok, "info": info}[self.level]
        out = [head(f"[{self.rule}] {self.message}")]
        if self.fix:
            out.append(f"       fix: {self.fix}")
        if self.since:
            out.append(f"       source: {self.since}")
        return "\n".join(out)

    def as_dict(self):
        return {"rule": self.rule, "level": self.level, "message": self.message,
                "fix": self.fix, "source": self.since}


def audit_object(portal: Portal, obj: dict, deep: bool = True) -> list:
    """Run the recurring-findings checklist against one portal object."""
    findings = []
    t = type_of(obj)
    acc = obj.get("accession") or acc_of(obj.get("@id", ""))

    # -- description ---------------------------------------------------------
    # Missing description was an open to-do for four consecutive meetings.
    if t in FILESET_TYPES:
        desc = (obj.get("description") or "").strip()
        if not desc:
            findings.append(Finding(
                "description", "fail",
                f"{t} {acc} has no description.",
                "1-2 sentences saying what was predicted and how. PATCH `description`.",
                "Oct 10 2025 / Jan 16 / Feb 20 2026 meetings"))
        elif len(desc) < 25:
            findings.append(Finding(
                "description", "warn",
                f"description is very short ({len(desc)} chars): {desc!r}",
                "The DACC asked for 1-2 full sentences.",
                "Feb 20 2026 meeting"))
        else:
            findings.append(Finding("description", "ok", "description present."))

    # -- input_file_sets on a prediction set --------------------------------
    if t == "PredictionSet":
        if not obj.get("input_file_sets"):
            findings.append(Finding(
                "input-file-sets", "fail",
                f"PredictionSet {acc} has no input_file_sets.",
                "Collect the accessions of every file set you consumed and PATCH "
                "`input_file_sets`. Try: igvf-sub fileset-files <curated-set> --field accession",
                "Feb 20 2026 meeting ('Unexpected input file sets' audit)"))
        else:
            findings.append(Finding("input-file-sets", "ok",
                                    f"input_file_sets: {len(obj['input_file_sets'])} linked."))
        if not obj.get("samples") and not obj.get("donors"):
            findings.append(Finding(
                "sample", "warn",
                "PredictionSet has neither samples nor donors.",
                "Each prediction set needs a sample, e.g. IGVFSM0002XTWQ (HepG2) "
                "or IGVFSM4381OIBT (K562).",
                "Aug 21 2026 meeting"))

    # -- file-level requirements --------------------------------------------
    if t in FILE_TYPES:
        if not obj.get("derived_from"):
            findings.append(Finding(
                "derived-from", "fail",
                f"{t} {acc} has no derived_from.",
                "List every input file's accession. If you filtered the inputs, that "
                "filtering must be described in the file format specification document.",
                "Feb 20 / May 15 / Jun 26 2026 meetings"))
        else:
            findings.append(Finding("derived-from", "ok",
                                    f"derived_from: {len(obj['derived_from'])} files."))

        # A missing genomic reference is a hard audit error on tabular files.
        refs = obj.get("reference_files") or []
        if t == "TabularFile" and not refs:
            findings.append(Finding(
                "reference-files", "fail",
                f"TabularFile {acc} has no reference_files -- this returns a portal audit error.",
                f"PATCH `reference_files` with a generic GRCh38 reference.",
                "Apr 10 2026 meeting"))
        elif refs:
            findings.append(Finding("reference-files", "ok",
                                    f"reference_files: {', '.join(acc_of(r) for r in refs)}"))
            if any(acc_of(r) == LEGACY_GRCH38_REF for r in refs):
                findings.append(Finding(
                    "ref-archived", "warn",
                    f"reference_files points at {LEGACY_GRCH38_REF}, which is now status=archived "
                    "on prod.",
                    "The DACC recommended this reference through Jun 2026, but it has since "
                    "been archived. Check what a currently-released sibling file uses and "
                    "consider re-pointing.",
                    "Verified live against the portal, Sep 2026"))

        if not obj.get("file_format_specifications"):
            findings.append(Finding(
                "format-spec", "fail",
                f"{t} {acc} has no file_format_specifications.",
                "Submit a Document with document_type='file format specification' stating "
                "0- vs 1-based coordinates, whether end is inclusive, and each column's "
                "meaning. Link it via `file_format_specifications` -- NOT via `documents`.",
                "Mar 20 / May 15 / Jun 26 / Aug 21 2026 meetings"))
        else:
            findings.append(Finding("format-spec", "ok", "file_format_specifications linked."))

        if obj.get("documents") and not obj.get("file_format_specifications"):
            findings.append(Finding(
                "spec-wrong-property", "fail",
                "A document is linked via `documents` but `file_format_specifications` is empty.",
                "The DACC was explicit: link the format spec via file_format_specifications, "
                "not documents.",
                "Jun 26 2026 meeting"))

        if not obj.get("analysis_step_version"):
            findings.append(Finding(
                "analysis-step-version", "warn",
                f"{t} {acc} has no analysis_step_version.",
                "Once the workflow object exists, link it via the file's "
                "`analysis_step_version`.",
                "Jun 26 / Aug 21 2026 meetings"))
        else:
            findings.append(Finding("analysis-step-version", "ok",
                                    "analysis_step_version linked."))

        # Gzip. Checked via the submitted filename, since the portal stores the
        # object not the bytes.
        name = obj.get("submitted_file_name") or ""
        href = obj.get("href") or ""
        if name and not (name.endswith(".gz") or href.endswith(".gz")):
            findings.append(Finding(
                "gzip", "fail",
                f"submitted_file_name does not look gzipped: {name}",
                "gzip the file and re-upload. This was raised at every meeting.",
                "Jan-Aug 2026 meetings, repeatedly"))

        # BED files must be bed3 and 0-based for element references.
        if obj.get("file_format") == "bed":
            if not obj.get("file_format_type"):
                findings.append(Finding(
                    "bed-subtype", "fail",
                    "file_format=bed but file_format_type is unset.",
                    "Set file_format_type, e.g. 'bed3' for a plain coordinate list.",
                    "Aug 21 2026 meeting"))
            if obj.get("content_type") == "element references" \
                    and obj.get("file_format_type") not in (None, "bed3"):
                findings.append(Finding(
                    "bed3", "warn",
                    f"content_type='element references' with file_format_type="
                    f"{obj.get('file_format_type')!r}; the DACC asked for bed3.",
                    "A coordinate list with no scores should be 0-based bed3.",
                    "Aug 21 2026 meeting"))

    # -- dead inputs ---------------------------------------------------------
    if deep:
        findings.extend(check_inputs(portal, obj))

    # -- release readiness ---------------------------------------------------
    status = obj.get("status")
    if status in ("in progress", "submitted"):
        findings.append(Finding(
            "status", "warn",
            f"status is '{status}' -- not released.",
            "The Sep 1 2026 deadline for Catalog v2.0 inclusion has passed; confirm with "
            "the DACC whether a late release still makes v2.0.",
            "Aug 21 2026 meeting (deadline Sep 1 2026)"))
    elif status:
        findings.append(Finding("status", "info", f"status: {status}"))

    return findings


def check_inputs(portal: Portal, obj: dict) -> list:
    """Flag inputs that are revoked/archived/replaced, and name replacements.

    This is the check that would have caught the Aug 2026 problem a month early:
    an input file was revoked for a reference-allele mismatch, which silently
    blocks release of everything derived from it.
    """
    findings = []
    linked = []
    for prop in ("derived_from", "input_file_sets", "reference_files", "files"):
        for ref in obj.get(prop) or []:
            linked.append((prop, acc_of(ref)))
    for prop, ident in linked:
        try:
            dep = portal.get(ident)
        except PortalError:
            findings.append(Finding("input-unreadable", "warn",
                                    f"{prop} -> {ident}: cannot read (permissions?)"))
            continue
        st = dep.get("status")
        if st in DEAD_STATUSES:
            sup = [acc_of(s) for s in dep.get("superseded_by") or []]
            detail = dep.get("revoke_detail") or ""
            msg = f"{prop} -> {ident} is status={st}."
            if detail:
                msg += f' Portal says: "{detail}"'
            fixtext = (f"Repoint {prop} to {', '.join(sup)}."
                       if sup else "Find the current replacement and repoint.")
            findings.append(Finding("dead-input", "fail", msg, fixtext,
                                    "Aug 21 2026 meeting"))
    if not any(f.rule == "dead-input" for f in findings) and linked:
        findings.append(Finding("dead-input", "ok",
                                f"all {len(linked)} linked objects are in a usable status."))
    return findings


# ---------------------------------------------------------------------------
# Local file validation
# ---------------------------------------------------------------------------

def is_gzip(path: str) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_maybe_gz(path: str):
    return gzip.open(path, "rt", newline="") if is_gzip(path) else open(path, "rt", newline="")


def md5_of(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def content_md5_of(path: str) -> str:
    """md5 of the DECOMPRESSED bytes -- the portal's content_md5sum."""
    h = hashlib.md5()
    opener = gzip.open if is_gzip(path) else open
    with opener(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sniff_delim(sample: str) -> str:
    counts = {d: sample.count(d) for d in ("\t", ",")}
    return max(counts, key=counts.get) if max(counts.values()) else "\t"


def validate_local(path: str, expect_bed3: bool = False) -> list:
    """Checks that fail a submission before the portal ever sees it."""
    findings = []
    if not os.path.exists(path):
        return [Finding("exists", "fail", f"no such file: {path}")]

    if is_gzip(path):
        findings.append(Finding("gzip", "ok", "file is gzipped."))
    else:
        findings.append(Finding(
            "gzip", "fail", f"{os.path.basename(path)} is not gzipped.",
            "Run: gzip -n " + path,
            "raised at every 2026 meeting"))

    findings.append(Finding("md5", "info", f"md5sum        = {md5_of(path)}"))
    findings.append(Finding("md5", "info", f"content_md5sum = {content_md5_of(path)}"))

    with open_maybe_gz(path) as fh:
        lines, comments = [], []
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            (comments if line.startswith("#") else lines).append(line)
            if len(lines) > 5000:
                break
    if not lines:
        return findings + [Finding("empty", "fail", "no data rows found.")]

    delim = sniff_delim(lines[0])
    ncol = len(lines[0].split(delim))
    findings.append(Finding("shape", "info",
                            f"{len(comments)} comment line(s), delimiter="
                            f"{'TAB' if delim == chr(9) else 'COMMA'}, {ncol} columns"))

    ragged = [i for i, l in enumerate(lines[1:], 2) if len(l.split(delim)) != ncol]
    if ragged:
        findings.append(Finding(
            "ragged", "fail",
            f"{len(ragged)} row(s) have a different column count than the first row "
            f"(e.g. line {ragged[0]}).",
            "Ragged rows fail portal validation."))
    else:
        findings.append(Finding("ragged", "ok", "all sampled rows have a consistent width."))

    if expect_bed3:
        findings.extend(_validate_bed3(lines, delim))
    return findings


def _validate_bed3(lines, delim) -> list:
    findings = []
    bad_type, neg, inverted, one_based_suspect = [], [], [], 0
    for i, line in enumerate(lines, 1):
        parts = line.split(delim)
        if len(parts) < 3:
            bad_type.append(i)
            continue
        try:
            start, end = int(parts[1]), int(parts[2])
        except ValueError:
            bad_type.append(i)
            continue
        if start < 0 or end < 0:
            neg.append(i)
        if end <= start:
            inverted.append(i)
        if end - start == 0:
            one_based_suspect += 1
    if bad_type:
        findings.append(Finding("bed-cols", "fail",
                                f"{len(bad_type)} row(s) lack 3 columns with integer start/end "
                                f"(e.g. line {bad_type[0]}).",
                                "bed3 is: chrom, start, end -- tab separated."))
    if neg:
        findings.append(Finding("bed-negative", "fail",
                                f"{len(neg)} row(s) have a negative coordinate."))
    if inverted:
        findings.append(Finding(
            "bed-halfopen", "fail",
            f"{len(inverted)} row(s) have end <= start (e.g. line {inverted[0]}).",
            "BED is 0-based half-open: a 1bp feature at position p is start=p-1, end=p, "
            "so end must be strictly greater than start.",
            "Mar 20 2026 meeting (0-based BED required)"))
    if one_based_suspect:
        findings.append(Finding(
            "bed-1based", "warn",
            f"{one_based_suspect} row(s) have end == start, which suggests 1-based "
            "coordinates written into a 0-based format.",
            "Convert to 0-based half-open before submitting."))
    if not (bad_type or neg or inverted or one_based_suspect):
        findings.append(Finding("bed3", "ok",
                                "coordinates look like valid 0-based half-open BED."))
    return findings


# ---------------------------------------------------------------------------
# Variant cross-check
# ---------------------------------------------------------------------------

def parse_variants(path, chrom_col=None, pos_col=None, id_col=None):
    """Map (chrom, pos) -> set of (ref, alt) allele pairs seen at that position.

    Alleles matter here: a reference-allele mismatch is corrected in place, so
    the position is unchanged while the alleles change. Keying on position alone
    would make exactly the dangerous case invisible.

    Handles the 'chr_pos_hg38_ref_alt' target_id used by the Sherwood base
    editing files and the explicit VariantChr/VariantStart/EffectAllele columns
    used by the FAVOR/CAMP prediction files.
    """
    with open_maybe_gz(path) as fh:
        rows = [l.rstrip("\r\n") for l in fh if l.strip() and not l.startswith("#")]
    if not rows:
        return {}, {"chrom_col": None, "pos_col": None, "id_col": None,
                    "fields": [], "n": 0, "with_alleles": 0}
    delim = sniff_delim(rows[0])
    reader = csv.DictReader(rows, delimiter=delim)
    fields = reader.fieldnames or []

    def pick(cands):
        for c in cands:
            for f in fields:
                if f.lower() == c:
                    return f
        return None

    chrom_col = chrom_col or pick(["variantchr", "chrom", "chr", "chromosome", "#chrom"])
    pos_col = pos_col or pick(["variantstart", "pos", "position", "start"])
    id_col = id_col or pick(["target_id", "variant_id", "vid"])
    ref_col = pick(["otherallele", "ref", "ref_allele", "reference_allele"])
    alt_col = pick(["effectallele", "alt", "alt_allele", "effect_allele"])

    posmap, with_alleles = {}, 0
    for row in reader:
        chrom = pos = ref = alt = None
        if chrom_col and pos_col and row.get(chrom_col) and row.get(pos_col) not in (None, ""):
            chrom = str(row[chrom_col]).replace("chr", "")
            try:
                pos = int(str(row[pos_col]).strip())
            except ValueError:
                pos = None
            if ref_col and alt_col:
                ref, alt = (row.get(ref_col) or "").upper(), (row.get(alt_col) or "").upper()
        elif id_col and row.get(id_col):
            parts = str(row[id_col]).split("_")
            if len(parts) >= 2:
                chrom = parts[0].replace("chr", "")
                try:
                    pos = int(parts[1])
                except ValueError:
                    pos = None
                # chrom_pos_build_ref_alt
                if len(parts) >= 5:
                    ref, alt = parts[-2].upper(), parts[-1].upper()
        if chrom is None or pos is None:
            continue
        key = (chrom, pos)
        posmap.setdefault(key, set())
        if ref and alt:
            posmap[key].add((ref, alt))
            with_alleles += 1
    return posmap, {"chrom_col": chrom_col, "pos_col": pos_col, "id_col": id_col,
                    "ref_col": ref_col, "alt_col": alt_col, "fields": fields,
                    "n": len(posmap), "with_alleles": with_alleles}


def normalise_pair(pair):
    """Treat a ref/alt swap as the same underlying edit for comparison purposes."""
    ref, alt = pair
    return tuple(sorted((ref, alt)))


def revised_positions(old_map, new_map, ignore_swaps=False):
    """Positions whose variant definition changed between two input revisions.

    Returns (dropped, allele_changed). `dropped` disappeared entirely;
    `allele_changed` kept the coordinate but changed alleles.
    """
    old_pos, new_pos = set(old_map), set(new_map)
    dropped = old_pos - new_pos
    allele_changed = set()
    for key in old_pos & new_pos:
        a, b = old_map[key], new_map[key]
        if not a or not b:
            continue                      # no allele info -- cannot compare
        if ignore_swaps:
            a = {normalise_pair(p) for p in a}
            b = {normalise_pair(p) for p in b}
        if a != b:
            allele_changed.add(key)
    return dropped, allele_changed


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def report(findings, as_json=False) -> int:
    """Print findings; return an exit code (1 if anything failed)."""
    if as_json:
        print(json.dumps([f.as_dict() for f in findings], indent=2))
    else:
        for f in findings:
            print(f.render())
        n_fail = sum(1 for f in findings if f.level == "fail")
        n_warn = sum(1 for f in findings if f.level == "warn")
        print()
        verdict = (_c("BLOCKED", "31;1") if n_fail
                   else _c("READY (with warnings)", "33;1") if n_warn
                   else _c("READY", "32;1"))
        print(f"{verdict}  {n_fail} failure(s), {n_warn} warning(s)")
    return 1 if any(f.level == "fail" for f in findings) else 0


def cmd_doctor(portal: Portal, args) -> int:
    print(f"igvf-sub {__version__}")
    print(f"mode:  {portal.mode}  ->  {portal.base}")
    print(f"keys:  {'present' if portal.authenticated else _c('MISSING', '31;1')}"
          f" (IGVF_API_KEY / IGVF_SECRET_KEY)")
    for var in ("IGVF_LAB", "IGVF_AWARD"):
        val = os.environ.get(var)
        print(f"{var.lower():6} {val if val else '(unset)'}")
    print()
    try:
        root = portal._request("GET", "/profiles", anon=True)
        print(ok(f"reachable -- {len(root)} profiles published"))
    except PortalError as exc:
        print(fail(f"host unreachable: {exc}"))
        return 1
    if portal.authenticated:
        try:
            me = portal._request("GET", "/session-properties")
            user = me.get("user", {})
            print(ok(f"authenticated as {user.get('title', '?')}"))
            labs = user.get("submits_for") or []
            if labs:
                print(info("submits_for: " + ", ".join(acc_of(l) for l in labs)))
            else:
                print(warn("this key can read but is not registered to submit for any lab."))
        except PortalError as exc:
            print(fail(f"credentials rejected on {portal.mode}: {exc}"))
            print(info("Keys are per-host: a prod key is not valid on staging, and "
                       "vice versa. Issue a key from this host's own profile page."))
            return 1
    else:
        print(warn("no keys -- only released (public) objects will be readable."))
    return 0


def cmd_get(portal: Portal, args) -> int:
    obj = portal.get(args.accession, frame=args.frame)
    if args.field:
        for f in args.field:
            print(f"{f}: {json.dumps(obj.get(f))}")
    else:
        print(json.dumps(obj, indent=2, sort_keys=True))
    return 0


def cmd_audit(portal: Portal, args) -> int:
    codes = 0
    for accession in args.accession:
        obj = portal.get(accession)
        if not args.json:
            print(_c(f"=== {accession}  ({type_of(obj)})", "1"))
            print(f"    {portal.ui_base}/{acc_of(obj.get('@id', accession))}/")
            print()
        findings = audit_object(portal, obj, deep=not args.shallow)
        if args.strict:
            findings = [f for f in findings if f.level in ("fail", "warn")]
        codes |= report(findings, args.json)
        # Recurse into member files -- a file set is only as releasable as its files.
        if args.recurse:
            for ref in obj.get("files") or []:
                child = portal.get(acc_of(ref))
                if not args.json:
                    print()
                    print(_c(f"--- member file {acc_of(ref)}  ({type_of(child)})", "1"))
                sub = audit_object(portal, child, deep=not args.shallow)
                if args.strict:
                    sub = [f for f in sub if f.level in ("fail", "warn")]
                codes |= report(sub, args.json)
    return codes


def cmd_inputs(portal: Portal, args) -> int:
    obj = portal.get(args.accession)
    print(_c(f"inputs of {args.accession} ({type_of(obj)})", "1"))
    print()
    rows = []
    for prop in ("derived_from", "input_file_sets", "reference_files", "files"):
        for ref in obj.get(prop) or []:
            ident = acc_of(ref)
            try:
                dep = portal.get(ident)
                st = dep.get("status", "?")
                sup = ", ".join(acc_of(s) for s in dep.get("superseded_by") or []) or "-"
                rows.append((prop, ident, st, sup, dep.get("revoke_detail") or ""))
            except PortalError:
                rows.append((prop, ident, "UNREADABLE", "-", ""))
    if not rows:
        print(info("no linked inputs."))
        return 0
    w = max(len(r[0]) for r in rows)
    for prop, ident, st, sup, detail in rows:
        flag = fail("") if st in DEAD_STATUSES or st == "UNREADABLE" else ok("")
        print(f"{flag} {prop:<{w}}  {ident:<16} {st:<12} superseded_by={sup}")
        if detail:
            print(f"       {detail}")
    dead = [r for r in rows if r[2] in DEAD_STATUSES]
    print()
    if dead:
        print(fail(f"{len(dead)} input(s) are unusable and will block release."))
        for prop, ident, st, sup, _ in dead:
            if sup != "-":
                print(f"       PATCH {prop}: {ident} -> {sup}")
        return 1
    print(ok("every input is in a usable status."))
    return 0


def cmd_fileset_files(portal: Portal, args) -> int:
    obj = portal.get(args.accession)
    refs = obj.get("files") or []
    out = []
    for ref in refs:
        ident = acc_of(ref)
        if args.field == "path":
            out.append(ref if isinstance(ref, str) else ref.get("@id"))
        else:
            out.append(ident)
    if args.filter_type:
        keep = []
        for ident in out:
            child = portal.get(acc_of(ident))
            if type_of(child) == args.filter_type:
                keep.append(ident)
        out = keep
    if args.plain:
        for x in out:
            print(x)
    else:
        # Paste-ready for the portal's JSON editing view.
        print(json.dumps(out, indent=2))
    print(f"\n{len(out)} file(s) in {args.accession}", file=sys.stderr)
    return 0


def cmd_validate(portal: Portal, args) -> int:
    findings = validate_local(args.path, expect_bed3=args.bed3)
    return report(findings, args.json)


def cmd_crosscheck(portal: Portal, args) -> int:
    """Are the variants revised in a corrected input still present in my file?

    The exact question the DACC asked on Aug 21 2026: an input file was revoked
    for reference-allele mismatches. Does the derived file need reuploading, or
    only a derived_from repoint?
    """
    tmp = args.workdir or "."
    os.makedirs(tmp, exist_ok=True)

    def resolve(spec):
        """Accession -> downloaded path; a local path passes straight through."""
        if os.path.exists(spec):
            return spec
        obj = portal.get(spec)
        ext = ".csv.gz" if obj.get("file_format") == "csv" else ".tsv.gz"
        dest = os.path.join(tmp, f"{acc_of(obj.get('@id', spec))}{ext}")
        if not os.path.exists(dest):
            portal.download(obj, dest)
        return dest

    old_path, new_path, mine_path = resolve(args.old), resolve(args.new), resolve(args.mine)
    old_map, old_meta = parse_variants(old_path)
    new_map, new_meta = parse_variants(new_path)
    mine_map, mine_meta = parse_variants(mine_path, args.chrom_col, args.pos_col)

    print(_c("column detection", "1"))
    for label, meta, path in (("old ", old_meta, old_path), ("new ", new_meta, new_path),
                              ("mine", mine_meta, mine_path)):
        print(f"  {label} {os.path.basename(path):<26} chrom={meta['chrom_col']} "
              f"pos={meta['pos_col']} id={meta['id_col']} "
              f"ref/alt={meta.get('ref_col')}/{meta.get('alt_col')}")
        print(f"       {meta['n']} distinct positions, {meta['with_alleles']} with alleles")
    print()

    if not old_meta["with_alleles"] or not new_meta["with_alleles"]:
        print(warn("could not read alleles from one of the input revisions; only whole-variant "
                   "removals can be detected, not in-place allele corrections."))

    dropped, changed = revised_positions(old_map, new_map, args.ignore_swaps)
    revised = dropped | changed
    print(_c("input revision", "1"))
    print(f"  positions dropped entirely:      {len(dropped)}")
    print(f"  positions with corrected alleles:{len(changed)}")
    print(f"  total revised positions:         {len(revised)}")
    print()

    # Control: if retained variants do not overlap my file at all, the
    # coordinate conventions differ and any negative result is meaningless.
    retained = (set(old_map) & set(new_map)) - changed
    control = sum(1 for k in retained if k in mine_map)
    print(_c("control", "1"))
    print(f"  unrevised input positions also in my file: {control} / {len(retained)}")
    if len(retained) and control == 0:
        print(fail("control FAILED -- zero overlap on unrevised variants, so the coordinate "
                   "conventions do not line up. The verdict below would be meaningless; "
                   "pass explicit --chrom-col/--pos-col."))
        return 2
    print(ok("control passed -- coordinates are comparable."))
    print()

    tol = args.tolerance
    impacted = []
    for chrom, pos in sorted(revised):
        for off in range(-tol, tol + 1):
            if (chrom, pos + off) in mine_map:
                kind = "dropped" if (chrom, pos) in dropped else "alleles corrected"
                impacted.append((chrom, pos, off, kind))
                break

    print(_c("verdict", "1"))
    if impacted:
        print(fail(f"{len(impacted)} revised input variant(s) ARE present in your file:"))
        for chrom, pos, off, kind in impacted:
            note = "" if off == 0 else f" (matched at offset {off:+d})"
            print(f"       chr{chrom}:{pos}{note}  [{kind}]")
            for pair in sorted(old_map.get((chrom, pos), set())):
                print(f"         was  ref/alt = {pair[0]}/{pair[1]}")
            for pair in sorted(new_map.get((chrom, pos), set())):
                print(f"         now  ref/alt = {pair[0]}/{pair[1]}")
        print()
        print("  -> Fix the alleles, re-upload the file, AND repoint derived_from.")
        return 1
    print(ok(f"none of the {len(revised)} revised input variants appear in your file."))
    print()
    print("  -> No reupload needed. Repoint derived_from to the corrected input, "
          "then the file can be released.")
    return 0


def cmd_patch(portal: Portal, args) -> int:
    payload = {}
    current = None

    # A REST PATCH REPLACES an array wholesale. Repointing one element of a
    # 25-entry derived_from by hand would silently drop the other 24, so build
    # the replacement array from the live value instead.
    for spec in args.repoint or []:
        if "=" not in spec or ":" not in spec.split("=", 1)[1]:
            die(f"--repoint expects PROP=OLD:NEW, got {spec!r}")
        prop, pair = spec.split("=", 1)
        old_acc, new_acc = pair.split(":", 1)
        if current is None:
            current = portal.get(args.accession)
        existing = current.get(prop)
        if not isinstance(existing, list):
            die(f"{prop} on {args.accession} is not an array (got {type(existing).__name__})")
        accs = [acc_of(x) for x in existing]
        if old_acc not in accs:
            die(f"{old_acc} is not in {prop} (current: {', '.join(accs)})")
        replaced = [new_acc if a == old_acc else a for a in accs]
        # de-duplicate while preserving order, in case NEW was already present
        seen, deduped = set(), []
        for a in replaced:
            if a not in seen:
                seen.add(a)
                deduped.append(a)
        payload[prop] = deduped
        print(info(f"{prop}: {len(accs)} -> {len(deduped)} entries, "
                   f"{old_acc} replaced by {new_acc}"))

    if args.json_payload:
        payload.update(json.loads(args.json_payload))
    for pair in args.set or []:
        if "=" not in pair:
            die(f"--set expects key=value, got {pair!r}")
        key, val = pair.split("=", 1)
        try:                       # let JSON through so arrays/bools work
            payload[key] = json.loads(val)
        except json.JSONDecodeError:
            payload[key] = val
    if not payload:
        die("nothing to patch: pass --set key=value or --json-payload")

    print(_c(f"PATCH {args.accession} on {portal.mode}", "1"))
    print(json.dumps(payload, indent=2))
    print()
    if not args.apply:
        print(warn("dry run -- nothing sent. Re-run with --apply to submit."))
        print(info("Equivalent igvf_utils call:"))
        rec = dict(payload, record_id=args.accession)
        print(f"       echo '{json.dumps(rec)}' > patch.json")
        print(f"       iu_register.py -m {portal.mode} -p {args.profile or '<profile>'} "
              f"-i patch.json --patch")
        return 0
    if portal.mode == "prod" and not args.yes:
        die("refusing to write to PRODUCTION without --yes "
            "(rehearse with -m staging first; sandbox is deprecated)")
    result = portal.patch(args.accession, payload)
    print(ok(f"patched {result.get('accession', args.accession)} -- "
             f"status={result.get('status')}"))
    return 0


def schema_required(prof: dict) -> "tuple[list, list]":
    """Required properties of a profile, including those inside a disjunction.

    Returns (always_required, alternative_groups).

    The IGVF schemas do not put every requirement in a top-level `required`.
    prediction_set has NO top-level required at all; it carries

        "oneOf": [{"required": ["lab","award","file_set_type","samples"]},
                  {"required": ["lab","award","file_set_type","donors"]}]

    -- so lab, award and file_set_type are always required, and you must
    supply samples OR donors. Reading only `required` reported "0 required"
    and emitted an empty template for the profile submitters use most. It
    also silently dropped `samples`, which is the exact property the Aug 2026
    meeting listed as a to-do ("Each prediction set needs a sample").

    allOf contributes unconditionally; anyOf and oneOf contribute their
    intersection unconditionally and the remainder as alternatives.
    """
    always = list(prof.get("required") or [])
    groups = []
    for key in ("allOf", "anyOf", "oneOf"):
        branches = prof.get(key) or []
        reqs = [list(b.get("required") or []) for b in branches
                 if isinstance(b, dict) and b.get("required")]
        if not reqs:
            continue
        if key == "allOf":
            for r in reqs:
                always.extend(x for x in r if x not in always)
            continue
        shared = set(reqs[0]).intersection(*[set(r) for r in reqs[1:]]) \
            if len(reqs) > 1 else set(reqs[0])
        always.extend(x for x in reqs[0] if x in shared and x not in always)
        alts = [[x for x in r if x not in shared] for r in reqs]
        alts = [a for a in alts if a]
        if len(alts) > 1:
            groups.append(alts)
        elif alts:
            always.extend(x for x in alts[0] if x not in always)
    return always, groups


def cmd_template(portal: Portal, args) -> int:
    prof = portal.profile(args.profile)
    required, alt_groups = schema_required(prof)
    props = prof.get("properties") or {}
    cols = list(required)
    # One column per alternative group, using the first branch, so the TSV is
    # submittable as-is; the note below says what the choice was.
    for alts in alt_groups:
        for name in alts[0]:
            if name not in cols:
                cols.append(name)
    for extra in args.include or []:
        if extra not in cols:
            cols.append(extra)
    print("\t".join(cols))
    print("\t".join(f"<{c}>" for c in cols))
    print(f"\n# profile: {args.profile}  ({len(props)} properties, "
          f"{len(required)} always required"
          f"{', ' + str(len(alt_groups)) + ' either/or group(s)' if alt_groups else ''})",
          file=sys.stderr)
    for alts in alt_groups:
        joined = "  OR  ".join(", ".join(a) for a in alts)
        print(f"#   either/or: {joined}", file=sys.stderr)
        print(f"#     -> this template uses {', '.join(alts[0])}; "
              f"delete it and use another branch instead if that fits your "
              f"submission.", file=sys.stderr)
    for c in cols:
        p = props.get(c, {})
        enum = p.get("enum") or (p.get("items") or {}).get("enum")
        line = f"#   {c}: {p.get('type', '?')}"
        if c in required:
            line += " [REQUIRED]"
        if enum:
            line += f"  enum({len(enum)})"
        print(line, file=sys.stderr)
        if enum and len(enum) <= 12:
            print(f"#       {enum}", file=sys.stderr)
    # iu_register.py's own mode flag: rehearse on staging. sandbox answers
    # HTTP 410 now, and this line used to recommend it.
    print("#\n# Submit with:\n#   iu_register.py -m staging -p "
          f"{args.profile} -i <this-file>.tsv -d   # dry run, then drop -d"
          "\n#   iu_register.py -m prod    -p "
          f"{args.profile} -i <this-file>.tsv -d   # only after staging passes",
          file=sys.stderr)
    return 0


def cmd_enum(portal: Portal, args) -> int:
    prof = portal.profile(args.profile)
    p = (prof.get("properties") or {}).get(args.property)
    if not p:
        die(f"{args.profile} has no property {args.property!r}")
    enum = p.get("enum") or (p.get("items") or {}).get("enum")
    if not enum:
        print(json.dumps(p, indent=2))
        return 0
    needle = (args.grep or "").lower()
    for v in enum:
        if needle in v.lower():
            print(v)
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

EPILOG = """\
examples:
  igvf-sub doctor
  igvf-sub audit IGVFDS2581EDPS --recurse
  igvf-sub inputs IGVFFI0856UJHP
  igvf-sub fileset-files IGVFDS5361GWAK           # paste-ready JSON array
  igvf-sub validate scores.bed.gz --bed3
  igvf-sub crosscheck --old IGVFFI7160EKDK --new IGVFFI1678CDBR --mine IGVFFI0856UJHP
  igvf-sub patch IGVFFI0856UJHP --set derived_from='["IGVFFI1678CDBR"]'
  igvf-sub template tabular_file --include derived_from reference_files
  igvf-sub enum tabular_file content_type --grep variant

Rehearse against staging before touching prod:  igvf-sub -m staging ...
(api.sandbox.igvf.org is deprecated and now returns HTTP 410.)
"""


def cmd_plan(portal: Portal, args) -> int:
    """Dependency-ordered submission plan from a draft JSON/YAML object list."""
    raw = Path(args.objects).read_text()
    try:
        objects = json.loads(raw)
    except json.JSONDecodeError:
        try:
            import yaml
            objects = yaml.safe_load(raw)
        except ImportError:
            die(f"{args.objects} is not JSON and PyYAML is not installed")
        except Exception as e:
            die(f"could not parse {args.objects}: {e}")
    if isinstance(objects, dict):
        objects = objects.get("objects") or [objects]
    out = plan_submission(objects)
    print(f"Submission plan — {len(out['steps'])} object(s), "
          f"{out['n_ready']} ready, {out['n_blocked']} blocked\n")
    for st in out["steps"]:
        mark = ok("ready") if st["ready"] else fail("BLOCKED")
        print(f"  {st['position']}. {st['type']:24} {mark}")
        if st.get("alias"):
            print(f"       alias: {st['alias']}")
        print(f"       {st['why_here']}")
        if st["missing_required"]:
            print(f"       missing: {', '.join(st['missing_required'])}")
    if out["unrecognised_types"]:
        print(warn(f"\n  not in the known submission order: "
                    f"{out['unrecognised_types']}"))
    print(f"\n{out['note']}")
    return 0 if out["n_blocked"] == 0 else 1


def cmd_check_revoked(portal: Portal, args) -> int:
    """Walk derived_from, flag dead inputs, and say which branch applies.

    The worked case this implements: a tabular file whose input was revoked
    for "ref allele mismatch to GRCh38, corrected in IGVFFI1678CDBR", with the
    offending alleles published as a .jsonl. Whether that needs a reupload or
    only a repoint depends on whether any of those alleles actually reached
    the derived file -- which is a question about data, not about metadata,
    and is why `crosscheck` exists.
    """
    obj = portal.get(args.accession)
    derived = obj.get("derived_from") or []
    print(f"=== {args.accession}  ({type_of(obj)})")
    print(f"  derived_from: {len(derived)} input(s)")
    dead = []
    for ref in derived:
        acc = acc_of(ref)
        try:
            up = portal.get(acc)
        except PortalError as e:
            print(warn(f"  {acc}: unreadable ({e})"))
            continue
        st = (up.get("status") or "").lower()
        if st in DEAD_STATUSES:
            note = (up.get("revoke_detail") or up.get("revoked_detail")
                     or up.get("description") or "")
            repl = None
            import re as _re
            m = _re.search(r"\b(IGVF[A-Z]{2}[0-9A-Z]{6,})\b", note or "")
            if m and m.group(1) != acc:
                repl = m.group(1)
            dead.append({"accession": acc, "status": st, "note": note,
                          "replacement": repl})
            print(fail(f"  {acc}: status={st}"))
            if note:
                print(f"       portal says: {note[:160]}")
            if repl:
                print(f"       replacement named in the notice: {repl}")
        else:
            print(ok(f"  {acc}: status={st}"))
    if not dead:
        print(ok("\n  no revoked or archived inputs — nothing to repair."))
        return 0
    print(f"\n  {len(dead)} dead input(s). Next step for each:")
    for d in dead:
        if d["replacement"]:
            print(f"    {d['accession']} -> {d['replacement']}")
            print(f"      Decide reupload vs repoint by checking whether the "
                  f"corrected records reached your file:")
            print(f"        igvfagent submit crosscheck --old {d['accession']} "
                  f"--new {d['replacement']} --mine {args.accession}")
            print(f"      If none did, repoint only:")
            print(f"        igvfagent submit patch {args.accession} "
                  f"--repoint derived_from={d['accession']}:{d['replacement']}")
        else:
            print(f"    {d['accession']}: the notice names no replacement; "
                  f"ask the DACC which accession supersedes it.")
    return 2


def cmd_status(portal: Portal, args) -> int:
    """A lab's file sets with release state and open audit counts."""
    lab = args.lab or os.environ.get("IGVF_LAB") or ""
    if not lab:
        die("no lab given and IGVF_LAB is unset")
    # Portal.search(item_type, limit=..., **filters) returns @graph already;
    # it does not take a params dict.
    rows = portal.search(args.item_type, limit=args.limit,
                          **{"lab.title": lab})
    print(f"=== {lab} — {len(rows)} {args.item_type}(s)")
    by_status = {}
    for r in rows:
        st = r.get("status", "?")
        by_status[st] = by_status.get(st, 0) + 1
        audits = r.get("audit") or {}
        n_audit = sum(len(v) for v in audits.values()) if isinstance(audits, dict) else 0
        flag = fail(f"{n_audit} audit") if n_audit else ok("clean")
        print(f"  {r.get('accession', '?'):16} {st:14} {flag}")
    print("\n  by status: " + ", ".join(f"{k}={v}" for k, v in
                                          sorted(by_status.items())))
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="igvfagent submit",
        description="Preflight, audit and submission helper for the IGVF data portal.",
        epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-m", "--mode", default=None,
                    help="prod | staging | hostname (default: $IGVF_MODE or prod); "
                         "'sandbox' is deprecated and maps to staging")
    ap.add_argument("-V", "--version", action="version", version=f"igvf-sub {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("doctor", help="check credentials and portal connectivity")
    p.set_defaults(fn=cmd_doctor)

    p = sub.add_parser("plan", help="dependency-ordered submission plan from a "
                                     "draft object list; sends nothing")
    p.add_argument("--objects", required=True,
                    help="JSON (or YAML) list of draft objects, each with a `type`")
    p.set_defaults(fn=cmd_plan)

    p = sub.add_parser("check-revoked",
                        help="walk derived_from and flag revoked/archived inputs")
    p.add_argument("accession")
    p.set_defaults(fn=cmd_check_revoked)

    p = sub.add_parser("status", help="a lab's file sets with release state and "
                                       "open audit counts")
    p.add_argument("--lab", default="", help="lab title; falls back to IGVF_LAB")
    p.add_argument("--type", dest="item_type", default="FileSet")
    p.add_argument("--limit", default="50")
    p.set_defaults(fn=cmd_status)

    p = sub.add_parser("get", help="fetch an object as JSON")
    p.add_argument("accession")
    p.add_argument("--frame", default="object", help="object | embedded | page (default: object)")
    p.add_argument("--field", nargs="*", help="print only these properties")
    p.set_defaults(fn=cmd_get)

    p = sub.add_parser("audit", help="run the recurring-DACC-findings checklist")
    p.add_argument("accession", nargs="+")
    p.add_argument("--recurse", action="store_true", help="also audit member files")
    p.add_argument("--shallow", action="store_true", help="skip fetching linked inputs")
    p.add_argument("--strict", action="store_true", help="hide passing checks")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_audit)

    p = sub.add_parser("inputs", help="list linked inputs and flag revoked/archived ones")
    p.add_argument("accession")
    p.set_defaults(fn=cmd_inputs)

    p = sub.add_parser("fileset-files", help="collect member file accessions from a file set")
    p.add_argument("accession")
    p.add_argument("--field", choices=["accession", "path"], default="accession")
    p.add_argument("--filter-type", help="keep only this @type, e.g. TabularFile")
    p.add_argument("--plain", action="store_true", help="one per line instead of JSON")
    p.set_defaults(fn=cmd_fileset_files)

    p = sub.add_parser("validate", help="validate a local file before upload")
    p.add_argument("path")
    p.add_argument("--bed3", action="store_true", help="apply 0-based half-open BED checks")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_validate)

    p = sub.add_parser("crosscheck",
                       help="do corrections in a revised input affect my derived file?")
    p.add_argument("--old", required=True, help="accession or path of the superseded input")
    p.add_argument("--new", required=True, help="accession or path of the corrected input")
    p.add_argument("--mine", required=True, help="accession or path of my derived file")
    p.add_argument("--tolerance", type=int, default=1,
                   help="positional slop when matching (default 1)")
    p.add_argument("--workdir", default=None, help="where to cache downloads")
    p.add_argument("--chrom-col", default=None, help="force my file's chromosome column")
    p.add_argument("--pos-col", default=None, help="force my file's position column")
    p.add_argument("--ignore-swaps", action="store_true",
                   help="treat a plain ref/alt swap as unchanged")
    p.set_defaults(fn=cmd_crosscheck)

    p = sub.add_parser("patch", help="PATCH an object (dry run unless --apply)")
    p.add_argument("accession")
    p.add_argument("--set", action="append", metavar="KEY=VALUE",
                   help="value is parsed as JSON when possible")
    p.add_argument("--json-payload", help="raw JSON object to merge in")
    p.add_argument("--repoint", action="append", metavar="PROP=OLD:NEW",
                   help="swap one accession inside an array property, preserving the rest")
    p.add_argument("--profile", help="profile name, for the printed iu_register equivalent")
    p.add_argument("--apply", action="store_true", help="actually send the request")
    p.add_argument("--yes", action="store_true", help="required to write to prod")
    p.set_defaults(fn=cmd_patch)

    p = sub.add_parser("template", help="emit a TSV header for iu_register from a profile")
    p.add_argument("profile", help="e.g. tabular_file, prediction_set, document")
    p.add_argument("--include", nargs="*", help="extra non-required columns to add")
    p.set_defaults(fn=cmd_template)

    p = sub.add_parser("enum", help="list a property's allowed values")
    p.add_argument("profile")
    p.add_argument("property")
    p.add_argument("--grep", help="substring filter")
    p.set_defaults(fn=cmd_enum)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    portal = Portal(mode=args.mode)
    try:
        return args.fn(portal, args)
    except PortalError as exc:
        die(str(exc))
    except KeyboardInterrupt:
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
