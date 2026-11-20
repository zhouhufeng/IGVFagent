#!/usr/bin/env python3
"""Run the IGVF Portal's own audits locally: on Portal objects, or on drafts before submitting.

The Portal (IGVF-DACC/igvfd) audits every object after it is saved; a
submitter used to learn the result only once the object was on the Portal,
often at a monthly meeting. This runs the same checkers, ported unchanged
into Scripts/igvfd_audit_rules.py (generated from igvfd/src/igvfd/audit at
the commit in IGVFD_COMMIT), against:

* objects already on the Portal (prod or staging), and
* draft objects (YAML/JSON/TSV, as `submit create` reads them) that are not
  submitted yet, with links to existing Portal objects followed for real.

What the server supplies to a checker is supplied here: `system['request']
.embed(path)` is a GET of the same path on the Portal API (or the draft with
that alias), `value` is the object in the frame the checker asked for, and the
type hierarchy is igvfd's own. Two Item checkers read the server's schema
registry and are ported by hand below (schema validation, linked statuses).
Checkers that need data a client cannot have (the ontology database) are
reported as skipped, never as passing.

    igvfagent submit -m staging portal-audit TSTFI69604689 IGVFDS2755YZNI --compare
    igvfagent submit -m staging portal-audit --drafts drafts.yaml
"""
from __future__ import annotations

import copy
import json
import posixpath
import re
import sys
import urllib.parse

try:
    from . import igvfd_audit_rules as R                       # type: ignore
    from .igvf_submission_skill import (Portal, PortalError, die, fail, info, ok, warn,  # type: ignore
                                        _load_objects, _coerce, _aliases, md5_of)
except ImportError:
    import igvfd_audit_rules as R                              # type: ignore
    from igvf_submission_skill import (Portal, PortalError, die, fail, info, ok, warn,  # type: ignore
                                       _load_objects, _coerce, _aliases, md5_of)

def _chrom_sizes(file_path):
    """igvfd's load_chrom_sizes_file reads src/igvfd/audit/_static/<assembly>.chrom.sizes;
    the same tables are embedded in the generated rules (CHROM_SIZES)."""
    return dict(R.CHROM_SIZES[str(file_path).rsplit("/", 1)[-1].split(".")[0]])


R.load_chrom_sizes_file = _chrom_sizes

LEVELS = {"ERROR": 60, "NOT_COMPLIANT": 50, "WARNING": 40, "INTERNAL_ACTION": 30}
LEVEL_NAMES = {v: k for k, v in LEVELS.items()}


class NeedsServerData(Exception):
    """A checker reached for data only the Portal server holds (its ontology)."""


class _Registry(dict):
    def __getitem__(self, key):
        raise NeedsServerData(f"registry[{key!r}] (Portal-side {key} data)")


# --------------------------------------------------------------- context ----

class AuditContext:
    """`request.embed` for the ported checkers: Portal GETs and local drafts, cached."""

    def __init__(self, portal: "Portal | None", drafts: "list[dict] | None" = None):
        self.portal = portal
        self.cache: "dict[str, object]" = {}
        self.drafts: "dict[str, dict]" = {}
        for d in drafts or []:
            for key in [d.get("@id")] + _aliases(d):
                if key:
                    self.drafts[key] = d
        self.n_get = 0

    def _draft_for(self, base: str) -> "dict | None":
        base = base.strip()
        if base in self.drafts:
            return self.drafts[base]
        parts = [p for p in base.split("/") if p]
        if parts and parts[-1] in self.drafts:
            return self.drafts[parts[-1]]
        return None

    def embed(self, *elements) -> dict:
        path = posixpath.join(*[str(e) for e in elements])
        if path in self.cache:
            hit = self.cache[path]
            if isinstance(hit, dict) and "__portal_error__" in hit:
                raise PortalError(hit["__portal_error__"])
            # A copy: igvfd checkers modify what they are given
            # (audit_input_for drains value['input_for'] while walking it).
            return copy.deepcopy(hit)
        base, sep, view = path.partition("@@")
        draft = self._draft_for(base)
        if draft is not None:
            out = _draft_view(draft, view)
        else:
            if self.portal is None:
                raise PortalError(f"no Portal to resolve {path}")
            p = path if path.startswith("/") else "/" + path
            raw_key = "raw:" + p
            if raw_key not in self.cache:
                self.n_get += 1
                try:
                    self.cache[raw_key] = self.portal._request("GET", p)
                except PortalError as e:
                    # Remembered like an answer, so a run (and a recorded
                    # fixture) sees the same refusal every time.
                    self.cache[path] = {"__portal_error__": str(e)}
                    raise
            out = self.cache[raw_key]
        if "skip_calculated=true" in view:
            out = _skip_calculated(out)
        self.cache[path] = out
        return copy.deepcopy(out)

    def get_frame(self, ident: str, frame: str) -> dict:
        """The object as an audit_checker frame gives it ('object', 'embedded', ...)."""
        draft = self._draft_for(ident)
        if draft is not None:
            return draft
        base = ident if ident.startswith("/") else f"/{ident}/"
        if not base.endswith("/"):
            base += "/"
        if frame.startswith("object"):
            q = frame[len("object"):]
            return self.embed(base + "@@object" + q)
        return self.embed(base)          # the default (page) frame embeds linked objects


class _Request:
    def __init__(self, ctx: AuditContext):
        self._ctx = ctx

    def embed(self, *elements, **kw):
        return self._ctx.embed(*elements)

    @property
    def root(self):
        # request.root.get_by_uuid(...) reads the server's own storage (which
        # S3 bucket a file sits in); a client cannot see it.
        raise NeedsServerData("request.root (Portal-side file storage)")


def _draft_view(d: dict, view: str) -> dict:
    name, _, query = view.partition("?")
    if name == "filtered_object":
        keep = urllib.parse.parse_qs(query).get("include", [])
        return {k: d[k] for k in keep if k in d}
    return d


# ------------------------------------------------- hand-ported Item checks ---

def _schema_links(schema: dict, prefix: str = "", skip=()) -> "list[str]":
    """snovault typeinfo.extract_schema_links: dotted paths of every linkTo property.

    The server reads its type's own schema file, which has no calculated
    properties; the served profile adds them, so they are skipped here.
    """
    out = []
    for name, prop in (schema.get("properties") or {}).items():
        if not prefix and name in skip:
            continue
        path = f"{prefix}{name}"
        target = prop["items"] if "items" in prop else prop
        if "properties" in target:
            out += _schema_links(target, path + ".")
        elif "linkTo" in target:
            out.append(path)
    return sorted(out)


def _skip_calculated(obj: dict) -> dict:
    """snovault's ?skip_calculated=true frame: stored properties only.

    @id, @type and uuid are calculated properties in snovault, so they go too;
    a checker on this frame takes the object's path from system['path'].
    The public API ignores the parameter, hence this.
    """
    if not isinstance(obj, dict) or "@type" not in obj:
        return obj
    calc = set(R.CALCULATED_PROPERTIES.get(obj["@type"][0], [])) | {"@id", "@type", "uuid"}
    return {k: v for k, v in obj.items() if k not in calc}


def audit_item_schema(value, system):
    """Ported by hand from igvfd/audit/item.py:audit_item_schema.

    The server validates the stored properties against the type's schema after
    upgrading them to the current schema version. A client sees the current
    version already (the Portal upgrades on read), so the upgrade step is not
    needed; the profile comes from /profiles/<type>.json.
    """
    schema = system.get("schema")
    if not schema:
        return
    try:
        import jsonschema
    except ImportError:
        return
    calculated = set(R.CALCULATED_PROPERTIES.get(value["@type"][0], []))
    props = {k: v for k, v in system.get("stored", value).items()
             if not k.startswith("@") and k not in calculated}
    props.setdefault("uuid", value.get("uuid", "00000000-0000-0000-0000-000000000000"))
    for k in ("uuid",):
        if not props.get(k):
            props.pop(k, None)
    sch = {k: v for k, v in schema.items() if k not in ("mixinProperties",)}
    try:
        validator = jsonschema.Draft202012Validator(sch)
    except Exception:
        validator = jsonschema.Draft7Validator(sch)
    for error in sorted(validator.iter_errors(props), key=lambda e: list(e.path)):
        category = "validation error"
        path = list(error.path)
        if path:
            category += ": " + "/".join(str(elem) for elem in path)
        detail = ("{} {} has schema error {}.".format(
            R.space_in_words(value["@type"][0]).capitalize(),
            R.audit_link(R.path_to_text(value["@id"]), value["@id"]),
            error.message))
        yield R.AuditFailure(category, detail, level="INTERNAL_ACTION")


def audit_item_status(value, system):
    """Ported by hand from igvfd/audit/item.py:audit_item_status.

    The server takes the linked paths from its type registry
    (context.type_info.schema_links); here they come from the profile's
    linkTo properties, which is what that registry is built from.
    """
    if "status" not in value:
        return
    if value["status"] in ["revoked", "disabled", "deleted", "replaced"]:
        return
    level = R.STATUS_HIERARCHY.get(value["status"], 50)
    if value["status"] in ["archived"]:
        level += 30
    request = system["request"]
    linked = set()
    calc = R.CALCULATED_PROPERTIES.get(value["@type"][0], [])
    for schema_path in _schema_links(system.get("schema") or {}, skip=calc):
        if schema_path == "supersedes":
            continue
        linked.update(x for x in R.simple_path_ids(value, schema_path) if isinstance(x, str))
    for path in sorted(linked):
        try:
            linked_value = request.embed(path + "@@filtered_object?include=@id&include=@type&include=uuid&include=status")
        except PortalError as e:
            # The server sees every object; these keys may not (an unreleased
            # object of another lab answers 403). Say so rather than pass it.
            if "HTTP 403" in str(e) or "HTTP 401" in str(e):
                yield R.AuditFailure("linked object not visible",
                                     f"{path} is linked but cannot be read with these keys, so its status "
                                     f"could not be checked.", level="UNVERIFIABLE")
            continue
        if "status" not in linked_value:
            continue
        if linked_value["status"] == "disabled":
            continue
        linked_level = R.STATUS_HIERARCHY.get(linked_value["status"], 50)
        if linked_value["status"] in ["archived"]:
            linked_level += 30
        if linked_value["status"] in ["deleted"]:
            linked_level += 100
        if linked_level == level:
            continue
        elif linked_level < level:
            detail = ("{} {} {} has {} subobject {} {}.".format(
                value["status"].capitalize(),
                R.space_in_words(value["@type"][0]).lower(),
                R.audit_link(R.path_to_text(value["@id"]), value["@id"]),
                linked_value["status"],
                R.space_in_words(linked_value["@type"][0]).lower(),
                R.audit_link(R.path_to_text(linked_value["@id"]), linked_value["@id"])))
            yield R.AuditFailure("mismatched status", detail, level="INTERNAL_ACTION")
        else:
            if "File" in value["@type"] and "FileSet" in linked_value["@type"]:
                detail = ("{} {} {} has {} subobject {} {} which is ahead in status.".format(
                    value["status"].capitalize(),
                    R.space_in_words(value["@type"][0]).lower(),
                    R.audit_link(R.path_to_text(value["@id"]), value["@id"]),
                    linked_value["status"],
                    R.space_in_words(linked_value["@type"][0]).lower(),
                    R.audit_link(R.path_to_text(linked_value["@id"]), linked_value["@id"])))
                yield R.AuditFailure("mismatched status", detail, level="INTERNAL_ACTION")


CHECKERS = list(R.CHECKERS) + [
    ("Item", "object", audit_item_schema, "item.py", None),
    ("Item", "object", audit_item_status, "item.py", None),
]


def _dispatch_table(fn):
    """The sub-checks of an igvfd *_object_dispatcher, so findings name the check that fired."""
    name = fn.__name__
    m = re.match(r"audit_(.+)_dispatcher(__\w+)?$", name)
    if not m:
        return None
    table = getattr(R, f"function_dispatcher_{m.group(1)}{m.group(2) or ''}", None)
    return table if isinstance(table, dict) else None


# ------------------------------------------------------------------ run ------

def type_lineage(obj: dict) -> "list[str]":
    t = obj.get("@type") or []
    if t:
        return list(t)
    name = R.ITEM_TYPES.get(str(obj.get("type", "")), str(obj.get("type", "")))
    return R.TYPE_LINEAGE.get(name, [name, "Item"])


def run_audits(obj: dict, ctx: AuditContext, schema: "dict | None" = None) -> "list[dict]":
    """Every igvfd checker that applies to `obj`, run as the Portal runs it.

    Returns dicts with category, level (name), detail, checker (the sub-check),
    dispatcher (the registered checker the Portal reports as `name`), path.
    A checker that needs server-only data or raises is returned as
    level 'SKIPPED' / 'CHECKER_ERROR' with the reason, never silently dropped.
    """
    lineage = type_lineage(obj)
    ident = obj.get("@id") or (_aliases(obj) or [""])[0]
    request = _Request(ctx)
    if schema is None and ctx.portal is not None:
        try:
            schema = ctx.portal._request("GET", f"/profiles/{lineage[0]}.json")
        except PortalError:
            schema = None
    frames: "dict[str, dict]" = {"object": obj}
    out = []
    for item_type, frame, fn, src, _cond in CHECKERS:
        if item_type not in lineage:
            continue
        if frame not in frames:
            try:
                frames[frame] = ctx.get_frame(ident, frame)
            except PortalError:
                frames[frame] = obj      # a draft or an in-memory object: one frame only
            if "skip_calculated=true" in frame:
                frames[frame] = _skip_calculated(frames[frame]) if "@type" in frames[frame] else frames[frame]
        value = frames[frame]
        system = {"request": request, "context": None, "registry": _Registry(),
                  "path": obj.get("@id"), "types": obj.get("@type") or lineage,
                  "schema": schema, "stored": frames.get("object?skip_calculated=true") or value}
        table = _dispatch_table(fn)
        subs = list(table.items()) if table else [(fn.__name__, fn)]
        for sub_name, sub in subs:
            try:
                for f in sub(copy.deepcopy(value), system):
                    lvl = f.level if isinstance(f.level, str) else LEVEL_NAMES.get(f.level, str(f.level))
                    out.append({"category": f.category, "level": lvl, "detail": f.detail,
                                "checker": sub_name, "dispatcher": fn.__name__.split("__")[0],
                                "file": src, "path": obj.get("@id")})
            except NeedsServerData as e:
                out.append({"category": "skipped", "level": "SKIPPED", "detail": str(e),
                            "checker": sub_name, "dispatcher": fn.__name__.split("__")[0],
                            "file": src, "path": value.get("@id")})
            except Exception as e:                              # noqa: BLE001
                out.append({"category": "checker error", "level": "CHECKER_ERROR",
                            "detail": f"{type(e).__name__}: {str(e)[:200]}", "checker": sub_name,
                            "dispatcher": fn.__name__.split("__")[0], "file": src,
                            "path": value.get("@id")})
    return out


def portal_audits(portal: Portal, ident: str) -> "list[dict]":
    """The Portal's own computed audits of one object (@@audit-self: not its linked objects')."""
    base = ident if ident.startswith("/") else f"/{ident}/"
    data = portal._request("GET", base.rstrip("/") + "/@@audit-self")
    return [{"category": a.get("category"), "level": a.get("level_name"), "detail": a.get("detail"),
             "dispatcher": a.get("name"), "path": a.get("path")} for a in data.get("audit") or []]


_LINK = re.compile(r"\[[^\]]*\]\([^)]*\)")


def _norm(detail: str) -> str:
    """Detail text for matching, with its list of linked objects in a fixed order.

    igvfd builds several lists from a Python set, so the Portal's own order of
    links in one message varies from run to run; the text around them does not.
    """
    text = re.sub(r"\s+", " ", str(detail or "")).strip()
    # Lists of property names (`a`, `b`, `c`) also come out of a set.
    text = re.sub(r"(`[^`]+`(?:, `[^`]+`)+)",
                  lambda m: ", ".join(sorted(m.group(1).split(", "))), text)
    links = sorted(_LINK.findall(text))
    return _LINK.sub("@", text) + " | " + " ".join(links)


def compare(ours: "list[dict]", theirs: "list[dict]") -> dict:
    """Match findings by (category, level, detail); a looser (category, level) match is reported too."""
    real = [o for o in ours if o["level"] in LEVELS]
    key = lambda a: (a["category"], a["level"], _norm(a["detail"]))
    loose = lambda a: (a["category"], a["level"])
    pool = [key(t) for t in theirs]
    matched, extra = [], []
    for o in real:
        if key(o) in pool:
            pool.remove(key(o))
            matched.append(o)
        else:
            extra.append(o)
    missed = [t for t in theirs if key(t) in pool and (pool.remove(key(t)) or True)]
    loose_pool = [loose(m) for m in missed]
    loose_hits = 0
    for e in extra:
        if loose(e) in loose_pool:
            loose_pool.remove(loose(e))
            loose_hits += 1
    return {"matched": matched, "extra": extra, "missed": missed, "category_level_matches": loose_hits}


# ------------------------------------------------------------------ CLI ------

def _print_findings(rows, title):
    print(f"\n=== {title}")
    if not rows:
        print(ok("no audits"))
    for r in sorted(rows, key=lambda r: -LEVELS.get(r["level"], 0)):
        mark = {"ERROR": fail, "NOT_COMPLIANT": fail, "WARNING": warn}.get(r["level"], info)
        print(mark(f"[{r['level']}] {r['category']} ({r['checker']}): {_norm(r['detail'])[:300]}"))


def _as_posted(d: dict, name: str, alias: str) -> dict:
    """A draft as the Portal would hold it right after `submit create` posts it.

    The same fills as create (lab/award from IGVF_LAB/IGVF_AWARD, md5sum and
    file_size of the local file), plus what the server sets on a new object
    (status in progress; upload_status pending on a file until the upload
    validates).
    """
    import os
    from pathlib import Path
    v = {k: x for k, x in d.items() if k not in ("type", "_local_file", "_attachment_file")}
    for k, env in (("lab", "IGVF_LAB"), ("award", "IGVF_AWARD")):
        if os.environ.get(env) and not v.get(k):
            v[k] = os.environ[env]
    local = d.get("_local_file") or d.get("submitted_file_name")
    lineage = R.TYPE_LINEAGE.get(name, [name, "Item"])
    if "File" in lineage:
        if local and Path(local).is_file():
            v.setdefault("md5sum", md5_of(local))
            v.setdefault("file_size", Path(local).stat().st_size)
        v.setdefault("upload_status", "pending")
    if d.get("_attachment_file") and not v.get("attachment"):
        v["attachment"] = {"download": Path(d["_attachment_file"]).name}
    v["@type"] = lineage
    v["@id"] = f"/{R.COLLECTIONS.get(name, name.lower())}/{alias}/"
    v.setdefault("status", "in progress")
    return v


def _normalise_links(portal, prepared) -> None:
    """Links as the Portal stores them: @id paths, not bare accessions or aliases.

    Drafts link by accession (IGVFFI0653VCGH) or by another draft's alias; on
    the Portal every link is stored as the target's @id, and the checkers
    build URLs and compare links on that form.
    """
    local = {}
    for v, _ in prepared:
        for k in _aliases(v):
            local[k] = v["@id"]
    resolved = {}

    def fix(x):
        if not isinstance(x, str) or x.startswith("/"):
            return x
        if x in local:
            return local[x]
        if x not in resolved:
            try:
                resolved[x] = portal._request("GET", f"/{x}/@@object").get("@id", x)
            except PortalError:
                resolved[x] = x
        return resolved[x]

    for v, schema in prepared:
        for path in _schema_links(schema or {}):
            if "." in path or path not in v:
                continue
            v[path] = [fix(x) for x in v[path]] if isinstance(v[path], list) else fix(v[path])


def _reverse_links(drafts: "list[dict]") -> None:
    """The calculated reverse links the audits read most (a file set's `files`), among the drafts."""
    by_key = {}
    for v in drafts:
        for k in [v["@id"]] + _aliases(v):
            by_key[k] = v
    for v in drafts:
        if "File" in v["@type"] and v.get("file_set") in by_key:
            fs = by_key[v["file_set"]]
            fs.setdefault("files", []).append(v["@id"])


def cmd_portal_audit(portal: Portal, args) -> int:
    rows_all = []
    if args.drafts:
        drafts = []
        for spec in args.drafts:
            drafts += _load_objects(spec, args.profile or "")
        prepared = []
        for d in drafts:
            name = R.ITEM_TYPES.get(str(d.get("type", "")), "")
            if not name:
                die(f"draft {(_aliases(d) or ['?'])[0]}: unknown type {d.get('type')!r}")
            try:
                schema = portal.profile(d["type"])
                _coerce(d, schema)
            except PortalError:
                schema = None
            alias = (_aliases(d) or ["draft"])[0]
            v = _as_posted(d, name, alias)
            prepared.append((v, schema))
        _normalise_links(portal, prepared)
        _reverse_links([v for v, _ in prepared])
        ctx = AuditContext(portal, [v for v, _ in prepared])
        for v, schema in prepared:
            rows = run_audits(v, ctx, schema=schema)
            rows_all += rows
            _print_findings(rows, f"{v['@type'][0]} {v['@id']} (draft)")
    for ident in args.accessions or []:
        ctx = AuditContext(portal)
        obj = ctx.get_frame(ident, "object")
        rows = run_audits(obj, ctx)
        rows_all += rows
        _print_findings(rows, f"{obj['@type'][0]} {obj['@id']}")
        if args.compare:
            cmp = compare(rows, portal_audits(portal, obj["@id"]))
            print(info(f"vs the Portal's own audits: {len(cmp['matched'])} match, "
                       f"{len(cmp['missed'])} only on the Portal, {len(cmp['extra'])} only here"))
            for m in cmp["missed"]:
                print(warn(f"  Portal only: [{m['level']}] {m['category']}: {_norm(m['detail'])[:200]}"))
            for e in cmp["extra"]:
                print(warn(f"  here only:   [{e['level']}] {e['category']}: {_norm(e['detail'])[:200]}"))
    if args.json:
        print(json.dumps(rows_all, indent=1, default=str))
    blocking = [r for r in rows_all if r["level"] in ("ERROR", "NOT_COMPLIANT")]
    print(f"\n{len(rows_all)} finding(s); {len(blocking)} ERROR/NOT_COMPLIANT "
          f"(igvfd {R.IGVFD_COMMIT})")
    return 1 if blocking else 0


def register(sub) -> None:
    p = sub.add_parser("portal-audit",
                       help="run the IGVF Portal's own audits (igvfd) on Portal objects or drafts")
    p.add_argument("accessions", nargs="*", help="objects on the Portal (accession, uuid or @id)")
    p.add_argument("--drafts", nargs="+", action="extend",
                   help="draft objects (JSON/YAML/TSV/CSV, as `create` reads them) to audit before submitting")
    p.add_argument("--profile", help="profile for sheet rows with no `type` column")
    p.add_argument("--compare", action="store_true",
                   help="also fetch each object's audits from the Portal and compare")
    p.add_argument("--json", action="store_true")
    p.set_defaults(fn=cmd_portal_audit)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(prog="igvfd-audits")
    ap.add_argument("-m", "--mode", default=None)
    sp = ap.add_subparsers(dest="cmd", required=True)
    register(sp)
    a = ap.parse_args()
    sys.exit(a.fn(Portal(mode=a.mode), a))
