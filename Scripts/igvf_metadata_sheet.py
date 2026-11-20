#!/usr/bin/env python3
"""IGVF-DACC igvf-metadata-submitter, as `igvfagent submit sheet-*` commands.

The DACC's metadata submitter (github.com/IGVF-DACC/igvf-metadata-submitter,
ported from commit 944f35b) is a Google Sheets add-on: one sheet per profile,
one row per object, a menu that GETs, POSTs, PATCHes and PUTs rows and uploads
files. It is what most labs submit with. This module ports its row <-> JSON
rules and its actions to local TSV / CSV / XLSX sheets, so a sheet made for the
add-on behaves the same here and a sheet made here opens in the add-on.

    igvfagent submit sheet-template tabular_file            # header + default row
    igvfagent submit sheet-get files.tsv                    # fill rows from the Portal
    igvfagent submit sheet-post files.tsv [--apply]         # POST rows, accessions written back
    igvfagent submit sheet-patch files.tsv --columns description [--apply]
    igvfagent submit sheet-patch files.tsv --all [--apply]
    igvfagent submit sheet-patch files.tsv --remove notes [--columns x] [--apply]
    igvfagent submit sheet-put files.tsv --admin [--apply]
    igvfagent submit sheet-upload files.tsv [--print-cmd | --apply]

Every write is a dry run unless --apply; prod also needs --yes. Results go to a
NEW sheet next to the input (never the input itself), printed as `Wrote: <path>`.
Like the add-on, these commands default to staging when no mode is given.

The sheet's profile comes from --profile, else from the file name (the add-on
keeps it in the sheet's developer metadata; a file `tabular_file.tsv` or
`tabular_file.anything.tsv` is the same profile).

Left out, because they only exist inside Google Sheets: header colours and
tooltips (written instead to a `<sheet>.legend.tsv` side file), dropdown
validation for enums (the legend lists the allowed values), the search dialog
for linked properties (the legend names the type each links to), Google Drive
attachment paths (an `attachment` cell here is `{"path": "<local file>"}`),
and the upload sidebar (`sheet-upload --apply` uploads directly).
"""
from __future__ import annotations

import csv
import datetime as _dt
import json
import mimetypes
import os
import re
import sys
from pathlib import Path

try:                                                     # inside the igvfagent package
    from .igvf_submission_skill import (Portal, PortalError, die, ok, warn, fail, info,
                                        _upload_bytes, md5_of, validate_local, _schema_errors)
except ImportError:                                      # standalone, from Scripts/
    from igvf_submission_skill import (Portal, PortalError, die, ok, warn, fail, info,  # type: ignore
                                       _upload_bytes, md5_of, validate_local, _schema_errors)

# ── constants (Metadata.js, Profile.js) ─────────────────────────────────────

SKIP = "#skip"
RESPONSE = "#response"
RESPONSE_TIME = "#response_time"
UPLOAD_ABSPATH = "#upload_abspath"
UPLOAD_STATUS = "#upload_status"
UPLOAD_CMD = "#upload_cmd"

# Column order in the header (Profile.js DEFAULT_PROP_PRIORITY).
DEFAULT_PROP_PRIORITY = [SKIP, RESPONSE, RESPONSE_TIME, UPLOAD_ABSPATH, UPLOAD_STATUS,
                         UPLOAD_CMD, "accession", "uuid", "name", "aliases", "award", "lab"]
# Which identifying property wins when a row has several (Profile.js).
IDENTIFYING_PROP_PRIORITY = ["accession", "uuid", "name"]
BIG_NUMBER_FOR_PRIORITY_SORTING = 1000
# "permission" values that make a property admin/system only (Profile.js).
ADMIN_OR_SYSTEM_PERMISSIONS = ["add_unvalidated", "edit_unvalidated", "expand", "impersonate",
                               "import_items", "index", "submit_for_any", "view_raw"]
SELECTED_PROP_KEYS_FOR_TOOLTIP = ["title", "description", "comment", "type", "readonly",
                                  "notSubmittable", "linkTo"]


# ── value helpers (Library.js) ──────────────────────────────────────────────

def to_boolean(val) -> bool:
    """Library.js toBoolean: "1", "true", "t", "o" (any case) are true."""
    return str(val).lower() in ("1", "true", "t", "o")


def _is_array_string(s: str) -> bool:
    """Library.js isArrayString."""
    t = s.strip()
    return t.startswith("[") and t.endswith("]")


def _is_json_string(s: str) -> bool:
    """Library.js isJsonString."""
    t = s.strip()
    return t.startswith("{") and t.endswith("}")


def camel_to_snake(name: str) -> str:
    """Library.js camelToSnake(uncapitalizeWord(name)): MeasurementSet -> measurement_set."""
    if not name:
        return name
    name = name[0].lower() + name[1:]
    return re.sub(r"[A-Z]", lambda m: "_" + m.group(0).lower(), name)


def now_string() -> str:
    """Library.js getCurrentLocalTimeString(""): local time, 'YYYY-MM-DD HH:MM:SS.mmm'."""
    return _dt.datetime.now().isoformat(sep=" ", timespec="milliseconds")


_NUM = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$")


def auto_parse(text):
    """What Google Sheets makes of a typed cell, before the add-on reads it.

    Sheets turns TRUE/FALSE into booleans and numeric text into numbers; dates
    the add-on reads back as their displayed text (rowToJson with
    bypassGoogleAutoParsing), which is the text itself here. A TSV/CSV cell is
    always text, so this reproduces that step. XLSX cells already carry types.
    """
    if not isinstance(text, str):
        if isinstance(text, (_dt.datetime, _dt.date)):
            return text.isoformat()
        return text
    if text.upper() in ("TRUE", "FALSE"):
        return text.upper() == "TRUE"
    if _NUM.match(text.strip()):
        v = float(text)
        return int(v) if v.is_integer() and not re.search(r"[.eE]", text) else v
    return text


# ── profile (Profile.js) ────────────────────────────────────────────────────

def _prop(profile, prop) -> dict:
    return (profile.get("properties") or {}).get(prop) or {}


def is_commented_prop(prop: str) -> bool:
    """Profile.js isCommentedProp."""
    return prop.startswith("#")


def is_searchable_prop(profile, prop) -> bool:
    """Profile.js isSearchableProp: linkTo on the property or on its items."""
    if not profile or prop.startswith("#") or prop not in (profile.get("properties") or {}):
        return False
    p = _prop(profile, prop)
    return "linkTo" in p or "linkTo" in (p.get("items") or {})


def link_to(profile, prop):
    p = _prop(profile, prop)
    return p.get("linkTo") or (p.get("items") or {}).get("linkTo")


def get_prop_type(profile, prop):
    """Profile.js getPropType."""
    if not profile or prop not in (profile.get("properties") or {}):
        return None
    return _prop(profile, prop).get("type")


def is_array_prop(profile, prop) -> bool:
    """Profile.js isArrayProp."""
    return bool(profile) and not prop.startswith("#") and get_prop_type(profile, prop) == "array"


def is_required_prop(profile, prop) -> bool:
    """Profile.js isRequiredProp: `required`, else recursively in `anyOf`.

    Exactly as the add-on: a profile whose requirements sit only under
    `oneOf`/`allOf` (prediction_set) reports none as required here. The
    stricter reading is igvf_submission_skill.schema_required.
    """
    if "required" in profile:
        return prop in (profile.get("required") or [])
    for sub in profile.get("anyOf") or []:
        if is_required_prop(sub, prop):
            return True
    return False


def is_identifying_prop(profile, prop) -> bool:
    """Profile.js isIdentifyingProp."""
    return prop in (profile.get("identifyingProperties") or [])


def is_readonly_prop(profile, prop) -> bool:
    """Profile.js isReadonlyProp."""
    return bool(_prop(profile, prop).get("readonly"))


def is_not_submittable_prop(profile, prop) -> bool:
    """Profile.js isNotSubmittableProp."""
    return bool(_prop(profile, prop).get("notSubmittable"))


def is_non_editable_prop(profile, prop) -> bool:
    """Profile.js isNonEditableProp."""
    return is_readonly_prop(profile, prop) or is_not_submittable_prop(profile, prop)


def has_do_not_submit_in_prop_comment(profile, prop) -> bool:
    """Profile.js hasDoNotSubmitInPropComment."""
    return str(_prop(profile, prop).get("comment") or "").lower().startswith("do not submit")


def is_admin_or_system_prop(profile, prop) -> bool:
    """Profile.js isAdminOrSystemProp (the add-on's permission list; IGVF's
    `admin_only` is not on it, so those properties are caught by readonly /
    "do not submit" instead, as in the add-on)."""
    return _prop(profile, prop).get("permission") in ADMIN_OR_SYSTEM_PERMISSIONS


def has_attachment(profile) -> bool:
    """Profile.js hasAttachment."""
    return bool(_prop(profile, "attachment").get("attachment"))


def is_gettable_prop(profile, prop, for_admin=False) -> bool:
    """Profile.js isGettableProp (also isPostableProp)."""
    if prop not in (profile.get("properties") or {}):
        return False
    if for_admin:
        return not is_not_submittable_prop(profile, prop)
    return (is_required_prop(profile, prop) or is_identifying_prop(profile, prop)
            or (not is_non_editable_prop(profile, prop)
                and not has_do_not_submit_in_prop_comment(profile, prop)
                and not is_admin_or_system_prop(profile, prop)))


is_postable_prop = is_gettable_prop


def get_default_for_prop(profile, prop):
    """Profile.js getDefaultForProp (None when there is no default)."""
    return _prop(profile, prop).get("default")


def get_identifying_prop_priority(prop) -> int:
    """Profile.js getIdentifyingPropPriority."""
    return IDENTIFYING_PROP_PRIORITY.index(prop) if prop in IDENTIFYING_PROP_PRIORITY \
        else BIG_NUMBER_FOR_PRIORITY_SORTING


def schema_version(profile):
    """Profile.js getProfileSchemaVersion."""
    return _prop(profile, "schema_version").get("default")


def type_cast_value(profile, prop, val, raw=None):
    """Profile.js typeCastValueByProfile: a number in a string-typed property
    becomes a string. Sheets would have turned "007" into 7 and back into "7";
    a text cell here keeps the text it was typed with (`raw`). A boolean in a
    string property is kept as typed text too (the add-on would send `true`
    and fail validation)."""
    if _prop(profile, prop).get("type") == "string" and isinstance(val, (int, float, bool)):
        if isinstance(raw, str):
            return raw
        return str(val) if not isinstance(val, bool) else ("TRUE" if val else "FALSE")
    return val


def type_cast_json(profile, obj, raw=None, keep_commented=False) -> dict:
    """Profile.js typeCastJsonValuesByProfile."""
    out = {}
    for prop, val in obj.items():
        if prop.startswith("#") and not keep_commented:
            continue
        out[prop] = type_cast_value(profile, prop, val, (raw or {}).get(prop))
    return out


def get_sorted_props(props, profile, priority=None) -> list:
    """Metadata.js getSortedProps: priority props, then required, then the
    rest in their own order, then commented (#) props."""
    priority = DEFAULT_PROP_PRIORITY if priority is None else priority
    out = []
    for p in list(priority) + list(profile.get("required") or []):
        if p in props and p not in out:
            out.append(p)
    for p in props:
        if not p.startswith("#") and p not in out:
            out.append(p)
    for p in props:
        if p.startswith("#") and p not in out:
            out.append(p)
    return out


def make_template(profile, for_admin=False) -> dict:
    """MetadataTemplate.js makeMetadataTemplateFromProfile."""
    out = {RESPONSE: None, RESPONSE_TIME: None}
    for prop in (profile.get("properties") or {}):
        if is_postable_prop(profile, prop, for_admin):
            out[prop] = get_default_for_prop(profile, prop)
    return out


def color_for_prop(profile, prop) -> str:
    """Profile.js getColorForProp, as the legend's role column."""
    if is_commented_prop(prop):
        return "commented"
    if is_required_prop(profile, prop):
        return "required (red)"
    if is_identifying_prop(profile, prop):
        return "identifying (blue)"
    if is_readonly_prop(profile, prop):
        return "read-only (gray)"
    if has_do_not_submit_in_prop_comment(profile, prop):
        return "do not submit (gray)"
    if is_not_submittable_prop(profile, prop):
        return "not submittable (gray)"
    return "submittable"


def tooltip_for_prop(profile, prop) -> str:
    """Profile.js makeTooltipForProp."""
    p = _prop(profile, prop)
    tip = "SEARCH AVAILABLE\n\n" if is_searchable_prop(profile, prop) else ""
    tip += "\n\n".join(f"* {k}\n{p[k]}" for k in SELECTED_PROP_KEYS_FOR_TOOLTIP if k in p)
    if "linkTo" in (p.get("items") or {}):
        tip += f"\n\n* linkTo\n{p['items']['linkTo']}"
    ex = p.get("submissionExample") or {}
    if "appscript" in ex:
        tip += f"\n\n* submissionExample (appscript)\n{ex['appscript']}"
    if "igvf_utils" in ex:
        tip += f"\n\n* submissionExample (igvf_utils)\n{ex['igvf_utils']}"
    dep = (profile.get("dependentSchemas") or {}).get(prop) or {}
    if "comment" in dep:
        tip += f"\n\n* dependency\n{dep['comment']}"
    return tip


# ── the sheet (Sheet.js) ────────────────────────────────────────────────────

class Sheet:
    """A header row and data rows, read from / written to TSV, CSV or XLSX.

    Cells hold what Google Sheets would hold: auto-parsed values (see
    auto_parse). `raw` keeps each cell's typed text, for type_cast_value.
    """

    def __init__(self, header=None, rows=None, raw=None, hidden=None, name="Sheet1"):
        self.header = list(header or [])
        self.rows = [list(r) for r in (rows or [])]
        self.raw = [list(r) for r in (raw or self.rows)]
        self.hidden = set(hidden or ())
        self.name = name

    # Sheet.js findColumnByHeaderValue (0-based here, None when absent)
    def col(self, prop):
        return self.header.index(prop) if prop in self.header else None

    def cell(self, r, prop):
        c = self.col(prop)
        if c is None or c >= len(self.rows[r]):
            return ""
        return self.rows[r][c]

    def raw_cell(self, r, prop):
        c = self.col(prop)
        if c is None or r >= len(self.raw) or c >= len(self.raw[r]):
            return ""
        return self.raw[r][c]

    def row_to_json(self, r, keep_commented=True) -> dict:
        """Sheet.js rowToJson: empty cells dropped, `[...]`/`{...}` text
        parsed as JSON, # columns dropped unless keep_commented."""
        out = {}
        for prop in self.header:
            if not prop:
                continue
            if prop.startswith("#") and not keep_commented:
                continue
            val = self.cell(r, prop)
            if val == "" or val is None:
                continue
            if isinstance(val, str) and (_is_json_string(val) or _is_array_string(val)):
                try:
                    val = json.loads(val)
                except ValueError as exc:
                    raise ValueError(f"column {prop!r}: not valid JSON ({exc.msg})") from None
            out[prop] = val
        return out

    def raw_row(self, r) -> dict:
        return {p: self.raw_cell(r, p) for p in self.header if p}

    def update_header(self, props) -> list:
        """Sheet.js updateHeaderWithArray: keep the header, append new props."""
        for p in props:
            if p not in self.header:
                self.header.append(p)
        for r in range(len(self.rows)):
            self.rows[r] += [""] * (len(self.header) - len(self.rows[r]))
            self.raw[r] += [""] * (len(self.header) - len(self.raw[r]))
        return self.header

    def write_json_to_row(self, obj: dict, r: int, props=None) -> None:
        """Sheet.js writeJsonToRow: the whole row is rewritten from `obj`;
        arrays and objects as JSON, null as empty, absent props empty."""
        self.update_header(props if props is not None else list(obj))
        cells = []
        for p in self.header:
            if p in obj:
                v = obj[p]
                if isinstance(v, (list, dict)):
                    v = json.dumps(v, separators=(",", ":"), ensure_ascii=False)
                elif v is None:
                    v = ""
                cells.append(v)
            else:
                cells.append("")
        while r >= len(self.rows):
            self.rows.append([""] * len(self.header))
            self.raw.append([""] * len(self.header))
        self.rows[r] = cells
        self.raw[r] = [_to_text(v) for v in cells]

    def add_json(self, obj, props=None):
        """Sheet.js addJsonToSheet."""
        self.write_json_to_row(obj, len(self.rows), props)

    def is_skipped(self, r) -> bool:
        """Hidden rows (XLSX) and #skip = 1 rows are not read or written."""
        return r in self.hidden or (SKIP in self.header and to_boolean(self.cell(r, SKIP)))

    # ── I/O ──
    @classmethod
    def read(cls, path, sheet_name=None) -> "Sheet":
        p = Path(_may_read(path, "sheet"))
        if p.suffix.lower() in (".xlsx", ".xlsm"):
            try:
                import openpyxl
            except ImportError:
                die("reading .xlsx needs openpyxl (pip install openpyxl), or save the sheet as TSV")
            wb = openpyxl.load_workbook(p, data_only=True)
            ws = wb[sheet_name] if sheet_name else wb.worksheets[0]
            data = [list(r) for r in ws.iter_rows(values_only=True)]
            hidden = {i - 2 for i, d in ws.row_dimensions.items() if d.hidden and i >= 2}
            if not data:
                return cls(name=ws.title)
            header = _trim([("" if h is None else str(h)) for h in data[0]])
            rows, raw = [], []
            for row in data[1:]:
                vals = [("" if v is None else v) for v in row][:len(header)]
                vals += [""] * (len(header) - len(vals))
                rows.append([auto_parse(v) for v in vals])
                raw.append([_to_text(v) for v in vals])
            return cls(header, rows, raw, hidden, ws.title)
        delim = "," if p.suffix.lower() == ".csv" else "\t"
        with open(p, newline="", encoding="utf-8-sig") as fh:
            data = list(csv.reader(fh, delimiter=delim))
        if not data:
            return cls()
        header = _trim(data[0])
        rows, raw = [], []
        for row in data[1:]:
            vals = row[:len(header)] + [""] * (len(header) - len(row))
            rows.append([auto_parse(v) for v in vals])
            raw.append(list(vals))
        return cls(header, rows, raw, set(), p.stem)

    def write(self, path) -> Path:
        p = Path(_may_write(path, "sheet"))
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.suffix.lower() in (".xlsx", ".xlsm"):
            import openpyxl
            wb = openpyxl.Workbook()
            ws = wb.active
            ws.title = self.name or "Sheet1"
            ws.append(self.header)
            for i, row in enumerate(self.rows):
                ws.append([_to_text(v) if isinstance(v, bool) else v for v in row])
                if i in self.hidden:
                    ws.row_dimensions[i + 2].hidden = True
            wb.save(p)
            return p
        delim = "," if p.suffix.lower() == ".csv" else "\t"
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter=delim, lineterminator="\n")
            w.writerow(self.header)
            for row in self.rows:
                w.writerow([_to_text(v) for v in row])
        return p


def _to_text(v) -> str:
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    if v is None:
        return ""
    if isinstance(v, (list, dict)):
        return json.dumps(v, separators=(",", ":"), ensure_ascii=False)
    return str(v)


def _trim(header) -> list:
    """Sheet.js getCellValuesInRow: up to the last non-empty header cell."""
    header = [str(h).strip() for h in header]
    while header and not header[-1]:
        header.pop()
    return header


# ── rows against the Portal (Metadata.js) ───────────────────────────────────

def find_identifying(sheet: Sheet, r: int, profile) -> "tuple":
    """Metadata.js findIdentifyingPropValColInRow: accession, then uuid, then
    name, then the profile's other identifying properties; an array one
    (aliases) gives its first element."""
    props = sorted(profile.get("identifyingProperties") or [], key=get_identifying_prop_priority)
    for prop in props:
        val = sheet.cell(r, prop)
        if val not in ("", None):
            if is_array_prop(profile, prop):
                if isinstance(val, str):
                    try:
                        val = json.loads(val)
                    except ValueError:
                        val = [val]
                val = val[0] if val else ""
                if not val:
                    continue
            return prop, str(val), sheet.col(prop)
    return None, None, None


def set_attachment(att) -> "dict | None":
    """Metadata.js setAttachment, with a local file in place of a Drive path."""
    import base64
    path = (att or {}).get("path") if isinstance(att, dict) else None
    if not path:
        raise ValueError('attachment must be {"path": "<local file>"}')
    f = Path(_may_read(str(Path(path).expanduser()), "attachment"))
    if not f.is_file():
        raise ValueError(f"attachment {path} not found")
    mime = mimetypes.guess_type(f.name)[0] or "application/octet-stream"
    if f.name.endswith(".gz"):
        mime = "application/gzip"          # the add-on maps application/x-gzip to this
    return {"download": f.name, "type": mime,
            "href": f"data:{mime};base64," + base64.b64encode(f.read_bytes()).decode()}


def object_url(profile_name, ident) -> str:
    """Metadata.js makeMetadataUrl (path part): /<profile>/<id>/."""
    return f"/{profile_name}/{ident}/"


def _http(portal: Portal, method, path, payload=None, params=None) -> "tuple[int, dict]":
    """One request; (status, JSON body), errors included, like muteHttpExceptions."""
    try:
        body = portal._request(method, path, params=params, payload=payload)
        return (201 if method == "POST" else 200), body
    except PortalError as exc:
        text = str(exc)
        m = re.search(r"HTTP (\d+)", text)
        code = int(m.group(1)) if m else 0
        detail = text.split("\n", 1)[1].strip() if "\n" in text else text
        try:
            return code, json.loads(detail.split("\n  Hint:")[0])
        except ValueError:
            return code, {"error": detail[:500]}


def get_metadata(portal, profile, profile_name, prop, ident, for_admin=False) -> dict:
    """Metadata.js getMetadataFromPortal: GET frame=object, gettable props only."""
    code, body = _http(portal, "GET", object_url(profile_name, ident), params={"frame": "object"})
    obj = {RESPONSE: f"GET,{code}", RESPONSE_TIME: now_string(), prop: ident}
    if code == 200:
        obj.update({k: v for k, v in body.items() if is_gettable_prop(profile, k, for_admin)})
    else:
        obj[RESPONSE] += "\n" + json.dumps(body, indent=2)
    return obj


# ── commands ────────────────────────────────────────────────────────────────

def _sheet_portal(portal: Portal, args) -> Portal:
    """The add-on's default endpoint is staging (Endpoint.js DEFAULT_ENDPOINT_WRITE)."""
    if args.mode is None and not os.environ.get("IGVF_MODE"):
        print(info("no mode given; using staging, as the metadata submitter does (-m prod for production)"))
        return Portal(mode="staging")
    return portal


def _profile_name(args) -> str:
    """The sheet's profile: --profile, else the file name's first part."""
    name = args.profile or Path(args.sheet).name.split(".")[0]
    return camel_to_snake(name) if name[:1].isupper() else name


def _load_profile(portal, name) -> dict:
    """Profile.js getProfile (the add-on's regex-escape fix is for JavaScript's
    RegExp; Python reads those patterns as written)."""
    try:
        prof = portal.profile(name)
    except PortalError as exc:
        die(f"no profile `{name}` on {portal.mode} (snake_case like tabular_file, or "
            f"CamelCase like TabularFile; never plural): {str(exc)[:160]}")
    if not prof.get("properties"):
        die(f"`{name}` on {portal.mode} is not a profile")
    return prof


def _out_path(args, action) -> Path:
    if getattr(args, "out", None):
        return Path(args.out)
    src = Path(args.sheet)
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    ext = "".join(src.suffixes[-1:]) or ".tsv"
    return src.with_name(f"{src.name[:-len(ext)] if ext else src.name}.{action}.{stamp}{ext}")


def _guard_write(portal, args):
    if args.apply and portal.mode == "prod" and not args.yes:
        die("refusing to write to PRODUCTION without --yes (rehearse with -m staging first)")


def write_legend(profile, props, path) -> Path:
    """The add-on's header colours, fonts, tooltips and dropdowns, as a table."""
    path = Path(_may_write(path, "legend"))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["property", "role", "type", "array", "searchable (links to)", "allowed values", "notes"])
        for p in props:
            if p.startswith("#"):
                w.writerow([p, "commented (not sent)", "", "", "", "", ""])
                continue
            q = _prop(profile, p)
            w.writerow([p, color_for_prop(profile, p), q.get("type", ""),
                        "yes" if is_array_prop(profile, p) else "",
                        link_to(profile, p) or "",
                        ", ".join(map(str, q.get("enum") or (q.get("items") or {}).get("enum") or [])),
                        tooltip_for_prop(profile, p).replace("\n", " ")])
    return path


def _skill():
    try:
        from . import igvf_submission_skill as sk  # type: ignore
    except ImportError:
        import igvf_submission_skill as sk  # type: ignore
    return sk


def _local_path(path: str, also: "Optional[str]" = None) -> str:
    """See igvf_submission_skill._local: a relative path is also looked for under the working directory and
    under ``also``."""
    try:
        return _skill()._local(path, also)
    except Exception:
        return path


def _may_read(path, what: str = "file") -> str:
    """See igvf_submission_skill._may_read: on a shared deployment, only the asker's own files."""
    return _skill()._may_read(path, what)


def _may_write(path, what: str = "output") -> str:
    return _skill()._may_write(path, what)


def cmd_sheet_template(portal: Portal, args) -> int:
    """MetadataTemplate.js addMetadataTemplateToSheet: header + one default row."""
    portal = _sheet_portal(portal, args)
    name = camel_to_snake(args.profile) if args.profile[:1].isupper() else args.profile
    profile = _load_profile(portal, name)
    tpl = make_template(profile, args.admin)
    props = get_sorted_props(list(tpl), profile)
    sh = Sheet(name=name)
    sh.add_json(tpl, props)
    if args.out:
        out = Path(args.out)
    else:   # not the current directory: the agent runs tools from the project root
        out = Path(os.environ.get("IGVF_PROJECT_ROOT") or ".") / "Data" / "Submitter" / "sheets" / f"{name}.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
    sh.write(out)
    legend = write_legend(profile, props, out.with_name(out.name.split(".")[0] + ".legend.tsv"))
    print(ok(f"{name} (schema version {schema_version(profile)}): {len(props)} columns, "
             f"{sum(is_required_prop(profile, p) for p in props)} required"))
    print(f"Wrote: {out}")
    print(f"Wrote: {legend}")
    return 0


def cmd_sheet_get(portal: Portal, args) -> int:
    """Metadata.js updateSheetWithMetadataFromPortal."""
    portal = _sheet_portal(portal, args)
    name = _profile_name(args)
    profile = _load_profile(portal, name)
    sh = Sheet.read(args.sheet, args.sheet_name)
    n = 0
    for r in range(len(sh.rows)):
        if sh.is_skipped(r):
            continue
        prop, ident, _ = find_identifying(sh, r, profile)
        if not ident:
            continue
        obj = get_metadata(portal, profile, name, prop, ident, args.admin)
        sh.write_json_to_row(obj, r, get_sorted_props(list(obj), profile))
        print((ok if obj[RESPONSE].startswith("GET,200") else fail)(f"row {r + 2}: {prop}={ident} "
                                                                   f"{obj[RESPONSE].splitlines()[0]}"))
        n += 1
    out = sh.write(_out_path(args, "get"))
    print(info(f"{n} row(s) fetched from {portal.mode}"))
    print(f"Wrote: {out}")
    return 0


def _row_payload(sh, r, profile) -> "tuple[dict, dict]":
    """(row JSON with # columns, typed payload without them), Metadata.js."""
    before = sh.row_to_json(r, keep_commented=True)
    payload = type_cast_json(profile, before, sh.raw_row(r))
    if has_attachment(profile) and payload.get("attachment"):
        payload["attachment"] = set_attachment(payload["attachment"])
    return before, payload


def _prop_errors(profile, payload) -> "list[str]":
    """Each sent property against its own schema (a PATCH carries a subset, so
    whole-object `required` errors do not apply)."""
    mini = {"type": "object", "properties": {k: _prop(profile, k) for k in payload
                                             if k in (profile.get("properties") or {})}}
    unknown = [f"{k}: not a property of {profile.get('title', 'this profile')}"
               for k in payload if k not in (profile.get("properties") or {})]
    return unknown + _schema_errors(mini, payload)


def _submit(portal: Portal, args, method) -> int:
    """Metadata.js submitSheetToPortal for POST, PATCH, PATCH_REMOVE and PUT."""
    portal = _sheet_portal(portal, args)
    _guard_write(portal, args)
    name = _profile_name(args)
    profile = _load_profile(portal, name)
    sh = Sheet.read(args.sheet, args.sheet_name)
    cols = [c.strip() for c in (getattr(args, "columns", "") or "").split(",") if c.strip()]
    remove = [c.strip() for c in (getattr(args, "remove", "") or "").split(",") if c.strip()]
    if method == "PATCH" and remove:
        method = "PATCH_REMOVE"
    if method == "PATCH" and not cols and not getattr(args, "all", False):
        die("sheet-patch needs --columns a,b (PATCH selected columns), --all (PATCH all "
            "columns) or --remove x,y (PATCH with property removal)")
    n = problems = 0
    for r in range(len(sh.rows)):
        if sh.is_skipped(r):
            continue
        try:
            before, json_ = _row_payload(sh, r, profile)
        except ValueError as exc:
            sh.write_json_to_row({**_safe_row(sh, r), RESPONSE: f"{method},error\n{exc}",
                                  RESPONSE_TIME: now_string()}, r)
            print(fail(f"row {r + 2}: {exc}"))
            problems += 1
            continue
        if method in ("PATCH", "PATCH_REMOVE") and cols:
            payload = {k: v for k, v in json_.items() if k in cols}
        elif method == "PATCH_REMOVE":
            payload = {}
        else:
            payload = dict(json_)
        if (method == "PATCH_REMOVE" and not payload and not remove) or \
                (method != "PATCH_REMOVE" and not payload):
            continue
        ident = prop = None
        if method != "POST":
            prop, ident, _ = find_identifying(sh, r, profile)
            if not ident:
                continue
        label = f"row {r + 2}" + (f" {prop}={ident}" if ident else "")
        tail = ""
        if method == "PATCH":
            tail = "\nSelected props: " + (",".join(cols) if cols else "ALL")
        elif method == "PATCH_REMOVE":
            tail = ("\nSelected props: " + (",".join(cols) if cols else "NONE")
                    + "\nRemove props: " + (",".join(remove) if remove else "NONE"))

        if not args.apply:
            errs = _schema_errors(profile, payload) if method in ("POST", "PUT") \
                else _prop_errors(profile, payload)
            problems += bool(errs)
            before[RESPONSE] = f"DRYRUN {method}{tail}\n" + ("ValidationSuccess" if not errs
                                                             else json.dumps(errs, indent=2))
            before[RESPONSE_TIME] = now_string()
            sh.write_json_to_row(before, r)
            print((fail if errs else ok)(f"{label}: {method} would send {len(payload)} "
                                         f"propert{'y' if len(payload) == 1 else 'ies'}"
                                         + (f" and remove {', '.join(remove)}" if remove else "")
                                         + (f"; {len(errs)} validation error(s)" if errs else "")))
            if method == "POST" and (json_.get("accession") or json_.get("uuid")):
                print(warn(f"       {label} already has an accession/uuid; POSTing it again "
                           f"makes a second object or a 409 conflict -- PATCH it instead"))
            for e in errs[:5]:
                print("       " + e)
            n += 1
            continue

        if method == "POST":
            code, body = _http(portal, "POST", f"/{name}/", payload)
        elif method == "PATCH_REMOVE":
            code, body = _http(portal, "GET", object_url(name, ident), params={"frame": "edit"})
            if code == 200:
                edited = dict(body)
                for k in remove:
                    edited.pop(k, None)
                edited.update(payload)
                code, body = _http(portal, "PUT", object_url(name, ident), edited)
        else:
            code, body = _http(portal, method, object_url(name, ident), payload)

        before[RESPONSE] = f"{method},{code}{tail}"
        before[RESPONSE_TIME] = now_string()
        if code == 201:
            got = (body.get("@graph") or [body])[0]
            for p in profile.get("identifyingProperties") or []:
                if not is_commented_prop(p):
                    before[p] = got.get(p)
            before[RESPONSE] += "\n" + json.dumps(body, indent=2)
        elif code != 200:
            if code == 422:
                before[RESPONSE] += ("\nIf error message is not helpful, run the same command "
                                     "without --apply to validate.\n")
            before[RESPONSE] += "\n" + json.dumps(body, indent=2)
            problems += 1
        sh.write_json_to_row(before, r)
        print((ok if code in (200, 201) else fail)(
            f"{label}: {method},{code}" + (f" -> {before.get('accession') or before.get('uuid')}"
                                           if code == 201 else "")))
        n += 1
    out = sh.write(_out_path(args, method.lower() + ("" if args.apply else "-dryrun")))
    if not args.apply:
        verb = {"POST": "POST", "PUT": "PUT", "PATCH": "PATCH",
                "PATCH_REMOVE": "PATCH with property removal"}[method]
        print(warn(f"dry run -- nothing sent ({n} row(s) checked). Re-run with --apply to "
                   f"{verb} on {portal.mode}."))
        if method in ("PATCH", "PATCH_REMOVE") and n and not problems:
            _save_patch_plan(portal, args, name, cols, remove, n)
    print(f"Wrote: {out}")
    return 1 if problems else 0


def _save_patch_plan(portal, args, profile: str, cols, remove, n: int) -> Path:
    """A clean PATCH dry run, kept so that what was checked is what gets sent: a copy of the sheet
    as validated, with the options, under Data/Submitter (the web app's Submit button runs plans
    from there only, as it does create's `Submission plan`)."""
    import shutil
    import time
    root = Path(os.environ.get("IGVF_PROJECT_ROOT") or ".")
    plan_dir = root / "Data" / "Submitter" / f"{time.strftime('%Y%m%d_%H%M%S')}_patch_{portal.mode}"
    plan_dir.mkdir(parents=True, exist_ok=True)
    sheet = plan_dir / Path(args.sheet).name
    shutil.copy2(args.sheet, sheet)
    plan = {"command": "sheet-patch", "sheet": str(sheet), "profile": profile,
            "sheet_name": getattr(args, "sheet_name", None), "columns": cols,
            "all": bool(getattr(args, "all", False)), "remove": remove,
            "admin": bool(getattr(args, "admin", False)), "mode": portal.mode, "rows": n}
    (plan_dir / "patch.json").write_text(json.dumps(plan, indent=1))
    try:
        from . import igvf_submission_skill as sk  # type: ignore
    except ImportError:
        import igvf_submission_skill as sk  # type: ignore
    sk.save_plan_meta(plan_dir / "patch.json", "patch", n, portal.mode)
    print(f"Patch plan: {plan_dir / 'patch.json'}  ({n} object(s), {portal.mode})")
    print(f"To send exactly this change with your own Portal keys: igvfagent submit send {plan_dir / 'patch.json'}"
          + (" --yes" if portal.mode == "prod" else "") + "  (on the web app: its Apply button)")
    return plan_dir / "patch.json"


def patch_plan_command(plan_path, mode: str) -> "list[str] | str":
    """The `submit sheet-patch ... --apply` arguments (after `igvfagent`) that send a saved plan, or
    why it cannot be sent: the sheet must be the copy the dry run checked, in the plan's own folder,
    and the plan must be for the Portal it is being sent to."""
    try:
        plan = json.loads(Path(plan_path).read_text())
        Path(plan["sheet"]).resolve().relative_to(Path(plan_path).resolve().parent)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return f"The patch plan cannot be read or does not hold its own sheet ({exc})."
    if plan.get("mode") != mode:
        return f"The patch plan was checked on {plan.get('mode')}, not {mode}."
    return patch_plan_argv(plan)


def patch_plan_argv(plan: dict) -> "list[str]":
    """The `submit sheet-patch` arguments that send a saved plan (after `igvfagent`)."""
    argv = ["submit", "sheet-patch", plan["sheet"], "--profile", plan["profile"], "-m", plan["mode"], "--apply"]
    if plan.get("sheet_name"):
        argv += ["--sheet-name", plan["sheet_name"]]
    if plan.get("columns"):
        argv += ["--columns", ",".join(plan["columns"])]
    if plan.get("all"):
        argv.append("--all")
    if plan.get("remove"):
        argv += ["--remove", ",".join(plan["remove"])]
    if plan.get("admin"):
        argv.append("--admin")
    if plan["mode"] == "prod":
        argv.append("--yes")
    return argv


def _safe_row(sh, r) -> dict:
    """A row's cells as they are, for writing an error back without parsing."""
    return {p: v for p, v in zip(sh.header, sh.rows[r]) if p and v not in ("", None)}


def cmd_sheet_post(portal, args) -> int:
    return _submit(portal, args, "POST")


def cmd_sheet_patch(portal, args) -> int:
    return _submit(portal, args, "PATCH")


def cmd_sheet_put(portal, args) -> int:
    if not args.admin:
        die("PUT replaces the whole object on the Portal and removes every property the row "
            "leaves out; it is for DACC admins. Pass --admin to confirm, and GET the rows with "
            "`sheet-get --admin` first.")
    return _submit(portal, args, "PUT")


def cmd_sheet_upload(portal: Portal, args) -> int:
    """Upload.js initUpload (the sidebar) and generateS3UploadCmd."""
    portal = _sheet_portal(portal, args)
    if args.apply and not args.print_cmd:
        _guard_write(portal, args)
    name = _profile_name(args)
    profile = _load_profile(portal, name)
    sh = Sheet.read(args.sheet, args.sheet_name)
    if not any(sh.col(p) is not None for p in profile.get("identifyingProperties") or []):
        die(f"no identifying column ({', '.join(profile.get('identifyingProperties') or [])}) "
            f"in the header; add one, filled for each row")
    need = UPLOAD_CMD if args.print_cmd else UPLOAD_STATUS
    if sh.col(UPLOAD_ABSPATH) is None or sh.col(need) is None:
        # Upload.js inserts the two columns at the left when they are missing.
        sh.header = [c for c in (UPLOAD_ABSPATH, need) if sh.col(c) is None] + sh.header
        for r in range(len(sh.rows)):
            pad = [""] * (len(sh.header) - len(sh.rows[r]))
            sh.rows[r] = pad + sh.rows[r]
            sh.raw[r] = pad + sh.raw[r]
        out = sh.write(_out_path(args, "upload-template"))
        print(warn(f"added {UPLOAD_ABSPATH} and {need}; fill {UPLOAD_ABSPATH} with each file's "
                   f"absolute path and run again on {out.name}"))
        print(f"Wrote: {out}")
        return 1
    n = problems = 0
    for r in range(len(sh.rows)):
        if sh.is_skipped(r):
            continue
        row = sh.row_to_json(r, keep_commented=True)
        prop, ident, _ = find_identifying(sh, r, profile)
        if not ident:
            row[UPLOAD_STATUS] = "Missing identifying value (e.g. accession, uuid)."
            sh.write_json_to_row(row, r)
            problems += 1
            continue
        code, rec = _http(portal, "GET", f"/files/{ident}/", params={"frame": "object"})
        status = rec.get("status") if code == 200 else None
        cerr = rec.get("content_error_detail") if status == "content error" else None
        if status == "released":
            row[UPLOAD_STATUS] = "error: uploading to a released accession is not allowed."
            sh.write_json_to_row(row, r)
            problems += 1
            continue
        path = row.get(UPLOAD_ABSPATH)
        if not path:
            row[UPLOAD_STATUS] = f"Missing {UPLOAD_ABSPATH}."
            sh.write_json_to_row(row, r)
            problems += 1
            continue
        # A relative entry is looked for from the working directory and from the sheet's own folder too.
        path = _may_read(_local_path(str(path), str(Path(args.sheet).resolve().parent)), "upload file")
        if not args.apply and not args.print_cmd:
            exists = Path(str(path)).is_file()
            row[UPLOAD_STATUS] = (f"status: {status}" + (f": {cerr}" if cerr else "") + "."
                                  + ("" if exists else f" local file not found: {path}"))
            sh.write_json_to_row(row, r)
            print((ok if exists else fail)(f"row {r + 2}: {ident} ({status}) <- {path}"))
            problems += not exists
            n += 1
            continue
        c2, up = _http(portal, "GET", f"/files/{ident}/@@upload", params={"frame": "object"})
        creds = ((up.get("@graph") or [{}])[0] or {}).get("upload_credentials") if c2 == 200 else None
        if not creds:
            # Credentials expire; igvf_utils asks for fresh ones with a POST.
            c3, up = _http(portal, "POST", f"/files/{ident}/@@upload", {})
            creds = ((up.get("@graph") or [{}])[0] or {}).get("upload_credentials") if c3 in (200, 201) else None
        if not creds:
            row[UPLOAD_STATUS] = "Failed to get upload credentials."
            sh.write_json_to_row(row, r)
            problems += 1
            continue
        row[UPLOAD_STATUS] = f"status: {status}" + (f": {cerr}" if cerr else "") + "."
        if args.print_cmd:
            bucket = creds["upload_url"].split("/")[2]
            key = "/".join(creds["upload_url"].split("/")[3:])
            row[UPLOAD_CMD] = (f'AWS_ACCESS_KEY_ID="{creds["access_key"]}" AWS_SECRET_ACCESS_KEY='
                               f'"{creds["secret_key"]}" AWS_SESSION_TOKEN="{creds["session_token"]}" '
                               f'aws s3api put-object --bucket "{bucket}" --key "{key}" --body "{path}"')
            print(ok(f"row {r + 2}: {ident} upload command written"))
        else:
            try:
                how = _upload_bytes(creds, str(path))
                row[UPLOAD_STATUS] += f" uploaded via {how}."
                print(ok(f"row {r + 2}: {ident} <- {Path(str(path)).name} uploaded via {how}"))
            except (PortalError, OSError) as exc:
                row[UPLOAD_STATUS] += f" upload failed: {exc}"
                print(fail(f"row {r + 2}: {ident}: {exc}"))
                problems += 1
        sh.write_json_to_row(row, r)
        n += 1
    out = sh.write(_out_path(args, "upload-cmd" if args.print_cmd else
                             ("upload" if args.apply else "upload-dryrun")))
    if args.print_cmd:
        print(warn("the commands carry temporary AWS credentials; do not share the sheet"))
    elif not args.apply:
        print(warn("dry run -- nothing uploaded. Re-run with --apply to upload, or --print-cmd "
                   "for `aws s3api` command lines"))
    print(f"Wrote: {out}")
    return 1 if problems else 0


# ── CLI ─────────────────────────────────────────────────────────────────────

def register(sub) -> None:
    """Add the sheet-* subcommands to `igvfagent submit`."""
    def common(p, writes=True):
        p.add_argument("sheet", help="TSV, CSV or XLSX sheet, one row per object")
        p.add_argument("--profile", help="profile (default: the file name's first part, e.g. "
                                         "tabular_file.tsv)")
        p.add_argument("--sheet-name", help="worksheet of an XLSX file (default: the first)")
        p.add_argument("--out", help="output sheet (default: next to the input, timestamped)")
        p.add_argument("--admin", action="store_true",
                       help="include admin-only properties (the add-on's ADMIN actions)")
        if writes:
            p.add_argument("--apply", action="store_true", help="actually send")
            p.add_argument("--yes", action="store_true", help="required to write to prod")

    p = sub.add_parser("sheet-template", help="metadata-submitter sheet for a profile: header + default row")
    p.add_argument("profile")
    p.add_argument("--out")
    p.add_argument("--admin", action="store_true")
    p.set_defaults(fn=cmd_sheet_template)

    p = sub.add_parser("sheet-get", help="fill each row from the Portal (GET by accession/uuid/alias)")
    common(p, writes=False)
    p.set_defaults(fn=cmd_sheet_get)

    p = sub.add_parser("sheet-post", help="POST each row; new accessions written back (dry run unless --apply)")
    common(p)
    p.set_defaults(fn=cmd_sheet_post)

    p = sub.add_parser("sheet-patch", help="PATCH rows: --columns, --all or --remove (dry run unless --apply)")
    common(p)
    p.add_argument("--columns", help="comma-separated columns to send (PATCH selected columns)")
    p.add_argument("--all", action="store_true", help="send every non-empty column (PATCH all columns)")
    p.add_argument("--remove", help="comma-separated properties to delete from the Portal object "
                                    "(GET frame=edit, delete, PUT); combine with --columns")
    p.set_defaults(fn=cmd_sheet_patch)

    p = sub.add_parser("sheet-put", help="PUT each row, replacing the object (admin only, needs --admin)")
    common(p)
    p.set_defaults(fn=cmd_sheet_put)

    p = sub.add_parser("sheet-upload", help="upload files named in #upload_abspath, or --print-cmd")
    common(p)
    p.add_argument("--print-cmd", action="store_true",
                   help="write `aws s3api put-object` commands to #upload_cmd instead of uploading")
    p.set_defaults(fn=cmd_sheet_upload)


if __name__ == "__main__":
    try:
        from .igvf_submission_skill import main
    except ImportError:
        from igvf_submission_skill import main  # type: ignore
    sys.exit(main())
