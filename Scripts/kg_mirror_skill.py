#!/usr/bin/env python3
"""IGVF Knowledge Graph local mirror skill.

Streams each Arango collection (except a configurable skip list) via the
read-only AQL cursor API and persists it as a Parquet shard on disk,
then registers the result as views in a DuckDB warehouse for offline
querying.

Why this exists:
    The full IGVF Catalog KG is ~2 TB on disk inside Arango. The two
    "planet-scale" collections — ``variants`` (~944 GB) and
    ``variants_variants`` (~531 GB) — would dominate any local mirror,
    and on this workstation we have ~1.3 TB free. By default this skill
    skips those two and mirrors everything else — roughly 470 GB
    uncompressed in Arango, compressing to ~80-120 GB on disk as ZSTD-
    compressed Parquet shards. The smaller collections download in
    seconds; the medium ones (genomic_elements, mm_variants) take
    minutes; only the coding_variants family takes hours.

Commands:
    kg-mirror inventory       — list collections with doc counts + bytes
    kg-mirror pull            — mirror a single collection
    kg-mirror pull-all        — mirror everything except the skip list
    kg-mirror register        — register parquet shards as DuckDB views
    kg-mirror verify          — sanity check the local mirror
    kg-mirror write-playbook  — emit the skill's markdown playbook

Storage layout:
    Data/Warehouse/KG/<collection>/<NNNN>.parquet  — per-batch shards
    Data/Warehouse/KG/_state/<collection>.json     — resume cursor
    Data/Warehouse/igvf_kg_mirror.duckdb           — DuckDB warehouse
                                                     with one view per
                                                     collection.

License: Apache-2.0. Uses only stdlib + pyarrow + duckdb (Apache-2 / BSD).
"""

from __future__ import annotations

import argparse
import base64
import csv
import json
import logging
import os
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
WAREHOUSE_DIR = DATA_DIR / "Warehouse"
KG_DIR = WAREHOUSE_DIR / "KG"
STATE_DIR = KG_DIR / "_state"
SKILL_DOC_DIR = DOCS_DIR / "Skills"
DUCKDB_PATH = WAREHOUSE_DIR / "igvf_kg_mirror.duckdb"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _endpoints import resolve as _resolve_endpoint  # noqa: E402

ARANGO_BASE = _resolve_endpoint("arango", "IGVF_ARANGO_BASE")

# Pathological collections that AQL cursor cannot stream cleanly.
# - genes_coding_variants_scores      ~64 GB  / 69 K docs   (bimodal: most
#                                        docs are 127 B but a few embed
#                                        ~80 MB per-variant score matrices
#                                        that time out the AQL cursor)
# - genes_coding_variants_scores_grp  ~50 GB  / 69 K docs   (same problem)
# Mirror these per-gene through the Catalog REST API instead.
# NOTE: `variants` (~944 GB / 1.87 B docs) and `variants_variants`
# (~531 GB / 5.93 B edges) ARE included by default — they are the heavy
# tail of the graph and the mirror is supposed to be complete. Expect
# multi-day pulls for those; the skill is resumable.
DEFAULT_SKIP = (
    "genes_coding_variants_scores",
    "genes_coding_variants_scores_grp",
)


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"kg_mirror_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return log_path


def safe_label(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in s)


def _require_pkg(name: str, hint: str) -> Any:
    try:
        return __import__(name)
    except Exception as exc:
        raise SystemExit(f"Missing dependency '{name}'. {hint}\nInstall: pip install {name}") from exc


# ---------------------------------------------------------------------------
# Arango client (basic-auth read-only)
# ---------------------------------------------------------------------------

class ArangoClient:
    def __init__(self, base: str, user: str | None, password: str | None):
        self.base = base.rstrip("/")
        self.user = user
        self.password = password
        self._auth_header: dict[str, str] = {}
        if user and password:
            token = base64.b64encode(f"{user}:{password}".encode()).decode()
            self._auth_header = {"Authorization": f"Basic {token}"}

    def _request(self, path: str, *, method: str = "GET", body: dict | None = None,
                  timeout: int = 120) -> Any:
        url = f"{self.base}{path}"
        data = None
        hdr = {"Accept": "application/json", **self._auth_header}
        if body is not None:
            data = json.dumps(body).encode()
            hdr["Content-Type"] = "application/json"
        req = urllib.request.Request(url, data=data, headers=hdr, method=method)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.load(resp)

    def list_collections(self) -> list[dict[str, Any]]:
        data = self._request("/_api/collection?excludeSystem=true")
        return data.get("result", [])

    def collection_count(self, name: str) -> int:
        return int(self._request(f"/_api/collection/{name}/count").get("count", 0))

    def collection_figures(self, name: str) -> dict[str, Any]:
        return self._request(f"/_api/collection/{name}/figures").get("figures", {})

    def open_cursor(self, query: str, *, batch_size: int = 5000,
                     bind_vars: dict | None = None) -> dict[str, Any]:
        body: dict[str, Any] = {
            "query": query,
            "batchSize": batch_size,
            "count": True,
        }
        if bind_vars:
            body["bindVars"] = bind_vars
        return self._request("/_api/cursor", method="POST", body=body)

    def next_batch(self, cursor_id: str) -> dict[str, Any]:
        return self._request(f"/_api/cursor/{cursor_id}", method="PUT")


def make_client() -> ArangoClient:
    return ArangoClient(
        ARANGO_BASE,
        os.environ.get("IGVF_ARANGO_USER", "guest"),
        os.environ.get("IGVF_ARANGO_PASSWORD"),
    )


# ---------------------------------------------------------------------------
# Inventory
# ---------------------------------------------------------------------------

def _collection_type_label(c: dict[str, Any]) -> str:
    return "edge" if c.get("type") == 3 else "doc"


def cmd_inventory(args: argparse.Namespace) -> int:
    setup_logging()
    client = make_client()
    cols = client.list_collections()
    rows = []
    for c in cols:
        name = c["name"]
        try:
            count = client.collection_count(name)
        except Exception:
            count = None
        try:
            fig = client.collection_figures(name)
            size = int(fig.get("documentsSize", 0))
        except Exception:
            size = 0
        rows.append({
            "type": _collection_type_label(c),
            "collection": name,
            "documents": count,
            "bytes": size,
            "skip_default": name in DEFAULT_SKIP,
        })
    rows.sort(key=lambda r: -(r["documents"] or 0))
    KG_DIR.mkdir(parents=True, exist_ok=True)
    out = KG_DIR / "_inventory.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    total_docs = sum((r["documents"] or 0) for r in rows)
    total_bytes = sum((r["bytes"] or 0) for r in rows)
    kept_bytes = sum((r["bytes"] or 0) for r in rows if not r["skip_default"])
    print(f"Inventory: {out}")
    print(f"{len(rows)} collections | {total_docs:,} docs | "
          f"{total_bytes/1e9:.1f} GB raw")
    print(f"Excluding default skip ({', '.join(DEFAULT_SKIP)}): "
          f"{kept_bytes/1e9:.1f} GB raw to mirror.")
    return 0


# ---------------------------------------------------------------------------
# Pull a single collection via streaming AQL cursor
# ---------------------------------------------------------------------------

def _state_path(collection: str) -> Path:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    return STATE_DIR / f"{collection}.json"


def _read_state(collection: str) -> dict[str, Any]:
    p = _state_path(collection)
    if not p.exists():
        return {"collection": collection, "rows_written": 0, "batches": 0,
                "last_completed_batch": -1, "status": "new"}
    return json.loads(p.read_text())


def _write_state(state: dict[str, Any]) -> None:
    _state_path(state["collection"]).write_text(json.dumps(state, indent=2))


def _collection_dir(collection: str) -> Path:
    p = KG_DIR / collection
    p.mkdir(parents=True, exist_ok=True)
    return p


def _flatten_row(row: Any) -> dict[str, Any]:
    """Flatten one Arango document; coerce non-scalar fields to JSON strings.

    Arango docs can mix scalars with nested dicts/lists in arbitrary ways.
    Parquet/pandas needs a consistent schema, so any non-scalar value is
    json.dumps()'d into a string column. The raw JSON is still recoverable
    later with json.loads().
    """
    if not isinstance(row, dict):
        return {"value": json.dumps(row)}
    out: dict[str, Any] = {}
    for k, v in row.items():
        if v is None or isinstance(v, (str, int, float, bool)):
            out[k] = v
        else:
            out[k] = json.dumps(v, sort_keys=True, default=str)
    return out


_INT64_MAX = (1 << 63) - 1
_INT64_MIN = -(1 << 63)


_FLOAT64_EXACT_INT = (1 << 53)  # 2^53 = max integer exactly representable in IEEE 754 double


def _coerce_column(values: list[Any]) -> list[Any]:
    """Pick a Parquet-safe representation for a column.

    Rules (apply in order; demote to string if any triggers):
      1. ANY int overflows int64 (i.e. |v| > 2^63 - 1).
      2. Column mixes (int, float). PyArrow promotes the column to
         float64, which can ONLY exactly represent integers up to
         |v| <= 2^53. If any int is bigger than that, the conversion
         loses precision and Arrow raises ArrowInvalid. So when a
         mixed-numeric column carries any |int| > 2^53, demote to str.
      3. Column mixes a non-numeric type with anything else (str ∪
         numeric, list ∪ scalar, etc.). PyArrow can't union those.
    Otherwise leave the values untouched and let PyArrow infer.
    """
    types: "set[type]" = set()
    has_overflow = False
    has_big_int_with_float = False
    has_float = False
    for v in values:
        if v is None:
            continue
        t = type(v)
        types.add(t)
        if t is int and not (_INT64_MIN <= v <= _INT64_MAX):
            has_overflow = True
        if t is float:
            has_float = True
    # Pre-scan for big ints that would lose precision when promoted to float64
    if has_float:
        for v in values:
            if isinstance(v, int) and not isinstance(v, bool) and abs(v) > _FLOAT64_EXACT_INT:
                has_big_int_with_float = True
                break
    numeric_only = types.issubset({bool, int, float})
    if has_overflow or has_big_int_with_float or (len(types) > 1 and not numeric_only):
        return [None if v is None else (v if isinstance(v, str) else str(v)) for v in values]
    return values


def _write_shard(rows: list[dict[str, Any]], collection: str, batch_idx: int) -> Path:
    pa = _require_pkg("pyarrow", "Required to write Parquet shards.")
    pq = __import__("pyarrow.parquet", fromlist=["parquet"])
    out = _collection_dir(collection) / f"{batch_idx:05d}.parquet"
    # Union schema across all rows so missing fields become nulls
    keys = {k for r in rows for k in r.keys()}
    cols: dict[str, list[Any]] = {k: [r.get(k) for r in rows] for k in keys}
    cols = {k: _coerce_column(v) for k, v in cols.items()}
    table = pa.Table.from_pydict(cols)
    pq.write_table(table, out, compression="zstd")
    return out


def _pick_initial_batch_size(client: ArangoClient, collection: str,
                              user_default: int) -> int:
    """Auto-size the AQL batch based on average doc size.

    Arango times out (HTTP 503/504) if we ask for too many bytes in one go.
    For wide-doc collections like ``genes_coding_variants_scores`` (~1 MB
    per doc), a 5000-row batch is ~5 GB and inevitably times out. Use the
    ratio ``documentsSize / count`` to target a single-request payload of
    at most ~64 MB.
    """
    try:
        fig = client.collection_figures(collection)
        size = int(fig.get("documentsSize", 0))
        count = client.collection_count(collection)
        if count > 0 and size > 0:
            avg = size / count
            target_bytes = 64 * 1024 * 1024  # 64 MB per response
            cap = max(1, int(target_bytes / max(avg, 1)))
            return max(50, min(user_default, cap))
    except Exception:
        pass
    return user_default


def _pull_collection(client: ArangoClient, collection: str, *,
                      batch_size: int = 5000,
                      max_rows: int | None = None,
                      restart: bool = False,
                      max_retries: int = 6,
                      grow_after_n_successes: int = 8) -> dict[str, Any]:
    state = _read_state(collection)
    if restart:
        state = {"collection": collection, "rows_written": 0, "batches": 0,
                  "last_completed_batch": -1, "status": "restart"}
    next_batch = state["last_completed_batch"] + 1
    skip = state["rows_written"]
    # Target batch size is what the user asked for; effective is what we are
    # currently using (may be smaller after a server-side timeout). Elastic
    # recovery: after N consecutive successful batches, double effective
    # until we reach target again.
    target_batch = batch_size
    effective_batch = _pick_initial_batch_size(client, collection, batch_size)
    if effective_batch != batch_size:
        logging.info("Auto-sized batch for %s: %d -> %d (avg doc size constraint)",
                      collection, batch_size, effective_batch)
    # Restore effective_batch from state if smaller than computed, so we
    # don't immediately re-trigger the same 5xx. But never let stale state
    # drop us BELOW the target's auto-sized floor on a fresh start.
    saved = state.get("effective_batch_size")
    if saved and isinstance(saved, int) and saved < effective_batch:
        effective_batch = saved
        logging.info("Resuming with reduced batch size from state: %d (target=%d)",
                      effective_batch, target_batch)
    consecutive_successes = 0
    logging.info("Pulling %s (resuming at row %d, batch %d, "
                  "effective_batch=%d, target=%d)",
                  collection, skip, next_batch, effective_batch, target_batch)

    started = time.time()
    rows_total = 0
    retries = 0
    while True:
        if max_rows is not None and rows_total >= max_rows:
            break
        page_size = effective_batch
        if max_rows is not None:
            page_size = min(page_size, max_rows - rows_total)
        aql = (f"FOR d IN {collection} "
                f"SORT d._key "
                f"LIMIT @skip, @count "
                f"RETURN d")
        try:
            cursor = client.open_cursor(aql, batch_size=page_size,
                                          bind_vars={"skip": skip, "count": page_size})
            retries = 0  # success — reset
        except urllib.error.HTTPError as exc:
            retries += 1
            if retries > max_retries:
                logging.error("Giving up on %s at skip=%d after %d retries (HTTP %s).",
                               collection, skip, max_retries, exc.code)
                state["status"] = f"failed_http_{exc.code}"
                _write_state(state)
                return state
            # Exponential backoff AND shrink the batch on 5xx (likely timeout).
            backoff = min(60.0, 2.0 ** retries)
            if exc.code in (502, 503, 504, 429) and effective_batch > 100:
                effective_batch = max(100, effective_batch // 2)
                consecutive_successes = 0
                logging.warning("HTTP %s for %s skip=%d retry=%d — shrinking batch to %d, sleeping %.1fs",
                                 exc.code, collection, skip, retries, effective_batch, backoff)
            else:
                logging.warning("HTTP %s for %s skip=%d retry=%d — sleeping %.1fs",
                                 exc.code, collection, skip, retries, backoff)
            time.sleep(backoff)
            continue
        except Exception as exc:
            retries += 1
            if retries > max_retries:
                logging.error("Giving up on %s at skip=%d after %d retries: %s",
                               collection, skip, max_retries, exc)
                state["status"] = f"failed_{type(exc).__name__}"
                _write_state(state)
                return state
            backoff = min(60.0, 2.0 ** retries)
            logging.warning("Error pulling %s skip=%d retry=%d: %s — sleeping %.1fs",
                             collection, skip, retries, exc, backoff)
            time.sleep(backoff)
            continue
        rows = cursor.get("result", [])
        if not rows:
            state["status"] = "done"
            _write_state(state)
            break
        flat = [_flatten_row(r) for r in rows]
        shard = _write_shard(flat, collection, next_batch)
        skip += len(rows)
        rows_total += len(rows)
        state["rows_written"] = skip
        state["batches"] += 1
        state["last_completed_batch"] = next_batch
        state["status"] = "in_progress"
        state["last_shard"] = str(shard)
        state["effective_batch_size"] = effective_batch
        _write_state(state)
        next_batch += 1
        elapsed = time.time() - started
        rate = rows_total / max(elapsed, 1e-3)
        logging.info("  batch %d: +%d rows -> %s (cumulative %d, %.1f rows/s)",
                      next_batch - 1, len(rows), shard.name, skip, rate)
        consecutive_successes += 1
        if consecutive_successes >= grow_after_n_successes and effective_batch < target_batch:
            new_batch = min(target_batch, effective_batch * 2)
            if new_batch != effective_batch:
                logging.info("  ↑ grew effective batch %d -> %d (after %d consecutive successes)",
                              effective_batch, new_batch, consecutive_successes)
                effective_batch = new_batch
                consecutive_successes = 0
        if len(rows) < page_size:
            state["status"] = "done"
            _write_state(state)
            break
    elapsed = time.time() - started
    state["last_pull_seconds"] = round(elapsed, 1)
    _write_state(state)
    return state


def cmd_pull(args: argparse.Namespace) -> int:
    setup_logging()
    client = make_client()
    state = _pull_collection(client, args.collection,
                              batch_size=args.batch_size,
                              max_rows=args.max_rows,
                              restart=args.restart)
    print(f"\n{args.collection}: status={state['status']}, "
          f"rows={state['rows_written']:,}, batches={state['batches']}, "
          f"elapsed={state.get('last_pull_seconds', 0):.1f}s")
    return 0


# ---------------------------------------------------------------------------
# Key-range pulls: the billion-row collections, in parallel
# ---------------------------------------------------------------------------
#
# `pull` pages with `SORT d._key LIMIT @skip, @count`. That is correct but
# the server walks every skipped row, so each page costs more than the last
# (15 s per 5,000 rows at 130 M rows in); the multi-billion-row collections
# would take years.
#
# The obvious fix, keyset paging (`FILTER d._key > @last SORT d._key`), LOSES
# ROWS on the IGVF catalog. It is a sharded cluster: a key-range FILTER is
# answered by each shard's primary index, which compares keys byte by byte,
# while the SORT merging the shards uses a collation ("--IcY..." sorts before
# "--IR_J...", yet the filter holds no key in between). A key that sorts
# after a page's last row but is byte-wise smaller than it is excluded from
# every later page.
#
# A key-range filter ON ITS OWN is consistent: [lo, hi) splits partition a
# collection exactly (counts on both sides of a cut always sum to the
# total), and byte order is Python's str order. So a collection is split
# into byte-wise ranges, and each range is streamed by ONE cursor with no
# SORT at all. A cursor cannot be resumed, so a range that fails part-way is
# discarded and pulled again; ranges are kept small (~25 M rows, well under
# an hour) so that costs little. Shards land in a private directory and are
# moved into place only once the row count matches the plan, so `register`
# never sees a partial range.

_CUT_CHARS = [chr(c) for c in range(0x21, 0x7F)]


def _range_filter(lo: str | None, hi: str | None) -> tuple[str, dict[str, Any]]:
    parts, bind = [], {}
    if lo is not None:
        parts.append("d._key >= @lo")
        bind["lo"] = lo
    if hi is not None:
        parts.append("d._key < @hi")
        bind["hi"] = hi
    return (f"FILTER {' AND '.join(parts)} " if parts else ""), bind


def _range_count(client: ArangoClient, collection: str, lo: str | None,
                 hi: str | None) -> int | None:
    """Rows in [lo, hi) by the primary index, or None if counting timed out."""
    where, bind = _range_filter(lo, hi)
    try:
        res = client.open_cursor(
            f"FOR d IN {collection} {where}COLLECT WITH COUNT INTO n RETURN n",
            batch_size=1, bind_vars=bind or None)
        return int(res["result"][0])
    except urllib.error.HTTPError as exc:
        if exc.code in (502, 503, 504):
            return None
        raise


def _range_edges(client: ArangoClient, collection: str, lo: str | None,
                 hi: str | None) -> tuple[str | None, str | None]:
    """Two keys from the range (byte-wise min/max of the ends of the sort).

    Only a hint for where to cut: the sort order is the collation's, not the
    filter's, so these are not guaranteed to be the byte-wise extremes.
    """
    where, bind = _range_filter(lo, hi)
    keys: list[str] = []
    for direction in ("", "DESC "):
        res = client.open_cursor(
            f"FOR d IN {collection} {where}SORT d._key {direction}LIMIT 1 RETURN d._key",
            batch_size=1, bind_vars=bind or None)
        keys += res["result"]
    return (min(keys), max(keys)) if keys else (None, None)


def _plan_ranges(client: ArangoClient, collection: str, target_rows: int, *,
                 max_depth: int = 16) -> list[dict[str, Any]]:
    """Contiguous byte-wise ranges of at most ~target_rows, measured by counts.

    For a range that is too big: cut at P + each printable character, where
    P is the common prefix of two keys inside it, count the pieces, recurse.
    Cuts are ordered in Python (byte order, which is the filter's order), so
    consecutive ranges share a boundary and every key lands in exactly one.
    """
    out: list[dict[str, Any]] = []

    def split(lo: str | None, hi: str | None, n: int | None, depth: int) -> None:
        if n is not None and (n <= target_rows or depth >= max_depth):
            out.append({"lo": lo, "hi": hi, "rows": n})
            return
        a, b = _range_edges(client, collection, lo, hi)
        if a is None:
            out.append({"lo": lo, "hi": hi, "rows": 0})
            return
        cuts: list[str] = []
        if a != b:
            i = 0
            while i < min(len(a), len(b)) and a[i] == b[i]:
                i += 1
            for prefix in (a[:i], a[:i + 1]):
                cuts = sorted({prefix + ch for ch in _CUT_CHARS
                               if (lo is None or prefix + ch > lo)
                               and (hi is None or prefix + ch < hi)})
                if cuts:
                    break
        if not cuts:
            # One key, or no cut falls inside: take it whole.
            out.append({"lo": lo, "hi": hi, "rows": n})
            return
        bounds = [lo] + cuts + [hi]
        logging.info("plan %s: depth %d, [%r, %r) has %s rows -> %d pieces",
                     collection, depth, lo, hi,
                     f"{n:,}" if n is not None else "?", len(bounds) - 1)
        for x, y in zip(bounds, bounds[1:]):
            t0 = time.time()
            m = _range_count(client, collection, x, y)
            if m or time.time() - t0 > 5:
                logging.info("  [%r, %r) = %s  (%.1fs)", x, y,
                             f"{m:,}" if m is not None else "count timed out",
                             time.time() - t0)
            split(x, y, m, depth + 1)

    split(None, None, client.collection_count(collection), 0)
    # Fold empty pieces into a neighbour, keeping the plan contiguous: range
    # i ends exactly where range i+1 begins, from the first key to the last.
    merged: list[dict[str, Any]] = []
    for r in out:
        if merged and (r["rows"] == 0 or merged[-1]["rows"] == 0):
            a, b = merged[-1]["rows"], r["rows"]
            merged[-1]["hi"] = r["hi"]
            merged[-1]["rows"] = None if a is None or b is None else a + b
        else:
            merged.append(dict(r))
    for i, r in enumerate(merged):
        r["tag"] = f"r{i:04d}"
    return merged


def _plan_path(collection: str) -> Path:
    return _state_path(f"{collection}__ranges")


def _stream_range(client: ArangoClient, collection: str, r: dict[str, Any], *,
                  batch_size: int, attempts: int = 4) -> dict[str, Any]:
    """Pull one planned range with a single streaming cursor; all or nothing."""
    tag = r["tag"]
    state_file = _state_path(f"{collection}__{tag}")
    if state_file.exists():
        state = json.loads(state_file.read_text())
        if state.get("status") == "done":
            return state
    final_dir = _collection_dir(collection)
    work = final_dir / f"_partial_{tag}"
    where, bind = _range_filter(r["lo"], r["hi"])
    batch = _pick_initial_batch_size(client, collection, batch_size)
    state: dict[str, Any] = {"collection": collection, "tag": tag, "lo": r["lo"],
                             "hi": r["hi"], "planned_rows": r["rows"]}
    for attempt in range(1, attempts + 1):
        # A cursor cannot be resumed: every attempt starts the range afresh.
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)
        started, n, idx, cursor_id = time.time(), 0, 0, None
        try:
            res = client._request("/_api/cursor", method="POST", timeout=600, body={
                "query": f"FOR d IN {collection} {where}RETURN d",
                "bindVars": bind, "batchSize": batch, "ttl": 1800,
                "options": {"stream": True}})
            while True:
                cursor_id = res.get("id")
                rows = res.get("result", [])
                if rows:
                    _write_rows(work / f"{idx:05d}.parquet", [_flatten_row(x) for x in rows])
                    n += len(rows)
                    idx += 1
                if not res.get("hasMore"):
                    break
                res = client._request(f"/_api/cursor/{cursor_id}", method="PUT", timeout=600)
        except Exception as exc:
            if cursor_id:
                try:
                    client._request(f"/_api/cursor/{cursor_id}", method="DELETE", timeout=60)
                except Exception:
                    pass
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (502, 503, 504) and batch > 1:
                batch = max(1, batch // 2)
            logging.warning("%s[%s] attempt %d failed after %d rows: %s (batch now %d)",
                            collection, tag, attempt, n, exc, batch)
            time.sleep(min(300.0, 30.0 * attempt))
            continue
        if r["rows"] is not None and n != r["rows"]:
            state.update(status="count_mismatch", rows=n)
            logging.error("%s[%s]: streamed %d rows, plan counted %d; keeping it out.",
                          collection, tag, n, r["rows"])
            _state_path(f"{collection}__{tag}").write_text(json.dumps(state, indent=2))
            return state
        for shard in sorted(work.glob("*.parquet")):
            shard.rename(final_dir / f"{tag}_{shard.name}")
        work.rmdir()
        state.update(status="done", rows=n, shards=idx,
                     seconds=round(time.time() - started, 1))
        state_file.write_text(json.dumps(state, indent=2))
        logging.info("%s[%s] done: %d rows in %d shards, %.0f rows/s", collection, tag,
                     n, idx, n / max(time.time() - started, 1e-3))
        return state
    state.update(status="failed")
    state_file.write_text(json.dumps(state, indent=2))
    return state


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    pa = _require_pkg("pyarrow", "Required to write Parquet shards.")
    pq = __import__("pyarrow.parquet", fromlist=["parquet"])
    keys = {k for row in rows for k in row.keys()}
    cols = {k: _coerce_column([row.get(k) for row in rows]) for k in keys}
    pq.write_table(pa.Table.from_pydict(cols), path, compression="zstd")


def cmd_plan_ranges(args: argparse.Namespace) -> int:
    setup_logging()
    client = make_client()
    p = _plan_path(args.collection)
    if p.exists() and not args.replan:
        plan = json.loads(p.read_text())
        print(f"Plan exists: {p} ({len(plan['ranges'])} ranges); --replan to redo it.")
        return 0
    t0 = time.time()
    ranges = _plan_ranges(client, args.collection, args.target_rows)
    total = client.collection_count(args.collection)
    p.write_text(json.dumps({"collection": args.collection, "target_rows": args.target_rows,
                             "planned_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
                             "total_rows": total, "ranges": ranges}, indent=2))
    counted = sum(r["rows"] or 0 for r in ranges)
    print(f"{args.collection}: {len(ranges)} ranges; {counted:,} of {total:,} rows counted; "
          f"largest {max((r['rows'] or 0) for r in ranges):,}; "
          f"{sum(r['rows'] is None for r in ranges)} uncounted  ({time.time() - t0:.0f}s)")
    print(f"Wrote {p}")
    return 0


def cmd_pull_planned(args: argparse.Namespace) -> int:
    setup_logging()
    p = _plan_path(args.collection)
    if not p.exists():
        raise SystemExit(f"No plan at {p}; run `kg-mirror plan-ranges --collection "
                         f"{args.collection}` first.")
    ranges = json.loads(p.read_text())["ranges"]
    if not 0 <= args.index < len(ranges):
        raise SystemExit(f"--index {args.index} out of range (plan has {len(ranges)}).")
    state = _stream_range(make_client(), args.collection, ranges[args.index],
                          batch_size=args.batch_size)
    print(f"\n{args.collection}[{state['tag']}]: status={state['status']}, "
          f"rows={state.get('rows', 0):,} (planned {state['planned_rows'] or 0:,})")
    return 0 if state["status"] == "done" else 1


# ---------------------------------------------------------------------------
# Pull-all orchestration
# ---------------------------------------------------------------------------

def _default_order(client: ArangoClient, skip: set[str]) -> list[str]:
    """Small collections first so the user sees progress quickly."""
    cols = client.list_collections()
    rows = []
    for c in cols:
        if c["name"] in skip:
            continue
        try:
            n = client.collection_count(c["name"])
        except Exception:
            n = float("inf")
        rows.append((n, c["name"]))
    rows.sort()
    return [name for _, name in rows]


def cmd_pull_all(args: argparse.Namespace) -> int:
    setup_logging()
    client = make_client()
    skip = set(args.skip.split(",")) if args.skip else set(DEFAULT_SKIP)
    if args.include_giants:
        skip = set()
    order = _default_order(client, skip)
    print(f"Mirroring {len(order)} collections "
          f"(skipping: {sorted(skip) or 'none'}).")
    summary = []
    for collection in order:
        if args.only and collection not in set(args.only.split(",")):
            continue
        if args.max_collection_bytes is not None:
            try:
                fig = client.collection_figures(collection)
                size = int(fig.get("documentsSize", 0))
                if size > args.max_collection_bytes:
                    logging.info("Skipping %s (%.1f GB > %s GB cap)",
                                  collection, size / 1e9,
                                  args.max_collection_bytes / 1e9)
                    continue
            except Exception:
                pass
        try:
            state = _pull_collection(client, collection,
                                      batch_size=args.batch_size,
                                      max_rows=args.max_rows,
                                      restart=args.restart)
            summary.append({"collection": collection, **state})
        except Exception as exc:
            logging.exception("Failed pulling %s: %s", collection, exc)
            summary.append({"collection": collection, "status": "error",
                              "error": str(exc)})
    # Write summary
    KG_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = KG_DIR / f"_summary_{time.strftime('%Y%m%d_%H%M%S')}.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nWrote summary: {summary_path}")
    return 0


# ---------------------------------------------------------------------------
# Register Parquet shards as DuckDB views
# ---------------------------------------------------------------------------

def cmd_register(_args: argparse.Namespace) -> int:
    duckdb = _require_pkg("duckdb", "Required to register the warehouse.")
    setup_logging()
    WAREHOUSE_DIR.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(DUCKDB_PATH))
    n_views = 0
    for sub in sorted(KG_DIR.iterdir()):
        if not sub.is_dir() or sub.name.startswith("_"):
            continue
        shards = sorted(sub.glob("*.parquet"))
        if not shards:
            continue
        view = f"kg_{sub.name}"
        glob = str(sub.resolve() / "*.parquet")
        con.execute(f"CREATE OR REPLACE VIEW {view} AS SELECT * FROM read_parquet('{glob}')")
        cnt = con.execute(f"SELECT count(*) FROM {view}").fetchone()[0]
        logging.info("View %s: %d rows over %d shards", view, cnt, len(shards))
        n_views += 1
    con.close()
    print(f"Registered {n_views} views in {DUCKDB_PATH}")
    return 0


def cmd_verify(_args: argparse.Namespace) -> int:
    duckdb = _require_pkg("duckdb", "Required to verify the warehouse.")
    setup_logging()
    if not DUCKDB_PATH.exists():
        raise SystemExit(f"No warehouse at {DUCKDB_PATH}; run `register` first.")
    con = duckdb.connect(str(DUCKDB_PATH), read_only=True)
    rows = con.execute(
        "SELECT view_name FROM duckdb_views() WHERE view_name LIKE 'kg_%' ORDER BY view_name"
    ).fetchall()
    print(f"{'view':45}  {'rows':>14}")
    print("-" * 65)
    for (view,) in rows:
        n = con.execute(f"SELECT count(*) FROM {view}").fetchone()[0]
        print(f"{view:45}  {n:>14,}")
    con.close()
    return 0


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------

def write_playbook() -> Path:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    path = SKILL_DOC_DIR / "KG_MIRROR_SKILL.md"
    path.write_text(
        """# Skill: Local IGVF Knowledge Graph Mirror

Streams Arango collections via the read-only AQL cursor API and persists
them as Parquet shards on disk, then registers a DuckDB warehouse with
one view per collection. Lets `igvfagent kg ...` and downstream skills
run **offline** against the cached copy.

## Defaults

- Skip list: `variants` (~944 GB doc) and `variants_variants` (~531 GB edge).
  These two account for ~75% of the upstream KG. Override with
  `--skip ""` or `--include-giants`.
- Storage:
  - `Data/Warehouse/KG/<collection>/<NNNN>.parquet` — ZSTD-compressed shards
  - `Data/Warehouse/KG/_state/<collection>.json` — resume cursor
  - `Data/Warehouse/igvf_kg_mirror.duckdb` — DuckDB with `kg_<collection>` views

## Commands

```bash
# 1. Inventory: list all collections with doc count + bytes
igvfagent kg-mirror inventory

# 2. Mirror everything except the two giants
igvfagent kg-mirror pull-all

# 3. Mirror one collection at a time
igvfagent kg-mirror pull --collection genes
igvfagent kg-mirror pull --collection coding_variants --batch-size 10000

# 4. Resume a long-running pull (state is on disk)
igvfagent kg-mirror pull --collection variants_proteins   # resumes automatically

# 5. Cap to small/medium collections only (under 10 GB each)
igvfagent kg-mirror pull-all --max-collection-bytes 10000000000

# 6. Register Parquet shards as DuckDB views
igvfagent kg-mirror register
igvfagent kg-mirror verify
```

## Querying the local mirror

```python
import duckdb
con = duckdb.connect("Data/Warehouse/igvf_kg_mirror.duckdb", read_only=True)
con.sql("SELECT count(*) FROM kg_genes").show()
con.sql("SELECT * FROM kg_proteins LIMIT 5").show()
# Edge collections: same naming convention — kg_genes_pathways, kg_variants_genes, ...
```

## Schema notes

- Every row is a flat dict. Scalar fields (str/int/float/bool) stay as
  scalars. Anything nested (list, dict) is `json.dumps`'d into a string
  column so the Parquet schema stays consistent across batches.
- `_id`, `_key`, `_from`, `_to` are preserved verbatim. Edge collections
  always carry `_from` and `_to` references like `genes/<gene_id>`.

## License

Apache-2.0. Uses only stdlib + pyarrow (Apache-2) + duckdb (MIT). The
Arango HTTP API is hit via basic-auth with the read-only `guest` account.
""",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="IGVF KG local mirror (Arango -> DuckDB/Parquet).")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("inventory", help="List collections with doc count + bytes.")

    p = sub.add_parser("pull", help="Mirror a single collection (resumable).")
    p.add_argument("--collection", required=True)
    p.add_argument("--batch-size", type=int, default=5000)
    p.add_argument("--max-rows", type=int, default=None,
                    help="Optional cap (for smoke testing).")
    p.add_argument("--restart", action="store_true",
                    help="Ignore on-disk state and start over.")

    p = sub.add_parser("plan-ranges",
                        help="Split a collection into byte-wise key ranges of about "
                             "--target-rows each (measured by counts), for parallel pulls.")
    p.add_argument("--collection", required=True)
    p.add_argument("--target-rows", type=int, default=25_000_000)
    p.add_argument("--replan", action="store_true",
                    help="Overwrite an existing plan (only before its ranges are pulled).")

    p = sub.add_parser("pull-planned",
                        help="Stream one range of a plan-ranges plan (all or nothing; "
                             "a finished range is skipped on re-run).")
    p.add_argument("--collection", required=True)
    p.add_argument("--index", type=int, required=True)
    p.add_argument("--batch-size", type=int, default=50000)

    p = sub.add_parser("pull-all", help="Mirror everything except the skip list.")
    p.add_argument("--skip", default=",".join(DEFAULT_SKIP),
                    help=f"Comma-separated collections to skip (default: {','.join(DEFAULT_SKIP)}).")
    p.add_argument("--only", default=None,
                    help="Comma-separated collections to include (overrides default order).")
    p.add_argument("--include-giants", action="store_true",
                    help="Override skip list and mirror variants/variants_variants too.")
    p.add_argument("--max-collection-bytes", type=int, default=None,
                    help="Per-collection upper bound on Arango documentsSize (bytes).")
    p.add_argument("--batch-size", type=int, default=5000)
    p.add_argument("--max-rows", type=int, default=None,
                    help="Optional cap per collection (for smoke testing).")
    p.add_argument("--restart", action="store_true")

    sub.add_parser("register", help="Register parquet shards as DuckDB views.")
    sub.add_parser("verify", help="Print row counts per view.")
    sub.add_parser("write-playbook", help="Write the skill's markdown playbook.")

    args = parser.parse_args(argv)
    if args.command == "inventory":
        return cmd_inventory(args)
    if args.command == "pull":
        return cmd_pull(args)
    if args.command == "plan-ranges":
        return cmd_plan_ranges(args)
    if args.command == "pull-planned":
        return cmd_pull_planned(args)
    if args.command == "pull-all":
        return cmd_pull_all(args)
    if args.command == "register":
        return cmd_register(args)
    if args.command == "verify":
        return cmd_verify(args)
    if args.command == "write-playbook":
        path = write_playbook()
        print(f"Wrote {path}")
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
