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

DEFAULT_SKIP = ("variants", "variants_variants")


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


def _write_shard(rows: list[dict[str, Any]], collection: str, batch_idx: int) -> Path:
    pa = _require_pkg("pyarrow", "Required to write Parquet shards.")
    pq = __import__("pyarrow.parquet", fromlist=["parquet"])
    out = _collection_dir(collection) / f"{batch_idx:05d}.parquet"
    # Union schema across all rows so missing fields become nulls
    keys = {k for r in rows for k in r.keys()}
    table = pa.Table.from_pylist([{k: r.get(k) for k in keys} for r in rows])
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
                      max_retries: int = 6) -> dict[str, Any]:
    state = _read_state(collection)
    if restart:
        state = {"collection": collection, "rows_written": 0, "batches": 0,
                  "last_completed_batch": -1, "status": "restart"}
    next_batch = state["last_completed_batch"] + 1
    skip = state["rows_written"]
    effective_batch = _pick_initial_batch_size(client, collection, batch_size)
    if effective_batch != batch_size:
        logging.info("Auto-sized batch for %s: %d -> %d (avg doc size constraint)",
                      collection, batch_size, effective_batch)
    logging.info("Pulling %s (resuming at row %d, batch %d, batch_size=%d)",
                  collection, skip, next_batch, effective_batch)

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
