"""kg-mirror key-range pulls get every row exactly once from the IGVF catalog.

`pull` pages with `SORT d._key LIMIT @skip, @count`: correct, but each page
costs more than the last, so the billion-row collections never finish. The
fix that looks obvious, keyset paging (`FILTER d._key > @last SORT d._key`),
loses rows there: the catalog is a sharded cluster whose key FILTER compares
bytes (each shard's primary index) while the SORT merging the shards uses a
collation, so a key sorted after a page's last row but byte-wise below it is
excluded from every later page.

`plan-ranges` + `pull-planned` rely only on what is consistent: a key-range
FILTER on its own. This checks, against a stand-in server that behaves the
same way (byte-order filter, collated sort, rows streamed in no order):
  - the plan is contiguous and every range is within the target size, even
    for keys shaped like `variants` (SPDI ids, random ids at the ends);
  - streaming every planned range yields every row exactly once;
  - a stream that dies part-way leaves nothing behind and the retry is exact;
  - a finished range is not pulled again, and a count mismatch is kept out;
  - the plan still covers every key where the optimizer's collation and the
    scan's byte order disagree about a range's bounds;
  - the stand-in really does lose rows under keyset paging, so the test is
    modelling the failure that matters.

    python3 Scripts/test_kg_mirror.py
"""
from __future__ import annotations

import io
import os
import random
import sys
import tempfile
import urllib.error
from pathlib import Path

TMP = tempfile.mkdtemp(prefix="kgmirror_test_")
os.environ["IGVF_PROJECT_ROOT"] = TMP
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kg_mirror_skill as M  # noqa: E402

M.time.sleep = lambda _s: None  # no real backoff in a test

FAILED: "list[str]" = []


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}"
          + ("" if ok else f"   got {got!r}, want {want!r}"))
    if not ok:
        FAILED.append(name)


def collated(k: str):
    # The merged SORT's order: not byte order ("Ic" before "IR").
    return (k.lower(), k)


class FakeCatalog:
    """Byte-order key filters, a collated SORT, streamed rows in no order."""

    def __init__(self, keys, *, die_after_batches: "int | None" = None):
        self.docs = [{"_key": k, "v": i} for i, k in enumerate(keys)]
        self.die_after = die_after_batches
        self.cursors: "dict[str, list]" = {}
        self.deleted: "list[str]" = []
        self.rng = random.Random(7)

    def _in(self, lo, hi):
        # The catalog's optimizer answers "nothing" when the collation puts lo
        # at or after hi, before scanning; the scan itself compares bytes.
        if lo is not None and hi is not None and collated(lo) >= collated(hi):
            return []
        return [d for d in self.docs
                if (lo is None or d["_key"] >= lo) and (hi is None or d["_key"] < hi)]

    def collection_figures(self, _name):
        return {"documentsSize": 100 * len(self.docs)}

    def collection_count(self, _name):
        return len(self.docs)

    def open_cursor(self, query, *, batch_size=5000, bind_vars=None):
        b = bind_vars or {}
        if query.startswith("FOR p IN @pairs"):
            return {"result": [collated(x) < collated(y) for x, y in b["pairs"]]}
        rows = self._in(b.get("lo"), b.get("hi"))
        if "COLLECT WITH COUNT" in query:
            return {"result": [len(rows)]}
        if query.endswith("LIMIT 1 RETURN d._key"):
            rows = sorted(rows, key=lambda d: collated(d["_key"]),
                          reverse=" DESC " in query)
            return {"result": [d["_key"] for d in rows[:1]]}
        if "FILTER d._key > @last" in query:        # keyset paging, for the control
            rows = [d for d in self.docs if d["_key"] > b["last"]] if b.get("last") else self.docs
            rows = sorted(rows, key=lambda d: collated(d["_key"]))
            return {"result": rows[:b["count"]]}
        raise AssertionError(f"unexpected query {query}")

    def _request(self, path, *, method="GET", body=None, timeout=120):
        if method == "POST" and path == "/_api/cursor":
            assert "SORT" not in body["query"], "range streams must not sort"
            b = body["bindVars"]
            rows = self._in(b.get("lo"), b.get("hi"))
            self.rng.shuffle(rows)
            cid = f"c{len(self.cursors)}"
            n = body["batchSize"]
            self.cursors[cid] = [rows[i:i + n] for i in range(0, len(rows), n)] or [[]]
            self.cursors[cid + ":served"] = 0
            return self._next(cid)
        if method == "PUT":
            return self._next(path.rsplit("/", 1)[1])
        if method == "DELETE":
            self.deleted.append(path.rsplit("/", 1)[1])
            return {}
        raise AssertionError(f"unexpected {method} {path}")

    def _next(self, cid):
        served = self.cursors[cid + ":served"]
        if self.die_after is not None and served >= self.die_after:
            self.die_after = None                    # die once, then behave
            raise urllib.error.HTTPError("x", 503, "Unavailable", {}, io.BytesIO())
        batches = self.cursors[cid]
        self.cursors[cid + ":served"] = served + 1
        return {"id": cid, "result": batches[served], "hasMore": served + 1 < len(batches)}


def mirrored(collection: str) -> "list[str]":
    pq = __import__("pyarrow.parquet", fromlist=["parquet"])
    keys: "list[str]" = []
    for shard in sorted((M.KG_DIR / collection).glob("*.parquet")):
        keys += pq.read_table(shard, columns=["_key"]).column("_key").to_pylist()
    return keys


def variants_like() -> "list[str]":
    # Mostly SPDI ids, with random-looking ids at both ends of the order:
    # the shape that defeats splitting on the first character.
    spdi = [f"NC_0000{c:02d}.11:{p}:A:G" for c in range(1, 23) for p in range(0, 400, 7)]
    edge = ["--IcYmh", "--IR_J5", "-0lxoka", "-0QGRJ2", "_x9", "zzu9yM", "zzNy6A", "Aq1", "b7Z"]
    return spdi + edge


def main() -> int:
    keys = variants_like()
    want = sorted(keys)

    print("\ncontrol: the stand-in loses rows under keyset paging, as the catalog does")
    fake, got, last = FakeCatalog(keys), [], None
    while True:
        page = fake.open_cursor("FOR d IN c FILTER d._key > @last SORT d._key LIMIT @count RETURN d",
                                bind_vars={"last": last, "count": 5})["result"]
        if not page:
            break
        got += [d["_key"] for d in page]
        last = page[-1]["_key"]
    check("keyset paging misses rows here", len(set(got)) < len(keys), True)

    print("\nplan")
    plan = M._plan_ranges(FakeCatalog(keys), "c1", target_rows=100)
    check("counts every row", sum(r["rows"] for r in plan), len(keys))
    check("no range above the target", max(r["rows"] for r in plan) <= 100, True)
    check("the SPDI bulk is split, not left as one range", len(plan) >= len(keys) // 100, True)
    check("contiguous: open at both ends, each range starts where the last ended",
          (plan[0]["lo"], plan[-1]["hi"], all(a["hi"] == b["lo"] for a, b in zip(plan, plan[1:]))),
          (None, None, True))

    print("\nplan where the collation and byte order disagree (the DZANK1 case)")
    # Byte-wise "DZANK1" lies in [DZ, D[), but the collation sorts "[" and
    # "a" before "Z", so the server answered 0 for that piece: 385,779
    # coding_variants keys fell in no range of the first plans.
    genes = [f"{g}_{i:03d}" for g in ("DYDC1", "DZANK1", "DZIP3", "Dab1", "E2F1", "EGFR") for i in range(40)]
    plan2 = M._plan_ranges(FakeCatalog(genes), "c9", target_rows=30)
    check("every key is counted by some range", sum(r["rows"] for r in plan2), len(genes))
    fake = FakeCatalog(genes)
    for r in plan2:
        M._stream_range(fake, "c9", r, batch_size=17)
    got = mirrored("c9")
    check("pulling the plan gets every key once", (sorted(got), len(got)), (sorted(genes), len(genes)))

    print("\nstream every planned range")
    fake = FakeCatalog(keys)
    states = [M._stream_range(fake, "c1", r, batch_size=17) for r in plan]
    check("every range done", {s["status"] for s in states}, {"done"})
    got = mirrored("c1")
    check("every row, exactly once", (sorted(got), len(got)), (want, len(want)))
    check("no partial directories left",
          [p.name for p in (M.KG_DIR / "c1").iterdir() if p.is_dir()], [])
    before = sorted(p.name for p in (M.KG_DIR / "c1").glob("*.parquet"))
    M._stream_range(FakeCatalog([]), "c1", plan[0], batch_size=17)
    check("a finished range is not pulled again",
          sorted(p.name for p in (M.KG_DIR / "c1").glob("*.parquet")), before)

    print("\na stream that dies part-way")
    fake = FakeCatalog(keys, die_after_batches=3)
    r = {"lo": None, "hi": None, "rows": len(keys), "tag": "r0000"}
    st = M._stream_range(fake, "c2", r, batch_size=50)
    check("retried to completion", st["status"], "done")
    got = mirrored("c2")
    check("no row from the dead attempt survives as a duplicate", (sorted(got), len(got)),
          (want, len(want)))
    check("the dead cursor was released", fake.deleted, ["c0"])

    print("\na job killed mid-range (SLURM time limit), then re-run")
    work = M._collection_dir("c4") / "_partial_r0000"
    work.mkdir(parents=True)
    M._write_rows(work / "00099.parquet", [{"_key": want[0], "v": -1}])
    st = M._stream_range(FakeCatalog(keys), "c4",
                         {"lo": None, "hi": None, "rows": len(keys), "tag": "r0000"},
                         batch_size=50)
    got = mirrored("c4")
    check("the killed job's leftover shard is not published", (sorted(got), len(got)),
          (want, len(want)))

    print("\na range whose rows do not match the plan")
    st = M._stream_range(FakeCatalog(keys), "c3",
                         {"lo": None, "hi": None, "rows": len(keys) + 1, "tag": "r0000"},
                         batch_size=50)
    check("marked, not published", (st["status"], mirrored("c3")), ("count_mismatch", []))

    print(f"\n{'FAILED: ' + ', '.join(FAILED) if FAILED else 'all passed'}")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
