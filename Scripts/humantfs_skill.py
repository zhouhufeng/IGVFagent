#!/usr/bin/env python3
"""Human Transcription Factors (HTF) database — local mirror + query.

Builds and queries a local SQLite mirror of the Lambert/Jolma/Hughes
**Human Transcription Factors** database v1.01
(https://humantfs.ccbr.utoronto.ca, Lambert et al., *Cell* 2018,
doi:10.1016/j.cell.2018.01.029).

The database is the canonical answer to "is this gene a transcription
factor, and what does it bind?": 2,765 assessed proteins of which
**1,639 are curated TFs**, each with its DNA-binding domain family,
binding mode, motif status and cross-references (Ensembl, HGNC,
EntrezGene, InterPro, PDB). It is the TF list the Tabula Sapiens 2.0
paper uses, and it underpins any cell-type-specificity or regulon
analysis this agent runs.

Subcommands
-----------
  build-db     Download v1.01 and build ``Data/HumanTFs/human_tfs.sqlite``.
  info         Row counts, DBD families, motif coverage, source versions.
  list         List TFs, filterable by DBD family / binding mode / motif.
  lookup       Everything the database knows about one or more genes.
  is-tf        Partition a pasted gene list into TFs and non-TFs.
  motifs       Motif records for a TF (source, ID, evidence).
  families     DBD family census.
  export       Write a TF list (symbols or Ensembl IDs) for another tool.
  write-playbook

Why a local mirror rather than live calls: the upstream site serves
static files with no query API, the whole database is ~2 MB, and every
downstream analysis needs to ask "is X a TF?" thousands of times. One
download, then offline forever.

License: Apache-2.0. The HTF database is redistributed by download at
run time, never vendored into this repository; its own terms and the
Lambert 2018 citation apply to the data.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import os
import sqlite3
import sys
import time
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Optional

ROOT = Path(
    os.environ.get("IGVF_PROJECT_ROOT")
    or Path(__file__).resolve().parents[1]
).resolve()
DATA_DIR = ROOT / "Data" / "HumanTFs"
DB_PATH = DATA_DIR / "human_tfs.sqlite"
DOCS_DIR = ROOT / "Docs"
LOG_DIR = DOCS_DIR / "Logs"
REPORT_DIR = DOCS_DIR / "HumanTFs"
SKILL_DOC_DIR = DOCS_DIR / "Skills"

HTF_VERSION = "v_1.01"
HTF_BASE = f"https://humantfs.ccbr.utoronto.ca/download/{HTF_VERSION}"
HTF_FILES = {
    "database": f"{HTF_BASE}/DatabaseExtract_{HTF_VERSION}.csv",
    "tf_names": f"{HTF_BASE}/TF_names_{HTF_VERSION}.txt",
    "tf_ensembl": f"{HTF_BASE}/TFs_Ensembl_{HTF_VERSION}.txt",
    "motifs": f"{HTF_BASE}/Human_TF_MotifList_{HTF_VERSION}.csv",
}
PWM_URL = f"{HTF_BASE}/PWMs.zip"

# The paper's own count, and a useful tripwire: if a rebuild yields a
# different number the upstream release changed under us.
EXPECTED_TF_COUNT = 1639


def setup_logging() -> Path:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    p = LOG_DIR / f"humantfs_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(p), logging.StreamHandler(sys.stdout)],
    )
    return p


def _get(url: str, *, timeout: int = 180) -> bytes:
    req = urllib.request.Request(
        url, headers={"User-Agent": "IGVFagent/humantfs"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def connect(*, read_only: bool = True) -> sqlite3.Connection:
    if not DB_PATH.exists():
        raise SystemExit(
            f"No HTF database at {DB_PATH}.\n"
            f"Build it first:  igvfagent humantfs build-db")
    con = sqlite3.connect(f"file:{DB_PATH}?mode=ro" if read_only
                          else str(DB_PATH), uri=read_only)
    con.row_factory = sqlite3.Row
    return con


# ---------------------------------------------------------------------------
# build-db
# ---------------------------------------------------------------------------

SCHEMA = """
CREATE TABLE IF NOT EXISTS tf (
    ensembl_id      TEXT PRIMARY KEY,
    symbol          TEXT,
    dbd             TEXT,
    is_tf           INTEGER NOT NULL,
    tf_assessment   TEXT,
    binding_mode    TEXT,
    motif_status    TEXT,
    entrez_id       TEXT,
    entrez_desc     TEXT,
    interpro        TEXT,
    pdb             TEXT,
    notes           TEXT
);
CREATE INDEX IF NOT EXISTS idx_tf_symbol ON tf(symbol);
CREATE INDEX IF NOT EXISTS idx_tf_is_tf  ON tf(is_tf);
CREATE INDEX IF NOT EXISTS idx_tf_dbd    ON tf(dbd);

CREATE TABLE IF NOT EXISTS motif (
    ensembl_id      TEXT,
    symbol          TEXT,
    motif_id        TEXT,
    source          TEXT,
    evidence        TEXT,
    best_motif      INTEGER
);
CREATE INDEX IF NOT EXISTS idx_motif_symbol ON motif(symbol);
CREATE INDEX IF NOT EXISTS idx_motif_ens    ON motif(ensembl_id);

CREATE TABLE IF NOT EXISTS pwm (
    motif_id        TEXT PRIMARY KEY,
    n_positions     INTEGER,
    matrix_json     TEXT
);

CREATE TABLE IF NOT EXISTS meta (
    key             TEXT PRIMARY KEY,
    value           TEXT
);
"""


def _col(row: "dict[str, Any]", *names: str) -> str:
    """First non-empty value among candidate column names.

    The upstream CSV's headers have drifted between releases (``HGNC
    symbol`` vs ``HGNC_symbol``), so match on several spellings rather
    than pinning one.
    """
    for n in names:
        v = row.get(n)
        if v not in (None, "", "NA", "None", "nan"):
            return str(v).strip()
    return ""


def cmd_build_db(args: argparse.Namespace) -> int:
    setup_logging()
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if DB_PATH.exists() and not args.force:
        raise SystemExit(
            f"{DB_PATH} already exists. Pass --force to rebuild.")

    logging.info("downloading HTF %s", HTF_VERSION)
    raw: "dict[str, bytes]" = {}
    for key, url in HTF_FILES.items():
        logging.info("  %s", url)
        raw[key] = _get(url)
    (DATA_DIR / f"DatabaseExtract_{HTF_VERSION}.csv").write_bytes(raw["database"])

    if DB_PATH.exists():
        DB_PATH.unlink()
    con = sqlite3.connect(str(DB_PATH))
    con.executescript(SCHEMA)

    # ── main table ────────────────────────────────────────────────────
    rows = list(csv.DictReader(io.StringIO(raw["database"].decode("utf8", "replace"))))
    n_tf = 0
    for r in rows:
        is_tf = 1 if _col(r, "Is TF?", "Is TF").lower().startswith("y") else 0
        n_tf += is_tf
        con.execute(
            "INSERT OR REPLACE INTO tf VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (
                _col(r, "Ensembl ID", "Ensembl_ID"),
                _col(r, "HGNC symbol", "HGNC_symbol", "Symbol"),
                _col(r, "DBD"),
                is_tf,
                _col(r, "TF assessment", "TF_assessment"),
                _col(r, "Binding mode", "Binding_mode"),
                _col(r, "Motif status", "Motif_status"),
                _col(r, "EntrezGene ID", "EntrezGene_ID"),
                _col(r, "EntrezGene Description", "EntrezGene_Description"),
                _col(r, "Interpro ID(s)", "Interpro_ID"),
                _col(r, "PDB ID", "PDB_ID"),
                _col(r, "Final Notes", "Final_Notes"),
            ),
        )

    # ── motifs ────────────────────────────────────────────────────────
    n_motif = 0
    try:
        mrows = list(csv.DictReader(
            io.StringIO(raw["motifs"].decode("utf8", "replace"))))
        for r in mrows:
            con.execute(
                "INSERT INTO motif VALUES (?,?,?,?,?,?)",
                (
                    _col(r, "Ensembl ID", "Ensembl_ID"),
                    _col(r, "HGNC symbol", "HGNC_symbol", "Symbol"),
                    _col(r, "CIS-BP ID", "Motif ID", "CISBP_ID", "Motif_ID"),
                    _col(r, "MSource_Identifier", "Source", "MSource_Type"),
                    _col(r, "TF_Status", "Evidence", "Motif_Type"),
                    1 if _col(r, "Best Motif(s)? (Figure 2A)",
                              "Best Motif").lower().startswith("t") else 0,
                ),
            )
            n_motif += 1
    except Exception as exc:  # pragma: no cover - upstream shape drift
        logging.warning("motif table skipped: %s", exc)

    # ── optional PWMs ─────────────────────────────────────────────────
    n_pwm = 0
    if args.with_pwms:
        logging.info("downloading PWMs")
        try:
            z = zipfile.ZipFile(io.BytesIO(_get(PWM_URL, timeout=600)))
            for name in z.namelist():
                if not name.lower().endswith(".txt") or name.endswith("/"):
                    continue
                mid = Path(name).stem
                mat = _parse_pwm(z.read(name).decode("utf8", "replace"))
                if mat:
                    con.execute(
                        "INSERT OR REPLACE INTO pwm VALUES (?,?,?)",
                        (mid, len(mat), json.dumps(mat)))
                    n_pwm += 1
        except Exception as exc:
            logging.warning("PWM download skipped: %s", exc)

    for k, v in (("version", HTF_VERSION),
                 ("built", time.strftime("%Y-%m-%d %H:%M:%S")),
                 ("source", "https://humantfs.ccbr.utoronto.ca"),
                 ("citation", "Lambert et al., Cell 2018, "
                              "doi:10.1016/j.cell.2018.01.029"),
                 ("n_rows", str(len(rows))), ("n_tf", str(n_tf)),
                 ("n_motif", str(n_motif)), ("n_pwm", str(n_pwm))):
        con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", (k, v))
    con.commit()
    con.close()

    print(f"Database:  {DB_PATH}")
    print(f"Rows:      {len(rows):,} assessed proteins")
    print(f"TFs:       {n_tf:,}")
    print(f"Motifs:    {n_motif:,}")
    if args.with_pwms:
        print(f"PWMs:      {n_pwm:,}")
    if n_tf != EXPECTED_TF_COUNT:
        print(f"\nNOTE: expected {EXPECTED_TF_COUNT} TFs for {HTF_VERSION} but "
              f"got {n_tf}. The upstream release may have changed; anything "
              f"citing the 1,639 figure should be re-checked.")
    return 0


def _parse_pwm(text: str) -> "list[list[float]]":
    """Parse a CIS-BP style PWM: header, then Pos A C G T rows."""
    out: "list[list[float]]" = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) != 5:
            continue
        try:
            out.append([float(x) for x in parts[1:]])
        except ValueError:
            continue  # header row
    return out


# ---------------------------------------------------------------------------
# Queries
# ---------------------------------------------------------------------------

def cmd_info(_args: argparse.Namespace) -> int:
    con = connect()
    meta = {r["key"]: r["value"] for r in con.execute("SELECT * FROM meta")}
    print(f"HTF database:  {DB_PATH}")
    for k in ("version", "built", "source", "citation"):
        if k in meta:
            print(f"  {k:10} {meta[k]}")
    n_all = con.execute("SELECT count(*) FROM tf").fetchone()[0]
    n_tf = con.execute("SELECT count(*) FROM tf WHERE is_tf=1").fetchone()[0]
    n_mot = con.execute("SELECT count(*) FROM motif").fetchone()[0]
    n_pwm = con.execute("SELECT count(*) FROM pwm").fetchone()[0]
    print(f"\n  assessed proteins  {n_all:,}")
    print(f"  curated TFs        {n_tf:,}")
    print(f"  motif records      {n_mot:,}")
    print(f"  PWMs               {n_pwm:,}")
    fam = con.execute(
        "SELECT dbd, count(*) n FROM tf WHERE is_tf=1 AND dbd<>'' "
        "GROUP BY dbd ORDER BY n DESC LIMIT 10").fetchall()
    n_fam = con.execute(
        "SELECT count(DISTINCT dbd) FROM tf WHERE is_tf=1 AND dbd<>''"
    ).fetchone()[0]
    print(f"\n  top DBD families ({n_fam} total):")
    for r in fam:
        print(f"    {r['dbd']:34} {r['n']:>5}")
    ms = con.execute(
        "SELECT motif_status, count(*) n FROM tf WHERE is_tf=1 "
        "GROUP BY motif_status ORDER BY n DESC").fetchall()
    print("\n  motif status:")
    for r in ms:
        print(f"    {(r['motif_status'] or '(blank)'):34} {r['n']:>5}")
    con.close()
    return 0


def _rows_out(rows: "list[sqlite3.Row]", cols: "list[str]",
              out: Optional[str], label: str) -> Optional[Path]:
    if not out:
        return None
    p = Path(out)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])
    return p


def cmd_list(args: argparse.Namespace) -> int:
    con = connect()
    q = "SELECT ensembl_id, symbol, dbd, binding_mode, motif_status FROM tf WHERE is_tf=1"
    params: "list[Any]" = []
    if args.dbd:
        q += " AND dbd LIKE ?"
        params.append(f"%{args.dbd}%")
    if args.binding_mode:
        q += " AND binding_mode LIKE ?"
        params.append(f"%{args.binding_mode}%")
    if args.with_motif:
        q += " AND motif_status LIKE '%Known motif%'"
    q += " ORDER BY symbol"
    rows = con.execute(q, params).fetchall()
    cols = ["symbol", "ensembl_id", "dbd", "binding_mode", "motif_status"]
    p = _rows_out(rows, cols, args.out, "tf_list")
    for r in rows[: args.limit]:
        print(f"  {r['symbol']:14} {r['ensembl_id']:18} {r['dbd']}")
    print(f"\n{len(rows):,} TF(s)" +
          (f"; showing {min(args.limit, len(rows))}" if len(rows) > args.limit else ""))
    if p:
        print(f"Wrote {p}")
    con.close()
    return 0


def _split_genes(text: str) -> "list[str]":
    import re
    return [g.strip().upper() for g in re.split(r"[,\s;]+", text or "") if g.strip()]


def cmd_lookup(args: argparse.Namespace) -> int:
    con = connect()
    genes = _split_genes(args.genes)
    if args.input:
        genes += _split_genes(Path(args.input).read_text())
    if not genes:
        raise SystemExit("give --genes or --input")
    for g in genes:
        r = con.execute(
            "SELECT * FROM tf WHERE upper(symbol)=? OR upper(ensembl_id)=?",
            (g, g)).fetchone()
        if not r:
            print(f"{g}: not in the HTF database (not assessed)")
            continue
        tag = "TF" if r["is_tf"] else "NOT a TF"
        print(f"\n{r['symbol']} ({r['ensembl_id']}) — {tag}")
        for k in ("dbd", "tf_assessment", "binding_mode", "motif_status",
                  "entrez_desc", "interpro", "pdb"):
            if r[k]:
                print(f"  {k:14} {r[k]}")
        mots = con.execute(
            "SELECT motif_id, source, evidence FROM motif WHERE upper(symbol)=?",
            (g,)).fetchall()
        if mots:
            print(f"  motifs         {len(mots)} "
                  f"(e.g. {', '.join(m['motif_id'] for m in mots[:3] if m['motif_id'])})")
    con.close()
    return 0


def cmd_is_tf(args: argparse.Namespace) -> int:
    """Partition a gene list into TFs and non-TFs."""
    con = connect()
    genes = _split_genes(args.genes)
    if args.input:
        genes += _split_genes(Path(args.input).read_text())
    if not genes:
        raise SystemExit("give --genes or --input")
    seen: "set[str]" = set()
    genes = [g for g in genes if not (g in seen or seen.add(g))]

    tfs, non, unknown = [], [], []
    for g in genes:
        r = con.execute(
            "SELECT symbol, is_tf, dbd FROM tf WHERE upper(symbol)=? "
            "OR upper(ensembl_id)=?", (g, g)).fetchone()
        if r is None:
            unknown.append(g)
        elif r["is_tf"]:
            tfs.append((r["symbol"], r["dbd"]))
        else:
            non.append(r["symbol"])
    con.close()

    print(f"Input genes:  {len(genes)}")
    print(f"  TFs:        {len(tfs)}")
    print(f"  not TFs:    {len(non)}")
    print(f"  unassessed: {len(unknown)}")
    if tfs:
        print("\nTFs:")
        for s, d in tfs[: args.limit]:
            print(f"  {s:14} {d}")
    if args.out:
        p = Path(args.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gene", "is_tf", "dbd"])
            for s, d in tfs:
                w.writerow([s, "yes", d])
            for s in non:
                w.writerow([s, "no", ""])
            for s in unknown:
                w.writerow([s, "unassessed", ""])
        print(f"\nWrote {p}")
    return 0


def cmd_motifs(args: argparse.Namespace) -> int:
    con = connect()
    rows = con.execute(
        "SELECT * FROM motif WHERE upper(symbol)=? ORDER BY best_motif DESC",
        (args.gene.upper(),)).fetchall()
    if not rows:
        print(f"No motif records for {args.gene}.")
        con.close()
        return 0
    print(f"{args.gene}: {len(rows)} motif record(s)")
    for r in rows[: args.limit]:
        star = "*" if r["best_motif"] else " "
        print(f" {star} {r['motif_id']:22} {r['source']:22} {r['evidence']}")
    n_pwm = con.execute(
        "SELECT count(*) FROM pwm WHERE motif_id IN "
        "(SELECT motif_id FROM motif WHERE upper(symbol)=?)",
        (args.gene.upper(),)).fetchone()[0]
    if n_pwm:
        print(f"  {n_pwm} of these have a PWM stored locally")
    con.close()
    return 0


def cmd_families(args: argparse.Namespace) -> int:
    con = connect()
    rows = con.execute(
        "SELECT dbd, count(*) n FROM tf WHERE is_tf=1 AND dbd<>'' "
        "GROUP BY dbd ORDER BY n DESC").fetchall()
    total = sum(r["n"] for r in rows)
    print(f"{len(rows)} DNA-binding-domain families over {total:,} TFs\n")
    for r in rows[: args.limit]:
        bar = "#" * max(1, round(60 * r["n"] / rows[0]["n"]))
        print(f"  {r['dbd'][:30]:30} {r['n']:>5}  {bar}")
    _rows_out(rows, ["dbd", "n"], args.out, "families")
    if args.out:
        print(f"\nWrote {args.out}")
    con.close()
    return 0


def cmd_export(args: argparse.Namespace) -> int:
    con = connect()
    field = "ensembl_id" if args.format == "ensembl" else "symbol"
    rows = con.execute(
        f"SELECT {field} v FROM tf WHERE is_tf=1 AND {field}<>'' "
        f"ORDER BY {field}").fetchall()
    vals = [r["v"] for r in rows]
    out = Path(args.out) if args.out else (DATA_DIR / f"tf_{args.format}.txt")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(vals) + "\n")
    print(f"Wrote {len(vals):,} TF {args.format}s to {out}")
    con.close()
    return 0


# ---------------------------------------------------------------------------
# Library API (used by the Tabula Sapiens skill)
# ---------------------------------------------------------------------------

def tf_symbols() -> "set[str]":
    """Every curated TF symbol, uppercased. Empty set if no DB is built."""
    if not DB_PATH.exists():
        return set()
    con = connect()
    out = {r[0].upper() for r in con.execute(
        "SELECT symbol FROM tf WHERE is_tf=1 AND symbol<>''")}
    con.close()
    return out


def tf_table() -> "dict[str, dict[str, Any]]":
    """Symbol -> full TF record, for annotating analysis output."""
    if not DB_PATH.exists():
        return {}
    con = connect()
    out = {r["symbol"].upper(): dict(r) for r in con.execute(
        "SELECT * FROM tf WHERE is_tf=1 AND symbol<>''")}
    con.close()
    return out


# ---------------------------------------------------------------------------
# Playbook
# ---------------------------------------------------------------------------

PLAYBOOK = """# Human Transcription Factors (HTF) database skill

Local SQLite mirror of the Lambert/Jolma/Hughes **Human Transcription
Factors** database v1.01 (https://humantfs.ccbr.utoronto.ca; Lambert et
al., *Cell* 2018, doi:10.1016/j.cell.2018.01.029).

2,765 assessed proteins, of which **1,639 are curated TFs**, each with
its DNA-binding-domain family, binding mode, motif status and
cross-references. This is the TF list the Tabula Sapiens 2.0 paper uses,
and the reference any regulon or cell-type-specificity analysis needs.

## Build it once

```bash
igvfagent humantfs build-db              # ~2 MB, seconds
igvfagent humantfs build-db --with-pwms  # also fetch position weight matrices
igvfagent humantfs info
```

Lands at `Data/HumanTFs/human_tfs.sqlite`. The build refuses to
overwrite without `--force`, and warns if the TF count is not 1,639 --
the number every downstream claim is calibrated to.

## Query it

```bash
igvfagent humantfs is-tf --genes "FOXP3, EOMES, ACTB, GAPDH, IKZF1"
igvfagent humantfs lookup --genes SATB2
igvfagent humantfs list --dbd "C2H2 ZF" --limit 20
igvfagent humantfs list --with-motif --out known_motif_tfs.tsv
igvfagent humantfs motifs --gene GATA1
igvfagent humantfs families
igvfagent humantfs export --format symbol --out tf_symbols.txt
igvfagent humantfs export --format ensembl
```

`is-tf` is the workhorse: paste any gene list and it partitions into
TFs, non-TFs and *unassessed* — the third bucket matters, because "not
in the database" is not the same as "not a TF".

## As a library

```python
from igvfagent.humantfs_skill import tf_symbols, tf_table
tfs = tf_symbols()          # {'AATF', 'ABL1', ...}
rec = tf_table()['GATA1']   # dbd, binding_mode, motif_status, ...
```

## Notes

**Zinc fingers dominate.** C2H2 ZF is by far the largest DBD family, and
a large share of those TFs are computationally predicted with no
validated motif — which is exactly why the Tabula Sapiens paper flags
ubiquitous zinc-finger TFs as an understudied group. `list --dbd "C2H2
ZF"` plus `--with-motif` gives you the split.

**Motif status is not binding evidence.** `Known motif` means a motif
exists in CIS-BP for that TF, from direct or inferred evidence; it says
nothing about whether the TF is active in your cells. Pair it with
regulon activity (`igvfagent tabula tf-regulons`) before claiming
function.

**The data is downloaded, not vendored.** This repository stores no HTF
content; `build-db` fetches it at run time. Cite Lambert 2018 and honour
the upstream terms when redistributing derived tables.
"""


def cmd_write_playbook(_a: argparse.Namespace) -> int:
    SKILL_DOC_DIR.mkdir(parents=True, exist_ok=True)
    p = SKILL_DOC_DIR / "HUMANTFS_SKILL.md"
    p.write_text(PLAYBOOK, encoding="utf-8")
    print(f"Wrote {p}")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: "list[str] | None" = None) -> int:
    ap = argparse.ArgumentParser(
        prog="humantfs",
        description="Human Transcription Factors database (Lambert 2018) — "
                    "local mirror and query.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("build-db", help="Download and build the local DB.")
    p.add_argument("--force", action="store_true")
    p.add_argument("--with-pwms", action="store_true",
                   help="Also download and store position weight matrices.")
    p.set_defaults(func=cmd_build_db)

    p = sub.add_parser("info", help="Database summary.")
    p.set_defaults(func=cmd_info)

    p = sub.add_parser("list", help="List TFs with optional filters.")
    p.add_argument("--dbd", default=None, help="DBD family substring.")
    p.add_argument("--binding-mode", default=None)
    p.add_argument("--with-motif", action="store_true")
    p.add_argument("--limit", type=int, default=40)
    p.add_argument("--out", default=None)
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("lookup", help="Full record for one or more genes.")
    p.add_argument("--genes", default=None)
    p.add_argument("--input", default=None)
    p.set_defaults(func=cmd_lookup)

    p = sub.add_parser("is-tf", help="Partition a gene list into TFs/non-TFs.")
    p.add_argument("--genes", default=None)
    p.add_argument("--input", default=None)
    p.add_argument("--limit", type=int, default=50)
    p.add_argument("--out", default=None)
    p.set_defaults(func=cmd_is_tf)

    p = sub.add_parser("motifs", help="Motif records for a TF.")
    p.add_argument("--gene", required=True)
    p.add_argument("--limit", type=int, default=30)
    p.set_defaults(func=cmd_motifs)

    p = sub.add_parser("families", help="DBD family census.")
    p.add_argument("--limit", type=int, default=40)
    p.add_argument("--out", default=None)
    p.set_defaults(func=cmd_families)

    p = sub.add_parser("export", help="Write the TF list for another tool.")
    p.add_argument("--format", choices=["symbol", "ensembl"], default="symbol")
    p.add_argument("--out", default=None)
    p.set_defaults(func=cmd_export)

    p = sub.add_parser("write-playbook", help="Emit the markdown playbook.")
    p.set_defaults(func=cmd_write_playbook)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
