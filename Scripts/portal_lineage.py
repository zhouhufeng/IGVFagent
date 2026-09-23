#!/usr/bin/env python3
"""Walk the IGVF Portal graph from any accession to what has already been computed.

The IGVF data model links raw data to its processed results through a handful
of properties (see the "example mapping" diagram of the Portal's data model):

    MeasurementSet --files--> SequenceFile (reads) --seqspecs--> ConfigurationFile (seqspec)
    MeasurementSet --input_for--> AnalysisSet (intermediate analysis)
    AnalysisSet --files--> MatrixFile / TabularFile / AlignmentFile (derived_from reads)
    AnalysisSet --input_for--> AnalysisSet (principal analysis) --files--> cell annotations
    AnalysisSet --documents--> Document (plate map, QC)
    MeasurementSet --related_measurement_sets--> the other modality of a multiome
    MeasurementSet --auxiliary_sets--> AuxiliarySet (MULTI-seq / hashing / guides)
    MeasurementSet --samples--> MultiplexedSample --barcode_map--> TabularFile
                              (barcode to sample mapping) --file_set--> CuratedSet (barcodes)
    Sample --file_sets--> every FileSet made from that pool
    File --quality_metrics--> *QualityMetric (the pipeline's own QC numbers)

and on the prediction side (the Portal's prediction data-model diagram):

    PredictionSet --input_file_sets--> ModelSet --input_file_sets--> experimental AnalysisSet
    PredictionSet --input_file_sets--> AnalysisSet / CuratedSet (the data it was run on)
    PredictionSet --files--> TabularFile (predictions) --derived_from--> ModelFile
    ModelSet --files--> ModelFile --derived_from--> TabularFile (training data, in an AnalysisSet)
    PredictionSet --large_scale_loci_list / large_scale_gene_list--> file (in a CuratedSet)
    PredictionSet / ModelSet --samples / donors--> (virtual) sample or donor

so from raw reads the walk reaches, through input_for, the predictions (e.g.
scE2G element-gene links) that were computed from them, and from a prediction
set it climbs back to its model, the model's training data and the samples.

Earlier IGVFagent code read one hop (MeasurementSet.input_for -> files) and
missed the rest: the principal analysis built on the intermediate one, the
multiome partner and its analyses, the hashing auxiliary set, the barcode to
sample map, and above all the QC metrics the uniform pipeline already
published. This module follows every one of those edges, in both directions,
from any starting accession (measurement / analysis / auxiliary / curated set,
file, sample), and records what it cannot open (403 = unreleased or needs
credentials) instead of silently dropping it.

On top of the graph it builds the answer a user needs before any download:
per analysis product (RNA matrix, ATAC fragments, peaks, cell annotations,
demultiplexing, alignments, predictions ...) the best already-processed file,
where it came from (uniform pipeline or lab), its size and access, the QC the
Portal holds for it, and how much raw data a from-scratch run would need.

Library use:
    from portal_lineage import walk, build_plan
    g = walk("IGVFDS9875NBZW")
    plan = build_plan(g)

CLI (via `igvfagent processed lineage ...`, see processed_first_skill.py).
"""
from __future__ import annotations

import json
import logging
import re
import time
import urllib.parse
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

log = logging.getLogger("portal_lineage")

ACCESSION_RE = re.compile(r"\b(IGVF(?:DS|FI|SM|DO|WF)[0-9]{4}[A-Z]{4}|TSTDS[0-9]+[A-Z]*|TSTFI[0-9]+[A-Z]*)\b")

FILESET_TYPES = {"MeasurementSet", "AnalysisSet", "AuxiliarySet", "CuratedSet", "PredictionSet",
                 "ConstructLibrarySet", "ModelSet", "FileSet", "PseudobulkSet", "SingleCellSet"}
# processed files whose own record is worth a request (QC links, provenance, access)
DETAIL_PRODUCTS = {"rna_matrix", "atac_matrix", "fragments", "alignments", "peaks", "guide_assignment",
                   "predictions", "model", "demultiplexing", "cell_annotations", "perturbation_effects"}
SAMPLE_TYPES = {"MultiplexedSample", "Tissue", "InVitroSystem", "PrimaryCell", "WholeOrganism",
                "TechnicalSample", "Sample"}
DONOR_TYPES = {"HumanDonor", "RodentDonor", "Donor"}

# content_type keyword -> product. Order matters: first match wins.
PRODUCT_RULES: "list[tuple[str, tuple[str, ...]]]" = [
    ("raw_reads", ("reads",)),
    ("seqspec", ("seqspec",)),
    ("onlist", ("onlist", "barcode whitelist", "barcode inclusion")),
    ("demultiplexing", ("barcode to sample", "cell hashing", "sample assignment", "demultiplex")),
    ("cell_annotations", ("cell annotation", "annotated cell", "cell type annotation", "cell metadata")),
    ("rna_matrix", ("cell by gene", "gene count matrix", "sparse gene count", "gene by cell", "transcriptome")),
    ("atac_matrix", ("cell by peak", "peak matrix", "cell by bin")),
    ("guide_assignment", ("guide assignment", "guide quantification", "guide count", "grna")),
    ("perturbation_effects", ("perturbation effect", "differential", "element quantification", "gene quantification",
                              "fold change", "screen results")),
    ("fragments", ("fragments",)),
    ("peaks", ("peaks",)),
    ("alignments", ("alignments",)),
    ("signal", ("signal", "coverage", "bigwig")),
    ("predictions", ("element to gene", "predictions", "variant effect", "variant to element", "regulatory element")),
    ("model", ("snp effect matrix", "model", "weights", "position weight", "pwm", "motif")),
    ("index", ("index",)),
    ("reference", ("genome reference", "transcriptome reference", "genome index")),
]
# Products a downstream analysis starts from, in the order a report lists them.
PRIMARY_PRODUCTS = ["rna_matrix", "atac_matrix", "fragments", "peaks", "cell_annotations", "demultiplexing",
                    "guide_assignment", "perturbation_effects", "predictions", "model", "loci_list", "signal",
                    "alignments"]
SET_PRODUCT = {"PredictionSet": "predictions", "ModelSet": "model"}
# per-product format preference overrides (a predictions TABLE beats a browser track)
PRODUCT_FORMAT_PREF = {"predictions": {"tsv": 0, "bed": 1, "bedpe": 2, "csv": 3, "bigwig": 8, "bigbed": 7},
                       "cell_annotations": {"tsv": 0, "csv": 1}}
FORMAT_PREF = {"h5ad": 0, "h5mu": 1, "hdf5": 2, "h5": 2, "mtx": 3, "tsv": 4, "bed": 4, "bedpe": 4, "csv": 5,
               "rds": 6, "tar": 7, "bigwig": 3, "bam": 8}
# Headline numbers per QC object type (all numeric fields are kept in the JSON).
QC_HEADLINE = [
    "n_reads", "n_pseudoaligned", "p_pseudoaligned", "percentage_reads_on_onlist", "num_barcodes_on_onlist",
    "total_umis", "median_umis_per_barcode", "mean_umis_per_barcode", "n_uniquely_mapped_reads", "n_mapped_reads",
    "pct_duplicates", "tss_enrichment", "frip", "n_barcodes_on_onlist", "n_corrected_barcodes", "n_cells",
    "median_genes_per_cell", "fraction_reads_in_cells", "mapping_rate",
]

Fetch = Callable[[str], "tuple[int, Any]"]


def live_fetch() -> Fetch:
    """Portal GET with the configured IGVF credentials (reads are via raw_data_pipeline.portal_json)."""
    try:
        from igvfagent.raw_data_pipeline import portal_json  # type: ignore
    except Exception:
        from raw_data_pipeline import portal_json  # type: ignore
    return portal_json


def _typ(obj: dict, path: str = "") -> str:
    t = obj.get("@type")
    if isinstance(t, list) and t:
        return str(t[0])
    if isinstance(t, str):
        return t
    seg = (path or obj.get("@id") or "").strip("/").split("/")[0]
    return {"measurement-sets": "MeasurementSet", "analysis-sets": "AnalysisSet", "auxiliary-sets": "AuxiliarySet",
            "curated-sets": "CuratedSet", "prediction-sets": "PredictionSet", "model-sets": "ModelSet",
            "construct-library-sets": "ConstructLibrarySet", "multiplexed-samples": "MultiplexedSample",
            "tissues": "Tissue", "in-vitro-systems": "InVitroSystem", "primary-cells": "PrimaryCell",
            "whole-organisms": "WholeOrganism", "human-donors": "HumanDonor", "rodent-donors": "RodentDonor",
            "documents": "Document", "workflows": "Workflow", "pseudobulk-sets": "PseudobulkSet",
            "single-cell-sets": "SingleCellSet"}.get(
        seg, "File" if seg.endswith("-files") else
        "".join(w.capitalize() for w in seg[:-1].split("-")) if seg.endswith("-sets") else seg)


def _ids(v) -> "list[str]":
    out = []
    for x in v or []:
        if isinstance(x, dict):
            if x.get("@id"):
                out.append(x["@id"])
        elif isinstance(x, str):
            out.append(x)
    return out


def _acc(path: str) -> str:
    m = ACCESSION_RE.search(path or "")
    return m.group(1) if m else (path or "").strip("/").split("/")[-1]


def classify_file(f: dict) -> str:
    ct = str(f.get("content_type") or "").lower()
    for product, keys in PRODUCT_RULES:
        if any(k in ct for k in keys):
            return product
    return "other_processed"


def is_qc(obj: dict) -> bool:
    t = obj.get("@type") or []
    return any("QualityMetric" in str(x) for x in (t if isinstance(t, list) else [t]))


class Graph:
    def __init__(self, root: str):
        self.root = root
        self.nodes: "dict[str, dict]" = {}
        self.edges: "list[tuple[str, str, str]]" = []
        self.blocked: "dict[str, int]" = {}
        self.requests = 0

    def add_edge(self, a: str, b: str, label: str) -> None:
        e = (a, b, label)
        if a and b and e not in self._edge_set:
            self._edge_set.add(e)
            self.edges.append(e)

    _edge_set: "set" = None  # type: ignore

    def to_json(self) -> dict:
        return {"root": self.root, "nodes": self.nodes, "edges": [list(e) for e in self.edges],
                "blocked": self.blocked, "requests": self.requests}


def _file_node(f: dict, parent_set: str = "") -> dict:
    return {
        "kind": "file", "type": _typ(f), "accession": f.get("accession") or _acc(f.get("@id", "")),
        "@id": f.get("@id"), "content_type": f.get("content_type") or "", "file_format": f.get("file_format") or "",
        "file_size": f.get("file_size"), "controlled_access": f.get("controlled_access"),
        "status": f.get("status") or "", "href": f.get("href") or "", "s3_uri": f.get("s3_uri") or "",
        "anvil_url": f.get("anvil_url") or "", "md5sum": f.get("md5sum") or "",
        "illumina_read_type": f.get("illumina_read_type") or "", "summary": f.get("summary") or "",
        "file_set": parent_set or (f.get("file_set", {}).get("@id") if isinstance(f.get("file_set"), dict)
                                   else f.get("file_set") or ""),
        "derived_from": _ids(f.get("derived_from")), "seqspecs": _ids(f.get("seqspecs")),
        "product": classify_file(f), "quality_metrics": _ids(f.get("quality_metrics")),
        "assembly": f.get("assembly") or "", "transcriptome_annotation": f.get("transcriptome_annotation") or "",
        "filtered": f.get("filtered"), "workflows": [w.get("name") or w.get("@id") if isinstance(w, dict) else w
                                                    for w in (f.get("workflows") or [])],
        "software": [s.get("summary") for s in ((f.get("analysis_step_version") or {}).get("software_versions") or [])
                     if isinstance(s, dict)] if isinstance(f.get("analysis_step_version"), dict) else [],
    }


# Walk directions. The graph is followed with a relation so the walk stays on
# the provenance of the starting object instead of spreading through shared
# pools: a principal analysis takes every sub-pool's intermediate analysis as
# input, and a multiplexed sample lists every file set made from the pool, so
# an undirected walk from one sub-pool reaches all of them (442 objects and
# 812 GB of "raw reads" from IGVFDS9875NBZW on 2026-09-23, with a matrix from
# a different sub-pool ranked first).
#   root    the starting object: up (inputs) and down (outputs)
#   family  its multiome partner / auxiliary sets: down only
#   down    built from root/family: further down only; a prediction set's
#           ModelSet input is followed up (model -> training data)
#   up      inputs of root: further up only
#   ref     samples, barcode maps, curated barcode sets, documents: recorded,
#           not expanded into sibling datasets
#   model_use  prediction sets that APPLY a model trained on this data to other
#           data (scE2G trained on K562 CRISPR, run on 100 other cell types):
#           recorded, not treated as results of this data
RELATION_RANK = {"root": 0, "family": 1, "down": 2, "up": 3, "model_use": 4, "ref": 5, "sibling": 6}
QC_PRODUCTS = {"rna_matrix", "atac_matrix", "fragments", "alignments", "peaks", "guide_assignment"}


def walk(accession: str, fetch: "Optional[Fetch]" = None, max_depth: int = 6, max_nodes: int = 300,
         fetch_qc: bool = True, reverse_search: bool = True, fanout: int = 25, siblings: bool = False,
         workers: int = 8) -> Graph:
    """Directed breadth-first walk of the Portal graph from `accession` (or an @id path)."""
    fetch = fetch or live_fetch()
    g = Graph(accession)
    g._edge_set = set()
    cache: "dict[str, tuple[int, Any]]" = {}

    def get(path: str, frame: str) -> "tuple[int, Any]":
        key = f"{path}|{frame}"
        if key not in cache:
            sep = "&" if "?" in path else "?"
            url = f"{path}{sep}format=json" + (f"&frame={frame}" if frame else "")
            g.requests += 1
            cache[key] = fetch(url)
        return cache[key]

    start = accession if accession.startswith("/") else f"/{accession}/"
    queue: "deque[tuple[str, int, str]]" = deque([(start, 0, "root")])
    seen: "dict[str, str]" = {}
    root_id = None
    relation: "dict[str, str]" = {}

    def enqueue(path: str, depth: int, rel: str) -> None:
        if not path or depth > max_depth or len(seen) >= max_nodes:
            return
        old = seen.get(path) or relation.get(path)
        if old is not None and RELATION_RANK.get(old, 9) <= RELATION_RANK.get(rel, 9):
            return
        queue.append((path, depth, rel))

    def capped(src: str, items: "list[str]", field: str) -> "list[str]":
        """Shared hubs (a curated training set feeding hundreds of models) are recorded, not all walked."""
        if len(items) > fanout:
            g.nodes.setdefault(src, {})
            g.nodes[src].setdefault("truncated", {})[field] = len(items)
            return items[:fanout]
        return items

    def stub(path: str, rel: str) -> None:
        if path not in g.nodes:
            g.nodes[path] = {"kind": "stub", "type": _typ({}, path), "accession": _acc(path), "@id": path, "relation": rel}

    def prefetch() -> None:
        """Fetch the pending frontier in parallel; processing stays sequential and deterministic."""
        todo = []
        for pth, _d, _r in list(queue)[:64]:
            fr = "embedded" if (_typ({}, pth) in FILESET_TYPES or pth == start) else "object"
            if f"{pth}|{fr}" not in cache and (pth, fr) not in todo:
                todo.append((pth, fr))
        if len(todo) < 2 or workers < 2:
            return
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for (pth, fr), res in zip(todo, ex.map(lambda x: _raw(x[0], x[1]), todo)):
                cache[f"{pth}|{fr}"] = res

    def _raw(path: str, frame: str):
        sep = "&" if "?" in path else "?"
        g.requests += 1
        return fetch(f"{path}{sep}format=json" + (f"&frame={frame}" if frame else ""))

    while queue:
        prefetch()
        path, depth, rel = queue.popleft()
        prev_rel = seen.get(path)
        if prev_rel is not None and RELATION_RANK.get(prev_rel, 9) <= RELATION_RANK.get(rel, 9):
            continue
        seen[path] = rel
        typ_guess = _typ({}, path)
        frame = "embedded" if (typ_guess in FILESET_TYPES or path == start) else "object"
        status, obj = get(path, frame)
        if status != 200 or not isinstance(obj, dict) or obj.get("status") == "error":
            code = status if status else 0
            g.blocked[path] = code
            g.nodes[path] = {"kind": "blocked", "type": typ_guess, "accession": _acc(path), "@id": path,
                             "http_status": code, "relation": rel}
            continue
        if obj.get("@graph") and not obj.get("@id"):
            obj = obj["@graph"][0]
        oid = obj.get("@id") or path
        if path == start:
            root_id = oid
            g.root = oid
        seen[oid] = rel
        relation[oid] = rel
        t = _typ(obj, oid)

        if t in FILESET_TYPES or "FileSet" in (obj.get("@type") or []):
            prev = g.nodes.get(oid) or {}
            node = {"kind": "fileset", "type": t, "accession": obj.get("accession") or _acc(oid), "@id": oid,
                    "relation": rel,
                    "file_set_type": obj.get("file_set_type") or "", "status": obj.get("status") or "",
                    "summary": (obj.get("summary") or "")[:300], "aliases": obj.get("aliases") or [],
                    "uniform_pipeline_status": obj.get("uniform_pipeline_status") or "",
                    "controlled_access": obj.get("controlled_access"),
                    "assay": ((obj.get("assay_term") or {}).get("term_name") if isinstance(obj.get("assay_term"), dict) else "")
                             or ", ".join(obj.get("assay_titles") or []),
                    "preferred_assay_titles": obj.get("preferred_assay_titles") or [],
                    "lab": (obj.get("lab") or {}).get("title") if isinstance(obj.get("lab"), dict) else obj.get("lab"),
                    "workflows": [{"name": w.get("name"), "accession": w.get("accession"),
                                   "uniform": w.get("uniform_pipeline"), "version": w.get("workflow_version")}
                                  for w in (obj.get("workflows") or []) if isinstance(w, dict)],
                    "pipeline_parameters": obj.get("pipeline_parameters"), "depth": depth,
                    "n_files": len(obj.get("files") or []), "doi": obj.get("doi") or "",
                    "multiome_size": obj.get("multiome_size"), "is_on_anvil": obj.get("is_on_anvil"),
                    "scope": obj.get("scope"), "model_name": obj.get("model_name"),
                    "model_version": obj.get("model_version"), "prediction_objects": obj.get("prediction_objects"),
                    "software_versions": [x.get("summary") or x.get("@id") if isinstance(x, dict) else x
                                          for x in obj.get("software_versions") or []],
                    "n_input_for": len(obj.get("input_for") or []),
                    "n_input_file_sets": len(obj.get("input_file_sets") or []),
                    "truncated": prev.get("truncated", {})}
            g.nodes[oid] = node
            expand_files = rel in ("root", "family", "down", "up")
            if rel == "model_use":
                # keep only whether it also used this data directly; its files are listed, not fetched
                node["data_inputs"] = [u for u in _ids(obj.get("input_file_sets")) if _typ({}, u) != "ModelSet"]
            for f in obj.get("files") or []:
                if isinstance(f, str):
                    f = {"@id": f}
                fid = f.get("@id")
                fn = _file_node(f, oid)
                fn["relation"] = rel
                if t in SET_PRODUCT and fn["product"] not in ("raw_reads", "index", "seqspec"):
                    fn["product"] = SET_PRODUCT[t]
                g.nodes.setdefault(fid, fn)
                g.nodes[fid]["relation"] = min((g.nodes[fid].get("relation") or rel, rel), key=lambda r: RELATION_RANK.get(r, 9))
                g.add_edge(oid, fid, "files")
                if expand_files and fn["product"] in DETAIL_PRODUCTS and fetch_qc and rel != "model_use":
                    enqueue(fid, depth + 1, rel)
                for sq in fn["seqspecs"]:
                    g.add_edge(fid, sq, "seqspecs")
            for fid in _ids(obj.get("onlist_files")):
                g.add_edge(oid, fid, "onlist_files")
                enqueue(fid, depth + 1, "ref")
            # downstream: what was built from this set. A model trained on this data is
            # downstream; the predictions it makes on OTHER data are "model_use".
            if rel in ("root", "family", "down"):
                for d in capped(oid, _ids(obj.get("input_for")), "input_for"):
                    g.add_edge(oid, d, "input_for")
                    enqueue(d, depth + 1, "model_use" if t == "ModelSet" and rel == "down" else "down")
                if reverse_search and not obj.get("input_for") and t != "CuratedSet":
                    st, res = get(f"/search/?type=FileSet&input_file_sets.@id={urllib.parse.quote(oid)}&limit=50", "")
                    for hit in (res or {}).get("@graph", []) if st == 200 else []:
                        if hit.get("@id") and hit["@id"] != oid:
                            g.add_edge(oid, hit["@id"], "input_for")
                            enqueue(hit["@id"], depth + 1, "down")
            elif obj.get("input_for"):
                node["n_input_for"] = len(obj.get("input_for") or [])
            # upstream: what this set was built from
            for u in _ids(obj.get("input_file_sets")):
                g.add_edge(u, oid, "input_for")
                ut = _typ({}, u)
                if rel in ("root", "up"):
                    enqueue(u, depth + 1, "up")
                elif rel == "down" and t in ("PredictionSet", "ModelSet") and ut in ("ModelSet", "AnalysisSet", "CuratedSet") \
                        and u not in relation:
                    # a prediction's model and the model's training data are its provenance
                    enqueue(u, depth + 1, "up" if ut == "ModelSet" or t == "ModelSet" else "sibling")
                    if ut not in ("ModelSet",) and t != "ModelSet":
                        stub(u, "sibling")
                else:
                    stub(u, "sibling")
            if rel in ("root",):
                for a in _ids(obj.get("auxiliary_sets")):
                    g.add_edge(oid, a, "auxiliary_sets")
                    enqueue(a, depth + 1, "family")
                for grp in obj.get("related_measurement_sets") or []:
                    for m in _ids((grp or {}).get("measurement_sets")):
                        g.add_edge(oid, m, "related_" + str((grp or {}).get("series_type") or "measurement_set"))
                        enqueue(m, depth + 1, "family")
                for m in _ids(obj.get("measurement_sets")):          # AuxiliarySet -> MeasurementSet
                    g.add_edge(m, oid, "auxiliary_sets")
                    enqueue(m, depth + 1, "family")
            else:
                for a in _ids(obj.get("auxiliary_sets")):
                    g.add_edge(oid, a, "auxiliary_sets")
                for grp in obj.get("related_measurement_sets") or []:
                    for m in _ids((grp or {}).get("measurement_sets")):
                        g.add_edge(oid, m, "related_" + str((grp or {}).get("series_type") or "measurement_set"))
            for c in _ids(obj.get("control_file_sets")):
                g.add_edge(oid, c, "control_file_sets")
                if rel == "root":
                    enqueue(c, depth + 1, "up")
            for fld in ("large_scale_loci_list", "large_scale_gene_list"):
                lf = obj.get(fld)
                lf = lf.get("@id") if isinstance(lf, dict) else lf
                if lf:
                    g.add_edge(oid, lf, fld)
                    g.nodes.setdefault(lf, {"kind": "file", "type": _typ({}, lf), "accession": _acc(lf), "@id": lf,
                                            "product": "loci_list", "content_type": fld.replace("_", " "),
                                            "file_set": "", "controlled_access": None, "relation": "ref"})
                    enqueue(lf, depth + 1, "ref")
            for fld in ("small_scale_loci_list", "small_scale_gene_list"):
                if obj.get(fld):
                    node[fld + "_n"] = len(obj.get(fld) or [])
            if rel in ("root", "family", "down", "up"):
                for sid in capped(oid, _ids(obj.get("samples")), "samples"):
                    g.add_edge(oid, sid, "samples")
                    enqueue(sid, depth + 1, "ref")
                for d in obj.get("documents") or []:
                    did = d.get("@id") if isinstance(d, dict) else d
                    g.add_edge(oid, did, "documents")
                    enqueue(did, depth + 1, "ref")
            for d in (obj.get("donors") or []) if rel in ("root", "down", "up") else []:
                if isinstance(d, str):
                    d = {"@id": d, "accession": _acc(d)}
                if isinstance(d, dict) and d.get("@id"):
                    g.nodes.setdefault(d["@id"], {"kind": "donor", "type": _typ(d, d["@id"]), "accession": d.get("accession"),
                                                  "@id": d["@id"], "sex": d.get("sex"), "taxa": d.get("taxa"),
                                                  "aliases": d.get("aliases") or [], "relation": "ref"})
                    g.add_edge(oid, d["@id"], "donors")
        elif is_qc(obj):
            nums = {k: v for k, v in obj.items() if isinstance(v, (int, float)) and not isinstance(v, bool)
                    and k not in ("schema_version",)}
            attach = {k: v.get("href") for k, v in obj.items() if isinstance(v, dict) and v.get("href")}
            g.nodes[oid] = {"kind": "qc", "type": t, "@id": oid, "accession": _acc(oid), "relation": rel,
                            "description": obj.get("description") or "", "metrics": nums, "attachments": attach,
                            "quality_metric_of": _ids(obj.get("quality_metric_of")), "status": obj.get("status")}
            for f in _ids(obj.get("quality_metric_of")):
                g.add_edge(f, oid, "quality_metrics")
        elif t in SAMPLE_TYPES or "Sample" in (obj.get("@type") or []):
            g.nodes[oid] = {"kind": "sample", "type": t, "@id": oid, "accession": obj.get("accession") or _acc(oid),
                            "relation": rel, "summary": (obj.get("summary") or "")[:300], "status": obj.get("status"),
                            "sample_terms": [x.get("term_name") if isinstance(x, dict) else x
                                             for x in obj.get("sample_terms") or []],
                            "n_multiplexed_samples": len(obj.get("multiplexed_samples") or []),
                            "n_donors": len(obj.get("donors") or []), "cellular_sub_pool": obj.get("cellular_sub_pool"),
                            "multiplexing_methods": obj.get("multiplexing_methods") or [],
                            "taxa": obj.get("taxa"), "n_file_sets": len(obj.get("file_sets") or [])}
            bm = obj.get("barcode_map")
            bm = bm.get("@id") if isinstance(bm, dict) else bm
            if bm:
                g.add_edge(oid, bm, "barcode_map")
                enqueue(bm, depth + 1, "ref")
            sample_is_root = rel == "root"
            for fs in capped(oid, _ids(obj.get("file_sets")), "file_sets"):
                g.add_edge(oid, fs, "file_sets")
                if sample_is_root:
                    enqueue(fs, depth + 1, "family")
                elif siblings:
                    enqueue(fs, depth + 1, "sibling")
                else:
                    stub(fs, "sibling")
            for ms in _ids(obj.get("multiplexed_samples")):
                g.add_edge(oid, ms, "multiplexed_samples")
        elif t == "Document":
            att = obj.get("attachment") or {}
            g.nodes[oid] = {"kind": "document", "type": t, "@id": oid, "accession": _acc(oid), "relation": rel,
                            "document_type": obj.get("document_type") or "", "description": obj.get("description") or "",
                            "attachment": (oid.rstrip("/") + "/" + att.get("href")) if isinstance(att, dict) and att.get("href") else ""}
        else:  # a file object
            fn = _file_node(obj)
            prev = g.nodes.get(oid) or {}
            fn["file_set"] = prev.get("file_set") or fn["file_set"]
            if prev.get("product") in ("predictions", "model", "loci_list"):
                fn["product"] = prev["product"]
            merged = {**prev, **{k: v for k, v in fn.items() if v not in (None, "", [])}}
            merged["relation"] = prev.get("relation") or rel
            g.nodes[oid] = merged
            if fn["file_set"]:
                g.add_edge(fn["file_set"], oid, "files")
                if rel == "root":                   # started from a file: its set is the real root
                    enqueue(fn["file_set"], depth + 1, "root")
                elif rel == "ref" and fn["product"] in ("demultiplexing", "onlist", "loci_list"):
                    enqueue(fn["file_set"], depth + 1, "ref")   # barcode map -> curated barcode set
            for d in fn["derived_from"]:
                g.add_edge(d, oid, "derived_from")
                if not d.startswith("/sequence-files/") and fn["product"] in ("predictions", "model", "loci_list"):
                    enqueue(d, depth + 1, "up")
            for sq in fn["seqspecs"]:
                g.add_edge(oid, sq, "seqspecs")
            qms = fn["quality_metrics"]
            if fetch_qc and not qms and fn["product"] in QC_PRODUCTS and rel in ("root", "family", "down"):
                st, res = get(f"/search/?type=QualityMetric&quality_metric_of={urllib.parse.quote(oid)}&limit=20", "object")
                qms = [h.get("@id") for h in (res or {}).get("@graph", [])] if st == 200 else []
                for h in (res or {}).get("@graph", []) if st == 200 else []:
                    cache[f"{h.get('@id')}|object"] = (200, h)
            for q in qms:
                g.add_edge(oid, q, "quality_metrics")
                if fetch_qc and rel in ("root", "family", "down"):
                    enqueue(q, depth + 1, rel)
    if root_id is None:
        g.root = start
    g.cap_reached = len(seen) >= max_nodes
    return g


# ---------------------------------------------------------------------------
# Plan: what to use instead of raw reads
# ---------------------------------------------------------------------------

def _gb(n) -> float:
    try:
        return float(n) / 1e9
    except (TypeError, ValueError):
        return 0.0


def _access(f: dict, have_credentials: bool) -> str:
    ca = f.get("controlled_access")
    if ca is True:
        return "controlled" + ("" if have_credentials else " (needs IGVF credentials + approved DUA)")
    if ca is False:
        return "public"
    return "not stated"


def build_plan(g: Graph, wants: "Optional[list[str]]" = None, have_credentials: bool = False) -> dict:
    nodes = g.nodes
    sets = {k: v for k, v in nodes.items() if v.get("kind") == "fileset"}
    files = {k: v for k, v in nodes.items() if v.get("kind") == "file"}
    qcs = {k: v for k, v in nodes.items() if v.get("kind") == "qc"}
    qc_by_file: "dict[str, list[dict]]" = defaultdict(list)
    for q in qcs.values():
        for f in q.get("quality_metric_of") or []:
            qc_by_file[f].append(q)
    for a, b, lab in g.edges:
        if lab == "quality_metrics" and b in qcs and qcs[b] not in qc_by_file[a]:
            qc_by_file[a].append(qcs[b])

    def set_rank(sid: str) -> "tuple[int, int]":
        s = sets.get(sid) or {}
        uniform = s.get("uniform_pipeline_status") == "completed" or any(w.get("uniform") for w in s.get("workflows") or [])
        principal = s.get("file_set_type") == "principal analysis"
        return (0 if uniform else 1, 0 if principal else 1)

    products: "dict[str, list[dict]]" = defaultdict(list)
    for fid, f in files.items():
        p = f.get("product") or "other_processed"
        if p in ("raw_reads", "seqspec", "onlist", "index", "reference"):
            continue
        s = sets.get(f.get("file_set") or "", {})
        rec = {
            "accession": f.get("accession"), "@id": fid, "content_type": f.get("content_type"),
            "file_format": f.get("file_format"), "size_gb": round(_gb(f.get("file_size")), 3),
            "access": _access(f, have_credentials), "controlled_access": f.get("controlled_access"),
            "status": f.get("status"), "href": f.get("href"), "file_set": s.get("accession") or _acc(f.get("file_set", "")),
            "file_set_type": s.get("file_set_type") or "", "uniform_pipeline_status": s.get("uniform_pipeline_status") or "",
            "workflows": [w.get("name") for w in s.get("workflows") or [] if w.get("name")],
            "software": f.get("software") or [], "annotation": f.get("transcriptome_annotation") or f.get("assembly") or "",
            "filtered": f.get("filtered"),
            "qc": [{"type": q["type"], "metrics": {k: q["metrics"][k] for k in QC_HEADLINE if k in q["metrics"]},
                    "all_metrics": q["metrics"], "attachments": q.get("attachments") or {}} for q in qc_by_file.get(fid, [])],
            "_rank": set_rank(f.get("file_set") or ""),
            "relation": f.get("relation") or (sets.get(f.get("file_set") or "", {}).get("relation")) or "ref",
        }
        rec["_product"] = p
        products[p].append(rec)

    def frank(r):
        blocked = 1 if (r["controlled_access"] is True and not have_credentials) else 0
        direct = 0 if r["relation"] in ("root", "family", "down") else 1
        fmt = str(r["file_format"]).lower()
        pref = PRODUCT_FORMAT_PREF.get(r.get("_product") or "", {}).get(fmt, FORMAT_PREF.get(fmt, 9))
        return (direct, r["_rank"][0], blocked, r["_rank"][1], pref, r["size_gb"] or 1e9)

    ranked = {}
    for p, lst in products.items():
        lst = sorted(lst, key=frank)
        for r in lst:
            r.pop("_rank", None)
            r.pop("_product", None)
        ranked[p] = lst
    raw = [f for f in files.values() if f.get("product") == "raw_reads"
           and (sets.get(f.get("file_set") or "", {}).get("relation") in ("root", "family", "up"))
           and sets.get(f.get("file_set") or "", {}).get("type") in ("MeasurementSet", "AuxiliarySet")]
    raw_by_set: "dict[str, dict]" = defaultdict(lambda: {"n": 0, "gb": 0.0, "controlled": 0})
    for f in raw:
        k = sets.get(f.get("file_set") or "", {}).get("accession") or _acc(f.get("file_set") or "")
        raw_by_set[k]["n"] += 1
        raw_by_set[k]["gb"] += _gb(f.get("file_size"))
        raw_by_set[k]["controlled"] += int(f.get("controlled_access") is True)
    order = [p for p in PRIMARY_PRODUCTS if p in ranked] + sorted(p for p in ranked if p not in PRIMARY_PRODUCTS)
    if wants:
        order = [p for p in order if p in wants]
    start_here = []
    for p in order:
        best = ranked[p][0]
        if best["relation"] not in ("root", "family", "down") and p not in ("model",):
            continue   # only inputs / references carry this product: not a result of this data
        start_here.append({"product": p, "relation": best["relation"], **{k: best[k] for k in ("accession", "content_type", "file_format", "size_gb",
                                                                 "access", "file_set", "file_set_type",
                                                                 "uniform_pipeline_status", "workflows", "href")},
                           "qc_headline": (best["qc"][0]["metrics"] if best["qc"] else {}),
                           "alternatives": [r["accession"] for r in ranked[p][1:4]]})
    siblings = sorted({v.get("accession") for v in nodes.values() if v.get("relation") == "sibling"
                       and v.get("kind") in ("stub", "fileset")} - {None})
    analysis_sets = [{"accession": s["accession"], "file_set_type": s.get("file_set_type"), "relation": s.get("relation"),
                      "uniform_pipeline_status": s.get("uniform_pipeline_status"), "status": s.get("status"),
                      "workflows": s.get("workflows"), "summary": s.get("summary"), "n_files": s.get("n_files")}
                     for s in sets.values() if s.get("type") == "AnalysisSet" and s.get("relation") != "sibling"]
    measurement_sets = [{"accession": s["accession"], "type": s["type"], "assay": s.get("assay"),
                         "file_set_type": s.get("file_set_type"), "status": s.get("status"),
                         "controlled_access": s.get("controlled_access")}
                        for s in sets.values() if s.get("type") in ("MeasurementSet", "AuxiliarySet")
                        and s.get("relation") in ("root", "family", "up")]
    blocked = [{"@id": k, "accession": _acc(k), "type": _typ({}, k), "http_status": v,
                "meaning": ("exists but not visible: unreleased, or needs credentials with access"
                            if v == 403 else "not found" if v == 404 else f"request failed ({v})")}
               for k, v in g.blocked.items()]
    raw_gb = sum(v["gb"] for v in raw_by_set.values())
    use_gb = sum(x["size_gb"] or 0 for x in start_here if x["product"] in ("rna_matrix", "atac_matrix", "fragments",
                                                                         "peaks", "cell_annotations", "demultiplexing"))
    configs = [{"accession": f.get("accession"), "product": f.get("product"), "content_type": f.get("content_type"),
                "file_set": sets.get(f.get("file_set") or "", {}).get("accession") or _acc(f.get("file_set") or ""),
                "href": f.get("href"), "size_gb": round(_gb(f.get("file_size")), 6), "access": _access(f, have_credentials)}
               for f in files.values() if f.get("product") in ("seqspec", "onlist")]
    samples = [v for v in nodes.values() if v.get("kind") == "sample"]
    documents = [v for v in nodes.values() if v.get("kind") == "document"]
    pred_models = [{"accession": s.get("accession"), "type": s.get("type"), "file_set_type": s.get("file_set_type"),
                    "scope": s.get("scope"), "model_name": s.get("model_name"), "model_version": s.get("model_version"),
                    "software_versions": s.get("software_versions") or [], "summary": s.get("summary"),
                    "status": s.get("status"), "n_files": s.get("n_files"),
                    "inputs": [_acc(a) for a, b, lab in g.edges if b == k and lab == "input_for"],
                    "relation": s.get("relation"),
                    "truncated": s.get("truncated") or {}}
                   for k, s in sets.items() if s.get("type") in ("PredictionSet", "ModelSet")]
    root = nodes.get(g.root, {})
    verdict = ("processed" if any(x["product"] in ("rna_matrix", "atac_matrix", "fragments", "peaks", "cell_annotations",
                                                   "perturbation_effects", "predictions") for x in start_here)
               else "raw_only" if raw else "nothing_found")
    return {"root": {"@id": g.root, "accession": root.get("accession") or _acc(g.root), "type": root.get("type"),
                     "summary": root.get("summary"), "assay": root.get("assay"),
                     "uniform_pipeline_status": root.get("uniform_pipeline_status")},
            "verdict": verdict, "start_here": start_here, "products": ranked, "analysis_sets": analysis_sets,
            "measurement_sets": measurement_sets, "raw_by_set": dict(raw_by_set), "raw_gb": round(raw_gb, 2),
            "processed_gb_to_download": round(use_gb, 2), "configs": configs, "samples": samples,
            "documents": documents, "blocked": blocked, "have_credentials": have_credentials,
            "predictions_and_models": pred_models, "siblings": siblings,
            "truncated": {s.get("accession"): s.get("truncated") for s in sets.values() if s.get("truncated")},
            "n_nodes": len(nodes), "n_edges": len(g.edges), "requests": g.requests,
            "node_cap_reached": bool(getattr(g, "cap_reached", False))}


# ---------------------------------------------------------------------------
# Report + figure
# ---------------------------------------------------------------------------

PRODUCT_LABEL = {"rna_matrix": "RNA cell x gene matrix", "atac_matrix": "ATAC cell x peak matrix",
                 "fragments": "ATAC fragments", "peaks": "peaks", "cell_annotations": "cell annotations",
                 "demultiplexing": "sample demultiplexing (hashing / barcode map)", "guide_assignment": "guide assignment",
                 "perturbation_effects": "perturbation effects", "predictions": "predictions", "signal": "signal tracks",
                 "alignments": "alignments (BAM)", "other_processed": "other processed files"}
NEXT_STEP = {
    "rna_matrix": "igvfagent sc-analyze pipeline --input <downloaded .h5ad>",
    "fragments": "igvfagent igvf-sc-pipeline tss-enrichment --fragments <downloaded fragments> ...",
    "atac_matrix": "igvfagent sc-analyze pipeline --input <downloaded matrix>",
    "cell_annotations": "join the annotation table to the matrix obs by barcode",
    "demultiplexing": "use as the barcode -> sample (donor) map for per-donor analyses",
}


def _fmt_metric(v) -> str:
    if isinstance(v, float):
        return f"{v:,.3f}" if abs(v) < 1000 else f"{v:,.0f}"
    if isinstance(v, int):
        return f"{v:,}"
    return str(v)


def write_report(plan: dict, g: Graph, out: Path, figure: "Optional[Path]" = None) -> Path:
    r = plan["root"]
    lines = [f"# What the IGVF Portal already holds for {r['accession']}", "",
             f"{r.get('type') or ''}: {r.get('summary') or ''}", ""]
    v = plan["verdict"]
    if v == "processed":
        lines += [f"**Processed results already exist.** Start from them instead of the raw reads: about "
                  f"{plan['processed_gb_to_download']:.1f} GB of processed files versus {plan['raw_gb']:.1f} GB of reads "
                  f"reachable from this accession.", ""]
    elif v == "raw_only":
        lines += [f"**Only raw data were found** ({plan['raw_gb']:.1f} GB of reads). No analysis set with processed "
                  f"outputs is reachable yet, so the pipeline has to be run from the reads.", ""]
    else:
        lines += ["**Nothing downloadable was reachable from this accession.**", ""]
    if plan["start_here"]:
        lines += ["## Start here", "", "| product | file | format | size | access | from | pipeline |",
                  "|---|---|---|---|---|---|---|"]
        for s in plan["start_here"]:
            pipe = "uniform, completed" if s["uniform_pipeline_status"] == "completed" else \
                (", ".join(s["workflows"]) or s["file_set_type"] or "lab-submitted")
            lines.append(f"| {PRODUCT_LABEL.get(s['product'], s['product'])} | {s['accession']} | {s['file_format']} | "
                         f"{s['size_gb']:.2f} GB | {s['access']} | {s['file_set']} ({s['file_set_type']}) | {pipe} |")
        lines.append("")
    qc_rows = []
    for p, lst in plan["products"].items():
        for rec in lst:
            for q in rec["qc"]:
                if q["metrics"]:
                    qc_rows.append((rec["accession"], p, q["type"], q["metrics"]))
    if qc_rows:
        lines += ["## QC the Portal already computed", "",
                  "These numbers come from the pipeline's QualityMetric objects; no reads were processed here.", ""]
        seen = set()
        for acc, p, typ, m in qc_rows:
            key = (typ, tuple(sorted(m.items())))
            if key in seen:
                continue
            seen.add(key)
            lines.append(f"- **{typ}** (on {acc}, {PRODUCT_LABEL.get(p, p)}): "
                         + "; ".join(f"{k} {_fmt_metric(val)}" for k, val in m.items()))
        lines.append("")
    if plan["analysis_sets"]:
        lines += ["## Analysis sets", "", "| accession | type | uniform pipeline | workflow | files |", "|---|---|---|---|---|"]
        for s in plan["analysis_sets"]:
            wf = ", ".join(f"{w.get('name')} {w.get('version') or ''}".strip() for w in s.get("workflows") or [])
            lines.append(f"| {s['accession']} | {s['file_set_type']} | {s['uniform_pipeline_status'] or '-'} | {wf or '-'} | {s['n_files']} |")
        lines.append("")
    direct_pm = [m for m in plan.get("predictions_and_models") or [] if m.get("relation") != "model_use"]
    use_pm = [m for m in plan.get("predictions_and_models") or [] if m.get("relation") == "model_use"]
    if use_pm:
        lines += [f"{len(use_pm)} prediction set(s) apply a model trained on this data to other data "
                  f"(not results of this dataset): " + ", ".join(m["accession"] for m in use_pm[:20])
                  + (" ..." if len(use_pm) > 20 else ""), ""]
    if direct_pm:
        lines += ["## Predictions and models built on this data", "",
                  "| set | type | kind | scope / model | inputs | summary |", "|---|---|---|---|---|---|"]
        for m in sorted(direct_pm, key=lambda x: (x["type"], x["accession"] or "")):
            what = m.get("scope") or " ".join(x for x in (m.get("model_name"), m.get("model_version")) if x) or ""
            lines.append(f"| {m['accession']} | {m['type']} | {m.get('file_set_type') or ''} | {what} | "
                         f"{', '.join(m['inputs'][:4])} | {(m.get('summary') or '')[:120]} |")
        lines.append("")
    if plan.get("truncated"):
        lines += ["Shared hubs were not walked in full: " + "; ".join(
            f"{k} {', '.join(f'{f} {n}' for f, n in v.items())}" for k, v in plan["truncated"].items())
            + " links (raise --fanout to follow more).", ""]
    if plan.get("siblings"):
        lines += [f"Sibling datasets from the same sample pool or used alongside this one (not walked, "
                  f"`--siblings` to include): {', '.join(plan['siblings'][:30])}"
                  + (" ..." if len(plan["siblings"]) > 30 else ""), ""]
    if plan["measurement_sets"]:
        lines += ["## Raw data", "", "| set | type | assay | reads | size | controlled |", "|---|---|---|---|---|---|"]
        for m in plan["measurement_sets"]:
            rb = plan["raw_by_set"].get(m["accession"], {"n": 0, "gb": 0, "controlled": 0})
            lines.append(f"| {m['accession']} | {m['type']} ({m.get('file_set_type')}) | {m.get('assay') or ''} | "
                         f"{rb['n']} | {rb['gb']:.1f} GB | {rb['controlled']} of {rb['n']} |")
        lines.append("")
    if plan["configs"]:
        lines += ["## Configuration (seqspec, onlists)", ""]
        for c in plan["configs"]:
            lines.append(f"- {c['accession']} {c['content_type']} of {c['file_set']} ({c['access']})")
        lines.append("")
    if plan["samples"]:
        lines += ["## Samples", ""]
        for s in plan["samples"]:
            lines.append(f"- {s['accession']} {s['type']}: {s.get('summary') or ''}")
        lines.append("")
    if plan["documents"]:
        lines += ["## Documents", ""]
        for d in plan["documents"]:
            lines.append(f"- {d['accession']} {d.get('document_type')}: {d.get('description') or ''}")
        lines.append("")
    if plan["blocked"]:
        lines += ["## Linked but not visible", ""]
        docs = [b for b in plan["blocked"] if b["type"] == "Document"]
        for b in plan["blocked"]:
            if b["type"] != "Document":
                lines.append(f"- {b['accession']} ({b['type']}): HTTP {b['http_status']}, {b['meaning']}")
        if docs:
            lines.append(f"- {len(docs)} document(s) (protocols, plate maps): HTTP 403, not released")
        lines += ["", "Supply IGVF credentials with access (IGVF_ACCESS_KEY / IGVF_SECRET_ACCESS_KEY, then "
                  "`igvfagent auth-check`) to see these.", ""]
    if plan["start_here"]:
        lines += ["## Next steps", "", "```bash",
                  f"igvfagent processed fetch {r['accession']} --want "
                  + ",".join(s["product"] for s in plan["start_here"][:3])]
        for s in plan["start_here"][:3]:
            if s["product"] in NEXT_STEP:
                lines.append(NEXT_STEP[s["product"]])
        lines += ["```", ""]
    lines += ["## Graph", "", f"{plan['n_nodes']} objects, {plan['n_edges']} links, {plan['requests']} Portal requests."
              + (" The object cap was reached, so distant links may be missing; rerun with a larger --max-nodes."
                 if plan.get("node_cap_reached") else ""), ""]
    if figure:
        lines += [f"![lineage]({figure.name})", ""]
    lines += ["```mermaid", mermaid(g), "```", ""]
    out.write_text("\n".join(lines))
    return out


def _short(n: dict) -> str:
    if n.get("kind") == "fileset":
        return f"{n['accession']}\\n{n['type']}\\n{n.get('file_set_type') or ''}"
    if n.get("kind") == "sample":
        return f"{n['accession']}\\n{n['type']}"
    if n.get("kind") == "blocked":
        return f"{n['accession']}\\n{n['type']} (HTTP {n.get('http_status')})"
    return n.get("accession") or n.get("@id", "")


def _drawable(g: Graph) -> "set[str]":
    """Provenance only: root/family/down/up sets and their files, plus samples, barcode maps and
    curated barcode sets hanging directly off the root or its family. Siblings become one count."""
    core = {k for k, n in g.nodes.items() if n.get("kind") == "fileset" and n.get("relation") in ("root", "family", "down", "up")}
    core |= {k for k, n in g.nodes.items() if n.get("kind") == "blocked" and n.get("relation") in ("down", "up", "family")}
    keep = set(core)
    for a, b, lab in g.edges:
        if a in core and lab == "files":
            keep.add(b)
    anchor = {k for k in core if g.nodes[k].get("relation") in ("root", "family")}
    for a, b, lab in g.edges:
        if a in anchor and lab in ("samples", "documents", "onlist_files", "large_scale_loci_list", "large_scale_gene_list"):
            keep.add(b)
    for a, b, lab in list(g.edges):
        if a in keep and lab == "barcode_map":
            keep.add(b)
            fs = (g.nodes.get(b) or {}).get("file_set")
            if fs:
                keep.add(fs)
    if any(n.get("kind") == "donor" for n in g.nodes.values()):
        keep |= {k for k, n in g.nodes.items() if n.get("kind") == "donor"}
    return keep


def collapse(g: Graph) -> "tuple[dict[str, dict], list[tuple[str, str, str]]]":
    """Group files into one node per (file set, product) so the picture stays readable."""
    group_of = {}
    groups: "dict[str, dict]" = {}
    keep = _drawable(g)
    n_sib = sum(1 for n in g.nodes.values() if n.get("relation") == "sibling" and n.get("kind") in ("stub", "fileset"))
    for k, n in g.nodes.items():
        if k not in keep:
            continue
        if n.get("kind") == "file":
            key = f"{n.get('file_set') or 'nofs'}::{n.get('product')}"
            group_of[k] = key
            gr = groups.setdefault(key, {"kind": "filegroup", "product": n.get("product"), "file_set": n.get("file_set"),
                                         "n": 0, "gb": 0.0, "formats": Counter()})
            gr["n"] += 1
            gr["gb"] += _gb(n.get("file_size"))
            gr["formats"][n.get("file_format")] += 1
        elif n.get("kind") == "qc":
            continue
        else:
            group_of[k] = k
            groups[k] = n
    edges = set()
    for a, b, lab in g.edges:
        if lab in ("derived_from", "quality_metrics", "seqspecs", "multiplexed_samples", "file_sets"):
            continue
        ga, gb = group_of.get(a), group_of.get(b)
        if ga and gb and ga != gb:
            edges.add((ga, gb, lab))
    if n_sib:
        groups["__siblings__"] = {"kind": "sibling_summary", "n": n_sib}
        root_samples = [k for k in groups if (g.nodes.get(k) or {}).get("kind") == "sample"]
        for sid in root_samples[:1]:
            edges.add((sid, "__siblings__", "file_sets"))
    return groups, sorted(edges)


def mermaid(g: Graph) -> str:
    groups, edges = collapse(g)
    donors = [k for k, n in groups.items() if n.get("kind") == "donor"]
    if len(donors) > 1:
        for k in donors:
            groups.pop(k)
        groups["__donors__"] = {"kind": "donor", "label": f"{len(donors)} donors"}
        edges = sorted({(a if a not in donors else "__donors__", b if b not in donors else "__donors__", lab)
                        for a, b, lab in edges})
    ids = {k: f"n{i}" for i, k in enumerate(groups)}
    out = ["flowchart LR"]
    for k, n in groups.items():
        if n.get("kind") == "sibling_summary":
            lab = f"{n['n']} sibling datasets (not walked)"
        elif n.get("kind") == "donor" and n.get("label"):
            lab = n["label"]
        elif n.get("kind") == "filegroup":
            lab = f"{n['n']} x {PRODUCT_LABEL.get(n['product'], n['product'])}<br/>{n['gb']:.1f} GB"
        else:
            lab = _short(n).replace("\\n", "<br/>")
        out.append(f'  {ids[k]}["{lab}"]')
    for a, b, lab in edges:
        if a in ids and b in ids:
            out.append(f"  {ids[a]} -->|{lab}| {ids[b]}")
    return "\n".join(out)


COLORS = {"MeasurementSet": "#f4a3a3", "AuxiliarySet": "#f4c9a3", "AnalysisSet": "#f2b6e8", "principal": "#9ef0e6",
          "CuratedSet": "#bfe8b0", "sample": "#f5b97a", "donor": "#cdb89a", "file": "#e3e3e8", "raw_reads": "#d6d6de",
          "blocked": "#ffffff", "document": "#f7d6d6", "PredictionSet": "#c6d8ff"}


def draw(g: Graph, path: Path) -> "Optional[Path]":
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyBboxPatch
    except Exception:
        return None
    groups, edges = collapse(g)

    def column(k: str, n: dict) -> int:
        kind = n.get("kind")
        if kind == "donor":
            return 0
        if kind == "sample":
            return 1
        if kind == "filegroup":
            fs = groups.get(n.get("file_set") or "", {})
            base = column(n.get("file_set"), fs) if fs and fs is not n else 3
            return base + (0 if n.get("product") in ("raw_reads", "seqspec", "onlist", "demultiplexing") else 1) \
                if fs.get("type") in ("MeasurementSet", "AuxiliarySet", "CuratedSet") else base + 1
        if kind == "document":
            return 7
        if kind == "sibling_summary":
            return 2
        t = n.get("type")
        if t in ("CuratedSet",):
            return 2
        if t in ("MeasurementSet", "AuxiliarySet"):
            return 3
        if t == "AnalysisSet":
            return 6 if n.get("file_set_type") == "principal analysis" else 4
        if kind == "blocked":
            return 6 if t == "AnalysisSet" else 3
        return 5
    cols: "dict[int, list[str]]" = defaultdict(list)
    donors = [k for k, n in groups.items() if n.get("kind") == "donor"]
    for k, n in groups.items():
        if n.get("kind") == "donor":
            continue
        cols[column(k, n)].append(k)
    if donors:
        cols[0] = ["__donors__"]
        groups["__donors__"] = {"kind": "donor", "label": f"{len(donors)} donors"}
        edges = [(a if a not in donors else "__donors__", b if b not in donors else "__donors__", lab) for a, b, lab in edges]
    pos = {}
    maxh = max((len(v) for v in cols.values()), default=1)
    for c, ks in cols.items():
        for i, k in enumerate(sorted(ks)):
            pos[k] = (c * 2.6, -(i + 0.5) * (maxh / max(1, len(ks))) * 1.3)
    fig, ax = plt.subplots(figsize=(max(10, 2.6 * (max(cols) + 1)), max(4, 1.3 * maxh + 1)))
    for a, b, lab in set(edges):
        if a in pos and b in pos and a != b:
            (x1, y1), (x2, y2) = pos[a], pos[b]
            ax.annotate("", xy=(x2 - 1.0, y2), xytext=(x1 + 1.0, y1),
                        arrowprops=dict(arrowstyle="-|>", color="#555555", lw=0.8, shrinkA=2, shrinkB=2))
            ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.08, lab, fontsize=6, color="#333333", ha="center")
    for k, (x, y) in pos.items():
        n = groups[k]
        if n.get("kind") == "filegroup":
            label = f"{n['n']} x {PRODUCT_LABEL.get(n['product'], n['product'])}\n{n['gb']:.1f} GB " + \
                    "/".join(str(f) for f in n["formats"] if f)
            col = COLORS["raw_reads"] if n["product"] == "raw_reads" else COLORS["file"]
        elif n.get("kind") == "sibling_summary":
            label, col = f"{n['n']} sibling datasets\nfrom the same pool\n(not walked)", "#f4f4f4"
        elif n.get("kind") == "donor":
            label, col = n.get("label") or "donor", COLORS["donor"]
        elif n.get("kind") == "sample":
            label, col = f"{n['accession']}\n{n['type']}", COLORS["sample"]
        elif n.get("kind") == "document":
            label, col = f"Document\n{n.get('document_type')}", COLORS["document"]
        elif n.get("kind") == "blocked":
            label, col = f"{n['accession']}\n{n['type']}\nHTTP {n.get('http_status')} (not visible)", COLORS["blocked"]
        else:
            t = n.get("type")
            label = f"{n['accession']}\n{t}\n{n.get('file_set_type') or ''}"
            if n.get("uniform_pipeline_status"):
                label += f"\nuniform: {n['uniform_pipeline_status']}"
            col = COLORS["principal"] if n.get("file_set_type") == "principal analysis" else COLORS.get(t, "#eeeeee")
        ls = "--" if n.get("kind") == "blocked" else "-"
        ax.add_patch(FancyBboxPatch((x - 1.0, y - 0.42), 2.0, 0.84, boxstyle="round,pad=0.02", fc=col, ec="#333333",
                                    lw=1.0, linestyle=ls))
        ax.text(x, y, label, ha="center", va="center", fontsize=6.5)
    ax.set_xlim(-1.3, max(p[0] for p in pos.values()) + 1.3 if pos else 1)
    ys = [p[1] for p in pos.values()] or [0]
    ax.set_ylim(min(ys) - 0.8, max(ys) + 0.8)
    ax.axis("off")
    ax.set_title(f"IGVF Portal lineage from {_acc(g.root)}", loc="left", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Offline fixture (the shape of IGVFDS9875NBZW on 2026-09-23) for tests
# ---------------------------------------------------------------------------

def fixture_portal() -> "tuple[Fetch, dict]":
    ms, atac, aux = "/measurement-sets/TSTDS0001RNA/", "/measurement-sets/TSTDS0002ATC/", "/auxiliary-sets/TSTDS0003AUX/"
    uni, lab, prin = "/analysis-sets/TSTDS0004UNI/", "/analysis-sets/TSTDS0005LAB/", "/analysis-sets/TSTDS0006PRN/"
    samp, bmap, cur = "/multiplexed-samples/TSTSM0001MUX/", "/tabular-files/TSTFI0020BMP/", "/curated-sets/TSTDS0007BAR/"

    def reads(i, fs, size=10_000_000_000, ctrl=True):
        return {"@id": f"/sequence-files/TSTFI{i:04d}RDS/", "accession": f"TSTFI{i:04d}RDS", "content_type": "reads",
                "file_format": "fastq", "file_size": size, "controlled_access": ctrl, "status": "released",
                "seqspecs": [f"/configuration-files/TSTFI{i:04d}SQS/"], "href": f"/sequence-files/TSTFI{i:04d}RDS/@@download/x"}

    def seqspec(i):
        return {"@id": f"/configuration-files/TSTFI{i:04d}SQS/", "accession": f"TSTFI{i:04d}SQS", "content_type": "seqspec",
                "file_format": "yaml", "file_size": 1400, "status": "released"}
    objs = {
        ms: {"@id": ms, "@type": ["MeasurementSet", "FileSet", "Item"], "accession": "TSTDS0001RNA", "status": "released",
             "file_set_type": "experimental data", "controlled_access": True, "summary": "snRNA-seq (10x multiome with MULTI-seq)",
             "assay_term": {"term_name": "single-nucleus RNA sequencing assay"},
             "files": [reads(1, ms), reads(2, ms), reads(3, ms, 1_000_000_000), seqspec(1)],
             "input_for": [{"@id": uni, "accession": "TSTDS0004UNI"}, {"@id": lab, "accession": "TSTDS0005LAB"}],
             "auxiliary_sets": [{"@id": aux, "accession": "TSTDS0003AUX"}],
             "related_measurement_sets": [{"series_type": "multiome", "measurement_sets": [{"@id": atac}]}],
             "samples": [{"@id": samp}], "onlist_files": ["/tabular-files/TSTFI0030ONL/"],
             "donors": [{"@id": "/human-donors/TSTDO0001AAA/", "accession": "TSTDO0001AAA", "sex": "male"},
                        {"@id": "/human-donors/TSTDO0002BBB/", "accession": "TSTDO0002BBB", "sex": "female"}]},
        atac: {"@id": atac, "@type": ["MeasurementSet", "FileSet", "Item"], "accession": "TSTDS0002ATC", "status": "released",
               "file_set_type": "experimental data", "assay_term": {"term_name": "single-nucleus ATAC-seq"},
               "files": [reads(11, atac, 20_000_000_000), reads(12, atac, 20_000_000_000)],
               "input_for": [{"@id": uni}, {"@id": lab}],
               "related_measurement_sets": [{"series_type": "multiome", "measurement_sets": [{"@id": ms}]}],
               "samples": [{"@id": samp}]},
        aux: {"@id": aux, "@type": ["AuxiliarySet", "FileSet", "Item"], "accession": "TSTDS0003AUX", "status": "released",
              "file_set_type": "lipid-conjugated oligo sequencing", "files": [reads(21, aux, 500_000_000)],
              "input_for": [lab], "measurement_sets": [ms]},
        uni: {"@id": uni, "@type": ["AnalysisSet", "FileSet", "Item"], "accession": "TSTDS0004UNI", "status": "released",
              "file_set_type": "intermediate analysis", "uniform_pipeline_status": "completed",
              "workflows": [{"name": "IGVF Single Cell Uniform Pipeline", "accession": "TSTWF1", "uniform_pipeline": True,
                             "workflow_version": "v1.1.0"}],
              "input_file_sets": [atac, ms],
              "files": [{"@id": "/matrix-files/TSTFI0101H5A/", "accession": "TSTFI0101H5A", "content_type": "cell by gene matrix",
                         "file_format": "h5ad", "file_size": 3_700_000_000, "status": "released",
                         "derived_from": ["/sequence-files/TSTFI0001RDS/"]},
                        {"@id": "/matrix-files/TSTFI0102TAR/", "accession": "TSTFI0102TAR",
                         "content_type": "kallisto cell by gene matrix", "file_format": "tar", "file_size": 13_000_000_000,
                         "status": "released"},
                        {"@id": "/tabular-files/TSTFI0103FRG/", "accession": "TSTFI0103FRG", "content_type": "fragments",
                         "file_format": "bed", "file_size": 5_600_000_000, "controlled_access": False, "status": "released"},
                        {"@id": "/alignment-files/TSTFI0104BAM/", "accession": "TSTFI0104BAM", "content_type": "alignments",
                         "file_format": "bam", "file_size": 38_000_000_000, "controlled_access": True, "status": "released"},
                        {"@id": "/index-files/TSTFI0105BAI/", "accession": "TSTFI0105BAI", "content_type": "index",
                         "file_format": "bai", "file_size": 7_000_000, "status": "released"}]},
        lab: {"@id": lab, "@type": ["AnalysisSet", "FileSet", "Item"], "accession": "TSTDS0005LAB", "status": "released",
              "file_set_type": "intermediate analysis", "workflows": [{"name": "Lab Cell Ranger", "uniform_pipeline": False}],
              "input_file_sets": [aux, ms, atac], "input_for": [prin, "/prediction-sets/TSTDS0008PRD/"],
              "files": [{"@id": "/matrix-files/TSTFI0201HDF/", "accession": "TSTFI0201HDF", "content_type": "cell by gene matrix",
                         "file_format": "hdf5", "file_size": 225_000_000, "status": "released"},
                        {"@id": "/tabular-files/TSTFI0202HSH/", "accession": "TSTFI0202HSH", "content_type": "cell hashing barcodes",
                         "file_format": "tsv", "file_size": 540_000, "controlled_access": False, "status": "released"},
                        {"@id": "/tabular-files/TSTFI0203PKS/", "accession": "TSTFI0203PKS", "content_type": "peaks",
                         "file_format": "bed", "file_size": 2_800_000, "controlled_access": True, "status": "released"}]},
        samp: {"@id": samp, "@type": ["MultiplexedSample", "Sample", "Item"], "accession": "TSTSM0001MUX", "status": "released",
               "summary": "multiplexed iPSC, 2 donors", "barcode_map": bmap, "file_sets": [ms, atac, aux, uni, lab, prin],
               "multiplexed_samples": ["/in-vitro-systems/A/", "/in-vitro-systems/B/"], "donors": ["x", "y"]},
        bmap: {"@id": bmap, "@type": ["TabularFile", "File", "Item"], "accession": "TSTFI0020BMP",
               "content_type": "barcode to sample mapping", "file_format": "tsv", "file_size": 233, "controlled_access": False,
               "status": "released", "file_set": cur},
        cur: {"@id": cur, "@type": ["CuratedSet", "FileSet", "Item"], "accession": "TSTDS0007BAR", "file_set_type": "barcodes",
              "status": "released", "files": [{"@id": bmap, "accession": "TSTFI0020BMP", "content_type": "barcode to sample mapping",
                                              "file_format": "tsv", "file_size": 233}]},
        "/tabular-files/TSTFI0030ONL/": {"@id": "/tabular-files/TSTFI0030ONL/", "@type": ["TabularFile", "File", "Item"],
                                          "accession": "TSTFI0030ONL", "content_type": "barcode onlist", "file_format": "tsv",
                                          "file_size": 50_000_000, "status": "released", "file_set": cur},
        "/matrix-files/TSTFI0101H5A/": {"@id": "/matrix-files/TSTFI0101H5A/", "@type": ["MatrixFile", "File", "Item"],
                                         "accession": "TSTFI0101H5A", "content_type": "cell by gene matrix", "file_format": "h5ad",
                                         "file_size": 3_700_000_000, "file_set": uni, "quality_metrics": ["/qm/rna1/"],
                                         "transcriptome_annotation": "GENCODE 43", "derived_from": ["/sequence-files/TSTFI0001RDS/"],
                                         "analysis_step_version": {"software_versions": [{"summary": "kallisto-bustools v1.1.0"}]}},
        "/tabular-files/TSTFI0103FRG/": {"@id": "/tabular-files/TSTFI0103FRG/", "@type": ["TabularFile", "File", "Item"],
                                          "accession": "TSTFI0103FRG", "content_type": "fragments", "file_format": "bed",
                                          "file_size": 5_600_000_000, "controlled_access": False, "file_set": uni},
        "/prediction-sets/TSTDS0008PRD/": {
            "@id": "/prediction-sets/TSTDS0008PRD/", "@type": ["PredictionSet", "FileSet", "Item"], "accession": "TSTDS0008PRD",
            "status": "released", "file_set_type": "element-gene links", "scope": "genome-wide",
            "summary": "element-gene links prediction using scE2G v1.2", "input_file_sets": ["/model-sets/TSTDS0009MOD/", lab],
            "large_scale_loci_list": "/reference-files/TSTFI0401LOC/", "samples": [samp],
            "files": [{"@id": "/tabular-files/TSTFI0301E2G/", "accession": "TSTFI0301E2G", "content_type": "thresholded element gene links",
                       "file_format": "tsv", "file_size": 120_000_000, "controlled_access": False, "status": "released",
                       "derived_from": ["/model-files/TSTFI0501MDL/"]}]},
        "/tabular-files/TSTFI0301E2G/": {"@id": "/tabular-files/TSTFI0301E2G/", "@type": ["TabularFile", "File", "Item"],
                                          "accession": "TSTFI0301E2G", "content_type": "thresholded element gene links",
                                          "file_format": "tsv", "file_size": 120_000_000, "file_set": "/prediction-sets/TSTDS0008PRD/",
                                          "derived_from": ["/model-files/TSTFI0501MDL/"]},
        "/model-sets/TSTDS0009MOD/": {
            "@id": "/model-sets/TSTDS0009MOD/", "@type": ["ModelSet", "FileSet", "Item"], "accession": "TSTDS0009MOD",
            "status": "released", "file_set_type": "element-gene links", "model_name": "scE2G", "model_version": "v1.2",
            "input_file_sets": ["/analysis-sets/TSTDS0010TRN/"], "input_for": ["/prediction-sets/TSTDS0008PRD/"],
            "files": [{"@id": "/model-files/TSTFI0501MDL/", "accession": "TSTFI0501MDL", "content_type": "model weights",
                       "file_format": "pkl", "file_size": 10_000, "derived_from": ["/tabular-files/TSTFI0601TRN/"]}]},
        "/model-files/TSTFI0501MDL/": {"@id": "/model-files/TSTFI0501MDL/", "@type": ["ModelFile", "File", "Item"],
                                        "accession": "TSTFI0501MDL", "content_type": "model weights", "file_format": "pkl",
                                        "file_size": 10_000, "file_set": "/model-sets/TSTDS0009MOD/",
                                        "derived_from": ["/tabular-files/TSTFI0601TRN/"]},
        "/tabular-files/TSTFI0601TRN/": {"@id": "/tabular-files/TSTFI0601TRN/", "@type": ["TabularFile", "File", "Item"],
                                          "accession": "TSTFI0601TRN", "content_type": "CRISPR screen training data",
                                          "file_format": "tsv", "file_size": 2_000_000, "file_set": "/analysis-sets/TSTDS0010TRN/"},
        "/analysis-sets/TSTDS0010TRN/": {"@id": "/analysis-sets/TSTDS0010TRN/", "@type": ["AnalysisSet", "FileSet", "Item"],
                                          "accession": "TSTDS0010TRN", "file_set_type": "principal analysis", "status": "released",
                                          "input_for": ["/model-sets/TSTDS0009MOD/"] + [f"/model-sets/TSTDS9{i:03d}HUB/" for i in range(40)],
                                          "files": [{"@id": "/tabular-files/TSTFI0601TRN/", "accession": "TSTFI0601TRN",
                                                     "content_type": "CRISPR screen training data", "file_format": "tsv"}]},
        "/reference-files/TSTFI0401LOC/": {"@id": "/reference-files/TSTFI0401LOC/", "@type": ["ReferenceFile", "File", "Item"],
                                            "accession": "TSTFI0401LOC", "content_type": "loci", "file_format": "bed",
                                            "file_size": 50_000, "file_set": cur},
        "/qm/rna1/": {"@id": "/qm/rna1/", "@type": ["SingleCellRnaSeqQualityMetric", "QualityMetric", "Item"],
                      "description": "RNAseq Kallisto Bustools QC metric", "n_reads": 767652391, "p_pseudoaligned": 74.9,
                      "percentage_reads_on_onlist": 93.37, "total_umis": 16775972,
                      "quality_metric_of": ["/matrix-files/TSTFI0101H5A/", "/matrix-files/TSTFI0102TAR/"],
                      "rnaseq_kb_info": {"href": "@@download/rnaseq_kb_info/kb_info.json"}},
        "/qm/atac1/": {"@id": "/qm/atac1/", "@type": ["SingleCellAtacSeqQualityMetric", "QualityMetric", "Item"],
                       "description": "ATACseq chromap fragments QC metric", "pct_duplicates": 29.3, "n_reads": 1839833072,
                       "quality_metric_of": ["/tabular-files/TSTFI0103FRG/"]},
    }
    calls: "list[str]" = []

    def fetch(url: str):
        calls.append(url)
        path = url.split("?")[0]
        q = urllib.parse.parse_qs(url.split("?", 1)[1]) if "?" in url else {}
        if path == "/search/":
            if q.get("type") == ["QualityMetric"]:
                of = q.get("quality_metric_of", [""])[0]
                hits = [o for o in objs.values() if of in (o.get("quality_metric_of") or [])]
                return (200, {"@graph": hits}) if hits else (404, {"@graph": []})
            if q.get("type") == ["FileSet"]:
                tgt = q.get("input_file_sets.@id", [""])[0]
                hits = [o for o in objs.values() if tgt in (o.get("input_file_sets") or [])]
                return (200, {"@graph": hits}) if hits else (404, {"@graph": []})
            return 404, None
        if path == prin:
            return 403, None
        if path.strip("/") in [o.get("accession") for o in objs.values()]:
            for o in objs.values():
                if o.get("accession") == path.strip("/"):
                    return 200, o
        if path in objs:
            return 200, objs[path]
        # a file with no own object: return a minimal object
        if path.split("/")[1].endswith("-files"):
            for o in objs.values():
                for f in o.get("files") or []:
                    if isinstance(f, dict) and f.get("@id") == path:
                        return 200, {**f, "@type": ["File"], "file_set": o["@id"]}
        return 404, None
    return fetch, {"calls": calls, "ids": {"ms": ms, "atac": atac, "aux": aux, "uni": uni, "lab": lab, "prin": prin,
                                         "samp": samp, "bmap": bmap, "cur": cur}}


def selftest(checks: "list", tmp: Path) -> None:
    def check(cond, msg):
        checks.append((bool(cond), msg))
        print(("  ok    " if cond else "  FAIL  ") + msg)
    fetch, info = fixture_portal()
    ids = info["ids"]
    g = walk("TSTDS0001RNA", fetch=fetch)
    kinds = Counter(n["kind"] for n in g.nodes.values())
    check(all(i in g.nodes for i in (ids["ms"], ids["atac"], ids["aux"], ids["uni"], ids["lab"], ids["samp"], ids["cur"])),
          "walk reaches the multiome partner, auxiliary set, both analysis sets, the multiplexed sample and the curated barcode set")
    check(ids["prin"] in g.blocked and g.blocked[ids["prin"]] == 403, "the unreleased principal analysis is recorded as blocked (403), not dropped")
    check((ids["lab"], ids["prin"], "input_for") in g.edges and (ids["samp"], ids["bmap"], "barcode_map") in g.edges,
          "edges input_for (intermediate -> principal) and barcode_map are recorded")
    check(kinds["qc"] == 2, f"both QC metric objects fetched ({kinds['qc']})")
    plan = build_plan(g)
    sh = {s["product"]: s for s in plan["start_here"]}
    check(plan["verdict"] == "processed", "verdict: processed results exist")
    check(sh["rna_matrix"]["accession"] == "TSTFI0101H5A", "RNA matrix: the uniform-pipeline h5ad beats the lab hdf5 and the 13 GB tar")
    check(sh["fragments"]["accession"] == "TSTFI0103FRG" and sh["fragments"]["access"] == "public", "fragments: public uniform file")
    check(sh["demultiplexing"]["accession"] in ("TSTFI0202HSH", "TSTFI0020BMP"), "demultiplexing product found (hashing barcodes / barcode map)")
    check(sh["rna_matrix"]["qc_headline"].get("p_pseudoaligned") == 74.9, "portal QC attached to the chosen matrix")
    check("alignments" in sh and "controlled" in sh["alignments"]["access"], "alignments listed with controlled access")
    check(abs(plan["raw_gb"] - 61.5) < 0.01, f"raw reads reachable: {plan['raw_gb']} GB across RNA, ATAC and MULTI-seq")
    check(plan["processed_gb_to_download"] < 10, f"processed start set is {plan['processed_gb_to_download']} GB")
    check(any(c["product"] == "onlist" for c in plan["configs"]) and any(c["product"] == "seqspec" for c in plan["configs"]),
          "seqspec and onlist configuration files listed")
    plan2 = build_plan(g, wants=["fragments"])
    check([s["product"] for s in plan2["start_here"]] == ["fragments"], "--want filters the start-here list")
    fig = draw(g, tmp / "lineage.png")
    rep = write_report(plan, g, tmp / "report.md", fig)
    txt = rep.read_text()
    check("Start here" in txt and "TSTFI0101H5A" in txt and "QC the Portal already computed" in txt and "HTTP 403" in txt,
          "report has start-here, portal QC and the not-visible section")
    check("flowchart LR" in txt, "report embeds a mermaid lineage graph")
    check(fig is None or fig.stat().st_size > 5000, "lineage figure drawn")
    # walking from a processed file climbs back to its sets
    g2 = walk("/matrix-files/TSTFI0101H5A/", fetch=fetch, max_depth=3)
    check(ids["uni"] in g2.nodes and ids["ms"] in g2.nodes, "walk from a file climbs to its analysis set and the raw measurement set")
    pm = {m["accession"]: m for m in plan["predictions_and_models"]}
    check("TSTDS0008PRD" in pm and "TSTDS0009MOD" in pm and "TSTDS0005LAB" in pm["TSTDS0008PRD"]["inputs"],
          "from raw reads the walk reaches the prediction set built on its analysis set, and the prediction set's model")
    check(sh.get("predictions", {}).get("accession") == "TSTFI0301E2G", "prediction output file typed as predictions by its set")
    g4 = walk("TSTDS0008PRD", fetch=fetch, max_depth=4, fanout=10)
    check("/tabular-files/TSTFI0601TRN/" in g4.nodes and "/analysis-sets/TSTDS0010TRN/" in g4.nodes,
          "from a prediction set: model file -> derived_from -> training data -> experimental analysis set")
    check(("/prediction-sets/TSTDS0008PRD/", "/reference-files/TSTFI0401LOC/", "large_scale_loci_list") in g4.edges,
          "large_scale_loci_list edge recorded")
    check(not any("HUB" in k for k in g4.nodes), "an upstream training set's other downstream models are not walked")
    g6 = walk("/analysis-sets/TSTDS0010TRN/", fetch=fetch, max_depth=4, fanout=50)
    pr = g6.nodes.get("/prediction-sets/TSTDS0008PRD/") or {}
    check(pr.get("relation") == "model_use", "predictions made by a model trained on this data elsewhere are tagged model_use")
    plan6 = build_plan(g6)
    check(not any(x["product"] == "predictions" for x in plan6["start_here"]),
          "model_use predictions are not offered as results of the training data")
    g5 = walk("/analysis-sets/TSTDS0010TRN/", fetch=fetch, max_depth=3, fanout=10)
    trn = g5.nodes["/analysis-sets/TSTDS0010TRN/"]
    n_hub = sum(1 for k in g5.nodes if "HUB" in k)
    check(trn.get("truncated", {}).get("input_for") == 41 and n_hub <= 10,
          f"hub fan-out capped when walking down from a shared set (41 recorded, {n_hub} walked)")
    check(not any(n.get("relation") == "sibling" and n.get("kind") == "fileset" for n in g.nodes.values()),
          "sibling datasets of the pool are recorded as stubs, not expanded")
    g3 = walk("TSTSM0001MUX", fetch=fetch, max_depth=2)
    check(ids["uni"] in g3.nodes and ids["bmap"] in g3.nodes, "walk from a multiplexed sample reaches its file sets and barcode map")
