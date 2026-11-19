"""Strict, relocatable feature loading and shared metrics for B/E clustering.

Native feature rows are joined by explicit row number and stable sample ID.
Family labels are evaluation sidecars, never input features or ordering keys.
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATASETS = {"B": "Dataset_B_Divergent", "E": "Dataset_E_MIBiG_BS1"}
ROUTES = {"esm2": "esm", "mac": "mac", "dna_full": "dna"}
ARTIFACTS = ("embeddings.npy", "bgc_ids.tsv", "metadata.json", "failures.tsv")


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")


def read_tsv(path, required=()):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if (not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames))
                or not set(required) <= set(reader.fieldnames)):
            raise ValueError(f"Missing/duplicate TSV columns: {path}")
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError(f"Malformed TSV row: {path}")
        return rows


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    if not fields:
        raise ValueError(f"Explicit columns required for empty TSV: {path}")
    with Path(path).open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def add_common_arguments(parser, *, default_output):
    parser.add_argument("--features-root", type=Path, default=ROOT / "data/processed/features_be_v1",
                        help="Root containing Dataset_*/{esm2,mac,dna_full}/ native exports")
    parser.add_argument("--cohort", type=Path, default=ROOT / "data/processed/pilot_be_v1/cohort.tsv")
    parser.add_argument("--datasets", nargs="+", choices=tuple(DATASETS), default=["B", "E"])
    parser.add_argument("--routes", nargs="+", choices=tuple(ROUTES), default=list(ROUTES))
    parser.add_argument("--output", type=Path, default=ROOT / default_output,
                        help="New output directory; existing results are never overwritten")
    parser.add_argument("--no-plots", action="store_true", help="Disable optional diagnostic figures")


def load_cohort(path, dataset_names):
    from pilot_features import validate_rows

    path = Path(path).expanduser().resolve()
    rows = read_tsv(path, ("dataset", "sample_id", "original_bgc_id", "archive_relpath",
                           "source_sha256", "gcf_raw", "tsv_row_ids_json"))
    metadata_path = path.parent / "metadata.json"
    if metadata_path.is_file():
        metadata = read_json(metadata_path)
        expected_hash = metadata.get("artifact_sha256", {}).get(path.name)
        if expected_hash is None or expected_hash != file_sha256(path):
            raise ValueError("Cohort checksum differs from its preparation metadata")
    rows = [row for row in rows if row["dataset"] in dataset_names]
    validate_rows(rows, dataset_names)
    if any(not row["gcf_raw"].strip() for row in rows):
        raise ValueError("Empty GCF label in cohort")
    return sorted(rows, key=lambda row: (row["dataset"], row["sample_id"]))


def cosine_distances(vectors):
    values = np.asarray(vectors, dtype=np.float64)
    if values.ndim != 2 or min(values.shape) < 1 or not np.isfinite(values).all():
        raise ValueError("Expected a finite nonempty vector matrix")
    norms = np.linalg.norm(values, axis=1)
    if np.any(norms <= 0) or not np.isfinite(norms).all():
        raise ValueError("Cosine distance is undefined for zero/invalid vectors")
    unit = values / norms[:, None]
    similarities = np.clip(unit @ unit.T, -1.0, 1.0)
    distances = np.clip(1.0 - similarities, 0.0, 2.0)
    distances = (distances + distances.T) / 2
    np.fill_diagonal(distances, 0.0)
    return unit, distances


def load_feature_bundle(directory, rows, dataset, route):
    """Load one complete export; accept TSV row reordering, never intersect IDs."""
    directory = Path(directory).expanduser().resolve()
    if route not in ROUTES or not rows:
        raise ValueError("Unknown route or empty cohort")
    rows = sorted(rows, key=lambda row: row["sample_id"])
    expected_ids = [row["sample_id"] for row in rows]
    if len(set(expected_ids)) != len(expected_ids):
        raise ValueError("Duplicate cohort IDs")
    for name in ARTIFACTS:
        if not (directory / name).is_file():
            raise ValueError(f"Missing native feature artifact: {directory / name}")
    metadata = read_json(directory / "metadata.json")
    if (metadata.get("format") != "bgc_vectors_v1" or metadata.get("status") != "complete"
            or metadata.get("route") != ROUTES[route] or metadata.get("excluded_count") != 0
            or metadata.get("upstream_excluded_count", 0) != 0 or metadata.get("error")):
        raise ValueError(f"Export is incomplete, excluded samples, or has the wrong route: {directory}")
    source = metadata.get("input", {})
    if source.get("input_count") != len(rows) or source.get("exclusions"):
        raise ValueError(f"Input coverage differs from the cohort: {directory}")
    if read_tsv(directory / "failures.tsv", ("bgc_id", "source", "stage", "error")):
        raise ValueError(f"Export has failure records: {directory}")
    id_rows = read_tsv(directory / "bgc_ids.tsv", ("row", "bgc_id", "token_count"))
    ids = [row["bgc_id"] for row in id_rows]
    if len(ids) != len(set(ids)) or set(ids) != set(expected_ids):
        missing, extra = set(expected_ids) - set(ids), set(ids) - set(expected_ids)
        raise ValueError(f"Feature IDs do not match full cohort: missing={sorted(missing)}, extra={sorted(extra)}")
    indices = [int(row["row"]) for row in id_rows]
    if sorted(indices) != list(range(len(rows))):
        raise ValueError("Feature row indices must be a unique, complete 0..N-1 range")
    if any(int(row["token_count"]) <= 0 for row in id_rows):
        raise ValueError("Invalid token count in feature IDs")
    array = np.load(directory / "embeddings.npy", allow_pickle=False)
    if (not isinstance(array, np.ndarray) or array.dtype != np.float32 or array.ndim != 2
            or array.shape[0] != len(rows) or array.shape[1] < 1 or not np.isfinite(array).all()):
        raise ValueError("Feature matrix must be finite FP32 [cohort_count, dimension]")
    if (metadata.get("count") != len(rows) or metadata.get("dimension") != array.shape[1]
            or metadata.get("dtype") != "float32" or metadata.get("normalization") != "l2"):
        raise ValueError("Feature matrix and metadata disagree")
    norms = np.linalg.norm(array.astype(np.float64), axis=1)
    if not np.allclose(norms, 1.0, rtol=1e-5, atol=1e-5):
        raise ValueError("Exported vectors are not L2 normalized")
    records = metadata.get("records", {})
    if set(records) != set(expected_ids):
        raise ValueError("Metadata record coverage differs from cohort")
    for row in rows:
        record = records[row["sample_id"]]
        if (record.get("bgc_id") != row["sample_id"]
                or record.get("source_sha256") != row["source_sha256"]):
            raise ValueError(f"Source ID/hash mismatch: {row['sample_id']}")
    if route == "dna_full":
        from pilot_feature_validation import _validate_dna
        _validate_dna(metadata, records, id_rows)
    index_by_id = {row["bgc_id"]: int(row["row"]) for row in id_rows}
    source_indices = [index_by_id[key] for key in expected_ids]
    unit, distances = cosine_distances(array[source_indices])
    cache_hash = (metadata.get("saved_cache", {}).get("sha256") if route == "esm2"
                  else source.get("cache_sha256") if route == "mac" else None)
    provenance = {
        "feature_directory": str(directory), "dataset": dataset, "route": route,
        "input_count": len(rows), "dimension": int(array.shape[1]), "coverage": 1.0,
        "artifact_sha256": {name: file_sha256(directory / name) for name in ARTIFACTS},
        "token_cache_sha256": cache_hash,
        "raw_norm_range": [float(norms.min()), float(norms.max())],
        "cosine_distance_quantiles": {
            str(q): float(np.quantile(distances[np.triu_indices(len(rows), 1)], q))
            for q in (0.0, 0.25, 0.5, 0.75, 1.0)
        } if len(rows) > 1 else {},
        "checkpoint": metadata.get("extractor", {}),
    }
    return {"dataset": dataset, "route": route, "ids": expected_ids, "rows": rows,
            "vectors": unit, "distances": distances, "labels": [r["gcf_raw"] for r in rows],
            "source_row_indices": source_indices, "provenance": provenance}


def load_inputs(args):
    if len(set(args.datasets)) != len(args.datasets) or len(set(args.routes)) != len(args.routes):
        raise ValueError("Repeated datasets/routes are not allowed")
    features = args.features_root.expanduser().resolve()
    dataset_names = [DATASETS[key] for key in args.datasets]
    rows = load_cohort(args.cohort, dataset_names)
    output = args.output.expanduser().resolve()
    if output.is_relative_to(features) or output.is_relative_to(args.cohort.expanduser().resolve().parent):
        raise ValueError("Clustering outputs must be separate from features and prepared inputs")
    state_path = features / "run.json"
    state = read_json(state_path) if state_path.is_file() else None
    bundles = []
    for dataset in dataset_names:
        subset = [row for row in rows if row["dataset"] == dataset]
        group = []
        for route in args.routes:
            bundle = load_feature_bundle(features / dataset / route, subset, dataset, route)
            if state is not None:
                step = state.get("steps", {}).get(dataset + "/" + route, {})
                expected = step.get("validation", {}).get("artifact_sha256")
                if step.get("status") != "complete" or expected != bundle["provenance"]["artifact_sha256"]:
                    raise ValueError(f"Export no longer matches saved extraction validation: {dataset}/{route}")
            bundle["provenance"]["extraction_run_verified"] = state is not None
            group.append(bundle)
        by_route = {bundle["route"]: bundle for bundle in group}
        if "esm2" in by_route and "mac" in by_route:
            left = by_route["esm2"]["provenance"]["token_cache_sha256"]
            right = by_route["mac"]["provenance"]["token_cache_sha256"]
            if not left or left != right:
                raise ValueError(f"ESM2 and MAC do not share a verified token cache: {dataset}")
        bundles.extend(group)
    return bundles


def canonical_labels(membership):
    values = np.asarray(membership)
    if values.ndim != 1 or not len(values):
        raise ValueError("Expected nonempty one-dimensional cluster labels")
    mapping = {}
    return np.array([mapping.setdefault(value, len(mapping)) for value in values.tolist()], dtype=np.int64)


def score_partition(labels_true, labels_pred, distances):
    from sklearn.metrics import (adjusted_mutual_info_score, adjusted_rand_score,
                                 homogeneity_completeness_v_measure, silhouette_score)

    labels_pred = np.asarray(labels_pred)
    count = len(labels_true)
    if labels_pred.shape != (count,) or not count or np.asarray(distances).shape != (count, count):
        raise ValueError("Partition, true labels and distances must cover the same samples")
    if not np.issubdtype(labels_pred.dtype, np.integer) or (labels_pred < 0).any():
        raise ValueError("These baselines must assign every sample; negative/noise labels are unsupported")
    _, sizes = np.unique(labels_pred, return_counts=True)
    h, c, v = homogeneity_completeness_v_measure(labels_true, labels_pred)
    silhouette = None
    reason = "undefined: requires 2 <= predicted clusters <= N-1"
    if 2 <= len(sizes) < count:
        silhouette = float(silhouette_score(distances, labels_pred, metric="precomputed"))
        reason = ""
    return {"sample_count": count, "n_clusters": len(sizes),
            "n_singleton_clusters": int(np.sum(sizes == 1)),
            "singleton_sample_fraction": float(np.sum(sizes == 1) / count),
            "coverage": 1.0, "true_gcf_count": len(set(labels_true)),
            "homogeneity": float(h), "completeness": float(c), "v_measure": float(v),
            "ari": float(adjusted_rand_score(labels_true, labels_pred)),
            "ami": float(adjusted_mutual_info_score(labels_true, labels_pred)),
            "silhouette_cosine": silhouette, "silhouette_reason": reason}


def create_output(path):
    output = Path(path).expanduser().resolve()
    if output.exists():
        raise ValueError(f"Output already exists; choose a new directory: {output}")
    output.mkdir(parents=True)
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".matplotlib"))
    write_json(output / "status.json", {"status": "running", "started_at_utc": datetime.now(timezone.utc).isoformat()})
    return output


def write_assignments(path, bundle, runs):
    def records():
        for run_id, membership in runs:
            if len(membership) != len(bundle["ids"]):
                raise ValueError("Incomplete cluster assignments")
            for index, label in enumerate(membership):
                row = bundle["rows"][index]
                yield {"run_id": run_id, "row": index, "sample_id": row["sample_id"],
                       "original_bgc_id": row["original_bgc_id"], "gcf_raw": row["gcf_raw"],
                       "cluster_id": int(label)}
    write_tsv(path, records(), ("run_id", "row", "sample_id", "original_bgc_id", "gcf_raw", "cluster_id"))


def source_provenance():
    sources = ("clustering_common.py", "cluster_average.py", "cluster_leiden.py",
               "pilot_features.py", "pilot_feature_validation.py")
    versions = {"python": sys.version.split()[0]}
    for name in ("numpy", "scipy", "scikit-learn", "igraph", "leidenalg", "matplotlib"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return {"versions": versions,
            "source_sha256": {name: file_sha256(Path(__file__).with_name(name)) for name in sources
                              if Path(__file__).with_name(name).is_file()}}


def finish_run(output, args, bundles, method, parameters, metrics_rows):
    write_tsv(output / "metrics.tsv", metrics_rows)
    for bundle in bundles:
        directory = output / bundle["dataset"] / bundle["route"]
        directory.mkdir(parents=True, exist_ok=True)
        write_tsv(directory / "samples.tsv", [
            {"row": index, "source_feature_row": bundle["source_row_indices"][index], **row}
            for index, row in enumerate(bundle["rows"])
        ])
        np.save(directory / "cosine_distances.npy", bundle["distances"], allow_pickle=False)
    provenance = source_provenance()
    artifact_hashes = {p.relative_to(output).as_posix(): file_sha256(p)
                       for p in sorted(output.rglob("*"))
                       if p.is_file() and not any(part.startswith(".") for part in p.relative_to(output).parts)
                       and p.name != "status.json"}
    write_json(output / "metadata.json", {
        "format": "bgc_clustering_baseline_v1", "status": "complete", "method": method,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "parameters": parameters, "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "clustering_space": "original embedding cosine; float64 normalization for cosine arithmetic only",
        "projection": None, "input_order": "ascending stable sample_id",
        "cohort": {"path": str(args.cohort.expanduser().resolve()), "sha256": file_sha256(args.cohort)},
        "inputs": [bundle["provenance"] for bundle in bundles],
        "label_policy": "raw Unicode GCF strings, dataset-local; evaluation only; singletons retained",
        "selection_policy": "all requested settings reported; no best test-label-selected configuration",
        "comparison_limits": "end-to-end embedding+clusterer baseline; no BiG-SCAPE run; E/MAC training overlap not resolved",
        "artifact_sha256": artifact_hashes, **provenance,
    })
    status = output / "status.json"
    status.write_text(json.dumps({"status": "complete", "finished_at_utc": datetime.now(timezone.utc).isoformat()}, indent=2) + "\n")
