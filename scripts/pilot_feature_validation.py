"""Validate native BGC_annotation exports against a fixed pilot manifest.

No torch checkpoint/cache deserialization is needed: cache provenance is checked
using file hashes. Validation never repairs, reorders or drops output rows.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

OUTPUT_FILES = ("embeddings.npy", "bgc_ids.tsv", "metadata.json", "failures.tsv")
NATIVE_ROUTES = {"esm2": "esm", "mac": "mac", "dna_full": "dna"}


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _read_tsv(path, required):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _require(reader.fieldnames is not None and set(required) <= set(reader.fieldnames),
                 f"{path.name}: missing required TSV columns {required}")
        _require(len(reader.fieldnames) == len(set(reader.fieldnames)),
                 f"{path.name}: duplicate TSV columns")
        rows = list(reader)
        _require(all(None not in row and all(value is not None for value in row.values())
                     for row in rows), f"{path.name}: malformed TSV row")
        return rows


def _positive_int(value):
    return type(value) is int and value > 0


def _validate_dna(metadata, records, id_rows):
    source = metadata["input"]
    encoder = metadata.get("extractor")
    _require(isinstance(encoder, dict), "DNA: missing extractor metadata")
    _require(metadata.get("sequence_mode") == "full"
             and encoder.get("sequence_mode") == "full", "DNA: expected full sequence mode")
    _require("length_filter" in source and source["length_filter"] is None,
             "DNA: a length filter is enabled or its policy is missing")
    _require(source.get("explicit_exclusions") == [], "DNA: explicit exclusions are not allowed")
    _require("max_tokens" in encoder and encoder["max_tokens"] is None
             and "max_content_tokens" in encoder and encoder["max_content_tokens"] is None,
             "DNA: full input token limits are not allowed")
    _require(encoder.get("config_length_enforced") is False
             and encoder.get("automatic_window_fallback") is False,
             "DNA: configuration length limit or automatic window fallback is enabled/unknown")
    _require(encoder.get("window_policy") == "full_sequence_single_forward",
             "DNA: expected a full sequence single forward policy")
    _require(metadata.get("window_batch_size") == 1,
             "DNA: full mode requires one sequence per forward")
    _require(metadata.get("length_excluded_count") == 0 and metadata.get("invalid_count") == 0
             and metadata.get("length_exclusions") == [], "DNA: excluded or invalid records")
    for item in id_rows:
        key = item["bgc_id"]
        record = records[key]
        count = record.get("token_count")
        bases = record.get("base_count")
        _require(record.get("sequence_mode") == "full"
                 and record.get("window_count") == 1 and record.get("forward_pass_count") == 1,
                 f"DNA {key}: expected exactly one full window and one forward")
        _require(_positive_int(count) and item["token_count"] == str(count)
                 and record.get("content_token_count") == count,
                 f"DNA {key}: inconsistent content token count")
        _require(_positive_int(bases) and record.get("sequence_length") == bases
                 and record.get("sequence_region") == [0, bases],
                 f"DNA {key}: incomplete whole-record sequence coverage")
        _require(record.get("forward_input_shapes") == [[1, count + 2]]
                 and record.get("total_input_tokens") == count + 2
                 and record.get("max_forward_tokens") == count + 2,
                 f"DNA {key}: forward shape does not cover all tokens")
        _require(record.get("window_content_token_counts") == [count]
                 and record.get("windows") == [{"token_start": 0, "token_end": count,
                                                "base_start": 0, "base_end": bases}],
                 f"DNA {key}: window does not cover the entire input")


def validate_export(export_dir, rows, route, *, cache_path=None):
    """Return hashes/counts after strict validation, or raise ``ValueError``.

    ``rows`` must be the complete manifest in ascending ``sample_id`` order,
    containing at least ``sample_id`` and ``source_sha256``. Native exports must
    use these IDs as their GenBank stems and include every manifest entry.
    """
    try:
        return _validate_export(Path(export_dir), rows, route, cache_path=cache_path)
    except (OSError, UnicodeError, EOFError, csv.Error, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot validate feature export {export_dir}: {exc}") from exc


def _validate_export(export_dir, rows, route, *, cache_path):
    import numpy as np

    _require(route in NATIVE_ROUTES, f"Unknown pilot route: {route!r}")
    _require(isinstance(rows, list) and bool(rows), "Expected a nonempty manifest row list")
    _require(all(isinstance(row, dict) and isinstance(row.get("sample_id"), str)
                 and row["sample_id"] and row["sample_id"].strip() == row["sample_id"]
                 and isinstance(row.get("source_sha256"), str)
                 and len(row["source_sha256"]) == 64
                 and all(char in "0123456789abcdef" for char in row["source_sha256"])
                 for row in rows), "Manifest rows require valid sample_id and source_sha256")
    expected_ids = [row["sample_id"] for row in rows]
    _require(len(set(expected_ids)) == len(expected_ids), "Manifest has duplicate sample IDs")
    _require(expected_ids == sorted(expected_ids), "Manifest sample IDs must be sorted")
    missing = [name for name in OUTPUT_FILES if not (export_dir / name).is_file()]
    _require(not missing, f"Missing export artifacts: {missing}")
    metadata = json.loads((export_dir / "metadata.json").read_text(encoding="utf-8"))
    _require(isinstance(metadata, dict), "metadata.json must contain an object")
    _require(metadata.get("format") == "bgc_vectors_v1", "Unsupported export format")
    _require(metadata.get("status") == "complete", "Export status is not complete")
    _require(metadata.get("route") == NATIVE_ROUTES[route], "Export route does not match requested route")
    _require(type(metadata.get("excluded_count")) is int and metadata["excluded_count"] == 0,
             "Export has exclusions or lacks excluded_count")
    upstream = metadata.get("upstream_excluded_count", 0)
    _require(type(upstream) is int and upstream == 0, "Export has upstream exclusions")
    _require(not metadata.get("error") and metadata.get("export_complete", True) is True,
             "Export contains failure metadata")
    source = metadata.get("input")
    _require(isinstance(source, dict), "Missing input provenance")
    _require(not source.get("exclusions") and not metadata.get("exclusions"),
             "Export provenance contains exclusions")
    _require(type(source.get("input_count")) is int and source["input_count"] == len(rows),
             "Input count does not match the full manifest")
    failures = _read_tsv(export_dir / "failures.tsv", ("bgc_id", "source", "stage", "error"))
    _require(not failures, "failures.tsv contains excluded/failed samples")
    id_rows = _read_tsv(export_dir / "bgc_ids.tsv", ("row", "bgc_id", "token_count"))
    ids = [row["bgc_id"] for row in id_rows]
    _require(len(set(ids)) == len(ids), "bgc_ids.tsv contains duplicate sample IDs")
    _require(ids == expected_ids, "bgc_ids.tsv IDs/order do not match the complete sorted manifest")
    for index, item in enumerate(id_rows):
        _require(item["row"] == str(index), "bgc_ids.tsv row indices do not match matrix order")
        _require(item["token_count"].isascii() and item["token_count"].isdecimal()
                 and int(item["token_count"]) > 0, f"Invalid token count for {item['bgc_id']}")
    array = np.load(export_dir / "embeddings.npy", allow_pickle=False)
    _require(isinstance(array, np.ndarray), "embeddings.npy must be a single array")
    _require(array.dtype == np.dtype("float32"), "embeddings.npy must use float32")
    _require(array.ndim == 2 and array.shape[0] == len(rows) and array.shape[1] > 0,
             "Embedding matrix shape does not match the manifest")
    _require(np.isfinite(array).all(), "Embeddings contain NaN or infinity")
    _require(np.allclose(np.linalg.norm(array.astype(np.float64), axis=1), 1.0,
                        atol=1e-5, rtol=1e-5), "Embeddings are not L2 normalized")
    _require(type(metadata.get("count")) is int and metadata["count"] == array.shape[0]
             and type(metadata.get("dimension")) is int and metadata["dimension"] == array.shape[1],
             "Metadata count/dimension does not match the matrix")
    _require(metadata.get("dtype") == "float32" and metadata.get("normalization") == "l2",
             "Metadata dtype/normalization does not match required vector format")
    records = metadata.get("records")
    _require(isinstance(records, dict) and set(records) == set(expected_ids),
             "Metadata records do not match the full manifest")
    for row in rows:
        key = row["sample_id"]
        record = records[key]
        _require(isinstance(record, dict) and record.get("source_sha256") == row["source_sha256"],
                 f"Source SHA256 does not match manifest for {key}")
        _require(record.get("bgc_id") == key, f"Record BGC ID does not match manifest for {key}")
    if route == "dna_full":
        _validate_dna(metadata, records, id_rows)
    summary = {"route": route, "count": int(array.shape[0]), "dimension": int(array.shape[1]),
               "artifact_sha256": {name: file_sha256(export_dir / name) for name in OUTPUT_FILES}}
    if cache_path is not None:
        _require(route in {"esm2", "mac"}, "Token cache validation applies only to protein routes")
        cache_path = Path(cache_path)
        _require(cache_path.is_file(), f"Missing token cache: {cache_path}")
        cache_hash = file_sha256(cache_path)
        provenance = metadata.get("saved_cache") if route == "esm2" else source
        field = "sha256" if route == "esm2" else "cache_sha256"
        _require(isinstance(provenance, dict) and provenance.get(field) == cache_hash,
                 f"{route}: cache SHA256 does not match export provenance")
        summary["cache_sha256"] = cache_hash
    return summary
