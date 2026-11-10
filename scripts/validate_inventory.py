#!/usr/bin/env python3
"""Independently verify raw relocation and inventory evidence without writing files.

This validates an inventory, not a final evaluation cohort. Parser failures,
multi-record files, and unresolved label relationships are valid retained entries.
Only the Python standard library is required; no encoders are loaded.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys


SCHEMA_VERSION = "curated_bgc_inventory_v1"
TABLES = {"source_files.tsv", "samples.tsv", "records.tsv", "label_rows.tsv",
          "duplicate_groups.tsv", "duplicate_members.tsv", "issues.tsv", "dataset_summary.tsv"}
GBK_SUFFIXES = {".gbk", ".gb", ".genbank"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def compact(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def text_sha256(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def file_sha256(path):
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def safe_path(root, relative):
    rel = Path(relative)
    require(not rel.is_absolute() and ".." not in rel.parts, f"Unsafe relative path: {relative}")
    path = root / rel
    for part in (path, *path.parents):
        require(not part.is_symlink(), f"Symlink is unsupported: {part}")
    return path


def rows(path, expected_columns=None):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"Missing header: {path}")
        require(len(set(reader.fieldnames)) == len(reader.fieldnames), f"Duplicate columns: {path}")
        if expected_columns is not None:
            require(reader.fieldnames == expected_columns, f"Column schema differs: {path.name}")
        result = list(reader)
    require(all(None not in row and None not in row.values() for row in result), f"Malformed TSV rows: {path}")
    return result


def unique_index(items, key, label):
    index = {item[key]: item for item in items}
    require(all(item[key] for item in items) and len(index) == len(items), f"Empty/duplicate {label}")
    return index


def json_array(row, field):
    value = json.loads(row[field])
    require(isinstance(value, list), f"Expected JSON array in {field}")
    return value


def is_true(value):
    require(value in {"True", "False"}, f"Invalid Boolean: {value}")
    return value == "True"


def verify_raw(workspace, source_root, sources):
    indexed = unique_index(sources, "source_relpath", "source path")
    unique_index(sources, "archive_relpath", "archive path")
    actual = set()
    for path in source_root.rglob("*"):
        require(not path.is_symlink(), f"Symlink in raw data: {path}")
        if path.is_dir():
            continue
        require(path.is_file(), f"Special file in raw data: {path}")
        actual.add(path.relative_to(workspace).as_posix())
    require(actual == set(indexed), "Raw file set differs from source_files.tsv")
    canonical = []
    for row in sources:
        path = safe_path(workspace, row["source_relpath"])
        require(path == safe_path(source_root, row["archive_relpath"]), "Source/archive path mapping differs")
        require(path.stat().st_size == int(row["size_bytes"]), f"Raw file size changed: {path}")
        require(file_sha256(path) == row["source_sha256"], f"Raw SHA256 changed: {path}")
        canonical.append({**row, "size_bytes": int(row["size_bytes"])})

    report = safe_path(workspace, "data/processed/raw_relocation")
    migration_meta = json.loads((report / "metadata.json").read_text(encoding="utf-8"))
    plan = json.loads((report / "plan.json").read_text(encoding="utf-8"))
    require(migration_meta["status"] == "verified", "Raw migration is not verified")
    require(file_sha256(report / "files.tsv") == migration_meta["files_tsv_sha256"], "Migration files.tsv hash differs")
    require(file_sha256(report / "plan.json") == migration_meta["plan_json_sha256"], "Migration plan hash differs")
    migration = rows(report / "files.tsv")
    migrated = unique_index(migration, "raw_path", "migration destination path")
    require(set(migrated) == set(indexed), "Migration and inventory source file sets differ")
    for path, row in indexed.items():
        before = migrated[path]
        require((row["source_sha256"], int(row["size_bytes"])) ==
                (before["sha256"], int(before["size_bytes"])), f"Migration/source provenance differs: {path}")
        require(row["archive_relpath"] == before["relative_path"], f"Migration relative path differs: {path}")
    require(len(migration) == migration_meta["file_count"] == plan["file_count"], "Migration file counts differ")
    require(sum(int(row["size_bytes"]) for row in sources) == migration_meta["total_bytes"] == plan["total_bytes"],
            "Migration byte counts differ")
    directories = sorted(["."] + [p.relative_to(source_root).as_posix() for p in source_root.rglob("*") if p.is_dir()])
    require(directories == plan["directories"], "Raw directory set differs from migration plan")
    return indexed, canonical


def verify_inventory(workspace, inventory):
    workspace = Path(workspace).absolute()
    inventory = Path(inventory)
    if not inventory.is_absolute():
        inventory = workspace / inventory
    inventory = safe_path(workspace, inventory.relative_to(workspace).as_posix())
    require(inventory.is_dir(), f"Missing inventory: {inventory}")
    metadata = json.loads((inventory / "metadata.json").read_text(encoding="utf-8"))
    require(metadata["schema_version"] == SCHEMA_VERSION and metadata["status"] == "complete", "Unsupported/incomplete inventory")
    require(metadata["output_dir"] == inventory.relative_to(workspace).as_posix(), "Metadata output_dir differs from actual inventory path")
    expected_outputs = TABLES | {"summary.md"}
    require(set(metadata["output_sha256"]) == expected_outputs, "Unexpected/missing output checksum entries")
    require({p.name for p in inventory.iterdir()} == expected_outputs | {"metadata.json"}, "Inventory file set differs")
    for name, checksum in metadata["output_sha256"].items():
        path = safe_path(inventory, name)
        require(path.is_file() and file_sha256(path) == checksum, f"Output SHA256 differs: {name}")
    require(set(metadata["table_columns"]) == TABLES, "Table schemas differ")
    tables = {name: rows(inventory / name, metadata["table_columns"][name]) for name in TABLES}
    source_root = safe_path(workspace, metadata["source_root"])
    source_index, canonical_sources = verify_raw(workspace, source_root, tables["source_files.tsv"])
    require(text_sha256(compact(canonical_sources)) == metadata["source_snapshot_sha256"], "Source snapshot digest differs")
    require(metadata["source_verified_unchanged_after_inventory"] is True, "Source was not verified unchanged")
    policies = metadata["policies"]
    require(policies["length_filter"] is None and policies["repairs"] == "none" and policies["exclusions"] == "none",
            "Inventory filtering/repair policy differs")

    samples = tables["samples.tsv"]
    records = tables["records.tsv"]
    labels = tables["label_rows.tsv"]
    sample_index = unique_index(samples, "sample_id", "sample ID")
    record_index = unique_index(records, "record_uid", "record UID")
    label_index = unique_index(labels, "row_id", "label row ID")
    sample_sources = unique_index(samples, "source_relpath", "sample source path")
    gbk_paths = {path for path in source_index if Path(path).suffix.lower() in GBK_SUFFIXES}
    require(set(sample_sources) == gbk_paths, "GBK/sample source sets differ")
    records_by_sample = defaultdict(list)
    for record in records:
        sid = record["sample_id"]
        require(sid in sample_index, f"Record references unknown sample: {sid}")
        require(record["dataset"] == sample_index[sid]["dataset"], "Record/sample dataset differs")
        require(record["record_uid"] == f"{sid}:record:{record['record_index']}", "Record UID/index differs")
        require(int(record["length_bp"]) >= 0, "Negative record length")
        records_by_sample[sid].append(record)
    for sample in samples:
        sid = sample["sample_id"]
        require(sid == "bgc_" + text_sha256(SCHEMA_VERSION + "\0" + sample["archive_relpath"])[:24], "Stable sample ID differs")
        source = source_index[sample["source_relpath"]]
        require(all(sample[key] == source[key] for key in ("source_sha256", "size_bytes", "archive_relpath")), "Sample source identity differs")
        require(not is_true(sample["excluded"]) and json_array(sample, "repairs_applied_json") == [], "Sample excluded or repaired")
        rr = records_by_sample[sid]
        require(len(rr) == int(sample["record_count"]), f"Record count differs: {sid}")
        require(sorted(int(r["record_index"]) for r in rr) == list(range(1, len(rr) + 1)), "Non-contiguous record indices")
        for record_key, sample_key in (("length_bp", "total_record_length_bp"), ("cds_count", "cds_count"),
                                      ("asdomain_count", "asdomain_count"), ("pfam_domain_count", "pfam_domain_count")):
            require(sum(int(r[record_key]) for r in rr) == int(sample[sample_key]), f"Record aggregate differs: {sample_key}")
        for field in ("tsv_row_ids_json", "candidate_tsv_row_ids_json"):
            for rid in json_array(sample, field):
                require(rid in label_index, f"Unknown TSV row reference: {rid}")
                label = label_index[rid]
                require(label["dataset"] == sample["dataset"] and label["BGC_raw"] == sample["original_bgc_id"], "Sample/label identity differs")
        linked = json_array(sample, "tsv_row_ids_json")
        require(sorted({label_index[rid]["GCF_raw"] for rid in linked}) == json_array(sample, "tsv_labels_raw_json"), "Sample label display differs")
        if sample["mapping_status"] in {"exact_path_label", "unique_basename_folder_mismatch"}:
            require(len(linked) == 1 and sid in json_array(label_index[linked[0]], "linked_sample_ids_json"), "Resolved mapping not reciprocal")
        for route in ("dna", "protein"):
            require(sample[route + "_parser_status"] in {"ok", "error", "not_checked"}, "Unknown parser status")
        if sample["protein_parser_status"] != "ok":
            require(not sample["protein_input_sha256"] and not sample["protein_token_count"] and not sample["protein_residue_coverage"],
                    "Failed/unchecked protein parser retained success-only fields")

    # Compare raw TSV rows, including whitespace and physical line numbers, to the inventory.
    expected_label_rows = {}
    for path in sorted(source_root.glob("Dataset_*/curated_GCF_assignments.tsv")):
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle, delimiter="\t")
            require(next(reader, None) == ["BGC", "GCF"], "Raw label header differs")
            previous_line = reader.line_num
            for cells in reader:
                start, previous_line = previous_line + 1, reader.line_num
                rid = f"{path.parent.name}:tsv:{start}"
                expected_label_rows[rid] = (path.relative_to(workspace).as_posix(), start, reader.line_num, cells)
    require(set(expected_label_rows) == set(label_index), "Raw TSV/inventory row sets differ")
    for label in labels:
        path, start, end, cells = expected_label_rows[label["row_id"]]
        require(label["source_relpath"] == path and label["source_sha256"] == source_index[path]["source_sha256"], "Label source provenance differs")
        require((int(label["line_number"]), int(label["line_end"])) == (start, end), "Label source line numbers differ")
        require(json_array(label, "raw_cells_json") == cells, "Raw label cells changed")
        require([label["BGC_raw"], label["GCF_raw"]] == (cells + ["", ""])[:2], "Raw label values changed")
        for field in ("candidate_sample_ids_json", "linked_sample_ids_json"):
            for sid in json_array(label, field):
                require(sid in sample_index, "Unknown label/sample reference")
                sample = sample_index[sid]
                require((sample["dataset"], sample["original_bgc_id"]) == (label["dataset"], label["BGC_raw"]), "Label/sample identity differs")
                if field == "linked_sample_ids_json":
                    require(label["row_id"] in json_array(sample, "tsv_row_ids_json"), "Label link not reciprocal")

    groups = tables["duplicate_groups.tsv"]
    group_index = unique_index(groups, "group_id", "duplicate group ID")
    grouped = defaultdict(list)
    member_keys = set()
    for member in tables["duplicate_members.tsv"]:
        gid, sid, rid = member["group_id"], member["sample_id"], member["record_uid"]
        require(gid in group_index and sid in sample_index, "Duplicate member references unknown group/sample")
        require(member["group_type"] == group_index[gid]["group_type"], "Duplicate group type differs")
        require(not rid or (rid in record_index and record_index[rid]["sample_id"] == sid), "Duplicate record reference differs")
        key = (gid, sid, rid)
        require(key not in member_keys, "Repeated duplicate membership")
        member_keys.add(key)
        grouped[gid].append(member)
    for group in groups:
        members = grouped[group["group_id"]]
        sample_ids = sorted({m["sample_id"] for m in members})
        require(len(members) == int(group["member_count"]) >= 2, "Duplicate member count differs")
        require(len(sample_ids) == int(group["sample_count"]) and sample_ids == json_array(group, "sample_ids_json"), "Duplicate sample count/list differs")
        dataset_labels = defaultdict(set)
        for sid in sample_ids:
            sample = sample_index[sid]
            dataset_labels[sample["dataset"]].update(json_array(sample, "tsv_labels_raw_json"))
        require(json.loads(group["labels_by_dataset_json"]) == {k: sorted(v) for k, v in dataset_labels.items()}, "Duplicate label evidence differs")
        require(int(group["dataset_count"]) == len(dataset_labels) and is_true(group["cross_dataset"]) == (len(dataset_labels) > 1), "Duplicate dataset count differs")
        require(is_true(group["within_dataset_label_conflict"]) == any(len(v) > 1 for v in dataset_labels.values()), "Duplicate label conflict flag differs")
    issues = tables["issues.tsv"]
    for issue in issues:
        require(not issue["sample_id"] or issue["sample_id"] in sample_index, "Issue references unknown sample")
        details = json.loads(issue["details_json"])
        if "group_id" in details:
            require(details["group_id"] in group_index, "Issue references unknown group")
        if "row_id" in details:
            require(details["row_id"] in label_index, "Issue references unknown label row")

    dataset_rows = tables["dataset_summary.tsv"]
    dataset_index = unique_index(dataset_rows, "dataset", "dataset summary ID")
    require(set(dataset_index) == {row["dataset"] for row in samples + labels}, "Dataset summary set differs")
    for dataset, summary in dataset_index.items():
        ss = [s for s in samples if s["dataset"] == dataset]
        rr = [r for r in records if r["dataset"] == dataset]
        ll = [l for l in labels if l["dataset"] == dataset]
        names = Counter(s["original_bgc_id"] for s in ss)
        families = Counter(l["GCF_raw"] for l in ll if l["row_status"] == "ok")
        expected = {"file_count": len(ss), "record_count": len(rr), "tsv_row_count": len(ll),
                    "distinct_bgc_ids": len(names), "raw_gcf_count": len(families),
                    "raw_singleton_gcf_count": sum(n == 1 for n in families.values()),
                    "exact_label_matches": sum(s["mapping_status"] == "exact_path_label" for s in ss),
                    "folder_label_mismatches": sum(s["mapping_status"] == "unique_basename_folder_mismatch" for s in ss),
                    "unresolved_file_mappings": sum(s["mapping_status"] in {"ambiguous", "missing"} for s in ss),
                    "unlinked_tsv_rows": sum(l["mapping_status"] != "linked" for l in ll),
                    "repeated_basename_groups": sum(n > 1 for n in names.values()),
                    "extra_basename_entries": sum(n - 1 for n in names.values()),
                    "genbank_parse_ok": sum(s["genbank_status"] == "ok" for s in ss),
                    "multi_record_files": sum(int(s["record_count"]) > 1 for s in ss),
                    "dna_parser_ok": sum(s["dna_parser_status"] == "ok" for s in ss),
                    "protein_parser_ok": sum(s["protein_parser_status"] == "ok" for s in ss),
                    "minimum_record_bp": min((int(r["length_bp"]) for r in rr), default=""),
                    "maximum_record_bp": max((int(r["length_bp"]) for r in rr), default=""),
                    "records_without_asdomains": sum(int(r["asdomain_count"]) == 0 for r in rr),
                    "excluded_count": 0, "repaired_count": 0}
        require(all(summary[key] == str(value) for key, value in expected.items()), f"Dataset counts differ: {dataset}")
    counts = {"source_files": len(source_index), "gbk_files": len(samples), "records": len(records),
              "tsv_rows": len(labels), "datasets": len(dataset_rows),
              "genbank_parse_ok": sum(s["genbank_status"] == "ok" for s in samples),
              "dna_parser_ok": sum(s["dna_parser_status"] == "ok" for s in samples),
              "protein_parser_ok": sum(s["protein_parser_status"] == "ok" for s in samples),
              "multi_record_files": sum(int(s["record_count"]) > 1 for s in samples),
              "label_mapping_status": dict(Counter(s["mapping_status"] for s in samples)),
              "issue_counts": dict(Counter(i["kind"] for i in issues)),
              "duplicate_group_counts": dict(Counter(g["group_type"] for g in groups)),
              "cross_dataset_duplicate_groups": dict(Counter(g["group_type"] for g in groups if is_true(g["cross_dataset"]))),
              "excluded_count": 0, "repaired_count": 0}
    require(counts == metadata["counts"], "Metadata counts differ from tables")
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--inventory", type=Path, default=Path("data/processed/inventory_v1"))
    args = parser.parse_args(argv)
    try:
        counts = verify_inventory(args.workspace, args.inventory)
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as exc:
        parser.exit(1, f"Inventory verification failed: {exc}\n")
    print(f"Verified: {counts['source_files']} unchanged raw files; {counts['gbk_files']} GBK entries; "
          f"{counts['records']} records; {counts['tsv_rows']} label rows; {counts['datasets']} datasets.")
    print("Output hashes, relocation provenance, references and counts agree. No data were modified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
