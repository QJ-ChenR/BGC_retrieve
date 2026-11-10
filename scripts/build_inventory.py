#!/usr/bin/env python3
"""Inventory curated GBKs without repairing, filtering, or selecting samples.

Each GBK file is a source entry, not necessarily an independent biological BGC.
Records and TSV rows retain their own identities. No sequence is concatenated,
translated, cropped, or encoded by a model. Outputs are published to a new directory.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import importlib
import json
import platform
from pathlib import Path
import sys
import tempfile
import warnings

import Bio
from Bio import SeqIO
from Bio.Seq import Seq


SCHEMA_VERSION = "curated_bgc_inventory_v1"
GBK_SUFFIXES = {".gbk", ".gb", ".genbank"}
DNA_ALPHABET = set("ACGTRYSWKMBDHVN")
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZUO")


def compact(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def digest(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def file_sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def sample_identifier(archive_relpath):
    return "bgc_" + digest(SCHEMA_VERSION + "\0" + archive_relpath)[:24]


def relative(path, workspace):
    return Path(path).relative_to(workspace).as_posix()


def snapshot(source_root, workspace):
    if source_root.is_symlink() or not source_root.is_dir():
        raise ValueError(f"Expected a real source directory: {source_root}")
    rows = []
    for path in sorted(source_root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Source symlinks are unsupported: {path}")
        if path.is_dir():
            continue
        if not path.is_file():
            raise ValueError(f"Source special files are unsupported: {path}")
        before = path.stat()
        checksum = file_sha256(path)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"Source changed while hashing: {path}")
        rows.append({"archive_relpath": path.relative_to(source_root).as_posix(),
                     "source_relpath": relative(path, workspace),
                     "size_bytes": before.st_size, "source_sha256": checksum})
    if not rows:
        raise ValueError("Source directory has no files")
    return rows


def write_tsv(path, rows, leading=()):
    fields = list(dict.fromkeys(leading))
    fields += sorted(set().union(*(row.keys() for row in rows)) - set(fields)) if rows else []
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return fields


def write_json(path, value):
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")


def issue(issues, kind, sample_id="", dataset="", details=None, severity="review"):
    issues.append({"kind": kind, "sample_id": sample_id, "dataset": dataset,
                   "severity": severity, "details_json": compact(details or {})})


def read_labels(source_root, workspace, sources, issues):
    rows = []
    for path in sorted(source_root.glob("Dataset_*/curated_GCF_assignments.tsv")):
        dataset = path.parent.name
        source = sources[path.relative_to(source_root).as_posix()]
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header != ["BGC", "GCF"]:
                raise ValueError(f"Unexpected TSV header at {path}: {header}")
            previous_line = reader.line_num
            for cells in reader:
                start_line = previous_line + 1
                previous_line = reader.line_num
                valid = len(cells) == 2 and all(cell.strip() for cell in cells)
                bgc, label = (cells + ["", ""])[:2]
                row = {"row_id": f"{dataset}:tsv:{start_line}", "dataset": dataset,
                       "source_relpath": relative(path, workspace),
                       "source_sha256": source["source_sha256"],
                       "line_number": start_line, "line_end": reader.line_num,
                       "BGC_raw": bgc, "GCF_raw": label,
                       "raw_cells_json": compact(cells),
                       "label_has_edge_whitespace": label != label.strip(),
                       "bgc_has_edge_whitespace": bgc != bgc.strip(),
                       "row_status": "ok" if valid else "invalid",
                       "candidate_sample_ids_json": "[]", "linked_sample_ids_json": "[]"}
                rows.append(row)
                if not valid:
                    issue(issues, "invalid_label_row", dataset=dataset,
                          details={"row_id": row["row_id"], "cells": cells})
                if row["label_has_edge_whitespace"] or row["bgc_has_edge_whitespace"]:
                    issue(issues, "raw_label_whitespace", dataset=dataset,
                          details={"row_id": row["row_id"], "BGC_raw": bgc, "GCF_raw": label},
                          severity="info")
    return rows


def associate_labels(samples, labels, issues):
    sample_names, label_names = defaultdict(list), defaultdict(list)
    for sample in samples:
        sample_names[(sample["dataset"], sample["original_bgc_id"])].append(sample)
    for label in labels:
        if label["row_status"] == "ok":
            label_names[(label["dataset"], label["BGC_raw"])].append(label)
    links = defaultdict(list)
    for sample in samples:
        key = (sample["dataset"], sample["original_bgc_id"])
        candidates, same_name = label_names[key], sample_names[key]
        exact = [row for row in candidates if row["GCF_raw"] == sample["folder_gcf_raw"]]
        exact_files = [row for row in same_name if row["folder_gcf_raw"] == sample["folder_gcf_raw"]]
        if len(exact) == len(exact_files) == 1:
            status, associated = "exact_path_label", exact
        elif len(same_name) == len(candidates) == 1:
            status, associated = "unique_basename_folder_mismatch", candidates
        else:
            status = "ambiguous" if candidates else "missing"
            associated = candidates  # Candidates are retained, NOT resolved labels.
        sample.update(mapping_status=status,
                      tsv_row_ids_json=compact([row["row_id"] for row in associated]),
                      tsv_labels_raw_json=compact(sorted({row["GCF_raw"] for row in associated})),
                      candidate_tsv_row_ids_json=compact([row["row_id"] for row in candidates]),
                      basename_file_count=len(same_name), basename_tsv_row_count=len(candidates),
                      basename_has_multiple_labels=len({row["GCF_raw"] for row in candidates}) > 1)
        if status in {"exact_path_label", "unique_basename_folder_mismatch"}:
            for row in associated:
                links[row["row_id"]].append(sample["sample_id"])
        if status != "exact_path_label":
            issue(issues, "label_mapping_" + status, sample["sample_id"], sample["dataset"],
                  {"folder_gcf_raw": sample["folder_gcf_raw"],
                   "candidate_rows": [row["row_id"] for row in candidates],
                   "candidate_labels": [row["GCF_raw"] for row in candidates]})
    for row in labels:
        candidates = sample_names[(row["dataset"], row["BGC_raw"])]
        row["candidate_sample_ids_json"] = compact([s["sample_id"] for s in candidates])
        row["linked_sample_ids_json"] = compact(links[row["row_id"]])
        row["mapping_status"] = ("linked" if links[row["row_id"]]
                                 else "unresolved" if candidates else "missing_file")
        if not links[row["row_id"]]:
            issue(issues, "label_row_" + row["mapping_status"], dataset=row["dataset"],
                  details={"row_id": row["row_id"]})


def record_inventory(record, index, sample):
    row = {"sample_id": sample["sample_id"], "dataset": sample["dataset"],
           "record_uid": f"{sample['sample_id']}:record:{index}", "record_index": index,
           "record_id_raw": record.id, "record_name_raw": record.name,
           "description_raw": record.description, "length_bp": len(record.seq),
           "organism_raw": record.annotations.get("organism", ""),
           "accessions_json": compact(record.annotations.get("accessions", [])),
           "taxonomy_json": compact(record.annotations.get("taxonomy", [])),
           "topology": record.annotations.get("topology", "unknown"),
           "sequence_sha256": "", "sequence_rc_canonical_sha256": "",
           "sequence_status": "ok", "sequence_error": "", "gc_fraction_acgt": "",
           "ambiguous_base_count": "", "ambiguous_base_fraction": "",
           "ambiguous_base_counts_json": "{}", "unsupported_dna_chars_json": "[]"}
    try:
        sequence = str(record.seq).upper()
        counts = Counter(sequence)
        acgt = sum(counts[base] for base in "ACGT")
        invalid = sorted(set(sequence) - DNA_ALPHABET)
        row.update(sequence_sha256=digest(sequence),
                   gc_fraction_acgt=(counts["G"] + counts["C"]) / acgt if acgt else "",
                   ambiguous_base_count=len(sequence) - acgt,
                   ambiguous_base_fraction=(len(sequence) - acgt) / len(sequence) if sequence else "",
                   ambiguous_base_counts_json=compact({k: v for k, v in counts.items() if k not in "ACGT"}),
                   unsupported_dna_chars_json=compact(invalid))
        if sequence and not invalid:
            reverse = str(Seq(sequence).reverse_complement())
            row["sequence_rc_canonical_sha256"] = digest(min(sequence, reverse))
        else:
            row["sequence_status"] = "invalid_alphabet" if invalid else "empty"
    except (ValueError, TypeError, UnicodeError) as exc:
        row.update(sequence_status="error", sequence_error=f"{type(exc).__name__}: {exc}")
    cds = [feature for feature in record.features if feature.type == "CDS"]
    annotations = [{"type": feature.type, "location": str(feature.location),
                    "qualifiers": feature.qualifiers} for feature in record.features]
    row["feature_annotation_sha256"] = digest(compact(annotations))
    row.update(cds_count=len(cds), cds_with_translation_count=0,
               cds_missing_translation_count=0, cds_multiple_translation_count=0,
               cds_terminal_stop_count=0, cds_internal_stop_count=0,
               cds_invalid_translation_count=0, cds_pseudo_count=0,
               biosynthetic_cds_count=0, biosynthetic_additional_cds_count=0,
               asdomain_count=sum(f.type == "aSDomain" for f in record.features),
               pfam_domain_count=sum(f.type == "PFAM_domain" for f in record.features))
    for feature in cds:
        values = feature.qualifiers.get("translation", [])
        if not values or not any(value.strip() for value in values):
            row["cds_missing_translation_count"] += 1
        else:
            row["cds_with_translation_count"] += 1
        row["cds_multiple_translation_count"] += len(values) > 1
        translations = ["".join(value.split()).upper() for value in values]
        row["cds_terminal_stop_count"] += any(value.endswith("*") for value in translations)
        row["cds_internal_stop_count"] += any("*" in value[:-1] for value in translations)
        row["cds_invalid_translation_count"] += any(set(value) - AA_ALPHABET for value in translations)
        row["cds_pseudo_count"] += any(key in feature.qualifiers for key in ("pseudo", "pseudogene"))
        row["biosynthetic_cds_count"] += "biosynthetic" in feature.qualifiers.get("gene_kind", [])
        row["biosynthetic_additional_cds_count"] += "biosynthetic-additional" in feature.qualifiers.get("gene_kind", [])
    taxids = {value.removeprefix("taxon:") for feature in record.features if feature.type == "source"
              for value in feature.qualifiers.get("db_xref", []) if value.startswith("taxon:")}
    row["taxon_ids_json"] = compact(sorted(taxids))
    return row


def load_parsers(encoder_project):
    if encoder_project is None:
        return {}, {"enabled": False}
    encoder_project = Path(encoder_project).resolve()
    names = ["data_preparation/__init__.py", "data_preparation/feature_inputs.py",
             "data_preparation/dna_inputs.py", "data_preparation/feature_cache.py"]
    provenance = {"enabled": True, "project": str(encoder_project),
                  "source_sha256": {name: file_sha256(encoder_project / name) for name in names}}
    sys.path.insert(0, str(encoder_project))
    try:
        protein = importlib.import_module("data_preparation.feature_inputs")
        dna = importlib.import_module("data_preparation.dna_inputs")
    finally:
        sys.path.pop(0)
    for module in (protein, dna):
        if not Path(module.__file__).resolve().is_relative_to(encoder_project):
            raise ValueError("A different data_preparation package is already imported")
    provenance.update(protein_parser_version=protein.PARSER_VERSION,
                      dna_parser_version=dna.DNA_INPUT_VERSION)
    return {"dna": dna.read_dna, "protein": protein.parse_bgc}, provenance


def audit_sample(sample, path, parsers, issues):
    sample.update(genbank_status="ok", genbank_error="", record_count=0,
                  records_complete=False, parse_warnings_json="[]",
                  ordered_file_dna_sha256="", single_record_dna_rc_sha256="",
                  protein_input_sha256="", total_record_length_bp=0,
                  cds_count=0, asdomain_count=0, pfam_domain_count=0,
                  protein_token_count="", protein_residue_coverage="",
                  review_status="inventory_only", repairs_applied_json="[]", excluded=False)
    rows = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            with path.open(encoding="utf-8") as handle:
                for index, record in enumerate(SeqIO.parse(handle, "genbank"), 1):
                    rows.append(record_inventory(record, index, sample))
            if not rows:
                raise ValueError("No GenBank records found")
            sample["records_complete"] = True
        except Exception as exc:
            sample.update(genbank_status="error", genbank_error=f"{type(exc).__name__}: {exc}")
            issue(issues, "genbank_parse_error", sample["sample_id"], sample["dataset"],
                  {"error": sample["genbank_error"], "records_recovered": len(rows)})
    sample["parse_warnings_json"] = compact(
        [{"category": w.category.__name__, "message": str(w.message)} for w in caught])
    for warning in caught:
        issue(issues, "genbank_warning", sample["sample_id"], sample["dataset"],
              {"category": warning.category.__name__, "message": str(warning.message)}, "info")
    sample["record_count"] = len(rows)
    sample["total_record_length_bp"] = sum(row["length_bp"] for row in rows)
    for key in ("cds_count", "asdomain_count", "pfam_domain_count"):
        sample[key] = sum(row[key] for row in rows)
    if sample["records_complete"] and all(row["sequence_status"] == "ok" for row in rows):
        # Hash a framed LIST of record hashes, not concatenated DNA.
        sample["ordered_file_dna_sha256"] = digest(compact([row["sequence_sha256"] for row in rows]))
        if len(rows) == 1:
            sample["single_record_dna_rc_sha256"] = rows[0]["sequence_rc_canonical_sha256"]
    if len(rows) > 1:
        issue(issues, "multi_record_file", sample["sample_id"], sample["dataset"],
              {"record_count": len(rows), "policy": "retain records; do not split labels or concatenate"})
    for row in rows:
        if row["sequence_status"] != "ok":
            issue(issues, "record_sequence_" + row["sequence_status"], sample["sample_id"], sample["dataset"],
                  {"record_uid": row["record_uid"], "error": row["sequence_error"],
                   "unsupported_chars": json.loads(row["unsupported_dna_chars_json"])})
    for route in ("dna", "protein"):
        sample[route + "_parser_status"] = "not_checked"
        sample[route + "_parser_error"] = ""
        sample[route + "_parser_warnings_json"] = "[]"
        if route not in parsers:
            continue
        with warnings.catch_warnings(record=True) as parser_warnings:
            warnings.simplefilter("always")
            try:
                result = parsers[route](path)
                if route == "protein":
                    tokens = result["tokens"]
                    token_count = len(tokens)
                    token_hash = digest(compact([t["sequence"] for t in tokens]))
                    parents = result["parents"]
                    denominator = sum(parent["protein_length"] for parent in parents)
                    coverage = (
                        sum(parent["represented_residues"] for parent in parents) / denominator
                        if denominator else "")
                    sample.update(protein_token_count=token_count, protein_input_sha256=token_hash,
                                  protein_residue_coverage=coverage)
                sample[route + "_parser_status"] = "ok"
            except Exception as exc:
                message = f"{type(exc).__name__}: {exc}".replace(str(path), sample["source_relpath"])
                sample[route + "_parser_status"] = "error"
                sample[route + "_parser_error"] = message
                issue(issues, route + "_parser_error", sample["sample_id"], sample["dataset"], {"error": message})
        sample[route + "_parser_warnings_json"] = compact(
            [{"category": w.category.__name__, "message": str(w.message)} for w in parser_warnings])
    return rows


def duplicate_inventory(samples, records, issues):
    by_id = {row["sample_id"]: row for row in samples}
    groups, members = [], []
    definitions = [("basename", samples, lambda r: compact([r["dataset"], r["original_bgc_id"]])),
                   ("file_bytes", samples, lambda r: r["source_sha256"]),
                   ("ordered_file_dna", samples, lambda r: r["ordered_file_dna_sha256"]),
                   ("single_record_dna_rc", samples, lambda r: r["single_record_dna_rc_sha256"]),
                   ("record_dna", records, lambda r: r["sequence_sha256"] if r["sequence_status"] == "ok" else ""),
                   ("record_dna_rc", records, lambda r: r["sequence_rc_canonical_sha256"]),
                   ("protein_input", samples, lambda r: r["protein_input_sha256"])]
    for kind, source_rows, get_key in definitions:
        buckets = defaultdict(list)
        for row in source_rows:
            key = get_key(row)
            if key:
                buckets[key].append(row)
        for key, rows in sorted(buckets.items()):
            if len(rows) < 2:
                continue
            ids = sorted({row["sample_id"] for row in rows})
            labels_by_dataset = defaultdict(set)
            for sid in ids:
                sample = by_id[sid]
                labels_by_dataset[sample["dataset"]].update(json.loads(sample["tsv_labels_raw_json"]))
            conflict = any(len(labels) > 1 for labels in labels_by_dataset.values())
            gid = f"{kind}_" + digest(key)[:24]
            groups.append({"group_id": gid, "group_type": kind, "key": key,
                           "member_count": len(rows), "sample_count": len(ids),
                           "dataset_count": len(labels_by_dataset),
                           "sample_ids_json": compact(ids),
                           "labels_by_dataset_json": compact({k: sorted(v) for k, v in labels_by_dataset.items()}),
                           "within_dataset_label_conflict": conflict,
                           "cross_dataset": len(labels_by_dataset) > 1})
            for row in rows:
                members.append({"group_id": gid, "group_type": kind,
                                "sample_id": row["sample_id"], "record_uid": row.get("record_uid", "")})
            if kind == "basename":
                issue(issues, "repeated_basename", dataset=rows[0]["dataset"],
                      details={"group_id": gid, "sample_ids": ids})
            elif kind in {"ordered_file_dna", "single_record_dna_rc"} and conflict:
                issue(issues, "sequence_label_conflict_candidate", details={"group_id": gid, "sample_ids": ids})
    return groups, members


def summarize(samples, records, labels, groups):
    datasets = sorted({row["dataset"] for row in samples + labels})
    result = []
    for dataset in datasets:
        ss = [row for row in samples if row["dataset"] == dataset]
        rr = [row for row in records if row["dataset"] == dataset]
        ll = [row for row in labels if row["dataset"] == dataset]
        names = Counter(row["original_bgc_id"] for row in ss)
        families = Counter(row["GCF_raw"] for row in ll if row["row_status"] == "ok")
        result.append({"dataset": dataset, "file_count": len(ss), "tsv_row_count": len(ll),
                       "distinct_bgc_ids": len(names), "raw_gcf_count": len(families),
                       "raw_singleton_gcf_count": sum(count == 1 for count in families.values()),
                       "exact_label_matches": sum(s["mapping_status"] == "exact_path_label" for s in ss),
                       "folder_label_mismatches": sum(s["mapping_status"] == "unique_basename_folder_mismatch" for s in ss),
                       "unresolved_file_mappings": sum(s["mapping_status"] in {"ambiguous", "missing"} for s in ss),
                       "unlinked_tsv_rows": sum(l["mapping_status"] != "linked" for l in ll),
                       "repeated_basename_groups": sum(count > 1 for count in names.values()),
                       "extra_basename_entries": sum(count - 1 for count in names.values()),
                       "genbank_parse_ok": sum(s["genbank_status"] == "ok" for s in ss),
                       "record_count": len(rr), "multi_record_files": sum(s["record_count"] > 1 for s in ss),
                       "dna_parser_ok": sum(s["dna_parser_status"] == "ok" for s in ss),
                       "protein_parser_ok": sum(s["protein_parser_status"] == "ok" for s in ss),
                       "minimum_record_bp": min((r["length_bp"] for r in rr), default=""),
                       "maximum_record_bp": max((r["length_bp"] for r in rr), default=""),
                       "records_without_asdomains": sum(r["asdomain_count"] == 0 for r in rr),
                       "excluded_count": 0, "repaired_count": 0})
    return result


def build_inventory(source_root, output_dir, workspace, encoder_project=None):
    source_root, output_dir, workspace = map(lambda p: Path(p).absolute(), (source_root, output_dir, workspace))
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(f"Refusing to overwrite inventory: {output_dir}")
    if source_root.is_relative_to(output_dir) or output_dir.is_relative_to(source_root):
        raise ValueError("Source and output directories must not contain one another")
    for directory in (source_root, output_dir.parent):
        for parent in (directory, *directory.parents):
            if parent.is_symlink():
                raise ValueError(f"Symlink directory component is unsupported: {parent}")
    source_rows = snapshot(source_root, workspace)
    sources = {row["archive_relpath"]: row for row in source_rows}
    issues, samples, records = [], [], []
    for source in source_rows:
        path = source_root / source["archive_relpath"]
        if path.suffix.lower() not in GBK_SUFFIXES:
            continue
        parts = Path(source["archive_relpath"]).parts
        if len(parts) != 4 or parts[1] != "GCF_assignments" or not parts[0].startswith("Dataset_"):
            raise ValueError(f"Unexpected GBK layout; refusing to ignore or reinterpret it: {path}")
        samples.append({"sample_id": sample_identifier(source["archive_relpath"]),
                        "dataset": parts[0], "original_bgc_id": path.stem,
                        "folder_gcf_raw": parts[2], **source})
    if not samples:
        raise ValueError("No GBK source entries found")
    if len({s["sample_id"] for s in samples}) != len(samples):
        raise ValueError("Sample ID collision")
    labels = read_labels(source_root, workspace, sources, issues)
    associate_labels(samples, labels, issues)
    parsers, parser_provenance = load_parsers(encoder_project)
    for index, sample in enumerate(samples, 1):
        records.extend(audit_sample(sample, source_root / sample["archive_relpath"], parsers, issues))
        if index % 250 == 0:
            print(f"Inventoried {index}/{len(samples)} GBK files", flush=True)
    groups, members = duplicate_inventory(samples, records, issues)
    dataset_rows = summarize(samples, records, labels, groups)
    # The snapshot covers TSVs and README as well as all GBKs.
    if snapshot(source_root, workspace) != source_rows:
        raise ValueError("Raw source changed during inventory; no completed output published")
    if parser_provenance["enabled"]:
        for name, checksum in parser_provenance["source_sha256"].items():
            if file_sha256(Path(parser_provenance["project"]) / name) != checksum:
                raise ValueError(f"Parser source changed during inventory: {name}")
    counts = {"source_files": len(source_rows), "gbk_files": len(samples),
              "records": len(records), "tsv_rows": len(labels), "datasets": len(dataset_rows),
              "genbank_parse_ok": sum(s["genbank_status"] == "ok" for s in samples),
              "dna_parser_ok": sum(s["dna_parser_status"] == "ok" for s in samples),
              "protein_parser_ok": sum(s["protein_parser_status"] == "ok" for s in samples),
              "multi_record_files": sum(s["record_count"] > 1 for s in samples),
              "label_mapping_status": dict(Counter(s["mapping_status"] for s in samples)),
              "issue_counts": dict(Counter(row["kind"] for row in issues)),
              "duplicate_group_counts": dict(Counter(row["group_type"] for row in groups)),
              "cross_dataset_duplicate_groups": dict(Counter(row["group_type"] for row in groups if row["cross_dataset"])),
              "excluded_count": 0, "repaired_count": 0}
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="." + output_dir.name + ".staging-", dir=output_dir.parent))
    tables = [("source_files.tsv", source_rows, ("source_relpath", "archive_relpath", "source_sha256", "size_bytes")),
              ("samples.tsv", samples, ("sample_id", "dataset", "source_relpath", "archive_relpath")),
              ("records.tsv", records, ("sample_id", "record_uid", "record_index", "dataset")),
              ("label_rows.tsv", labels, ("row_id", "dataset", "source_relpath", "line_number", "BGC_raw", "GCF_raw")),
              ("duplicate_groups.tsv", groups, ("group_id", "group_type", "member_count", "sample_count", "within_dataset_label_conflict")),
              ("duplicate_members.tsv", members, ("group_id", "group_type", "sample_id", "record_uid")),
              ("issues.tsv", issues, ("kind", "severity", "sample_id", "dataset", "details_json")),
              ("dataset_summary.tsv", dataset_rows, ("dataset", "file_count", "tsv_row_count", "record_count"))]
    schemas = {name: write_tsv(staging / name, rows, leading) for name, rows, leading in tables}
    summary = ["# 原始数据样本清单（尚未筛选评测队列）", "",
               f"已保留 {len(samples)} 个 GBK 文件条目、{len(records)} 个 GenBank records、{len(labels)} 条原始标签行。",
               "没有长度过滤、去重删除、标签修正、翻译推导或模型编码；文件条目不等于独立生物学样本。", "",
               "| 数据集 | GBK | TSV 行 | records | DNA 兼容 | 蛋白兼容 | 多 record 文件 |",
               "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    summary += [f"| {r['dataset']} | {r['file_count']} | {r['tsv_row_count']} | {r['record_count']} | "
                f"{r['dna_parser_ok']} | {r['protein_parser_ok']} | {r['multi_record_files']} |" for r in dataset_rows]
    summary += ["", "兼容列仅为现有编码器的单文件 CPU 解析检查，未运行批量发现或模型。未启用检查时，计数为 0，状态为 not_checked。",
                "基础 GBK 读取成功与现有编码器兼容是两个不同的状态。失败条目仍留在清单内。", "",
                "## 标签与重复", "",
                "标签关联状态：`" + compact(counts["label_mapping_status"]) + "`。",
                "重复组数量：`" + compact(counts["duplicate_group_counts"]) + "`。",
                "跨数据集重复组数量：`" + compact(counts["cross_dataset_duplicate_groups"]) + "`。", "",
                "精确文件—TSV 对应不代表生物学标签冲突已解决。unique_basename_folder_mismatch 仅关联唯一名称，保留目录与 TSV 的不同原值。",
                "ambiguous 行保留候选 TSV 行和标签，不把它们认定为已验证的多标签真值。标签空白（含 NBSP）不被修改。",
                "不同数据集的 GCF 标签各自独立；within_dataset_label_conflict 仅提示同一数据集内不同标签，仍需人工核实。",
                "record 级重复不等于整个文件／BGC 重复；反向互补分组只用于审计，未翻转输入。未开展近同源聚类。", "",
                "## 追溯", "",
                "samples.tsv 通过 sample_id 连接 records.tsv、duplicate_members.tsv 和 issues.tsv；通过 tsv_row_ids_json 连接 label_rows.tsv。",
                "source_files.tsv 包含全部原始文件的路径、字节数和 SHA256。原始目录移动前后的映射见 ../raw_relocation/files.tsv。",
                "样本 ID 是 archive_relpath 的带版本哈希，独立于 raw 的绝对路径；源文件内容变化由 source_sha256 检出。",
                "ordered_file_dna_sha256 哈希的是有序 record 序列哈希列表，没有拼接多条序列。",
                "protein_input_sha256 只在现有蛋白解析器成功时记录其有序输入 token 序列，不代表模型输出或完整 CDS 序列。",
                "metadata.json 记录脚本、解析器源码、版本、策略和各输出校验和。", ""]
    (staging / "summary.md").write_text("\n".join(summary), encoding="utf-8")
    metadata = {"schema_version": SCHEMA_VERSION, "status": "complete",
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "source_root": relative(source_root, workspace),
                "output_dir": relative(output_dir, workspace),
                "script_sha256": file_sha256(Path(__file__)), "python_version": platform.python_version(),
                "biopython_version": Bio.__version__, "encoder_parsers": parser_provenance,
                "source_snapshot_sha256": digest(compact(source_rows)),
                "source_verified_unchanged_after_inventory": True,
                "policies": {"sample_unit": "source GBK file; not confirmed biological independence",
                             "record_index": "one-based, original file order",
                             "sample_id": "bgc_ + first 24 hex of sha256(schema_version + NUL + archive_relpath)",
                             "label_matching": "exact raw dataset/basename/folder first; unique basename mismatch explicit; otherwise candidates",
                             "label_normalization": "none; raw Unicode whitespace preserved",
                             "dna_hash": "uppercase as-stored sequence; SHA256 UTF-8",
                             "rc_hash": "SHA256 lexical minimum of uppercase sequence and reverse complement; audit only",
                             "gc_denominator": "A+C+G+T only",
                             "feature_annotation_hash": "ordered feature type/location/qualifiers as canonical JSON",
                             "length_filter": None, "repairs": "none", "exclusions": "none",
                             "deduplication": "groups reported; no entries removed", "homology_analysis": "not performed",
                             "evaluation_cohort": "not finalized; inventory is not an encoder input directory"},
                "counts": counts, "table_columns": schemas,
                "output_sha256": {path.name: file_sha256(path) for path in sorted(staging.iterdir())}}
    write_json(staging / "metadata.json", metadata)
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(f"Output appeared while running; staging retained at {staging}")
    staging.rename(output_dir)
    return counts


def main(argv=None):
    default_workspace = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=default_workspace)
    parser.add_argument("--source", type=Path, default=Path("data/raw/curated_datasets"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/inventory_v1"))
    parser.add_argument("--encoder-project", type=Path, help="Enable CPU-only single-file parser checks")
    args = parser.parse_args(argv)
    workspace = args.workspace.absolute()
    source = args.source if args.source.is_absolute() else workspace / args.source
    output = args.output if args.output.is_absolute() else workspace / args.output
    try:
        counts = build_inventory(source, output, workspace, args.encoder_project)
    except (ValueError, OSError, ImportError, csv.Error) as exc:
        parser.exit(1, f"Inventory failed: {exc}\n")
    print(json.dumps(counts, indent=2, ensure_ascii=False))
    print(f"Published inventory: {output}")


if __name__ == "__main__":
    main()
