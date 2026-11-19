#!/usr/bin/env python3
"""Publish sequence-proven G/H coordinate repairs under processed; never edit raw.

Preserves original basenames, family folders, labels, DNA and translations.
Only proven domain endpoint numbers are patched in the original GBK bytes.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import sys
import tempfile
import warnings

import Bio
from Bio import SeqIO
from Bio.SeqFeature import SimpleLocation

from build_inventory import (SCHEMA_VERSION, compact, digest, file_sha256,
                             load_parsers, sample_identifier, snapshot, write_json, write_tsv)
from domain_coordinate_repair import apply_byte_repairs, plan_repairs
from pilot_features import read_tsv, safe_child

ROOT = Path(__file__).resolve().parents[1]
DATASETS = {"Dataset_G_Planomonospora": 11, "Dataset_H_Reclassification": 90}


def resolve(workspace, path):
    return (workspace / path).expanduser().resolve()


def parsed_record(data):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        record = SeqIO.read(io.StringIO(data.decode("utf-8")), "genbank")
    return record, [{"category": w.category.__name__, "message": str(w.message)} for w in caught]


def check_preserved(before, after, changes):
    """Reject any semantic difference beyond the explicitly supported endpoints."""
    for name in ("id", "name", "description", "seq", "annotations", "dbxrefs", "letter_annotations"):
        if getattr(before, name) != getattr(after, name):
            raise ValueError(f"Repair changed record {name}")
    if len(before.features) != len(after.features):
        raise ValueError("Repair changed feature count")
    by_index = {r["feature_index"]: r for r in changes}
    if len(by_index) != len(changes):
        raise ValueError("Multiple repairs for one feature")
    for index, (old, new) in enumerate(zip(before.features, after.features)):
        qualifiers = dict(old.qualifiers)
        location = old.location
        if index in by_index:
            change = by_index[index]
            location = SimpleLocation(change["new_genomic_start"], change["new_genomic_end"],
                                      strand=change["strand"])
            if change["new_protein_end"] is not None:
                qualifiers["protein_end"] = [str(change["new_protein_end"])]
        if (old.type != new.type or old.id != new.id or qualifiers != new.qualifiers
                or location != new.location):
            raise ValueError(f"Unexpected annotation change at feature {index}")


def parse_protein(path, parser):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result = parser(path)
            status, error = "ok", ""
        except ValueError as exc:
            result, status, error = None, "error", str(exc).replace(str(path), path.name)
    return result, status, error, [{"category": w.category.__name__, "message": str(w.message)} for w in caught]


def token_digest(result):
    return digest(compact([t["sequence"] for t in result["tokens"]])) if result else ""


def load_inventory(workspace, raw, inventory):
    metadata = json.loads((inventory / "metadata.json").read_text(encoding="utf-8"))
    if metadata.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported inventory schema")
    hashes = {}
    for name in ("samples.tsv", "source_files.tsv"):
        hashes[name] = file_sha256(inventory / name)
        if metadata["output_sha256"].get(name) != hashes[name]:
            raise ValueError(f"Inventory changed: {name}")
    rows = [r for r in read_tsv(inventory / "samples.tsv") if r["dataset"] in DATASETS]
    if Counter(r["dataset"] for r in rows) != Counter(DATASETS):
        raise ValueError("G/H inventory sample counts differ from expected 11/90")
    if len({r["sample_id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate inventory sample ID")
    for row in rows:
        labels = json.loads(row["tsv_labels_raw_json"])
        if (row["mapping_status"] != "exact_path_label" or row["record_count"] != "1"
                or row["excluded"] != "False" or row["repairs_applied_json"] != "[]"
                or len(labels) != 1 or labels[0] != row["folder_gcf_raw"]
                or sample_identifier(row["archive_relpath"]) != row["sample_id"]):
            raise ValueError(f"Unresolved identity, label or record policy: {row['sample_id']}")
        row["gcf_raw"] = labels[0]
        path = safe_child(raw, row["archive_relpath"])
        if path.stem != row["original_bgc_id"] or path.parts[-3] != "GCF_assignments":
            raise ValueError(f"Unexpected raw layout: {path}")
    for dataset in DATASETS:
        subset = [r for r in rows if r["dataset"] == dataset]
        if len({r["original_bgc_id"] for r in subset}) != len(subset):
            raise ValueError(f"Duplicate original IDs in {dataset}")
        labels = read_tsv(raw / dataset / "curated_GCF_assignments.tsv")
        if sorted((r["BGC"], r["GCF"]) for r in labels) != sorted(
                (r["original_bgc_id"], r["gcf_raw"]) for r in subset):
            raise ValueError(f"Raw labels disagree with inventory: {dataset}")
    before = snapshot(raw, workspace)
    sources = {r["archive_relpath"]: r for r in read_tsv(inventory / "source_files.tsv")}
    if {r["archive_relpath"] for r in before} != set(sources):
        raise ValueError("Raw file set differs from source inventory")
    for row in before:
        old = sources[row["archive_relpath"]]
        if (row["source_sha256"] != old["source_sha256"]
                or row["size_bytes"] != int(old["size_bytes"])):
            raise ValueError(f"Raw file differs from inventory: {row['archive_relpath']}")
    return sorted(rows, key=lambda r: (r["dataset"], r["archive_relpath"])), before, hashes


def summary_text(dataset_rows, repairs, review, raw_count):
    protein_edits = sum(r["new_protein_end"] is not None for r in repairs)
    lines = ["# G/H 数据集坐标修复", "",
             "原始数据保留在 data/raw/curated_datasets；本目录是单独生成的修复副本。",
             "目录、原始 GBK 文件名、GCF 标签、DNA、CDS/domain 翻译及其他注释保持不变。", "",
             "| 数据集 | GBK 数 | 修改文件数 | 修复 domain 数 | 修复前蛋白解析通过 | 修复后蛋白解析通过 | 修复后 DNA 解析通过 |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for r in dataset_rows:
        lines.append(f"| {r['dataset']} | {r['sample_count']} | {r['changed_files']} | {r['repaired_domains']} | "
                     f"{r['protein_ok_before']} | {r['protein_ok_after']} | {r['dna_ok_after']} |")
    lines += ["", f"共修改 {len(repairs)} 个 domain 的基因组末端，其中 {protein_edits} 个同步修改 protein_end；",
              "其余原本没有蛋白坐标的注释只修正基因组区间，没有新增 protein_start/protein_end。", "",
              "## 修复依据", "",
              "- 唯一、同链且包含 domain 的父 CDS；已有标识符关系无冲突。",
              "- domain 翻译与父 CDS 的对应蛋白区间完全匹配。",
              "- 原 domain DNA 翻译等于已存 domain 翻译去掉最后一个残基。",
              "- 沿编码方向延伸一个密码子后，DNA 翻译与已存 domain 翻译完全相同，且新坐标仍在父 CDS 内。",
              "- 正链延长 genomic end 3 bp；负链减小 genomic start 3 bp。涉及旧闭区间 protein_end 时只加 1。",
              "- 只替换原 GBK 字节中的坐标数字，未用 GenBank 重写器重排文件。PFAM_domain 保持原类型。", "",
              "repairs.tsv 的 genomic_start/end 使用 Biopython 的零基、右端不包含坐标；原 GBK 文本为一基闭区间。",
              "protein_start/end 的原值及修复值单独记录，不能把所有注释统一 end += 1。", "",
              "## 验证与待核查项", "",
              f"已核对 raw 全部 {raw_count} 个文件的修复前后 SHA256，与原始清单也一致。",
              "输出 GBK 已重新解析，逐 feature 验证仅发生列出的坐标变化；所有 DNA 和翻译保持一致。",
              "全部样本通过当前 BGC_annotation 蛋白与 DNA 解析器；未做长度过滤、去重或样本剔除。",
              f"另有 {len(review)} 条 domain 注释需人工核查，保留原值，详见 [review_items.tsv](review_items.tsv)。",
              "解析兼容性不表示每条生物学注释都已完备，也不表示模型推理已完成。", "",
              "## 文件与后续使用", "",
              "- Dataset_G_Planomonospora/GCF_assignments/ 与 Dataset_H_Reclassification/GCF_assignments/：完整修复副本。",
              "- 各数据集 curated_GCF_assignments.tsv：与 raw 逐字节一致，仍按原始 BGC 文件名匹配。",
              "- [samples.tsv](samples.tsv)：稳定 sample_id、原始 ID/标签、raw/processed hash 和逐样本验证。",
              "- [repairs.tsv](repairs.tsv)：逐 domain 修改和序列证据；[domain_checks.tsv](domain_checks.tsv)：全部 domain 检查。",
              "- [raw_integrity.tsv](raw_integrity.tsv)：原始文件保护记录；[metadata.json](metadata.json)：来源、版本和输出校验值。", "",
              "现有 pilot_features.py 与聚类入口仍限定 B/E；本次仅准备 G/H 修复数据，没有扩展这些入口或运行模型。",
              "后续提取应使用本目录 GBK，并将 processed_sha256 作为实际模型输入哈希，raw_source_sha256 保留为来源证据。", ""]
    return "\n".join(lines)


def run(args):
    workspace = args.workspace.expanduser().resolve()
    raw, inventory, output, encoder = [resolve(workspace, p) for p in
                                       (args.raw_root, args.inventory, args.output, args.encoder_project)]
    processed_root = (workspace / "data/processed").resolve()
    if (output.exists() or output == processed_root or not output.is_relative_to(processed_root)
            or output.is_relative_to(raw) or output.is_relative_to(inventory)):
        raise ValueError("Choose a new subdirectory of data/processed, separate from raw and inventory")
    if not raw.is_relative_to((workspace / "data/raw").resolve()):
        raise ValueError("Expected raw source under this workspace's data/raw")
    samples, before, inventory_hashes = load_inventory(workspace, raw, inventory)
    parsers, parser_provenance = load_parsers(encoder)
    protein_module = sys.modules[parsers["protein"].__module__]
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".gh-repair-", dir=output.parent) as temporary:
        stage = Path(temporary) / "result"
        stage.mkdir()
        results, repairs, checks = [], [], []
        for sample in samples:
            source = safe_child(raw, sample["archive_relpath"])
            raw_bytes = source.read_bytes()
            if file_sha256(source) != sample["source_sha256"]:
                raise ValueError(f"Raw sample changed during repair: {source}")
            original, raw_warnings = parsed_record(raw_bytes)
            changes, domain_checks = plan_repairs(original, protein_module)
            repaired_bytes = apply_byte_repairs(raw_bytes, original, changes)
            repaired, after_warnings = parsed_record(repaired_bytes)
            check_preserved(original, repaired, changes)
            # A second application must propose no additional changes.
            remaining, _ = plan_repairs(repaired, protein_module)
            if remaining:
                raise ValueError(f"Repair is not idempotent: {source}")
            dest = safe_child(stage, sample["archive_relpath"])
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(repaired_bytes)
            original_protein, pre_status, pre_error, pre_warnings = parse_protein(source, parsers["protein"])
            repaired_protein, post_status, post_error, post_warnings = parse_protein(dest, parsers["protein"])
            if post_status != "ok":
                raise ValueError(f"Processed protein validation failed: {sample['archive_relpath']}: {post_error}")
            original_dna, _ = parsers["dna"](source)
            repaired_dna, _ = parsers["dna"](dest)
            if original_dna != repaired_dna:
                raise ValueError(f"DNA changed: {source}")
            if original_protein and token_digest(original_protein) != token_digest(repaired_protein):
                raise ValueError(f"Previously valid protein token sequences changed: {source}")
            base = {"dataset": sample["dataset"], "sample_id": sample["sample_id"],
                    "original_bgc_id": sample["original_bgc_id"], "archive_relpath": sample["archive_relpath"]}
            repairs.extend({**base, **change} for change in changes)
            checks.extend({**base, **check} for check in domain_checks)
            results.append({
                **base, "gcf_raw": sample["gcf_raw"], "tsv_row_ids_json": sample["tsv_row_ids_json"],
                "raw_source_relpath": source.relative_to(workspace).as_posix(),
                "raw_source_sha256": sample["source_sha256"],
                "processed_relpath": sample["archive_relpath"], "processed_sha256": file_sha256(dest),
                "changed": raw_bytes != repaired_bytes, "repaired_domain_count": len(changes),
                "protein_end_edits": sum(c["new_protein_end"] is not None for c in changes),
                "record_id": original.id, "record_count": 1,
                "length_bp": len(original.seq), "dna_sha256": digest(original_dna),
                "dna_unchanged": True, "translations_unchanged": True, "unlisted_annotations_unchanged": True,
                "protein_parser_status_before": pre_status, "protein_parser_error_before": pre_error,
                "protein_parser_status_after": post_status, "dna_parser_status_before": "ok", "dna_parser_status_after": "ok",
                "protein_input_sha256_before": token_digest(original_protein),
                "protein_input_sha256_after": token_digest(repaired_protein),
                "protein_token_count": len(repaired_protein["tokens"]),
                "warnings_before_json": compact(raw_warnings + pre_warnings),
                "warnings_after_json": compact(after_warnings + post_warnings),
            })
        label_hashes = {}
        for dataset in DATASETS:
            rel = Path(dataset) / "curated_GCF_assignments.tsv"
            (stage / rel).write_bytes((raw / rel).read_bytes())
            label_hashes[rel.as_posix()] = file_sha256(stage / rel)
        after = snapshot(raw, workspace)
        if before != after:
            raise ValueError("Raw source changed during repair; processed results not published")
        for name, expected in parser_provenance["source_sha256"].items():
            if file_sha256(encoder / name) != expected:
                raise ValueError(f"Parser source changed during repair: {name}")
        for name, expected in inventory_hashes.items():
            if file_sha256(inventory / name) != expected:
                raise ValueError(f"Inventory changed during repair: {name}")
        dataset_rows = []
        for dataset in DATASETS:
            rows = [r for r in results if r["dataset"] == dataset]
            changes = [r for r in repairs if r["dataset"] == dataset]
            dataset_rows.append({"dataset": dataset, "sample_count": len(rows),
                                 "gcf_count": len({r["gcf_raw"] for r in rows}),
                                 "changed_files": sum(r["changed"] for r in rows),
                                 "repaired_domains": len(changes),
                                 "protein_end_edits": sum(r["new_protein_end"] is not None for r in changes),
                                 "protein_ok_before": sum(r["protein_parser_status_before"] == "ok" for r in rows),
                                 "protein_ok_after": len(rows), "dna_ok_after": len(rows)})
        review = [r for r in checks if r["status"] == "needs_review"]
        write_tsv(stage / "samples.tsv", results, ("sample_id", "dataset", "original_bgc_id", "gcf_raw"))
        write_tsv(stage / "repairs.tsv", repairs, ("sample_id", "dataset", "feature_index"))
        write_tsv(stage / "domain_checks.tsv", checks, ("sample_id", "dataset", "feature_index", "status"))
        write_tsv(stage / "review_items.tsv", review, ("sample_id", "dataset", "feature_index", "status", "reason"))
        write_tsv(stage / "dataset_summary.tsv", dataset_rows, ("dataset", "sample_count"))
        write_tsv(stage / "raw_integrity.tsv", [
            {"archive_relpath": r["archive_relpath"], "size_bytes": r["size_bytes"],
             "sha256_before": r["source_sha256"], "sha256_after": s["source_sha256"], "unchanged": True}
            for r, s in zip(before, after)], ("archive_relpath", "sha256_before", "sha256_after", "unchanged"))
        (stage / "summary.md").write_text(summary_text(dataset_rows, repairs, review, len(before)), encoding="utf-8")
        artifacts = {p.relative_to(stage).as_posix(): file_sha256(p) for p in sorted(stage.rglob("*")) if p.is_file()}
        write_json(stage / "metadata.json", {
            "format": "bgc_gh_coordinate_repairs_v1", "status": "complete", "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "raw_root": str(raw), "output": str(output), "inventory": str(inventory),
            "inventory_sha256": inventory_hashes, "raw_label_sha256": label_hashes,
            "counts": dataset_rows, "domain_checks": len(checks), "unresolved_review_items": len(review),
            "raw_file_count_verified": len(before), "raw_unchanged": True,
            "biopython_version": Bio.__version__, "python_version": sys.version.split()[0],
            "encoder_parsers": parser_provenance,
            "source_sha256": {name: file_sha256(Path(__file__).with_name(name)) for name in
                              ("repair_gh.py", "domain_coordinate_repair.py", "build_inventory.py", "pilot_features.py")},
            "policy": {"samples": "all G/H inventory entries; no filtering or deduplication",
                       "labels": "raw Unicode values and original labels copied byte-for-byte",
                       "changes": "sequence-proven domain C-terminal coordinates only; raw byte patches",
                       "raw_protection": "all raw files match inventory and pre/post snapshots",
                       "sequence_changes": False, "length_filter": None, "model_inference": False,
                       "identity": "retain original basenames and inventory sample IDs in manifest; distinct raw and processed hashes",
                       "publish": "only after all processed samples pass existing protein and DNA parsers; stage then rename"},
            "artifact_sha256": artifacts,
        })
        if output.exists():
            raise ValueError("Output appeared during repair; refusing to replace it")
        stage.rename(output)
    print(f"Published {len(results)} G/H GBKs with {len(repairs)} domain repairs to {output}")
    return output


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=ROOT)
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw/curated_datasets"))
    parser.add_argument("--inventory", type=Path, default=Path("data/processed/inventory_v1"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/repaired_gh_v1"))
    parser.add_argument("--encoder-project", type=Path, default=Path("../BGC_annotation"))
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run(args)
    except (ValueError, OSError, KeyError, ImportError) as exc:
        parser.exit(1, f"G/H repair failed: {exc}\n")


if __name__ == "__main__":
    main()
