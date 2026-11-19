#!/usr/bin/env python3
"""Restore existing Geneious aSDomain annotations in a new, audited H snapshot.

No domain discovery, boundary inference, sequence edits, deduplication or model
inference. Unresolvable annotations remain unchanged and block sample readiness.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
import tempfile

import Bio

from build_inventory import compact, digest, file_sha256, load_parsers, snapshot, write_json, write_tsv
from geneious_asdomain_restore import apply_reconstruction_bytes, plan_reconstruction
from pilot_features import read_tsv, safe_child
from repair_gh import load_inventory, parse_protein, parsed_record, token_digest

ROOT = Path(__file__).resolve().parents[1]
DATASET = "Dataset_H_Reclassification"
FORMAT = "bgc_h_asdomain_reconstruction_v1"


def verify_artifacts(root, hashes):
    if not hashes:
        raise ValueError("Empty source artifact manifest")
    for relative, expected in hashes.items():
        path = safe_child(root, relative)
        if not path.is_file() or file_sha256(path) != expected:
            raise ValueError(f"Source artifact mismatch: {relative}")


def validate_source(source, inventory_rows):
    metadata = json.loads((source / "metadata.json").read_text(encoding="utf-8"))
    if (metadata.get("format") != "bgc_gh_coordinate_repairs_v1"
            or metadata.get("status") != "complete" or metadata.get("raw_unchanged") is not True):
        raise ValueError("Expected a completed, integrity-checked G/H coordinate repair snapshot")
    verify_artifacts(source, metadata["artifact_sha256"])
    rows = [r for r in read_tsv(source / "samples.tsv") if r["dataset"] == DATASET]
    expected = {r["sample_id"]: r for r in inventory_rows if r["dataset"] == DATASET}
    if len(rows) != len(expected) or {r["sample_id"] for r in rows} != set(expected):
        raise ValueError("H sample membership differs from inventory")
    for row in rows:
        original = expected[row["sample_id"]]
        for name in ("archive_relpath", "original_bgc_id", "gcf_raw", "tsv_row_ids_json"):
            if row[name] != original[name]:
                raise ValueError(f"Source identity/label differs from inventory: {name}")
        if (row["raw_source_sha256"] != original["source_sha256"]
                or row["processed_relpath"] != original["archive_relpath"]
                or metadata["artifact_sha256"].get(row["processed_relpath"]) != row["processed_sha256"]):
            raise ValueError("Source raw/processed hash provenance mismatch")
        if (row["record_count"] != "1" or row["protein_parser_status_after"] != "ok"
                or row["dna_parser_status_after"] != "ok"
                or any(row[k] != "True" for k in (
                    "dna_unchanged", "translations_unchanged", "unlisted_annotations_unchanged"))):
            raise ValueError("Source sample has unresolved coordinate repair validation")
    label_rel = f"{DATASET}/curated_GCF_assignments.tsv"
    labels = read_tsv(source / label_rel)
    if sorted((r["BGC"], r["GCF"]) for r in labels) != sorted(
            (r["original_bgc_id"], r["gcf_raw"]) for r in rows):
        raise ValueError("Source H labels do not match samples")
    return sorted(rows, key=lambda r: r["archive_relpath"]), metadata


def tsv_rows(rows):
    """Store nested evidence as JSON, never Python repr."""
    return [{(key + "_json" if isinstance(value, (dict, list, tuple)) else key):
             compact(value) if isinstance(value, (dict, list, tuple)) else value
             for key, value in row.items()} for row in rows]


def token_stats(parsed):
    tokens = parsed["tokens"]
    return {"protein_token_count": len(tokens),
            "domain_token_count": sum(t["source"] == "domain" for t in tokens),
            "whole_protein_token_count": sum(t["source"] == "protein" for t in tokens),
            "max_token_aa": max(len(t["sequence"]) for t in tokens),
            "tokens_over_1021aa": sum(len(t["sequence"]) > 1021 for t in tokens)}


def check_token_transition(before, after, changes):
    """Old domains and unaffected proteins survive exactly; each restored domain is used."""
    before_tokens = {t["feature_index"]: t for t in before["tokens"]}
    after_tokens = {t["feature_index"]: t for t in after["tokens"]}
    parents = {c["parent_feature_index"] for c in changes}
    expected = (set(before_tokens) - parents) | {c["feature_index"] for c in changes}
    if set(after_tokens) != expected or len(after_tokens) != len(after["tokens"]):
        raise ValueError("Unexpected protein token coverage after reconstruction")
    for index, old in before_tokens.items():
        if index in parents:
            continue
        new = after_tokens[index]
        for key in ("source", "sequence", "protein_start", "protein_end", "parent_identifiers"):
            if new[key] != old[key]:
                raise ValueError(f"Unrelated existing token changed: feature {index}/{key}")
    for change in changes:
        new = after_tokens[change["feature_index"]]
        if (new["source"] != "domain" or new["sequence"] != change["translation"]
                or new["protein_start"] != change["protein_start"]
                or new["protein_end"] != change["protein_end"]):
            raise ValueError("Restored domain disagrees with current encoder input")


def summary_text(counts, review, inherited_count):
    lines = ["# H 现有 aSDomain 恢复", "",
             f"状态：**{counts['status']}**。保留全部 {counts['sample_count']} 个 H 样本和原 GCF 标签。",
             "本目录以 repaired_gh_v1 为输入；raw、坐标修复 v1 和旧特征/队列均未覆盖。", "",
             "| 检查 | 数量 |", "| --- | ---: |",
             f"| 存在 Geneious aSDomain 的样本 | {counts['samples_with_hidden_domains']} |",
             f"| 原隐藏 aSDomain | {counts['hidden_domains_before']} |",
             f"| 成功恢复 | {counts['restored_domains']} |",
             f"| 待复核、保持原注释 | {counts['unresolved_asdomain_annotations']} |",
             f"| 修正旧父 locus_tag 的恢复条目 | {counts['parent_id_corrections']} |",
             f"| 保留重复位置组 / 额外条目 | {counts['duplicate_groups']} / {counts['duplicate_extra_entries']} |",
             f"| 现有蛋白 / DNA 解析通过 | {counts['protein_parser_ok']} / {counts['dna_parser_ok']} |",
             f"| aSDomain 恢复检查通过的样本 | {counts['asdomain_ready_samples']} |",
             f"| >1021 aa 输入单元（前 → 后） | {counts['long_tokens_before']} → {counts['long_tokens_after']} |", "",
             "只恢复带有精确 `Geneious type: aSDomain` 标记的 misc_feature；要求唯一同链父 CDS、整密码子区间、",
             "基因组翻译与父 CDS 翻译切片完全相等。保留原基因组边界和重复条目，补齐 translation、",
             "protein_start/end，必要时显式修正父 locus_tag。所有 DNA、CDS 翻译、原生 aSDomain、",
             "PFAM_domain 和其他非目标 feature 保持不变。没有新做 Pfam 切分或模型推理。", "",
             "protein_start/end 和日志中的 genomic_start/end 均为零基、右端不包含坐标。",
             "这种恢复证明原边界与当前序列内部一致，不独立证明原始 HMM 边界的生物学正确性。", ""]
    if review:
        lines += ["## 仍待处理", "",
                  "以下样本保留在清单内，但标记为 `needs_review`；整个 H 尚不能称为完整恢复的基准输入。",
                  "解析通过不等于注释完整，不能自动删掉这些样本或混入旧 H 特征后比较性能。", ""]
        for row in review:
            lines.append(f"- {row['original_bgc_id']} feature {row['feature_index']}: {row.get('reason', '')}")
        lines += ["", "原始 COMMENT（含人工编辑记录）保存在 review_items.tsv 的 record_comment / record_structured_comment_json 列。", ""]
        if any(r["original_bgc_id"] == "UOC915" for r in review):
            lines += ["UOC915 的原始编辑记录注明在 17584–17585 插入 GG，将小基因合并成 NRPS1gene。",
                      "未恢复的两条重复注释横跨该插入，其 1154 bp 区间及相对父 CDS 的 3115 nt 起点均不整帧；",
                      "当前父 CDS 翻译与 DNA 一致，不能靠边界取整恢复这两条旧注释。", ""]
    lines += [f"坐标修复阶段遗留的 {inherited_count} 条 H 复核项另存 inherited_review_items.tsv；",
              "它们不计入本次隐藏 aSDomain 数，且本次未修改。", "",
              "## 可追溯输出", "",
              "- samples.tsv：样本身份、标签、raw/source/processed 哈希、前后输入统计和逐样本状态。",
              "- reconstructions.tsv / domain_checks.tsv：逐条变化、坐标、序列与父 CDS 校验证据。",
              "- parent_id_mappings.tsv / duplicate_annotations.tsv：父 ID 修正和保留的重复注释。",
              "- review_items.tsv：没有足够恢复证据的注释；不做边界取整或猜测。",
              "- source_provenance/：输入快照的元数据、修复记录和清单；metadata.json：所有产物哈希。",
              "- raw_integrity.tsv / source_integrity.tsv：输入内容未改变的校验记录。", "",
              "旧 pilot_gh_v1 的 H 尚未采用此次恢复，不能用于新的 H 基准；本次未创建模型运行队列。",
              "现有 prepare-repaired 仅接受旧坐标修复 schema，本恢复快照不可冒充该 schema。", ""]
    return "\n".join(lines)


def reconstruct(args):
    workspace = Path(args.workspace).expanduser().resolve()
    resolve = lambda path: (workspace / path).expanduser().resolve()
    source, output = resolve(args.source), resolve(args.output)
    raw, inventory, encoder = resolve(args.raw), resolve(args.inventory), resolve(args.encoder_project)
    processed = workspace / "data/processed"
    if (not output.is_relative_to(processed) or output == processed
            or output.is_relative_to(source) or source.is_relative_to(output)
            or output.exists()):
        raise ValueError("Output must be a new processed directory outside the input snapshot")
    if not raw.is_relative_to(workspace / "data/raw"):
        raise ValueError("Raw input must remain under workspace data/raw")
    inventory_rows, raw_before, inventory_hashes = load_inventory(workspace, raw, inventory)
    samples, _ = validate_source(source, inventory_rows)
    source_before = snapshot(source, workspace)
    parsers, parser_provenance = load_parsers(encoder)
    code_hashes = {name: file_sha256(Path(__file__).parent / name) for name in (
        "reconstruct_h.py", "geneious_asdomain_restore.py", "repair_gh.py", "build_inventory.py", "pilot_features.py")}
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".h-asdomain-", dir=output.parent) as temporary:
        stage = Path(temporary) / "result"
        stage.mkdir()
        results, reconstructions, checks = [], [], []
        for sample in samples:
            original_path = safe_child(source, sample["processed_relpath"])
            original_bytes = original_path.read_bytes()
            original, original_warnings = parsed_record(original_bytes)
            changes, domain_checks = plan_reconstruction(original)
            restored_bytes = apply_reconstruction_bytes(original_bytes, original, changes)
            restored, restored_warnings = parsed_record(restored_bytes)
            remaining, remaining_checks = plan_reconstruction(restored)
            if remaining or len(remaining_checks) != len(domain_checks) - len(changes):
                raise ValueError(f"Reconstruction is not idempotent: {original_path}")
            dest = safe_child(stage, sample["archive_relpath"])
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(restored_bytes)
            before, before_status, before_error, bw = parse_protein(original_path, parsers["protein"])
            after, after_status, after_error, aw = parse_protein(dest, parsers["protein"])
            if before_status != "ok" or after_status != "ok":
                raise ValueError(f"Protein parser failed: {sample['original_bgc_id']}: {before_error}; {after_error}")
            if token_digest(before) != sample["protein_input_sha256_after"]:
                raise ValueError("Source protein input differs from coordinate repair evidence")
            check_token_transition(before, after, changes)
            before_dna, _ = parsers["dna"](original_path)
            after_dna, _ = parsers["dna"](dest)
            if before_dna != after_dna or digest(after_dna) != sample["dna_sha256"]:
                raise ValueError("DNA changed relative to source evidence")
            base = {key: sample[key] for key in (
                "sample_id", "dataset", "original_bgc_id", "archive_relpath")}
            reconstructions.extend({**base, **change} for change in changes)
            for check in domain_checks:
                row = {**base, **check}
                if check["status"] == "needs_review":
                    row["record_comment"] = original.annotations.get("comment", "")
                    row["record_structured_comment"] = original.annotations.get("structured_comment", {})
                checks.append(row)
            unresolved = sum(c["status"] == "needs_review" for c in domain_checks)
            results.append({**base, "gcf_raw": sample["gcf_raw"], "tsv_row_ids_json": sample["tsv_row_ids_json"],
                "raw_source_relpath": sample["raw_source_relpath"], "raw_source_sha256": sample["raw_source_sha256"],
                "source_relpath": sample["processed_relpath"], "source_sha256": sample["processed_sha256"],
                "processed_relpath": sample["archive_relpath"], "processed_sha256": file_sha256(dest),
                "changed": original_bytes != restored_bytes, "hidden_asdomain_count_before": len(domain_checks),
                "restored_domain_count": len(changes), "unresolved_asdomain_count": unresolved,
                "asdomain_status": "needs_review" if unresolved else "complete",
                "ready_for_encoding": not unresolved, "record_count": 1, "record_id": original.id,
                "length_bp": len(original.seq), "dna_sha256": digest(after_dna),
                "dna_unchanged": True, "cds_translations_unchanged": True, "non_target_features_unchanged": True,
                "protein_parser_status_after": after_status, "dna_parser_status_after": "ok",
                "protein_input_sha256_before": token_digest(before), "protein_input_sha256_after": token_digest(after),
                **{key + "_before": value for key, value in token_stats(before).items()},
                **{key + "_after": value for key, value in token_stats(after).items()},
                "warnings_before_json": compact(original_warnings + bw),
                "warnings_after_json": compact(restored_warnings + aw)})
        label_rel = f"{DATASET}/curated_GCF_assignments.tsv"
        shutil.copyfile(source / label_rel, stage / label_rel)
        provenance = stage / "source_provenance"
        provenance.mkdir()
        for name in ("metadata.json", "samples.tsv", "repairs.tsv", "domain_checks.tsv", "review_items.tsv"):
            shutil.copyfile(source / name, provenance / name)
        inherited = [r for r in read_tsv(source / "review_items.tsv") if r["dataset"] == DATASET]
        review = [r for r in checks if r["status"] == "needs_review"]
        mappings = [r for r in reconstructions if r["old_qualifiers"].get("locus_tag") != r["new_qualifiers"].get("locus_tag")]
        duplicates = [r for r in checks if r.get("duplicate_group")]
        groups = Counter((r["sample_id"], r["duplicate_group"]) for r in duplicates)
        counts = {"sample_count": len(results), "gcf_count": len({r["gcf_raw"] for r in results}),
                  "changed_files": sum(r["changed"] for r in results),
                  "samples_with_hidden_domains": sum(r["hidden_asdomain_count_before"] > 0 for r in results),
                  "hidden_domains_before": len(checks), "restored_domains": len(reconstructions),
                  "unresolved_asdomain_annotations": len(review), "parent_id_corrections": len(mappings),
                  "duplicate_groups": len(groups), "duplicate_extra_entries": sum(n - 1 for n in groups.values()),
                  "protein_parser_ok": len(results), "dna_parser_ok": len(results),
                  "asdomain_ready_samples": sum(r["ready_for_encoding"] for r in results),
                  "long_tokens_before": sum(r["tokens_over_1021aa_before"] for r in results),
                  "long_tokens_after": sum(r["tokens_over_1021aa_after"] for r in results),
                  "status": "needs_review" if review else "complete"}
        for name, rows, leading in (
                ("samples", results, ("sample_id", "dataset", "original_bgc_id", "gcf_raw", "asdomain_status")),
                ("reconstructions", reconstructions, ("sample_id", "original_bgc_id", "feature_index")),
                ("domain_checks", checks, ("sample_id", "original_bgc_id", "feature_index", "status", "reason")),
                ("review_items", review, ("sample_id", "original_bgc_id", "feature_index", "status", "reason")),
                ("parent_id_mappings", mappings, ("sample_id", "original_bgc_id", "feature_index")),
                ("duplicate_annotations", duplicates, ("sample_id", "original_bgc_id", "feature_index", "duplicate_group")),
                ("inherited_review_items", inherited, ("sample_id", "dataset", "feature_index", "reason"))):
            write_tsv(stage / f"{name}.tsv", tsv_rows(rows), leading)
        if snapshot(raw, workspace) != raw_before or snapshot(source, workspace) != source_before:
            raise ValueError("Raw or coordinate-repaired input changed during reconstruction")
        for root, hashes in ((inventory, inventory_hashes), (encoder, parser_provenance["source_sha256"]),
                             (Path(__file__).parent, code_hashes)):
            verify_artifacts(root, hashes)
        for name, rows in (("raw_integrity", raw_before), ("source_integrity", source_before)):
            write_tsv(stage / f"{name}.tsv", [{**r, "unchanged": True} for r in rows], ("archive_relpath",))
        (stage / "summary.md").write_text(summary_text(counts, review, len(inherited)), encoding="utf-8")
        write_json(stage / "metadata.json", {
            "format": FORMAT, "status": counts["status"], "snapshot_complete": True,
            "ready_for_full_h_encoding": not review, "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "counts": counts, "source": str(source), "output": str(output), "raw_root": str(raw),
            "raw_unchanged": True, "source_unchanged": True, "raw_file_count_verified": len(raw_before),
            "source_metadata_sha256": file_sha256(source / "metadata.json"),
            "inventory_sha256": inventory_hashes, "encoder_parsers": parser_provenance,
            "source_sha256": code_hashes, "biopython_version": Bio.__version__, "python_version": sys.version.split()[0],
            "inherited_review_items": len(inherited),
            "policy": {"existing_asdomain_only": True, "pfam_splitting": False, "boundary_changes": False,
                       "sequence_changes": False, "sample_filtering": False, "annotation_deduplication": False,
                       "model_inference": False, "unresolved": "preserve annotation and mark sample needs_review",
                       "coordinates": "zero-based half-open", "publish": "new snapshot; readiness reported separately"},
            "artifact_sha256": {p.relative_to(stage).as_posix(): file_sha256(p)
                                for p in sorted(stage.rglob("*")) if p.is_file()}})
        if output.exists():
            raise ValueError("Output appeared during reconstruction; refusing overwrite")
        stage.rename(output)
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", default=str(ROOT))
    parser.add_argument("--source", default="data/processed/repaired_gh_v1")
    parser.add_argument("--output", default="data/processed/reconstructed_h_v1")
    parser.add_argument("--raw", default="data/raw/curated_datasets")
    parser.add_argument("--inventory", default="data/processed/inventory_v1")
    parser.add_argument("--encoder-project", default="../BGC_annotation")
    args = parser.parse_args(argv)
    try:
        counts = reconstruct(args)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"Reconstruction failed: {exc}\n")
    print(json.dumps(counts, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
