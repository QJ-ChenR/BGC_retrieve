"""Adapt the audited G/H repair snapshot to portable, stable-ID encoder inputs.

This module is called by pilot_features with that module injected as ``pilot``.
Only CPU input parsers run here; checkpoints and GPU inference are not loaded.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import warnings

REPAIR_SCHEMA = "bgc_gh_coordinate_repairs_v1"
PROVENANCE_FILES = ("metadata.json", "samples.tsv", "repairs.tsv", "review_items.tsv")
EXTRA_COHORT_FIELDS = ("raw_source_sha256", "processed_source_relpath", "repair_count")
CHECK_FIELDS = (
    "dataset", "sample_id", "source_sha256", "protein_parser_status", "dna_parser_status",
    "length_bp", "n_base_count", "ambiguous_base_count", "cds_count", "domain_token_count",
    "whole_protein_token_count", "protein_token_count", "max_token_residues",
    "tokens_over_1021_residues", "esm_1021_window_count", "dna_sha256",
    "protein_input_sha256", "warnings_json",
)


def _compact(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _digest(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _check_hashes(root, hashes, pilot, description):
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError(f"Missing {description} hashes")
    for relative, expected in hashes.items():
        path = pilot.safe_child(root, relative)
        if not path.is_file() or pilot.sha256(path) != expected:
            raise ValueError(f"{description} checksum mismatch: {relative}")


def _repair_metadata(metadata):
    if metadata.get("format") != REPAIR_SCHEMA or metadata.get("status") != "complete":
        raise ValueError("Repair snapshot is incomplete or has the wrong format")
    if metadata.get("raw_unchanged") is not True:
        raise ValueError("Repair snapshot lacks raw-preservation verification")
    if not metadata.get("raw_root") or not isinstance(metadata.get("artifact_sha256"), dict):
        raise ValueError("Repair snapshot lacks source/artifact provenance")


def _selected_rows(root, metadata, dataset_names, pilot):
    """Validate identity/labels inside a repair snapshot, or its embedded copy."""
    _repair_metadata(metadata)
    samples = pilot.read_tsv(root / "samples.tsv")
    if len(samples) != len({r["sample_id"] for r in samples}):
        raise ValueError("Duplicate sample ID in repair provenance")
    selected = [r for r in samples if r["dataset"] in dataset_names]
    pilot.validate_rows(selected, dataset_names)
    counts = {r["dataset"]: r["sample_count"] for r in metadata["counts"]}
    repairs = pilot.read_tsv(root / "repairs.tsv")
    repair_counts = Counter(r["sample_id"] for r in repairs)
    sample_ids = {r["sample_id"] for r in samples}
    if set(repair_counts) - sample_ids:
        raise ValueError("Repair log contains unknown sample IDs")
    for row in selected:
        if (row["record_count"] != "1" or row["protein_parser_status_after"] != "ok"
                or row["dna_parser_status_after"] != "ok"
                or any(row[field] != "True" for field in
                       ("dna_unchanged", "translations_unchanged", "unlisted_annotations_unchanged"))):
            raise ValueError(f"Repaired sample no longer meets audited policy: {row['sample_id']}")
        if row["processed_relpath"] != row["archive_relpath"]:
            raise ValueError(f"Repaired path disagrees with archive identity: {row['sample_id']}")
        for field in ("raw_source_sha256", "processed_sha256", "protein_input_sha256_after", "dna_sha256"):
            value = row[field]
            if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                raise ValueError(f"Invalid {field}: {row['sample_id']}")
        if int(row["repaired_domain_count"]) != repair_counts[row["sample_id"]]:
            raise ValueError(f"Repair count disagrees with log: {row['sample_id']}")
        if metadata["artifact_sha256"].get(row["processed_relpath"]) != row["processed_sha256"]:
            raise ValueError(f"Processed hash disagrees with repair manifest: {row['sample_id']}")
    for dataset in dataset_names:
        subset = [r for r in selected if r["dataset"] == dataset]
        if counts.get(dataset) != len(subset):
            raise ValueError(f"Repair metadata sample count mismatch: {dataset}")
        rel = dataset + "/curated_GCF_assignments.tsv"
        path = pilot.safe_child(root, rel)
        expected = metadata["raw_label_sha256"].get(rel)
        if (not path.is_file() or pilot.sha256(path) != expected
                or metadata["artifact_sha256"].get(rel) != expected):
            raise ValueError(f"Original labels disagree with repair hashes: {dataset}")
        labels = pilot.read_tsv(path)
        if sorted((r["BGC"], r["GCF"]) for r in labels) != sorted(
                (r["original_bgc_id"], r["gcf_raw"]) for r in subset):
            raise ValueError(f"Original labels disagree with repaired samples: {dataset}")
    return sorted(selected, key=lambda r: (r["dataset"], r["sample_id"]))


def _cohort_row(sample, pilot):
    row = {key: sample[key] for key in pilot.COHORT_FIELDS if key != "source_sha256"}
    return {**row, "source_sha256": sample["processed_sha256"],
            "raw_source_sha256": sample["raw_source_sha256"],
            "processed_source_relpath": sample["processed_relpath"],
            "repair_count": sample["repaired_domain_count"]}


def _input_check(path, sample, parsers, pilot):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        protein = parsers["protein"](path)
        dna, dna_metadata = parsers["dna"](path)
    tokens = protein["tokens"]
    sequences = [t["sequence"] for t in tokens]
    token_hash = _digest(_compact(sequences))
    if not sequences or any(not s for s in sequences):
        raise ValueError(f"No usable protein tokens: {sample['sample_id']}")
    if (protein["source_sha256"] != sample["processed_sha256"]
            or dna_metadata["source_sha256"] != sample["processed_sha256"]
            or pilot.sha256(path) != sample["processed_sha256"]):
        raise ValueError(f"Parser input bytes changed: {sample['sample_id']}")
    if (protein["bgc_id"] != sample["sample_id"] or dna_metadata["bgc_id"] != sample["sample_id"]):
        raise ValueError(f"Parser changed stable sample ID: {sample['sample_id']}")
    if token_hash != sample["protein_input_sha256_after"]:
        raise ValueError(f"Protein token sequences changed since repair: {sample['sample_id']}")
    if len(tokens) != int(sample["protein_token_count"]):
        raise ValueError(f"Protein token count changed since repair: {sample['sample_id']}")
    if len(dna) != int(sample["length_bp"]) or _digest(dna) != sample["dna_sha256"]:
        raise ValueError(f"DNA sequence changed since repair: {sample['sample_id']}")
    if {t["source"] for t in tokens} - {"domain", "protein"}:
        raise ValueError(f"Unexpected protein token source: {sample['sample_id']}")
    return {
        "dataset": sample["dataset"], "sample_id": sample["sample_id"],
        "source_sha256": sample["processed_sha256"], "protein_parser_status": "ok",
        "dna_parser_status": "ok", "length_bp": len(dna), "n_base_count": dna.count("N"),
        "ambiguous_base_count": sum(c not in "ACGT" for c in dna),
        "cds_count": len(protein["parents"]),
        "domain_token_count": sum(t["source"] == "domain" for t in tokens),
        "whole_protein_token_count": sum(t["source"] == "protein" for t in tokens),
        "protein_token_count": len(tokens), "max_token_residues": max(map(len, sequences)),
        "tokens_over_1021_residues": sum(len(s) > 1021 for s in sequences),
        "esm_1021_window_count": sum((len(s) + 1020) // 1021 for s in sequences),
        "dna_sha256": _digest(dna), "protein_input_sha256": token_hash,
        "warnings_json": _compact([{"category": w.category.__name__, "message": str(w.message)}
                                    for w in caught]),
    }


def validate_repaired_cohort(prepared, metadata, rows, pilot):
    """Validate embedded repair provenance without accessing old absolute paths."""
    prepared = Path(prepared)
    provenance = prepared / "repair_provenance"
    dataset_names = metadata["datasets"]
    if set(dataset_names) - {pilot.DATASETS[k][0] for k in ("G", "H") if k in pilot.DATASETS}:
        raise ValueError("Repaired pilot only supports G/H")
    required = {"cohort.tsv", "input_checks.tsv", "summary.md"}
    required.update("repair_provenance/" + name for name in PROVENANCE_FILES)
    for dataset in dataset_names:
        required.update((dataset + "/curated_GCF_assignments.tsv",
                         "repair_provenance/" + dataset + "/curated_GCF_assignments.tsv"))
    if required - set(metadata["artifact_sha256"]):
        raise ValueError("Prepared manifest lacks required repair provenance artifacts")
    _check_hashes(prepared, metadata["artifact_sha256"], pilot, "Prepared artifact")
    repair_meta = pilot.read_json(provenance / "metadata.json")
    if pilot.sha256(provenance / "metadata.json") != metadata["repair_metadata_sha256"]:
        raise ValueError("Embedded repair metadata checksum mismatch")
    if metadata["raw_root"] != repair_meta["raw_root"]:
        raise ValueError("Raw provenance root changed")
    for name in PROVENANCE_FILES[1:]:
        if pilot.sha256(provenance / name) != repair_meta["artifact_sha256"].get(name):
            raise ValueError(f"Embedded provenance disagrees with original repair hashes: {name}")
    selected = _selected_rows(provenance, repair_meta, dataset_names, pilot)
    expected = [_cohort_row(r, pilot) for r in selected]
    keys = (*pilot.COHORT_FIELDS, *EXTRA_COHORT_FIELDS)
    canonical = lambda items: sorted(tuple(r[key] for key in keys) for r in items)
    if canonical(rows) != canonical(expected):
        raise ValueError("Prepared cohort disagrees with embedded repair provenance")
    if metadata["counts"] != dict(Counter(r["dataset"] for r in rows)):
        raise ValueError("Prepared metadata counts disagree with cohort")
    checks = pilot.read_tsv(prepared / "input_checks.tsv")
    if (len(checks) != len(rows) or len({c["sample_id"] for c in checks}) != len(rows)
            or {c["sample_id"] for c in checks} != {r["sample_id"] for r in rows}):
        raise ValueError("Input checks do not cover the full cohort")
    source_by_id = {r["sample_id"]: r for r in selected}
    for check in checks:
        sample = source_by_id[check["sample_id"]]
        if (check["dataset"] != sample["dataset"] or check["protein_parser_status"] != "ok"
                or check["dna_parser_status"] != "ok" or check["source_sha256"] != sample["processed_sha256"]
                or check["protein_input_sha256"] != sample["protein_input_sha256_after"]
                or check["dna_sha256"] != sample["dna_sha256"]
                or check["length_bp"] != sample["length_bp"]
                or check["protein_token_count"] != sample["protein_token_count"]):
            raise ValueError(f"Prepared parser check disagrees with repair: {check['sample_id']}")
    for dataset in dataset_names:
        actual = pilot.read_tsv(prepared / dataset / "curated_GCF_assignments.tsv")
        if sorted((r["BGC"], r["GCF"]) for r in actual) != sorted(
                (r["sample_id"], r["gcf_raw"]) for r in rows if r["dataset"] == dataset):
            raise ValueError(f"Prepared labels disagree with cohort: {dataset}")


def _summary(rows, checks, review, dataset_names):
    lines = ["# G/H encoder 输入准备", "", "全部样本使用修复后的 GBK 字节，稳定 sample_id 作为文件名；原始 ID、GCF 标签及 raw/processed 哈希保留在 cohort.tsv。",
             "", "| 数据集 | BGC | GCF | 修改过的 GBK | 蛋白解析通过 | DNA 解析通过 | 最长 DNA/bp | 最长蛋白 token/aa | ESM 1021 窗口总数 |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for dataset in dataset_names:
        subset = [r for r in rows if r["dataset"] == dataset]
        checked = [c for c in checks if c["dataset"] == dataset]
        lines.append(f"| {dataset} | {len(subset)} | {len({r['gcf_raw'] for r in subset})} | "
                     f"{sum(int(r['repair_count']) > 0 for r in subset)} | {len(checked)} | {len(checked)} | "
                     f"{max(c['length_bp'] for c in checked)} | {max(c['max_token_residues'] for c in checked)} | "
                     f"{sum(c['esm_1021_window_count'] for c in checked)} |")
    lines.extend(["", "本次只执行现有 CPU 蛋白/DNA 解析器：蛋白 token 序列及完整 DNA 与修复快照一致。",
                  "ESM2 与 BGC-MAC 共用蛋白 token；超过 1021 aa 的 token 由现有 ESM 接口分窗并按残基数加权。",
                  "DNA 使用完整原始序列，不做长度过滤、截断或自动改为 windowed；未加载 tokenizer/checkpoint，DNA token 数和显存需求需服务器推理时记录。",
                  f"本批保留 {len(review)} 条待人工复核注释，见 repair_provenance/review_items.tsv；无样本剔除。",
                  "未运行 ESM2、BGC-MAC、DNA encoder 推理或聚类，解析通过不代表 GPU 推理已通过。", "",
                  "模型输入与全部必要来源清单都在本目录；复制整个目录到服务器即可校验，不依赖旧 raw/repaired 绝对路径。",
                  "source_sha256 对应模型实际读取的修复文件；raw_source_sha256 仅记录原始备份来源。", ""])
    return "\n".join(lines)


def prepare_repaired(args, pilot):
    # Extraction/validation can load a portable cohort with only stdlib imports.
    # Biopython and the local encoder parsers are needed only during preparation.
    import build_inventory

    source, output, encoder = [Path(p).expanduser().resolve() for p in
                               (args.repaired_root, args.output, args.encoder_project)]
    metadata_path = source / "metadata.json"
    source_metadata_hash = pilot.sha256(metadata_path)
    source_metadata = pilot.read_json(metadata_path)
    _repair_metadata(source_metadata)
    raw = Path(source_metadata["raw_root"]).expanduser().resolve()
    if (output.exists() or output.is_relative_to(source) or output.is_relative_to(raw)
            or source.is_relative_to(output) or raw.is_relative_to(output)):
        raise ValueError("Choose a new output directory separate from repaired inputs and raw data")
    if not args.datasets or set(args.datasets) - {"G", "H"}:
        raise ValueError("prepare-repaired only supports G/H")
    dataset_names = [pilot.DATASETS[key][0] for key in args.datasets]
    _check_hashes(source, source_metadata["artifact_sha256"], pilot, "Repair artifact")
    required_source = set(PROVENANCE_FILES[1:])
    if required_source - set(source_metadata["artifact_sha256"]):
        raise ValueError("Repair snapshot lacks required provenance artifacts")
    samples = _selected_rows(source, source_metadata, dataset_names, pilot)
    for dataset in dataset_names:
        actual = {p.relative_to(source).as_posix()
                  for p in (source / dataset / "GCF_assignments").rglob("*")
                  if p.is_file() and p.suffix.lower() in {".gbk", ".gb", ".genbank"}}
        expected = {r["processed_relpath"] for r in samples if r["dataset"] == dataset}
        if actual != expected:
            raise ValueError(f"Repaired GBK file set differs from sample manifest: {dataset}")
    rows = [_cohort_row(r, pilot) for r in samples]
    pilot.validate_rows(rows, dataset_names)
    parsers, parser_provenance = build_inventory.load_parsers(encoder)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".pilot-gh-", dir=output.parent) as temporary:
        stage = Path(temporary) / "prepared"
        stage.mkdir()
        checks = []
        for sample in samples:
            original = pilot.safe_child(source, sample["processed_relpath"])
            dest = stage / sample["dataset"] / "inputs" / (sample["sample_id"] + ".gbk")
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original, dest)
            if pilot.sha256(dest) != sample["processed_sha256"]:
                raise ValueError(f"Processed source changed during copy: {sample['sample_id']}")
            checks.append(_input_check(dest, sample, parsers, pilot))
        pilot.write_tsv(stage / "cohort.tsv", rows, (*pilot.COHORT_FIELDS, *EXTRA_COHORT_FIELDS))
        pilot.write_tsv(stage / "input_checks.tsv", checks, CHECK_FIELDS)
        provenance_files = list(PROVENANCE_FILES)
        for dataset in dataset_names:
            label_rel = dataset + "/curated_GCF_assignments.tsv"
            provenance_files.append(label_rel)
            pilot.write_tsv(stage / label_rel,
                            [{"BGC": r["sample_id"], "GCF": r["gcf_raw"]}
                             for r in rows if r["dataset"] == dataset], ("BGC", "GCF"))
        for relative in provenance_files:
            dest = stage / "repair_provenance" / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(pilot.safe_child(source, relative), dest)
        review = [r for r in pilot.read_tsv(source / "review_items.tsv") if r["dataset"] in dataset_names]
        (stage / "summary.md").write_text(_summary(rows, checks, review, dataset_names), encoding="utf-8")
        artifacts = {p.relative_to(stage).as_posix(): pilot.sha256(p)
                     for p in sorted(stage.rglob("*")) if p.is_file()}
        prepared_metadata = {
            "format": pilot.REPAIRED_SCHEMA, "status": "complete", "created_at_utc": pilot.now(),
            "datasets": dataset_names, "counts": dict(Counter(r["dataset"] for r in rows)),
            "raw_root": source_metadata["raw_root"], "source_root": str(source),
            "repair_metadata_sha256": source_metadata_hash, "encoder_parsers": parser_provenance,
            "artifact_sha256": artifacts,
            "script_sha256": {"prepare_repaired_features.py": pilot.sha256(__file__),
                              "pilot_features.py": pilot.sha256(Path(pilot.__file__)),
                              "build_inventory.py": pilot.sha256(Path(build_inventory.__file__))},
            "policy": {"identity": "inventory stable sample_id", "labels": "unchanged raw Unicode",
                       "source_sha256": "actual processed GBK bytes", "raw_source_sha256": "original backup provenance",
                       "length_filter": None, "deduplicated": False, "sample_exclusions": [],
                       "protein_window_residues": 1021, "dna_sequence_mode": "full", "model_inference": False},
        }
        pilot.write_json(stage / "metadata.json", prepared_metadata)
        validate_repaired_cohort(stage, prepared_metadata, rows, pilot)
        _check_hashes(source, source_metadata["artifact_sha256"], pilot, "Repair source after preparation")
        if pilot.sha256(metadata_path) != source_metadata_hash:
            raise ValueError("Repair metadata changed during preparation")
        _check_hashes(encoder, parser_provenance["source_sha256"], pilot, "Encoder parser source")
        if output.exists():
            raise ValueError("Output appeared during preparation; refusing to replace it")
        stage.rename(output)
    print(f"Prepared {len(rows)} repaired G/H BGCs: {output}")
    return output
