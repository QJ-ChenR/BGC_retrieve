"""Portable encoder cohorts from reconstructed H plus optional repaired G.

Preparation runs CPU parsers only. A complete queue may retain unresolved
aSDomain annotations: their sample IDs and counts remain explicit in every
manifest, and no sample is silently removed to obtain a complete annotation set.
"""

from collections import Counter
import json
from pathlib import Path
import shutil
import tempfile

from prepare_repaired_features import (
    CHECK_FIELDS, PROVENANCE_FILES as REPAIR_FILES, _check_hashes, _compact,
    _digest, _input_check, _repair_metadata, _selected_rows,
)

H = "Dataset_H_Reclassification"
G = "Dataset_G_Planomonospora"
SOURCE_SCHEMA = "bgc_h_asdomain_reconstruction_v1"
OLD_FILES = ("metadata.json", "samples.tsv", "repairs.tsv", "domain_checks.tsv", "review_items.tsv")
H_FILES = ("metadata.json", "samples.tsv", "reconstructions.tsv", "domain_checks.tsv",
           "review_items.tsv", "parent_id_mappings.tsv", "duplicate_annotations.tsv",
           "inherited_review_items.tsv", H + "/curated_GCF_assignments.tsv",
           *("source_provenance/" + name for name in OLD_FILES))
EXTRA_COHORT_FIELDS = (
    "raw_source_sha256", "processed_source_relpath", "source_kind", "coordinate_source_sha256",
    "repair_count", "restored_domain_count", "unresolved_asdomain_count", "inherited_review_count",
    "annotation_status", "annotation_complete", "review_policy",
)
REVIEW_FIELDS = ("annotation_status", "annotation_complete", "unresolved_asdomain_count",
                 "unresolved_sample_ids", "review_policy")
INPUT_CHECK_FIELDS = (*CHECK_FIELDS, "annotation_status", "annotation_complete",
                      "unresolved_asdomain_count", "review_policy")


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha(value):
    return isinstance(value, str) and len(value) == 64 and not set(value) - set("0123456789abcdef")


def _read_embedded(root, names, pilot):
    metadata = pilot.read_json(root / "metadata.json")
    hashes = metadata.get("artifact_sha256", {})
    _require(isinstance(hashes, dict) and hashes, "Missing source artifact hashes")
    required = set(names) - {"metadata.json"}
    _require(not required - set(hashes), "Source manifest lacks required provenance artifacts")
    _check_hashes(root, {name: hashes[name] for name in required}, pilot, "Embedded source provenance")
    return metadata


def _identity(row, sample):
    for key in ("sample_id", "dataset", "original_bgc_id", "archive_relpath"):
        _require(row.get(key) == sample.get(key), f"Annotation log identity differs: {key}")


def _entries(rows, samples, description):
    result = {}
    for row in rows:
        _require(row.get("sample_id") in samples, f"Unknown sample in {description}")
        _identity(row, samples[row["sample_id"]])
        index = int(row["feature_index"])
        _require(index >= 0, f"Invalid feature index in {description}")
        key = row["sample_id"], index
        _require(key not in result, f"Duplicate annotation in {description}")
        result[key] = row
    return result


def _same_log_subset(actual, expected, description):
    _require(set(actual) == set(expected), f"{description} membership differs from domain checks")
    for key, row in actual.items():
        # Tables have different unions of columns; an absent value is blank.
        keys = set(row) | set(expected[key])
        _require(all(row.get(field, "") == expected[key].get(field, "") for field in keys),
                 f"{description} evidence differs from domain checks")


def _h_rows(root, pilot):
    """Validate H's reconstruction and its embedded coordinate-repair ancestry."""
    metadata = _read_embedded(root, H_FILES, pilot)
    _require(metadata.get("format") == SOURCE_SCHEMA and metadata.get("snapshot_complete") is True
             and metadata.get("raw_unchanged") is True and metadata.get("source_unchanged") is True,
             "Expected a complete, integrity-checked H reconstruction snapshot")
    policy = metadata.get("policy", {})
    _require(policy.get("existing_asdomain_only") is True and all(policy.get(key) is False for key in
             ("pfam_splitting", "boundary_changes", "sequence_changes", "sample_filtering", "annotation_deduplication")),
             "H reconstruction policy differs from the supported evidence policy")
    rows = pilot.read_tsv(root / "samples.tsv")
    pilot.validate_rows(rows, [H])
    samples = {r["sample_id"]: r for r in rows}
    old_root = root / "source_provenance"
    old = _read_embedded(old_root, OLD_FILES, pilot)
    _repair_metadata(old)
    _require(pilot.sha256(old_root / "metadata.json") == metadata["source_metadata_sha256"],
             "Coordinate-repair metadata differs from reconstruction provenance")
    _require(metadata["raw_root"] == old["raw_root"], "Reconstruction raw provenance root differs")
    previous_rows = pilot.read_tsv(old_root / "samples.tsv")
    _require(len(previous_rows) == len({r["sample_id"] for r in previous_rows}), "Duplicate source sample ID")
    previous = {r["sample_id"]: r for r in previous_rows if r["dataset"] == H}
    _require(set(previous) == set(samples), "H samples differ from coordinate-repair source")
    pilot.validate_rows([{**r, "source_sha256": r["processed_sha256"]} for r in previous.values()], [H])
    old_repairs = pilot.read_tsv(old_root / "repairs.tsv")
    repair_counts = Counter(r["sample_id"] for r in old_repairs)
    inherited = [r for r in pilot.read_tsv(old_root / "review_items.tsv") if r["dataset"] == H]
    inherited_saved = pilot.read_tsv(root / "inherited_review_items.tsv")
    _require(sorted(map(_compact, inherited)) == sorted(map(_compact, inherited_saved)),
             "Inherited coordinate review evidence changed")
    inherited_counts = Counter(r["sample_id"] for r in inherited)
    for row in rows:
        before = previous[row["sample_id"]]
        for key in ("archive_relpath", "original_bgc_id", "dataset", "gcf_raw", "tsv_row_ids_json",
                    "raw_source_relpath", "raw_source_sha256", "record_count", "length_bp", "dna_sha256"):
            _require(row[key] == before[key], f"H identity/sequence provenance differs: {key}")
        _require(row["source_relpath"] == before["processed_relpath"] == row["archive_relpath"]
                 and row["source_sha256"] == before["processed_sha256"]
                 and old["artifact_sha256"].get(before["processed_relpath"]) == before["processed_sha256"],
                 "Coordinate-repaired source hash/path differs")
        _require(row["protein_input_sha256_before"] == before["protein_input_sha256_after"]
                 and row["protein_token_count_before"] == before["protein_token_count"],
                 "H input tokens disagree with coordinate-repair source")
        _require(int(before["repaired_domain_count"]) == repair_counts[row["sample_id"]],
                 "Inherited coordinate repair count disagrees with log")
        _require(row["record_count"] == "1" and row["protein_parser_status_after"] == "ok"
                 and row["dna_parser_status_after"] == "ok" and all(row[key] == "True" for key in
                 ("dna_unchanged", "cds_translations_unchanged", "non_target_features_unchanged")),
                 "H sample does not satisfy reconstruction preservation checks")
        _require(row["processed_relpath"] == row["archive_relpath"]
                 and metadata["artifact_sha256"].get(row["processed_relpath"]) == row["processed_sha256"],
                 "H processed input hash/path differs from source manifest")
        for key in ("raw_source_sha256", "source_sha256", "processed_sha256", "dna_sha256", "protein_input_sha256_after"):
            _require(_sha(row[key]), f"Invalid H {key}")
    checks = _entries(pilot.read_tsv(root / "domain_checks.tsv"), samples, "domain checks")
    _require(all(r["status"] in {"restorable", "needs_review"} for r in checks.values()),
             "Unexpected H domain check status")
    changes = _entries(pilot.read_tsv(root / "reconstructions.tsv"), samples, "reconstructions")
    review = _entries(pilot.read_tsv(root / "review_items.tsv"), samples, "review items")
    _same_log_subset(changes, {key: r for key, r in checks.items() if r["status"] == "restorable"}, "Reconstruction")
    _same_log_subset(review, {key: r for key, r in checks.items() if r["status"] == "needs_review"}, "Review")
    for change in changes.values():
        peptide = change["translation"]
        qualifiers = json.loads(change["new_qualifiers_json"])
        _require(change["feature_type"] == "misc_feature" and change["new_feature_type"] == "aSDomain"
                 and change["parent_slice_matches"] == change["dna_translation_matches"] == "True"
                 and _digest(peptide) == change["translation_sha256"] == change["dna_translation_sha256"]
                 and int(change["protein_end"]) - int(change["protein_start"]) == len(peptide)
                 and qualifiers.get("translation") == [peptide]
                 and qualifiers.get("protein_start") == [change["protein_start"]]
                 and qualifiers.get("protein_end") == [change["protein_end"]],
                 "Reconstruction lacks consistent peptide/coordinate evidence")
    mappings = _entries(pilot.read_tsv(root / "parent_id_mappings.tsv"), samples, "parent mappings")
    expected_mappings = {key: r for key, r in changes.items()
                         if json.loads(r["old_qualifiers_json"]).get("locus_tag") != json.loads(r["new_qualifiers_json"]).get("locus_tag")}
    _same_log_subset(mappings, expected_mappings, "Parent mapping")
    duplicates = _entries(pilot.read_tsv(root / "duplicate_annotations.tsv"), samples, "duplicate annotations")
    _same_log_subset(duplicates, {key: r for key, r in checks.items() if r.get("duplicate_group")}, "Duplicate annotation")
    groups = Counter((r["sample_id"], r["duplicate_group"]) for r in duplicates.values())
    _require(all(int(r["duplicate_group_size"]) == groups[r["sample_id"], r["duplicate_group"]]
                 for r in duplicates.values()), "Duplicate group size differs from retained annotations")
    per_check = Counter(r["sample_id"] for r in checks.values())
    per_change = Counter(r["sample_id"] for r in changes.values())
    per_review = Counter(r["sample_id"] for r in review.values())
    for row in rows:
        sid = row["sample_id"]
        _require(int(row["hidden_asdomain_count_before"]) == per_check[sid]
                 and int(row["restored_domain_count"]) == per_change[sid]
                 and int(row["unresolved_asdomain_count"]) == per_review[sid], "H sample restoration/review counts differ from logs")
        _require(row["asdomain_status"] == ("needs_review" if per_review[sid] else "complete")
                 and row["ready_for_encoding"] == str(not per_review[sid]), "H sample annotation readiness disagrees with review log")
        _require(row["changed"] == str(row["source_sha256"] != row["processed_sha256"]), "H changed flag disagrees with GBK hashes")
    status = "needs_review" if review else "complete"
    counts = metadata["counts"]
    expected_counts = {"sample_count": len(rows), "gcf_count": len({r["gcf_raw"] for r in rows}),
                       "restored_domains": len(changes), "hidden_domains_before": len(checks),
                       "unresolved_asdomain_annotations": len(review), "parent_id_corrections": len(mappings),
                       "asdomain_ready_samples": sum(not per_review[r["sample_id"]] for r in rows),
                       "duplicate_groups": len(groups), "duplicate_extra_entries": sum(n - 1 for n in groups.values()),
                       "protein_parser_ok": len(rows), "dna_parser_ok": len(rows), "status": status}
    _require(all(counts.get(key) == value for key, value in expected_counts.items()), "H metadata counts differ from logs")
    _require(metadata["status"] == status and metadata["ready_for_full_h_encoding"] is (not review)
             and metadata["inherited_review_items"] == len(inherited), "H snapshot readiness differs from review evidence")
    label_rel = H + "/curated_GCF_assignments.tsv"
    label_hash = pilot.sha256(root / label_rel)
    _require(label_hash == old["artifact_sha256"].get(label_rel) == old["raw_label_sha256"].get(label_rel),
             "H labels differ from original coordinate-repair labels")
    labels = pilot.read_tsv(root / label_rel)
    _require(sorted((r["BGC"], r["GCF"]) for r in labels) == sorted((r["original_bgc_id"], r["gcf_raw"]) for r in rows),
             "H labels differ from sample identities")
    return rows, metadata, repair_counts, inherited_counts


def _cohort_row(sample, pilot, *, reconstructed, repair_count, inherited_count):
    row = {key: sample[key] for key in pilot.COHORT_FIELDS if key != "source_sha256"}
    unresolved = int(sample["unresolved_asdomain_count"]) if reconstructed else 0
    return {**row, "source_sha256": sample["processed_sha256"], "raw_source_sha256": sample["raw_source_sha256"],
            "processed_source_relpath": sample["processed_relpath"],
            "source_kind": "asdomain_reconstructed" if reconstructed else "coordinate_repaired",
            "coordinate_source_sha256": sample["source_sha256"] if reconstructed else sample["processed_sha256"],
            "repair_count": str(repair_count),
            "restored_domain_count": sample["restored_domain_count"] if reconstructed else "0",
            "unresolved_asdomain_count": str(unresolved), "inherited_review_count": str(inherited_count),
            "annotation_status": "needs_review" if unresolved else "complete",
            "annotation_complete": str(not unresolved), "review_policy": "retain"}


def _review_metadata(rows):
    unresolved = sum(int(row["unresolved_asdomain_count"]) for row in rows)
    return {"annotation_status": "needs_review" if unresolved else "complete", "annotation_complete": not unresolved,
            "unresolved_asdomain_count": unresolved,
            "unresolved_sample_ids": sorted(row["sample_id"] for row in rows if int(row["unresolved_asdomain_count"])),
            "review_policy": "retain"}


def _g_rows(root, pilot):
    metadata = _read_embedded(root, (*REPAIR_FILES, G + "/curated_GCF_assignments.tsv"), pilot)
    rows = _selected_rows(root, metadata, [G], pilot)
    inherited = Counter(r["sample_id"] for r in pilot.read_tsv(root / "review_items.tsv") if r["dataset"] == G)
    return rows, metadata, inherited


def _check_g_ancestry(rows, hroot, pilot):
    """Compare G content to H's source snapshot, without equating machine paths."""
    ancestor_root = hroot / "source_provenance"
    ancestor_meta = pilot.read_json(ancestor_root / "metadata.json")
    ancestors = [r for r in pilot.read_tsv(ancestor_root / "samples.tsv") if r["dataset"] == G]
    pilot.validate_rows([{**r, "source_sha256": r["processed_sha256"]} for r in ancestors], [G])
    by_id = {r["sample_id"]: r for r in ancestors}
    counts = Counter(r["sample_id"] for r in pilot.read_tsv(ancestor_root / "repairs.tsv"))
    fields = ("sample_id", "dataset", "archive_relpath", "original_bgc_id", "gcf_raw", "tsv_row_ids_json",
              "raw_source_sha256", "processed_sha256", "dna_sha256", "protein_input_sha256_after",
              "repaired_domain_count", "record_count", "length_bp", "protein_token_count")
    _require({r["sample_id"] for r in rows} == set(by_id), "G membership differs from H's coordinate-repair ancestry")
    for row in rows:
        ancestor = by_id[row["sample_id"]]
        _require(ancestor_meta["artifact_sha256"].get(ancestor["processed_relpath"]) == ancestor["processed_sha256"]
                 and int(ancestor["repaired_domain_count"]) == counts[row["sample_id"]],
                 "H's ancestral G sample disagrees with its repair manifest")
        _require(all(row[key] == ancestor[key] for key in fields),
                 "G content differs from H's coordinate-repair ancestry")


def validate_reconstructed_cohort(prepared, metadata, rows, pilot):
    """Validate the portable cohort using embedded provenance, without CPU/GPU imports."""
    prepared = Path(prepared)
    _require(metadata.get("format") == pilot.RECONSTRUCTED_SCHEMA and metadata.get("status") == "complete",
             "Reconstructed cohort is incomplete or has the wrong format")
    datasets = metadata["datasets"]
    _require(H in datasets and not set(datasets) - {G, H}, "Reconstructed cohorts require H and optionally G")
    pilot.validate_rows(rows, datasets)
    required = {"cohort.tsv", "input_checks.tsv", "summary.md"}
    required.update("reconstruction_provenance/" + name for name in H_FILES)
    required.update(dataset + "/curated_GCF_assignments.tsv" for dataset in datasets)
    if G in datasets:
        required.update("repair_provenance/" + name for name in (*REPAIR_FILES, G + "/curated_GCF_assignments.tsv"))
    _require(not required - set(metadata["artifact_sha256"]), "Prepared manifest lacks required reconstruction provenance")
    _check_hashes(prepared, metadata["artifact_sha256"], pilot, "Prepared artifact")
    hroot = prepared / "reconstruction_provenance"
    samples, source_meta, repair_counts, inherited = _h_rows(hroot, pilot)
    _require(pilot.sha256(hroot / "metadata.json") == metadata["reconstruction_metadata_sha256"], "Embedded reconstruction metadata changed")
    _require(metadata["raw_root"] == source_meta["raw_root"], "Prepared raw provenance root differs")
    raw_roots = {H: source_meta["raw_root"]}
    expected = [_cohort_row(s, pilot, reconstructed=True, repair_count=repair_counts[s["sample_id"]],
                            inherited_count=inherited[s["sample_id"]]) for s in samples]
    source_samples = {s["sample_id"]: {**s, "protein_token_count": s["protein_token_count_after"]} for s in samples}
    if G in datasets:
        groot = prepared / "repair_provenance"
        gsamples, gmeta, ginherited = _g_rows(groot, pilot)
        _require(pilot.sha256(groot / "metadata.json") == metadata["repair_metadata_sha256"], "Embedded G repair metadata changed")
        _check_g_ancestry(gsamples, hroot, pilot)
        raw_roots[G] = gmeta["raw_root"]
        expected.extend(_cohort_row(s, pilot, reconstructed=False, repair_count=s["repaired_domain_count"],
                                    inherited_count=ginherited[s["sample_id"]]) for s in gsamples)
        source_samples.update({s["sample_id"]: s for s in gsamples})
    fields = (*pilot.COHORT_FIELDS, *EXTRA_COHORT_FIELDS)
    canonical = lambda values: sorted(tuple(str(r[key]) for key in fields) for r in values)
    _require(canonical(rows) == canonical(expected), "Prepared cohort differs from reconstruction/repair provenance")
    _require(metadata.get("raw_roots") == raw_roots, "Prepared raw provenance roots differ from their source records")
    _require(metadata["counts"] == dict(Counter(r["dataset"] for r in rows)), "Prepared sample counts differ")
    _require(all(metadata.get(key) == value for key, value in _review_metadata(rows).items()), "Prepared annotation/review policy differs from retained cohort")
    checks = pilot.read_tsv(prepared / "input_checks.tsv")
    _require(len(checks) == len(rows) and len({r["sample_id"] for r in checks}) == len(rows)
             and {r["sample_id"] for r in checks} == set(source_samples), "Input checks do not cover the full cohort")
    cohort = {r["sample_id"]: r for r in rows}
    for check in checks:
        source = source_samples[check["sample_id"]]
        row = cohort[check["sample_id"]]
        expected_check = {"dataset": source["dataset"], "source_sha256": source["processed_sha256"],
                          "protein_parser_status": "ok", "dna_parser_status": "ok",
                          "protein_input_sha256": source["protein_input_sha256_after"],
                          "dna_sha256": source["dna_sha256"], "length_bp": source["length_bp"],
                          "protein_token_count": source["protein_token_count"],
                          **{key: row[key] for key in ("annotation_status", "annotation_complete", "unresolved_asdomain_count", "review_policy")}}
        _require(all(str(check.get(key)) == str(value) for key, value in expected_check.items()),
                 f"Prepared parser check differs from source: {check['sample_id']}")
    for dataset in datasets:
        labels = pilot.read_tsv(prepared / dataset / "curated_GCF_assignments.tsv")
        subset = [r for r in rows if r["dataset"] == dataset]
        _require(sorted((r["BGC"], r["GCF"]) for r in labels) == sorted((r["sample_id"], r["gcf_raw"]) for r in subset),
                 "Prepared labels differ from stable-ID cohort")
        inputs = prepared / dataset / "inputs"
        _require({p.name for p in inputs.iterdir()} == {r["sample_id"] + ".gbk" for r in subset}, "Prepared GBK membership differs")
        for row in subset:
            _require(pilot.sha256(inputs / (row["sample_id"] + ".gbk")) == row["source_sha256"], "Prepared GBK bytes changed")


def _summary(rows, checks):
    state = _review_metadata(rows)
    return "\n".join([
        "# G/H 恢复后 encoder 输入队列", "", "队列状态：complete。仅表示身份、标签、来源和 CPU 解析检查完成。",
        f"aSDomain 注释状态：{state['annotation_status']}；保留未解决注释 {state['unresolved_asdomain_count']} 条。",
        f"全部保留 {len(rows)} 个样本，待复核样本：{', '.join(state['unresolved_sample_ids']) or '无'}。", "",
        "本队列采用 review_policy=retain，已知注释缺口随清单及运行结果保留；不自动删掉 UOC915 或其他样本。",
        "annotation_status 仅指本次隐藏 aSDomain 的恢复；旧坐标修复的 PFAM 等复核项记录为 inherited_review_count。",
        "H 来自恢复后的 GBK，G（如选用）来自原坐标修复快照；raw 和旧队列未修改。",
        "ESM2 与 BGC-MAC 共用新蛋白 token；DNA 为完整序列，不做长度过滤。未运行模型推理。",
        f"CPU 解析：{len(checks)}/{len(rows)}；超过 1021 aa 的 token 数：{sum(int(c['tokens_over_1021_residues']) for c in checks)}。", "",
        "cohort.tsv 保留原始 ID、GCF、稳定 sample_id、实际 GBK/raw/坐标修复哈希及恢复/复核计数。",
        "input_checks.tsv 来自实际 CPU 蛋白和 DNA 解析；全部必要来源清单嵌入 reconstruction_provenance/ 和可选 repair_provenance/。",
        "整个队列目录可以迁移；验证不读取旧绝对路径。新输出必须使用新缓存，不能复用旧 H ESM/MAC token 缓存。", "",
    ])


def prepare_reconstructed(args, pilot):
    """Prepare H or G/H; keep unresolved annotations and all original samples."""
    import build_inventory

    _require(args.datasets and "H" in args.datasets and not set(args.datasets) - {"G", "H"}
             and len(args.datasets) == len(set(args.datasets)), "Select H or G,H for prepare-reconstructed")
    datasets = [pilot.DATASETS[key][0] for key in args.datasets]
    source = Path(args.reconstructed_root).expanduser().resolve()
    output = Path(args.output).expanduser().resolve()
    encoder = Path(args.encoder_project).expanduser().resolve()
    hrows, hmeta, repair_counts, inherited = _h_rows(source, pilot)
    raw_roots = {H: hmeta["raw_root"]}
    sources = {H: source}
    metadata_by_source = {H: hmeta}
    samples = list(hrows)
    cohort = [_cohort_row(s, pilot, reconstructed=True, repair_count=repair_counts[s["sample_id"]],
                           inherited_count=inherited[s["sample_id"]]) for s in hrows]
    source_files = {H: H_FILES}
    if G in datasets:
        _require(args.repaired_root is not None, "G requires --repaired-root")
        gsource = Path(args.repaired_root).expanduser().resolve()
        grows, gmeta, ginherited = _g_rows(gsource, pilot)
        _check_g_ancestry(grows, source, pilot)
        raw_roots[G] = gmeta["raw_root"]
        sources[G], metadata_by_source[G] = gsource, gmeta
        source_files[G] = (*REPAIR_FILES, G + "/curated_GCF_assignments.tsv")
        samples.extend(grows)
        cohort.extend(_cohort_row(s, pilot, reconstructed=False, repair_count=s["repaired_domain_count"],
                                  inherited_count=ginherited[s["sample_id"]]) for s in grows)
    protected = [*sources.values(), (Path(pilot.ROOT) / "data/raw").resolve(),
                 *(Path(path).expanduser().resolve() for path in raw_roots.values())]
    _require(not output.exists() and all(not output.is_relative_to(path) and not path.is_relative_to(output) for path in protected),
             "Choose a new output directory separate from reconstruction, repair and raw inputs")
    initial_metadata_hashes = {}
    for dataset, root in sources.items():
        initial_metadata_hashes[dataset] = pilot.sha256(root / "metadata.json")
        _check_hashes(root, metadata_by_source[dataset]["artifact_sha256"], pilot, "Source artifact")
        selected = [r for r in samples if r["dataset"] == dataset]
        actual = {p.relative_to(root).as_posix() for p in (root / dataset / "GCF_assignments").rglob("*")
                  if p.is_file() and p.suffix.lower() in {".gbk", ".gb", ".genbank"}}
        _require(actual == {r["processed_relpath"] for r in selected}, "Source GBK membership differs from sample manifest")
    cohort.sort(key=lambda r: (r["dataset"], r["sample_id"]))
    pilot.validate_rows(cohort, datasets)
    parsers, parser_provenance = build_inventory.load_parsers(encoder)
    code_hashes = {"prepare_reconstructed_features.py": pilot.sha256(__file__),
                   "prepare_repaired_features.py": pilot.sha256(Path(__file__).with_name("prepare_repaired_features.py")),
                   "pilot_features.py": pilot.sha256(pilot.__file__), "build_inventory.py": pilot.sha256(build_inventory.__file__)}
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".pilot-reconstructed-", dir=output.parent) as temporary:
        stage = Path(temporary) / "prepared"
        stage.mkdir()
        checks = []
        by_id = {r["sample_id"]: r for r in cohort}
        for sample in sorted(samples, key=lambda r: (r["dataset"], r["sample_id"])):
            source_path = pilot.safe_child(sources[sample["dataset"]], sample["processed_relpath"])
            target = stage / sample["dataset"] / "inputs" / (sample["sample_id"] + ".gbk")
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_path, target)
            normalized = {**sample, "protein_token_count": sample["protein_token_count_after"]} if sample["dataset"] == H else sample
            checked = _input_check(target, normalized, parsers, pilot)
            row = by_id[sample["sample_id"]]
            checked.update({key: row[key] for key in ("annotation_status", "annotation_complete", "unresolved_asdomain_count", "review_policy")})
            checks.append(checked)
        pilot.write_tsv(stage / "cohort.tsv", cohort, (*pilot.COHORT_FIELDS, *EXTRA_COHORT_FIELDS))
        pilot.write_tsv(stage / "input_checks.tsv", checks, INPUT_CHECK_FIELDS)
        for dataset, root in sources.items():
            folder = "reconstruction_provenance" if dataset == H else "repair_provenance"
            for relative in source_files[dataset]:
                target = stage / folder / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(pilot.safe_child(root, relative), target)
            pilot.write_tsv(stage / dataset / "curated_GCF_assignments.tsv",
                            [{"BGC": r["sample_id"], "GCF": r["gcf_raw"]} for r in cohort if r["dataset"] == dataset], ("BGC", "GCF"))
        (stage / "summary.md").write_text(_summary(cohort, checks), encoding="utf-8")
        metadata = {"format": pilot.RECONSTRUCTED_SCHEMA, "status": "complete", "created_at_utc": pilot.now(),
                    "datasets": datasets, "counts": dict(Counter(r["dataset"] for r in cohort)),
                    **_review_metadata(cohort), "raw_root": hmeta["raw_root"], "source_root": str(source),
                    "raw_roots": raw_roots,
                    "source_roots": {dataset: str(root) for dataset, root in sources.items()},
                    "reconstruction_metadata_sha256": initial_metadata_hashes[H],
                    "encoder_parsers": parser_provenance, "script_sha256": code_hashes,
                    "policy": {"identity": "inventory stable sample_id", "source_sha256": "actual processed GBK bytes",
                               "raw_source_sha256": "original backup provenance", "review_policy": "retain",
                               "annotation_status_scope": "existing Geneious aSDomain reconstruction",
                               "sample_exclusions": [], "deduplicated": False, "length_filter": None,
                               "protein_window_residues": 1021, "dna_sequence_mode": "full", "model_inference": False},
                    "artifact_sha256": {p.relative_to(stage).as_posix(): pilot.sha256(p) for p in sorted(stage.rglob("*")) if p.is_file()}}
        if G in sources:
            metadata["repair_metadata_sha256"] = initial_metadata_hashes[G]
        pilot.write_json(stage / "metadata.json", metadata)
        validate_reconstructed_cohort(stage, metadata, cohort, pilot)
        for dataset, root in sources.items():
            _check_hashes(root, metadata_by_source[dataset]["artifact_sha256"], pilot, "Source after preparation")
            _require(pilot.sha256(root / "metadata.json") == initial_metadata_hashes[dataset], "Source metadata changed during preparation")
        _check_hashes(encoder, parser_provenance["source_sha256"], pilot, "Encoder parser source")
        _check_hashes(Path(__file__).parent, code_hashes, pilot, "Preparation code")
        _require(not output.exists(), "Output appeared during preparation; refusing overwrite")
        stage.rename(output)
    print(f"Prepared {len(cohort)} reconstructed G/H BGCs: {output}; annotation_status={metadata['annotation_status']}; review_policy=retain")
    return output
