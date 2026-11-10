#!/usr/bin/env python3
"""Prepare and extract the fixed B/E pilot using the existing BGC_annotation CLIs."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
DATASETS = {"B": ("Dataset_B_Divergent", 39), "E": ("Dataset_E_MIBiG_BS1", 281)}
ROUTES = ("esm2", "mac", "dna_full")
INVENTORY_SCHEMA = "curated_bgc_inventory_v1"
PREPARED_SCHEMA = "bgc_be_pilot_inputs_v1"
COHORT_FIELDS = ("dataset", "sample_id", "original_bgc_id", "archive_relpath",
                 "source_sha256", "gcf_raw", "tsv_row_ids_json")
ENCODER_FILES = (
    "experiment/extract_features.py", "experiment/extract_dna_features.py",
    "data_preparation/feature_cache.py", "data_preparation/feature_inputs.py",
    "data_preparation/frozen_esm.py", "data_preparation/dna_inputs.py",
    "data_preparation/frozen_dna.py", "model/features.py", "model/encoders.py",
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    """Replace only our own small status file, with atomic publication."""
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")
    temp.replace(path)


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_child(root, relative):
    rel = PurePosixPath(relative)
    if rel.is_absolute() or ".." in rel.parts or not rel.parts:
        raise ValueError(f"Unsafe relative path: {relative}")
    path = root.joinpath(*rel.parts)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Source escapes its root: {relative}")
    return path


def validate_rows(rows, dataset_names):
    if len(dataset_names) != len(set(dataset_names)):
        raise ValueError("Repeated dataset selection")
    expected_counts = dict(DATASETS.values())
    if not dataset_names or set(dataset_names) - set(expected_counts):
        raise ValueError("Only the audited B/E pilot datasets are supported")
    if {r["dataset"] for r in rows} != set(dataset_names):
        raise ValueError("Unexpected dataset membership in cohort")
    ids = [r["sample_id"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate sample IDs in cohort")
    for dataset in dataset_names:
        subset = [r for r in rows if r["dataset"] == dataset]
        if len(subset) != expected_counts[dataset]:
            raise ValueError(f"{dataset}: expected {expected_counts[dataset]} samples, got {len(subset)}")
        if len({r["original_bgc_id"] for r in subset}) != len(subset):
            raise ValueError(f"{dataset}: duplicate original BGC names")
    for row in rows:
        rel = PurePosixPath(row["archive_relpath"])
        expected_id = "bgc_" + hashlib.sha256(
            (INVENTORY_SCHEMA + "\0" + row["archive_relpath"]).encode("utf-8")
        ).hexdigest()[:24]
        if row["sample_id"] != expected_id or rel.parts[0] != row["dataset"]:
            raise ValueError(f"Invalid archive identity: {row['sample_id']}")
        if not row["gcf_raw"] or rel.stem != row["original_bgc_id"]:
            raise ValueError(f"Invalid original ID/label: {row['sample_id']}")


def prepare(args):
    raw, inventory, output = (p.expanduser().resolve() for p in
                              (args.raw_root, args.inventory, args.output))
    if output.exists() or output.is_relative_to(raw):
        raise ValueError("Choose a new output directory outside raw data")
    metadata = read_json(inventory / "metadata.json")
    if metadata.get("schema_version") != INVENTORY_SCHEMA:
        raise ValueError("Unsupported inventory schema")
    inventory_hashes = {}
    for name in ("samples.tsv", "source_files.tsv"):
        inventory_hashes[name] = sha256(inventory / name)
        if inventory_hashes[name] != metadata["output_sha256"].get(name):
            raise ValueError(f"Inventory checksum mismatch: {name}")
    dataset_names = [DATASETS[key][0] for key in args.datasets]
    samples = [r for r in read_tsv(inventory / "samples.tsv") if r["dataset"] in dataset_names]
    sources = {r["archive_relpath"]: r for r in read_tsv(inventory / "source_files.tsv")}
    rows = []
    for sample in samples:
        if (sample["mapping_status"] != "exact_path_label" or sample["record_count"] != "1"
                or sample["dna_parser_status"] != "ok" or sample["protein_parser_status"] != "ok"
                or sample["excluded"] != "False" or sample["repairs_applied_json"] != "[]"):
            raise ValueError(f"Pilot entry no longer meets audited policy: {sample['sample_id']}")
        labels = json.loads(sample["tsv_labels_raw_json"])
        if len(labels) != 1 or labels[0] != sample["folder_gcf_raw"]:
            raise ValueError(f"Unresolved label: {sample['sample_id']}")
        source = safe_child(raw, sample["archive_relpath"])
        if sha256(source) != sample["source_sha256"]:
            raise ValueError(f"Server GBK differs from inventory: {source}")
        rows.append({**{key: sample[key] for key in COHORT_FIELDS if key != "gcf_raw"},
                     "gcf_raw": labels[0]})
    rows.sort(key=lambda r: (r["dataset"], r["sample_id"]))
    validate_rows(rows, dataset_names)
    label_hashes = {}
    for dataset in dataset_names:
        subset = [r for r in rows if r["dataset"] == dataset]
        actual_paths = {p.relative_to(raw).as_posix()
                        for p in (raw / dataset / "GCF_assignments").rglob("*")
                        if p.is_file() and p.suffix.lower() in {".gbk", ".gb", ".genbank"}}
        if actual_paths != {r["archive_relpath"] for r in subset}:
            raise ValueError(f"Raw GBK file set differs from inventory: {dataset}")
        rel = dataset + "/curated_GCF_assignments.tsv"
        label_hashes[rel] = sha256(raw / rel)
        if label_hashes[rel] != sources[rel]["source_sha256"]:
            raise ValueError(f"Server labels differ from inventory: {rel}")
        labels = read_tsv(raw / rel)
        if len(labels) != len(subset) or sorted((r["BGC"], r["GCF"]) for r in labels) != sorted(
            (r["original_bgc_id"], r["gcf_raw"]) for r in subset
        ):
            raise ValueError(f"Raw labels disagree with pilot cohort: {dataset}")
    output.parent.mkdir(parents=True, exist_ok=True)
    # Build a private staging directory, then publish the complete cohort at once.
    with tempfile.TemporaryDirectory(prefix=".pilot-", dir=output.parent) as tmp:
        stage = Path(tmp) / "prepared"
        stage.mkdir()
        for row in rows:
            dest = stage / row["dataset"] / "inputs" / (row["sample_id"] + ".gbk")
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(safe_child(raw, row["archive_relpath"]), dest)
            if sha256(dest) != row["source_sha256"]:
                raise ValueError(f"Source changed during copy: {row['archive_relpath']}")
        write_tsv(stage / "cohort.tsv", rows, COHORT_FIELDS)
        for dataset in dataset_names:
            write_tsv(stage / dataset / "curated_GCF_assignments.tsv",
                      [{"BGC": r["sample_id"], "GCF": r["gcf_raw"]}
                       for r in rows if r["dataset"] == dataset], ("BGC", "GCF"))
        artifact_hashes = {p.relative_to(stage).as_posix(): sha256(p)
                           for p in stage.rglob("*.tsv")}
        write_json(stage / "metadata.json", {
            "format": PREPARED_SCHEMA, "status": "complete", "created_at_utc": now(),
            "datasets": dataset_names, "counts": {d: sum(r["dataset"] == d for r in rows)
                                                  for d in dataset_names},
            "raw_root": str(raw), "inventory": str(inventory),
            "inventory_sha256": inventory_hashes, "raw_label_sha256": label_hashes,
            "artifact_sha256": artifact_hashes, "script_sha256": sha256(__file__),
            "policy": {"identity": "inventory sample_id; byte-identical GBK copies",
                       "labels": "raw Unicode retained; separate from encoder inputs",
                       "length_filter": None, "deduplicated": False, "repairs": []},
        })
        stage.rename(output)
    print(f"Prepared {len(rows)} BGCs: {output}")


def load_prepared(path):
    path = Path(path).expanduser().resolve()
    metadata = read_json(path / "metadata.json")
    if metadata.get("format") != PREPARED_SCHEMA or metadata.get("status") != "complete":
        raise ValueError("Prepared cohort is incomplete or has the wrong format")
    for name, digest in metadata["artifact_sha256"].items():
        if sha256(safe_child(path, name)) != digest:
            raise ValueError(f"Prepared manifest changed: {name}")
    rows = read_tsv(path / "cohort.tsv")
    validate_rows(rows, metadata["datasets"])
    for dataset in metadata["datasets"]:
        subset = [r for r in rows if r["dataset"] == dataset]
        inputs = path / dataset / "inputs"
        if {p.name for p in inputs.iterdir()} != {r["sample_id"] + ".gbk" for r in subset}:
            raise ValueError(f"Prepared input file set changed: {dataset}")
        for row in subset:
            if sha256(inputs / (row["sample_id"] + ".gbk")) != row["source_sha256"]:
                raise ValueError(f"Prepared input changed: {row['sample_id']}")
    return path, metadata, rows


def file_signature(path):
    path = Path(path).expanduser().resolve()
    if path.is_file():
        return {"path": str(path), "sha256": sha256(path)}
    if path.is_dir():
        files = {p.relative_to(path).as_posix(): sha256(p) for p in sorted(path.rglob("*"))
                 if p.is_file() and p.suffix.lower() in {".json", ".txt", ".model", ".py", ".bin", ".safetensors"}}
        if files:
            return {"path": str(path), "files": files}
    raise ValueError(f"Missing or empty model/config/tokenizer path: {path}")


def extraction_config(args, prepared, prep_metadata):
    routes = [r for r in ROUTES if r in args.routes or (r == "esm2" and "mac" in args.routes)]
    project = args.encoder_project.expanduser().resolve()
    checkpoints = {}
    for field, needed in (("esm_checkpoint", "esm2" in routes),
                          ("mac_checkpoint", "mac" in routes), ("dna_checkpoint", "dna_full" in routes),
                          ("dna_tokenizer", "dna_full" in routes and args.dna_tokenizer is not None),
                          ("dna_config", "dna_full" in routes and args.dna_config is not None)):
        value = getattr(args, field)
        if needed:
            if value is None:
                raise ValueError(f"Required option: --{field.replace('_', '-')}")
            checkpoints[field] = file_signature(value)
    if args.esm_batch_size <= 0 or args.batch_size <= 0 or args.mac_member < 0:
        raise ValueError("Batch sizes must be positive and MAC member nonnegative")
    if args.dna_attention_backend == "flash_attention_2" and args.dna_dtype == "float32" and "dna_full" in routes:
        raise ValueError("FlashAttention-2 requires bfloat16 or float16")
    return {
        "prepared": str(prepared), "cohort_sha256": prep_metadata["artifact_sha256"]["cohort.tsv"],
        "routes": routes, "encoder_project": str(project), "python": args.python,
        "checkpoints": checkpoints, "device": args.device, "batch_size": args.batch_size,
        "esm_batch_size": args.esm_batch_size, "esm_pooling": "bos", "max_residues": 1021,
        "mac_member": args.mac_member, "trusted_input": args.trusted_input,
        "dna_dtype": args.dna_dtype, "dna_attention_backend": args.dna_attention_backend,
        "sequence_mode": "full", "max_bp": None, "max_tokens": None,
        "source_sha256": {name: sha256(project / name) for name in ENCODER_FILES},
        "runner_sha256": sha256(__file__),
        "validator_sha256": sha256(Path(__file__).with_name("pilot_feature_validation.py")),
    }


def command_for(config, dataset, route, output):
    base = [config["python"], "-u", "-m"]
    inputs = Path(config["prepared"]) / dataset / "inputs"
    route_dir = output / dataset / route
    cache = output / dataset / "esm2" / "tokens.pth"
    checkpoints = config["checkpoints"]
    if route == "dna_full":
        command = base + ["experiment.extract_dna_features", "--gbk", str(inputs),
                          "--checkpoint", checkpoints["dna_checkpoint"]["path"],
                          "--sequence-mode", "full", "--window-batch-size", "1",
                          "--dtype", config["dna_dtype"],
                          "--attention-backend", config["dna_attention_backend"]]
        for field, flag in (("dna_tokenizer", "--tokenizer"), ("dna_config", "--config")):
            if field in checkpoints:
                command += [flag, checkpoints[field]["path"]]
    else:
        command = base + ["experiment.extract_features", "--route", "esm" if route == "esm2" else "mac",
                          "--batch-size", str(config["batch_size"])]
        if route == "esm2":
            command += ["--gbk", str(inputs), "--esm-checkpoint", checkpoints["esm_checkpoint"]["path"],
                        "--esm-pooling", "bos", "--max-residues", "1021",
                        "--esm-batch-size", str(config["esm_batch_size"]), "--save-cache", str(cache)]
        else:
            command += ["--esm-cache", str(cache), "--checkpoint", checkpoints["mac_checkpoint"]["path"],
                        "--member", str(config["mac_member"])]
        if config["trusted_input"]:
            command.append("--trusted-input")
    return command + ["--device", config["device"], "--output", str(route_dir)]


def check_step(output, dataset, route, rows, saved=None):
    from pilot_feature_validation import validate_export
    cache = output / dataset / "esm2" / "tokens.pth" if route in {"esm2", "mac"} else None
    result = validate_export(output / dataset / route, rows, route, cache_path=cache)
    if saved is not None and result != saved:
        raise ValueError(f"Previously validated artifacts changed: {dataset}/{route}")
    return result


def extract(args):
    prepared, prep_metadata, rows = load_prepared(args.prepared)
    output = args.output.expanduser().resolve()
    if output.is_relative_to(prepared) or output.is_relative_to(Path(prep_metadata["raw_root"])):
        raise ValueError("Feature output must be separate from prepared inputs and raw data")
    config = extraction_config(args, prepared, prep_metadata)
    state_path = output / "run.json"
    if output.exists():
        if not args.resume or not state_path.is_file():
            raise ValueError("Output already exists; use --resume for the same run or choose a new output")
        state = read_json(state_path)
        if state.get("config") != config:
            raise ValueError("Run configuration/code/checkpoint changed; choose a new output directory")
    else:
        output.mkdir(parents=True)
        state = {"format": "bgc_be_feature_run_v1", "status": "running", "created_at_utc": now(),
                 "config": config, "datasets": prep_metadata["datasets"], "steps": {}}
        write_json(state_path, state)
    # A lock prevents simultaneous resume commands from sharing GPU output paths.
    lock = output / ".running.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise ValueError(f"Run lock exists: {lock}; check for a running process before removing a stale lock") from exc
    with os.fdopen(descriptor, "w") as handle:
        handle.write(f"pid={os.getpid()} time={now()}\n")
    try:
        # Reload only after holding the lock; another process may just have finished.
        state = read_json(state_path)
        if state.get("config") != config:
            raise ValueError("Run configuration changed while acquiring the lock")
        state["status"] = "running"
        state.pop("error", None)
        write_json(state_path, state)
        for dataset in prep_metadata["datasets"]:
            subset = sorted((r for r in rows if r["dataset"] == dataset), key=lambda r: r["sample_id"])
            for route in config["routes"]:
                key = dataset + "/" + route
                previous = state["steps"].get(key, {})
                if previous.get("status") == "complete":
                    check_step(output, dataset, route, subset, previous["validation"])
                    print(f"Verified; reusing {key}", flush=True)
                    continue
                route_dir = output / dataset / route
                if route_dir.exists():
                    archive = output / "_failed" / dataset
                    archive.mkdir(parents=True, exist_ok=True)
                    archive = archive / (route + "_" + str(time.time_ns()))
                    route_dir.rename(archive)
                    state.setdefault("archived_attempts", []).append(str(archive.relative_to(output)))
                route_dir.mkdir(parents=True)
                command = command_for(config, dataset, route, output)
                step = {"status": "running", "started_at_utc": now(), "command": command}
                state["steps"][key] = step
                write_json(state_path, state)
                print(f"Running {key} ({len(subset)} samples)\n{shlex.join(command)}", flush=True)
                start = time.monotonic()
                try:
                    with (route_dir / "execution.log").open("w", encoding="utf-8") as log:
                        result = subprocess.run(command, cwd=config["encoder_project"], stdout=log,
                                                stderr=subprocess.STDOUT, check=False)
                    step["returncode"] = result.returncode
                    if result.returncode:
                        raise RuntimeError(f"{key} failed (exit {result.returncode}); see {route_dir / 'execution.log'}")
                    step["validation"] = check_step(output, dataset, route, subset)
                    step["status"] = "complete"
                    print(f"Validated {key}: {step['validation']['count']} vectors", flush=True)
                except BaseException as exc:
                    step["status"] = "failed"
                    step["error"] = f"{type(exc).__name__}: {exc}"
                    raise
                finally:
                    step["elapsed_seconds"] = time.monotonic() - start
                    step["finished_at_utc"] = now()
                    write_json(state_path, state)
        state["status"] = "complete"
        state["finished_at_utc"] = now()
    except BaseException as exc:
        state["status"] = "failed"
        state["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        write_json(state_path, state)
        lock.unlink()
    print(f"All requested routes validated: {state_path}")


def validate(args):
    prepared, prep_metadata, rows = load_prepared(args.prepared)
    output = args.output.expanduser().resolve()
    state = read_json(output / "run.json")
    if (state.get("format") != "bgc_be_feature_run_v1" or state.get("status") != "complete"
            or state["datasets"] != prep_metadata["datasets"]
            or state["config"]["cohort_sha256"] != sha256(prepared / "cohort.tsv")):
        raise ValueError("Run is incomplete or belongs to a different cohort")
    for dataset in prep_metadata["datasets"]:
        subset = sorted((r for r in rows if r["dataset"] == dataset), key=lambda r: r["sample_id"])
        for route in state["config"]["routes"]:
            step = state["steps"][dataset + "/" + route]
            if step["status"] != "complete":
                raise ValueError(f"Incomplete step: {dataset}/{route}")
            result = check_step(output, dataset, route, subset, step["validation"])
            print(f"{dataset}/{route}: {result['count']} x {result['dimension']}, validated")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    prep = commands.add_parser("prepare", help="Copy the complete audited B/E cohort; no model dependencies")
    prep.add_argument("--raw-root", type=Path, default=ROOT / "data/raw/curated_datasets")
    prep.add_argument("--inventory", type=Path, default=ROOT / "data/processed/inventory_v1")
    prep.add_argument("--output", type=Path, default=ROOT / "data/processed/pilot_be_v1")
    prep.add_argument("--datasets", nargs="+", choices=tuple(DATASETS), default=["B", "E"])
    prep.set_defaults(handler=prepare)
    run = commands.add_parser("extract", help="Run frozen ESM2, MAC and full-DNA extraction")
    run.add_argument("--prepared", type=Path, required=True)
    run.add_argument("--encoder-project", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--routes", nargs="+", choices=ROUTES, default=list(ROUTES))
    for name in ("esm-checkpoint", "mac-checkpoint", "dna-checkpoint", "dna-tokenizer", "dna-config"):
        run.add_argument("--" + name, type=Path)
    run.add_argument("--python", default=sys.executable, help="Interpreter for the BGC_annotation environment")
    run.add_argument("--device", default="cuda:0")
    run.add_argument("--esm-batch-size", type=int, default=8)
    run.add_argument("--batch-size", type=int, default=16)
    run.add_argument("--mac-member", type=int, default=0)
    run.add_argument("--dna-dtype", choices=("bfloat16", "float16", "float32"), default="bfloat16")
    run.add_argument("--dna-attention-backend", choices=("flash_attention_2", "sdpa", "eager"),
                     default="flash_attention_2")
    run.add_argument("--trusted-input", action="store_true", help="Permit legacy loading only for trusted protein checkpoints")
    run.add_argument("--resume", action="store_true", help="Verify completed stages; archive and restart failed stages")
    run.set_defaults(handler=extract)
    check = commands.add_parser("validate", help="Validate every requested output against the complete cohort")
    check.add_argument("--prepared", type=Path, required=True)
    check.add_argument("--output", type=Path, required=True)
    check.set_defaults(handler=validate)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.handler(args)
    except (ValueError, OSError, KeyError, RuntimeError, TypeError) as exc:
        parser.exit(1, f"Pilot workflow failed: {exc}\n")


if __name__ == "__main__":
    main()
