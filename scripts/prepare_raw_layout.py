#!/usr/bin/env python3
"""Move the original curated data once, preserving and verifying every file.

Only this script relocates raw data. It never edits GenBank or label contents.
An interrupted run leaves files.tsv and plan.json as recovery evidence and is
not resumed automatically. A completed run can be rerun to verify its evidence.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import stat
import sys


SCHEMA_VERSION = "raw_relocation_v1"
SOURCE = Path("data/curated_datasets")
DESTINATION = Path("data/raw/curated_datasets")
REPORT = Path("data/processed/raw_relocation")
FILE_FIELDS = ("original_path", "raw_path", "relative_path", "size_bytes", "sha256")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_directory(path: Path) -> None:
    if path.is_symlink() or not path.is_dir():
        raise ValueError(f"Expected a real directory, not a symlink: {path}")


def snapshot_tree(root: Path) -> tuple[list[dict], list[str]]:
    """Return all regular files and directories; reject links/special entries."""
    _require_directory(root)
    files = []
    directories = ["."]
    for current, dirnames, filenames in os.walk(root, followlinks=False):
        current = Path(current)
        dirnames.sort()
        filenames.sort()
        for name in dirnames:
            child = current / name
            if not stat.S_ISDIR(child.lstat().st_mode):
                raise ValueError(f"Symlink or special directory entry is unsupported: {child}")
            directories.append(child.relative_to(root).as_posix())
        for name in filenames:
            child = current / name
            before = child.lstat()
            if not stat.S_ISREG(before.st_mode):
                raise ValueError(f"Symlink or special file entry is unsupported: {child}")
            checksum = sha256_file(child)
            after = child.lstat()
            if (before.st_size, before.st_mtime_ns, before.st_ino) != (
                after.st_size, after.st_mtime_ns, after.st_ino
            ):
                raise ValueError(f"File changed while hashing: {child}")
            files.append({"relative_path": child.relative_to(root).as_posix(),
                          "size_bytes": before.st_size, "sha256": checksum})
    return sorted(files, key=lambda row: row["relative_path"]), sorted(directories)


def tree_sha256(files: list[dict], directories: list[str]) -> str:
    content = json.dumps({"files": files, "directories": directories},
                         ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def _write_json_new(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _validate_parent_paths(workspace: Path) -> None:
    _require_directory(workspace)
    _require_directory(workspace / "data")
    for relative in (Path("data/raw"), Path("data/processed")):
        path = workspace / relative
        if path.exists() or path.is_symlink():
            _require_directory(path)


def _verify_completed(workspace: Path) -> dict:
    report = workspace / REPORT
    _require_directory(report)
    metadata_path = report / "metadata.json"
    if not metadata_path.is_file() or metadata_path.is_symlink():
        raise ValueError(f"Incomplete relocation evidence at {report}; inspect plan.json "
                         "and both source/destination before any manual recovery.")
    for name in ("metadata.json", "plan.json", "files.tsv"):
        path = report / name
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Missing or non-regular relocation evidence: {path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    plan = json.loads((report / "plan.json").read_text(encoding="utf-8"))
    if metadata.get("status") != "verified" or plan.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Relocation evidence is not a completed, supported report.")
    if plan.get("source") != SOURCE.as_posix() or plan.get("destination") != DESTINATION.as_posix():
        raise ValueError("Relocation evidence names unexpected source/destination paths.")
    if (workspace / SOURCE).exists() or (workspace / SOURCE).is_symlink():
        raise ValueError("Original source still exists alongside a completed relocation report.")
    if sha256_file(report / "files.tsv") != metadata.get("files_tsv_sha256"):
        raise ValueError("Relocation file inventory checksum changed.")
    if sha256_file(report / "plan.json") != metadata.get("plan_json_sha256"):
        raise ValueError("Relocation plan checksum changed.")
    with (report / "files.tsv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(FILE_FIELDS):
            raise ValueError("Unexpected relocation file inventory columns.")
        expected = []
        for row in reader:
            rel = row["relative_path"]
            if row["original_path"] != (SOURCE / rel).as_posix() or row["raw_path"] != (DESTINATION / rel).as_posix():
                raise ValueError("Relocation path mapping is inconsistent.")
            expected.append({"relative_path": rel, "size_bytes": int(row["size_bytes"]),
                             "sha256": row["sha256"]})
    actual, directories = snapshot_tree(workspace / DESTINATION)
    if actual != expected or directories != plan["directories"]:
        raise ValueError("Raw data differ from the recorded relocation inventory.")
    if tree_sha256(actual, directories) != plan["tree_sha256"]:
        raise ValueError("Raw tree checksum differs from relocation plan.")
    return {**metadata, "already_verified": True}


def prepare_layout(workspace: Path) -> dict:
    workspace = Path(workspace).absolute()
    _validate_parent_paths(workspace)
    source, destination, report = (workspace / path for path in (SOURCE, DESTINATION, REPORT))
    if report.exists() or report.is_symlink():
        return _verify_completed(workspace)
    if destination.exists() or destination.is_symlink():
        raise ValueError(f"Refusing to merge into or overwrite existing destination: {destination}")
    before, directories = snapshot_tree(source)
    if not before:
        raise ValueError(f"Refusing to relocate an empty source tree: {source}")

    (workspace / "data/raw").mkdir(exist_ok=True)
    (workspace / "data/processed").mkdir(exist_ok=True)
    if source.stat().st_dev != destination.parent.stat().st_dev:
        raise ValueError("Source and destination are on different filesystems; atomic rename required.")
    report.mkdir()  # Never overwrite an existing report, including an interrupted run.
    files_path = report / "files.tsv"
    with files_path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FILE_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in before:
            writer.writerow({"original_path": (SOURCE / row["relative_path"]).as_posix(),
                             "raw_path": (DESTINATION / row["relative_path"]).as_posix(), **row})
        handle.flush()
        os.fsync(handle.fileno())
    plan = {"schema_version": SCHEMA_VERSION, "created_at_utc": _utc_now(),
            "source": SOURCE.as_posix(), "destination": DESTINATION.as_posix(),
            "file_count": len(before), "total_bytes": sum(row["size_bytes"] for row in before),
            "directories": directories, "tree_sha256": tree_sha256(before, directories),
            "script_sha256": sha256_file(Path(__file__)), "python_version": sys.version,
            "operation": "same-filesystem directory rename; no content changes"}
    _write_json_new(report / "plan.json", plan)
    if destination.exists() or destination.is_symlink():
        raise ValueError(f"Destination appeared during preparation; refusing move: {destination}")
    source.rename(destination)
    after, after_directories = snapshot_tree(destination)
    if before != after or directories != after_directories:
        raise ValueError(f"Post-move verification failed; preserve raw data and evidence in {report}")
    metadata = {"schema_version": SCHEMA_VERSION, "status": "verified",
                "verified_at_utc": _utc_now(), "file_count": len(after),
                "total_bytes": plan["total_bytes"], "tree_sha256": plan["tree_sha256"],
                "files_tsv_sha256": sha256_file(files_path),
                "plan_json_sha256": sha256_file(report / "plan.json")}
    _write_json_new(report / "metadata.json", metadata)
    return {**metadata, "already_verified": False}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args(argv)
    try:
        result = prepare_layout(args.workspace)
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        parser.exit(1, f"Raw relocation failed: {exc}\n")
    action = "Already relocated; verified again" if result["already_verified"] else "Relocated and verified"
    print(f"{action}: {result['file_count']} files, {result['total_bytes']} bytes.")
    print(f"Evidence: {args.workspace / REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
