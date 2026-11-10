"""Integrity and refusal checks for relocation, using temporary data only."""

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_raw_layout.py"
SPEC = importlib.util.spec_from_file_location("prepare_raw_layout", SCRIPT)
layout = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(layout)


class RawLayoutTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.workspace = Path(self.temp.name)
        self.source = self.workspace / layout.SOURCE
        family = self.source / "Dataset_A" / "family\u00a0"
        family.mkdir(parents=True)
        (family / "sample.gbk").write_bytes(b"LOCUS example\nORIGIN\n acgt\n//\n")
        (self.source / "README.md").write_bytes(b"source information\n")
        (self.source / "labels.tsv").write_bytes(b"BGC\tGCF\nsample\tfamily\xc2\xa0\n")
        (self.source / "empty_folder").mkdir()

    def test_preserves_all_bytes_paths_and_empty_directories_and_reruns(self):
        before = layout.snapshot_tree(self.source)
        result = layout.prepare_layout(self.workspace)
        self.assertFalse(result["already_verified"])
        self.assertFalse(self.source.exists())
        self.assertEqual(layout.snapshot_tree(self.workspace / layout.DESTINATION), before)
        with (self.workspace / layout.REPORT / "files.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 3)
        self.assertTrue(any(row["relative_path"].endswith("labels.tsv") for row in rows))
        self.assertTrue(all(row["original_path"].startswith("data/curated_datasets/") for row in rows))
        self.assertTrue(layout.prepare_layout(self.workspace)["already_verified"])

    def test_refuses_existing_destination_without_touching_source(self):
        destination = self.workspace / layout.DESTINATION
        destination.mkdir(parents=True)
        (destination / "existing").write_text("do not overwrite")
        with self.assertRaisesRegex(ValueError, "existing destination"):
            layout.prepare_layout(self.workspace)
        self.assertTrue(self.source.exists())
        self.assertEqual((destination / "existing").read_text(), "do not overwrite")
        self.assertFalse((self.workspace / layout.REPORT).exists())

    def test_refuses_symlink_input_and_parent(self):
        link = self.source / "link"
        link.symlink_to(self.source / "README.md")
        with self.assertRaisesRegex(ValueError, "Symlink"):
            layout.prepare_layout(self.workspace)
        link.unlink()
        (self.workspace / "data/raw").symlink_to(self.source, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "symlink"):
            layout.prepare_layout(self.workspace)
        self.assertTrue(self.source.exists())

    def test_tampered_raw_file_rejected_on_verification(self):
        layout.prepare_layout(self.workspace)
        (self.workspace / layout.DESTINATION / "README.md").write_text("changed")
        with self.assertRaisesRegex(ValueError, "Raw data differ"):
            layout.prepare_layout(self.workspace)

    def test_partial_move_failure_keeps_plan_and_source(self):
        with patch.object(Path, "rename", side_effect=OSError("simulated rename failure")):
            with self.assertRaisesRegex(OSError, "simulated"):
                layout.prepare_layout(self.workspace)
        report = self.workspace / layout.REPORT
        self.assertTrue(self.source.exists())
        self.assertTrue((report / "files.tsv").is_file())
        self.assertEqual(json.loads((report / "plan.json").read_text())["file_count"], 3)
        self.assertFalse((report / "metadata.json").exists())
        with self.assertRaisesRegex(ValueError, "Incomplete relocation"):
            layout.prepare_layout(self.workspace)

    def test_changed_tree_after_move_does_not_publish_verified_metadata(self):
        real_snapshot = layout.snapshot_tree
        def changed_snapshot(root):
            rows, directories = real_snapshot(root)
            if root == self.workspace / layout.DESTINATION:
                rows[0] = {**rows[0], "sha256": "wrong"}
            return rows, directories
        with patch.object(layout, "snapshot_tree", side_effect=changed_snapshot):
            with self.assertRaisesRegex(ValueError, "Post-move verification"):
                layout.prepare_layout(self.workspace)
        self.assertTrue((self.workspace / layout.DESTINATION).is_dir())
        self.assertFalse((self.workspace / layout.REPORT / "metadata.json").exists())
        self.assertTrue((self.workspace / layout.REPORT / "plan.json").exists())


if __name__ == "__main__":
    unittest.main()
