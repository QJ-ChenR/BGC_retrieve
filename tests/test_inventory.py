"""Small, synthetic provenance and coverage checks for the inventory builder."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_inventory import build_inventory


EXPECTED_OUTPUTS = {
    "samples.tsv", "records.tsv", "label_rows.tsv", "duplicate_groups.tsv",
    "issues.tsv", "source_files.tsv", "dataset_summary.tsv", "metadata.json",
    "summary.md",
}


class InventoryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="bgc-inventory-test-")
        self.addCleanup(temporary.cleanup)
        self.workspace = Path(temporary.name)
        self.source = self.workspace / "data" / "raw" / "curated_datasets"
        self.source.mkdir(parents=True)
        self.output = self.workspace / "data" / "processed" / "inventory"

    def dataset(self, labels, name="Dataset_A_Test"):
        directory = self.source / name
        directory.mkdir(parents=True)
        with (directory / "curated_GCF_assignments.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(("BGC", "GCF"))
            writer.writerows(labels)
        return directory

    @staticmethod
    def record(sequence="ATGAAAGCTTAA", record_id="record_1", translation="MKA"):
        record = SeqRecord(Seq(sequence), id=record_id, name=record_id,
                           description="synthetic inventory fixture")
        record.annotations = {"molecule_type": "DNA", "organism": "Test organism"}
        qualifiers = {"locus_tag": [record_id + "_cds"]}
        if translation is not None:
            qualifiers["translation"] = [translation]
        record.features = [SeqFeature(
            SimpleLocation(0, min(9, len(sequence)), strand=1),
            type="CDS", qualifiers=qualifiers,
        )]
        return record

    def gbk(self, dataset, family, basename, records=None):
        path = dataset / "GCF_assignments" / family / (basename + ".gbk")
        path.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records if records is not None else [self.record()], path, "genbank")
        return path

    @staticmethod
    def rows(output, filename):
        with (output / filename).open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    @staticmethod
    def source_hashes(source):
        return {
            path.relative_to(source).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in source.rglob("*") if path.is_file()
        }

    def run_inventory(self, *, source=None, output=None, encoder_project=None):
        source = source if source is not None else self.source
        output = output if output is not None else self.output
        before = self.source_hashes(source)
        summary = build_inventory(
            source_root=source, output_dir=output, workspace=self.workspace,
            encoder_project=encoder_project,
        )
        self.assertIsInstance(summary, dict)
        self.assertEqual(before, self.source_hashes(source), "Inventory changed source bytes")
        self.assertTrue(EXPECTED_OUTPUTS.issubset({path.name for path in output.iterdir()}))
        with (output / "metadata.json").open(encoding="utf-8") as handle:
            self.assertIsInstance(json.load(handle), dict)
        return output

    @staticmethod
    def is_true(value):
        return value.lower() in {"true", "1", "yes"}

    def test_same_basename_different_dna_keeps_two_source_identities(self):
        dataset = self.dataset([("shared", "family_a"), ("shared", "family_b")])
        self.gbk(dataset, "family_a", "shared")
        self.gbk(dataset, "family_b", "shared", [self.record("ATGCCCGGGTAA")])
        output = self.run_inventory()

        samples = self.rows(output, "samples.tsv")
        self.assertEqual(len(samples), 2)
        self.assertEqual(len({row["sample_id"] for row in samples}), 2)
        self.assertEqual(len({row["source_sha256"] for row in samples}), 2)
        self.assertEqual({row["mapping_status"] for row in samples}, {"exact_path_label"})
        self.assertEqual({row["dna_parser_status"] for row in samples}, {"not_checked"})
        self.assertEqual({row["protein_parser_status"] for row in samples}, {"not_checked"})
        groups = self.rows(output, "duplicate_groups.tsv")
        basename_groups = [row for row in groups if row["group_type"] == "basename"]
        self.assertEqual(len(basename_groups), 1)
        self.assertEqual(int(basename_groups[0]["sample_count"]), 2)
        self.assertFalse(any(row["group_type"] == "ordered_file_dna" for row in groups))

    def test_identical_dna_with_conflicting_labels_is_retained_and_flagged(self):
        dataset = self.dataset([("copy_a", "family_a"), ("copy_b", "family_b")])
        self.gbk(dataset, "family_a", "copy_a", [self.record(record_id="source_a")])
        self.gbk(dataset, "family_b", "copy_b", [self.record(record_id="source_b")])
        output = self.run_inventory()

        samples = self.rows(output, "samples.tsv")
        self.assertEqual(len(samples), 2)
        self.assertEqual({tuple(json.loads(row["tsv_labels_raw_json"])) for row in samples},
                         {("family_a",), ("family_b",)})
        self.assertEqual(len({row["source_sha256"] for row in samples}), 2)
        groups = [row for row in self.rows(output, "duplicate_groups.tsv")
                  if row["group_type"] == "ordered_file_dna"]
        self.assertEqual(len(groups), 1)
        self.assertEqual(int(groups[0]["sample_count"]), 2)
        self.assertTrue(self.is_true(groups[0]["within_dataset_label_conflict"]))

    def test_nonbreaking_space_and_tsv_source_lines_are_preserved(self):
        raw_label = "holomycin\u00a0"
        dataset = self.dataset([("one", raw_label)])
        self.gbk(dataset, raw_label, "one")
        output = self.run_inventory()

        label_rows = self.rows(output, "label_rows.tsv")
        self.assertEqual(len(label_rows), 1)
        self.assertEqual(label_rows[0]["BGC_raw"], "one")
        self.assertEqual(label_rows[0]["GCF_raw"], raw_label)
        self.assertEqual(int(label_rows[0]["line_number"]), 2)
        sample = self.rows(output, "samples.tsv")[0]
        self.assertEqual(sample["folder_gcf_raw"], raw_label)
        self.assertEqual(sample["mapping_status"], "exact_path_label")
        self.assertEqual(json.loads(sample["tsv_labels_raw_json"]), [raw_label])
        self.assertEqual(json.loads(sample["tsv_row_ids_json"]), [label_rows[0]["row_id"]])

    def test_unique_basename_folder_mismatch_preserves_both_labels(self):
        dataset = self.dataset([("one", "other_11")])
        self.gbk(dataset, "other_11.2", "one")
        output = self.run_inventory()

        sample = self.rows(output, "samples.tsv")[0]
        self.assertEqual(sample["mapping_status"], "unique_basename_folder_mismatch")
        self.assertEqual(sample["folder_gcf_raw"], "other_11.2")
        self.assertEqual(json.loads(sample["tsv_labels_raw_json"]), ["other_11"])
        self.assertTrue(self.rows(output, "issues.tsv"))

    def test_multirecord_file_keeps_record_boundaries(self):
        dataset = self.dataset([("multi", "family"), ("joined", "family")])
        first = "ATGAAAGCTTAA"
        second = "ATGCCCGGGAAATAA"
        self.gbk(dataset, "family", "multi", [
            self.record(first, "part_1"), self.record(second, "part_2"),
        ])
        self.gbk(dataset, "family", "joined", [self.record(first + second, "joined")])
        output = self.run_inventory()

        samples = self.rows(output, "samples.tsv")
        multi = next(row for row in samples if row["archive_relpath"].endswith("/multi.gbk"))
        self.assertEqual(int(multi["record_count"]), 2)
        self.assertEqual(multi["genbank_status"], "ok")
        records = [row for row in self.rows(output, "records.tsv")
                   if row["sample_id"] == multi["sample_id"]]
        records.sort(key=lambda row: int(row["record_index"]))
        self.assertEqual([int(row["length_bp"]) for row in records], [len(first), len(second)])
        self.assertEqual(len({row["record_uid"] for row in records}), 2)
        self.assertFalse(any(row["group_type"] == "ordered_file_dna"
                             for row in self.rows(output, "duplicate_groups.tsv")))

    def test_dna_longer_than_40000_bp_is_included_without_truncation(self):
        dataset = self.dataset([("long", "family")])
        sequence = "ATG" + "A" * 39998
        self.gbk(dataset, "family", "long", [self.record(sequence)])
        output = self.run_inventory()

        sample = self.rows(output, "samples.tsv")[0]
        self.assertEqual(sample["genbank_status"], "ok")
        record = self.rows(output, "records.tsv")[0]
        self.assertEqual(int(record["length_bp"]), 40001)
        self.assertEqual(record["sequence_sha256"], hashlib.sha256(sequence.encode()).hexdigest())

    def test_reverse_complements_share_only_orientation_invariant_dna_groups(self):
        dataset = self.dataset([("forward", "family"), ("reverse", "family")])
        sequence = "ATGAAAGCTTAA"
        reverse = str(Seq(sequence).reverse_complement())
        self.assertNotEqual(sequence, reverse)
        self.gbk(dataset, "family", "forward", [self.record(sequence)])
        self.gbk(dataset, "family", "reverse", [self.record(reverse)])
        output = self.run_inventory()

        records = self.rows(output, "records.tsv")
        self.assertEqual(len({row["sequence_sha256"] for row in records}), 2)
        self.assertEqual(len({row["sequence_rc_canonical_sha256"] for row in records}), 1)
        groups = self.rows(output, "duplicate_groups.tsv")
        group_types = {row["group_type"] for row in groups}
        self.assertTrue({"single_record_dna_rc", "record_dna_rc"}.issubset(group_types))
        self.assertNotIn("record_dna", group_types)
        self.assertNotIn("ordered_file_dna", group_types)

    def test_malformed_gbk_remains_in_file_and_label_inventory(self):
        dataset = self.dataset([("broken", "family")])
        path = dataset / "GCF_assignments" / "family" / "broken.gbk"
        path.parent.mkdir(parents=True)
        path.write_text("This is not a GenBank record.\n", encoding="utf-8")
        output = self.run_inventory()

        samples = self.rows(output, "samples.tsv")
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["genbank_status"], "error")
        self.assertEqual(int(samples[0]["record_count"]), 0)
        self.assertEqual(samples[0]["source_sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
        self.assertEqual(len(self.rows(output, "label_rows.tsv")), 1)
        self.assertFalse(self.rows(output, "records.tsv"))
        self.assertTrue(self.rows(output, "issues.tsv"))
        self.assertGreaterEqual(len(self.rows(output, "source_files.tsv")), 2)

    def test_repeated_tsv_rows_are_not_overwritten_or_arbitrarily_selected(self):
        dataset = self.dataset([("one", "family"), ("one", "family")])
        self.gbk(dataset, "family", "one")
        output = self.run_inventory()

        labels = self.rows(output, "label_rows.tsv")
        self.assertEqual(len(labels), 2)
        self.assertEqual({int(row["line_number"]) for row in labels}, {2, 3})
        self.assertEqual(len({row["row_id"] for row in labels}), 2)
        sample = self.rows(output, "samples.tsv")[0]
        self.assertEqual(sample["mapping_status"], "ambiguous")
        self.assertEqual(set(json.loads(sample["tsv_row_ids_json"])),
                         {row["row_id"] for row in labels})

    def test_dataset_local_labels_do_not_create_cross_dataset_conflicts(self):
        dataset_a = self.dataset([("one", "family_a")], name="Dataset_A_Test")
        dataset_b = self.dataset([("two", "family_b")], name="Dataset_B_Test")
        self.gbk(dataset_a, "family_a", "one")
        self.gbk(dataset_b, "family_b", "two")
        output = self.run_inventory()

        groups = [row for row in self.rows(output, "duplicate_groups.tsv")
                  if row["group_type"] == "ordered_file_dna"]
        self.assertEqual(len(groups), 1)
        self.assertEqual(int(groups[0]["sample_count"]), 2)
        self.assertFalse(self.is_true(groups[0]["within_dataset_label_conflict"]))
        self.assertEqual(len(self.rows(output, "dataset_summary.tsv")), 2)

    def test_sample_id_survives_root_move_and_content_change(self):
        dataset = self.dataset([("one", "family")])
        self.gbk(dataset, "family", "one")
        output = self.run_inventory()
        initial = self.rows(output, "samples.tsv")[0]

        moved_source = self.workspace / "other_raw_prefix" / "curated_datasets"
        shutil.copytree(self.source, moved_source)
        moved_output = self.run_inventory(source=moved_source, output=self.workspace / "moved_output")
        moved = self.rows(moved_output, "samples.tsv")[0]
        self.assertEqual(initial["sample_id"], moved["sample_id"])
        self.assertEqual(initial["archive_relpath"], moved["archive_relpath"])
        self.assertEqual(initial["source_sha256"], moved["source_sha256"])

        changed_path = next(moved_source.rglob("*.gbk"))
        SeqIO.write([self.record("ATGCCCGGGTAA")], changed_path, "genbank")
        changed_output = self.run_inventory(source=moved_source, output=self.workspace / "changed_output")
        changed = self.rows(changed_output, "samples.tsv")[0]
        self.assertEqual(initial["sample_id"], changed["sample_id"])
        self.assertNotEqual(initial["source_sha256"], changed["source_sha256"])

    def test_existing_output_is_rejected_without_overwriting_artifacts(self):
        dataset = self.dataset([("one", "family")])
        self.gbk(dataset, "family", "one")
        self.run_inventory()
        before = self.source_hashes(self.output)
        with self.assertRaises((ValueError, FileExistsError)):
            build_inventory(self.source, self.output, self.workspace)
        self.assertEqual(before, self.source_hashes(self.output))

    def test_cpu_parser_checks_preserve_missing_translation_and_multirecord_inputs(self):
        encoder_project = Path(os.environ.get(
            "BGC_ENCODER_PROJECT", str(PROJECT_ROOT.parent / "BGC_annotation")
        ))
        if not (encoder_project / "data_preparation" / "feature_inputs.py").is_file():
            self.skipTest("Adjacent encoder project absent; set BGC_ENCODER_PROJECT for CPU parser checks")
        dataset = self.dataset([("good", "family"), ("missing", "family"), ("multi", "family")])
        self.gbk(dataset, "family", "good")
        self.gbk(dataset, "family", "missing", [self.record(translation=None)])
        self.gbk(dataset, "family", "multi", [self.record(record_id="a"), self.record(record_id="b")])
        output = self.run_inventory(encoder_project=encoder_project)

        samples = {Path(row["archive_relpath"]).stem: row
                   for row in self.rows(output, "samples.tsv")}
        self.assertEqual(len(samples), 3)
        for row in samples.values():
            self.assertEqual(row["genbank_status"], "ok")
        self.assertEqual(samples["good"]["dna_parser_status"], "ok")
        self.assertEqual(samples["good"]["protein_parser_status"], "ok")
        self.assertEqual(samples["missing"]["dna_parser_status"], "ok")
        self.assertEqual(samples["missing"]["protein_parser_status"], "error")
        self.assertEqual(samples["multi"]["dna_parser_status"], "error")
        self.assertEqual(samples["multi"]["protein_parser_status"], "error")
        missing_records = [row for row in self.rows(output, "records.tsv")
                           if row["sample_id"] == samples["missing"]["sample_id"]]
        self.assertEqual(len(missing_records), 1)
        self.assertEqual(int(missing_records[0]["cds_count"]), 1)
        self.assertEqual(int(missing_records[0]["cds_missing_translation_count"]), 1)

    def test_terminal_and_internal_stops_remain_distinct(self):
        dataset = self.dataset([("terminal", "family"), ("internal", "family"), ("both", "family")])
        for name, translation in (("terminal", "MKA*"), ("internal", "MK*A"), ("both", "MKA**")):
            self.gbk(dataset, "family", name, [self.record(record_id=name, translation=translation)])
        output = self.run_inventory()
        records = {row["record_id_raw"]: row for row in self.rows(output, "records.tsv")}
        self.assertEqual(records["terminal"]["cds_terminal_stop_count"], "1")
        self.assertEqual(records["terminal"]["cds_internal_stop_count"], "0")
        self.assertEqual(records["internal"]["cds_terminal_stop_count"], "0")
        self.assertEqual(records["internal"]["cds_internal_stop_count"], "1")
        self.assertEqual(records["both"]["cds_terminal_stop_count"], "1")
        self.assertEqual(records["both"]["cds_internal_stop_count"], "1")

    def test_failed_parser_metadata_does_not_create_protein_duplicate_groups(self):
        dataset = self.dataset([("first", "family"), ("second", "family")])
        for name in ("first", "second"):
            self.gbk(dataset, "family", name)
        def incompatible_parser(path):
            return {"tokens": [{"sequence": "MKA"}], "parents": [{}]}
        with patch("scripts.build_inventory.load_parsers",
                   return_value=({"protein": incompatible_parser}, {"enabled": False})):
            output = self.run_inventory()
        for row in self.rows(output, "samples.tsv"):
            self.assertEqual(row["protein_parser_status"], "error")
            self.assertEqual(row["protein_input_sha256"], "")
            self.assertEqual(row["protein_token_count"], "")
        self.assertFalse(any(row["group_type"] == "protein_input"
                             for row in self.rows(output, "duplicate_groups.tsv")))


if __name__ == "__main__":
    unittest.main()
