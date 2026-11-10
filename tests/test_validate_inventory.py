"""Verify that the independent validator detects tampering and broken references."""

import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.build_inventory import build_inventory
from scripts.prepare_raw_layout import prepare_layout
from scripts.validate_inventory import verify_inventory


class InventoryValidationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="inventory-validation-")
        self.addCleanup(temporary.cleanup)
        self.workspace = Path(temporary.name)
        source = self.workspace / "data/curated_datasets"
        dataset = source / "Dataset_A_Test"
        family = dataset / "GCF_assignments/family_a"
        family.mkdir(parents=True)
        (source / "README.md").write_text("Original source reference\n")
        (dataset / "curated_GCF_assignments.tsv").write_text("BGC\tGCF\none\tfamily_a\nbroken\tfamily_a\n")
        record = SeqRecord(Seq("ATGACCTAA"), id="record_1", name="record_1", description="fixture")
        record.annotations["molecule_type"] = "DNA"
        SeqIO.write([record, record], family / "one.gbk", "genbank")
        (family / "broken.gbk").write_text("Invalid GenBank content\n")
        prepare_layout(self.workspace)
        self.source = self.workspace / "data/raw/curated_datasets"
        self.output = self.workspace / "data/processed/inventory_v1"
        build_inventory(self.source, self.output, self.workspace)

    def metadata(self):
        return json.loads((self.output / "metadata.json").read_text())

    def update_metadata(self, metadata):
        (self.output / "metadata.json").write_text(json.dumps(metadata))

    def test_valid_inventory_accepts_retained_parse_failures_and_multirecord(self):
        counts = verify_inventory(self.workspace, self.output)
        self.assertEqual(counts["gbk_files"], 2)
        self.assertEqual(counts["records"], 2)
        self.assertEqual(counts["genbank_parse_ok"], 1)
        self.assertEqual(counts["multi_record_files"], 1)

    def test_output_tamper_is_detected(self):
        with (self.output / "summary.md").open("a") as handle:
            handle.write("Unexpected edit\n")
        with self.assertRaisesRegex(ValueError, "Output SHA256"):
            verify_inventory(self.workspace, self.output)

    def test_raw_tamper_is_detected(self):
        with (self.source / "README.md").open("a") as handle:
            handle.write("Unexpected edit\n")
        with self.assertRaisesRegex(ValueError, "Raw file size|Raw SHA256"):
            verify_inventory(self.workspace, self.output)

    def test_unknown_record_reference_is_detected_even_with_updated_checksum(self):
        path = self.output / "records.tsv"
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields, data = reader.fieldnames, list(reader)
        data[0]["sample_id"] = "not_a_source_sample"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(data)
        metadata = self.metadata()
        metadata["output_sha256"]["records.tsv"] = hashlib.sha256(path.read_bytes()).hexdigest()
        self.update_metadata(metadata)
        with self.assertRaisesRegex(ValueError, "unknown sample"):
            verify_inventory(self.workspace, self.output)

    def test_metadata_count_mismatch_is_detected(self):
        metadata = self.metadata()
        metadata["counts"]["gbk_files"] += 1
        self.update_metadata(metadata)
        with self.assertRaisesRegex(ValueError, "Metadata counts"):
            verify_inventory(self.workspace, self.output)


if __name__ == "__main__":
    unittest.main()
