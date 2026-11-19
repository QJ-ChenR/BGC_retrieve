"""Guard raw protection, transactional publication and annotation preservation."""

from copy import deepcopy
from io import StringIO
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import repair_gh as workflow


class RepairGHWorkflowTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gh-workflow-test-")
        self.addCleanup(temporary.cleanup)
        self.workspace = Path(temporary.name)
        self.args = workflow.build_parser().parse_args([
            "--workspace", str(self.workspace),
            "--encoder-project", str(ROOT.parent / "BGC_annotation"),
        ])
        self.output = self.workspace / self.args.output

    def fixture(self):
        record = SeqRecord(Seq("AAAATGGCTTGTGATGAATTTGGTCCC"), id="fixture",
                           name="fixture", description="unchanged description")
        record.annotations["molecule_type"] = "DNA"
        record.features = [
            SeqFeature(SimpleLocation(3, 24, strand=1), type="CDS", qualifiers={
                "locus_tag": ["gene1"], "translation": ["MACDEFG"], "transl_table": ["11"]}),
            SeqFeature(SimpleLocation(9, 15, strand=1), type="aSDomain", qualifiers={
                "locus_tag": ["gene1"], "aSDomain": ["domain"], "translation": ["CDE"],
                "protein_start": ["2"], "protein_end": ["4"]}),
        ]
        stream = StringIO()
        SeqIO.write(record, stream, "genbank")
        return record, stream.getvalue().encode("utf-8")

    def test_existing_output_is_rejected_before_input_processing(self):
        self.output.mkdir(parents=True)
        sentinel = self.output / "keep.txt"
        sentinel.write_bytes(b"keep existing result")
        with mock.patch.object(workflow, "load_inventory") as loader:
            with self.assertRaisesRegex(ValueError, "new subdirectory"):
                workflow.run(self.args)
            loader.assert_not_called()
        self.assertEqual(sentinel.read_bytes(), b"keep existing result")

    def test_raw_or_external_output_is_rejected_before_input_processing(self):
        for output in (Path("data/raw/new-repair"), Path("outside-processed"), Path("data/processed")):
            with self.subTest(output=output):
                self.args.output = output
                with mock.patch.object(workflow, "load_inventory") as loader:
                    with self.assertRaisesRegex(ValueError, "new subdirectory"):
                        workflow.run(self.args)
                    loader.assert_not_called()
                self.assertFalse((self.workspace / output).exists())

    def test_processed_parser_failure_does_not_publish_or_modify_raw(self):
        _, raw_bytes = self.fixture()
        relative = "Dataset_G_Planomonospora/GCF_assignments/family/fixture.gbk"
        raw = self.workspace / self.args.raw_root / relative
        raw.parent.mkdir(parents=True)
        raw.write_bytes(raw_bytes)
        sample = {"dataset": "Dataset_G_Planomonospora", "archive_relpath": relative,
                  "source_sha256": workflow.file_sha256(raw)}
        with mock.patch.object(workflow, "load_inventory", return_value=([sample], [], {})), \
                mock.patch.object(workflow, "parse_protein", side_effect=[
                    (None, "error", "original inclusive endpoint", []),
                    (None, "error", "injected processed-parser failure", []),
                ]) as parser:
            with self.assertRaisesRegex(ValueError, "Processed protein validation failed"):
                workflow.run(self.args)
        self.assertEqual(parser.call_count, 2)
        self.assertEqual(raw.read_bytes(), raw_bytes)
        self.assertFalse(self.output.exists())
        self.assertEqual(list(self.output.parent.glob(".gh-repair-*")), [])

    def test_preservation_check_rejects_unlisted_dna_change(self):
        original, _ = self.fixture()
        changed = deepcopy(original)
        changed.seq = Seq("C") + changed.seq[1:]
        with self.assertRaisesRegex(ValueError, "record seq"):
            workflow.check_preserved(original, changed, [])

    def test_preservation_check_rejects_translation_and_other_qualifier_changes(self):
        original, _ = self.fixture()
        workflow.check_preserved(original, deepcopy(original), [])
        for index, qualifier, value in ((0, "translation", "MACDEFA"),
                                        (1, "translation", "CDF"),
                                        (0, "product", "unrequested annotation")):
            with self.subTest(index=index, qualifier=qualifier):
                changed = deepcopy(original)
                changed.features[index].qualifiers[qualifier] = [value]
                with self.assertRaisesRegex(ValueError, "Unexpected annotation change"):
                    workflow.check_preserved(original, changed, [])


if __name__ == "__main__":
    unittest.main()
