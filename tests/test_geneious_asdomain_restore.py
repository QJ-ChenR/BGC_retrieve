"""CPU-only reconstruction evidence and minimal-byte-edit regression tests."""

from copy import deepcopy
from io import StringIO
from pathlib import Path
import sys
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from geneious_asdomain_restore import apply_reconstruction_bytes, plan_reconstruction


def encode_record(record, newline="\n"):
    stream = StringIO()
    SeqIO.write(record, stream, "genbank")
    raw = stream.getvalue().replace("\n", newline).encode()
    return raw, SeqIO.read(StringIO(raw.decode()), "genbank")


def fixture(strand=1, codon_start=1, dna="ATGGCTGAATTTGGTTAA", peptide="MAEFG", table=1):
    coding = "T" * (codon_start - 1) + dna
    sequence = Seq(coding)
    if strand == -1:
        sequence = sequence.reverse_complement()
    parent = SeqFeature(SimpleLocation(0, len(sequence), strand=strand), type="CDS", qualifiers={
        "locus_tag": ["correct_parent"], "translation": [peptide],
        "codon_start": [str(codon_start)], "transl_table": [str(table)],
    })
    start, end = codon_start - 1 + 3, codon_start - 1 + 12
    if strand == -1:
        start, end = len(sequence) - end, len(sequence) - start
    domain = SeqFeature(SimpleLocation(start, end, strand=strand), type="misc_feature", qualifiers={
        "note": ["Geneious type: aSDomain", "original evidence"],
        "locus_tag": ["correct_parent"], "standard_name": ["original_name"],
    })
    record = SeqRecord(sequence, id="TEST.1", name="TEST", description="evidence fixture",
                       annotations={"molecule_type": "DNA", "comment": "Preserve original comment."},
                       features=[parent, domain])
    return record


class ReconstructionTests(unittest.TestCase):
    def test_positive_and_negative_strand_with_codon_start(self):
        for strand in (-1, 1):
            for codon_start in (1, 2, 3):
                with self.subTest(strand=strand, codon_start=codon_start):
                    record = fixture(strand, codon_start)
                    old = deepcopy(record)
                    changes, checks = plan_reconstruction(record)
                    self.assertEqual(len(changes), 1)
                    change = changes[0]
                    self.assertEqual((change["protein_start"], change["protein_end"]), (1, 4))
                    self.assertEqual(change["translation"], "AEF")
                    self.assertEqual(change["status"], "restorable")
                    self.assertEqual(change["new_qualifiers"]["note"], old.features[1].qualifiers["note"])
                    self.assertNotIn("aSDomain", change["new_qualifiers"])
                    self.assertNotIn("score", change["new_qualifiers"])
                    self.assertEqual(record.features, old.features)
                    self.assertEqual(checks[0]["dna_translation_matches"], True)

    def test_translation_table_is_used(self):
        changes, _ = plan_reconstruction(fixture(dna="ATGGCTTGATTTGGTTAA", peptide="MAWFG", table=4))
        self.assertEqual(changes[0]["translation"], "AWF")

    def test_locus_conflict_maps_to_sequence_verified_parent(self):
        record = fixture()
        record.features[1].qualifiers["locus_tag"] = ["stale_parent"]
        changes, _ = plan_reconstruction(record)
        mapping = changes[0]["identifier_mapping"]
        self.assertEqual(mapping["old_locus_tag"], "stale_parent")
        self.assertEqual(mapping["resolved_locus_tag"], "correct_parent")
        self.assertTrue(mapping["changed"])
        self.assertEqual(changes[0]["new_qualifiers"]["locus_tag"], ["correct_parent"])
        self.assertEqual(record.features[1].qualifiers["locus_tag"], ["stale_parent"])

    def test_only_exact_geneious_marker_is_restored(self):
        record = fixture()
        for marker in ("Geneious type: PFAM_domain", "Geneious type: aSModule",
                       "Geneious type: aSDomain extra", "some Geneious type: aSDomain"):
            other = deepcopy(record.features[1])
            other.qualifiers["note"] = [marker]
            record.features.append(other)
        recognized = deepcopy(record.features[1])
        recognized.type = "aSDomain"
        record.features.append(recognized)
        changes, checks = plan_reconstruction(record)
        self.assertEqual([c["feature_index"] for c in changes], [1])
        self.assertEqual(len(checks), 1)

    def assert_review(self, record, reason):
        changes, checks = plan_reconstruction(record)
        self.assertEqual(changes, [])
        self.assertEqual(checks[0]["status"], "needs_review")
        self.assertIn(reason, checks[0]["reason"])
        return checks

    def test_nontriplet_and_out_of_frame_are_not_rounded(self):
        for start, end in ((3, 13), (4, 13), (4, 12)):
            record = fixture()
            record.features[1].location = SimpleLocation(start, end, strand=1)
            self.assert_review(record, "not complete codons")
            self.assertEqual(record.features[1].location, SimpleLocation(start, end, strand=1))

    def test_parent_sequence_disagreement_is_rejected(self):
        record = fixture()
        record.features[0].qualifiers["translation"] = ["MWWFG"]
        self.assert_review(record, "does not equal")

    def test_stop_codons_are_rejected(self):
        self.assert_review(fixture(dna="ATGGCTTGATTTGGTTAA", peptide="MAWFG"), "stop codons")

    def test_missing_ambiguous_and_wrong_strand_parents(self):
        record = fixture()
        record.features[0].type = "gene"
        self.assert_review(record, "found 0")
        record = fixture()
        record.features.append(deepcopy(record.features[0]))
        self.assert_review(record, "found 2")
        record = fixture()
        record.features[0].location = SimpleLocation(0, len(record.seq), strand=-1)
        self.assert_review(record, "found 0")

    def test_fuzzy_and_compound_locations_need_review(self):
        record = fixture()
        record.features[1].location = SimpleLocation(BeforePosition(3), 12, strand=1)
        self.assert_review(record, "domain location")
        record = fixture()
        record.features[1].location = CompoundLocation([SimpleLocation(3, 6, strand=1),
                                                        SimpleLocation(9, 12, strand=1)])
        self.assert_review(record, "domain location")
        record = fixture()
        record.features[0].location = CompoundLocation([SimpleLocation(0, 6, strand=1),
                                                        SimpleLocation(9, 18, strand=1)])
        self.assert_review(record, "parent location")

    def test_duplicate_annotations_are_retained_and_reported(self):
        record = fixture()
        duplicate = deepcopy(record.features[1])
        del duplicate.qualifiers["standard_name"]
        record.features.append(duplicate)
        raw, parsed = encode_record(record)
        changes, checks = plan_reconstruction(parsed)
        self.assertEqual(len(changes), 2)
        self.assertEqual(checks[0]["duplicate_group"], checks[1]["duplicate_group"])
        self.assertEqual([c["duplicate_group_size"] for c in checks], [2, 2])
        output = apply_reconstruction_bytes(raw, parsed, changes)
        restored = SeqIO.read(StringIO(output.decode()), "genbank")
        self.assertEqual([f.type for f in restored.features], ["CDS", "aSDomain", "aSDomain"])
        self.assertIn("standard_name", restored.features[1].qualifiers)
        self.assertNotIn("standard_name", restored.features[2].qualifiers)

    def test_existing_conflicting_qualifiers_require_review(self):
        for key, value in (("translation", "WWW"), ("protein_start", "2"),
                           ("protein_end", "5"), ("protein_id", "conflicting")):
            with self.subTest(key=key):
                record = fixture()
                record.features[1].qualifiers[key] = [value]
                self.assert_review(record, key)

    def test_invalid_parent_metadata_and_empty_domains_require_review(self):
        for key, values in (("codon_start", ["4"]), ("transl_table", ["not_a_table"]),
                            ("translation", [""]), ("locus_tag", ["one", "two"])):
            with self.subTest(key=key):
                record = fixture()
                record.features[0].qualifiers[key] = values
                changes, checks = plan_reconstruction(record)
                self.assertEqual(changes, [])
                self.assertEqual(checks[0]["status"], "needs_review")
        record = fixture()
        record.features[1].location = SimpleLocation(3, 3, strand=1)
        self.assert_review(record, "empty")

    def test_unrelated_bytes_and_target_notes_are_preserved(self):
        for newline in ("\n", "\r\n"):
            with self.subTest(newline=repr(newline)):
                record = fixture()
                record.features[1].qualifiers["locus_tag"] = ["stale_parent"]
                pfam = deepcopy(record.features[1])
                pfam.qualifiers["note"] = ["Geneious type: PFAM_domain"]
                record.features.append(pfam)
                raw, parsed = encode_record(record, newline)
                changes, _ = plan_reconstruction(parsed)
                output = apply_reconstruction_bytes(raw, parsed, changes)
                original_head = raw[:raw.index(b"     misc_feature")]
                original_tail = raw[raw.index(b"     misc_feature", raw.index(b"     misc_feature") + 1):]
                self.assertTrue(output.startswith(original_head))
                self.assertTrue(output.endswith(original_tail))
                self.assertIn(b'/standard_name="original_name"' + newline.encode(), output)
                self.assertIn(b'/note="original evidence"' + newline.encode(), output)
                self.assertEqual(output.count(b"\r\n"), output.count(b"\n") if newline == "\r\n" else 0)
                restored = SeqIO.read(StringIO(output.decode()), "genbank")
                self.assertEqual(restored.features[1].location, parsed.features[1].location)
                self.assertEqual(restored.features[1].qualifiers["locus_tag"], ["correct_parent"])
                self.assertEqual(restored.features[2], parsed.features[2])

    def test_idempotence_and_failed_candidates_remain_unchanged(self):
        record = fixture()
        bad = deepcopy(record.features[1])
        bad.location = SimpleLocation(4, 13, strand=1)
        record.features.append(bad)
        raw, parsed = encode_record(record)
        changes, _ = plan_reconstruction(parsed)
        output = apply_reconstruction_bytes(raw, parsed, changes)
        restored = SeqIO.read(StringIO(output.decode()), "genbank")
        next_changes, next_checks = plan_reconstruction(restored)
        self.assertEqual(next_changes, [])
        self.assertEqual(len(next_checks), 1)
        self.assertEqual(next_checks[0]["feature_index"], 2)
        self.assertEqual(next_checks[0]["status"], "needs_review")
        self.assertEqual(restored.features[2], parsed.features[2])
        self.assertEqual(apply_reconstruction_bytes(output, restored, next_changes), output)

    def test_stale_or_tampered_plans_and_source_records_are_rejected(self):
        raw, parsed = encode_record(fixture())
        changes, _ = plan_reconstruction(parsed)
        altered = deepcopy(changes)
        altered[0]["new_qualifiers"]["translation"] = ["WWW"]
        with self.assertRaisesRegex(ValueError, "unproven"):
            apply_reconstruction_bytes(raw, parsed, altered)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            apply_reconstruction_bytes(raw, parsed, changes * 2)
        altered_record = deepcopy(parsed)
        altered_record.features[0].qualifiers["translation"] = ["MWWFG"]
        with self.assertRaisesRegex(ValueError, "differs"):
            apply_reconstruction_bytes(raw, altered_record, changes)

    def test_long_translations_are_wrapped_and_round_trip(self):
        record = fixture(dna="ATG" + "GCT" * 250 + "TAA", peptide="M" + "A" * 250)
        record.features[1].location = SimpleLocation(3, 753, strand=1)
        raw, parsed = encode_record(record)
        changes, _ = plan_reconstruction(parsed)
        output = apply_reconstruction_bytes(raw, parsed, changes)
        restored = SeqIO.read(StringIO(output.decode()), "genbank")
        self.assertEqual(restored.features[1].qualifiers["translation"], ["A" * 250])
        self.assertTrue(all(len(line) <= 80 for line in output.splitlines()))


if __name__ == "__main__":
    unittest.main()
