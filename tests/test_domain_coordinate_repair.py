"""Check exact peptide proof, both strands, and narrowly scoped byte patches."""

from copy import deepcopy
import importlib.util
from io import StringIO
from pathlib import Path
import sys
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from domain_coordinate_repair import apply_byte_repairs, plan_repairs, GENOMIC_ONLY_RULE

spec = importlib.util.spec_from_file_location("repair_test_encoder_parser", ROOT.parent / "BGC_annotation/data_preparation/feature_inputs.py")
protein = importlib.util.module_from_spec(spec)
spec.loader.exec_module(protein)


class DomainCoordinateRepairTests(unittest.TestCase):
    def fixture(self, strand=1, missing_coordinates=False, valid=False, feature_type="aSDomain"):
        coding = Seq("ATGGCTTGTGATGAATTTGGT")  # MACDEFG, domain CDE = aa [2:5)
        sequence = Seq("AAA") + (coding if strand == 1 else coding.reverse_complement()) + Seq("CCC")
        record = SeqRecord(sequence, id="fixture", name="fixture", description="unaltered annotation")
        record.annotations["molecule_type"] = "DNA"
        cds = SeqFeature(SimpleLocation(3, 24, strand=strand), type="CDS",
                         qualifiers={"locus_tag": ["gene1"], "translation": ["MACDEFG"], "transl_table": ["11"]})
        coordinates = (9, 18 if valid else 15) if strand == 1 else (9 if valid else 12, 18)
        qualifiers = {"locus_tag": ["gene1"], "aSDomain": ["test_domain"], "translation": ["CDE"]}
        if not missing_coordinates:
            qualifiers.update(protein_start=["2"], protein_end=["5" if valid else "4"])
        domain = SeqFeature(SimpleLocation(*coordinates, strand=strand), type=feature_type, qualifiers=qualifiers)
        record.features = [cds, domain]
        return record

    def serialize(self, record, crlf=False):
        stream = StringIO()
        SeqIO.write(record, stream, "genbank")
        raw = stream.getvalue().encode("utf-8")
        return raw.replace(b"\n", b"\r\n") if crlf else raw

    def test_positive_joint_repair_changes_only_two_numeric_tokens(self):
        record = self.fixture()
        raw = self.serialize(record)
        before = deepcopy(record)
        changes, checks = plan_repairs(record, protein)
        self.assertEqual(len(changes), 1)
        change = changes[0]
        self.assertEqual(change["new_protein_end"], 5)
        self.assertEqual(change["new_genomic_end"], 18)
        self.assertEqual(change["domain_sha256"], change["repaired_dna_translation_sha256"])
        fixed = apply_byte_repairs(raw, record, changes)
        expected = raw.replace(b"aSDomain        10..15", b"aSDomain        10..18").replace(b'/protein_end="4"', b'/protein_end="5"')
        self.assertEqual(fixed, expected)
        self.assertEqual(str(record.features[1].location), str(before.features[1].location))
        self.assertEqual(record.features[1].qualifiers, before.features[1].qualifiers)
        parsed = SeqIO.read(StringIO(fixed.decode()), "genbank")
        self.assertEqual(parsed.seq, record.seq)
        self.assertEqual(str(parsed.features[1].extract(parsed.seq).translate(table=11)), "CDE")
        self.assertEqual(checks[0]["status"], "repaired")

    def test_negative_strand_uses_start_minus_three_and_preserves_crlf(self):
        record = self.fixture(strand=-1)
        raw = self.serialize(record, crlf=True)
        changes, _ = plan_repairs(record, protein)
        fixed = apply_byte_repairs(raw, record, changes)
        expected = raw.replace(b"complement(13..18)", b"complement(10..18)").replace(b'/protein_end="4"', b'/protein_end="5"')
        self.assertEqual(fixed, expected)
        self.assertEqual(fixed.count(b"\r\n"), raw.count(b"\r\n"))
        self.assertEqual(changes[0]["new_genomic_start"], 9)

    def test_missing_coordinates_only_repairs_genomic_token_and_infers_unique_slice(self):
        record = self.fixture(missing_coordinates=True)
        raw = self.serialize(record)
        changes, _ = plan_repairs(record, protein)
        self.assertEqual(changes[0]["rule"], GENOMIC_ONLY_RULE)
        self.assertIsNone(changes[0]["new_protein_end"])
        self.assertEqual(changes[0]["resolved_protein_end"], 5)
        fixed = apply_byte_repairs(raw, record, changes)
        self.assertEqual(fixed, raw.replace(b"aSDomain        10..15", b"aSDomain        10..18"))
        self.assertNotIn(b"/protein_end", fixed)

    def test_valid_domains_unchanged_and_pfam_type_preserved_when_repaired(self):
        record = self.fixture(valid=True)
        raw = self.serialize(record)
        changes, checks = plan_repairs(record, protein)
        self.assertEqual(changes, [])
        self.assertEqual(checks[0]["status"], "unchanged")
        self.assertEqual(apply_byte_repairs(raw, record, changes), raw)
        record = self.fixture(feature_type="PFAM_domain")
        changes, _ = plan_repairs(record, protein)
        self.assertEqual(changes[0]["feature_type"], "PFAM_domain")
        fixed = apply_byte_repairs(self.serialize(record), record, changes)
        parsed = SeqIO.read(StringIO(fixed.decode()), "genbank")
        self.assertEqual(parsed.features[1].type, "PFAM_domain")

    def test_ambiguous_parent_rejected_even_when_alias_identifies_one(self):
        record = self.fixture()
        other = deepcopy(record.features[0])
        other.qualifiers["locus_tag"] = ["gene2"]
        record.features.append(other)
        with self.assertRaisesRegex(ValueError, "unique geometric"):
            plan_repairs(record, protein)

    def test_absent_parent_pfam_logged_but_asdomain_fails(self):
        record = self.fixture(feature_type="PFAM_domain")
        record.features[1].location = SimpleLocation(0, 3, strand=1)
        changes, checks = plan_repairs(record, protein)
        self.assertEqual(changes, [])
        self.assertEqual(checks[0]["status"], "needs_review")
        record.features[1].type = "aSDomain"
        with self.assertRaisesRegex(ValueError, "expected one parent"):
            plan_repairs(record, protein)

    def test_incorrect_parent_slice_or_dna_evidence_is_rejected(self):
        record = self.fixture()
        record.features[1].qualifiers["translation"] = ["CDQ"]
        with self.assertRaisesRegex(ValueError, "parent-sequence"):
            plan_repairs(record, protein)
        record = self.fixture()
        record.seq = record.seq[:15] + Seq("CAA") + record.seq[18:]
        with self.assertRaisesRegex(ValueError, "genomic translation"):
            plan_repairs(record, protein)

    def test_unquoted_end_supported_but_multiline_location_or_mismatched_plan_rejected(self):
        record = self.fixture()
        raw = self.serialize(record).replace(b'/protein_end="4"', b'/protein_end=4')
        changes, _ = plan_repairs(record, protein)
        self.assertIn(b"/protein_end=5", apply_byte_repairs(raw, record, changes))
        with self.assertRaisesRegex(ValueError, "location text"):
            apply_byte_repairs(raw.replace(b"aSDomain        10..15", b"aSDomain        10..\n                     15"), record, changes)
        changes[0]["new_genomic_end"] = 19
        with self.assertRaisesRegex(ValueError, "genomic extension"):
            apply_byte_repairs(raw, record, changes)


if __name__ == "__main__":
    unittest.main()
