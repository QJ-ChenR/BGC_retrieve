"""Publication/provenance regression tests with small in-memory CPU parser fixtures."""
from argparse import Namespace
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from Bio.SeqFeature import SimpleLocation

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import reconstruct_h as workflow
from build_inventory import compact, digest, file_sha256, snapshot, write_json, write_tsv
from pilot_features import read_json, read_tsv
from test_geneious_asdomain_restore import encode_record, fixture


def protein_parser(path):
    record, _ = workflow.parsed_record(path.read_bytes())
    tokens = []
    for parent_index, parent in enumerate(record.features):
        if parent.type != "CDS":
            continue
        domains = [(i, f) for i, f in enumerate(record.features) if f.type == "aSDomain"
                   and f.qualifiers.get("locus_tag") == parent.qualifiers.get("locus_tag")]
        for index, feature in domains or [(parent_index, parent)]:
            seq = feature.qualifiers["translation"][0]
            tokens.append({"feature_index": index, "source": "domain" if domains else "protein",
                           "sequence": seq, "protein_start": int(feature.qualifiers.get("protein_start", [0])[0]),
                           "protein_end": int(feature.qualifiers.get("protein_end", [len(seq)])[0]),
                           "parent_identifiers": {"locus_tag": parent.qualifiers["locus_tag"][0]}})
    return {"tokens": tokens}


def dna_parser(path):
    record, _ = workflow.parsed_record(path.read_bytes())
    return str(record.seq), {}


class ReconstructionWorkflowTests(unittest.TestCase):
    def setup_workspace(self, root, unresolved=False):
        raw = root / "data/raw/curated_datasets"
        source = root / "data/processed/repaired_gh_v1"
        inventory = root / "data/processed/inventory_v1"
        rel = f"{workflow.DATASET}/GCF_assignments/A/fixture.gbk"
        record = fixture()
        if unresolved:
            record.features[1].location = SimpleLocation(4, 12, strand=1)
            record.annotations["structured_comment"] = {"Editing": {"Note": "Original editing evidence"}}
        content, _ = encode_record(record)
        for directory in (raw, source):
            (directory / rel).parent.mkdir(parents=True)
            (directory / rel).write_bytes(content)
            write_tsv(directory / workflow.DATASET / "curated_GCF_assignments.tsv", [{"BGC": "fixture", "GCF": "A"}])
        inventory.mkdir()
        for name in ("samples.tsv", "source_files.tsv"):
            (inventory / name).write_text("fixture\n")
        checksum = file_sha256(source / rel)
        base = {"sample_id": "fixture_id", "dataset": workflow.DATASET, "archive_relpath": rel,
                "original_bgc_id": "fixture", "gcf_raw": "A", "tsv_row_ids_json": compact(["row:1"])}
        sample = {**base, "raw_source_sha256": checksum, "raw_source_relpath": f"data/raw/curated_datasets/{rel}",
                  "processed_relpath": rel, "processed_sha256": checksum, "record_count": "1",
                  "protein_parser_status_after": "ok", "dna_parser_status_after": "ok",
                  "dna_unchanged": "True", "translations_unchanged": "True", "unlisted_annotations_unchanged": "True",
                  "protein_input_sha256_after": workflow.token_digest(protein_parser(source / rel)),
                  "dna_sha256": digest(str(record.seq))}
        write_tsv(source / "samples.tsv", [sample])
        for name in ("repairs.tsv", "domain_checks.tsv", "review_items.tsv"):
            write_tsv(source / name, [], ("dataset",))
        write_json(source / "metadata.json", {
            "format": "bgc_gh_coordinate_repairs_v1", "status": "complete", "raw_unchanged": True,
            "artifact_sha256": {p.relative_to(source).as_posix(): file_sha256(p)
                                for p in source.rglob("*") if p.is_file()}})
        args = Namespace(workspace=str(root), source=str(source), output=str(root / "data/processed/restored"),
                         raw=str(raw), inventory=str(inventory), encoder_project=str(root / "encoder"))
        inventory_result = ([{**base, "source_sha256": checksum}], snapshot(raw, root),
                            {name: file_sha256(inventory / name) for name in ("samples.tsv", "source_files.tsv")})
        return args, inventory_result, content

    def run_fixture(self, args, inventory_result):
        parsers = {"protein": protein_parser, "dna": dna_parser}
        with patch.object(workflow, "load_inventory", return_value=inventory_result), \
                patch.object(workflow, "load_parsers", return_value=(parsers, {"source_sha256": {"fixture": "unused"}})), \
                patch.object(workflow, "verify_artifacts", wraps=workflow.verify_artifacts) as verify:
            # A real unchanged parser source is still hash-checked by orchestration.
            encoder = Path(args.encoder_project)
            encoder.mkdir()
            (encoder / "fixture").write_text("CPU parser fixture\n")
            with patch.object(workflow, "load_parsers", return_value=(
                    parsers, {"source_sha256": {"fixture": file_sha256(encoder / "fixture")}})):
                result = workflow.reconstruct(args)
            self.assertGreaterEqual(verify.call_count, 4)
            return result

    def test_complete_publication_preserves_sources_and_hashes_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            args, evidence, content = self.setup_workspace(root)
            source_before = snapshot(Path(args.source), root)
            counts = self.run_fixture(args, evidence)
            output = Path(args.output)
            self.assertEqual(counts["restored_domains"], 1)
            self.assertEqual(counts["status"], "complete")
            self.assertEqual(snapshot(Path(args.source), root), source_before)
            self.assertEqual(snapshot(Path(args.raw), root), evidence[1])
            metadata = read_json(output / "metadata.json")
            self.assertTrue(metadata["ready_for_full_h_encoding"])
            workflow.verify_artifacts(output, metadata["artifact_sha256"])
            sample = read_tsv(output / "samples.tsv")[0]
            self.assertEqual(sample["source_sha256"], sample["raw_source_sha256"])
            self.assertNotEqual(sample["source_sha256"], sample["processed_sha256"])
            self.assertEqual(sample["ready_for_encoding"], "True")
            with self.assertRaisesRegex(ValueError, "new processed directory"):
                workflow.reconstruct(args)

    def test_unresolved_domain_preserved_and_readiness_blocked_despite_parsing(self):
        with tempfile.TemporaryDirectory() as temp:
            args, evidence, content = self.setup_workspace(Path(temp), unresolved=True)
            counts = self.run_fixture(args, evidence)
            output = Path(args.output)
            metadata = read_json(output / "metadata.json")
            self.assertEqual(counts["status"], "needs_review")
            self.assertTrue(metadata["snapshot_complete"])
            self.assertFalse(metadata["ready_for_full_h_encoding"])
            sample = read_tsv(output / "samples.tsv")[0]
            self.assertEqual(sample["protein_parser_status_after"], "ok")
            self.assertEqual(sample["ready_for_encoding"], "False")
            self.assertEqual((output / sample["processed_relpath"]).read_bytes(), content)
            review = read_tsv(output / "review_items.tsv")
            self.assertEqual(len(review), 1)
            self.assertIn("Original editing evidence", review[0]["record_structured_comment_json"])

    def test_source_tampering_refuses_publication(self):
        with tempfile.TemporaryDirectory() as temp:
            args, evidence, content = self.setup_workspace(Path(temp))
            source_file = next(Path(args.source).rglob("*.gbk"))
            source_file.write_bytes(content.replace(b"original_name", b"changed_name"))
            with patch.object(workflow, "load_inventory", return_value=evidence):
                with self.assertRaisesRegex(ValueError, "artifact mismatch"):
                    workflow.reconstruct(args)
            self.assertFalse(Path(args.output).exists())

    def test_label_conflicts_refused_even_when_source_manifest_is_consistent(self):
        with tempfile.TemporaryDirectory() as temp:
            args, evidence, _ = self.setup_workspace(Path(temp))
            wrong_inventory = deepcopy(evidence[0])
            wrong_inventory[0]["gcf_raw"] = "B"
            with self.assertRaisesRegex(ValueError, "identity/label"):
                workflow.validate_source(Path(args.source), wrong_inventory)

    def test_output_cannot_be_raw_or_inside_source(self):
        with tempfile.TemporaryDirectory() as temp:
            args, _, _ = self.setup_workspace(Path(temp))
            for dest in (Path(args.raw) / "new", Path(args.source) / "new", Path(temp) / "data/processed"):
                args.output = str(dest)
                with self.assertRaisesRegex(ValueError, "new processed directory"):
                    workflow.reconstruct(args)

    def test_token_validation_rejects_unrelated_change_and_silent_drop(self):
        token = {"feature_index": 0, "source": "protein", "sequence": "MAEFG", "protein_start": 0,
                 "protein_end": 5, "parent_identifiers": {"locus_tag": "parent"}}
        before = {"tokens": [token]}
        with self.assertRaisesRegex(ValueError, "coverage"):
            workflow.check_token_transition(before, {"tokens": []}, [])
        after = deepcopy(before)
        after["tokens"][0]["sequence"] = "MAEFA"
        with self.assertRaisesRegex(ValueError, "existing token changed"):
            workflow.check_token_transition(before, after, [])


if __name__ == "__main__":
    unittest.main()
