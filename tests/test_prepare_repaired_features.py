"""Portable G/H inputs and provenance guards; actual CPU parsers, no model loading."""

import argparse
from copy import deepcopy
import io
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
import pilot_features as pilot
import prepare_repaired_features as repaired
import build_inventory


class RepairedFeatureInputsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gh-input-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "repaired"
        self.raw = self.root / "raw"
        self.output = self.root / "prepared"
        self.encoder = ROOT.parent / "BGC_annotation"
        self.datasets = {"G": ("Dataset_G_Planomonospora", 1),
                         "H": ("Dataset_H_Reclassification", 1)}
        patch = mock.patch.dict(pilot.DATASETS, self.datasets, clear=True)
        patch.start()
        self.addCleanup(patch.stop)
        self.samples = []
        self.repairs = []
        self.parsers, _ = build_inventory.load_parsers(self.encoder)
        for index, (key, (dataset, _)) in enumerate(self.datasets.items()):
            # G models a coordinate-repaired record; H is byte-identical to raw.
            family = "family\u00a0" if index == 0 else "阴性单例"
            original = "original_" + key
            relative = f"{dataset}/GCF_assignments/{family}/{original}.gbk"
            sample_id = "bgc_" + pilot.hashlib.sha256(
                (pilot.INVENTORY_SCHEMA + "\0" + relative).encode("utf-8")).hexdigest()[:24]
            source = self.source / relative
            raw = self.raw / relative
            source.parent.mkdir(parents=True)
            raw.parent.mkdir(parents=True)
            record = self.record(original)
            self.write_record(source, record)
            before = deepcopy(record)
            if index == 0:
                before.features[1].location = SimpleLocation(9, 15, strand=1)
                before.features[1].qualifiers["protein_end"] = ["4"]
                self.repairs.append({"sample_id": sample_id, "dataset": dataset, "feature_index": "1"})
            self.write_record(raw, before)
            protein = self.parsers["protein"](source)
            dna, _ = self.parsers["dna"](source)
            sample = {
                "dataset": dataset, "sample_id": sample_id, "original_bgc_id": original,
                "archive_relpath": relative, "gcf_raw": family,
                "tsv_row_ids_json": repaired._compact([f"{dataset}:tsv:2"]),
                "raw_source_relpath": "raw/" + relative,
                "raw_source_sha256": pilot.sha256(raw), "processed_relpath": relative,
                "processed_sha256": pilot.sha256(source), "record_count": "1",
                "repaired_domain_count": "1" if index == 0 else "0", "length_bp": str(len(dna)),
                "dna_sha256": repaired._digest(dna),
                "dna_unchanged": "True", "translations_unchanged": "True",
                "unlisted_annotations_unchanged": "True", "protein_parser_status_after": "ok",
                "dna_parser_status_after": "ok",
                "protein_input_sha256_after": repaired._digest(repaired._compact(
                    [t["sequence"] for t in protein["tokens"]])),
                "protein_token_count": str(len(protein["tokens"])),
            }
            self.samples.append(sample)
            labels = [{"BGC": original, "GCF": family}]
            pilot.write_tsv(self.source / dataset / "curated_GCF_assignments.tsv", labels, ("BGC", "GCF"))
            pilot.write_tsv(self.raw / dataset / "curated_GCF_assignments.tsv", labels, ("BGC", "GCF"))
        pilot.write_tsv(self.source / "samples.tsv", self.samples, tuple(self.samples[0]))
        pilot.write_tsv(self.source / "repairs.tsv", self.repairs, tuple(self.repairs[0]))
        pilot.write_tsv(self.source / "review_items.tsv", [], ("sample_id", "dataset", "reason"))
        # Extra repair artifact is verified even though it is not copied into pilot inputs.
        (self.source / "domain_checks.tsv").write_text("status\nok\n", encoding="utf-8")
        self.metadata = {
            "format": repaired.REPAIR_SCHEMA, "status": "complete", "raw_root": str(self.raw),
            "raw_unchanged": True,
            "counts": [{"dataset": name, "sample_count": count} for name, count in self.datasets.values()],
            "raw_label_sha256": {name + "/curated_GCF_assignments.tsv": pilot.sha256(
                self.source / name / "curated_GCF_assignments.tsv") for name, _ in self.datasets.values()},
            "artifact_sha256": {p.relative_to(self.source).as_posix(): pilot.sha256(p)
                                for p in self.source.rglob("*") if p.is_file()},
        }
        pilot.write_json(self.source / "metadata.json", self.metadata)
        self.args = argparse.Namespace(repaired_root=self.source, output=self.output,
                                       datasets=["G", "H"], encoder_project=self.encoder)

    @staticmethod
    def record(name):
        # One 3 aa aSDomain plus a full 2043 aa CDS exercise the 1021 aa window policy.
        sequence = "AAAATGGCTTGTGATGAATTTGGTCCC" + "GCT" * 2043 + "NRY"
        record = SeqRecord(Seq(sequence), id=name, name=name, description="fixture")
        record.annotations["molecule_type"] = "DNA"
        record.features = [
            SeqFeature(SimpleLocation(3, 24, strand=1), type="CDS", qualifiers={
                "locus_tag": ["gene1"], "translation": ["MACDEFG"], "transl_table": ["11"]}),
            SeqFeature(SimpleLocation(9, 18, strand=1), type="aSDomain", qualifiers={
                "locus_tag": ["gene1"], "aSDomain": ["domain"], "translation": ["CDE"],
                "protein_start": ["2"], "protein_end": ["5"]}),
            SeqFeature(SimpleLocation(27, 27 + 2043 * 3, strand=1), type="CDS", qualifiers={
                "locus_tag": ["gene2"], "translation": ["A" * 2043], "transl_table": ["11"]}),
        ]
        return record

    @staticmethod
    def write_record(path, record):
        stream = io.StringIO()
        SeqIO.write(record, stream, "genbank")
        path.write_text(stream.getvalue(), encoding="utf-8")

    @staticmethod
    def replace_tsv(path, rows):
        path.unlink()
        pilot.write_tsv(path, rows, tuple(rows[0]))

    @staticmethod
    def snapshot(path):
        return {p.relative_to(path).as_posix(): pilot.sha256(p) for p in path.rglob("*") if p.is_file()}

    def update_output_hash(self, relative):
        metadata = pilot.read_json(self.output / "metadata.json")
        metadata["artifact_sha256"][relative] = pilot.sha256(self.output / relative)
        pilot.write_json(self.output / "metadata.json", metadata)

    def test_processed_bytes_stable_ids_unicode_and_actual_parser_checks(self):
        before = self.snapshot(self.source), self.snapshot(self.raw)
        pilot.prepare_repaired(self.args)
        _, metadata, rows = pilot.load_prepared(self.output)
        self.assertEqual(metadata["format"], pilot.REPAIRED_SCHEMA)
        self.assertEqual(metadata["counts"], dict(self.datasets.values()))
        self.assertEqual({r["sample_id"] for r in rows}, {r["sample_id"] for r in self.samples})
        self.assertEqual({r["gcf_raw"] for r in rows}, {"family\u00a0", "阴性单例"})
        changed = next(r for r in rows if r["repair_count"] == "1")
        self.assertNotEqual(changed["source_sha256"], changed["raw_source_sha256"])
        for row in rows:
            copied = self.output / row["dataset"] / "inputs" / (row["sample_id"] + ".gbk")
            self.assertEqual(copied.read_bytes(), (self.source / row["processed_source_relpath"]).read_bytes())
            self.assertEqual(pilot.sha256(copied), row["source_sha256"])
        for check in pilot.read_tsv(self.output / "input_checks.tsv"):
            self.assertEqual(check["cds_count"], "2")
            self.assertEqual(check["domain_token_count"], "1")
            self.assertEqual(check["whole_protein_token_count"], "1")
            self.assertEqual(check["protein_token_count"], "2")
            self.assertEqual(check["max_token_residues"], "2043")
            self.assertEqual(check["esm_1021_window_count"], "4")
            self.assertEqual(check["n_base_count"], "1")
            self.assertEqual(check["ambiguous_base_count"], "3")
        self.assertEqual(before, (self.snapshot(self.source), self.snapshot(self.raw)))
        with self.assertRaisesRegex(ValueError, "new output"):
            pilot.prepare_repaired(self.args)

    def test_relocation_needs_no_original_repair_raw_or_encoder_paths(self):
        pilot.prepare_repaired(self.args)
        self.source.rename(self.root / "source-unavailable")
        self.raw.rename(self.root / "raw-unavailable")
        relocated = self.root / "server" / "inputs"
        relocated.parent.mkdir()
        self.output.rename(relocated)
        with mock.patch.object(build_inventory, "load_parsers", side_effect=AssertionError("no parser imports")):
            _, metadata, rows = pilot.load_prepared(relocated)
        self.assertEqual(len(rows), 2)
        self.assertEqual(metadata["source_root"], str(self.source))
        self.assertFalse(Path(metadata["source_root"]).exists())
        self.assertFalse(Path(metadata["raw_root"]).exists())

    def test_changed_or_missing_source_and_all_artifact_hashes_rejected(self):
        paths = [self.source / self.samples[0]["processed_relpath"],
                 self.source / self.samples[1]["processed_relpath"],
                 self.source / "domain_checks.tsv"]
        for path in paths:
            with self.subTest(path=path):
                old = path.read_bytes()
                path.write_bytes(old + b"changed\n")
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    pilot.prepare_repaired(self.args)
                self.assertFalse(self.output.exists())
                path.unlink()
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    pilot.prepare_repaired(self.args)
                path.write_bytes(old)

    def test_labels_cannot_disagree_even_if_manifest_hash_is_refreshed(self):
        dataset = self.samples[0]["dataset"]
        relative = dataset + "/curated_GCF_assignments.tsv"
        path = self.source / relative
        labels = pilot.read_tsv(path)
        labels[0]["GCF"] = "family"  # Removing the nonbreaking space is a label change.
        self.replace_tsv(path, labels)
        self.metadata["artifact_sha256"][relative] = pilot.sha256(path)
        self.metadata["raw_label_sha256"][relative] = pilot.sha256(path)
        pilot.write_json(self.source / "metadata.json", self.metadata)
        with self.assertRaisesRegex(ValueError, "labels disagree"):
            pilot.prepare_repaired(self.args)
        self.assertFalse(self.output.exists())

    def test_raw_hash_cannot_masquerade_as_processed_input(self):
        pilot.prepare_repaired(self.args)
        path = self.output / "cohort.tsv"
        rows = pilot.read_tsv(path)
        rows[0]["source_sha256"] = rows[0]["raw_source_sha256"]
        self.replace_tsv(path, rows)
        self.update_output_hash("cohort.tsv")
        with self.assertRaisesRegex(ValueError, "embedded repair provenance"):
            pilot.load_prepared(self.output)

    def test_missing_or_tampered_embedded_provenance_and_inputs_rejected(self):
        pilot.prepare_repaired(self.args)
        paths = [self.output / "repair_provenance" / name for name in repaired.PROVENANCE_FILES]
        paths.append(next(self.output.rglob("*.gbk")))
        for path in paths:
            with self.subTest(path=path):
                old = path.read_bytes()
                path.write_bytes(old + b"changed\n")
                with self.assertRaises((ValueError, OSError)):
                    pilot.load_prepared(self.output)
                path.unlink()
                with self.assertRaises((ValueError, OSError)):
                    pilot.load_prepared(self.output)
                path.write_bytes(old)
        # Updating the outer hash cannot untie samples from the original repair metadata.
        path = self.output / "repair_provenance" / "samples.tsv"
        rows = pilot.read_tsv(path)
        rows[0]["gcf_raw"] = "different"
        self.replace_tsv(path, rows)
        self.update_output_hash("repair_provenance/samples.tsv")
        with self.assertRaisesRegex(ValueError, "original repair hashes"):
            pilot.load_prepared(self.output)

    def test_prepared_label_tampering_and_missing_required_hash_rejected(self):
        pilot.prepare_repaired(self.args)
        relative = self.samples[0]["dataset"] + "/curated_GCF_assignments.tsv"
        path = self.output / relative
        rows = pilot.read_tsv(path)
        rows[0]["GCF"] = "different"
        self.replace_tsv(path, rows)
        self.update_output_hash(relative)
        with self.assertRaisesRegex(ValueError, "Prepared labels disagree"):
            pilot.load_prepared(self.output)
        metadata = pilot.read_json(self.output / "metadata.json")
        del metadata["artifact_sha256"]["repair_provenance/repairs.tsv"]
        pilot.write_json(self.output / "metadata.json", metadata)
        with self.assertRaisesRegex(ValueError, "required repair provenance"):
            pilot.load_prepared(self.output)

    def test_output_protection_and_incomplete_or_unsafe_source_metadata(self):
        for output in (self.raw / "new", self.source / "new", self.root):
            with self.subTest(output=output):
                self.args.output = output
                with self.assertRaisesRegex(ValueError, "new output"):
                    pilot.prepare_repaired(self.args)
        self.args.output = self.output
        self.metadata["status"] = "running"
        pilot.write_json(self.source / "metadata.json", self.metadata)
        with self.assertRaisesRegex(ValueError, "incomplete"):
            pilot.prepare_repaired(self.args)
        self.metadata["status"] = "complete"
        self.metadata["artifact_sha256"]["../raw/escape"] = "0" * 64
        pilot.write_json(self.source / "metadata.json", self.metadata)
        with self.assertRaisesRegex(ValueError, "Unsafe relative path"):
            pilot.prepare_repaired(self.args)

    def test_parser_drift_failure_never_publishes_or_changes_sources(self):
        before = self.snapshot(self.source), self.snapshot(self.raw)
        bad = {**self.parsers, "protein": mock.Mock(side_effect=ValueError("injected parser failure"))}
        with mock.patch.object(build_inventory, "load_parsers", return_value=(bad, {})):
            with self.assertRaisesRegex(ValueError, "injected parser failure"):
                pilot.prepare_repaired(self.args)
        self.assertFalse(self.output.exists())
        self.assertEqual(list(self.root.glob(".pilot-gh-*")), [])
        self.assertEqual(before, (self.snapshot(self.source), self.snapshot(self.raw)))
        actual_parser = self.parsers["protein"]

        def changed(path):
            result = actual_parser(path)
            result["tokens"][0]["sequence"] = "AAA"
            return result

        with mock.patch.object(build_inventory, "load_parsers", return_value=(
                {**self.parsers, "protein": changed}, {})):
            with self.assertRaisesRegex(ValueError, "token sequences changed"):
                pilot.prepare_repaired(self.args)
        self.assertFalse(self.output.exists())

    def test_old_schema_cannot_masquerade_as_gh(self):
        pilot.prepare_repaired(self.args)
        metadata = pilot.read_json(self.output / "metadata.json")
        metadata["format"] = pilot.PREPARED_SCHEMA
        pilot.write_json(self.output / "metadata.json", metadata)
        with self.assertRaisesRegex(ValueError, "source policy"):
            pilot.load_prepared(self.output)

    def test_gh_extract_and_validate_use_repaired_run_schema_and_all_routes(self):
        pilot.prepare_repaired(self.args)
        esm, mac, dna = [self.root / name for name in ("esm.pt", "mac.pt", "dna")]
        esm.write_bytes(b"esm fixture")
        mac.write_bytes(b"mac fixture")
        dna.mkdir()
        (dna / "config.json").write_text("{}")
        args = pilot.build_parser().parse_args([
            "extract", "--prepared", str(self.output), "--output", str(self.root / "features"),
            "--encoder-project", str(self.encoder), "--esm-checkpoint", str(esm),
            "--mac-checkpoint", str(mac), "--dna-checkpoint", str(dna), "--trusted-input",
        ])
        commands = []

        def child(command, **kwargs):
            commands.append(command)
            return argparse.Namespace(returncode=0)

        def check(output, dataset, route, rows, saved=None):
            result = {"route": route, "count": len(rows), "dimension": 2}
            self.assertEqual(len(rows), 1)
            if saved is not None:
                self.assertEqual(saved, result)
            return result

        with mock.patch.object(pilot.subprocess, "run", side_effect=child), \
                mock.patch.object(pilot, "check_step", side_effect=check):
            pilot.extract(args)
            state = pilot.read_json(args.output / "run.json")
            self.assertEqual(state["format"], "bgc_repaired_feature_run_v1")
            self.assertEqual(state["status"], "complete")
            self.assertEqual(len(commands), 6)
            self.assertEqual(len(state["steps"]), 6)
            pilot.validate(argparse.Namespace(prepared=self.output, output=args.output))
            args.resume = True
            pilot.extract(args)
            self.assertEqual(len(commands), 6)
            state["format"] = "bgc_be_feature_run_v1"
            pilot.write_json(args.output / "run.json", state)
            with self.assertRaisesRegex(ValueError, "different cohort"):
                pilot.validate(argparse.Namespace(prepared=self.output, output=args.output))
        for command in commands:
            if "experiment.extract_dna_features" in command:
                self.assertEqual(command[command.index("--sequence-mode") + 1], "full")
                self.assertNotIn("--max-bp", command)
                self.assertNotIn("--max-tokens", command)
            else:
                self.assertIn("--trusted-input", command)


if __name__ == "__main__":
    unittest.main()
