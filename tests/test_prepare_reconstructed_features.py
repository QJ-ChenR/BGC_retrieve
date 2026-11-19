"""Server preparation, portable provenance and unresolved-annotation retention.

Fixtures are restored with the real CPU restoration pipeline in temporary
directories. No repository data, checkpoints or model inference are required.
"""

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
import build_inventory
import pilot_features as pilot
import prepare_repaired_features as repaired
import reconstruct_h


class ReconstructedFeatureInputsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="reconstructed-input-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.raw = self.root / "data/raw/curated_datasets"
        self.coordinate = self.root / "data/processed/repaired"
        self.recovered = self.root / "data/processed/reconstructed"
        self.inventory = self.root / "data/processed/inventory"
        self.output = self.root / "data/processed/prepared"
        self.encoder = ROOT.parent / "BGC_annotation"
        self.datasets = {"G": ("Dataset_G_Planomonospora", 1),
                         "H": ("Dataset_H_Reclassification", 2)}
        patch = mock.patch.dict(pilot.DATASETS, self.datasets, clear=True)
        patch.start()
        self.addCleanup(patch.stop)
        self.parsers, _ = build_inventory.load_parsers(self.encoder)
        self.make_coordinate_snapshot()
        self.inventory.mkdir()
        for name in ("samples.tsv", "source_files.tsv"):
            (self.inventory / name).write_text("fixture\n", encoding="utf-8")
        inventory_result = (
            [{**row, "source_sha256": row["raw_source_sha256"]} for row in self.coordinate_rows],
            build_inventory.snapshot(self.raw, self.root),
            {name: pilot.sha256(self.inventory / name) for name in ("samples.tsv", "source_files.tsv")},
        )
        args = argparse.Namespace(workspace=self.root, source=self.coordinate, output=self.recovered,
                                  raw=self.raw, inventory=self.inventory, encoder_project=self.encoder)
        with mock.patch.object(reconstruct_h, "load_inventory", return_value=inventory_result):
            reconstruct_h.reconstruct(args)
        self.args = argparse.Namespace(reconstructed_root=self.recovered, repaired_root=self.coordinate,
                                       datasets=["G", "H"], output=self.output, encoder_project=self.encoder)

    @staticmethod
    def record(name):
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
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(stream.getvalue(), encoding="utf-8")

    @staticmethod
    def write_table(path, rows, fields=None):
        if path.exists():
            path.unlink()
        if fields is None:
            fields = tuple(dict.fromkeys(key for row in rows for key in row))
        pilot.write_tsv(path, rows, fields)

    @staticmethod
    def snapshot(path):
        return {p.relative_to(path).as_posix(): pilot.sha256(p) for p in path.rglob("*") if p.is_file()}

    def make_coordinate_snapshot(self):
        self.coordinate_rows, repairs = [], []
        for key, original, family, unresolved in (
                ("G", "original_G", "family\u00a0", False),
                ("H", "original_H_ready", "家族 α", False),
                ("H", "original_H_review", "家族 β", True)):
            dataset = self.datasets[key][0]
            rel = f"{dataset}/GCF_assignments/{family}/{original}.gbk"
            sample_id = build_inventory.sample_identifier(rel)
            record = self.record(original)
            if key == "H":
                domain = record.features[1]
                domain.type = "misc_feature"
                domain.qualifiers = {"note": ["Geneious type: aSDomain"],
                                     "locus_tag": ["stale_gene1" if not unresolved else "gene1"],
                                     "standard_name": ["existing_domain"]}
                if unresolved:
                    bad = deepcopy(domain)
                    bad.location = SimpleLocation(10, 18, strand=1)
                    record.features.append(bad)
                    record.annotations["structured_comment"] = {"Editing": {"Note": "Original insertion evidence"}}
            self.write_record(self.coordinate / rel, record)
            raw_record = deepcopy(record)
            if key == "G":
                raw_record.features[1].location = SimpleLocation(9, 15, strand=1)
                raw_record.features[1].qualifiers["protein_end"] = ["4"]
                repairs.append({"sample_id": sample_id, "dataset": dataset, "feature_index": "1"})
            self.write_record(self.raw / rel, raw_record)
            protein = self.parsers["protein"](self.coordinate / rel)
            dna, _ = self.parsers["dna"](self.coordinate / rel)
            self.coordinate_rows.append({
                "dataset": dataset, "sample_id": sample_id, "original_bgc_id": original,
                "archive_relpath": rel, "gcf_raw": family,
                "tsv_row_ids_json": repaired._compact([f"{dataset}:tsv:{len(self.coordinate_rows) + 2}"]),
                "raw_source_relpath": f"data/raw/curated_datasets/{rel}",
                "raw_source_sha256": pilot.sha256(self.raw / rel), "processed_relpath": rel,
                "processed_sha256": pilot.sha256(self.coordinate / rel), "record_count": "1",
                "repaired_domain_count": "1" if key == "G" else "0", "length_bp": str(len(dna)),
                "dna_sha256": repaired._digest(dna), "dna_unchanged": "True",
                "translations_unchanged": "True", "unlisted_annotations_unchanged": "True",
                "protein_parser_status_after": "ok", "dna_parser_status_after": "ok",
                "protein_input_sha256_after": reconstruct_h.token_digest(protein),
                "protein_token_count": str(len(protein["tokens"])),
            })
        for dataset, _ in self.datasets.values():
            labels = [{"BGC": r["original_bgc_id"], "GCF": r["gcf_raw"]}
                      for r in self.coordinate_rows if r["dataset"] == dataset]
            for base in (self.coordinate, self.raw):
                self.write_table(base / dataset / "curated_GCF_assignments.tsv", labels)
        self.write_table(self.coordinate / "samples.tsv", self.coordinate_rows)
        self.write_table(self.coordinate / "repairs.tsv", repairs)
        self.write_table(self.coordinate / "review_items.tsv", [], ("sample_id", "dataset", "reason"))
        self.write_table(self.coordinate / "domain_checks.tsv", [], ("sample_id", "dataset", "status"))
        pilot.write_json(self.coordinate / "metadata.json", {
            "format": repaired.REPAIR_SCHEMA, "status": "complete", "raw_root": str(self.raw),
            "raw_unchanged": True,
            "counts": [{"dataset": name, "sample_count": count} for name, count in self.datasets.values()],
            "raw_label_sha256": {name + "/curated_GCF_assignments.tsv": pilot.sha256(
                self.coordinate / name / "curated_GCF_assignments.tsv") for name, _ in self.datasets.values()},
            "artifact_sha256": self.snapshot(self.coordinate),
        })

    def refresh_hash(self, root, relative):
        metadata = pilot.read_json(root / "metadata.json")
        metadata["artifact_sha256"][relative] = pilot.sha256(root / relative)
        pilot.write_json(root / "metadata.json", metadata)

    def prepare(self):
        return pilot.prepare_reconstructed(self.args)

    def reconstruct_after_copying_coordinate_snapshot_to_server(self, historical_raw_root=None):
        """Run real reconstruction on current raw data with historical source paths."""
        historical_raw_root = str(historical_raw_root or "/old/client/BGC_retrieve/data/raw/curated_datasets")
        coordinate_metadata = pilot.read_json(self.coordinate / "metadata.json")
        coordinate_metadata["raw_root"] = historical_raw_root
        pilot.write_json(self.coordinate / "metadata.json", coordinate_metadata)
        self.recovered = self.root / "data/processed/reconstructed_on_server"
        inventory_result = (
            [{**row, "source_sha256": row["raw_source_sha256"]} for row in self.coordinate_rows],
            build_inventory.snapshot(self.raw, self.root),
            {name: pilot.sha256(self.inventory / name) for name in ("samples.tsv", "source_files.tsv")},
        )
        args = argparse.Namespace(workspace=self.root, source=self.coordinate, output=self.recovered,
                                  raw=self.raw, inventory=self.inventory, encoder_project=self.encoder)
        with mock.patch.object(reconstruct_h, "load_inventory", return_value=inventory_result):
            reconstruct_h.reconstruct(args)
        self.args.reconstructed_root = self.recovered
        current = pilot.read_json(self.recovered / "metadata.json")
        previous = pilot.read_json(self.recovered / "source_provenance/metadata.json")
        self.assertEqual(current["raw_root"], str(self.raw))
        self.assertEqual(previous["raw_root"], historical_raw_root)
        self.assertNotEqual(current["raw_root"], previous["raw_root"])
        return historical_raw_root

    def test_reconstruction_rerun_on_server_accepts_historical_raw_root_and_remains_portable(self):
        historical_raw_root = self.reconstruct_after_copying_coordinate_snapshot_to_server()
        expected = {row["sample_id"]: row for row in self.coordinate_rows}
        expected.update({row["sample_id"]: row for row in pilot.read_tsv(self.recovered / "samples.tsv")})
        queues = []
        for selection in (["G", "H"], ["H"]):
            with self.subTest(datasets=selection):
                self.args.datasets = selection
                self.args.repaired_root = self.coordinate if "G" in selection else None
                output = self.root / "data/processed" / ("prepared_" + "".join(selection))
                self.args.output = output
                self.prepare()
                _, metadata, rows = pilot.load_prepared(output)
                datasets = {self.datasets[key][0] for key in selection}
                expected_ids = {sid for sid, row in expected.items() if row["dataset"] in datasets}
                self.assertEqual({row["sample_id"] for row in rows}, expected_ids)
                self.assertEqual(metadata["status"], "complete")
                self.assertEqual(metadata["annotation_status"], "needs_review")
                self.assertEqual(metadata["unresolved_asdomain_count"], 1)
                self.assertFalse(metadata["annotation_complete"])
                roots = {self.datasets["H"][0]: str(self.raw)}
                if "G" in selection:
                    roots[self.datasets["G"][0]] = historical_raw_root
                self.assertEqual(metadata["raw_roots"], roots)
                for row in rows:
                    original = expected[row["sample_id"]]
                    for key in ("original_bgc_id", "gcf_raw", "raw_source_sha256"):
                        self.assertEqual(row[key], original[key])
                    self.assertEqual(row["source_sha256"], original["processed_sha256"])
                queues.append((output, rows, metadata["unresolved_sample_ids"]))
        for source in (self.coordinate, self.recovered, self.raw):
            source.rename(self.root / (source.name + "-unavailable"))
        for output, original_rows, review_ids in queues:
            relocated = self.root / "deployed" / output.name
            relocated.parent.mkdir(exist_ok=True)
            output.rename(relocated)
            with mock.patch.object(build_inventory, "load_parsers", side_effect=AssertionError("no parser imports")):
                _, metadata, rows = pilot.load_prepared(relocated)
            self.assertEqual(rows, original_rows)
            self.assertEqual(metadata["unresolved_sample_ids"], review_ids)
            self.assertEqual(metadata["unresolved_asdomain_count"], 1)

    def test_server_rerun_rejects_conflicting_ancestral_raw_hash_with_refreshed_manifests(self):
        self.reconstruct_after_copying_coordinate_snapshot_to_server()
        self.args.datasets = ["H"]
        self.args.repaired_root = None
        ancestor = self.recovered / "source_provenance"
        rows = pilot.read_tsv(ancestor / "samples.tsv")
        h_sample = next(row for row in rows if row["dataset"] == self.datasets["H"][0])
        h_sample["raw_source_sha256"] = "0" * 64
        self.write_table(ancestor / "samples.tsv", rows)
        self.refresh_hash(ancestor, "samples.tsv")
        metadata = pilot.read_json(self.recovered / "metadata.json")
        metadata["source_metadata_sha256"] = pilot.sha256(ancestor / "metadata.json")
        for relative in ("source_provenance/samples.tsv", "source_provenance/metadata.json"):
            metadata["artifact_sha256"][relative] = pilot.sha256(self.recovered / relative)
        pilot.write_json(self.recovered / "metadata.json", metadata)
        with self.assertRaisesRegex(ValueError, "identity/sequence provenance differs: raw_source_sha256"):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_server_rerun_still_checks_ancestral_metadata_hash(self):
        self.reconstruct_after_copying_coordinate_snapshot_to_server()
        metadata = pilot.read_json(self.recovered / "metadata.json")
        metadata["source_metadata_sha256"] = "0" * 64
        pilot.write_json(self.recovered / "metadata.json", metadata)
        with self.assertRaisesRegex(ValueError, "Coordinate-repair metadata differs"):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_server_rerun_rejects_changed_gbk_even_when_outer_artifact_hash_is_refreshed(self):
        self.reconstruct_after_copying_coordinate_snapshot_to_server()
        sample = pilot.read_tsv(self.recovered / "samples.tsv")[0]
        path = self.recovered / sample["processed_relpath"]
        path.write_bytes(path.read_bytes() + b"\n")
        self.refresh_hash(self.recovered, sample["processed_relpath"])
        with self.assertRaisesRegex(ValueError, "H processed input hash/path differs"):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_h_only_server_rerun_also_protects_historical_ancestral_raw_directory(self):
        historical = self.root / "old-project/data/raw/curated_datasets"
        historical.mkdir(parents=True)
        (historical / "backup.txt").write_bytes(b"historical source must remain unchanged\n")
        self.reconstruct_after_copying_coordinate_snapshot_to_server(historical)
        self.args.datasets = ["H"]
        self.args.repaired_root = None
        self.args.output = historical / "new_queue"
        protected_sources = (self.coordinate, self.recovered, self.raw, historical)
        before = [self.snapshot(path) for path in protected_sources]
        with self.assertRaisesRegex(ValueError, "new output"):
            self.prepare()
        self.assertFalse(self.args.output.exists())
        self.assertEqual(before, [self.snapshot(path) for path in protected_sources])

    def test_h_and_g_use_distinct_sources_and_preserve_reviewed_sample(self):
        before = [self.snapshot(root) for root in (self.raw, self.coordinate, self.recovered)]
        self.prepare()
        _, metadata, rows = pilot.load_prepared(self.output)
        self.assertEqual(metadata["format"], pilot.RECONSTRUCTED_SCHEMA)
        self.assertEqual(metadata["status"], "complete")
        self.assertEqual(metadata["counts"], dict(self.datasets.values()))
        self.assertEqual(metadata["annotation_status"], "needs_review")
        self.assertFalse(metadata["annotation_complete"])
        self.assertEqual(metadata["unresolved_asdomain_count"], 1)
        self.assertEqual(metadata["review_policy"], "retain")
        self.assertEqual(len(rows), 3)
        review = next(row for row in rows if row["original_bgc_id"] == "original_H_review")
        self.assertEqual(metadata["unresolved_sample_ids"], [review["sample_id"]])
        self.assertEqual(review["annotation_status"], "needs_review")
        self.assertEqual(review["unresolved_asdomain_count"], "1")
        self.assertEqual({r["gcf_raw"] for r in rows}, {"family\u00a0", "家族 α", "家族 β"})
        self.assertEqual(len({r["source_kind"] for r in rows}), 2)
        for row in rows:
            source = self.coordinate if row["dataset"] == self.datasets["G"][0] else self.recovered
            copied = self.output / row["dataset"] / "inputs" / (row["sample_id"] + ".gbk")
            self.assertEqual(copied.read_bytes(), (source / row["processed_source_relpath"]).read_bytes())
            self.assertEqual(pilot.sha256(copied), row["source_sha256"])
            self.assertNotEqual(row["source_sha256"], row["raw_source_sha256"])
        checks = pilot.read_tsv(self.output / "input_checks.tsv")
        self.assertEqual(len(checks), 3)
        self.assertTrue(all(row["protein_token_count"] == "2" for row in checks))
        self.assertEqual(before, [self.snapshot(root) for root in (self.raw, self.coordinate, self.recovered)])
        with self.assertRaises(ValueError):
            self.prepare()

    def test_h_only_does_not_need_coordinate_source_directory(self):
        self.args.datasets = ["H"]
        self.args.repaired_root = None
        self.coordinate.rename(self.root / "coordinate-source-unavailable")
        self.raw.rename(self.root / "raw-unavailable")
        self.prepare()
        _, metadata, rows = pilot.load_prepared(self.output)
        self.assertEqual(metadata["counts"], {self.datasets["H"][0]: 2})
        self.assertEqual(len(rows), 2)
        self.assertTrue(any(r["annotation_status"] == "needs_review" for r in rows))

    def test_relocated_queue_loads_without_original_paths_or_parser_imports(self):
        self.prepare()
        for source in (self.coordinate, self.recovered, self.raw):
            source.rename(self.root / (source.name + "-unavailable"))
        relocated = self.root / "server/prepared"
        relocated.parent.mkdir()
        self.output.rename(relocated)
        with mock.patch.object(build_inventory, "load_parsers", side_effect=AssertionError("CPU parsers must not load")):
            _, metadata, rows = pilot.load_prepared(relocated)
        self.assertEqual(len(rows), 3)
        self.assertEqual(metadata["annotation_status"], "needs_review")

    def test_equivalent_g_snapshot_from_another_server_preserves_both_raw_roots(self):
        coordinate_metadata = pilot.read_json(self.coordinate / "metadata.json")
        other_raw_root = "/srv/other-machine/BGC_retrieve/data/raw/curated_datasets"
        coordinate_metadata["raw_root"] = other_raw_root
        pilot.write_json(self.coordinate / "metadata.json", coordinate_metadata)
        self.prepare()
        for source in (self.coordinate, self.recovered, self.raw):
            source.rename(self.root / (source.name + "-unavailable"))
        with mock.patch.object(build_inventory, "load_parsers", side_effect=AssertionError("no parser imports")):
            _, metadata, rows = pilot.load_prepared(self.output)
        self.assertEqual(len(rows), 3)
        self.assertEqual(metadata["raw_roots"], {
            self.datasets["H"][0]: str(self.raw), self.datasets["G"][0]: other_raw_root,
        })
        self.assertEqual(metadata["raw_root"], str(self.raw))
        self.assertFalse(self.raw.exists())

    def test_g_content_must_match_h_ancestry_even_with_refreshed_g_hashes(self):
        samples_path = self.coordinate / "samples.tsv"
        original = samples_path.read_bytes()
        for field in ("protein_input_sha256_after", "dna_sha256"):
            with self.subTest(field=field):
                rows = pilot.read_tsv(samples_path)
                row = next(r for r in rows if r["dataset"] == self.datasets["G"][0])
                row[field] = "0" * 64
                self.write_table(samples_path, rows)
                self.refresh_hash(self.coordinate, "samples.tsv")
                # The G snapshot's own files/hash manifest are consistent; it is
                # specifically the cross-snapshot ancestry comparison that fails.
                with self.assertRaisesRegex(ValueError, "G content differs.*ancestry"):
                    self.prepare()
                self.assertFalse(self.output.exists())
                samples_path.write_bytes(original)
                self.refresh_hash(self.coordinate, "samples.tsv")

    def test_current_workspace_raw_is_protected_even_with_unavailable_old_source(self):
        self.args.datasets = ["H"]
        self.args.repaired_root = None
        self.coordinate.rename(self.root / "coordinate-source-unavailable")
        self.raw.rename(self.root / "raw-unavailable")
        self.args.output = self.root / "data/raw/new_should_reject"
        with mock.patch.object(pilot, "ROOT", self.root):
            with self.assertRaisesRegex(ValueError, "new output"):
                self.prepare()
        self.assertFalse(self.args.output.exists())

    def test_unlisted_or_changed_source_artifacts_block_preparation(self):
        sample = pilot.read_tsv(self.recovered / "samples.tsv")[0]
        paths = [self.recovered / sample["processed_relpath"], self.recovered / "reconstructions.tsv",
                 self.coordinate / "domain_checks.tsv"]
        for path in paths:
            with self.subTest(path=path):
                old = path.read_bytes()
                path.write_bytes(old + b"tampered\n")
                with self.assertRaises(ValueError):
                    self.prepare()
                self.assertFalse(self.output.exists())
                path.write_bytes(old)
        extra = self.recovered / self.datasets["H"][0] / "GCF_assignments/extra.gbk"
        extra.write_bytes((self.recovered / sample["processed_relpath"]).read_bytes())
        with self.assertRaises(ValueError):
            self.prepare()

    def test_labels_cannot_change_even_after_outer_hash_update(self):
        relative = self.datasets["H"][0] + "/curated_GCF_assignments.tsv"
        path = self.recovered / relative
        rows = pilot.read_tsv(path)
        rows[0]["GCF"] = "different family"
        self.write_table(path, rows)
        self.refresh_hash(self.recovered, relative)
        with self.assertRaises(ValueError):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_missing_review_rows_and_incorrect_reconstruction_counts_rejected(self):
        review_path = self.recovered / "review_items.tsv"
        original = review_path.read_bytes()
        self.write_table(review_path, [], ("sample_id", "dataset", "feature_index", "status", "reason"))
        self.refresh_hash(self.recovered, "review_items.tsv")
        with self.assertRaises(ValueError):
            self.prepare()
        review_path.write_bytes(original)
        self.refresh_hash(self.recovered, "review_items.tsv")
        rows = pilot.read_tsv(self.recovered / "samples.tsv")
        rows[0]["restored_domain_count"] = "200"
        self.write_table(self.recovered / "samples.tsv", rows)
        self.refresh_hash(self.recovered, "samples.tsv")
        with self.assertRaises(ValueError):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_prepared_metadata_cannot_hide_unresolved_annotations(self):
        self.prepare()
        metadata_path = self.output / "metadata.json"
        original = metadata_path.read_bytes()
        for field, false_value in (("annotation_status", "complete"), ("annotation_complete", True),
                                   ("unresolved_asdomain_count", 0), ("unresolved_sample_ids", []),
                                   ("review_policy", "drop")):
            with self.subTest(field=field):
                metadata = pilot.read_json(metadata_path)
                metadata[field] = false_value
                pilot.write_json(metadata_path, metadata)
                with self.assertRaises(ValueError):
                    pilot.load_prepared(self.output)
                metadata_path.write_bytes(original)

    def test_cohort_hash_or_annotation_status_cannot_override_embedded_provenance(self):
        self.prepare()
        path = self.output / "cohort.tsv"
        original = path.read_bytes()
        for field, false_value in (("source_sha256", "0" * 64), ("annotation_status", "complete")):
            with self.subTest(field=field):
                rows = pilot.read_tsv(path)
                row = next(r for r in rows if r["original_bgc_id"] == "original_H_review")
                row[field] = false_value
                self.write_table(path, rows)
                self.refresh_hash(self.output, "cohort.tsv")
                with self.assertRaises(ValueError):
                    pilot.load_prepared(self.output)
                path.write_bytes(original)
                self.refresh_hash(self.output, "cohort.tsv")

    def test_nested_provenance_tampering_is_detected_after_outer_hash_update(self):
        self.prepare()
        path = self.output / "reconstruction_provenance/source_provenance/samples.tsv"
        self.assertTrue(path.is_file())
        rows = pilot.read_tsv(path)
        rows[0]["gcf_raw"] = "tampered"
        self.write_table(path, rows)
        self.refresh_hash(self.output, path.relative_to(self.output).as_posix())
        with self.assertRaises(ValueError):
            pilot.load_prepared(self.output)

    def test_parser_failure_does_not_publish_or_mutate_sources(self):
        before = [self.snapshot(root) for root in (self.raw, self.coordinate, self.recovered)]
        parsers, provenance = build_inventory.load_parsers(self.encoder)
        bad = {**parsers, "protein": mock.Mock(side_effect=ValueError("injected parser failure"))}
        with mock.patch.object(build_inventory, "load_parsers", return_value=(bad, provenance)):
            with self.assertRaisesRegex(ValueError, "injected parser failure"):
                self.prepare()
        self.assertFalse(self.output.exists())
        self.assertEqual(before, [self.snapshot(root) for root in (self.raw, self.coordinate, self.recovered)])


if __name__ == "__main__":
    unittest.main()
