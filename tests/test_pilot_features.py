"""Input provenance and subprocess orchestration tests; no model/GPU loading."""

import argparse
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import pilot_features as pilot


class PilotWorkflowTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="be-pilot-test-")
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.raw = self.root / "raw"
        self.inventory = self.root / "inventory"
        self.prepared = self.root / "prepared"
        self.output = self.root / "features"
        self.inventory.mkdir()
        self.dataset = "Dataset_B_Divergent"
        patch = mock.patch.dict(pilot.DATASETS, {"B": (self.dataset, 2)}, clear=True)
        patch.start()
        self.addCleanup(patch.stop)
        self.samples, sources = [], []
        for index, family in enumerate(("family\u00a0", "neg_singleton")):
            original_id = f"original_{index}"
            rel = f"{self.dataset}/GCF_assignments/{family}/{original_id}.gbk"
            path = self.raw / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(f"synthetic GenBank bytes {index}\n".encode())
            sample_id = "bgc_" + pilot.hashlib.sha256(
                (pilot.INVENTORY_SCHEMA + "\0" + rel).encode()).hexdigest()[:24]
            self.samples.append({
                "dataset": self.dataset, "sample_id": sample_id, "original_bgc_id": original_id,
                "archive_relpath": rel, "source_sha256": pilot.sha256(path),
                "tsv_row_ids_json": json.dumps([f"{self.dataset}:tsv:{index + 2}"]),
                "tsv_labels_raw_json": json.dumps([family]), "folder_gcf_raw": family,
                "mapping_status": "exact_path_label", "record_count": "1",
                "dna_parser_status": "ok", "protein_parser_status": "ok",
                "excluded": "False", "repairs_applied_json": "[]",
            })
            sources.append({"archive_relpath": rel, "source_sha256": pilot.sha256(path)})
        labels = self.raw / self.dataset / "curated_GCF_assignments.tsv"
        pilot.write_tsv(labels, [{"BGC": r["original_bgc_id"], "GCF": r["folder_gcf_raw"]}
                                 for r in self.samples], ("BGC", "GCF"))
        sources.append({"archive_relpath": labels.relative_to(self.raw).as_posix(),
                        "source_sha256": pilot.sha256(labels)})
        pilot.write_tsv(self.inventory / "samples.tsv", self.samples, tuple(self.samples[0]))
        pilot.write_tsv(self.inventory / "source_files.tsv", sources, tuple(sources[0]))
        pilot.write_json(self.inventory / "metadata.json", {
            "schema_version": pilot.INVENTORY_SCHEMA,
            "output_sha256": {name: pilot.sha256(self.inventory / name)
                              for name in ("samples.tsv", "source_files.tsv")},
        })
        self.prepare_args = argparse.Namespace(raw_root=self.raw, inventory=self.inventory,
                                               output=self.prepared, datasets=["B"])

    def extraction_args(self, *extra):
        project = self.root / "encoder"
        for name in pilot.ENCODER_FILES:
            file = project / name
            file.parent.mkdir(parents=True, exist_ok=True)
            file.write_text("# fixture\n")
        esm, mac, dna = (self.root / name for name in ("esm.pt", "mac.ckpt", "dna"))
        esm.write_bytes(b"fixture esm")
        mac.write_bytes(b"fixture mac")
        dna.mkdir(exist_ok=True)
        (dna / "config.json").write_text("{}")
        return pilot.build_parser().parse_args([
            "extract", "--prepared", str(self.prepared), "--encoder-project", str(project),
            "--output", str(self.output), "--esm-checkpoint", str(esm),
            "--mac-checkpoint", str(mac), "--dna-checkpoint", str(dna), *extra,
        ])

    def test_prepare_keeps_bytes_stable_ids_and_unicode_labels(self):
        pilot.prepare(self.prepare_args)
        _, metadata, rows = pilot.load_prepared(self.prepared)
        self.assertEqual(metadata["counts"], {self.dataset: 2})
        self.assertIn("family\u00a0", [r["gcf_raw"] for r in rows])
        for row in rows:
            copy = self.prepared / self.dataset / "inputs" / (row["sample_id"] + ".gbk")
            self.assertEqual(copy.read_bytes(), (self.raw / row["archive_relpath"]).read_bytes())
        with self.assertRaisesRegex(ValueError, "new output"):
            pilot.prepare(self.prepare_args)

    def test_changed_source_or_labels_fail_before_publish(self):
        for path in (self.raw / self.samples[0]["archive_relpath"],
                     self.raw / self.dataset / "curated_GCF_assignments.tsv"):
            original = path.read_bytes()
            path.write_bytes(original + b"changed\n")
            with self.assertRaisesRegex(ValueError, "differs|differ"):
                pilot.prepare(self.prepare_args)
            self.assertFalse(self.prepared.exists())
            path.write_bytes(original)

    def test_extra_raw_file_and_changed_prepared_input_rejected(self):
        extra = self.raw / self.dataset / "GCF_assignments" / "extra.gbk"
        extra.write_bytes(b"not in inventory")
        with self.assertRaisesRegex(ValueError, "file set"):
            pilot.prepare(self.prepare_args)
        extra.unlink()
        pilot.prepare(self.prepare_args)
        staged = next((self.prepared / self.dataset / "inputs").iterdir())
        staged.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "Prepared input changed"):
            pilot.load_prepared(self.prepared)

    def test_commands_share_cache_and_full_dna_has_no_limits(self):
        pilot.prepare(self.prepare_args)
        args = self.extraction_args("--routes", "mac", "dna_full")
        prepared, metadata, _ = pilot.load_prepared(self.prepared)
        config = pilot.extraction_config(args, prepared, metadata)
        self.assertEqual(config["routes"], list(pilot.ROUTES))
        commands = {route: pilot.command_for(config, self.dataset, route, self.output)
                    for route in pilot.ROUTES}
        esm, mac, dna = (commands[r] for r in pilot.ROUTES)
        self.assertEqual(esm[esm.index("--save-cache") + 1], mac[mac.index("--esm-cache") + 1])
        self.assertEqual(dna[dna.index("--sequence-mode") + 1], "full")
        self.assertEqual(dna[dna.index("--dtype") + 1], "bfloat16")
        for flag in ("--skip-invalid", "--max-bp", "--max-tokens", "--ids-from"):
            self.assertNotIn(flag, dna)
        self.assertNotIn("--trusted-input", esm)

    def test_resume_reuses_complete_stages_and_archives_failure(self):
        pilot.prepare(self.prepare_args)
        args = self.extraction_args()
        attempts = []

        def child(command, **kwargs):
            route = Path(command[command.index("--output") + 1]).name
            attempts.append(route)
            return argparse.Namespace(returncode=1 if attempts == ["esm2", "mac"] else 0)

        def checked(output, dataset, route, rows, saved=None):
            result = {"route": route, "count": len(rows), "dimension": 2}
            if saved is not None:
                self.assertEqual(saved, result)
            return result

        with mock.patch.object(pilot.subprocess, "run", side_effect=child), \
                mock.patch.object(pilot, "check_step", side_effect=checked):
            with self.assertRaisesRegex(RuntimeError, "mac failed"):
                pilot.extract(args)
            state = pilot.read_json(self.output / "run.json")
            self.assertEqual(state["status"], "failed")
            self.assertEqual(state["steps"][self.dataset + "/esm2"]["status"], "complete")
            self.assertFalse((self.output / ".running.lock").exists())
            args.resume = True
            pilot.extract(args)
            self.assertEqual(attempts, ["esm2", "mac", "mac", "dna_full"])
            state = pilot.read_json(self.output / "run.json")
            self.assertEqual(state["status"], "complete")
            self.assertEqual(len(state["archived_attempts"]), 1)
            self.assertTrue((self.output / state["archived_attempts"][0] / "execution.log").exists())
            pilot.extract(args)
            self.assertEqual(len(attempts), 4)
            args.mac_member = 1
            with self.assertRaisesRegex(ValueError, "configuration/code/checkpoint changed"):
                pilot.extract(args)

    def test_existing_lock_does_not_mutate_run_status(self):
        pilot.prepare(self.prepare_args)
        args = self.extraction_args("--routes", "dna_full")
        prepared, metadata, _ = pilot.load_prepared(self.prepared)
        config = pilot.extraction_config(args, prepared, metadata)
        self.output.mkdir()
        pilot.write_json(self.output / "run.json", {"config": config, "status": "complete"})
        (self.output / ".running.lock").write_text("another process")
        args.resume = True
        with self.assertRaisesRegex(ValueError, "Run lock exists"):
            pilot.extract(args)
        self.assertEqual(pilot.read_json(self.output / "run.json")["status"], "complete")


if __name__ == "__main__":
    unittest.main()
