"""Server runner integration with recorded annotation gaps; no model inference."""
import argparse
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import pilot_features as pilot


class ReconstructedRunTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="reconstructed-run-test-")
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.prepared = self.root / "prepared"
        self.prepared.mkdir()
        self.output = self.root / "features"
        self.project = self.root / "encoder"
        for relative in pilot.ENCODER_FILES:
            path = self.project / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# CPU test fixture\n")
        self.datasets = [pilot.DATASETS[key][0] for key in ("G", "H")]
        self.rows = [{"sample_id": "g", "dataset": self.datasets[0], "annotation_status": "complete"},
                     {"sample_id": "h", "dataset": self.datasets[1], "annotation_status": "needs_review"}]
        pilot.write_tsv(self.prepared / "cohort.tsv", self.rows, tuple(self.rows[0]))
        self.metadata = {
            "format": pilot.RECONSTRUCTED_SCHEMA, "status": "complete", "datasets": self.datasets,
            "raw_root": str(self.root / "raw"), "source_root": str(self.root / "restored"),
            "raw_roots": {"G": str(self.root / "server_raw"), "H": str(self.root / "raw")},
            "source_roots": {"G": str(self.root / "coordinate"), "H": str(self.root / "restored")},
            "artifact_sha256": {"cohort.tsv": pilot.sha256(self.prepared / "cohort.tsv")},
            "annotation_status": "needs_review", "annotation_complete": False,
            "unresolved_asdomain_count": 2, "unresolved_sample_ids": ["h"], "review_policy": "retain",
        }
        pilot.write_json(self.prepared / "metadata.json", self.metadata)
        checkpoint = self.root / "fixture.pt"
        checkpoint.write_bytes(b"No real weights are loaded")
        self.args = pilot.build_parser().parse_args([
            "extract", "--prepared", str(self.prepared), "--output", str(self.output),
            "--encoder-project", str(self.project), "--esm-checkpoint", str(checkpoint),
            "--mac-checkpoint", str(checkpoint), "--dna-checkpoint", str(checkpoint),
        ])
        load = mock.patch.object(pilot, "load_prepared", side_effect=lambda path: (
            self.prepared, deepcopy(self.metadata), deepcopy(self.rows)))
        load.start()
        self.addCleanup(load.stop)
        self.commands = []

    def child(self, command, **kwargs):
        self.commands.append(command)
        self.assertEqual(kwargs["cwd"], str(self.project))
        return argparse.Namespace(returncode=0)

    def check(self, output, dataset, route, rows, saved=None):
        # A needs_review input stays in every route's coverage check.
        self.assertEqual(len(rows), 1)
        if dataset == self.datasets[1]:
            self.assertEqual(rows[0]["sample_id"], "h")
            self.assertEqual(rows[0]["annotation_status"], "needs_review")
        result = {"route": route, "count": len(rows), "dimension": 2}
        if saved is not None:
            self.assertEqual(saved, result)
        return result

    def launch(self):
        with mock.patch.object(pilot.subprocess, "run", side_effect=self.child), \
                mock.patch.object(pilot, "check_step", side_effect=self.check):
            pilot.extract(self.args)

    def test_all_routes_keep_review_context_and_share_the_new_cache(self):
        self.launch()
        state = pilot.read_json(self.output / "run.json")
        self.assertEqual(state["format"], "bgc_reconstructed_feature_run_v1")
        self.assertEqual(state["status"], "complete")
        self.assertFalse(state["annotation_complete"])
        for key in pilot.ANNOTATION_FIELDS:
            self.assertEqual(state[key], self.metadata[key])
            self.assertEqual(state["config"][key], self.metadata[key])
        self.assertEqual(len(self.commands), 6)
        self.assertEqual(len(state["steps"]), 6)
        for dataset in self.datasets:
            esm = state["steps"][dataset + "/esm2"]["command"]
            mac = state["steps"][dataset + "/mac"]["command"]
            dna = state["steps"][dataset + "/dna_full"]["command"]
            cache = str(self.output / dataset / "esm2/tokens.pth")
            self.assertEqual(esm[esm.index("--save-cache") + 1], cache)
            self.assertEqual(mac[mac.index("--esm-cache") + 1], cache)
            self.assertEqual(dna[dna.index("--sequence-mode") + 1], "full")
            self.assertNotIn("--skip-invalid", esm)
            self.assertNotIn("--skip-invalid", dna)
        with mock.patch.object(pilot, "check_step", side_effect=self.check):
            pilot.validate(argparse.Namespace(prepared=self.prepared, output=self.output))
        self.args.resume = True
        self.launch()
        self.assertEqual(len(self.commands), 6)

    def test_run_cannot_drop_or_relabel_annotation_review_during_resume_or_validation(self):
        self.launch()
        original = pilot.read_json(self.output / "run.json")
        for field in pilot.ANNOTATION_FIELDS:
            with self.subTest(field=field):
                state = deepcopy(original)
                state.pop(field)
                pilot.write_json(self.output / "run.json", state)
                with self.assertRaisesRegex(ValueError, "annotation review"):
                    pilot.validate(argparse.Namespace(prepared=self.prepared, output=self.output))
                self.args.resume = True
                with self.assertRaisesRegex(ValueError, "annotation review"):
                    self.launch()
        state = deepcopy(original)
        state["config"]["annotation_complete"] = True
        pilot.write_json(self.output / "run.json", state)
        with self.assertRaisesRegex(ValueError, "annotation review"):
            self.launch()

    def test_prepared_metadata_drift_rejects_resume(self):
        self.launch()
        metadata = {**self.metadata, "unexpected_change": True}
        pilot.write_json(self.prepared / "metadata.json", metadata)
        self.args.resume = True
        with self.assertRaisesRegex(ValueError, "annotation review/provenance"):
            self.launch()

    def test_feature_output_cannot_write_inside_either_source_snapshot(self):
        for source in (*self.metadata["source_roots"].values(), *self.metadata["raw_roots"].values()):
            with self.subTest(source=source):
                self.args.output = Path(source) / "features"
                with self.assertRaisesRegex(ValueError, "separate"):
                    self.launch()
                self.assertFalse(self.args.output.exists())
        with mock.patch.object(pilot, "ROOT", self.root / "server_project"):
            self.args.output = pilot.ROOT / "data/raw/feature_output"
            with self.assertRaisesRegex(ValueError, "separate"):
                self.launch()
            self.assertFalse(self.args.output.exists())

    def test_cli_exposes_recovered_h_and_optional_g_without_changing_old_defaults(self):
        args = pilot.build_parser().parse_args(["prepare-reconstructed", "--datasets", "H"])
        self.assertEqual(args.datasets, ["H"])
        self.assertEqual(args.handler, pilot.prepare_reconstructed)
        args = pilot.build_parser().parse_args(["prepare-reconstructed"])
        self.assertEqual(args.datasets, ["G", "H"])
        self.assertEqual(args.output.name, "pilot_gh_reconstructed_v1")
        self.assertEqual(pilot.feature_run_format({"format": pilot.PREPARED_SCHEMA}), "bgc_be_feature_run_v1")
        self.assertEqual(pilot.feature_run_format({"format": pilot.REPAIRED_SCHEMA}), "bgc_repaired_feature_run_v1")


if __name__ == "__main__":
    unittest.main()
