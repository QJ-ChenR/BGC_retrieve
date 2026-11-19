"""Exercise PCA/average output, diagnostic selection, and full-cohort failures."""

from collections import defaultdict
import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import clustering_common as common
import cluster_pca_average as pipeline


class PCAAveragePipelineTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="pca-average-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.cohort = self.root / "cohort.tsv"
        self.cohort.write_text("synthetic cohort provenance fixture\n", encoding="utf-8")

    def bundle(self, values):
        unit, distances = common.cosine_distances(values)
        count = len(unit)
        rows = [{"sample_id": f"bgc_{i:02d}", "original_bgc_id": f"original_{i}",
                 "gcf_raw": "family\u00a0" if i < count // 2 else "second",
                 "dataset": "Dataset_B_Divergent", "source_sha256": f"{i:064x}"}
                for i in range(count)]
        return {"dataset": "Dataset_B_Divergent", "route": "esm2", "rows": rows,
                "ids": [r["sample_id"] for r in rows], "labels": [r["gcf_raw"] for r in rows],
                "vectors": unit, "distances": distances,
                "source_row_indices": list(reversed(range(count))),
                "provenance": {"dataset": "Dataset_B_Divergent", "route": "esm2",
                               "artifact_sha256": {}, "input_count": count}}

    def args(self, dimensions=(2, 99)):
        args = pipeline.build_parser().parse_args([
            "--cohort", str(self.cohort), "--output", str(self.root / "results"),
            "--datasets", "B", "--routes", "esm2", "--no-plots",
            "--dimensions", *map(str, dimensions), "--variance-targets", "0.95",
        ])
        args.reference_average = None
        return args

    def run_bundle(self, bundle, args):
        with mock.patch.object(pipeline.common, "load_inputs", return_value=[bundle]), contextlib.redirect_stdout(io.StringIO()):
            return pipeline.run(args)

    def test_complete_pipeline_preserves_ids_raw_geometry_and_all_condition_coverage(self):
        bundle = self.bundle(np.random.default_rng(17).normal(size=(9, 4)))
        output = self.run_bundle(bundle, self.args())
        directory = output / bundle["dataset"] / bundle["route"]
        ids = common.read_tsv(directory / "bgc_ids.tsv")
        self.assertEqual([r["bgc_id"] for r in ids], bundle["ids"])
        self.assertEqual([int(r["source_feature_row"]) for r in ids], bundle["source_row_indices"])
        np.testing.assert_array_equal(np.load(directory / "raw/cosine_distances.npy"), bundle["distances"])
        np.testing.assert_array_equal(np.load(directory / "raw/unit_vectors.npy"), bundle["vectors"])

        matched = common.read_tsv(output / "matched_metrics.tsv")
        self.assertIn("centered_only", {r["condition"] for r in matched})
        self.assertIn("pca_0002", {r["condition"] for r in matched})
        full_rank = next(r for r in matched if "full_rank" in r["selectors"].split(","))
        np.testing.assert_allclose(np.load(directory / full_rank["condition"] / "cosine_distances.npy"),
                                   np.load(directory / "centered_only/cosine_distances.npy"), atol=1e-12)
        with np.load(directory / "transform.npz") as transform:
            projected = (bundle["vectors"] - transform["mean"]) @ transform["components"].T
            expected_unit, _ = common.cosine_distances(projected)
            np.testing.assert_allclose(np.load(directory / full_rank["condition"] / "unit_vectors.npy"),
                                       expected_unit, atol=1e-12)

        for selected in matched:
            self.assertEqual(int(selected["sample_count"]), len(bundle["ids"]))
            self.assertEqual(float(selected["coverage"]), 1)
            assignments = common.read_tsv(directory / selected["condition"] / "clusters.tsv")
            by_run = defaultdict(list)
            for row in assignments:
                by_run[row["run_id"]].append(row)
            for members in by_run.values():
                self.assertEqual([r["sample_id"] for r in members], bundle["ids"])
                self.assertEqual([r["gcf_raw"] for r in members], bundle["labels"])
            metrics = common.read_tsv(directory / selected["condition"] / "metrics.tsv")
            expected = min(metrics, key=lambda r: (abs(int(r["n_clusters"]) - int(r["true_gcf_count"])),
                                                  float(r["threshold"])))
            self.assertEqual(selected["run_id"], expected["run_id"])

        metadata = common.read_json(output / "metadata.json")
        self.assertEqual(metadata["status"], "complete")
        self.assertEqual(metadata["raw_reference_checks"][0]["status"], "unavailable")
        self.assertEqual(metadata["parameters"]["leiden_runs"], 0)
        for relative, digest in metadata["artifact_sha256"].items():
            self.assertEqual(common.file_sha256(output / relative), digest)

    def test_oversized_request_logged_without_mislabeled_clamping(self):
        bundle = self.bundle(np.random.default_rng(3).normal(size=(6, 3)))
        output = self.run_bundle(bundle, self.args(dimensions=(2, 99)))
        skipped = common.read_tsv(output / "skipped.tsv")
        oversized = next(r for r in skipped if r["selector"] == "dim:99")
        self.assertEqual(oversized["status"], "skipped")
        self.assertEqual(int(oversized["requested_components"]), 99)
        self.assertEqual(int(oversized["maximum_supported"]), 3)
        self.assertFalse((output / bundle["dataset"] / bundle["route"] / "pca_0099").exists())

    def test_zero_projected_rows_fail_whole_condition_preserving_complete_raw_run(self):
        bundle = self.bundle([[1., 0.], [-1., 0.], [0., 1.], [0., -1.]])
        output = self.run_bundle(bundle, self.args(dimensions=(1, 2)))
        metadata = common.read_json(output / "metadata.json")
        self.assertEqual(metadata["status"], "complete_with_condition_failures")
        conditions = common.read_tsv(output / "conditions.tsv")
        failed = next(r for r in conditions if r["condition"] == "pca_0001")
        self.assertEqual(failed["status"], "failed")
        self.assertIn("no rows were dropped", failed["reason"])
        matched = common.read_tsv(output / "matched_metrics.tsv")
        self.assertNotIn("pca_0001", {r["condition"] for r in matched})
        self.assertTrue(all(int(r["sample_count"]) == 4 and float(r["coverage"]) == 1 for r in matched))
        self.assertEqual(next(r for r in matched if r["condition"] == "raw")["sample_count"], "4")
        self.assertFalse((output / bundle["dataset"] / bundle["route"] / "pca_0001").exists())

    def test_nearest_count_then_lower_threshold_selection_not_max_v(self):
        rows = [{"n_clusters": 9, "true_gcf_count": 2, "threshold": 0., "v_measure": .99},
                {"n_clusters": 2, "true_gcf_count": 2, "threshold": .5, "v_measure": .4},
                {"n_clusters": 2, "true_gcf_count": 2, "threshold": .6, "v_measure": .8}]
        selected = pipeline.choose_matched_cut(rows)
        self.assertEqual(selected["threshold"], .5)
        self.assertEqual(selected["v_measure"], .4)

    def test_existing_output_is_rejected_and_original_report_preserved(self):
        bundle = self.bundle(np.random.default_rng(5).normal(size=(6, 3)))
        args = self.args(dimensions=(2,))
        output = self.run_bundle(bundle, args)
        original = (output / "summary.md").read_bytes()
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.run_bundle(bundle, args)
        self.assertEqual((output / "summary.md").read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
