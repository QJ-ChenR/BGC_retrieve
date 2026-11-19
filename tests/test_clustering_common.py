"""Feature/label identity and numerical safeguards shared by both baselines."""

import argparse
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import clustering_common as common


class ClusteringInputTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="clustering-input-test-")
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.dataset = "Dataset_B_Divergent"
        self.directory = self.root / self.dataset / "esm2"
        self.directory.mkdir(parents=True)
        self.rows = [
            {"sample_id": "bgc_a", "source_sha256": "a" * 64,
             "dataset": self.dataset, "gcf_raw": "family\u00a0", "original_bgc_id": "original_a"},
            {"sample_id": "bgc_b", "source_sha256": "b" * 64,
             "dataset": self.dataset, "gcf_raw": "neg_b", "original_bgc_id": "original_b"},
        ]
        # Matrix rows differ from canonical IDs, and TSV physical order differs from matrix order.
        np.save(self.directory / "embeddings.npy", np.array([[0, 1], [1, 0]], dtype=np.float32))
        self.save_tsv("bgc_ids.tsv", ["row", "bgc_id", "token_count"],
                      [[1, "bgc_a", 4], [0, "bgc_b", 5]])
        self.save_tsv("failures.tsv", ["bgc_id", "source", "stage", "error"], [])
        self.metadata = {
            "format": "bgc_vectors_v1", "status": "complete", "route": "esm",
            "count": 2, "dimension": 2, "dtype": "float32", "normalization": "l2",
            "excluded_count": 0, "upstream_excluded_count": 0, "input": {"input_count": 2},
            "saved_cache": {"sha256": "c" * 64},
            "records": {row["sample_id"]: {"bgc_id": row["sample_id"],
                                           "source_sha256": row["source_sha256"],
                                           "source": "/server/path/not/present/locally.gbk"}
                        for row in self.rows},
        }
        self.save_metadata()

    def save_tsv(self, name, fields, records):
        with (self.directory / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(fields)
            writer.writerows(records)

    def save_metadata(self):
        (self.directory / "metadata.json").write_text(json.dumps(self.metadata), encoding="utf-8")

    def load(self):
        return common.load_feature_bundle(self.directory, list(reversed(self.rows)), self.dataset, "esm2")

    def test_align_by_explicit_row_then_sample_id_preserve_raw_unicode(self):
        bundle = self.load()
        self.assertEqual(bundle["ids"], ["bgc_a", "bgc_b"])
        self.assertEqual(bundle["source_row_indices"], [1, 0])
        np.testing.assert_array_equal(bundle["vectors"], np.eye(2))
        self.assertEqual(bundle["labels"], ["family\u00a0", "neg_b"])

    def test_missing_id_duplicate_id_and_duplicate_matrix_row_fail(self):
        cases = [ [[1, "bgc_a", 4]], [[0, "bgc_a", 4], [1, "bgc_a", 5]],
                  [[0, "bgc_a", 4], [0, "bgc_b", 5]] ]
        for rows in cases:
            with self.subTest(rows=rows):
                self.save_tsv("bgc_ids.tsv", ["row", "bgc_id", "token_count"], rows)
                with self.assertRaises(ValueError):
                    self.load()

    def test_partial_status_or_failure_records_are_not_silently_intersected(self):
        self.metadata["status"] = "complete_with_exclusions"
        self.save_metadata()
        with self.assertRaisesRegex(ValueError, "incomplete"):
            self.load()
        self.metadata["status"] = "complete"
        self.save_metadata()
        self.save_tsv("failures.tsv", ["bgc_id", "source", "stage", "error"],
                      [["bgc_a", "path", "encoding", "OOM"]])
        with self.assertRaisesRegex(ValueError, "failure records"):
            self.load()

    def test_hash_mismatch_rejected(self):
        self.metadata["records"]["bgc_a"]["source_sha256"] = "0" * 64
        self.save_metadata()
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            self.load()

    def test_nan_zero_or_unnormalized_vectors_rejected(self):
        for matrix in (np.array([[1, 0], [np.nan, 0]]), np.zeros((2, 2)), 2 * np.eye(2)):
            with self.subTest(matrix=matrix.tolist()):
                np.save(self.directory / "embeddings.npy", matrix.astype(np.float32))
                with self.assertRaises(ValueError):
                    self.load()

    def test_extraction_run_hash_mismatch_rejected(self):
        (self.root / "run.json").write_text(json.dumps({"steps": {
            self.dataset + "/esm2": {"status": "complete", "validation": {"artifact_sha256": {}}}
        }}))
        args = argparse.Namespace(features_root=self.root, datasets=["B"], routes=["esm2"],
                                  cohort=self.root / "elsewhere" / "cohort.tsv",
                                  output=self.root.parent / (self.root.name + "-results"))
        with mock.patch.object(common, "load_cohort", return_value=self.rows):
            with self.assertRaisesRegex(ValueError, "extraction validation"):
                common.load_inputs(args)


class ClusteringMathTests(unittest.TestCase):
    def test_cosine_geometry_preserved_without_projection(self):
        unit, distance = common.cosine_distances([[2, 0], [0, 3], [-4, 0]])
        np.testing.assert_allclose(np.linalg.norm(unit, axis=1), 1)
        np.testing.assert_array_equal(distance, [[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        with self.assertRaisesRegex(ValueError, "zero"):
            common.cosine_distances([[0, 0], [1, 1]])

    def test_perfect_partition_and_undefined_silhouettes(self):
        _, distance = common.cosine_distances([[1, 0], [1, .01], [0, 1], [.01, 1]])
        result = common.score_partition(["a", "a", "b", "b"], [0, 0, 1, 1], distance)
        self.assertEqual(result["ari"], 1)
        self.assertEqual(result["v_measure"], 1)
        self.assertEqual(result["coverage"], 1)
        self.assertGreater(result["silhouette_cosine"], .9)
        for membership in ([0, 0, 0, 0], [0, 1, 2, 3]):
            result = common.score_partition(["a", "a", "b", "b"], membership, distance)
            self.assertIsNone(result["silhouette_cosine"])
        with self.assertRaisesRegex(ValueError, "every sample"):
            common.score_partition(["a", "a", "b", "b"], [-1, 0, 1, 2], distance)


if __name__ == "__main__":
    unittest.main()
