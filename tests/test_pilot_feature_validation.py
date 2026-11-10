"""Reject incomplete or mismatched server feature exports without GPU/torch."""

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.pilot_feature_validation import file_sha256, validate_export


class PilotFeatureValidationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="pilot-feature-validation-")
        self.addCleanup(temporary.cleanup)
        self.output = Path(temporary.name)
        self.rows = [{"sample_id": "bgc_a", "source_sha256": "a" * 64},
                     {"sample_id": "bgc_b", "source_sha256": "b" * 64}]
        self.cache = self.output / "tokens.pth"
        self.cache.write_bytes(b"not deserialized by validation")
        np.save(self.output / "embeddings.npy", np.eye(2, dtype=np.float32))
        self.write_tsv("bgc_ids.tsv", ["row", "bgc_id", "token_count"],
                       [[0, "bgc_a", 3], [1, "bgc_b", 3]])
        self.write_tsv("failures.tsv", ["bgc_id", "source", "stage", "error"], [])
        self.metadata = {
            "format": "bgc_vectors_v1", "status": "complete", "route": "esm",
            "count": 2, "dimension": 2, "dtype": "float32", "normalization": "l2",
            "excluded_count": 0, "upstream_excluded_count": 0,
            "input": {"input_count": 2},
            "records": {row["sample_id"]: {"bgc_id": row["sample_id"],
                                           "source_sha256": row["source_sha256"]}
                        for row in self.rows},
            "saved_cache": {"sha256": file_sha256(self.cache)},
        }
        self.save_metadata()

    def write_tsv(self, name, columns, rows):
        with (self.output / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns)
            writer.writerows(rows)

    def save_metadata(self):
        (self.output / "metadata.json").write_text(json.dumps(self.metadata), encoding="utf-8")

    def validate(self, route="esm2", **kwargs):
        return validate_export(self.output, self.rows, route, **kwargs)

    def dna_metadata(self):
        self.metadata.update({
            "route": "dna", "sequence_mode": "full", "window_batch_size": 1,
            "length_exclusions": [], "length_excluded_count": 0, "invalid_count": 0,
            "input": {"input_count": 2, "length_filter": None, "explicit_exclusions": []},
            "extractor": {"sequence_mode": "full", "max_tokens": None,
                          "max_content_tokens": None, "config_length_enforced": False,
                          "automatic_window_fallback": False,
                          "window_policy": "full_sequence_single_forward"},
        })
        for record in self.metadata["records"].values():
            record.update({"sequence_mode": "full", "window_count": 1, "forward_pass_count": 1,
                           "token_count": 3, "content_token_count": 3, "base_count": 9,
                           "sequence_length": 9, "sequence_region": [0, 9],
                           "forward_input_shapes": [[1, 5]], "total_input_tokens": 5,
                           "max_forward_tokens": 5, "window_content_token_counts": [3],
                           "windows": [{"token_start": 0, "token_end": 3,
                                        "base_start": 0, "base_end": 9}]})
        self.save_metadata()

    def test_valid_esm_and_mac_share_hashed_cache(self):
        esm = self.validate(cache_path=self.cache)
        self.assertEqual((esm["count"], esm["dimension"]), (2, 2))
        self.assertEqual(set(esm["artifact_sha256"]),
                         {"embeddings.npy", "bgc_ids.tsv", "metadata.json", "failures.tsv"})
        self.metadata["route"] = "mac"
        self.metadata["input"]["cache_sha256"] = file_sha256(self.cache)
        self.save_metadata()
        mac = self.validate("mac", cache_path=self.cache)
        self.assertEqual(esm["cache_sha256"], mac["cache_sha256"])

    def test_valid_full_dna(self):
        self.dna_metadata()
        self.assertEqual(self.validate("dna_full")["count"], 2)

    def test_missing_artifact(self):
        (self.output / "failures.tsv").unlink()
        with self.assertRaisesRegex(ValueError, "Missing export artifacts"):
            self.validate()

    def test_missing_id_is_not_a_complete_export(self):
        self.write_tsv("bgc_ids.tsv", ["row", "bgc_id", "token_count"], [[0, "bgc_a", 3]])
        with self.assertRaisesRegex(ValueError, "complete sorted manifest"):
            self.validate()

    def test_wrong_order_or_duplicate_ids(self):
        for ids in [("bgc_b", "bgc_a"), ("bgc_a", "bgc_a")]:
            with self.subTest(ids=ids):
                self.write_tsv("bgc_ids.tsv", ["row", "bgc_id", "token_count"],
                               [[i, key, 3] for i, key in enumerate(ids)])
                with self.assertRaises(ValueError):
                    self.validate()

    def test_wrong_row_index(self):
        self.write_tsv("bgc_ids.tsv", ["row", "bgc_id", "token_count"],
                       [[1, "bgc_a", 3], [0, "bgc_b", 3]])
        with self.assertRaisesRegex(ValueError, "row indices"):
            self.validate()

    def test_fake_complete_with_failures_or_exclusions(self):
        for field in ("excluded_count", "upstream_excluded_count"):
            with self.subTest(field=field):
                self.metadata[field] = 1
                self.save_metadata()
                with self.assertRaisesRegex(ValueError, "exclusions"):
                    self.validate()
                self.metadata[field] = 0
        self.save_metadata()
        self.write_tsv("failures.tsv", ["bgc_id", "source", "stage", "error"],
                       [["bgc_b", "raw.gbk", "encoding", "failed"]])
        with self.assertRaisesRegex(ValueError, "failures.tsv contains"):
            self.validate()

    def test_nan_wrong_dtype_or_zero_norm(self):
        arrays = [np.array([[np.nan, 0], [0, 1]], dtype=np.float32),
                  np.eye(2, dtype=np.float64), np.zeros((2, 2), dtype=np.float32)]
        for array in arrays:
            with self.subTest(dtype=array.dtype, array=array.tolist()):
                np.save(self.output / "embeddings.npy", array)
                with self.assertRaises(ValueError):
                    self.validate()

    def test_cache_mismatch_for_both_protein_routes(self):
        self.cache.write_bytes(b"changed cache")
        with self.assertRaisesRegex(ValueError, "cache SHA256"):
            self.validate(cache_path=self.cache)
        self.metadata["route"] = "mac"
        self.metadata["input"]["cache_sha256"] = "0" * 64
        self.save_metadata()
        with self.assertRaisesRegex(ValueError, "cache SHA256"):
            self.validate("mac", cache_path=self.cache)

    def test_source_hash_or_record_set_mismatch(self):
        self.metadata["records"]["bgc_a"]["source_sha256"] = "0" * 64
        self.save_metadata()
        with self.assertRaisesRegex(ValueError, "Source SHA256"):
            self.validate()
        del self.metadata["records"]["bgc_a"]
        self.save_metadata()
        with self.assertRaisesRegex(ValueError, "Metadata records"):
            self.validate()

    def test_dna_rejects_length_filter_and_token_limit(self):
        self.dna_metadata()
        self.metadata["input"]["length_filter"] = {"max_bp": 40000}
        self.save_metadata()
        with self.assertRaisesRegex(ValueError, "length filter"):
            self.validate("dna_full")
        self.metadata["input"]["length_filter"] = None
        self.metadata["extractor"]["max_tokens"] = 2500
        self.save_metadata()
        with self.assertRaisesRegex(ValueError, "token limits"):
            self.validate("dna_full")

    def test_dna_rejects_windowed_or_truncated_forward(self):
        self.dna_metadata()
        self.metadata["records"]["bgc_b"]["forward_pass_count"] = 2
        self.save_metadata()
        with self.assertRaisesRegex(ValueError, "one full window"):
            self.validate("dna_full")
        self.metadata["records"]["bgc_b"]["forward_pass_count"] = 1
        self.metadata["records"]["bgc_b"]["forward_input_shapes"] = [[1, 4]]
        self.save_metadata()
        with self.assertRaisesRegex(ValueError, "forward shape"):
            self.validate("dna_full")


if __name__ == "__main__":
    unittest.main()
