"""Guard the comparison against selecting parameters by test V or best seed."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from summarize_vmeasure import aggregate_leiden, select_summaries, singleton_controls


class SummarySelectionTests(unittest.TestCase):
    def average_row(self, clusters, threshold, v):
        return {"dataset": "B", "route": "esm2", "run_id": str(threshold),
                "sample_count": 39, "true_gcf_count": 23, "n_clusters": clusters,
                "threshold": threshold, "v_measure": v, "ari": .1,
                "homogeneity": v, "completeness": v}

    def leiden_row(self, resolution, seed, clusters, v):
        return {"dataset": "B", "route": "esm2", "k": 10, "resolution": resolution,
                "seed": seed, "sample_count": 39, "true_gcf_count": 23,
                "n_clusters": clusters, "v_measure": v, "ari": .1,
                "homogeneity": v, "completeness": v}

    def test_average_prefers_family_count_and_threshold_not_high_v(self):
        rows = [self.average_row(39, 0, .95), self.average_row(23, .1, .7),
                self.average_row(23, .2, .8)]
        selected, _, oracle = select_summaries(rows, [])
        self.assertEqual(selected[0]["threshold"], .1)
        self.assertEqual(selected[0]["v_measure"], .7)
        self.assertEqual(oracle[0]["v_measure_mean"], .95)

    def test_leiden_uses_mean_absolute_count_error_and_all_seeds(self):
        rows = [self.leiden_row(1, 0, 13, .95), self.leiden_row(1, 1, 33, .9),
                self.leiden_row(2, 0, 22, .5), self.leiden_row(2, 1, 22, .7)]
        grouped = aggregate_leiden(rows)
        _, selected, _ = select_summaries([], grouped)
        self.assertEqual(selected[0]["resolution"], 2)
        self.assertAlmostEqual(selected[0]["v_measure_mean"], .6)
        self.assertAlmostEqual(selected[0]["v_measure_std"], .1)
        self.assertEqual(selected[0]["seed_count"], 2)

    def test_insufficient_resolution_range_is_explicit(self):
        rows = [self.leiden_row(1, 0, 4, .4), self.leiden_row(2, 0, 8, .5)]
        _, selected, _ = select_summaries([], aggregate_leiden(rows))
        self.assertFalse(selected[0]["target_bracketed_in_grid"])
        self.assertEqual(selected[0]["grid_n_clusters_max"], 8)

    def test_singleton_control_has_zero_ari_despite_positive_v(self):
        controls = singleton_controls([{"dataset": "B", "gcf_raw": x} for x in ["a", "a", "b", "b"]])
        self.assertGreater(controls[0]["v_measure"], .6)
        self.assertEqual(controls[0]["ari"], 0)


if __name__ == "__main__":
    unittest.main()
