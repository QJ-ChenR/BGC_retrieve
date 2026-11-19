"""Graph semantics, retained isolates and stochastic repeatability."""

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from cluster_leiden import build_knn_graph, run_leiden, seed_stability


class KnnGraphTests(unittest.TestCase):
    def test_identical_vectors_exclude_self_and_use_index_ties(self):
        graph = build_knn_graph(np.zeros((4, 4)), 1)
        self.assertEqual(graph["neighbors"].tolist(), [[1], [0], [0], [0]])
        self.assertEqual(graph["edges"], [(0, 1), (0, 2), (0, 3)])
        self.assertEqual(graph["weights"], [1, 1, 1])
        self.assertEqual(graph["stats"]["directed_neighbor_count"], 4)
        self.assertEqual(graph["stats"]["edge_count"], 3)

    def test_union_keeps_reciprocal_weight_once_and_mutual_removes_one_way(self):
        distances = np.array([[0, 0.1, 0.2], [0.1, 0, 0.8], [0.2, 0.8, 0]])
        union = build_knn_graph(distances, 1, "union")
        mutual = build_knn_graph(distances, 1, "mutual")
        self.assertEqual(union["edges"], [(0, 1), (0, 2)])
        np.testing.assert_allclose(union["weights"], [0.9, 0.8])
        self.assertEqual(mutual["edges"], [(0, 1)])
        self.assertEqual(mutual["isolated_vertices"], [2])
        self.assertEqual(mutual["stats"]["connected_component_count"], 2)

    def test_nonpositive_edges_omitted_but_every_vertex_retained(self):
        distances = np.array([[0, 0.1, 1.2], [0.1, 0, 1.4], [1.2, 1.4, 0]])
        graph = build_knn_graph(distances, 2)
        self.assertEqual(graph["n"], 3)
        self.assertEqual(graph["edges"], [(0, 1)])
        self.assertEqual(graph["isolated_vertices"], [2])
        self.assertEqual(graph["stats"]["nonpositive_candidate_edges_removed"], 2)

    def test_empty_graph_keeps_all_samples_as_singletons(self):
        graph = build_knn_graph(np.ones((4, 4)) - np.eye(4), 2)
        membership, details = run_leiden(graph, resolution=1.0, seed=0)
        self.assertEqual(membership.tolist(), list(range(4)))
        self.assertEqual(details["quality_reason"], "no_positive_edges")
        self.assertIsNone(details["quality"])
        self.assertEqual(graph["stats"]["connected_component_count"], 4)

    def test_invalid_matrix_and_neighbor_sizes_rejected(self):
        valid = np.ones((3, 3)) - np.eye(3)
        for k in (0, 3, 100, True):
            with self.subTest(k=k), self.assertRaises(ValueError):
                build_knn_graph(valid, k)
        for invalid in (np.array([[0, np.nan], [np.nan, 0]]),
                        np.array([[0, -1], [-1, 0]]),
                        np.array([[0, 0.2], [0.3, 0]]), np.eye(3)):
            with self.subTest(matrix=invalid), self.assertRaises(ValueError):
                build_knn_graph(invalid, 1)


@unittest.skipUnless(importlib.util.find_spec("igraph") and importlib.util.find_spec("leidenalg"),
                     "Leiden dependencies are not installed")
class LeidenTests(unittest.TestCase):
    def test_reproducible_partition_and_isolate_singleton(self):
        distances = np.full((7, 7), 1.2)
        distances[:3, :3] = 0.05
        distances[3:6, 3:6] = 0.05
        np.fill_diagonal(distances, 0)
        graph = build_knn_graph(distances, 2)
        first, details = run_leiden(graph, 1.0, 42)
        second, _ = run_leiden(graph, 1.0, 42)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(len(first), 7)
        self.assertEqual(len(set(first[:3])), 1)
        self.assertEqual(len(set(first[3:6])), 1)
        self.assertNotEqual(first[0], first[3])
        self.assertNotIn(first[6], first[:6])
        self.assertEqual(details["isolates_as_singletons"], 1)
        self.assertTrue(np.isfinite(details["quality"]))


@unittest.skipUnless(importlib.util.find_spec("sklearn"), "scikit-learn is not installed")
class StabilityTests(unittest.TestCase):
    def test_stability_is_partition_invariant_and_reports_single_seed(self):
        pairs, summary = seed_stability({1: [0, 0, 1], 7: [4, 4, 2]})
        self.assertEqual(pairs, [{"seed_a": 1, "seed_b": 7, "pairwise_ari": 1.0}])
        self.assertEqual(summary["pairwise_ari_mean"], 1.0)
        self.assertEqual(summary["pairwise_ari_std"], 0.0)
        pairs, summary = seed_stability({1: [0, 0, 1]})
        self.assertFalse(pairs)
        self.assertIsNone(summary["pairwise_ari_mean"])
        self.assertEqual(summary["stability_reason"], "requires_at_least_two_seeds")


if __name__ == "__main__":
    unittest.main()
