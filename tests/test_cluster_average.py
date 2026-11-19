"""Meaningful checks on the cosine/UPGMA hierarchy and threshold semantics."""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from cluster_average import build_average_linkage, cut_average_linkage, threshold_grid


class AverageLinkageTests(unittest.TestCase):
    def test_uses_precomputed_distance_not_distance_rows_as_observations(self):
        distances = np.array([[0, .2, .4], [.2, 0, .6], [.4, .6, 0]])
        tree = build_average_linkage(distances)
        self.assertEqual(set(tree[0, :2]), {0, 1})
        np.testing.assert_allclose(tree[:, 2], [.2, .5])
        np.testing.assert_array_equal(tree[:, 3], [2, 3])

    def test_average_is_size_weighted_upgma_not_wpgma(self):
        distances = np.array([[0, .1, .2, .9], [.1, 0, .4, .8],
                              [.2, .4, 0, .6], [.9, .8, .6, 0]])
        tree = build_average_linkage(distances)
        np.testing.assert_allclose(tree[:, 2], [.1, .3, (.9 + .8 + .6) / 3])
        np.testing.assert_array_equal(tree[:, 3], [2, 3, 4])

    def test_actual_cosine_vectors_recover_expected_pairs(self):
        vectors = np.array([[1., 0.], [.99, .1], [-1., 0.], [-.99, -.1]])
        vectors /= np.linalg.norm(vectors, axis=1, keepdims=True)
        distances = np.clip(1 - vectors @ vectors.T, 0, 2)
        np.fill_diagonal(distances, 0)
        tree = build_average_linkage(distances)
        labels = cut_average_linkage(tree, .02)
        np.testing.assert_array_equal(labels, [0, 0, 1, 1])
        np.testing.assert_array_equal(cut_average_linkage(tree, 2), [0, 0, 0, 0])

    def test_cut_includes_equal_height_and_does_not_force_requested_k(self):
        distances = np.ones((4, 4)) - np.eye(4)
        tree = build_average_linkage(distances)
        np.testing.assert_array_equal(threshold_grid(tree), [0, 1])
        np.testing.assert_array_equal(cut_average_linkage(tree, np.nextafter(1., 0)), [0, 1, 2, 3])
        np.testing.assert_array_equal(cut_average_linkage(tree, 1), [0, 0, 0, 0])

    def test_identical_vectors_have_one_attainable_zero_distance_partition(self):
        tree = build_average_linkage(np.zeros((4, 4)))
        np.testing.assert_array_equal(threshold_grid(tree), [0])
        np.testing.assert_array_equal(cut_average_linkage(tree, 0), [0, 0, 0, 0])

    def test_explicit_thresholds_sorted_deduplicated_without_label_use(self):
        tree = build_average_linkage(np.array([[0, .4], [.4, 0]]))
        np.testing.assert_array_equal(threshold_grid(tree, [2, .1, .1, 0]), [0, .1, 2])
        for values in ([], [-.1], [2.01], [float("nan")], [float("inf")]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                threshold_grid(tree, values)

    def test_invalid_distance_matrix_is_rejected(self):
        cases = [np.array([[0.]]), np.zeros((2, 3)), np.array([[0, -1], [-1, 0]]),
                 np.array([[0, 3], [3, 0]]), np.array([[0, .1], [.2, 0]]),
                 np.array([[1, .1], [.1, 0]]), np.array([[0, np.nan], [np.nan, 0]])]
        for distances in cases:
            with self.subTest(distances=distances), self.assertRaises(ValueError):
                build_average_linkage(distances)


if __name__ == "__main__":
    unittest.main()
