"""Numerical and cohort-preserving guarantees for PCA/cosine experiments."""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from clustering_common import cosine_distances
from pca_transform import fit_pca, normalized_projection, plan_dimensions


class PCATransformTests(unittest.TestCase):
    def setUp(self):
        self.values = np.array([[2., 1., 4.], [3., -2., 5.], [-1., 3., 2.],
                                [4., 0., 1.], [0., -1., 7.]])

    def test_full_rank_preserves_centered_cosine_not_original_cosine(self):
        fitted = fit_pca(self.values)
        _, full_distances = normalized_projection(fitted["scores"])
        _, centered_distances = cosine_distances(self.values - self.values.mean(axis=0))
        _, raw_distances = cosine_distances(self.values)
        np.testing.assert_allclose(full_distances, centered_distances, atol=1e-12)
        self.assertGreater(np.max(np.abs(full_distances - raw_distances)), .1)

    def test_saved_components_transform_and_reconstruct_without_whitening(self):
        fitted = fit_pca(self.values)
        centered = self.values - fitted["mean"]
        np.testing.assert_allclose(centered @ fitted["components"].T, fitted["scores"], atol=1e-12)
        np.testing.assert_allclose(fitted["scores"] @ fitted["components"], centered, atol=1e-12)
        np.testing.assert_allclose(np.var(fitted["scores"], axis=0, ddof=1),
                                   fitted["explained_variance"], atol=1e-12)
        np.testing.assert_allclose(fit_pca(self.values)["components"], fitted["components"], atol=0)

    def test_rank_is_limited_by_samples_features_and_numerical_dependencies(self):
        rng = np.random.default_rng(11)
        for shape in ((5, 20), (20, 5)):
            fitted = fit_pca(rng.normal(size=shape))
            self.assertEqual(fitted["theoretical_max_rank"], min(shape[0] - 1, shape[1]))
            self.assertEqual(fitted["numerical_rank"], min(shape[0] - 1, shape[1]))
        dependent = np.arange(4.)[:, None] * np.array([[1., 2., -3.]])
        self.assertEqual(fit_pca(dependent)["numerical_rank"], 1)

    def test_oversized_dimension_is_logged_not_clamped(self):
        fitted = fit_pca(self.values)
        plans, skipped = plan_dimensions(fitted, [2, 32, 32], [])
        self.assertEqual([p["n_components"] for p in plans], [2, 3])
        self.assertEqual(len(skipped), 1)
        self.assertEqual(skipped[0]["requested_components"], 32)
        self.assertEqual(skipped[0]["maximum_supported"], 3)
        self.assertEqual(plans[-1]["selectors"], ["full_rank"])

    def test_variance_targets_choose_minimal_dimension_and_share_aliases(self):
        fitted = {"numerical_rank": 3, "explained_variance_ratio": np.array([.6, .3, .1])}
        plans, skipped = plan_dimensions(fitted, [1, 2], [.5, .8, .95, 1.0])
        self.assertFalse(skipped)
        self.assertEqual([p["n_components"] for p in plans], [1, 2, 3])
        self.assertEqual(plans[0]["selectors"], ["dim:1", "variance:0.5"])
        self.assertEqual(plans[1]["selectors"], ["dim:2", "variance:0.8"])
        self.assertIn("full_rank", plans[2]["selectors"])
        self.assertAlmostEqual(plans[2]["retained_variance"], 1)

    def test_sample_at_mean_invalidates_projection_without_exclusion(self):
        fitted = fit_pca(np.array([[1., 2.], [2., 4.], [3., 6.]]))
        with self.assertRaisesRegex(ValueError, "rows \\[1\\].*no rows were dropped"):
            normalized_projection(fitted["scores"])
        with self.assertRaisesRegex(ValueError, "near-zero"):
            normalized_projection(np.array([[1., 0.], [1e-18, 1e-18]]))
        unit, _ = normalized_projection(np.array([[1., 2.], [3., 4.]]) * 1e-12)
        np.testing.assert_allclose(np.linalg.norm(unit, axis=1), 1)

    def test_invalid_inputs_and_dimension_requests_are_rejected(self):
        for values in ([1, 2], [[1, 2]], np.ones((4, 3)), [[1., 2.], [np.nan, 3.]]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                fit_pca(values)
        fitted = fit_pca(self.values)
        for dimensions in ([0], [-1], [1.5], [True]):
            with self.subTest(dimensions=dimensions), self.assertRaises(ValueError):
                plan_dimensions(fitted, dimensions, [])
        for targets in ([0], [1.1], [float("nan")], [True]):
            with self.subTest(targets=targets), self.assertRaises(ValueError):
                plan_dimensions(fitted, [], targets)


if __name__ == "__main__":
    unittest.main()
