"""Label-free, centered PCA and explicit dimension planning for cosine analysis.

Fit each dataset/representation separately. Inputs are used as supplied: there is
no feature scaling, whitening, or row normalization inside ``fit_pca``. Cosine
normalization happens only after projection and cannot silently drop samples.
"""

from __future__ import annotations

from numbers import Integral

import numpy as np

from clustering_common import cosine_distances


def fit_pca(vectors):
    """Fit deterministic full-SVD PCA, retaining only the numerical rank.

    Returned components are rows; transform another array with
    ``(array - fit['mean']) @ fit['components'].T``. Component signs are fixed by
    making the largest absolute loading positive. Repeated singular values can
    still admit different bases across LAPACK implementations, which preserves
    full-subspace geometry but can affect a cut inside a repeated eigenspace.
    """
    values = np.asarray(vectors, dtype=np.float64)
    if (values.ndim != 2 or values.shape[0] < 2 or values.shape[1] < 1
            or not np.isfinite(values).all()):
        raise ValueError("PCA requires a finite N x D matrix with N >= 2 and D >= 1")
    if np.all(values == values[0]):
        raise ValueError("PCA is undefined for an all-constant sample matrix")
    sample_count, input_dimension = values.shape
    mean = values.mean(axis=0)
    centered = values - mean
    if not np.isfinite(centered).all():
        raise ValueError("PCA centering produced nonfinite values")
    _, singular, components = np.linalg.svd(centered, full_matrices=False)
    theoretical_max_rank = min(sample_count - 1, input_dimension)
    rank_tolerance = float(max(centered.shape) * np.finfo(np.float64).eps * singular[0])
    numerical_rank = min(int(np.count_nonzero(singular > rank_tolerance)), theoretical_max_rank)
    if numerical_rank == 0:
        raise ValueError("PCA has zero numerical rank after centering")
    all_variance = np.square(singular) / (sample_count - 1)
    total_variance = float(all_variance.sum())
    if not np.isfinite(total_variance) or total_variance <= 0:
        raise ValueError("PCA total variance must be positive and finite")
    components = components[:numerical_rank].copy()
    pivots = np.argmax(np.abs(components), axis=1)
    signs = np.sign(components[np.arange(numerical_rank), pivots])
    components *= signs[:, None]
    scores = centered @ components.T
    return {
        "mean": mean,
        "components": components,
        "singular_values": singular[:numerical_rank].copy(),
        "explained_variance": all_variance[:numerical_rank].copy(),
        "explained_variance_ratio": all_variance[:numerical_rank] / total_variance,
        "scores": scores,
        "theoretical_max_rank": theoretical_max_rank,
        "numerical_rank": numerical_rank,
        "rank_tolerance": rank_tolerance,
        "input_dimension": input_dimension,
        "sample_count": sample_count,
    }


def plan_dimensions(fit, dimensions, variance_targets):
    """Resolve requested dimensions without silently clamping impossible ones.

    All selectors resolving to the same number of components share one plan.
    A full-numerical-rank control is always included. Selection uses the input
    covariance only, never family labels or clustering metrics.
    """
    rank = int(fit["numerical_rank"])
    ratios = np.asarray(fit["explained_variance_ratio"], dtype=np.float64)
    if (rank < 1 or ratios.shape != (rank,) or not np.isfinite(ratios).all()
            or np.any(ratios < 0) or not np.isclose(ratios.sum(), 1.0, rtol=1e-10, atol=1e-12)):
        raise ValueError("Invalid PCA rank or explained-variance ratios")
    cumulative = np.minimum(np.cumsum(ratios), 1.0)
    by_dimension, skips = {}, []

    def add(count, selector):
        selectors = by_dimension.setdefault(count, [])
        if selector not in selectors:
            selectors.append(selector)

    for count in dimensions:
        if isinstance(count, bool) or not isinstance(count, Integral) or count < 1:
            raise ValueError("Requested PCA dimensions must be positive integers")
        count = int(count)
        selector = f"dim:{count}"
        if count > rank:
            item = {"selector": selector, "requested_components": count,
                    "maximum_supported": rank,
                    "reason": "requested dimension exceeds the centered numerical rank"}
            if item not in skips:
                skips.append(item)
        else:
            add(count, selector)
    for target in variance_targets:
        if isinstance(target, bool):
            raise ValueError("Explained-variance targets must satisfy 0 < target <= 1")
        target = float(target)
        if not np.isfinite(target) or not 0 < target <= 1:
            raise ValueError("Explained-variance targets must satisfy 0 < target <= 1")
        # The numerical-rank tail can leave sum(ratios) infinitesimally below 1.
        count = rank if target == 1 else min(int(np.searchsorted(cumulative, target, side="left")) + 1, rank)
        add(count, f"variance:{target:g}")
    add(rank, "full_rank")
    plans = [{"condition": f"pca_{count:04d}", "n_components": count,
              "selectors": selectors, "retained_variance": float(cumulative[count - 1])}
             for count, selectors in sorted(by_dimension.items())]
    return plans, skips


def normalized_projection(scores):
    """Return normalized projected rows and cosine distances, or reject all.

    A row is numerically zero if its norm is no greater than
    ``eps_float64 * max(N, D) * max_row_norm``. This scale-relative check prevents
    numerical residue at the origin from becoming an arbitrary unit direction.
    Any offending row invalidates the condition; no sample is excluded.
    """
    values = np.asarray(scores, dtype=np.float64)
    if values.ndim != 2 or min(values.shape) < 1 or not np.isfinite(values).all():
        raise ValueError("Projected scores must be a finite nonempty N x D matrix")
    norms = np.linalg.norm(values, axis=1)
    if not np.isfinite(norms).all():
        raise ValueError("Projected row norms must be finite")
    tolerance = float(np.finfo(np.float64).eps * max(values.shape) * norms.max())
    bad_rows = np.flatnonzero(norms <= tolerance)
    if len(bad_rows):
        raise ValueError(f"Cosine is undefined for zero/near-zero projected rows {bad_rows.tolist()}; "
                         f"norm tolerance={tolerance:.17g}; no rows were dropped")
    return cosine_distances(values)
