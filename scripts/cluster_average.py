#!/usr/bin/env python3
"""Cluster native BGC exports using original cosine distance and average linkage.

The hierarchy and distance cuts use vectors only. Curated family labels are used
after clustering for evaluation, never to choose a cut or the number of clusters.
"""

from __future__ import annotations

import argparse

import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform

import clustering_common as common


METHOD = "cosine_average_linkage"


def build_average_linkage(distances):
    """Build UPGMA from a *condensed* original cosine distance matrix.

    Passing the square matrix directly to ``linkage`` would instead treat its
    rows as observation vectors, yielding a different and unintended distance.
    """
    distances = np.asarray(distances, dtype=np.float64)
    if distances.ndim != 2 or distances.shape[0] != distances.shape[1] or len(distances) < 2:
        raise ValueError("Average linkage requires an N x N distance matrix with N >= 2")
    if not np.isfinite(distances).all() or np.any(distances < 0) or np.any(distances > 2):
        raise ValueError("Cosine distances must be finite and lie in [0, 2]")
    if not np.allclose(distances, distances.T, atol=1e-12, rtol=0):
        raise ValueError("Distance matrix must be symmetric")
    if not np.allclose(np.diag(distances), 0, atol=1e-12, rtol=0):
        raise ValueError("Distance matrix must have a zero diagonal")
    condensed = squareform(distances, checks=False)
    return linkage(condensed, method="average", optimal_ordering=False)


def threshold_grid(tree, requested=None):
    """Use all distinct merge heights by default, plus the zero-distance cut.

    This is an input-dependent diagnostic sweep, not test-label tuning. Tied
    heights are cut together, so not every possible cluster count is attainable.
    """
    values = np.concatenate(([0.0], np.asarray(tree)[:, 2])) if requested is None else np.asarray(requested)
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError("Thresholds must be a nonempty finite one-dimensional list")
    if np.any(values < 0) or np.any(values > 2):
        raise ValueError("Cosine distance thresholds must lie in [0, 2]")
    return np.unique(values)


def cut_average_linkage(tree, threshold):
    if not np.isfinite(threshold) or not 0 <= threshold <= 2:
        raise ValueError("Cosine distance threshold must be finite and lie in [0, 2]")
    return common.canonical_labels(fcluster(tree, t=float(threshold), criterion="distance"))


def plot_outputs(directory, bundle, tree, metrics):
    """Save publication-independent diagnostic figures; no projection is used."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    count = len(bundle["ids"])
    fig, axis = plt.subplots(figsize=(16 if count > 60 else 13, 7))
    dendrogram(
        tree, ax=axis,
        labels=[row["original_bgc_id"] for row in bundle["rows"]],
        no_labels=count > 60, leaf_rotation=90, leaf_font_size=7,
        link_color_func=lambda _: "#334155",
    )
    axis.set_ylabel("Average-linkage cosine distance")
    axis.set_xlabel("BGC leaves (complete mapping: leaf_order.tsv)")
    axis.set_title(f"{bundle['dataset']} / {bundle['route']} — original cosine + average linkage")
    fig.tight_layout()
    for extension in ("svg", "png"):
        fig.savefig(directory / f"dendrogram.{extension}", dpi=180)
    plt.close(fig)

    thresholds = [row["threshold"] for row in metrics]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for key, label in (("v_measure", "V-measure"), ("ari", "ARI"),
                       ("homogeneity", "Homogeneity"), ("completeness", "Completeness")):
        axes[0].plot(thresholds, [row[key] for row in metrics], label=label)
    axes[0].set_xlabel("Cosine distance cut")
    axes[0].set_ylabel("Curated-GCF agreement")
    axes[0].legend(frameon=False)
    axes[1].step(thresholds, [row["n_clusters"] for row in metrics], where="post", label="Clusters")
    axes[1].step(thresholds, [row["n_singleton_clusters"] for row in metrics], where="post", label="Singleton clusters")
    axes[1].set_xlabel("Cosine distance cut")
    axes[1].set_ylabel("Count")
    axes[1].legend(frameon=False)
    fig.suptitle(f"{bundle['dataset']} / {bundle['route']} — diagnostic sweep; no selected optimum")
    fig.tight_layout()
    for extension in ("svg", "png"):
        fig.savefig(directory / f"threshold_curve.{extension}", dpi=180)
    plt.close(fig)


def run(args):
    # Validate explicit thresholds and all input routes before creating outputs.
    if args.thresholds is not None:
        threshold_grid(np.empty((0, 4)), args.thresholds)
    bundles = common.load_inputs(args)
    if any(len(bundle["ids"]) < 2 for bundle in bundles):
        raise ValueError("Average linkage needs at least two BGCs per dataset")
    output = common.create_output(args.output)
    all_metrics = []
    for bundle in bundles:
        directory = output / bundle["dataset"] / bundle["route"]
        directory.mkdir(parents=True)
        tree = build_average_linkage(bundle["distances"])
        np.save(directory / "linkage.npy", tree, allow_pickle=False)
        common.write_tsv(
            directory / "linkage_rows.tsv",
            [{"merge_index": index, "new_node": len(bundle["ids"]) + index,
              "left_node": int(row[0]), "right_node": int(row[1]),
              "distance": float(row[2]), "member_count": int(row[3])}
             for index, row in enumerate(tree)],
        )
        # Tree node indices < N refer to the canonical input order, not leaf order.
        common.write_tsv(directory / "input_order.tsv", [
            {"row": index, "sample_id": key, "original_bgc_id": bundle["rows"][index]["original_bgc_id"]}
            for index, key in enumerate(bundle["ids"])
        ])
        common.write_tsv(directory / "leaf_order.tsv", [
            {"leaf_position": position, "row": int(index), "sample_id": bundle["ids"][index],
             "original_bgc_id": bundle["rows"][index]["original_bgc_id"],
             "gcf_raw": bundle["labels"][index]}
            for position, index in enumerate(leaves_list(tree))
        ])
        runs, metrics = [], []
        for index, threshold in enumerate(threshold_grid(tree, args.thresholds)):
            run_id = f"cut_{index:04d}"
            membership = cut_average_linkage(tree, float(threshold))
            runs.append((run_id, membership))
            metrics.append({
                "dataset": bundle["dataset"], "route": bundle["route"], "method": METHOD,
                "run_id": run_id, "threshold": float(threshold),
                **common.score_partition(bundle["labels"], membership, bundle["distances"]),
            })
        common.write_assignments(directory / "clusters.tsv", bundle, runs)
        common.write_tsv(directory / "metrics.tsv", metrics)
        if not args.no_plots:
            plot_outputs(directory, bundle, tree, metrics)
        all_metrics.extend(metrics)
        print(f"{bundle['dataset']}/{bundle['route']}: {len(bundle['ids'])} BGCs, "
              f"{len(metrics)} distance cuts", flush=True)
    common.finish_run(output, args, bundles, METHOD, {
        "distance": "cosine in original embedding space",
        "linkage": "average (UPGMA), scipy condensed precomputed distances",
        "optimal_leaf_ordering": False,
        "input_order": "ascending stable sample_id; no family-label ordering",
        "threshold_policy": "all distinct merge heights plus zero" if args.thresholds is None else "explicit user thresholds",
        "requested_thresholds": args.thresholds,
        "cut_rule": "scipy fcluster criterion=distance; include merges at height <= threshold",
        "ties": "equal-height merges are cut together; not every cluster count is attainable",
        "merge_tie_breaking": "SciPy implementation on canonical sample-ID order; library version recorded",
        "label_use": "evaluation only; no cutoff or optimum selected from labels",
        "projection": None,
        "references": [
            "https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html",
            "https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html",
        ],
    }, all_metrics)
    return output


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    common.add_common_arguments(parser, default_output="data/processed/clustering_average_v1")
    parser.add_argument("--thresholds", type=float, nargs="+",
                        help="Cosine-distance cuts in [0,2]; default: zero and all distinct merge heights")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        output = run(args)
    except (ValueError, OSError, KeyError, RuntimeError, ImportError) as exc:
        parser.exit(1, f"Average-linkage analysis failed: {exc}\n")
    print(f"Saved average-linkage analysis to {output}")


if __name__ == "__main__":
    main()
