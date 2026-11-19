#!/usr/bin/env python3
"""Cluster BGC vectors with exact cosine kNN graphs and weighted Leiden.

No PCA/UMAP or curated labels enter graph construction or community detection.
The graph keeps every sample, with each isolated vertex reported as a singleton.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import math
from pathlib import Path

import numpy as np


def build_knn_graph(distances, k, symmetrize="union"):
    """Return a reproducible undirected positive-cosine weighted kNN graph.

    Rows must already be in canonical sample-ID order. Self-neighbors are never
    candidates, even when multiple rows have identical vectors. Equal distances
    are broken by row index. Reciprocal neighbor choices yield just one edge.
    """
    distances = np.asarray(distances, dtype=np.float64)
    if (distances.ndim != 2 or distances.shape[0] != distances.shape[1]
            or distances.shape[0] < 2 or not np.isfinite(distances).all()):
        raise ValueError("Expected a finite square cosine distance matrix for at least two samples")
    if (np.any(distances < 0) or np.any(distances > 2)
            or not np.allclose(distances, distances.T, atol=1e-12, rtol=0)
            or not np.allclose(np.diag(distances), 0, atol=1e-12, rtol=0)):
        raise ValueError("Cosine distances must be symmetric, in [0, 2], with a zero diagonal")
    n = distances.shape[0]
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or not 1 <= k < n:
        raise ValueError(f"k must satisfy 1 <= k < sample_count ({n}); no automatic clipping")
    if symmetrize not in {"union", "mutual"}:
        raise ValueError("symmetrize must be union or mutual")
    indices = np.arange(n)
    neighbors = np.empty((n, k), dtype=np.int64)
    for i in range(n):
        candidates = indices[indices != i]
        order = np.lexsort((candidates, distances[i, candidates]))
        neighbors[i] = candidates[order[:k]]
    directed = {(i, int(j)) for i in range(n) for j in neighbors[i]}
    candidates = sorted({(min(i, j), max(i, j)) for i, j in directed
                         if symmetrize == "union" or (j, i) in directed})
    edges, weights = [], []
    for i, j in candidates:
        weight = max(0.0, float(1.0 - distances[i, j]))
        if weight > 0:
            edges.append((i, j))
            weights.append(weight)

    # Union-find gives graph diagnostics without requiring igraph at import time.
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    degrees = np.zeros(n, dtype=np.int64)
    for i, j in edges:
        parent[find(j)] = find(i)
        degrees[i] += 1
        degrees[j] += 1
    isolated = np.flatnonzero(degrees == 0).tolist()
    stats = {
        "vertex_count": n,
        "k": int(k),
        "symmetrize": symmetrize,
        "directed_neighbor_count": n * int(k),
        "undirected_candidate_count": len(candidates),
        "edge_count": len(edges),
        "nonpositive_candidate_edges_removed": len(candidates) - len(edges),
        "connected_component_count": len({find(i) for i in range(n)}),
        "isolated_vertex_count": len(isolated),
        "minimum_degree": int(degrees.min()),
        "maximum_degree": int(degrees.max()),
        "mean_degree": float(degrees.mean()),
        "density": float(2 * len(edges) / (n * (n - 1))),
        "total_edge_weight": float(sum(weights)),
    }
    return {"n": n, "neighbors": neighbors, "edges": edges, "weights": weights,
            "isolated_vertices": isolated, "stats": stats}


def _canonical_membership(membership):
    remapping = {}
    return np.asarray([remapping.setdefault(int(label), len(remapping)) for label in membership],
                      dtype=np.int64)


def run_leiden(graph, resolution, seed, n_iterations=-1):
    """Return full-cohort memberships and RB configuration quality.

    Zero-edge graphs have a deterministic singleton partition and undefined
    modularity quality. Isolates are also singleton in nonempty graphs, avoiding
    arbitrary grouping of vertices that have no positive-similarity edges.
    """
    if not math.isfinite(resolution) or resolution <= 0:
        raise ValueError("resolution must be finite and positive")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or not 0 <= seed <= 2**31 - 1:
        raise ValueError("seed must be an integer between 0 and 2**31 - 1")
    if (isinstance(n_iterations, bool) or not isinstance(n_iterations, int)
            or n_iterations == 0 or n_iterations < -1):
        raise ValueError("n_iterations must be -1 (until convergence) or a positive integer")
    if not graph["edges"]:
        return np.arange(graph["n"], dtype=np.int64), {
            "quality": None, "quality_reason": "no_positive_edges",
            "isolates_as_singletons": len(graph["isolated_vertices"]),
        }
    try:
        import igraph as ig
        import leidenalg
    except ImportError as exc:
        raise RuntimeError("Leiden requires igraph and leidenalg; install requirements-clustering.txt") from exc
    network = ig.Graph(n=graph["n"], edges=graph["edges"], directed=False)
    network.es["weight"] = graph["weights"]
    partition = leidenalg.find_partition(
        network, leidenalg.RBConfigurationVertexPartition, weights="weight",
        resolution_parameter=float(resolution), seed=int(seed), n_iterations=n_iterations,
    )
    membership = list(partition.membership)
    next_label = max(membership, default=-1) + 1
    for index in graph["isolated_vertices"]:
        membership[index] = next_label
        next_label += 1
    membership = _canonical_membership(membership)
    # The same objective, now explicitly evaluated after the isolate policy.
    quality = float(leidenalg.RBConfigurationVertexPartition(
        network, initial_membership=membership.tolist(), weights="weight",
        resolution_parameter=float(resolution),
    ).quality())
    return membership, {"quality": quality, "quality_reason": "",
                        "isolates_as_singletons": len(graph["isolated_vertices"])}


def seed_stability(memberships):
    """Return pairwise ARIs without using curated labels or selecting a seed."""
    from sklearn.metrics import adjusted_rand_score

    pairs = [{"seed_a": int(a), "seed_b": int(b),
              "pairwise_ari": float(adjusted_rand_score(memberships[a], memberships[b]))}
             for a, b in combinations(sorted(memberships), 2)]
    scores = [row["pairwise_ari"] for row in pairs]
    summary = {"seed_count": len(memberships), "seed_pair_count": len(pairs),
               "pairwise_ari_mean": float(np.mean(scores)) if scores else None,
               "pairwise_ari_std": float(np.std(scores)) if scores else None,
               "pairwise_ari_min": float(min(scores)) if scores else None,
               "stability_reason": "" if scores else "requires_at_least_two_seeds"}
    return pairs, summary


def _positive_unique(values, name):
    if len(values) != len(set(values)) or any(not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError(f"{name} must contain distinct positive finite values")


def run(args):
    import clustering_common as common

    _positive_unique(args.neighbors, "neighbors")
    _positive_unique(args.resolutions, "resolutions")
    if len(set(args.seeds)) != len(args.seeds) or any(not 0 <= seed <= 2**31 - 1 for seed in args.seeds):
        raise ValueError("seeds must be distinct integers between 0 and 2**31 - 1")
    if args.iterations == 0 or args.iterations < -1:
        raise ValueError("iterations must be -1 or a positive integer")
    try:
        import igraph  # noqa: F401 -- fail before creating output if dependencies are missing
        import leidenalg  # noqa: F401
    except ImportError as exc:
        raise RuntimeError("Install requirements-clustering.txt before running Leiden") from exc
    bundles = common.load_inputs(args)
    for bundle in bundles:
        if max(args.neighbors) >= len(bundle["ids"]):
            raise ValueError(f"Every k must be below {bundle['dataset']} sample count {len(bundle['ids'])}")
    output = common.create_output(args.output)
    metrics_rows, stability_rows, stability_summary = [], [], []
    for bundle in bundles:
        target = output / bundle["dataset"] / bundle["route"]
        target.mkdir(parents=True)
        assignments = []
        for k in args.neighbors:
            graph = build_knn_graph(bundle["distances"], k, args.symmetrize)
            graph_dir = target / "graphs" / f"k_{k}"
            graph_dir.mkdir(parents=True)
            common.write_json(graph_dir / "graph.json", graph["stats"])
            common.write_tsv(graph_dir / "edges.tsv", [
                {"source_row": i, "target_row": j, "source_id": bundle["ids"][i],
                 "target_id": bundle["ids"][j], "cosine_distance": float(bundle["distances"][i, j]),
                 "weight": weight}
                for (i, j), weight in zip(graph["edges"], graph["weights"])
            ], fields=("source_row", "target_row", "source_id", "target_id", "cosine_distance", "weight"))
            common.write_tsv(graph_dir / "neighbors.tsv", [
                {"source_row": i, "neighbor_rank": rank + 1, "target_row": int(j),
                 "source_id": bundle["ids"][i], "target_id": bundle["ids"][int(j)],
                 "cosine_distance": float(bundle["distances"][i, j])}
                for i in range(graph["n"]) for rank, j in enumerate(graph["neighbors"][i])
            ])
            for resolution_index, resolution in enumerate(args.resolutions):
                memberships, group_metrics = {}, []
                for seed in args.seeds:
                    membership, quality = run_leiden(graph, resolution, seed, args.iterations)
                    run_id = f"k{k}_r{resolution_index:03d}_s{seed}"
                    row = {"dataset": bundle["dataset"], "route": bundle["route"], "run_id": run_id,
                           "k": k, "resolution": resolution, "seed": seed,
                           "symmetrize": args.symmetrize, "n_iterations": args.iterations,
                           **graph["stats"], "rb_configuration_quality": quality["quality"],
                           "quality_reason": quality["quality_reason"],
                           **common.score_partition(bundle["labels"], membership, bundle["distances"])}
                    metrics_rows.append(row)
                    group_metrics.append(row)
                    assignments.append((run_id, membership))
                    memberships[seed] = membership
                pairs, stability = seed_stability(memberships)
                prefix = {"dataset": bundle["dataset"], "route": bundle["route"],
                          "k": k, "resolution": resolution, "symmetrize": args.symmetrize}
                stability_rows.extend({**prefix, **pair} for pair in pairs)
                summary = {**prefix, **stability}
                for metric in ("n_clusters", "v_measure", "homogeneity", "completeness", "ari", "ami",
                               "n_singleton_clusters", "silhouette_cosine"):
                    values = [row[metric] for row in group_metrics if row.get(metric) is not None]
                    summary[metric + "_mean"] = float(np.mean(values)) if values else None
                    summary[metric + "_std"] = float(np.std(values)) if values else None
                stability_summary.append(summary)
        common.write_assignments(target / "clusters.tsv", bundle, assignments)
        print(f"Leiden {bundle['dataset']}/{bundle['route']}: {len(assignments)} configurations, "
              f"{len(bundle['ids'])} samples retained", flush=True)
    common.write_tsv(output / "stability.tsv", stability_rows,
                     fields=("dataset", "route", "k", "resolution", "symmetrize", "seed_a", "seed_b", "pairwise_ari"))
    common.write_tsv(output / "seed_summary.tsv", stability_summary)
    parameters = {
        "neighbors": args.neighbors, "resolutions": args.resolutions, "seeds": args.seeds,
        "symmetrize": args.symmetrize, "n_iterations": args.iterations,
        "metric": "cosine", "clustering_space": "original normalized embedding space",
        "neighbor_search": "exact; self excluded; ties by sorted sample ID",
        "edge_weight": "max(0, cosine_similarity); zero-weight edges omitted; reciprocal edges counted once",
        "objective": "RBConfigurationVertexPartition (weighted resolution modularity)",
        "isolate_policy": "retain every isolated vertex as its own singleton",
        "parameter_selection": "none; export all requested configurations and seeds",
        "stability": "pairwise adjusted Rand index across seeds at fixed k and resolution; population std",
        "reference": "https://leidenalg.readthedocs.io/en/stable/reference.html",
    }
    common.finish_run(output, args, bundles, "cosine_knn_leiden", parameters, metrics_rows)
    return output


def build_parser():
    import clustering_common as common

    parser = argparse.ArgumentParser(description=__doc__)
    common.add_common_arguments(parser, default_output="data/processed/clustering_leiden_v1")
    parser.add_argument("--neighbors", nargs="+", type=int, default=[5, 10, 15], help="Exact kNN sizes; each must be below N")
    parser.add_argument("--resolutions", nargs="+", type=float, default=[0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
    parser.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    parser.add_argument("--symmetrize", choices=("union", "mutual"), default="union")
    parser.add_argument("--iterations", type=int, default=-1, help="-1 runs until no improvement; otherwise positive")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run(args)
    except (ValueError, RuntimeError, OSError, ImportError, KeyError) as exc:
        parser.exit(1, f"Leiden clustering failed: {exc}\n")


if __name__ == "__main__":
    main()
