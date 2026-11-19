#!/usr/bin/env python3
"""Explore centered PCA + cosine + average linkage on complete native B/E exports.

PCA is fitted separately to each unlabeled dataset/route. All requested dimensions
are reported; true GCF counts select diagnostic tree cuts, never PCA dimensions.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

import clustering_common as common
from cluster_average import build_average_linkage, cut_average_linkage, threshold_grid
from pca_transform import fit_pca, normalized_projection, plan_dimensions


def choose_matched_cut(metrics):
    return min(metrics, key=lambda r: (abs(r["n_clusters"] - r["true_gcf_count"]), r["threshold"]))


def analyze_condition(directory, bundle, condition, unit, distances):
    """Reuse the baseline hierarchy, cut semantics and complete-cohort metrics."""
    directory.mkdir()
    np.save(directory / "unit_vectors.npy", unit, allow_pickle=False)
    np.save(directory / "cosine_distances.npy", distances, allow_pickle=False)
    tree = build_average_linkage(distances)
    np.save(directory / "linkage.npy", tree, allow_pickle=False)
    runs, metrics = [], []
    for index, threshold in enumerate(threshold_grid(tree)):
        run_id = f"cut_{index:04d}"
        membership = cut_average_linkage(tree, float(threshold))
        runs.append((run_id, membership))
        metrics.append({
            "dataset": bundle["dataset"], "route": bundle["route"],
            "condition": condition["condition"], "method": "cosine_average_linkage",
            "n_components": condition["n_components"],
            "run_id": run_id, "threshold": float(threshold),
            **common.score_partition(bundle["labels"], membership, distances),
        })
    common.write_assignments(directory / "clusters.tsv", bundle, runs)
    common.write_tsv(directory / "metrics.tsv", metrics)
    chosen = choose_matched_cut(metrics)
    common.write_json(directory / "metadata.json", {
        **condition, "status": "complete", "array_row_mapping": "../bgc_ids.tsv",
        "clustering_space": condition["condition"] + " L2 normalized cosine",
        "cut_selection": "closest true GCF count; tie: smaller threshold; diagnostic only",
        "selected_run_id": chosen["run_id"], "selected_threshold": chosen["threshold"],
        "sample_count": len(bundle["ids"]), "coverage": 1.0,
    })
    return metrics, {
        **chosen, "selectors": ",".join(condition["selectors"]),
        "retained_variance": condition["retained_variance"],
        "family_count_error": abs(chosen["n_clusters"] - chosen["true_gcf_count"]),
    }


def verify_raw_reference(path, bundle, metrics):
    if path is None or not Path(path).exists():
        return {"status": "unavailable", "path": str(path) if path else None}
    path = Path(path).resolve()
    meta = common.read_json(path / "metadata.json")
    if meta.get("status") != "complete" or meta.get("method") != "cosine_average_linkage":
        raise ValueError("Reference average-linkage run is not complete or has wrong method")
    relative = Path(bundle["dataset"]) / bundle["route"] / "metrics.tsv"
    for item in (Path("metrics.tsv"), relative, relative.parent / "samples.tsv",
                 relative.parent / "cosine_distances.npy"):
        digest = meta["artifact_sha256"].get(item.as_posix())
        if digest is None or common.file_sha256(path / item) != digest:
            raise ValueError(f"Reference artifact hash mismatch: {item}")
    source = next((r for r in meta["inputs"] if r["dataset"] == bundle["dataset"]
                   and r["route"] == bundle["route"]), None)
    if source is None or source["artifact_sha256"] != bundle["provenance"]["artifact_sha256"]:
        raise ValueError("Reference average run used different feature artifacts")
    samples = common.read_tsv(path / relative.parent / "samples.tsv")
    if ([r["sample_id"] for r in samples] != bundle["ids"]
            or [r["gcf_raw"] for r in samples] != bundle["labels"]):
        raise ValueError("Reference average run used different sample ordering or labels")
    np.testing.assert_array_equal(np.load(path / relative.parent / "cosine_distances.npy"), bundle["distances"])
    reference = common.read_tsv(path / relative)
    # Reference may use an explicit threshold grid. Compare the selected nearest-K
    # diagnostic only when it is the original all-height sweep.
    if meta["parameters"].get("requested_thresholds") is not None:
        raise ValueError("Raw reference must use the all-merge-height threshold sweep")
    if len(reference) != len(metrics):
        raise ValueError("Raw-reference cut count differs")
    for old, new in zip(reference, metrics):
        for key in ("threshold", "n_clusters", "v_measure", "ari", "homogeneity", "completeness"):
            if not np.isclose(float(old[key]), new[key], atol=1e-12, rtol=0):
                raise ValueError(f"Raw-reference {key} mismatch")
    return {"status": "verified", "path": str(path), "partitions": len(metrics),
            "metadata_sha256": common.file_sha256(path / "metadata.json"),
            "metric_tolerance": 1e-12, "cosine_distances_bitwise_equal": True}


def make_report(matched, skipped):
    lines = ["# PCA + average-linkage 比较", "",
             "每个数据集和模型分别对完整无标签队列拟合 PCA；所有条件均使用 cosine + average-linkage。",
             "PCA 只中心化，不按特征标准化、不 whitening；投影后逐行 L2 归一化。",
             "raw 为原始特征；centered_only 只减均值；full_rank 为保留全部数值有效主成分的对照。",
             "方差保留率相对于中心化后的总方差，计算于投影后 L2 归一化之前；它不是性能指标。", "",
             "每个条件扫描全部不同合并高度和零阈值，选择簇数最接近真实 GCF 数的切点，并列取较小阈值。",
             "维数按预设网格/方差目标决定，不按 V-measure 选最佳维数。这是传导式探索和真实家族数辅助诊断，不能当作独立测试成绩。",
             "PCA 保留全秩时应与 centered_only 的 cosine 距离相符；与 raw 的差别体现中心化的影响。", "",
             "## 预设 95% 方差条件及对照", "",
             "| 数据集 | 表征 | raw V | 仅中心化 V | PCA95 维数 | PCA95 V | Δraw | Δ中心化 | PCA95 ARI |",
             "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    keys = sorted({(r["dataset"], r["route"]) for r in matched})
    for dataset, route in keys:
        rows = [r for r in matched if (r["dataset"], r["route"]) == (dataset, route)]
        raw = next(r for r in rows if r["condition"] == "raw")
        centered = next((r for r in rows if r["condition"] == "centered_only"), None)
        pca = next((r for r in rows if "variance:0.95" in r["selectors"].split(",")), None)
        center_text = f"{centered['v_measure']:.4f}" if centered else "不可用"
        pca_cells = (f"{pca['n_components']} | {pca['v_measure']:.4f} | {pca['delta_v_raw']:+.4f} | "
                     f"{pca['delta_v_centered']:+.4f} | {pca['ari']:.4f}"
                     if pca and centered else "— | — | — | — | —")
        lines.append(f"| {dataset} | {route} | {raw['v_measure']:.4f} | {center_text} | {pca_cells} |")
    lines += ["", "95% 是运行前预设的展示点，不代表扫描中最佳维数；未请求或失败时显示空值。", "",
              "## 所有完成条件", "",
              "| 数据集 | 表征 | 条件 / 请求来源 | 维数 | 保留方差 | 簇数 / 真值 | V | ARI | Δraw | Δ中心化 |",
              "| --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |"]
    for r in matched:
        variance = "—" if r["retained_variance"] is None else f"{r['retained_variance']:.4f}"
        delta = "—" if r["delta_v_centered"] is None else f"{r['delta_v_centered']:+.4f}"
        lines.append(f"| {r['dataset']} | {r['route']} | {r['condition']} / {r['selectors']} | "
                     f"{r['n_components']} | {variance} | {r['n_clusters']} / {r['true_gcf_count']} | "
                     f"{r['v_measure']:.4f} | {r['ari']:.4f} | {r['delta_v_raw']:+.4f} | {delta} |")
    lines += ["", "## 覆盖与解释", "",
              f"共有 {len(matched)} 个完成条件；{len(skipped)} 个跳过/失败记录详见 [skipped.tsv](skipped.tsv)。",
              "超过数值有效秩的固定维数明确跳过，不自动裁剪后冒用请求维数。相同维数的多个请求合并计算，并保留所有 selectors。",
              "投影产生零/近零向量时，该条件整体失败并记录原因，不删除样本后评分。",
              "V-measure 无随机一致性校正，应同时看 ARI、簇数、homogeneity 和 completeness。",
              "MAC 的训练重叠尚未排查；本实验不改变编码器、不重算 Leiden、不更新 BiG-SCAPE 外部参考。", "",
              "[选定切点完整指标](matched_metrics.tsv) · [所有切点评分](metrics.tsv) · [条件清单](conditions.tsv) · [来源和校验记录](metadata.json)",
              "", "所有向量、距离矩阵和树的样本行号均对应同 dataset/route 下的 bgc_ids.tsv / samples.tsv。",
              "transform.npz 保存均值和主成分；unit_vectors.npy 是供本分析使用的 FP64 单位向量，不冒充原生编码器导出。", ""]
    return "\n".join(lines)


def plot_comparison(output, matched):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    datasets = sorted({r["dataset"] for r in matched})
    routes = list(dict.fromkeys(r["route"] for r in matched))
    for metric, title in (("v_measure", "V-measure"), ("ari", "ARI")):
        fig, axes = plt.subplots(len(datasets), len(routes), squeeze=False,
                                 figsize=(5 * len(routes), 3.7 * len(datasets)))
        for i, dataset in enumerate(datasets):
            for j, route in enumerate(routes):
                ax = axes[i, j]
                rows = [r for r in matched if (r["dataset"], r["route"]) == (dataset, route)]
                for name, color in (("raw", "#64748b"), ("centered_only", "#c56b15")):
                    control = next((r for r in rows if r["condition"] == name), None)
                    if control:
                        ax.axhline(control[metric], color=color, linestyle="--", label=name)
                projected = sorted((r for r in rows if r["condition"].startswith("pca_")),
                                   key=lambda r: r["n_components"])
                if projected:
                    ax.plot([r["n_components"] for r in projected], [r[metric] for r in projected],
                            "o-", color="#2563eb", label="PCA (all requested dimensions)")
                    ax.set_xscale("log", base=2)
                ax.set_title(f"{dataset.split('_')[1]} / {route}")
                ax.set_xlabel("Retained components (log2 scale)")
                ax.set_ylabel(title)
                ax.grid(alpha=.2)
                ax.legend(fontsize=8)
        fig.suptitle("Cosine + average linkage; cuts matched to true family count; dimensions not selected by score")
        fig.tight_layout()
        for suffix in ("png", "svg"):
            fig.savefig(output / f"{metric}_by_dimension.{suffix}", dpi=180)
        plt.close(fig)


def run(args):
    # Validate selectors before reading all inputs or creating any result directory.
    if not args.dimensions or any(d < 1 for d in args.dimensions):
        raise ValueError("Dimensions must be positive integers")
    if not args.variance_targets or any(not np.isfinite(t) or not 0 < t <= 1 for t in args.variance_targets):
        raise ValueError("Variance targets must be finite and in (0, 1]")
    bundles = common.load_inputs(args)
    if any(len(b["ids"]) < 2 for b in bundles):
        raise ValueError("PCA and average linkage need at least two samples")
    output = common.create_output(args.output)
    all_metrics, matched, conditions, skipped, fit_records, references = [], [], [], [], [], []
    for bundle in bundles:
        dataset, route = bundle["dataset"], bundle["route"]
        directory = output / dataset / route
        directory.mkdir(parents=True)
        common.write_tsv(directory / "bgc_ids.tsv", [
            {"row": i, "bgc_id": key, "source_feature_row": bundle["source_row_indices"][i]}
            for i, key in enumerate(bundle["ids"])])
        common.write_tsv(directory / "samples.tsv", [
            {"row": i, "source_feature_row": bundle["source_row_indices"][i], **row}
            for i, row in enumerate(bundle["rows"])])
        raw = {"condition": "raw", "n_components": bundle["vectors"].shape[1],
               "selectors": ["raw"], "retained_variance": None}
        metrics, selected = analyze_condition(directory / "raw", bundle, raw, bundle["vectors"], bundle["distances"])
        all_metrics.extend(metrics)
        local_matched = [selected]
        conditions.append({"dataset": dataset, "route": route, **raw, "selectors": "raw", "status": "complete"})
        references.append({"dataset": dataset, "route": route,
                           **verify_raw_reference(args.reference_average, bundle, metrics)})
        try:
            fit = fit_pca(bundle["vectors"])
        except ValueError as exc:
            skipped.append({"dataset": dataset, "route": route, "selector": "all_centered_pca",
                            "status": "failed", "reason": str(exc)})
            fit = None
        if fit is not None:
            np.savez(directory / "transform.npz", **{key: fit[key] for key in
                     ("mean", "components", "singular_values", "explained_variance", "explained_variance_ratio")})
            fit_records.append({"dataset": dataset, "route": route, **{k: v for k, v in fit.items()
                                if k not in ("mean", "components", "singular_values", "explained_variance",
                                             "explained_variance_ratio", "scores")}})
            common.write_tsv(directory / "explained_variance.tsv", [
                {"component": i + 1, "singular_value": float(fit["singular_values"][i]),
                 "explained_variance": float(fit["explained_variance"][i]),
                 "explained_variance_ratio": float(ratio),
                 "cumulative_ratio": float(np.sum(fit["explained_variance_ratio"][:i + 1]))}
                for i, ratio in enumerate(fit["explained_variance_ratio"])])
            plans, omitted = plan_dimensions(fit, args.dimensions, args.variance_targets)
            skipped.extend({"dataset": dataset, "route": route, "status": "skipped", **r} for r in omitted)
            plans.insert(0, {"condition": "centered_only", "n_components": bundle["vectors"].shape[1],
                             "selectors": ["centered_only"], "retained_variance": 1.0})
            centered_distances = None
            for plan in plans:
                scores = (bundle["vectors"] - fit["mean"] if plan["condition"] == "centered_only"
                          else fit["scores"][:, :plan["n_components"]])
                try:
                    unit, distances = normalized_projection(scores)
                except ValueError as exc:
                    record = {"dataset": dataset, "route": route, **plan,
                              "selectors": ",".join(plan["selectors"]), "status": "failed", "reason": str(exc)}
                    conditions.append(record)
                    skipped.append({**record, "selector": record["selectors"]})
                    continue
                if plan["condition"] == "centered_only":
                    centered_distances = distances
                elif "full_rank" in plan["selectors"] and centered_distances is not None:
                    np.testing.assert_allclose(distances, centered_distances, atol=1e-10, rtol=0)
                metrics, selected = analyze_condition(directory / plan["condition"], bundle, plan, unit, distances)
                all_metrics.extend(metrics)
                local_matched.append(selected)
                conditions.append({"dataset": dataset, "route": route, **plan,
                                   "selectors": ",".join(plan["selectors"]), "status": "complete"})
            print(f"{dataset}/{route}: PCA rank={fit['numerical_rank']}, {len(local_matched)} completed conditions", flush=True)
        baseline = local_matched[0]
        centered = next((r for r in local_matched if r["condition"] == "centered_only"), None)
        for row in local_matched:
            row["delta_v_raw"] = row["v_measure"] - baseline["v_measure"]
            row["delta_ari_raw"] = row["ari"] - baseline["ari"]
            row["delta_v_centered"] = row["v_measure"] - centered["v_measure"] if centered else None
            row["delta_ari_centered"] = row["ari"] - centered["ari"] if centered else None
        matched.extend(local_matched)
    common.write_tsv(output / "metrics.tsv", all_metrics)
    common.write_tsv(output / "matched_metrics.tsv", matched)
    common.write_tsv(output / "conditions.tsv", conditions)
    skip_fields = list(dict.fromkeys(["dataset", "route", "selector", "status", "reason"] +
                                    [key for row in skipped for key in row]))
    common.write_tsv(output / "skipped.tsv", skipped, skip_fields)
    (output / "summary.md").write_text(make_report(matched, skipped), encoding="utf-8")
    if not args.no_plots:
        plot_comparison(output, matched)
    provenance = common.source_provenance()
    for name in ("cluster_pca_average.py", "pca_transform.py"):
        provenance["source_sha256"][name] = common.file_sha256(Path(__file__).with_name(name))
    status = "complete_with_condition_failures" if any(r["status"] == "failed" for r in skipped) else "complete"
    common.write_json(output / "metadata.json", {
        "format": "bgc_pca_average_v1", "status": status,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "cohort": {"path": str(args.cohort.resolve()), "sha256": common.file_sha256(args.cohort)},
        "inputs": [b["provenance"] for b in bundles], "pca_fits": fit_records,
        "raw_reference_checks": references, "completed_conditions": len(matched),
        "scored_partitions": len(all_metrics), "skipped_or_failed_requests": len(skipped),
        "parameters": {
            "fit_scope": "separate full unlabeled dataset/route; transductive",
            "input": "original FP32 exported unit vectors renormalized in float64, as raw baseline",
            "pca": "float64 centered full SVD; no whitening or feature standardization; numerical rank capped at min(N-1,D)",
            "postprojection": "row L2 normalization, then cosine distances",
            "distance_cut_selection": "closest true GCF count, tie lower threshold; all merge heights plus zero",
            "dimension_selection": "declared dimensions and variance targets; full rank control; no V optimization",
            "rank_exceeded": "skip request explicitly; never silently clip",
            "zero_vector_policy": "entire condition fails explicitly; never drop samples",
            "leiden_runs": 0,
        },
        "artifact_sha256": {p.relative_to(output).as_posix(): common.file_sha256(p)
                            for p in sorted(output.rglob("*")) if p.is_file() and p.name != "status.json"
                            and not any(part.startswith(".") for part in p.relative_to(output).parts)},
        **provenance,
    })
    (output / "status.json").write_text(json.dumps({"status": status}, indent=2) + "\n", encoding="utf-8")
    return output


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    common.add_common_arguments(parser, default_output="data/processed/clustering_pca_average_v1")
    parser.add_argument("--dimensions", type=int, nargs="+", default=[2, 4, 8, 16, 32, 64, 128, 256])
    parser.add_argument("--variance-targets", type=float, nargs="+", default=[.9, .95, .99])
    parser.add_argument("--reference-average", type=Path,
                        default=common.ROOT / "data/processed/clustering_average_v1",
                        help="Verify raw baseline if this prior run exists; absent paths are recorded")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        output = run(args)
    except (ValueError, OSError, KeyError, RuntimeError, ImportError, AssertionError) as exc:
        parser.exit(1, f"PCA average-linkage analysis failed: {exc}\n")
    print(f"Saved PCA average-linkage analysis to {output}")


if __name__ == "__main__":
    main()
