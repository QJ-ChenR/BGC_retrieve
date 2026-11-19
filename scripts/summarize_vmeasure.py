#!/usr/bin/env python3
"""Recompute and summarize existing B/E clusterings; do not rerun clustering.

The main comparison is a label-informed family-count diagnostic, not a held-out
parameter-selection score. Every Leiden seed is retained; k is not optimized.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import statistics

import numpy as np
from sklearn.metrics import adjusted_rand_score, homogeneity_completeness_v_measure

from clustering_common import file_sha256, read_json, read_tsv, write_json, write_tsv

ROOT = Path(__file__).resolve().parents[1]
PAPER = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12936217/"
DISPLAY = {"esm2": "ESM2", "mac": "BGC-MAC", "dna_full": "DNA full"}


def recompute(directory, cohort):
    metadata = read_json(directory / "metadata.json")
    if metadata.get("status") != "complete":
        raise ValueError(f"Incomplete clustering output: {directory}")
    for relative, expected in metadata["artifact_sha256"].items():
        if file_sha256(directory / relative) != expected:
            raise ValueError(f"Clustering artifact changed: {relative}")
    rows = read_tsv(directory / "metrics.tsv")
    groups = defaultdict(list)
    for row in rows:
        groups[row["dataset"], row["route"]].append(row)
    count = 0
    for (dataset, route), metrics in groups.items():
        expected = sorted((r for r in cohort if r["dataset"] == dataset), key=lambda r: r["sample_id"])
        expected_ids = [r["sample_id"] for r in expected]
        truth = [r["gcf_raw"] for r in expected]
        members = defaultdict(list)
        for member in read_tsv(directory / dataset / route / "clusters.tsv"):
            members[member["run_id"]].append(member)
        if set(members) != {r["run_id"] for r in metrics} or len(metrics) != len(members):
            raise ValueError("Metrics and cluster assignment runs disagree")
        for metric in metrics:
            assignment = sorted(members[metric["run_id"]], key=lambda r: int(r["row"]))
            if ([r["sample_id"] for r in assignment] != expected_ids
                    or [int(r["row"]) for r in assignment] != list(range(len(expected)))
                    or [r["gcf_raw"] for r in assignment] != truth):
                raise ValueError("Cluster assignments do not cover the unchanged cohort")
            predicted = [int(r["cluster_id"]) for r in assignment]
            h, c, v = homogeneity_completeness_v_measure(truth, predicted)
            ari = adjusted_rand_score(truth, predicted)
            for key, value in (("homogeneity", h), ("completeness", c), ("v_measure", v), ("ari", ari)):
                if abs(float(metric[key]) - value) > 1e-12:
                    raise ValueError(f"Recomputed {key} disagrees: {dataset}/{route}/{metric['run_id']}")
                metric[key] = float(value)
            for key, value in (("sample_count", len(expected)), ("n_clusters", len(set(predicted))),
                               ("true_gcf_count", len(set(truth)))):
                if int(metric[key]) != value:
                    raise ValueError(f"Count mismatch: {key}")
                metric[key] = value
            if float(metric["coverage"]) != 1:
                raise ValueError("Expected full-cohort coverage")
            count += 1
    return rows, {"path": str(directory), "metadata_sha256": file_sha256(directory / "metadata.json"),
                  "recomputed_partitions": count, "maximum_allowed_metric_error": 1e-12}


def aggregate_leiden(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row["dataset"], row["route"], int(row["k"]), float(row["resolution"])].append(row)
    result = []
    for (dataset, route, k, resolution), group in sorted(groups.items()):
        seeds = [int(r["seed"]) for r in group]
        if len(set(seeds)) != len(seeds):
            raise ValueError("Duplicate Leiden seed in one configuration")
        target = group[0]["true_gcf_count"]
        summary = {"dataset": dataset, "route": route, "method": "leiden", "k": k,
                   "resolution": resolution, "seed_count": len(group),
                   "seeds": ",".join(map(str, sorted(seeds))), "sample_count": group[0]["sample_count"],
                   "true_gcf_count": target,
                   "n_clusters_min": min(r["n_clusters"] for r in group),
                   "n_clusters_max": max(r["n_clusters"] for r in group),
                   "mean_absolute_family_count_error": statistics.mean(abs(r["n_clusters"] - target) for r in group)}
        for metric in ("v_measure", "ari", "homogeneity", "completeness", "n_clusters"):
            values = [r[metric] for r in group]
            summary[metric + "_mean"] = statistics.mean(values)
            summary[metric + "_std"] = statistics.pstdev(values)
        result.append(summary)
    return result


def select_summaries(average, leiden):
    average_groups, leiden_groups = defaultdict(list), defaultdict(list)
    for row in average:
        average_groups[row["dataset"], row["route"]].append(row)
    for row in leiden:
        leiden_groups[row["dataset"], row["route"], row["k"]].append(row)
    average_selected, leiden_selected, oracle = [], [], []
    for (dataset, route), group in sorted(average_groups.items()):
        chosen = min(group, key=lambda r: (abs(r["n_clusters"] - r["true_gcf_count"]), float(r["threshold"])))
        keep = {key: chosen[key] for key in ("dataset", "route", "run_id", "sample_count", "true_gcf_count",
                                             "n_clusters", "threshold", "v_measure", "ari", "homogeneity", "completeness")}
        keep.update(method="average", family_count_error=abs(chosen["n_clusters"] - chosen["true_gcf_count"]),
                    selection="closest true GCF count; tie: smallest threshold; label-informed")
        average_selected.append(keep)
        best = min(group, key=lambda r: (-r["v_measure"], float(r["threshold"])))
        oracle.append({"dataset": dataset, "route": route, "method": "average",
                       "v_measure_mean": best["v_measure"], "v_measure_std": 0,
                       "n_clusters_mean": best["n_clusters"], "threshold": best["threshold"],
                       "ari_mean": best["ari"], "selection": "max same-label V; optimistic diagnostic ONLY"})
    for (_, _, _), group in sorted(leiden_groups.items()):
        chosen = min(group, key=lambda r: (r["mean_absolute_family_count_error"], r["resolution"]))
        low = min(r["n_clusters_min"] for r in group)
        high = max(r["n_clusters_max"] for r in group)
        leiden_selected.append({**chosen,
                                "target_bracketed_in_grid": low <= chosen["true_gcf_count"] <= high,
                                "grid_n_clusters_min": low, "grid_n_clusters_max": high,
                                "selection": "fixed k; min mean(abs(K-GCF)) over seeds; tie: smaller resolution; label-informed"})
    groups = defaultdict(list)
    for row in leiden:
        groups[row["dataset"], row["route"]].append(row)
    for (dataset, route), group in sorted(groups.items()):
        best = min(group, key=lambda r: (-r["v_measure_mean"], r["k"], r["resolution"]))
        oracle.append({"dataset": dataset, "route": route, "method": "leiden",
                       **{key: best[key] for key in ("v_measure_mean", "v_measure_std", "n_clusters_mean",
                                                     "ari_mean", "k", "resolution", "seed_count")},
                       "selection": "max same-label mean V across k/resolution; all seeds retained; optimistic diagnostic ONLY"})
    return average_selected, leiden_selected, oracle


def singleton_controls(cohort):
    rows = []
    for dataset in sorted({r["dataset"] for r in cohort}):
        truth = [r["gcf_raw"] for r in cohort if r["dataset"] == dataset]
        h, c, v = homogeneity_completeness_v_measure(truth, range(len(truth)))
        rows.append({"dataset": dataset, "sample_count": len(truth), "true_gcf_count": len(set(truth)),
                     "n_clusters": len(truth), "homogeneity": float(h), "completeness": float(c),
                     "v_measure": float(v), "ari": float(adjusted_rand_score(truth, range(len(truth))))})
    return rows


def report(averages, leiden, controls, published):
    lines = ["# B/E V-measure 比较", "",
             "B：39 个 BGC / 23 个 GCF；E：281 个 BGC / 90 个 GCF。所有本地划分均保留全部样本。",
             "已从成员表和原始 cohort 标签重算全部 1500 个划分的 V-measure、ARI，核对误差不超过 1e-12。", "",
             "## 选参口径", "",
             "Average 选择预测簇数最接近真实 GCF 数的切点，并列取较小 threshold；六组均精确匹配。",
             "Leiden 不跨 k 挑最高分：每个 k 分别选择五个种子的平均 |预测簇数−真实 GCF 数| 最小的 resolution，并列取较小 resolution。",
             "下方简表固定显示 k=10（原网格的中间值，未按 V-measure 选择）；完整 k=5/10/15 表在 leiden_by_k.tsv。",
             "均值±标准差是同一设置的五个随机种子的总体标准差，不是独立样本置信区间。", "",
             "这是使用真实家族数量的诊断比较，不是独立开发集选参或泛化成绩。BiG-SCAPE 论文也采用最近家族数的 cutoff 口径。",
             f"来源：[BiG-SCAPE 2.0 论文 Fig.2 与 Methods]({PAPER})。", "",
             "## 本地路线主表", "",
             "| 方法 | 表征 | B V-measure | B 预测簇数 | E V-measure | E 预测簇数 |",
             "| --- | --- | ---: | ---: | ---: | ---: |"]
    for method, entries in (("Average", averages), ("Leiden k=10", [r for r in leiden if r["k"] == 10])):
        for route in DISPLAY:
            b, e = [next(r for r in entries if r["route"] == route and r["dataset"].startswith(f"Dataset_{letter}_"))
                    for letter in ("B", "E")]
            if method == "Average":
                cells = [f"{b['v_measure']:.4f}", str(b["n_clusters"]), f"{e['v_measure']:.4f}", str(e["n_clusters"])]
            else:
                cells = [f"{b['v_measure_mean']:.4f} ± {b['v_measure_std']:.4f}", f"{b['n_clusters_min']}–{b['n_clusters_max']}",
                         f"{e['v_measure_mean']:.4f} ± {e['v_measure_std']:.4f}", f"{e['n_clusters_min']}–{e['n_clusters_max']}"]
            lines.append(f"| {method} | {DISPLAY[route]} | " + " | ".join(cells) + " |")
    lines += ["", "**E 的 Leiden 当前网格没有接近 90 簇：全网格最多 50 簇。** 因而本表中的 Leiden 是已有网格的最近点，",
              "不应据此断言 Leiden 算法弱于 Average 或 BiG-SCAPE；应先补足 resolution 覆盖，再按固定规则比较。",
              "B 也并非每个 k 都能接近 23 簇。完整表保留所有实际簇数及差值。", "", "## BiG-SCAPE 公开分数", ""]
    if published:
        lines += ["以下数值来自作者公开的 benchmark 原始输出，不是本机重新运行。每一模式单独列出，不取跨模式最高分当默认性能。",
                  "cutoff 和预测簇数按作者最近家族数规则提取。日志样本数为 B=39/E=281，但尚未逐条核验有效 ID。",
                  "部分公开混淆矩阵的样本/家族计数不完整（例如 E / 1.1.9 global 矩阵计数为 253，而日志为 281）；",
                  "因此保留作者 Summary 的分数，不用不完整矩阵重算后替换，也不将日志计数视为完整评分覆盖的证明。",
                  "详细来源及验证限制见 sources/provenance.json 和 bigscape_published.tsv。", "",
                  "| 版本 | alignment | extension | B V-measure | B 簇数 / cutoff | E V-measure | E 簇数 / cutoff |",
                  "| --- | --- | --- | ---: | --- | ---: | --- |"]
        groups = defaultdict(dict)
        for row in published:
            groups[row["version"], row["alignment"], row["extension"]][row["dataset"]] = row
        for key, group in sorted(groups.items()):
            values = []
            for letter in ("B", "E"):
                row = next((r for d, r in group.items() if d.startswith(f"Dataset_{letter}_")), None)
                values.extend([f"{float(row['v_measure']):.4f}",
                               f"{row.get('n_clusters') or '未知'} / {row.get('cutoff') or '未知'}"] if row else ["—", "—"])
            lines.append("| " + " | ".join([*key, *values]) + " |")
        urls = sorted({r["source_url"] for r in published if r.get("source_url")})
        lines += ["", "来源：" + "；".join(f"[作者公开 benchmark 输出]({url})" for url in urls) + "。",
                  "归档中的 summary 数值存在小数位取舍，精度与重算核对情况见来源明细；不把打印精度当作额外证据。"]
    else:
        lines.append("未导入可核实的 BiG-SCAPE 数值，不用图片估读值或其他数据集分数代替。")
    lines += ["", "## 为什么不把扫描最大 V 当作主结果", "",
              "V-measure 没有随机一致性校正。即使每个样本各自成簇，也能得到较高分：", "",
              "| 数据集 | 全 singleton V | ARI |", "| --- | ---: | ---: |"]
    for row in controls:
        lines.append(f"| {row['dataset']} | {row['v_measure']:.4f} | {row['ari']:.4f} |")
    lines += ["", "因此另存 oracle_diagnostic.tsv 中的扫描最大值只表示使用同一批标签选参后的乐观诊断，",
              "既不是独立泛化成绩，也不是算法性能的数学上界。主表同时保留家族数，详细 TSV 还保留 ARI、homogeneity 和 completeness。", "",
              "## 当前证据的解释", "",
              "在精确匹配家族数的 Average 结果中，两个数据集均为 BGC-MAC > ESM2 > DNA full。",
              "固定 k=10 的现有 Leiden 最近点均值也呈同一顺序，但不代表已完成分辨率匹配。",
              "E 来自 MIBiG，MAC 的训练重叠仍未排查；B/E 之间也存在序列重叠，不能将二者当作独立训练/测试划分。",
              "与论文分数的并列仅为外部参考；要声称严格同样本优劣，需核对对应 BiG-SCAPE 输出 ID 和标签后重算。", "",
              "## 文件与复现", "",
              "- [比较图](comparison.png) / [SVG](comparison.svg)",
              "- [Average 切点、簇数及指标](average_matched.tsv)",
              "- [Leiden 所有 k 的选定设置](leiden_by_k.tsv) / [所有设置的种子汇总](leiden_all_settings.tsv)",
              "- [全 singleton 对照](singleton_controls.tsv) / [同标签最大值诊断](oracle_diagnostic.tsv)",
              "- [运行与输出校验记录](metadata.json)", "",
              "在项目根目录、安装 requirements-clustering.txt 的环境中运行；已有输出时指定新的目录：", "",
              "```bash", "python scripts/summarize_vmeasure.py \\",
              "  --published data/processed/vmeasure_comparison_v1/sources/bigscape_published.tsv \\",
              "  --output data/processed/vmeasure_comparison_v2", "```", ""]
    return "\n".join(lines)


def plot_summary(output, averages, leiden, published):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    colors = ["#2563eb", "#e07a1f", "#13856c"]
    for axis, letter in zip(axes, ("B", "E")):
        for index, (route, color) in enumerate(zip(DISPLAY, colors)):
            a = next(r for r in averages if r["dataset"].startswith(f"Dataset_{letter}_") and r["route"] == route)
            l = next(r for r in leiden if r["dataset"].startswith(f"Dataset_{letter}_") and r["route"] == route and r["k"] == 10)
            axis.bar(index - .18, a["v_measure"], .34, color=color, label="Average: matched GCF count" if index == 0 else None)
            axis.bar(index + .18, l["v_measure_mean"], .34, color=color, alpha=.45,
                     yerr=l["v_measure_std"], capsize=3, label="Leiden: k=10, nearest existing setting" if index == 0 else None)
            axis.text(index - .18, a["v_measure"] + .018, f"K={a['n_clusters']}", ha="center", fontsize=8)
            axis.text(index + .18, l["v_measure_mean"] + .045, f"K={l['n_clusters_min']}–{l['n_clusters_max']}", ha="center", fontsize=8)
        scores = [float(r["v_measure"]) for r in published if r["dataset"].startswith(f"Dataset_{letter}_") and r["version"].startswith("2")]
        if scores:
            axis.axhspan(min(scores), max(scores), color="#777777", alpha=.16, label="Published BiG-SCAPE 2 modes (external)")
        axis.set_xticks(range(3), [DISPLAY[r] for r in DISPLAY])
        axis.set_ylim(0, 1.06)
        axis.set_title(f"Dataset {letter}: " + ("39 BGCs / 23 GCFs" if letter == "B" else "281 BGCs / 90 GCFs"))
        axis.set_ylabel("V-measure")
        axis.grid(axis="y", alpha=.2)
    axes[0].legend(loc="lower left", fontsize=8)
    fig.suptitle("Family-count-informed diagnostic; E Leiden grid does not reach 90 families", fontsize=11)
    fig.tight_layout()
    for suffix in ("png", "svg"):
        fig.savefig(output / f"comparison.{suffix}", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--average", type=Path, default=ROOT / "data/processed/clustering_average_v1")
    parser.add_argument("--leiden", type=Path, default=ROOT / "data/processed/clustering_leiden_v1")
    parser.add_argument("--cohort", type=Path, default=ROOT / "data/processed/pilot_be_v1/cohort.tsv")
    parser.add_argument("--output", type=Path, default=ROOT / "data/processed/vmeasure_comparison_v1")
    parser.add_argument("--published", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if (output / "summary.md").exists() or (output / "metadata.json").exists():
        parser.error("Summary already exists; choose a new --output directory")
    output.mkdir(parents=True, exist_ok=True)
    cohort = read_tsv(args.cohort)
    average, a_validation = recompute(args.average, cohort)
    leiden, l_validation = recompute(args.leiden, cohort)
    aggregated = aggregate_leiden(leiden)
    selected_a, selected_l, oracle = select_summaries(average, aggregated)
    controls = singleton_controls(cohort)
    published = read_tsv(args.published) if args.published else []
    write_tsv(output / "average_matched.tsv", selected_a)
    write_tsv(output / "leiden_by_k.tsv", selected_l)
    write_tsv(output / "leiden_all_settings.tsv", aggregated)
    write_tsv(output / "oracle_diagnostic.tsv", oracle)
    write_tsv(output / "singleton_controls.tsv", controls)
    if published:
        write_tsv(output / "bigscape_published.tsv", published)
    (output / "summary.md").write_text(report(selected_a, selected_l, controls, published), encoding="utf-8")
    plot_summary(output, selected_a, selected_l, published)
    write_json(output / "metadata.json", {
        "status": "complete", "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": file_sha256(__file__), "cohort_sha256": file_sha256(args.cohort),
        "validation": [a_validation, l_validation],
        "published_source": None if args.published is None else {"path": str(args.published.resolve()), "sha256": file_sha256(args.published)},
        "primary_leiden_k": 10, "seed_aggregation": "mean and population standard deviation; no seed selection",
        "selection_uses_true_family_count": True, "new_clustering_runs": False,
        "artifact_sha256": {p.name: file_sha256(p) for p in output.iterdir() if p.is_file()},
    })
    print(f"Recomputed {len(average) + len(leiden)} partitions; saved comparison to {output}")


if __name__ == "__main__":
    main()
