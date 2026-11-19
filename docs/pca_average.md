# PCA 后的 cosine＋average-linkage 比较

入口为 [`scripts/cluster_pca_average.py`](../scripts/cluster_pca_average.py)，
数值变换辅助函数为 [`scripts/pca_transform.py`](../scripts/pca_transform.py)。
这次实验比较 PCA 对 B/E × ESM2、BGC-MAC、DNA full 的影响，只运行 average-linkage。
使用现有服务器特征，不需要重新调用模型或使用 GPU。

默认实验已于 2026-11-12 完成：[汇总报告](../data/processed/clustering_pca_average_v1/summary.md)、
[各条件对照](../data/processed/clustering_pca_average_v1/matched_metrics.tsv)、
[V-measure 维数曲线](../data/processed/clustering_pca_average_v1/v_measure_by_dimension.png)。
共完成 72 个条件、12,488 个距离切点的评分；9 个超出 B 有效秩的维数请求被明确跳过，
没有条件失败。全部划分保留 B=39、E=281 个样本。

六组 raw 对照均通过已有 average 结果的核验：距离矩阵逐位一致，指标在 `1e-12` 内一致。
另从成员表独立重算了全部 72 个选定切点的 H/C/V/ARI，结果一致；
全部有效主成分与仅中心化的 cosine 距离最大误差低于 `3e-15`。
27 项相关测试通过，覆盖 PCA 数值、实验流程、共享输入与 average 层次树。

## 运行

在项目根目录、安装 `requirements-clustering.txt` 的 Python 环境中执行：

```bash
python scripts/cluster_pca_average.py
```

默认读取 `data/processed/features_be_v1/` 和
`data/processed/pilot_be_v1/cohort.tsv`，写入新的
`data/processed/clustering_pca_average_v1/`。本地默认结果已存在，重跑时应按下面的命令
指定新的 `--output`；已有目录不会覆盖。
默认固定维数为 `2 4 8 16 32 64 128 256`，方差保留率为 `0.90 0.95 0.99`，
同时加入原始向量、仅中心化和全部有效主成分的对照。

指定输入位置、维数和输出目录：

```bash
python scripts/cluster_pca_average.py \
  --features-root data/processed/features_be_v1 \
  --cohort data/processed/pilot_be_v1/cohort.tsv \
  --datasets B E \
  --routes esm2 mac dna_full \
  --dimensions 2 4 8 16 32 64 128 256 \
  --variance-targets 0.90 0.95 0.99 \
  --reference-average data/processed/clustering_average_v1 \
  --output data/processed/clustering_pca_average_v2
```

`--no-plots` 关闭图形输出。原始向量对照仍会在本次实验中计算；
`--reference-average` 指向已有 average 结果以核对原始基线，默认使用
`data/processed/clustering_average_v1/`。参考目录缺失时明确记录未核对，
不会用其他输入或另一种聚类方法代替。

## 变换与对照

每个数据集、每条表征路线单独拟合 PCA，B 和 E 不合并拟合。
流程为：已验证的 L2 特征 → float64 数值计算 → 减去该组均值 → full SVD PCA
→ 选定主成分 → 每行 L2 归一化 → cosine 距离 → average-linkage（UPGMA）。
不做逐特征标准化或 whitening。

| 条件 | 输入到 cosine 的向量 | 用途 |
| --- | --- | --- |
| `raw` | 原始 L2 特征 | 重现已有原始空间基线 |
| `centered_only` | 原始特征减去该组均值后重新 L2 | 判断中心化本身的影响 |
| `pca_*` | 前若干主成分坐标，重新 L2 | 判断截断维度后的效果 |
| `full_rank` 选择器对应的 PCA 条件 | 所有数值上有效的主成分，重新 L2 | 应保持仅中心化向量的 cosine 几何，作为实现对照 |

PCA 的中心化会改变 cosine。因而“PCA 比 raw 好”不自动说明丢弃维度有益；
同时比较相对 `raw` 和相对 `centered_only` 的变化，才能区分这两种影响。
全部有效主成分的条件可能降低数组列数，但不会有意丢弃该组样本中的有效变化。

中心化后的理论最大秩为 `min(N−1, D)`；实际使用 SVD 的数值秩。
B 的 39 个样本至多有 38 个有效主成分；E 至多有 280 个，且受输入维数约束
（例如 MAC 的 256 维）。超过有效秩的固定维数请求会记录为跳过，不悄悄改成较小维数。
方差选择器取达到目标累计解释方差的最少主成分。多个选择器指向同一维数时，
只计算一个条件并保留所有选择器，避免把相同设置当作重复证据。

解释方差和保留率均在投影后重新 L2 **之前**计算，分母是中心化数据的总方差。
它们描述向量变化的保留程度，不是 GCF 信息保留率，也不保证保留率越高分数越好。

## 样本与评分口径

继续使用显式行号与稳定 sample_id 对齐全部特征，保留原始 Unicode GCF 标签。
成功条件必须包含该数据集的全部样本。若中心化或投影使某行变为零/近零向量，
cosine 无法可靠计算，该条件会明确失败；不会删除该样本后继续评分。
原始数据、已有特征和已有聚类结果保持不变。

每个有效条件建立完整 average 层次树，扫描零阈值与所有不同的合并高度。
汇总取预测簇数最接近真实 GCF 数的切点，平局取较小阈值，保留实际簇数及差值。
这是沿用原基线的家族数辅助诊断，不按 V-measure 最大值选切点。
同高度合并一起发生，因此不保证每个目标簇数都可达到。

所有维数都报告 V-measure、ARI、homogeneity、completeness、覆盖和簇数。
95% 方差保留率只是预先固定的简表展示条件，不从这些标签上挑选最佳维数。
维数曲线用于探索；若据此选出了设置，需要在独立开发/评测安排下验证。

PCA 在当前整组无标签样本上拟合，属于对整组样本共同处理的 transductive 分析。
家族标签不参与 PCA 或层次树构建，但参与评分和汇总切点选择；这不是独立测试集泛化成绩。
若以后评估新样本，必须使用训练/开发集拟合并保存的固定变换，不能分别拟合后把坐标当作相同空间。
MAC 在 MIBiG 的训练重叠等既有解释限制仍然适用。

## 输出

首先查看输出根目录中的 `summary.md` 和 `matched_metrics.tsv`。

| 文件 | 内容 |
| --- | --- |
| `summary.md` | 原始、仅中心化、95% 方差条件简表，以及各设置结果和解释限制 |
| `matched_metrics.tsv` | 每个条件按最近真实家族数选出的切点、指标及相对对照的差值 |
| `metrics.tsv` | 所有有效条件、所有距离切点的完整评分 |
| `conditions.tsv` | 条件、选择器、维数、方差保留率及状态 |
| `skipped.tsv` | 超出有效秩等未执行请求及原因 |
| `metadata.json` | 参数、输入与脚本来源、验证信息、输出校验值 |
| `v_measure_by_dimension.png/.svg` | 各维数 V-measure 与原始/仅中心化对照 |
| `ari_by_dimension.png/.svg` | 各维数 ARI 与原始/仅中心化对照 |
| `Dataset_*/route/transform.npz` | 均值、主成分、奇异值、解释方差等可追溯变换参数 |
| `Dataset_*/route/` 的样本映射 | 原始特征行号、稳定 ID、原始标签与来源 |
| 各条件目录的 `unit_vectors.npy`、`cosine_distances.npy`、`linkage.npy` | 实际聚类使用的向量、距离和树 |
| 各条件目录的 `clusters.tsv`、`metrics.tsv` | 每个切点的完整成员和指标 |

投影结果是本实验的分析产物，不能仅凭文件名当作原始模型导出的 `bgc_vectors_v1` 特征。
读取时应使用本次保存的行号/ID 映射和变换元数据。
