# B/E 的两个原始 cosine 聚类基线

这两个脚本读取已经导出的 BGC 向量，在 CPU 上分别运行层次聚类和图社区聚类。
不需要重新加载 ESM2、BGC-MAC、DNA encoder，也不需要 GPU。
默认对 B、E 各自的 `esm2`、`mac`、`dna_full` 三条路线独立分析，保留 B 的
39 个样本和 E 的 281 个样本；两个数据集不合并成一个家族标签空间。

| 脚本 | 计算方式 | 主要用途 |
| --- | --- | --- |
| `scripts/cluster_average.py` | 原始 cosine 距离 → average-linkage → 距离阈值切树 | 查看全局家族层次及不同合并尺度 |
| `scripts/cluster_leiden.py` | 原始 cosine 距离 → 精确 kNN 图 → 加权 Leiden | 查看局部邻域社区及其参数、随机种子稳定性 |

两个脚本都不使用 PCA、UMAP 或二维投影；GCF 标签只在聚类完成后用于评分。
它们不会根据标签自动选择最佳阈值、最佳 resolution、最佳种子或指定真实家族数。

## 环境与输入

以下命令在项目 `BGC_retrieve` 根目录运行。建议新建 Python 3.12 环境，
`requirements-clustering.txt` 固定了本次使用的直接依赖版本：

```bash
python3.12 -m venv .venv-clustering
.venv-clustering/bin/python -m pip install -r requirements-clustering.txt
```

实际特征默认位于 `data/processed/features_be_v1/`，不是 `data/feature/`：

```text
data/processed/features_be_v1/
  run.json
  Dataset_B_Divergent/
    esm2/
    mac/
    dna_full/
  Dataset_E_MIBiG_BS1/
    esm2/
    mac/
    dna_full/
```

每条路线读取原生导出的四个文件：

- `embeddings.npy`：有限、非零、L2 归一化的 FP32 `[样本数, 维度]` 数组。
- `bgc_ids.tsv`：用 `row` 对应矩阵行，用 `bgc_id` 对应稳定 `sample_id`。
- `metadata.json`：必须为 `complete`，路线、维度、输入数及每个样本的源文件哈希需一致。
- `failures.tsv`：要求没有失败或排除的样本行。

标签和原始名称默认来自 `data/processed/pilot_be_v1/cohort.tsv`。
加载器先按 TSV 中的明确行号取向量，再按稳定 `sample_id` 排序；不依赖 TSV
物理行顺序、目录遍历顺序或 GCF 标签顺序。标签保留原始 Unicode，包括 B 中
的末尾不换行空格（NBSP），不会通过 `.strip()` 改写标签。

加载时要求选定路线覆盖完整清单，不会自动取成功 ID 的交集。若特征根目录含有
`run.json`，还会核对提取阶段保存的四个产物哈希；同时选择 ESM2 和 MAC 时，
会核对两者元数据中的 ESM token 缓存哈希一致。聚类无需加载 `tokens.pth`。
元数据中的服务器绝对路径仅保留为来源信息，不要求这些旧路径在本机存在。

cosine 使用 float64 运算，并重新归一化已验证的单位向量以控制数值误差；这是
cosine 计算的一部分，不是另一个学习到的特征空间。结果保存实际使用的距离矩阵。

## Average-linkage

运行默认完整分析：

```bash
.venv-clustering/bin/python scripts/cluster_average.py
```

默认输出到 `data/processed/clustering_average_v1/`。使用 SciPy average linkage
（UPGMA），传入 cosine 距离的 condensed 数组，而不是把方形距离矩阵误当作
新的特征矩阵。平均连接的距离是两个簇之间所有样本对距离的平均值。
实现定义见 [SciPy linkage 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)。

默认扫描 **0 加上树中所有不同的合并高度**。当前 ESM2 的 cosine 距离尺度很小，
统一使用 0.1、0.2 等较大阈值容易只得到一个簇；按实际合并高度扫描能保留全部
可观察的层次变化。这是仅依赖输入向量的诊断扫描，并没有使用 GCF 标签选择切点。
等高合并会一起发生，因此不保证能切出任意指定簇数。切树规则为
`fcluster(..., criterion="distance")`，合并高度不超过阈值的分支归入同簇。
参见 [SciPy fcluster 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html)。

也可显式指定阈值、子集及新输出目录，例如：

```bash
.venv-clustering/bin/python scripts/cluster_average.py \
  --datasets B \
  --routes esm2 mac dna_full \
  --thresholds 0 0.001 0.003 0.01 0.03 0.1 0.3 0.6 \
  --output data/processed/clustering_average_custom_v1
```

该列表只是调用示例，不是推荐的最终评测阈值。阈值范围必须在 `[0, 2]`。

每个 `Dataset_*/<route>/` 中包括：

| 产物 | 含义 |
| --- | --- |
| `linkage.npy`、`linkage_rows.tsv` | 完整层次树，保留 SciPy 节点索引、合并距离与成员数 |
| `input_order.tsv` | 树中原始叶节点索引与输入样本的映射 |
| `leaf_order.tsv` | 树图显示的叶子顺序，与输入行顺序分开记录 |
| `clusters.tsv` | 所有切点的样本归属，通过 `run_id` 连接指标表 |
| `metrics.tsv` | 本路线所有阈值对应的聚类指标 |
| `dendrogram.svg/png` | 完整树图；B 标注原始 BGC 名称，E 隐藏密集叶标签，可查叶子映射表 |
| `threshold_curve.svg/png` | 指标、簇数和 singleton 数随阈值变化的曲线 |

加 `--no-plots` 可关闭绘图，树、成员表及数值指标仍会输出。

## kNN + Leiden

运行默认完整分析：

```bash
.venv-clustering/bin/python scripts/cluster_leiden.py
```

默认输出到 `data/processed/clustering_leiden_v1/`，参数网格为：

- `--neighbors 5 10 15`
- `--resolutions 0.25 0.5 1 2 4 8`
- `--seeds 0 1 2 3 4`
- `--symmetrize union`
- `--iterations -1`：迭代到一轮没有进一步改进。

每个数据集、每条路线共 90 个参数与种子组合，全部保留。图构建规则为：

1. 在原始 cosine 距离中精确查找 k 个邻居，明确排除自身；同距离按稳定样本
   ID 排序后的行号决定顺序。要求 `1 <= k < N`，不会自动缩小过大的 k。
2. 默认使用无向 union 图：任意方向选中邻居就保留该候选边。互相选中的邻居
   仍只形成一条边，不把其权重加倍。`--symmetrize mutual` 可改为双方都选中才保留。
3. 边权为正 cosine 相似度 `max(0, 1 - distance)`，删除权重为 0 的候选边并记录数量。
4. 使用 `RBConfigurationVertexPartition` 的加权 resolution modularity 目标。
   所有样本都保留；孤立点各自成为 singleton，全图没有正边时也输出全 singleton。

目标与接口定义见 [Leiden 官方文档](https://leidenalg.readthedocs.io/en/stable/reference.html)。
`rb_configuration_quality` 是未归一化的优化目标值，不能跨不同 k、不同图或
不同表示直接当作生物学优劣评分。

较小的显式参数扫描示例：

```bash
.venv-clustering/bin/python scripts/cluster_leiden.py \
  --datasets B E \
  --neighbors 5 10 \
  --resolutions 0.5 1 2 \
  --seeds 0 1 2 \
  --output data/processed/clustering_leiden_custom_v1
```

每条路线保存 `clusters.tsv`；`graphs/k_<k>/` 保存 `neighbors.tsv`、`edges.tsv`
和 `graph.json`，可核查邻居、权重、连通分量及孤立点。
根目录的 `metrics.tsv` 保存逐参数、逐种子结果；`seed_summary.tsv` 保存固定
`(dataset, route, k, resolution)` 下指标均值和总体标准差；`stability.tsv`
保存同一设置不同种子之间的两两 ARI。不会挑选得分最高的种子。
仅运行一个种子时，稳定性字段为空并注明至少需要两个种子。

## 两个基线共有的输出与评分

两个输出根目录都包含 `metadata.json`、`status.json`、汇总 `metrics.tsv`；
每条路线还包含 `samples.tsv` 和 `cosine_distances.npy`。元数据记录参数、
依赖版本、代码哈希、输入哈希、距离分位数及产物哈希。

| 指标 | 解释 |
| --- | --- |
| `homogeneity` | 一个预测簇是否混入不同真实 GCF |
| `completeness` | 同一真实 GCF 是否被拆到多个预测簇 |
| `v_measure` | homogeneity 与 completeness 的调和平均；未校正随机一致性 |
| `ari`、`ami` | 经随机一致性校正的划分一致性指标 |
| `silhouette_cosine` | 原始 cosine 空间内的轮廓系数；不使用投影后的坐标 |
| `n_clusters`、`n_singleton_clusters` | 预测簇数及只有一个成员的预测簇数 |
| `singleton_sample_fraction` | 被分为 singleton 的样本比例 |
| `coverage` | 当前固定清单的保留比例；两个基线完整运行时为 1 |

只有一个簇或所有样本各自成簇时，silhouette 未定义，表中留空并写出原因，
不会填 0。singleton 不会从 GCF 评分或候选集合中删除。

如需改变特征或清单位置，两脚本都支持：

```bash
  --features-root /path/to/features_be_v1 \
  --cohort /path/to/pilot_be_v1/cohort.tsv
```

两个脚本均拒绝覆盖已有输出，也没有阶段恢复参数。修改参数或重跑时使用新的
`--output`，例如 `data/processed/clustering_average_v2`。运行中断可能留下
`status.json` 的 `running` 状态；只有根目录 `metadata.json` 和 `status.json`
都写明 `complete` 时，才把该目录作为完整结果使用。不要将未完成的成员表当作全队列结果。

## 比较边界

参数曲线用于检查层次、社区和稳定性，不代表可以把同一测试标签上的最高分
作为独立泛化性能。若使用 B/E 标签选择阈值、k 或 resolution，后续报告应明确
这是开发集诊断；可迁移的评测需要在独立开发数据上确定选择规则。

E 来自 MIBiG，BGC-MAC 可能训练接触过其中样本；B/E 之间也有跨数据集重复。
当前分别保留原始队列，不能直接把 B/E 当作互不重叠的训练/测试划分。
GCF 一致性也不等价于连续的产物化学相似度。

本步骤比较的是向量与聚类方法的组合，尚未在这些脚本中调用 BiG-SCAPE。
实际运行是否完整、所用依赖版本和参数，应以对应输出的元数据和日志为准。

## 本地验证记录（2026-11-11）

在 Python 3.12 隔离环境中，66 项测试全部通过。两个脚本均已按默认参数完成
六组真实特征的运行，输出位于上述两个 `*_v1` 目录：

| 基线 | B 每路线配置数 | E 每路线配置数 | 三路线合计 |
| --- | ---: | ---: | ---: |
| Average-linkage | 39 个切点 | 281 个切点 | 960 |
| kNN + Leiden | 90 个参数/种子组合 | 90 个参数/种子组合 | 540 |

各配置保留 B 全部 39 个及 E 全部 281 个样本。已独立核对成员表的 ID、行号、
原始标签、簇数，以及产物/代码哈希；两个基线使用的六份 cosine 距离矩阵
逐元素完全一致。Average 的树图与阈值曲线也已生成。这些检查验证运行和输出
一致性，不代表已经确定了最佳参数或证明某种模型更优。
