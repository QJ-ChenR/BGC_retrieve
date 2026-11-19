# BGC 特征与人工 GCF 聚类比较

已建立九个数据集的可追溯清单。B、E 的 ESM2、BGC-MAC 和 DNA 特征已经在
服务器提取并传回本地，位于 `data/processed/features_be_v1/`。
当前实现两个原始空间基线：cosine＋average-linkage、cosine＋kNN＋Leiden。
已汇总其 V-measure，并导入作者公开的 BiG-SCAPE 分数作为外部参考。
当前新增 PCA 后的 cosine＋average-linkage 比较，只测试 average-linkage。
G/H 的 domain 坐标修复入口见 [修复说明](docs/repair_gh.md)，修复副本写入
`data/processed/repaired_gh_v1/`，原始数据保留在 `data/raw/`。
G/H 全部 101 个输入此前通过了 CPU 解析，但 H 的部分 `aSDomain` 被 Geneious
保存成 `misc_feature`，没有进入蛋白/domain 输入；旧 `pilot_gh_v1/` 的解析通过
不代表 domain 注释已完整恢复。现已在新副本中恢复 **2363/2365 条**，
超过 1021 aa 的 token 从 **156 降为 0**（最长 1008 aa）；`UOC915` 的两条仍待复核，
整套 H 状态为 `needs_review`。见 [恢复说明](docs/reconstruct_h.md)；
暂勿用旧 H 队列启动新的三路线比较。
服务器桥接入口为 `pilot_features.py prepare-reconstructed`，保留全部 90 个 H 样本，
并将 `UOC915` 的已知注释缺失带入新队列与运行记录。新队列和模型提取均由用户
在服务器运行，命令见 [G/H 操作说明](docs/gh_features.md)；本地尚未发布新队列。
使用 Pfam 坐标重新分割的方案暂时搁置；尚未运行 G/H 模型推理或聚类。
其他数据集仍需进一步准备；本地 BiG-SCAPE 运行、微调和融合尚未开展。

聚类入口与说明：

- [cluster_average.py](scripts/cluster_average.py)：完整层次树、距离阈值扫描及树图。
- [cluster_leiden.py](scripts/cluster_leiden.py)：kNN 图、Leiden 社区及随机种子稳定性。
- [聚类使用说明](docs/clustering_baselines.md)：数据格式、参数、指标和输出。
- [PCA＋average-linkage](docs/pca_average.md)：固定维数/方差保留率扫描，以及原始、仅中心化、全部有效主成分对照。

在本项目根目录、已安装 `requirements-clustering.txt` 的环境中运行：

```bash
python scripts/cluster_average.py
python scripts/cluster_leiden.py
```

默认分析 B/E × ESM2/MAC/DNA 六组特征，输出到
`data/processed/clustering_average_v1/`、`data/processed/clustering_leiden_v1/`。
本地已完成这两个默认运行，全部样本保留，66 项测试通过。
已有输出时用 `--output` 指定新的目录。两脚本均按显式行号与 sample_id 对齐，
要求完整样本覆盖；家族标签只参与评分，不用于选参数或构建距离/图。

V-measure 对照见 **[汇总报告](data/processed/vmeasure_comparison_v1/summary.md)**、
[比较图](data/processed/vmeasure_comparison_v1/comparison.png)；入口为
[summarize_vmeasure.py](scripts/summarize_vmeasure.py)。汇总从成员表重算并核对了
1500 个划分的 homogeneity、completeness、V-measure 和 ARI，误差不超过 `1e-12`；
B 保留 39 个样本 / 23 个家族，E 保留 281 个样本 / 90 个家族。

汇总采用真实家族数辅助的诊断口径：Average 取预测簇数最近的切点，并列取较小阈值；
Leiden 对每个 k 分别最小化五个种子的平均 `|预测簇数−真实家族数|`，并列取较小 resolution。
简表固定显示原网格中间值 k=10，完整保留 k=5/10/15，均未按最高 V-measure 挑选。
精确匹配家族数的 Average 中，BGC-MAC 在 B/E 分别为 **0.8700 / 0.8133**；
E 的 Leiden 现有网格最多只有 50 簇，尚未覆盖 90 个家族的目标分辨率。
另有 4 项汇总选择规则测试通过。

BiG-SCAPE 参考来自 [作者 Zenodo 归档](https://zenodo.org/records/16838820)，
含 1.x 的 2 种模式、2.x 的 7 种模式，共 18 条 B/E 结果。本地没有重跑 BiG-SCAPE；
公开日志的名义样本数相符，但有效 ID 未逐条核验，且部分导出混淆矩阵计数不完整，
因此暂不能将并列分数解释为严格同样本比较。来源和核验记录保存在 `sources/`。

复算汇总时指定新目录（脚本拒绝覆盖已有报告）：

```bash
python scripts/summarize_vmeasure.py \
  --published data/processed/vmeasure_comparison_v1/sources/bigscape_published.tsv \
  --output data/processed/vmeasure_comparison_v2
```

PCA 实验使用相同特征与队列，在已安装聚类依赖的环境中运行：

```bash
python scripts/cluster_pca_average.py
```

默认扫描 `2/4/8/16/32/64/128/256` 维和 `90%/95%/99%` 方差保留率，
超过有效秩的请求会明确跳过。输出到新的 `data/processed/clustering_pca_average_v1/`。
PCA 按数据集/路线单独拟合，不做 whitening；同时报告仅中心化对照，
用来区分中心化与丢弃主成分的影响。所有条件沿用最近真实家族数的切点规则，
不按最高 V-measure 选择维数。命令、产物和解释限制见 [PCA 实验说明](docs/pca_average.md)。
默认运行已完成 72 个条件、12,488 个距离切点；9 个超过 B 有效秩的请求明确跳过。
全部划分保留完整队列，六组原始空间对照与已有 average 结果核验一致，27 项相关测试通过。
结果见 [PCA 汇总报告](data/processed/clustering_pca_average_v1/summary.md) 和
[维数曲线](data/processed/clustering_pca_average_v1/v_measure_by_dimension.png)。
本地重跑请用 `--output` 指定新目录。

服务器操作见 **[B/E 特征提取说明](docs/pilot_features.md)**，入口为
`python scripts/pilot_features.py prepare / extract / validate` 三个子命令。
本地已准备的 [pilot_be_v1/cohort.tsv](data/processed/pilot_be_v1/cohort.tsv)
包含 B 的 39 个和 E 的 281 个条目；输入 GBK 使用稳定 sample_id 命名，内容与原始文件逐字节一致。

G/H 数据修复使用 [`repair_gh.py`](scripts/repair_gh.py)。它只在序列和父 CDS
共同支持时修正 domain 端点，并保留逐条证据及修改前后哈希。输出保留原始文件名和
家族目录。运行和输出细节见
[G/H 坐标修复](docs/repair_gh.md)。

G/H 修复副本已完成：G 的蛋白解析由 8/11 提升到 **11/11**，H 由 75/90
提升到 **90/90**，DNA 解析全部通过。20 个文件共修正 884 条 domain 坐标，
原始 2941 个文件的哈希前后一致。详见
[修复结果](data/processed/repaired_gh_v1/summary.md) 和
[逐条修改记录](data/processed/repaired_gh_v1/repairs.tsv)。H 的两条无父 CDS
Pfam 注释保留原样并标记待核查，不影响当前蛋白解析。

`pilot_features.py prepare-repaired` 此前将 G/H 全部 101 个坐标修复副本准备为稳定 ID 队列，
位于 [pilot_gh_v1](data/processed/pilot_gh_v1/summary.md)，metadata 状态为 `complete`，
再次加载校验通过。11 项 G/H 准备测试与 18 项原有提取/校验测试全部通过。
这证明了来源校验和 CPU 解析可运行，但当时没有恢复 H 中隐藏的 Geneious aSDomain。
旧队列作为历史快照保留，不覆盖其文件或 metadata；H 新比较应使用恢复注释后重新
准备的版本。队列中的 `source_sha256`
对应修复后文件，`raw_source_sha256` 单独追溯原始备份；嵌入来源清单后可整目录
同步到服务器。新运行仍使用 ESM2 BOS/1021-aa 窗口策略，MAC 复用同批新缓存，
DNA 使用无长度过滤的 full/FA2/BF16。命令、输入规模和已验证边界见
[G/H 服务器操作](docs/gh_features.md)。当前聚类与评分脚本仍仅注册 B/E，待 G/H
特征返回后再适配 average/PCA 与评分，不新增 Leiden 测试。

H 的恢复入口为 [`reconstruct_h.py`](scripts/reconstruct_h.py)，默认从
`repaired_gh_v1/` 读取 H 的全部 90 个副本，写到新的
`data/processed/reconstructed_h_v1/`，保留 G、raw 和既有 processed 快照。
它只恢复带有确切 `Geneious type: aSDomain` 标记、且父 CDS 与阅读框/翻译证据
支持的旧注释；保留重复条目并记录旧父 ID 映射，不新增 Pfam 扫描或分割。
实际结果见 [恢复摘要](data/processed/reconstructed_h_v1/summary.md)：90 个样本全部
保留，89 个通过 aSDomain 恢复检查，`UOC915` 待复核；无超长 token 不代表其缺失
注释已解决。本次 22 项测试通过，未生成新模型输入队列或特征。
逐条证据和后续输入要求见 [H aSDomain 恢复说明](docs/reconstruct_h.md)。

新增 `pilot_features.py prepare-reconstructed` 从 `reconstructed_h_v1/` 读取 H，
并可从 `repaired_gh_v1/` 合入 G，按稳定 ID 准备完整队列。随后现有 `extract`、
`validate` 直接使用新输出。建议服务器目录分别为
`pilot_gh_reconstructed_v1/` 和 `features_gh_reconstructed_v1/`。
机械准备的 `complete` 与注释的 `needs_review` 分开记录；H 不自动排除
`UOC915`，也不把未恢复的两条注释改标为完整。此桥接已在 `/tmp` 用真实 G/H
101 个样本及 H-only 90 个样本通过 CPU 准备、加载与迁移验证；真实输入哈希不变，
未加载模型。上述正式队列和模型结果尚待用户在服务器创建。

## 目录布局

```text
data/
  raw/
    curated_datasets/          # 原始下载内容，整目录迁移，文件内容保持不变
      README.md
      Dataset_*/
        curated_GCF_assignments.tsv
        GCF_assignments/<原始家族目录>/*.gbk
  processed/
    raw_relocation/            # 迁移前后路径、逐文件校验和及验证记录
    inventory_v1/              # 当前样本清单
    pilot_be_v1/               # B/E 固定输入、身份映射与原始标签
    features_be_v1/            # 服务器返回的 B/E × ESM2/MAC/DNA 特征
    clustering_average_v1/     # 原始 cosine average-linkage 分析
    clustering_leiden_v1/      # 原始 cosine kNN Leiden 分析
    vmeasure_comparison_v1/    # 指标复核、选择规则、比较图和 BiG-SCAPE 来源
    clustering_pca_average_v1/ # PCA＋average 的默认输出位置
    repaired_gh_v1/            # G/H 坐标修复副本、逐条证据及 raw 完整性校验
    reconstructed_h_v1/        # H Geneious aSDomain 恢复副本；使用前检查 metadata 状态
    pilot_gh_v1/               # 历史队列，H 尚未恢复 Geneious aSDomain；勿用于新 H baseline
    pilot_gh_reconstructed_v1/ # 新 G/H 队列建议位置，待服务器 prepare-reconstructed 创建
    features_gh_reconstructed_v1/ # 新 G/H 三路线特征建议位置，尚未运行
    _history/                 # 首次生成、复核前的派生清单备份
scripts/
  prepare_raw_layout.py
  build_inventory.py
  validate_inventory.py
  pilot_features.py           # B/E 与修复后 G/H 准备、三路线提取、恢复及校验
  prepare_repaired_features.py # G/H 修复来源验证、CPU 输入检查与可迁移队列
  prepare_reconstructed_features.py # H 恢复来源桥接，保留待复核状态，可合入 G
  pilot_feature_validation.py # ID、向量、来源和 full-DNA 覆盖校验
  clustering_common.py        # 共享特征读取、ID/来源校验、cosine 与评分
  cluster_average.py
  cluster_leiden.py
  summarize_vmeasure.py       # 从已有成员表复核指标并汇总，不重跑聚类
  cluster_pca_average.py      # PCA 维数扫描及 average-linkage 对照
  pca_transform.py            # full SVD、数值秩与方差保留率选择
  repair_gh.py                # G/H 修复编排与前后验证
  domain_coordinate_repair.py # 有序列证据的 domain 端点修复
  reconstruct_h.py            # H 既有 Geneious aSDomain 恢复与逐条审计
  geneious_asdomain_restore.py # CDS/阅读框/翻译证据与最小文本恢复
tests/
```

原来的 `data/curated_datasets/` 已移至 `data/raw/curated_datasets/`，没有保留第二份原始目录。
迁移覆盖全部 **2941 个文件、631205225 字节**，包括 2911 个 GBK、9 个标签 TSV、
10 个 HTML、10 个 JSON 和 1 个 README；迁移前后逐文件 SHA256、大小和路径集合均已核对。
所有后续派生数据放在 `data/processed/`；不要改写 `data/raw/` 中的序列、注释或标签。

## 当前清单

首先查看 [汇总](data/processed/inventory_v1/summary.md) 和
[逐数据集统计](data/processed/inventory_v1/dataset_summary.tsv)。

| 文件 | 内容 |
| --- | --- |
| `samples.tsv` | 2911 个 GBK 文件条目：稳定 sample_id、来源路径、文件 hash、目录/TSV 原始标签、关联状态、record 数、解析器兼容性和错误 |
| `records.tsv` | 2934 个 record：原始顺序和 ID、DNA hash、反向互补规范化 hash、长度、GC、模糊碱基、物种/分类标识、CDS/domain 与翻译诊断 |
| `label_rows.tsv` | 2911 条 TSV 原始行：文件 hash、行号、原始 BGC/GCF 字符串以及对应/候选文件 ID |
| `source_files.tsv` | 全部 2941 个原始文件的当前路径、归档内路径、字节数和 SHA256 |
| `duplicate_groups.tsv` | 按不同证据定义的重复组与同一数据集内标签冲突候选；不合并样本 |
| `duplicate_members.tsv` | 重复组与 sample_id / record_uid 的对应关系 |
| `issues.tsv` | 名称重复、标签差异、多 record 和解析失败等待核查事项；一条样本可能出现多次 |
| `dataset_summary.tsv` | 每个数据集的文件/标签/record 数、兼容率和关联情况 |
| `metadata.json` | 版本、策略、输入快照、脚本/解析器源码 hash、依赖版本及输出校验和 |

文件条目并不等于独立 BGC。每个多 record GBK 在 `samples.tsv` 中仍只有一条，
其内部 records 逐条写入 `records.tsv`。不会拼接序列或将文件标签分配为多个独立训练样本。

### 身份与标签关联

- `sample_id` 是 `schema_version + NUL + archive_relpath` 的 SHA256 前 24 个十六进制字符，
  加 `bgc_` 前缀。`archive_relpath` 相对于 `curated_datasets/`，移动工作区不会改变 ID。
  它标识来源文件条目；同名文件不会覆盖，不根据序列去重。源文件字节变化由独立的 `source_sha256` 检出。
- `record_uid = sample_id:record:<1-based index>`，保留原始文件中的 record 顺序。
- `row_id` 标识数据集及 TSV 原始行号；`tsv_row_ids_json` 可连接回原始标签行。
- 原始标签保留 Unicode 空白，包括 Dataset B 的 9 条尾随 NBSP 标签。没有 `.strip()` 后重写原值。
- `exact_path_label`：数据集、原始 basename 和家族目录与 TSV 唯一精确对应，共 2904 条。
- `unique_basename_folder_mismatch`：文件和 TSV 各有一个同名条目，但目录与标签不同，共 7 条。
  当前均为 I 的目录 `other_11.2` 与 TSV `other_11`。这只建立来源关联，不认定两标签是生物学同义词。
- 通用脚本还支持 `ambiguous` 和 `missing`：保留候选行，绝不选择第一/最后一个标签。
  当前这两种状态均为 0，但生物学重复和多标签冲突仍未解决。

因此，**文件—TSV 对应成功，不代表同名/同序列条目的单标签真值已确定**。

### 哈希与重复组的含义

| group_type | 分组证据 |
| --- | --- |
| `basename` | 同一数据集内文件 stem 相同；不表示序列相同 |
| `file_bytes` | 原始 GBK 字节 SHA256 相同 |
| `ordered_file_dna` | 文件内有序 record 序列 hash 列表相同；有边界，不是 DNA 拼接 hash |
| `single_record_dna_rc` | 单 record 文件的 DNA 相同或互为反向互补 |
| `record_dna` | 单个 record 的正向 DNA 相同 |
| `record_dna_rc` | 单个 record 的 DNA 相同或互为反向互补 |
| `protein_input` | 现有蛋白解析器成功后生成的有序 token 序列完全相同 |

DNA hash 基于大写的原始序列；反向互补规范化只影响审计分组，没有翻转实际输入。
GC 分母为 A/C/G/T 总数；模糊碱基保留并另外计数。
`protein_input` 不包括被现有 domain 切分遗漏的残基，因此不表示所有 CDS 或完整 DNA 相同。
`feature_annotation_sha256` 记录有序 feature 的类型、位置及 qualifiers，不是全文件 hash。

当前有 105 个数据集内重名组、146 个正向 record 重复组和 150 个含反向互补的 record 重复组；
其中分别有 21、22 个 record 序列组跨数据集。不同分组类型会重叠，不能相加当作独立重复数。
`within_dataset_label_conflict` 仅提示同一数据集内的不同标签，跨数据集标签不当作统一家族词表。
尚未做序列相似性搜索或近同源分组。

### 兼容性与覆盖

全部 2911 个 GBK 能通过基础多 record 读取。现有 DNA 单文件解析器通过 2893 个，
18 个多 record 文件失败；现有蛋白单文件解析器通过 2694 个，217 个失败。
这些是 CPU 解析检查，未调用模型，也没有绕过验证或修复注释。

失败条目仍在清单中，`excluded=False`、`repairs_applied_json=[]`。
DNA 不设 bp/token 长度上限；最长 record 为 281670 bp，已完整登记。
本清单不是可直接交给编码器的最终输入目录，也不是已完成去重和标签裁定的评测队列。

## 复现与校验

本次构建使用 Python 3.12.14、Biopython 1.88。迁移和验证脚本只依赖标准库；
构建清单需要 Biopython，依赖记录在 `requirements-inventory.txt`。

在项目根目录执行。迁移脚本对已完成的迁移只会重新验证，拒绝覆盖或合并现有目录：

```bash
python scripts/prepare_raw_layout.py --workspace .
python scripts/validate_inventory.py --workspace .
```

清单构建拒绝覆盖已有输出。复核或调整策略时创建新的版本目录，例如：

```bash
python scripts/build_inventory.py \
  --workspace . \
  --encoder-project ../BGC_annotation \
  --output data/processed/inventory_v2
python scripts/validate_inventory.py --workspace . --inventory data/processed/inventory_v2
```

`--encoder-project` 仅加载 `read_dna` 和 `parse_bgc` 两个单文件解析 API，
不会加载模型权重或使用 GPU。省略它时，兼容状态为 `not_checked`，不能把通过数 0 解释为全部失败。
构建前后重新核对全部原始文件，出现变动时不发布完成结果。
迁移报告留有原始路径 `original_path` 到当前 `raw_path` 的映射；中断迁移不会自动覆盖或猜测恢复。

本地实际使用的解释器为 `/home/qiju/miniforge3/envs/BGC/bin/python`。
测试命令：

```bash
python -B -m unittest discover -s tests -v
```

其他七个数据集仍应根据 `issues.tsv`、重复组和兼容性错误逐类确定处理政策。
B/E 当前按原始队列分别提取，保留 singleton 和跨数据集重复；不把两者当作互不重叠的训练/测试划分。
这两个聚类基线已进入实现阶段；其他聚类方案、ESM2 微调和 DNA/蛋白融合仍是未来路线，详见 `AGENTS.md`。
