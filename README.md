# BGC 人工 GCF 数据集准备

已建立九个数据集的可追溯清单。当前选择 B、E 跑通 ESM2、BGC-MAC 和 DNA
特征提取；服务器调用脚本与固定输入队列已准备，真实模型提取由用户在服务器运行。
尚未执行聚类、BiG-SCAPE 比较或模型训练，其他数据集仍需进一步准备。

服务器操作见 **[B/E 特征提取说明](docs/pilot_features.md)**，入口为
`python scripts/pilot_features.py prepare / extract / validate` 三个子命令。
本地已准备的 [pilot_be_v1/cohort.tsv](data/processed/pilot_be_v1/cohort.tsv)
包含 B 的 39 个和 E 的 281 个条目；输入 GBK 使用稳定 sample_id 命名，内容与原始文件逐字节一致。

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
    _history/                 # 首次生成、复核前的派生清单备份
scripts/
  prepare_raw_layout.py
  build_inventory.py
  validate_inventory.py
  pilot_features.py           # B/E 准备、三路线提取、阶段恢复与校验
  pilot_feature_validation.py # ID、向量、来源和 full-DNA 覆盖校验
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
多聚类、ESM2 微调和 DNA/蛋白融合仍属于未来研究路线，详见 `AGENTS.md`。
