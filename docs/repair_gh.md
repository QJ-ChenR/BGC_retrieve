# G/H domain 坐标修复

入口为 [`scripts/repair_gh.py`](../scripts/repair_gh.py)，坐标证据与最小文本修改由
[`scripts/domain_coordinate_repair.py`](../scripts/domain_coordinate_repair.py) 实现。
修复对象是 G/H 中有相同坐标错误证据的 `aSDomain` 和 `PFAM_domain` 注释，
原始 DNA、翻译序列和 GCF 标签不重写。

## 运行

在项目根目录、已安装 `requirements-inventory.txt` 的环境中运行：

```bash
python scripts/repair_gh.py \
  --workspace . \
  --encoder-project ../BGC_annotation \
  --output data/processed/repaired_gh_v1
```

只处理 `Dataset_G_Planomonospora` 和 `Dataset_H_Reclassification`。
输出必须位于工作区的 `data/processed/`，且必须是新目录。重跑时换用
`--output data/processed/repaired_gh_v2`，不覆盖原始数据或已有修复结果。
脚本先在临时目录完成修复和验证，再发布完整目录；验证失败不会发布完成结果。
验证调用现有编码器项目中的 CPU 解析函数，不加载模型或使用 GPU。

## 修复依据

部分旧 domain 的结束坐标漏掉了最后一个氨基酸及其密码子。
这次逐条检查以下证据，不对整个数据集统一执行 `end += 1`：

1. domain 有唯一、同链且包含它的父 CDS，标识符与父子关系无歧义。
2. domain 注释翻译与父 CDS 的对应蛋白区间完全相同。
3. 原 domain 核酸区间的翻译等于 domain 翻译去掉最后一个氨基酸。
4. 沿蛋白 C 端延长一个密码子后，核酸翻译与完整 domain 翻译完全相同，
   且新位置仍在父 CDS 内，并与蛋白区间、链方向和 CDS 阅读框一致。

有显式旧蛋白坐标的条目，在上述证据通过后将 `protein_end` 改为原值加 1，
同时延长对应的核酸端点：正链 `end + 3`，负链 `start - 3`。
另有少量 G 条目没有 `protein_start`/`protein_end`，但翻译在父 CDS 中只有一个匹配；
这类条目只修复有相同证据的核酸端点，不新增原本不存在的蛋白坐标。

脚本只替换原 GBK 中必要的数字，不通过重新序列化整份 GenBank 文件改写格式。
原本正确的 domain 不修改；`PFAM_domain` 保持原类型，不转换为 `aSDomain`，
也不会因此进入现有 ESM domain token 输入。
DNA、CDS/domain 的 translation、其他注释、文件名、家族目录和标签均保留。
无法证明的修复不执行：影响当前蛋白解析的错误会阻止发布；解析器不使用的孤立
`PFAM_domain` 保留原样并标记 `needs_review`，不猜坐标、不删除 feature 或样本。

## 输出与追溯

```text
data/processed/repaired_gh_v1/
  Dataset_G_Planomonospora/
    curated_GCF_assignments.tsv
    GCF_assignments/<family>/<original_name>.gbk
  Dataset_H_Reclassification/
    curated_GCF_assignments.tsv
    GCF_assignments/<family>/<original_name>.gbk
  samples.tsv
  domain_checks.tsv
  repairs.tsv
  review_items.tsv
  raw_integrity.tsv
  dataset_summary.tsv
  metadata.json
  summary.md
```

| 文件 | 内容 |
| --- | --- |
| `samples.tsv` | 稳定 sample_id、原始 BGC ID/路径/标签、raw 和 processed 哈希、修复数量及前后解析状态 |
| `domain_checks.tsv` | 检查过的 domain 与父 CDS、坐标和序列匹配证据 |
| `repairs.tsv` | 每个修改 feature 的旧/新坐标、规则和序列证据 |
| `review_items.tsv` | 证据不足、保留原样的待核查注释 |
| `raw_integrity.tsv` | 原始文件在本次处理前后的完整性校验 |
| `dataset_summary.tsv` | 按数据集统计输入数、修改数及解析通过数 |
| `metadata.json` | 输入清单、脚本/解析器版本与哈希、验证结果、输出校验值 |
| `summary.md` | 结果摘要及使用限制 |

修复副本保留 G/H 全部 101 个文件条目，不按 DNA 长度过滤、去重或重分配标签。
原始 `data/raw/` 和 `inventory_v1/` 保持原状；清单中旧的解析失败记录仍描述原始文件，
修复后的状态应查看本目录的 `samples.tsv`。

## 当前状态与后续使用

2026-11-19 已完成实际运行，结果位于
[`repaired_gh_v1/summary.md`](../data/processed/repaired_gh_v1/summary.md)：

| 数据集 | 保留 GBK 数 | 修改文件数 | 修复 domain 数 | 蛋白解析：修复前 → 后 |
| --- | ---: | ---: | ---: | --- |
| G | 11 | 5 | 266 | 8/11 → 11/11 |
| H | 90 | 15 | 618 | 75/90 → 90/90 |

DNA 解析全部通过。共修复 668 条 aSDomain 与 216 条 PFAM_domain 的核酸末端，
其中 866 条同步修改显式 protein_end，另 18 条只修正核酸坐标。
其余 81 个 GBK 和两个标签 TSV 与 raw 逐字节一致。
全部 2941 个原始文件的哈希与清单一致，修复前后也完全一致；未改写原始清单。
输出经过逐 feature 保留性验证和再次修复不产生变化的检查。
13 项核心修复与编排测试通过，包括正负链、文本保留、证据不足时拒绝修复、
禁止覆盖以及验证失败时不发布部分结果。

审计中 H 的 `ApSA.gbk` 有两条 `PFAM_domain` 无法对应到注释中的父 CDS，
暂不修改并保留复核记录。当前蛋白解析器不使用 `PFAM_domain`，因此它们不阻止
蛋白输入解析；解析通过不代表全部生物学注释问题已经解决。

后续编码通过 `scripts/pilot_features.py prepare-repaired` 读取本目录，生成
`data/processed/pilot_gh_v1/`：输入内容保持为修复副本，文件名改为稳定 sample_id，
并嵌入修复来源清单。再使用相同脚本的 `extract`、`validate` 调用 ESM2、BGC-MAC
和 DNA。不要直接把本目录作为 `--prepared`，因为修复快照与编码队列的格式不同。
完整命令及接口检查见 [G/H 特征提取](gh_features.md)。
`pilot_gh_v1/` 已包含全部 101 个样本并通过加载校验，见
[编码输入准备摘要](../data/processed/pilot_gh_v1/summary.md)。
目前只执行修复、CPU 输入验证与编排适配，没有运行 G/H 模型推理、聚类或 BiG-SCAPE。
