# G/H 的 ESM2、BGC-MAC 与 DNA 特征提取

现有 `BGC_annotation` 编码器可解析修复后的全部 G/H 输入，无需修改模型或解析器。
本项目的 `pilot_features.py` 已增加 `prepare-repaired`，将坐标修复副本转换为与
B/E 相同的稳定 ID 输入布局；后续沿用 `extract` 和 `validate`。
本机检查只使用 CPU 解析器，尚未执行 G/H 模型推理或得到聚类性能。
本机九个编码器相关源码文件的哈希与 B/E 已完成运行记录一致；扩展后的校验入口
也重新通过了已有 B/E 六组导出的检查。
固定队列已发布到 `data/processed/pilot_gh_v1/`，metadata 状态为 `complete`，
G 的 11 个和 H 的 90 个样本全部通过再次加载校验；见
[准备摘要](../data/processed/pilot_gh_v1/summary.md) 和
[逐样本输入检查](../data/processed/pilot_gh_v1/input_checks.tsv)。
11 项 G/H 准备测试和 18 项原有提取/校验测试全部通过。

## 接口检查结果

| 项目 | G：Planomonospora | H：Reclassification |
| --- | ---: | ---: |
| BGC / 人工 GCF | 11 / 3 | 90 / 10 |
| 蛋白和 DNA 解析通过 | 11 / 11 | 90 / 90 |
| DNA 长度范围 / bp | 15,484–76,182 | 43,475–99,572 |
| 每个 BGC 的蛋白/domain token 数 | 15–71 | 13–92 |
| 最长蛋白/domain token / aa | 945 | 9,508 |
| 超过 1021 aa 的 token 数 | 0 | 156（涉及 42 个 BGC） |
| ESM2 1021-aa 窗口总数 | 372 | 5,196 |
| 保留的 N 碱基数 | 0 | 246 |

- **ESM2**：长蛋白/domain 由现有接口按最多 1021 aa 的不重叠窗口完整覆盖，
  窗口使用 BOS 向量，再按残基数加权还原为一个原始 token 向量。
  H 的 9,508 aa token 可以进入该流程，但它不是整条蛋白一次全上下文编码。
- **BGC-MAC**：复用本次 ESM2 导出的 token cache。分窗后的向量会先合并，
  因此 H 的 MAC 输入最多仍为每个 BGC 92 个 token。MAC 的位置编码按实际
  token 长度计算，不存在 G/H 数据集名称白名单。
- **DNA**：使用 full 模式、不设 bp/token 上限、不截断、无自动 windowed 回退；
  N 碱基保留，实际 tokenizer 的 UNK 统计由服务器导出记录。
  H 最长 99,572 bp，低于此前 B/E 已成功的最长 185,250 bp 输入；这提供资源
  经验参考，实际 DNA token 数和显存使用仍需在同一服务器环境中验证。

H 的 `ApSA.gbk` 仍有两条无父 CDS 的 `PFAM_domain` 待核查。现有蛋白解析器
不使用这种 feature，因此它们不阻止编码；检查通过不代表全部生物学注释已解决。

## 准备可迁移的固定队列

本地 `pilot_gh_v1/` 已准备完成，可直接整目录同步到服务器。下面是复现命令；
已有输出时需换用新的 `--output`，例如 `pilot_gh_v2`。
从本项目根目录运行，环境需有 `requirements-inventory.txt` 中的 CPU 依赖：

```bash
python scripts/pilot_features.py prepare-repaired \
  --repaired-root data/processed/repaired_gh_v1 \
  --encoder-project ../BGC_annotation \
  --datasets G H \
  --output data/processed/pilot_gh_v1
```

该命令只接受完整的 G 和/或 H 队列。若输出已存在，使用新版本目录，不覆盖已有
结果。`--repaired-root` 指向包含 `metadata.json` 和两个数据集目录的修复快照；
不能指向 raw 或其中某个 `GCF_assignments/`。旧的 `prepare` 仍用于 B/E 原始输入。

准备步骤核对修复快照的产物哈希、完整文件集合、逐条修复记录和标签，再复制 GBK，
执行实际 CPU 蛋白/DNA 解析并核对序列。保留全部 101 个样本、singleton 和原始标签，
不按长度筛选、去重或进一步修改注释。模型读取目录只包含按稳定 `sample_id` 命名的 GBK。

```text
data/processed/pilot_gh_v1/
  cohort.tsv
  input_checks.tsv
  metadata.json
  summary.md
  repair_provenance/
    metadata.json
    samples.tsv
    repairs.tsv
    review_items.tsv
    Dataset_*/curated_GCF_assignments.tsv
  Dataset_G_Planomonospora/
    inputs/bgc_<stable_hash>.gbk
    curated_GCF_assignments.tsv
  Dataset_H_Reclassification/
    inputs/bgc_<stable_hash>.gbk
    curated_GCF_assignments.tsv
```

`cohort.tsv` 保留稳定 ID、原始 ID/路径/GCF，以及两种哈希：

- `source_sha256`：模型实际读取的**修复后 GBK** 字节哈希；与编码器 metadata 核验。
- `raw_source_sha256`：原始备份字节哈希，用于追溯；修复过的文件与前者不同。

`input_checks.tsv` 保存每个样本的完整 DNA、蛋白输入序列哈希、解析状态、长度、
CDS/token 数和预计 ESM 窗口数。`repair_provenance/` 嵌入来源快照和原始 ID 标签。
**复制完整的 `pilot_gh_v1/` 到服务器即可继续**；加载时不依赖本机原来的 raw/repaired
绝对路径。不要单独复制 GBK 而遗漏来源清单。

## 在服务器提取三条路线

同步本项目完整 `scripts/` 和准备好的 `data/processed/pilot_gh_v1/`。
沿用完成 B/E 的服务器环境、编码器源码和三个 checkpoint。下列路径来自已返回的
`features_be_v1/run.json`；如果服务器目录已移动，只调整相应路径，并保持权重与配置一致。

在服务器 `/home/s252613/BGC_retrieve` 根目录执行：

```bash
/home/s252613/miniforge3/envs/natural_product/bin/python scripts/pilot_features.py extract \
  --prepared data/processed/pilot_gh_v1 \
  --encoder-project /home/s252613/BGC_encoder/BGC_annotation \
  --output data/processed/features_gh_v1 \
  --esm-checkpoint /home/s252613/BGC_encoder/BGC_annotation/data/esm2_t33_650M_UR50D.pt \
  --mac-checkpoint /home/s252613/BGC_encoder/BGC_annotation/ckpt/BGC_MAC/MAC_default/MAC_default.ckpt \
  --dna-checkpoint /home/s252613/checkpoint-13920 \
  --routes esm2 mac dna_full \
  --device cuda:0 \
  --esm-batch-size 8 \
  --batch-size 16 \
  --mac-member 0 \
  --dna-dtype bfloat16 \
  --dna-attention-backend flash_attention_2 \
  --trusted-input
```

`--trusted-input` 沿用 B/E 对可信旧 ESM2/MAC checkpoint 的加载设置。ESM2 固定
BOS、1021-aa 窗口；MAC 使用 member 0；DNA 使用 full、FA2/BF16。
提取配置会记录 checkpoint、tokenizer/config 和编码器源码哈希，便于核对与 B/E
一致性。实际维度以权重和导出 metadata 为准；同一组权重预计为 1280 / 256 / 768。

提取完成后校验：

```bash
/home/s252613/miniforge3/envs/natural_product/bin/python scripts/pilot_features.py validate \
  --prepared data/processed/pilot_gh_v1 \
  --output data/processed/features_gh_v1
```

每个数据集分别产生 `esm2/`、`mac/`、`dna_full/`。每条路线都必须有 G 的 **11 行**
或 H 的 **90 行**，原生文件仍为 `embeddings.npy`、`bgc_ids.tsv`、`metadata.json`、
`failures.tsv`，另有 `execution.log`。`esm2/tokens.pth` 由 MAC 复用；顶层
`run.json` 保存六个阶段的配置、状态和完整覆盖验证结果。
不能把部分成功的 ID 交集当作完整数据集结果。

中断后使用完全相同的命令追加 `--resume`。脚本验证并复用已完成阶段，失败阶段
归档后重跑；不是逐样本恢复。更改代码、checkpoint 或运行配置时使用新输出目录。
输出、恢复和校验的详细规则见 [B/E 操作说明](pilot_features.md)。

## 后续性能比较

服务器运行后，将整个 `features_gh_v1/` 传回本地，保留 `run.json`、metadata 和
ESM cache 的来源证据。当前这里只完成输入接口适配，尚无 G/H embedding 或 V-measure。
已有 average/PCA 和 V-measure 汇总脚本仍按 B/E 注册，不能仅替换目录就假定能分析 G/H。
特征返回后需扩展其数据集选择与固定队列校验，再沿用原始 cosine＋average-linkage、
PCA＋average-linkage 和同一评分口径；本阶段不增加 Leiden 测试。

比较 BiG-SCAPE 时也要核验实际样本覆盖和坐标版本。公开作者结果仍是外部参考，
不能因为名义样本数一致就称为严格同样本复现。人工 GCF 是本次聚类的辅助评价标签，
分数不等同于产品功能或独立泛化性能。
