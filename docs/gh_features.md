# G/H 的 ESM2、BGC-MAC 与 DNA 特征提取

当前服务器入口为 `pilot_features.py prepare-reconstructed` → `extract` → `validate`。
G 使用坐标修复副本，H 使用已恢复已有 aSDomain 的副本；创建新的
`pilot_gh_reconstructed_v1/` 和 `features_gh_reconstructed_v1/`，避免混入旧 H 输入。
本地只实现和验证代码，**尚未发布这套新编码队列，也没有运行 G/H 模型推理**；
下面命令由用户在服务器执行。

H 已恢复 2363/2365 条 Geneious aSDomain，保留全部 90 个样本，`UOC915` 的两条
旧注释仍为 `needs_review`。新桥接保留该样本并携带这一已知注释缺口，
默认保留完整 90 样本队列，并在结果中注明已知注释缺口；不会删除该样本或把待复核状态改成正常。
这次不新增 Pfam 分割、不重写 GBK，也不把机械准备完成解释为所有注释均已修复。
恢复证据见 [H aSDomain 恢复说明](reconstruct_h.md)。

## 当前输入规模

| 项目 | G：坐标修复后 | H：aSDomain 恢复后 |
| --- | ---: | ---: |
| BGC / 人工 GCF | 11 / 3 | 90 / 10 |
| 蛋白和 DNA 解析通过 | 11 / 11 | 90 / 90 |
| DNA 长度范围 / bp | 15,484–76,182 | 43,475–99,572 |
| 每个 BGC 的蛋白/domain token 数 | 15–71 | 48–123 |
| token 总数 | 372 | 6791 |
| 最长 token / aa | 945 | 1008 |
| 超过 1021 aa 的 token 数 | 0 | 0 |
| 保留的 N 碱基数 | 0 | 246 |

H 的 6791 个 token 包括 4718 个 aSDomain 和 2073 条未分割蛋白。恢复前的
`pilot_gh_v1/` 只使用了 4836 个 token，包含 156 条超长整蛋白，最长 9508 aa；
它遗漏了 Geneious 保存为 `misc_feature` 的旧 aSDomain。旧队列保持不变，
不要用于本次 H 运行，也不要将新 GBK 覆盖进去而继续使用旧哈希。
H 中无超长 token 不代表 `UOC915` 两条未恢复注释造成的输入缺失已经解决。

ESM2 继续使用 BOS 和最多 1021 aa 的不重叠窗口；当前 G/H token 均不需要分窗。
MAC 复用同一次提取的 ESM2 cache，H 每个 BGC 最多 123 个 token。
DNA 使用 full/FA2/BF16，不设 DNA 长度过滤、bp/token 上限或自动 windowed 回退。
H 的 `ApSA.gbk` 两条无父 CDS 的 `PFAM_domain` 仍保留；当前蛋白解析器不使用
该 feature 类型，其来源复核记录随恢复副本保留。

## 在服务器准备新队列

同步本项目完整 `scripts/` 和以下**完整快照目录**到服务器：

- `data/processed/repaired_gh_v1/`：G 的来源，以及坐标修复证据。
- `data/processed/reconstructed_h_v1/`：H 的来源，包含恢复后的 90 个 GBK、标签、
  审计 TSV、metadata 和嵌入的来源记录。

同步已有恢复快照即可，不必重新运行恢复脚本。若另行选择在服务器重新执行
`reconstruct_h.py`，则还需完整 raw、inventory 和坐标修复来源；其路径参数见
[恢复说明](reconstruct_h.md)。只准备 H 时仅需完整 H 恢复快照。

使用已有 B/E 服务器环境，需有 `requirements-inventory.txt` 中的 CPU 依赖。
以下路径沿用 B/E 的服务器配置；如目录已移动，只调整相应路径。
从服务器 `/home/s252613/BGC_retrieve` 根目录执行：

```bash
/home/s252613/miniforge3/envs/natural_product/bin/python scripts/pilot_features.py prepare-reconstructed \
  --reconstructed-root data/processed/reconstructed_h_v1 \
  --repaired-root data/processed/repaired_gh_v1 \
  --datasets G H \
  --encoder-project /home/s252613/BGC_encoder/BGC_annotation \
  --output data/processed/pilot_gh_reconstructed_v1
```

只准备 H 时使用 `--datasets H` 并省略 `--repaired-root`，不需要读取 G 的来源；后续 `extract` / `validate`
会按队列实际包含的数据集运行。输出必须是新目录；重跑请换版本，不覆盖已有队列。
`--reconstructed-root` 指向含 metadata 的 H 恢复快照根目录，不能指向其
`GCF_assignments/`、原始 raw 或旧 `pilot_gh_v1/`。

准备步骤核验来源产物哈希、样本身份、标签、恢复记录与待复核状态，复制实际 GBK
字节到稳定 `sample_id` 文件名，并重新执行 CPU 蛋白/DNA 解析与输入检查。
保留 G 的 11 个和 H 的 90 个条目、singleton、重复注释及原始标签；不按长度过滤，
不再修复注释。来源审计随队列嵌入，准备完成后可整目录迁移，不依赖本机来源绝对路径。

队列保存 `cohort.tsv`、`input_checks.tsv`、`metadata.json`、`summary.md`、来源记录，
以及各数据集的 `inputs/bgc_<stable_hash>.gbk` 和原始标签 TSV。
H 的恢复审计嵌入 `reconstruction_provenance/`，其中包含原有
`source_provenance/` 坐标修复链；选择 G 时另嵌入 `repair_provenance/`。
`source_sha256` 对应模型实际读取的 GBK，`raw_source_sha256` 单独追溯原始备份。
不要单独复制 GBK 而遗漏清单和来源文件。

`metadata.status: complete` 表示机械准备与校验完成。H 的 `needs_review`
状态及 `UOC915` 两条未恢复注释会保留在队列 metadata、成员表
和后续 `run.json` 中，不能由准备成功覆盖。检查 `summary.md` 时同时核对完整
样本覆盖和注释限制，而不是只看 `complete`。

新队列 schema 为 `bgc_reconstructed_pilot_inputs_v1`，提取运行 schema 为
`bgc_reconstructed_feature_run_v1`。prepared metadata、提取配置和运行记录均携带：

| 字段 | 当前 G/H 或 H-only 队列含义 |
| --- | --- |
| `annotation_status` | `needs_review`，仍有已知注释缺口 |
| `annotation_complete` | `false` |
| `unresolved_asdomain_count` | `2` |
| `unresolved_sample_ids` | `UOC915` 对应的稳定 sample_id |
| `review_policy` | `retain`，保留待复核样本 |

`annotation_status` / `annotation_complete` 描述本次 Geneious aSDomain 恢复状态。
坐标修复阶段遗留的 PFAM 等复核项另记为 `inherited_review_count`，不合并进这两条
未恢复 aSDomain，也不因其未被蛋白解析器使用而丢弃追溯记录。

## 提取三条路线

在同一服务器目录执行，使用与 B/E 一致的权重和配置：

```bash
/home/s252613/miniforge3/envs/natural_product/bin/python scripts/pilot_features.py extract \
  --prepared data/processed/pilot_gh_reconstructed_v1 \
  --encoder-project /home/s252613/BGC_encoder/BGC_annotation \
  --output data/processed/features_gh_reconstructed_v1 \
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

`--trusted-input` 沿用 B/E 对可信旧 ESM2/MAC checkpoint 的加载设置。
MAC 使用 member 0，并读取这次 ESM2 生成的 `esm2/tokens.pth`；不要复用恢复前的
H ESM2 cache。提取配置记录 checkpoint、tokenizer/config 和编码器源码哈希。
同组权重预期维度为 ESM2 1280 / MAC 256 / DNA 768，实际以导出 metadata 为准。

提取完成后校验：

```bash
/home/s252613/miniforge3/envs/natural_product/bin/python scripts/pilot_features.py validate \
  --prepared data/processed/pilot_gh_reconstructed_v1 \
  --output data/processed/features_gh_reconstructed_v1
```

每个数据集分别产生 `esm2/`、`mac/`、`dna_full/`；每条路线必须包含 G 的 **11 行**
或 H 的 **90 行**，包含 `UOC915`。原生文件仍为 `embeddings.npy`、`bgc_ids.tsv`、
`metadata.json`、`failures.tsv`，另有 `execution.log`；顶层 `run.json` 保存各阶段
配置、状态、覆盖验证和已知注释限制。不得把部分成功的 ID 交集当作完整数据集结果。

中断后使用完全相同的提取命令追加 `--resume`。脚本验证并复用已完成阶段，失败阶段
归档后重跑；不是逐样本恢复。更换输入快照、代码、checkpoint 或配置时使用新输出目录。
输出、恢复和校验的通用规则见 [B/E 操作说明](pilot_features.md)。

## 当前验证与后续性能比较

本地已在 `/tmp` 对真实输入完成 CPU 准备与加载：G 的 11 个、H 的 90 个样本全部
保留，token 总数分别为 372 / 6791，超过 1021 aa 的 token 均为 0，两条待复核
注释的状态完整保留。另已单独验证 H-only 的全部 90 个样本。
验证期间未导入 `torch`、未加载模型，真实输入文件哈希保持不变。

迁移验证在临时副本中将 G 的历史 `raw_root` 改为模拟服务器路径，搬移已准备队列并改名其
原输入副本，再限制 JSON/TSV 读取和哈希检查只能访问新队列；加载仍通过。
这验证了来源身份依赖内容哈希、已准备队列不依赖旧绝对路径。

上述均为临时 CPU 验证；没有在 `data/processed/` 创建正式新队列，也没有运行
服务器模型。旧 `pilot_gh_v1/` 的 101 个输入通过 CPU 解析是历史
结果，不代表新队列已在服务器准备成功，也不证明 H 全部 domain 注释已修复。
已有 B/E 的六组真实特征导出已再次通过校验，继续使用原有格式与入口。
合并测试全部通过 **47 项**：新增准备 13 项、新增运行编排 5 项、原 G/H 准备
11 项、原提取及原生导出校验 18 项。模型子进程在编排测试中使用替身；真实输入
只执行 CPU 解析。这些检查验证接口、来源与流程，不是 GPU 模型推理测试。

服务器完成后，将整个 `features_gh_reconstructed_v1/` 传回本地，保留 `run.json`、
metadata 和 ESM cache 的来源证据。当前尚无 G/H embedding 或 V-measure。
已有 average/PCA 和 V-measure 汇总脚本仍按 B/E 注册，不能仅替换目录就假定能分析 G/H。
特征返回后需扩展其数据集选择与固定队列校验，再沿用原始 cosine＋average-linkage、
PCA＋average-linkage 和同一评分口径；本阶段不增加 Leiden 测试。

H 的性能报告必须带上 `UOC915` 注释已知不完整的限制，不把它解释为完整恢复版本的
最终性能；不要自动删除该样本来提高分数。比较 BiG-SCAPE 时也要核验实际样本覆盖
和坐标版本，公开作者结果仍是外部参考。人工 GCF 是聚类的辅助评价标签，分数不等同
于产品功能或独立泛化性能。
