# H 既有 Geneious aSDomain 恢复

本步骤恢复 H 中被 Geneious 保存成 `misc_feature` 的已有 aSDomain，使其能被
当前蛋白输入解析器识别。沿用 aSDomain 分割策略，**暂不进行 Pfam 坐标分割**。
入口为 [`scripts/reconstruct_h.py`](../scripts/reconstruct_h.py)，核心证据检查与
文本恢复由 [`scripts/geneious_asdomain_restore.py`](../scripts/geneious_asdomain_restore.py) 实现。

## 为什么还需要恢复

此前的 [`repaired_gh_v1/`](../data/processed/repaired_gh_v1/summary.md) 修正了原生
domain 的端点错误，使 G/H 全部 101 个输入通过 CPU 解析。但 H 另有 2365 条
`misc_feature` 带有确切的 `/note="Geneious type: aSDomain"`，分布在 45/90 个 BGC。
这些条目缺少 `translation`、`protein_start`、`protein_end`，现有解析器只使用
类型为 `aSDomain` 的 feature，因此没有把它们用于蛋白分割。

旧 `pilot_gh_v1/` 中 H 的 156 条超长 token 均为整条 CDS，涉及 42 个 BGC，
最长 9508 aa；这些 CDS 均有尚未恢复的 Geneious aSDomain。
因此旧统计不能解释成“按已有 aSDomain 分割后仍超过 ESM2 上限”。
旧队列的 `complete` 仅描述当时的来源校验与解析准备已完成，不能证明 domain
注释完整。保留该历史快照，暂勿据此启动新的 H 三路线比较。

## 运行与数据边界

在项目根目录、已安装 `requirements-inventory.txt` 的环境中运行：

```bash
python scripts/reconstruct_h.py \
  --source data/processed/repaired_gh_v1 \
  --encoder-project ../BGC_annotation \
  --output data/processed/reconstructed_h_v1
```

默认来源是 `data/processed/repaired_gh_v1/`，保留之前已有序列证据的坐标修复。
输出是 H 的全部 90 个 GBK 副本及标签，仍使用原始文件名与家族目录。
G、`data/raw/`、原始清单和已有 processed 快照均不改写。
重跑时使用新的输出版本，不覆盖旧目录。这里只进行注释恢复和 CPU 输入验证，
不加载 ESM2、BGC-MAC 或 DNA 模型，也不运行聚类。
`--workspace`、`--raw`、`--inventory` 可指定工作区和追溯来源路径。
在服务器使用已生成副本时，同步完整恢复快照即可；只有重跑本恢复命令才需要
完整 raw、inventory 与坐标修复快照。

## 恢复依据

只考虑原类型为 `misc_feature`、且具有确切 Geneious aSDomain 标记的条目，
不把 `PFAM_domain` 或其他普通 `misc_feature` 转成 aSDomain。逐条核验：

1. 原核酸区间对应唯一、同链且包含它的父 CDS。
2. 区间相对父 CDS 的位置满足链方向、`codon_start` 和完整密码子的阅读框。
3. domain 核酸翻译与父 CDS 对应蛋白区间完全一致。
4. 生成的 `protein_start` / `protein_end` 为从零开始、左闭右开的蛋白坐标，
   与基因组位置和恢复出的 `translation` 一致。

通过后恢复 feature 类型、翻译和蛋白坐标；如旧 `locus_tag` 与证据支持的父 CDS
不一致，必须记录旧值、新父 ID 和匹配证据。不得只因 ID 相似就选择父 CDS。
缺失的 domain 家族名、HMM 分数和 Pfam ID 不推测填写；保留原有 `standard_name`
等注释。序列一致证明旧区间可按当前 CDS 解释，不等同于重新验证 HMM 的生物学边界。

实际运行中 2272 条直接满足证据，另 91 条通过有证据的旧父 ID 映射恢复。
`UOC915.gbk` 的两条注释共享 1154-bp 区间，且相对父 CDS 的起点不在阅读框内。
原文件 `COMMENT` 记录在位置 17584–17585 插入过 `GG`，落在该旧 domain 区间内；
当前 CDS 的 DNA 与翻译相符，但旧 domain 边界不符合当前阅读框。九个原始数据集
中也没有可用的完全相同父蛋白及可信 domain 边界，无法据此转移坐标。
这两条保留原样并列入复核，不能四舍五入坐标、统一加一个密码子或复制同名
domain 的位置。
逐条证据保存在输出的 `domain_checks.tsv` 和 `review_items.tsv`。

此前 [坐标修复](repair_gh.md) 能证明“末端少一个密码子”，是因为它有独立的旧
domain 翻译证据。本批 Geneious 注释没有这份翻译，不能直接套用该修复规则。
证据不足的条目保持原样，列为 `needs_review`，不删除其所在样本。

## 保留规则与追溯

- 保留全部 H 样本和 GCF 标签，不按 DNA 长度过滤，不重分配家族。
- 保留所有原有重复注释。实际记录 187 个相同父 CDS/位置的重复组、190 个额外
  条目，涉及 40 个 BGC；本次只报告，不去重、不调整其输入权重。
- DNA、CDS 翻译、record ID、已有原生 aSDomain/PFAM 注释以及其他 feature 不改写。
  使用最小文本变更，不整体重新序列化 GenBank。
- 保存逐条候选的原始类型、位置、父 ID、证据、动作和待复核原因；同时保留
  样本身份、源文件与输出文件哈希，以及 CPU 解析前后的 token/长度检查。
- 若仍有未解决条目，可发布包含全部 90 个样本、状态为 `needs_review` 的可检查
  副本，但不将其称为完整的 benchmark 输入。

## 实际结果

2026-11-19 已发布 [`reconstructed_h_v1/summary.md`](../data/processed/reconstructed_h_v1/summary.md)。
状态为 **`needs_review`**，全部 90 个 H 样本 / 10 个 GCF 保留。

| 检查 | 结果 |
| --- | ---: |
| 修改的 GBK | 45 / 90 |
| Geneious aSDomain 候选 | 2365 |
| 成功恢复 | 2363 |
| 待复核、保留原样 | 2，均在 `UOC915` |
| 有证据的父 `locus_tag` 修正 | 91 |
| 保留的重复组 / 额外条目 | 187 / 190 |
| 蛋白与 DNA 解析通过 | 均为 90 / 90 |
| aSDomain 恢复检查通过的样本 | 89 / 90 |
| token 总数，恢复前 → 后 | 4836 → 6791 |
| 超过 1021 aa 的 token，恢复前 → 后 | 156 → 0 |
| 最长 token，恢复前 → 后 | 9508 → 1008 aa |

恢复后 6791 个 token 中，4718 个来自识别到的 aSDomain，2073 个为没有识别到
aSDomain 的整条蛋白。其余 45 个 GBK 与来源副本逐字节一致。

**没有超长 token 不代表 UOC915 的两条注释已经解决。** 现有解析器仍会跳过这两条
`misc_feature`；完整 H 的 domain 输入仍存在已知缺失，不能仅凭解析通过或长度
检查通过就宣布 benchmark 就绪。另有坐标修复阶段遗留的两条 `ApSA` 无父 CDS
PFAM 注释，原样保留在 `inherited_review_items.tsv`，不计入本次 2365 条候选。

全部输出 GBK 已复读，验证恢复幂等、非目标 feature 保留和预期的 token 变化；
另已独立核验全部 106 个产物文件与当前脚本哈希。
本次 16 项核心测试与 6 项编排测试通过。全部 2941 个 raw 文件和来源快照的哈希
前后一致。未运行模型、未生成新 prepared 队列，也没有新的 G/H 聚类分数。

## 输出文件与状态

```text
data/processed/reconstructed_h_v1/
  Dataset_H_Reclassification/
    curated_GCF_assignments.tsv
    GCF_assignments/<family>/<original_name>.gbk
  samples.tsv
  reconstructions.tsv
  domain_checks.tsv
  review_items.tsv
  parent_id_mappings.tsv
  duplicate_annotations.tsv
  inherited_review_items.tsv
  raw_integrity.tsv
  source_integrity.tsv
  source_provenance/
    metadata.json
    samples.tsv
    repairs.tsv
    domain_checks.tsv
    review_items.tsv
  metadata.json
  summary.md
```

| 文件 | 用途 |
| --- | --- |
| `samples.tsv` | 全部 H 样本身份、标签、输入/输出哈希及解析前后统计 |
| `reconstructions.tsv` | 实际恢复条目和修改证据 |
| `domain_checks.tsv` | 全部候选的父 CDS、阅读框、核酸翻译与蛋白区间检查 |
| `review_items.tsv` | 本次未恢复、保持原样的条目及原因 |
| `parent_id_mappings.tsv` | 旧 `locus_tag` 到证据支持父 CDS 的映射 |
| `duplicate_annotations.tsv` | 保留的重复注释及其分组 |
| `inherited_review_items.tsv` | 之前坐标修复快照中保留的待核查事项 |
| `raw_integrity.tsv` / `source_integrity.tsv` | raw 与来源快照未改变的校验 |
| `source_provenance/` | 嵌入坐标修复快照的来源记录 |
| `metadata.json` / `summary.md` | 参数、策略、代码/产物哈希、完成状态与摘要 |

输出 schema 为 `bgc_h_asdomain_reconstruction_v1`。如果任一目标 aSDomain 尚未恢复，
`metadata.status` 为 `needs_review`，`ready_for_full_h_encoding` 为 `false`。
`snapshot_complete: true` 只表示全部副本与审计文件已写完并通过快照检查，
不表示注释全部恢复。后续服务器桥接默认保留全部 90 个 H 样本并携带
已知不完整状态；运行成功不能将该状态改标为完整恢复。

## 后续编码

恢复快照保留原始目录结构，**不是** `pilot_features.py extract` 的 `--prepared`
目录。现有 `prepare-repaired` 面向坐标修复快照格式，不能把新恢复快照直接塞入
该入口，或把新 GBK 覆盖到 `pilot_gh_v1/` 而继续使用旧哈希。
新的 `pilot_features.py prepare-reconstructed` 接受本恢复 schema，核验审计与
待复核状态，再生成带有完整恢复来源的新稳定 ID 队列。它可合入坐标修复后的 G，
也支持单独 H；`prepare-repaired` 保持原有坐标修复输入用途。

从服务器项目根目录运行：

```bash
python scripts/pilot_features.py prepare-reconstructed \
  --reconstructed-root data/processed/reconstructed_h_v1 \
  --repaired-root data/processed/repaired_gh_v1 \
  --datasets G H \
  --encoder-project /home/s252613/BGC_encoder/BGC_annotation \
  --output data/processed/pilot_gh_reconstructed_v1
```

准备时保留 H 全部 90 个样本和两条未恢复注释，复制 GBK 字节而不再次改注释。
新队列 metadata、成员表和后续运行记录保留 `UOC915` 的 `needs_review` 状态。
队列 `status: complete` 只表示机械准备完成，与原始恢复快照
的 `ready_for_full_h_encoding: false` 并不矛盾；后者仍说明 H 有已知注释缺失。
prepared schema 为 `bgc_reconstructed_pilot_inputs_v1`，提取运行 schema 为
`bgc_reconstructed_feature_run_v1`；`annotation_complete: false`、两条未恢复注释、
对应稳定 ID 和 `review_policy: retain` 随队列与运行记录保留。
本地在 `/tmp` 用真实 G/H 全部 101 个样本完成了 CPU 准备与加载验证，H-only
90 个样本也通过；未加载模型，未改变真实输入哈希。尚未发布正式新队列，
由用户在服务器执行上述命令；验证范围见 [G/H 操作说明](gh_features.md)。

新的输入须重新统计 domain token、残留超长 token 和 ESM2 窗口数，使用新队列
哈希重建 ESM2 token cache，再由 BGC-MAC 读取同批 cache。DNA 序列应保持不变，
但新运行仍应记录与实际输入匹配的来源证据。模型配置继续参考
[G/H 特征提取说明](gh_features.md)，本步骤不宣称已有新的模型结果或 V-measure。
