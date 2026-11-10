# B、E 数据集的服务器特征提取

本阶段为 B（Divergent，39 个文件）和 E（MIBiG_BS1，281 个文件）分别导出
ESM2、BGC-MAC 和完整 DNA 的冻结特征。每个 BGC 对应一行向量；两个数据集
分别保留，不合并 GCF 标签。本脚本不执行聚类、BiG-SCAPE 或模型训练。

入口为 `scripts/pilot_features.py`，依次使用 `prepare`、`extract`、`validate`。
下面命令从服务器的 `BIG-SCAPE` 根目录运行；请将 `/server/...` 换成实际路径。

## 服务器需要的文件与环境

同步以下内容：

- 本项目的 `scripts/` 和 `data/processed/inventory_v1/`。
- 与本次脚本匹配的 `BGC_annotation` 源码，包含
  `experiment/extract_features.py`、`experiment/extract_dna_features.py`
  及它们依赖的 `data_preparation/`、`model/`。
- B、E 原始数据。服务器不必复制其他七个数据集，但 `--raw-root` 下必须保留
  `Dataset_B_Divergent/`、`Dataset_E_MIBiG_BS1/` 两层数据集目录及其内部结构。
  该参数指向包含这两个目录的 `curated_datasets`，不是其中的 `GCF_assignments`。
- 本地 ESM2 `.pt` 权重、训练好的 BGC-MAC `.ckpt`、DNA 模型的 Hugging Face
  checkpoint 目录。仅复制配置文件不够；脚本不会自动下载权重。

在服务器已匹配 CUDA 的 PyTorch 环境中安装项目依赖：

```bash
python -m pip install -e '/server/BGC_annotation[dna]'
python -m pip install fair-esm rotary-embedding-torch
```

当前 `BGC_annotation/setup.py` 的 `dna` extra 包含 `torch>=2.4`、
`transformers==5.1.0`、`numpy` 和 `biopython`；蛋白支路还需 `fair-esm`、
`rotary-embedding-torch`。DNA 默认使用 Flash Attention 2，它不在上述 extra
中，需要服务器已有与 CUDA、PyTorch 匹配的 `flash-attn` 安装。沿用已验证的
服务器环境，避免为此替换 GPU 驱动或任意升级 PyTorch。

默认使用启动脚本的同一个 Python 解释器调用模型。若模型环境独立，可指定
`--python /server/env/bin/python`（使用绝对路径）；运行本脚本的环境仍需 NumPy
完成导出校验。`prepare` 本身仅使用 Python 标准库，不加载模型。

ESM2 checkpoint 应与 MAC 训练时使用的 ESM2 对应；常见配置为
`esm2_t33_650M_UR50D.pt`。MAC 选择一个 ensemble member，默认索引 `0`。
DNA 权重目录默认同时提供 config 和 tokenizer；若二者单独存放，使用下文的
`--dna-config` 和 `--dna-tokenizer`，并确保它们来自同一模型。

## 1. 准备固定样本队列

```bash
python scripts/pilot_features.py prepare \
  --raw-root /server/curated_datasets \
  --inventory data/processed/inventory_v1 \
  --output data/processed/pilot_be_v1
```

默认 `--datasets B E`，也可先使用 `--datasets B` 准备较小的试跑队列。
准备步骤使用库存清单中的稳定 `sample_id`，核对服务器原始文件与清单一致，
将 GBK 复制到新队列并以 `sample_id` 命名。文件内容保持不变，原始 ID、标签
和路径保留在对应清单中。这样模型读取的文件名不会包含 GCF 标签，也不会因
不同数据集的同名文件发生 ID 冲突。

保留全部 B、E 条目，包括 singleton；没有长度过滤、标签修正、序列修复或
去重删除。B 中含不换行空格的原始标签也保留。原始数据仍在 `raw`，准备结果
和后续特征均放在 `processed`。如果服务器原始文件与清单不一致，应核查数据
版本，不能跳过来源校验来继续提取。

准备结果的结构为：

```text
data/processed/pilot_be_v1/
  cohort.tsv                          # sample_id ↔ 原始 ID、路径、SHA256、GCF
  metadata.json                       # 来源、队列规则和清单校验和
  Dataset_B_Divergent/
    inputs/bgc_<stable_hash>.gbk       # 每个样本一个内容不变的副本
    curated_GCF_assignments.tsv        # BGC 列改为 sample_id，GCF 保留原值
  Dataset_E_MIBiG_BS1/
    inputs/bgc_<stable_hash>.gbk
    curated_GCF_assignments.tsv
```

## 2. 提取三个模型的特征

下面的 `--trusted-input` 仅适用于你确认来源可信的 ESM2/MAC 权重；它允许
读取旧 checkpoint 中的 Python pickle 对象。若不需要此兼容，可去掉该参数；
受限加载失败时不会自动切换为不受限加载。

```bash
python scripts/pilot_features.py extract \
  --prepared data/processed/pilot_be_v1 \
  --encoder-project /server/BGC_annotation \
  --output data/processed/features_be_v1 \
  --esm-checkpoint /server/models/esm2_t33_650M_UR50D.pt \
  --mac-checkpoint /server/models/MAC_default.ckpt \
  --dna-checkpoint /server/models/dna_checkpoint \
  --device cuda:0 \
  --esm-batch-size 8 \
  --batch-size 16 \
  --mac-member 0 \
  --trusted-input
```

默认运行 `--routes esm2 mac dna_full`；数据集由 `--prepared` 的队列决定，
`extract` 不另接收 `--datasets`。三条路线的输入必须来自这份固定队列。

| 路线 | 计算方式 | 一般维度 |
| --- | --- | ---: |
| `esm2` | ESM2 domain/蛋白 token 向量 → BGC 内 masked mean → L2 归一化 | 1280 |
| `mac` | 同一份 ESM2 token 向量 → 冻结的训练后 MAC BGCEncoder → masked mean → L2 | 256 |
| `dna_full` | 完整 BGC DNA 一次 backbone forward → 内容 BPE token mean → L2 | 768 |

实际维度以加载的 checkpoint 和输出 metadata 为准，不按表中常见值强行改变。

蛋白支路沿用现有解析：有 `aSDomain` 的 CDS 使用其 domain，没有 domain 的
CDS 使用完整翻译。ESM2 固定使用 BOS 读出；超过 1021 aa 的蛋白/domain
按不重叠窗口完整覆盖，再按窗口残基数加权。`--esm-batch-size` 控制 ESM2
窗口批量，`--batch-size` 控制 BGC 聚合批量，两者可以分别调整。

MAC 复用本次 ESM2 阶段新生成的 `esm2/tokens.pth`，不重新计算 ESM2，也不
接受来源不明的外部旧缓存。即使只指定 `--routes mac`，也会自动先生成所需的
ESM2 阶段。MAC 使用实际训练权重，缺失或不兼容时会报错。

DNA 固定为 **full 模式、每次一条完整 BGC**；默认使用 **bfloat16 和
Flash Attention 2**。
不设置 bp/token 长度上限，不按旧文档中的 40 kb 规则筛选；不截断，也不在
OOM 后自动切成 windowed。保留原始方向及模糊碱基，pooling 不包含特殊 token
和 padding。任意长度输入接口不代表任意长度都能装入显存；资源失败会明确
记录并使阶段失败，不把失败条目静默排除。

必要时可显式设置 `--dna-dtype bfloat16|float16|float32` 和
`--dna-attention-backend flash_attention_2|sdpa|eager`。FA2 不支持 float32；
改用其他后端也不会改变 full 输入策略。长序列 SDPA/eager 可能占用更多显存，
不同精度/后端应使用新输出目录并保留配置，不混作同一轮结果。

如 DNA config/tokenizer 与权重分开，向提取命令追加：

```bash
  --dna-tokenizer /server/models/dna_tokenizer \
  --dna-config /server/models/dna_config/config.json
```

也可用新输出目录只运行 DNA：

```bash
python scripts/pilot_features.py extract \
  --prepared data/processed/pilot_be_v1 \
  --encoder-project /server/BGC_annotation \
  --output data/processed/features_be_dna_v1 \
  --routes dna_full \
  --dna-checkpoint /server/models/dna_checkpoint \
  --device cuda:0
```

## 3. 输出与校验

默认全部路线完成后的主要目录为：

```text
data/processed/features_be_v1/
  run.json
  Dataset_B_Divergent/
    esm2/                 # 含本次生成的 tokens.pth
    mac/
    dna_full/
  Dataset_E_MIBiG_BS1/
    esm2/
    mac/
    dna_full/
```

每个路线目录保留原生提取器的四类产物：

- `embeddings.npy`：float32 二维矩阵，每行一个 L2 归一化的 BGC 向量。
- `bgc_ids.tsv`：`row` 对应矩阵行号，`bgc_id` 对应准备阶段的稳定 `sample_id`。
  `token_count` 的单位在蛋白和 DNA 支路不同，不应混用。
- `metadata.json`：模型来源、维度、pooling、源代码/输入哈希、运行参数等。
- `failures.tsv`：明确的失败记录；完整成功时仍可保留只有表头的文件。

每个路线还有 `execution.log`，记录原生提取器的标准输出和错误；运行时可查看
该文件了解进度。顶层 `run.json` 保存实际命令、各阶段状态、耗时和验证结果。

原始 BGC 名称和 GCF 通过准备阶段的 `cohort.tsv` 关联。不要按目录遍历顺序或标签行号
去对齐向量，也不要把一个数据集的家族名称当作另一个数据集的同义标签。

提取后运行：

```bash
python scripts/pilot_features.py validate \
  --prepared data/processed/pilot_be_v1 \
  --output data/processed/features_be_v1
```

校验范围由 `run.json` 中请求的阶段决定，包括自动加入的依赖阶段。完整队列
应分别得到 B 的 39 行和 E 的 281 行；检查对应 ID、数组形状、有限数值、L2
归一化、成功状态和来源信息。只成功一部分的结果不能作为完整 B/E 特征集，
也不能靠取三个模型成功 ID 的交集来掩盖缺失。

## 中断与重新运行

重复同一条 `extract` 命令并追加 `--resume`，可以复用已完成且通过当前配置与
产物校验的阶段。中断或失败阶段会归档到 `_failed`，然后从头重跑该阶段。
这不是逐样本续跑：例如 DNA 在第 100 条失败，重试该 DNA 阶段仍从第一条开始。
ESM2/MAC 已完整完成的阶段通过校验后可以跳过。

恢复时保留原命令的输入、模型、路线和参数。改变 checkpoint、MAC member、
设备、批量或其他运行配置，应选择新的 `--output`，防止将不同配置混在一起。
不要手动复制 `metadata.json` 或零散向量文件来伪造一个完成阶段。

脚本使用 `.running.lock` 阻止同一输出目录被同时恢复。如果进程被强制终止或
服务器重启而留下锁，先确认对应任务和模型子进程都已停止，再移除这个残留锁，
之后使用原命令加 `--resume`。正常完成或捕获到异常时脚本会自行清理锁。

## 结果边界

本地脚本与 CPU 检查不代替服务器真实权重、CUDA 和 Flash Attention 2 的运行
验证。第一次运行可先准备 B，再扩展到 B/E 的完整队列。

E 来自 MIBiG，BGC-MAC 的训练数据可能与之重叠；冻结模型不能消除训练暴露。
B 与 E 也存在跨数据集序列重复。因此这轮适合跑通特征导出和对照流程；如果
后续用 B 调参、E 测试或进行 fine-tuning，需要另行审计重叠，不能直接将结果
称为独立泛化测试。这里只导出向量，不据此判定哪种模型或聚类方法更优。
