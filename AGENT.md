# AGENT.md — DAH7PS 变构起源演化动力学重建（V4.1 SOP rev4 执行指南）

> **本仓库的“唯一真相源（Single Source of Truth）”**：`dah7ps_v4_final_sop_4.md`  
> 本文件（AGENT.md）的作用是把 SOP 变成**可执行的、可验收的**任务清单与工程契约，确保自动化执行时不偏航、不“偷懒跳步”、不引入不可追溯的临时决定。

---

## 0. 总原则（必须遵守）

1. **严格按 SOP 的 Phase 顺序执行**：除非某 Phase 明确标注“可选/备选”，否则不得跳过或重排。
2. **所有阈值/模型/随机性必须可追溯**：统一写入 `meta/params.json`；随机抽样必须写明 seed，并保存抽样清单（ids）。
3. **任何 QC 断言失败 = 立刻停止并回溯修复**：不得带病推进。
4. **不覆盖结果**：每次正式跑用独立 run 目录，例如 `results/run_YYYYMMDD/`；或使用 SOP 既定目录，但必须保留旧文件备份（加后缀 `.bak`）。
5. **不手工拼接祖先序列**：AF3 输入的全长祖先序列只能来自 Phase 4.4 的“亚型内全长嵌套 ASR”输出（见 SOP 的 ⚠ [CHECK-02]）。
6. **核心 MSA 禁止插入列膨胀**：核心映射必须走 `Stockholm → 剥离 Insert 列 → AFA`（SOP Phase 3.6）。禁止用 `hmmalign --outformat afa` 直接生成核心比对。
7. **DCA 门槛是硬门槛**：主分析要求 `Meff/L ≥ 3.0`（理想 ≥ 5.0）。不满足则该模块/联合 DCA 禁跑、禁入 ICDC（SOP ⚠ [CHECK-03]）。
8. **祖先 holo 条件必须门控**：对祖先节点必须先 Apo，再口袋/对接验证后才允许做 Holo（SOP ⚠ [CHECK-07]）。

---

## 1. 目录与文件契约

执行过程中必须维持以下目录结构（SOP Phase 0.3）：

- `data/`：原始输入与外部下载（UniProt/NR/UniRef、PDB、ligand、KDOPS 外群等）
- `meta/`：参数、软件版本、外部模型文件（`meta/models/`）
- `scripts/`：所有可执行脚本（必须提供 `--help`、参数检查、可复现输出）
- `results/`：所有中间产物与最终产物（按 SOP 目录分层）

最小必备文件：
- `dah7ps_v4_final_sop_4.md`
- `meta/params.json`
- `results/meta/software_versions.tsv`
- `results/meta/model_files.tsv`（若使用外部模型文件，如 3Di）

---

## 2. 执行总览（按 Phase）

下面的“Done 条件”是验收标准；未满足不得进入下一 Phase。

### Phase 0：环境与可复现性

**任务**
- 创建并激活 conda/mamba 环境（hmmer, mafft, clipkit, iqtree2, foldmason, seqkit, plmc, pastml, gromacs 等）
- 记录软件版本到 `results/meta/software_versions.tsv`
- 下载并锁定外部模型文件（例如 3Di `Q_3Di_models.nex`）到 `meta/models/`，写入 sha256 到 `results/meta/model_files.tsv`

**Done 条件**
- `software_versions.tsv` 存在且包含关键软件版本
- `model_files.tsv` 存在且包含 `Q_3Di_models.nex` 的 sha256（若用 3Di）

---

### Phase 1：数据挖掘（KDOPS 反向过滤 + 全库扫描）

**任务**
- 用 DAH7PS 核心 HMM 与 KDOPS HMM 进行双向评分（KDOPS 反向过滤）
- 生成候选 DAH7PS 序列集合与 domtblout 命中表

**Done 条件**
- `results/01_mining/` 下存在候选序列 fasta + `hits_*.domtbl`
- QC1（长度/去冗余）能读入这些文件

---

### Phase 2：质量控制与去冗余（QC1）

**任务**
- 过滤异常长度序列、去除低复杂度/明显错误序列
- CD-HIT 80% 生成 `nr80_*.fasta`（按亚型/分组可分文件）

**Done 条件**
- `results/02_qc/nr80_*.fasta` 产出
- `results/02_qc/qc_length_report.md` 产出

---

### Phase 3：结构感知核心 MSA（FoldMason 骨架 → HMM 映射）+ 模块注释

#### Phase 3.1–3.5：FoldMason 骨架与核心列定义
**Done 条件**
- `results/03_msa_core/skeleton_core_aa.fa`
- `results/03_msa_core/skeleton_3di.fa`
- `results/03_msa_core/core_columns.mask`

#### Phase 3.6：全量核心映射（关键！）
**任务**
- `extract_core_domains.py` 提取 core-only 序列片段，并输出 `core_domain_coords.tsv`
- `hmmalign` 输出 **Stockholm**，再用 `esl-alimask --rf-is-mask` 剥离 Insert 列，转为 `core_global_matchonly.afa`
- ⚠ 必须启用 Hit Stitching（SOP ⚠ [CHECK-06]）

**Done 条件**
- `results/03_msa_core/core_global_matchonly.afa` 存在
- 核心列数稳定（约 400–600），且显著小于 V3.1（>3000）的膨胀水平
- QC2 报告 `results/03_msa_core/qc_core_alignment.md` 产出

#### Phase 3.7：双版本修剪（树 vs ASR/DCA）
**Done 条件**
- `core_tree.afa` 与 `core_asr.afa` 产出

#### Phase 3.8：模块注释与模块 MSA
**任务**
- `annotate_modules.py` 输出 presence/absence 表
- `extract_module_seqs.py` 必须输出 `${module}_domain_coords.tsv`
- 构建模块 MSA（`ACT_msa.afa` 等）

**Done 条件**
- `module_presence_absence_strict.tsv`
- `ACT_msa.afa` 等模块 MSA 产出

#### Phase 3.9：Profile-anchored Stitching（生成全长亚型 MSA）
**任务**
- 仅对同一架构亚型（例如 Type Iβ-ACT）生成 `msa_full_Ib_v4.afa`
- linker 允许自由对齐；core/module 列必须继承原 MSA，不得漂移

**Done 条件**
- `results/03_msa_core/msa_full_Ib_v4.afa`
- `results/03_msa_full/msa_full_Ib_column_map.tsv`

---

### Phase 4：系统发育与 ASR（QC3 关键：根稳定性）

**任务**
- 外群定根：MFP + 至少一个 site-heterogeneous 模型（LG+C20+F+G / EX_EHO+F+G）
- 生成根稳定性报告 `QC3_root_stability.md`
- 嵌套 ASR（局部子树 + 全长亚型比对）：对 AF3 输入祖先生成连续全长序列

**Done 条件**
- `CoreTree_rooted_MFP.treefile`、`CoreTree_rooted_LGC20.treefile` 至少存在
- `QC3_root_stability.md` 通过（根位置一致或已声明不确定性并给出敏感性方案）
- `ASR_Ib_local.*` 产出（用于 AF3 输入）

---

### Phase 5：结构预测与 MD（QC4 关键：门控）

**任务**
- 每个祖先节点先做 Apo AF3
- 口袋检测 + 对接门控后才允许 holo
- MD 前必须过结构 QC（SOP ⚠ [CHECK-05]）
- 运行 2×2 因子矩阵（holo 条件可能被门控取消）

**Done 条件**
- `results/05_struct_md/pocket_gating/<node>_gating.md`
- `results/05_struct_md/qc_dynamics.md`

---

### Phase 6：DCA × 动力学融合（ICDC）

**任务**
- 核心 DCA（`core_dca.afa`）
- 模块 DCA（Meff/L ≥ 3 才跑）
- 联合 DCA（`Ib_ACT_joint_dca.afa`）提取跨域耦联 `Ib_ACT_crossdomain_top200.tsv`
- with/without 模块比较：必须 Meff 匹配下采样 + Z-score（SOP 6.2.5）
- ICDC 融合（坐标映射必须先通过 SOP ⚠ [CHECK-04]）

**Done 条件**
- `results/06_icdc/icdc_crosslayer_paths.tsv`
- `results/06_icdc/icdc_crossdomain_couplings.tsv`（若跑联合 DCA）
- `results/06_icdc/coordinate_map.tsv`（断言通过）

---

## 3. 脚本最小接口约定（必须实现）

若仓库中缺少 SOP 引用脚本，AGENT 必须按以下“接口契约”实现（保持可测试、可复现）：

- 所有脚本必须支持：
  - `--help`
  - 输入文件存在性检查
  - 输出目录自动创建
  - 失败时非 0 退出码
- 关键脚本与输出：
  - `extract_core_domains.py` → `all_core_only.fasta` + `core_domain_coords.tsv`（含 hit stitching 日志）
  - `prepare_dca_input.py` → `*_dca.afa` + `*_dca_stats.tsv`（含 Meff/L）
  - `stitch_full_length_msa.py` → `msa_full_*.afa` + `column_map.tsv`（断言 core 列不变）
  - `qc_root_stability.py` → `QC3_root_stability.md`
  - `coordinate_mapper.py` → `coordinate_map.tsv`（断言一致）

---

## 4. 常见故障与处理策略（快速索引）

- **核心 MSA 列数突然暴涨（>1000）**：99% 是 Insert 列未剥离或误用了 `hmmalign --outformat afa`。回到 Phase 3.6 修复。
- **Type II 序列核心被切半**：未做 hit stitching（SOP ⚠ [CHECK-06]）。
- **根位置在 MFP 与 LGC20 不一致**：典型 LBA 风险；必须声明根不确定性，并在两种根假设下重复 PastML/关键结论敏感性分析（SOP QC3）。
- **模块 DCA 看起来“很漂亮”但 Meff/L < 3**：一律视为噪声，禁入 ICDC（SOP ⚠ [CHECK-03]）。
- **祖先 holo 结构出现“神奇新口袋”**：未做 Apo-first 门控（SOP ⚠ [CHECK-07]），撤回 holo 条件，先做口袋/对接验证。
- **坐标映射对不上**：优先检查 PDB insertion code、ClipKIT 修剪导致列号变化、GROMACS 拓扑编号偏移（SOP ⚠ [CHECK-04]）。

---

## 5. 交付物清单（最终验收）

- `results/03_msa_core/qc_core_alignment.md`
- `results/04_phylogeny_asr/QC3_root_stability.md`
- `results/05_struct_md/qc_dynamics.md`
- `results/06_icdc/icdc_core_network.graphml`
- `results/06_icdc/icdc_module_network.graphml`
- `results/06_icdc/icdc_crosslayer_paths.tsv`
- （若跑联合 DCA）`results/06_icdc/icdc_crossdomain_couplings.tsv`

