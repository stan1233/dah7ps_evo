# AGENT.md — DAH7PS 变构起源演化动力学重建（V5.0 SOP rev5 执行指南）

> **本仓库的"唯一真相源（Single Source of Truth）"**：`PLAN.md`（V5.0 SOP rev5）
> 本文件（AGENT.md）的作用是把 SOP 变成**可执行的、可验收的**任务清单与工程契约，确保自动化执行时不偏航、不"偷懒跳步"、不引入不可追溯的临时决定。

---

## 0. 总原则（必须遵守）

1. **严格按 SOP 的 Phase 顺序执行**：除非某 Phase 明确标注"可选/备选"，否则不得跳过或重排。
2. **所有阈值/模型/随机性必须可追溯**：统一写入 `meta/params.json`；随机抽样必须写明 seed，并保存抽样清单（ids）。
3. **任何 QC 断言失败 = 立刻停止并回溯修复**：不得带病推进。
4. **不覆盖结果**：每次正式跑用独立 run 目录，例如 `results/run_YYYYMMDD/`；或使用 SOP 既定目录，但必须保留旧文件备份（加后缀 `.bak`）。
5. **不手工拼接祖先序列**：AF3 输入的全长祖先序列只能来自 Phase 4.4 的"亚型内全长嵌套 ASR"输出（见 SOP ⚠ [CHECK-02]）。
6. **核心 MSA 禁止插入列膨胀**：核心映射必须走 `Stockholm → 剥离 Insert 列 → AFA`（SOP Phase 3.6）。禁止用 `hmmalign --outformat afa` 直接生成核心比对。
7. **DCA 门槛是硬门槛**：主分析（核心层）要求 `Meff/L ≥ 3.0`（理想 ≥ 5.0）。模块层和联合跨域 DCA 因深度不足已移入"可选探索"，不作为论文主线结论（SOP ⚠ [CHECK-03]；V5.0 Phase 6）。
8. **装配体状态不可默认（V5.0 硬约束）**：不得将 oligomeric state 写死为"四聚体"。每个祖先节点的功能装配体必须通过 Phase 5.0 Assembly Adjudication 显式判定（文献扫描 + dimer/tetramer 平行 AF3 + PISA 评分），产出 `assembly_adjudication.tsv`。所有下游结构预测与 MD 的拷贝数读取该表（SOP ⚠ [CHECK-08]）。
9. **实验记录必须同步更新**：每次执行任何分析操作，必须将精确命令、参数、输出文件、结果摘要和时间戳记录到 `log.md`。禁止"先跑完再补记录"。
10. **CPU 核心全利用**：本机 28 核，工具设为 **20 线程**。`mafft --thread 20`（禁止 `--thread -1`）、`hmmsearch/hmmalign --cpu 20`、`iqtree -T 20`、`mmseqs --threads 20`、`cd-hit -T 20`、Python `max_workers=20`。
11. **树–比对 tip 集严格一致（V5.0 硬约束）**：IQ-TREE `-te` 要求树与比对的 tip 集完全匹配。Phase 4.3 ASR 前必须从定根树中 prune 掉 KDOPS 外群 tips → `CoreTree_rooted_ingroup.treefile`，并用 `assert_tip_match.py` 断言。不匹配 = 终止。
12. **V5.0 全部 Apo-only**：Phase 5 不做 Holo 预测，消除"现代配体幻觉"风险。V4.1 的 CHECK-07 已不再适用。
13. **ICDC 仅为 Discussion 展望**：不作为论文主线定量结论。核心 DCA × 祖先结构 × 有限 MD 的一致性仅作定性描述（SOP Phase 7.3）。

---

## 1. 目录与文件契约

执行过程中必须维持以下目录结构（SOP Phase 0.3）：

- `data/`：原始输入与外部下载
- `meta/`：参数、软件版本、外部模型文件（`meta/models/`）
- `scripts/`：所有可执行脚本（必须提供 `--help`、参数检查、可复现输出）
- `results/`：所有中间产物与最终产物
  - `results/03_msa_core/` — 核心 MSA 及其 QC
  - `results/03_msa_modules/` — 模块注释矩阵、模块序列、模块 MSA
  - `results/03_msa_full/` — 亚型内全长缝合 MSA + linker + column_map ⚠ 全长 MSA 禁止放入 `03_msa_core/`
  - `results/04_phylogeny_asr/` — 核心树、ASR 产出、PastML 结果
  - `results/05_struct_valid/` — 装配体判定、AF3 结构、MD 轨迹（V5.0 仅 Apo 验证性）
  - `results/06_dca/` — 核心层 DCA + 可选探索性 DCA
  - **禁止**跨目录混放

**ACT 低 prevalence 策略决策（2026-03-03）：** ACT strict = 47 seqs（L=142），Meff/L ≈ 0.2–0.3。ACT DCA 排除出主线证据链，仅作探索性附录。

最小必备文件：`PLAN.md`（V5.0）、`meta/params.json`、`results/meta/software_versions.tsv`、`results/meta/model_files.tsv`

---

## 2. 执行总览（按 Phase）

下面的"Done 条件"是验收标准；未满足不得进入下一 Phase。

### Phase 0：环境与可复现性 ✅

**Done 条件** — `software_versions.tsv` + `params.json` + `model_files.tsv` 均存在

### Phase 1：数据挖掘 ✅

**Done 条件** — `results/01_mining/` 候选序列 fasta + `hits_*.domtbl` + QC1 报告

### Phase 2：质量控制与去冗余 ✅

**Done 条件** — `results/02_qc/nr80_*.fasta` + `qc_length_report.md` + seeds60 + stepping stones

### Phase 3：结构感知核心 MSA + 模块注释

#### Phase 3.1–3.8 ✅

**Done 条件**
- `panel_candidates.tsv` + `panel_manifest.tsv`
- `skeleton_core_aa.fa` + `core_columns.mask`
- `core_global_matchonly.afa`（9,393 × 521）
- `core_tree.afa`（436 cols）+ `core_asr.afa`（472 cols）
- `module_presence_absence_strict.tsv` + 5 模块 MSA
- QC2 报告

#### Phase 3.9：Profile-anchored Stitching [待执行]

**任务**
- 仅对同一架构亚型（如 Type Iβ-ACT）生成全长缝合 MSA
- linker 允许自由对齐；core/module 列必须继承原 MSA，不得漂移

**Done 条件**
- `results/03_msa_full/msa_full_Ib_v4.afa`（⚠ 路径在 `03_msa_full/` 而非 `03_msa_core/`）
- `results/03_msa_full/msa_full_Ib_column_map.tsv`
- Core 段列数与 `core_asr.afa` 完全一致（断言通过）

---

### Phase 4：系统发育与分层 ASR [待执行]

**任务**
- 外群定根：MFP + 至少一个 site-heterogeneous 模型（LG+C20+F+G / EX_EHO+F+G）
- **⚠ Phase 4.3 树–比对 tip 集一致性（V5.0 硬约束）**：定根树包含 KDOPS 外群，但 `core_asr.afa` 不包含。ASR 前必须 prune 外群 → `CoreTree_rooted_ingroup.treefile`，并用 `assert_tip_match.py` 断言 tip 集严格一致。选择 prune 方案（而非合并外群到 ASR 比对），因为外群序列在 ASR 中只会引入噪声。
- 嵌套 ASR（局部子树 + 全长亚型比对，输入来自 `results/03_msa_full/`）

**Done 条件**
- `CoreTree_rooted_MFP.treefile` + `CoreTree_rooted_LGC20.treefile` 至少存在
- `CoreTree_rooted_ingroup.treefile`（pruned，用于 ASR）
- `QC3_root_stability.md` 通过
- `ASR_Ib_local.*` 产出

---

### Phase 5：关键祖先节点的结构验证 [待执行，V5.0 精简版]

> **V5.0 范围：** 2–4 个关键节点 × Apo-only × native-oligomer（由 5.0 判定）× 有限验证性 MD。**不做 Holo。不默认四聚体。**

**任务**
- **5.0 装配体判定（Assembly Adjudication）[CHECK-08]**：对每个祖先节点，先做文献/结构注释扫描 → dimer + tetramer 平行 AF3 → PISA 界面评分 → 产出 `assembly_adjudication.tsv`。后续步骤的拷贝数全部读取该表。
- 5.1 候选祖先节点选择（Pre-gain / Post-gain，2–4 个）
- 5.2 Apo 结构预测（AF3 拷贝数由 5.0 判定 + ESMFold 交叉验证）
- 5.2b 结构 QC 门控 [CHECK-05]：ipTM ≥ 0.6 → 能量最小化 → 10 ns 筛选
- 5.3 验证性 MD：native-oligomer-apo ≥200 ns × 2 rep（主）+ 可选 alternative-oligomer 10–20 ns 排除测试
- 5.4 有限动力学读出（RMSD/RMSF + 界面面积 + Fpocket 口袋拓扑比较，不做完整 DCCM/网络流）

**Done 条件**
- `results/05_struct_valid/assembly_adjudication.tsv`（每个节点有装配体判定与证据链）
- `results/05_struct_valid/qc_struct_validation.md`（含装配体判定证据、AF3 指标、MD 收敛性、界面稳定性）

---

### Phase 6：核心层共进化分析（DCA） [待执行，V5.0 聚焦版]

> **V5.0 范围：** 仅核心层 DCA 进入主线。模块 DCA、联合跨域 DCA → 可选探索，不入论文主线结论。

**任务**
- 6.1 核心层 DCA 输入准备（`core_asr.afa` → gap 过滤 → `core_dca.afa`，Meff/L 门控）
- 6.2 plmc 执行
- 6.3 显著性评估：top-L 接触验证（1KFL/1RZM/3NV8）+ 功能位点富集 + 跨亚型保守 vs 特异耦联
- 6.4 可选探索（不入主线）：模块层 DCA、联合跨域 DCA、Meff 匹配比较

**Done 条件**
- `results/06_dca/core_dca.afa` + `core_dca_stats.tsv`（Meff/L ≥ 5）
- `results/06_dca/core_significant_couplings.tsv`
- `results/06_dca/qc_core_dca.md`

---

### Phase 7：论文写作蓝图与跨证据一致性展望 [待执行，V5.0 新增]

**任务**
- 7.1 论文主线叙事锁定（结构感知 MSA → 定根树 + 模块 ASR → 核心 DCA → 祖先结构验证）
- 7.2 核心图表清单（Fig 1–6 + 补充图）
- 7.3 ICDC 在 Discussion 中的定位：倒数第二段"跨证据一致性展望"，定性展示核心 DCA × 祖先结构 × 有限 MD 的一致性信号，明确声明非正式 ICDC 融合，提出全因子 MD → DCCM → 正式融合作为未来方向

**Done 条件**
- 论文主线叙事文档 + Fig 1–6 齐备
- ICDC 定位写法确定

---

## 3. 脚本最小接口约定（必须实现）

若仓库中缺少 SOP 引用脚本，AGENT 必须按以下"接口契约"实现：

- 所有脚本必须支持 `--help`、输入文件存在性检查、输出目录自动创建、失败时非 0 退出码
- 关键脚本与输出：
  - `extract_core_domains.py` → `all_core_only.fasta` + `core_domain_coords.tsv`（含 hit stitching） ✅
  - `prepare_dca_input.py` → `*_dca.afa` + `*_dca_stats.tsv`（含 Meff/L）
  - `stitch_full_length_msa.py` → `results/03_msa_full/msa_full_*.afa` + `column_map.tsv`（断言 core 列不变）
  - `prune_tree.py` → `CoreTree_rooted_ingroup.treefile`（V5.0 新增）
  - `assert_tip_match.py` → 断言树与比对 tip 集相同（V5.0 新增）
  - `adjudicate_assembly.py` → `assembly_adjudication.tsv`（V5.0 新增）
  - `qc_root_stability.py` → `QC3_root_stability.md`
  - `coordinate_mapper.py` → `coordinate_map.tsv`（断言一致）

---

## 4. 常见故障与处理策略（快速索引）

- **核心 MSA 列数暴涨（>1000）**：Insert 列未剥离或误用 `hmmalign --outformat afa`。回到 Phase 3.6。
- **Type II 序列核心被切半**：未做 hit stitching（[CHECK-06]）。
- **根位置在 MFP 与 LGC20 不一致**：LBA 风险；声明根不确定性，两种根假设下重复敏感性分析（QC3）。
- **ASR 报错 tip 集不匹配**：定根树包含 KDOPS 外群但 core_asr.afa 不包含。必须先 prune → `CoreTree_rooted_ingroup.treefile`。
- **模块 DCA 看起来"很漂亮"但 Meff/L < 3**：一律视为噪声，禁入主线（[CHECK-03]）。
- **祖先结构预测的装配体界面崩溃**：可能选错了 oligomeric state。检查 `assembly_adjudication.tsv`，考虑用 alternative oligomer 重新预测（[CHECK-08]）。
- **坐标映射对不上**：检查 PDB insertion code、ClipKIT 修剪列号变化、GROMACS 拓扑编号偏移。
- **全长缝合 MSA 找不到**：检查路径是否在 `results/03_msa_full/` 而非 `results/03_msa_core/`。

---

## 5. 交付物清单（最终验收）

- `results/03_msa_core/qc_core_alignment.md` ✅
- `results/03_msa_modules/boundary_robustness.md` ✅
- `results/04_phylogeny_asr/QC3_root_stability.md`
- `results/05_struct_valid/assembly_adjudication.tsv`
- `results/05_struct_valid/qc_struct_validation.md`
- `results/06_dca/qc_core_dca.md`
- `results/06_dca/core_significant_couplings.tsv`
