# TASKS.md — DAH7PS V5.1 项目任务追踪

> 基于 `PLAN.md`（V5.1 SOP rev6）提取，状态根据 `log.md` 实验记录同步更新。
> 所有数字以 `results/meta/metrics_manifest.tsv` 为准。
>
> 状态标记：`[x]` 完成 / `[/]` 进行中 / `[ ]` 待做 / `[-]` 跳过/延后

---

## Phase 0：环境与可复现性 ✅

- [x] 0.1 创建 conda 环境 `dah7ps_v4`（python 3.11, hmmer, mafft, iqtree, mmseqs2, foldmason, seqkit, cd-hit）
- [x] 0.1 安装 pip 依赖（clipkit, pastml）
- [x] 0.2 记录软件版本 → `results/meta/software_versions.tsv`
- [x] 0.3 建立 SOP 目录结构（data/, meta/, scripts/, results/01-06, workflow/）
- [x] 0.4 初始化参数文件 → `meta/params.json`
- [x] 0.2.1 初始化模型文件记录 → `results/meta/model_files.tsv`
- [x] 0.2.1 下载 3Di 模型（Q.3Di.AF/Q.3Di.LLM）→ `meta/models/`（sha256 已记录）
- [x] 0.5 **[V5.1]** 初始化文档同步：创建 `results/meta/metrics_manifest.tsv` + `results/meta/progress_snapshot.md`
- [x] 0.6 **[V5.1]** 创建 `results/04_phylogeny_asr/root_scenarios.tsv`
- [x] 0.7 **[V5.1]** 创建 `results/04_phylogeny_asr/node_selection_registry.tsv`
- [x] 0.8 **[V5.1]** 更新 `meta/params.json` → 新增 `rooting` / `manuscript` / `metrics` 参数块

**Done 条件：** ✅ `software_versions.tsv` + `params.json` + `model_files.tsv` 均存在

---

## Phase 1：数据挖掘（KDOPS 反向过滤 + 全库扫描） ✅

### V3.1 已完成（需评估是否满足 V4.1 要求）

- [x] 1.1 种子比对 + HMM 构建（Ia/Ib/II）→ `results/01_mining/model_*.hmm`
- [x] 1.2 hmmsearch 扫描 UniRef90 → `results/01_mining/domhits_*.tbl`
- [x] 1.2 序列提取 → `results/01_mining/raw_full_*.fasta`（Ia=10,102 / Ib=16,211 / II=9,862）
- [x] 1.3 KDOPS 反向过滤（仅 Iβ）→ `results/01_mining/raw_full_Ib_clean.fasta`

### V4.1 新增要求

- [x] 1.1 **扩充 Type II 种子**（1 → 14 条，覆盖细菌 + 植物多分支）
- [x] 1.1 扩充 Type Ia（3 → 13 条）/ Type Ib（3 → 6 条）种子多样性
- [x] 1.3 准备 KDOPS 外群种子集（12 条，覆盖 9 个分类群）→ `data/seeds/kdops_outgroup.fasta`
- [x] 1.1 用扩充后种子重建 HMM → 重新 hmmsearch（Ia=10,071 / Ib=15,608 / II=18,529）
- [x] 1.3 KDOPS 反向过滤 V4.1（Ib: 15,608 → 7,869 clean）→ `results/01_mining/hits_Ib_clean.fasta`
- [x] QC1 产出 `results/01_mining/qc_mining_report.md`

### V5.1 解释边界

> KDOPS 在 Phase 1 中是成功的负选择工具；但在 Phase 4 中不再自动等价于唯一可信的深根。

---

## Phase 2：质量控制与去冗余 ✅

### Pre-Phase 2 门控（计划外追加）

- [x] Gate A：三亚型 hits 交集统计 → Ia∩II = 8,561（85%）
- [x] Gate A-2：Best-hit Ia vs II 竞争打分归属 → 全部归属 Ia（100% HIGH）
- [x] Gate A-2：生成互斥 FASTA（Ia=10,071 / Ib=7,869 / II=9,968）
- [x] Gate B：Ib 边缘 KDOPS 序列隔离（4 条标记 → `kdops_borderline_ids.txt`）

### V4.1 执行

- [x] 2.1 长度 + HMM 覆盖度三箱过滤（`qc_length_coverage.py`）
  - Ia/Ib: cov_mode=best, II: cov_mode=merged (multi-domain stitching)
  - PASS 合计 24,202（Ia=9,204 / Ib=6,401 / II=8,597）
  - Type II rescued by stitching: 2,093 条
- [x] 2.2 CD-HIT 80% 去冗余 → `nr80_*.fasta`（Ia=3,521 / Ib=3,073 / II=3,079, 总计 9,673）
- [x] 2.3 CD-HIT 60% 种子提取 → `seeds60_*.fasta`（Ia=581 / Ib=648 / II=649, 总计 1,878）
- [x] 2.4a MMseqs2 40% 跨亚型 stepping stones → 258 条（Ia=104 / Ib=46 / II=108, 零跨亚型混簇）
- [x] 2.4b 二次聚类实验（25%→183, 20%→181）→ 确认需 Phase 3.1 结构优先级筛选
- [x] QC2 产出 `results/02_qc/qc_length_report.md`

**Done 条件：** ✅ 全部完成

---

## Phase 3：结构感知核心 MSA + 模块注释 ✅

### 3.1 结构面板构建 ✅

- [x] 3.1A-1 `panel_candidates.tsv`（258 条 backbone 全量评估）+ `panel_manifest.tsv`（30 条）
- [x] 线路 A 决策：PDB 作为外置锚点，不占 30 条配额
- [x] 3.1A-2 PDB 锚点下载 (1KFL, 1RZM, 3NV8, 5CKV, 2B7O)
- [x] 3.1A-3 AFDB 下载（30 条，core pLDDT ≥ 70）
- [x] 最终面板 = 30 AFDB + 5 PDB = 35 结构

### 3.2-3.4 FoldMason 骨架与核心列 ✅

- [x] 3.2 FoldMason easy-msa → 46 seqs × 1966 cols 骨架
- [/] 3.3 FoldMason refinemsa — segfault（msa2lddt 替代）
- [x] 3.4 逐列 LDDT 核心列界定 → core_len=521, lddt_min=0.1814

### 3.6 全量核心映射 ✅

- [x] 3.6 核心 HMM → `core_global.hmm`（L=521）
- [x] 3.6 `extract_core_domains.py`（含 hit stitching + pad=20）→ 9,393 seqs
- [x] 3.6 hmmalign → Stockholm → esl-alimask → `core_global_matchonly.afa`（9,393 × 521）

### 3.7 双版本修剪 ✅

- [x] 3.7 ClipKIT kpic-smart-gap → `core_tree.afa`（436 cols）
- [x] 3.7 Minimal trim gap>0.95 → `core_asr.afa`（472 cols）
- [-] 3.7 FoldMason msa2lddt 结构复核 — ✅ Average LDDT = 0.2638

### 3.8 模块注释 ✅

- [x] 3.8 模块 HMM 库构建
- [x] 3.8 `annotate_modules.py` → strict/relaxed 双矩阵（9,393 rows）
- [x] 3.8 Strict⊆Relaxed 一致性验证通过
- [x] 3.8 模块 MSA → ACT(47×142), CM(408×266), α2β3(172×582), N_ext(3,130×8,153), C_tail(360×3,617)

### 3.9 Profile-anchored Stitching ✅

- [x] 3.9.1 定义亚型与架构子集 → `results/03_msa_full/Ib_ACT.ids`（47 条）
- [x] 3.9.2 提取 linker 片段 → `results/03_msa_full/Ib_ACT_linkers.fasta`（C-flank，min=94,median=299,max=348 aa）
- [x] 3.9.3 亚型内 linker E-INS-i 对齐 → `results/03_msa_full/Ib_ACT_linkers_einsi.afa`（47×426）
- [x] 3.9.4 `stitch_full_length_msa.py` → `results/03_msa_full/msa_full_Ib_v4.afa`（47×1040）+ `results/03_msa_full/msa_full_Ib_column_map.tsv`
  - ✅ QC2b 断言通过：core 段 472 列与 core_asr.afa 完全一致
  - 列分段：core(1–472) | ACT(473–614) | C-flank(615–1040)

---

## Phase 4：系统发育与分层 ASR [进行中，V5.1 重点修订]

> **V5.1 规则：计算可以继续，叙事必须分层。working tree ≠ narrative tree。**

### 4.1 核心树推断与 root scenario 管理

- [x] 4.1 KDOPS 外群并入核心比对（hmmalign sto → strip insert）
- [x] 4.1 全局树 baseline（MFP）→ `CoreTree_rooted_MFP.treefile`（Q.PFAM+F+R10, 1001 iter, LogL=-2835499.68）
- [ ] 4.1 全局树 LBA 抗性（LG+C20+F+G）→ `CoreTree_rooted_LGC20.treefile`
- [ ] 4.1 **[V5.1]** 补跑 S1/S2 支持度版本（`-alrt 1000 -B 1000`）（若 full-tree 太贵则在 focal subtree 上补做）
- [x] 4.1 **[V5.1]** 构建 S3 ingroup-only midpoint rooted tree → `CoreTree_rooted_midpoint_ingroup.treefile` ✅
- [x] 4.1 **[V5.1]** 构建 S4 ingroup-only MAD rooted tree → `CoreTree_rooted_MAD_ingroup.treefile` ✅ (rho=0.163617, split=[3060,6333])
- [ ] 4.1 **[V5.1]** 可选 S5：reduced representative set nonreversible/rootstrap
- [ ] 4.1 **[V5.1]** 维护 `results/04_phylogeny_asr/root_scenarios.tsv`

### 4.2 结构系统发育交叉验证 [V5.1 升级为高优先级]

- [x] 4.2 AA 树 vs 3Di 树交叉验证 → `tree_comparison.md` + `tree_comparison.tsv` ✅ (nRF=0.7442, 三大亚型单系一致, QC3-YELLOW)
- [x] 4.2 **[V5.1]** 三个判读问题：Q1 通过（亚型单系一致） / Q2 不适用（骨架树无 KDOPS） / Q3 无致命冲突（高 nRF 为 AA vs 3Di 预期行为）

### 4.3 核心氨基酸 ASR（可继续计算，解释分层）

- [x] 4.3 **Prune KDOPS → `CoreTree_rooted_ingroup.treefile`**（12 KDOPS pruned → 9,393 tips, bifurcating root）⚠ KDOPS polyphyletic in ML tree
- [x] 4.3 **`assert_tip_match.py` 断言 pruned tree tips == core_asr.afa tips** ✅ PASS 9,393/9,393
- [🏃] 4.3 核心氨基酸 ASR（用 pruned ingroup tree，Q.PFAM+F+R10, PID 23826）
- [ ] 4.3 **[V5.1]** 至少补跑 1 个 alternative scenario 的 ASR（优先 S2 或 S4）
- [ ] 4.3 **[V5.1]** 产出 `asr_node_summary.tsv`

### 4.4 嵌套 ASR [立即并行推进]

- [ ] 4.4 嵌套 ASR（Iβ-ACT exemplar，输入 `msa_full_Ib_v4.afa`）
- [ ] 4.4 **[V5.1]** focal Pre-gain / Post-gain 节点导出：ML 序列 + AltAll 序列 + PP 概要

### 4.5 Gap 祖先态重建

- [ ] 4.5 Gap 祖先态重建（方案 A/B/C 分层策略）
- [ ] 4.5 **[V5.1]** 所有 gap ASR 使用 pruned ingroup rooted scenario trees

### 4.6 模块 trait ASR（PastML）[V5.1 关键加强]

- [ ] 4.6 模块获得/丢失离散性状 ASR（PastML）
- [ ] 4.6 **[V5.1]** 必须运行 strict/relaxed × 至少 2 root scenarios（S1, S2 或 S4）
- [ ] 4.6 **[V5.1]** 输出 `module_origin_stability.tsv`（robust / semi-robust / root-sensitive / annotation-sensitive / unresolved）

### 4.7–4.8 AltAll 与模块层 ASR

- [ ] 4.7 AltAll 系综采样（PP₁=0.80, PP₂=0.20）— **[V5.1]** 只对 Phase 5 候选 + gain/loss 边界 + exemplar 关键祖先节点
- [ ] 4.8 模块层局部 ASR（探索性，不作为主线 gating）

### QC3：系统发育 / ASR 统一 gate [V5.1 强制执行]

- [ ] QC3 **[V5.1]** root scenario 列表与对应模型/方法
- [ ] QC3 **[V5.1]** UFBoot + SH-aLRT 支持度记录
- [ ] QC3 **[V5.1]** KDOPS subset sensitivity
- [ ] QC3 **[V5.1]** midpoint / MAD 对照
- [ ] QC3 **[V5.1]** AA tree vs 3Di tree 对照
- [ ] QC3 **[V5.1]** strict vs relaxed 模块事件稳定性
- [ ] QC3 **[V5.1]** focal node 多场景同向判断
- [ ] QC3 **[V5.1]** Phase 5 推荐：Go / Conditional Go / Hold
- [ ] QC3 产出 `results/04_phylogeny_asr/qc_phylogeny_asr.md` + `.tsv`

---

## Phase 5：关键祖先节点的结构验证 [待执行，V5.1 更严格门控]

> **V5.1：2–4 节点 × Apo-only × native-oligomer。不做 Holo。节点必须通过 QC3 gate + 四重门控。**

- [ ] 5.0 **装配体判定（Assembly Adjudication）[CHECK-08]**
  - 文献/PDB/PISA 装配体注释扫描
  - dimer + tetramer 平行 AF3 预测（均 Apo）
  - PISA 界面评分比较
  - 产出 `results/05_struct_valid/assembly_adjudication.tsv`
- [ ] 5.1 **[V5.1]** 候选祖先节点选择（四重门控：UFBoot ≥ 95 + root stability ≥ 2 scenarios + strict/relaxed stability + ASR interpretability）
  - Tier-1：主文主图（2 节点）
  - Tier-2：补充/reserve（1–2 节点）
  - 产出 `results/04_phylogeny_asr/node_selection_registry.tsv`
- [ ] 5.2 Apo 结构预测（AF3 拷贝数由 5.0 判定 + ESMFold 交叉验证）
- [ ] 5.2b 结构 QC 门控 [CHECK-05]：ipTM ≥ 0.6 + 能量最小化 + 10 ns 筛选
- [ ] 5.3 验证性 MD：native-oligomer-apo ≥200 ns × 2 rep（Tier-1 only）
- [ ] 5.3 可选：alternative-oligomer-apo 10–20 ns × 1 rep（排除性测试）
- [ ] 5.4 有限动力学读出（RMSD/RMSF + 界面面积 + Fpocket）— 只比较 root-robust tier-1 节点
- [ ] QC4 产出 `results/05_struct_valid/qc_struct_validation.md`

---

## Phase 6：核心层共进化分析（DCA） [待执行，V5.1 聚焦版]

> **V5.1：仅核心层 DCA 进入主线。模块/联合 DCA → 可选探索。**

### 主线分析

- [ ] 6.1 核心层 DCA 输入准备（`core_asr.afa` → gap 过滤 → `core_dca.afa`，Meff/L 门控 ≥ 3.0）[可立即并行]
- [ ] 6.2 plmc 执行 → `core_couplings.txt`
- [ ] 6.3 显著性评估：top-L 接触验证（1KFL/1RZM/3NV8）
- [ ] 6.3 功能位点富集 + 跨亚型保守 vs 特异耦联
- [ ] QC5 产出 `results/06_dca/qc_core_dca.md`

### 可选探索（不入论文主线结论）

- [ ] 6.4.1 模块层 DCA（探索性，ACT Meff/L ≈ 0.2–0.3，深度不足）
- [ ] 6.4.2 联合跨域 DCA（探索性）
- [ ] 6.4.3 模块获得前后耦联变化比较（Meff 匹配下采样，探索性）

---

## Phase 7：论文写作蓝图 [待执行，V5.1 重写主叙事]

> **V5.1 写作核心：从"唯一深根故事"改为"root-robust 主叙事"。**

- [ ] 7.1 论文主线叙事锁定（root-robust claims → main text; root-sensitive → Supplement/Discussion）
- [ ] 7.2 核心图表 Fig 1–6 齐备
  - Fig 3 **[V5.1]** root scenario 比较 + stable/unstable module events
  - Fig 4 **[V5.1]** Iβ-ACT nested full-length ASR exemplar
- [ ] 7.3 ICDC 在 Discussion 中的定位（"跨证据一致性展望"，定性描述，非正式融合）
- [ ] 7.4 **[V5.1]** Claim tiers 落地（root_robust / root_sensitive / exploratory）

---

## 脚本开发

### 已完成
- [x] `scripts/extract_core_domains.py`（含 hit stitching）
- [x] `scripts/define_core_columns.py`（LDDT 拐点法）
- [x] `scripts/annotate_modules.py`（strict/relaxed 双矩阵）
- [x] `scripts/extract_module_seqs.py`（C-tail + 按 matrix 逐模块）
- [x] `scripts/minimal_trim.py`
- [x] `scripts/extract_struct_subset.py`
- [x] `scripts/extract_linkers.py`
- [x] `scripts/stitch_full_length_msa.py`（含 core 列不变断言，输出到 `results/03_msa_full/`）
- [x] `scripts/select_sequences.py`
- [x] `scripts/prune_tree.py`
- [x] `scripts/assert_tip_match.py`

### 高优先补齐 [V5.1]
- [ ] `scripts/qc_root_stability.py`
- [ ] `scripts/compare_trees.py`
- [ ] `scripts/adjudicate_assembly.py`
- [ ] `scripts/prepare_dca_input.py`（含 Meff/L 门控）
- [ ] `scripts/compute_meff.py`
- [ ] `scripts/dca_significance.py`
- [ ] `scripts/coordinate_mapper.py`
- [ ] `scripts/build_metrics_manifest.py` **[V5.1 新增]**
- [ ] `scripts/summarize_module_stability.py` **[V5.1 新增]**

---

## 最终交付物

- [x] `results/03_msa_core/qc_core_alignment.md`
- [x] `results/03_msa_modules/boundary_robustness.md`
- [ ] `results/04_phylogeny_asr/qc_phylogeny_asr.md` **[V5.1 替代旧 QC3_root_stability.md]**
- [ ] `results/04_phylogeny_asr/root_scenarios.tsv` **[V5.1 新增]**
- [ ] `results/04_phylogeny_asr/module_origin_stability.tsv` **[V5.1 新增]**
- [ ] `results/04_phylogeny_asr/node_selection_registry.tsv` **[V5.1 新增]**
- [ ] `results/05_struct_valid/assembly_adjudication.tsv`
- [ ] `results/05_struct_valid/qc_struct_validation.md`
- [ ] `results/06_dca/core_significant_couplings.tsv`
- [ ] `results/06_dca/qc_core_dca.md`
- [x] `results/meta/metrics_manifest.tsv` **[V5.1 新增]** ✅ 2026-03-18 初始化（39 指标）
- [x] `results/meta/progress_snapshot.md` **[V5.1 新增]** ✅ 2026-03-18 初始化
