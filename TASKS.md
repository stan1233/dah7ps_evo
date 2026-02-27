# TASKS.md — DAH7PS V4.1 项目任务追踪

> 基于 `PLAN.md`（V4.1 SOP rev4）提取，状态根据 `log.md` 实验记录同步更新。
>
> 状态标记：`[x]` 完成 / `[/]` 进行中 / `[ ]` 待做 / `[-]` 跳过/延后

---

## Phase 0：环境与可复现性

- [x] 0.1 创建 conda 环境 `dah7ps_v4`（python 3.11, hmmer, mafft, iqtree, mmseqs2, foldmason, seqkit, cd-hit）
- [x] 0.1 安装 pip 依赖（clipkit, pastml）
- [x] 0.2 记录软件版本 → `results/meta/software_versions.tsv`
- [x] 0.3 建立 SOP 目录结构（data/, meta/, scripts/, results/01-06, workflow/）
- [x] 0.4 初始化参数文件 → `meta/params.json`
- [x] 0.2.1 初始化模型文件记录 → `results/meta/model_files.tsv`
- [x] 0.2.1 下载 3Di 模型（Q.3Di.AF/Q.3Di.LLM）→ `meta/models/`（用户手动完成，sha256 已记录）

**Done 条件：** ✅ `software_versions.tsv` + `params.json` + `model_files.tsv` 均存在

---

## Phase 1：数据挖掘（KDOPS 反向过滤 + 全库扫描）

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

---

## Phase 2：质量控制与去冗余

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
- [x] 2.3 CD-HIT 60% 种子提取 → `seeds60_*.fasta`（Ia=581 / Ib=648 / II=649, 总计 1,878, Gate B 已排除）
- [x] 2.4a MMseqs2 40% 跨亚型 stepping stones → Coverage Backbone 258 条（Ia=104 / Ib=46 / II=108, 零跨亚型混簇）
- [x] 2.4b 二次聚类实验（25%→183, 20%→181）→ 确认需 Phase 3.1 结构优先级筛选
- [x] QC2 产出 `results/02_qc/qc_length_report.md`

### Stepping Stones 两层定义（决策记录）

- **Coverage Backbone（258 条）**：`stepping_stones_rep_seq.fasta`，用于覆盖度验证
- **Structure Panel（目标 20-40 条）**：Phase 3.1 按 PDB > AFDB > ESMFold 优先级筛选
- SOP 的 20-50 约束重新绑定到 Structure Panel

**Done 条件：** ✅ `nr80_*.fasta` + `qc_length_report.md` + seeds60 + stepping stones 均产出

---

## Phase 3：结构感知核心 MSA + 模块注释

### 3.1 结构面板构建（Selection Contract 已锁定）

**Selection Contract**: N=30 (20–40), quota Ia=12/Ib=5/II=13 (±1; floors Ia≥8/Ib≥4/II≥8)
Priority: PDB > AFDB (core pLDDT≥70) > ESMFold (core pLDDT≥70). Core-region confidence, not full-length.

- [x] 3.1A-1 生成 `panel_candidates.tsv`（258 条 backbone 全量评估：PDB=0, AFDB=154, needs_ESMFold=124）
- [x] 3.1A-1 Selection 验证 → `panel_manifest.tsv`（30 条：Ia=12/Ib=5/II=13，全 AFDB）
- [x] **决策已定**：选择线路 A（PDB 作为外置锚点，不占 30 条配额）
- [x] 3.1A-2 PDB 锚点结构下载 (1KFL, 1RZM, 3NV8, 5CKV, 2B7O) → `data/structures/panel_dah7ps/PDB-*.cif`
- [x] 3.1A-3 AFDB 结构下载（30 条，core pLDDT ≥ 70，0 errors）→ `data/structures/panel_dah7ps/AF-*.pdb`
- [-] 3.1A-4 ESMFold 补缺口 → 不需要（30 AFDB 全部满足质量门槛）
- [x] 3.1A-5 最终面板 = 30 AFDB + 5 PDB 外置锚点 = 35 结构（在 SOP 20-40 范围内）
- [ ] 3.1B 下载 KDOPS 外群结构 → `data/structures/panel_kdops/`（Phase 4.1 再做）

### 3.2-3.4 FoldMason 骨架与核心列

- [x] 3.2 FoldMason easy-msa → `results/03_msa_core/skeleton_aa.fa` (46 seqs = 30 AFDB + 16 PDB chains) + `skeleton_3di.fa` + `skeleton.html` (5.6M)
- [/] 3.3 FoldMason refinemsa — ⚠ segfault（`createdb` 过滤 24 短链导致 DB-MSA 序列数不匹配；`msa2lddt` + `msa2lddtreport` 已跑成功作为替代质量评估）
- [ ] 3.4 逐列 LDDT 核心列界定（`define_core_columns.py`）→ `skeleton_core_aa.fa` + `core_columns.mask`

### 3.5-3.6 亚型内骨架与全量映射

- [ ] 3.5 亚型内种子 E-INS-i 骨架 → `results/03_msa_core/seeds60_*_einsi.afa`
- [ ] 3.6 核心 HMM 构建 → `core_global.hmm`
- [ ] 3.6 `extract_core_domains.py`（含 hit stitching [CHECK-06]）→ `all_core_only.fasta` + `core_domain_coords.tsv`
- [ ] 3.6 hmmalign → Stockholm → esl-alimask 剥离 Insert → `core_global_matchonly.afa`
- [ ] QC2 产出 `results/03_msa_core/qc_core_alignment.md`

### 3.7 双版本修剪

- [ ] 3.7 ClipKIT kpic-smart-gap → `core_tree.afa`
- [ ] 3.7 Minimal trim → `core_asr.afa`
- [ ] 3.7 FoldMason msa2lddt 结构复核

### 3.8 模块注释

- [ ] 3.8 `annotate_modules.py` → `module_presence_absence_strict.tsv` + `_relaxed.tsv`
- [ ] 3.8 模块边界稳健性 QC（三证据一致性 + 双版本矩阵）
- [ ] 3.8 `extract_module_seqs.py` → 各模块序列 + `*_domain_coords.tsv`
- [ ] 3.8 模块 MSA（MAFFT E-INS-i）→ `ACT_msa.afa` 等

### 3.9 Profile-anchored Stitching

- [ ] 3.9.1 定义亚型与架构子集
- [ ] 3.9.2 提取 linker 片段
- [ ] 3.9.3 亚型内 linker E-INS-i 对齐
- [ ] 3.9.4 `stitch_full_length_msa.py` → `msa_full_Ib_v4.afa` + `column_map.tsv`（断言 core 列不变）

---

## Phase 4：系统发育与分层 ASR

- [ ] 4.1 KDOPS 外群并入核心比对（hmmalign sto → strip insert）
- [ ] 4.1 全局树 baseline（MFP）→ `CoreTree_rooted_MFP.treefile`
- [ ] 4.1 全局树 LBA 抗性（LG+C20+F+G）→ `CoreTree_rooted_LGC20.treefile`
- [ ] 4.1 QC3 根稳定性报告 → `QC3_root_stability.md`
- [ ] 4.2 AA 树 vs 3Di 树交叉验证（需 3Di 模型文件）
- [ ] 4.3 核心氨基酸 ASR
- [ ] 4.4 嵌套 ASR（亚型内全长比对，输入来自 Phase 3.9）
- [ ] 4.5 Gap 祖先态重建（方案 A/B/C 分层策略）
- [ ] 4.6 模块获得/丢失离散性状 ASR（PastML，严格/宽松双版本敏感性）
- [ ] 4.7 AltAll 系综采样（PP₁=0.80, PP₂=0.20）
- [ ] 4.8 模块层局部 ASR

---

## Phase 5：结构预测与多聚体 MD

- [ ] 5.1 候选祖先集合定义（Pre-gain / Post-gain / 对照）
- [ ] 5.2 Apo AF3 四聚体预测（每个祖先节点）
- [ ] 5.2a Apo-first 门控 [CHECK-07]：口袋检测 → 对接 → holo 权限
- [ ] 5.2b MD 前结构 QC [CHECK-05]：ipTM/PAE → 能量最小化 → 10ns 筛选
- [ ] 5.3 2×2 因子矩阵 MD（四聚体-apo, 四聚体-effector, 单体对照）
- [ ] 5.4 动力学读出（RMSD/RMSF, DCCM, 网络流）
- [ ] QC4 产出 `results/05_struct_md/qc_dynamics.md`

---

## Phase 6：ICDC — 分层 DCA × 动力学网络整合

- [ ] 6.1.1 核心层 DCA 输入准备（`prepare_dca_input.py`，Meff/L ≥ 3.0 门控）
- [ ] 6.1.2 模块层 DCA 输入（Meff/L < 3 自动跳过 [CHECK-03]）
- [ ] 6.1.3 亚型内联合 DCA 输入
- [ ] 6.2.1 核心层 plmc DCA
- [ ] 6.2.2 模块层 DCA
- [ ] 6.2.3 联合 DCA → 跨域 core↔module 耦联 `Ib_ACT_crossdomain_top200.tsv`
- [ ] 6.2.4 显著性评估（top_L, Cβ-Cβ < 8Å）
- [ ] 6.2.5 模块获得前后耦联变化（Meff 匹配下采样 + Z-score）
- [ ] 6.3 ICDC 融合（五套坐标系映射验证 `coordinate_mapper.py` [CHECK-04]）
- [ ] 6.3 产出 `icdc_core_network.graphml` + `icdc_crosslayer_paths.tsv`

---

## 脚本开发（SOP 引用但尚未实现）

- [ ] `scripts/extract_core_domains.py`（含 hit stitching）
- [ ] `scripts/define_core_columns.py`（LDDT 拐点法）
- [ ] `scripts/annotate_modules.py`（多证据融合 + 双版本矩阵）
- [ ] `scripts/extract_module_seqs.py`
- [ ] `scripts/extract_linkers.py`
- [ ] `scripts/stitch_full_length_msa.py`（含 core 列不变断言）
- [ ] `scripts/prepare_dca_input.py`（含 Meff/L 门控）
- [ ] `scripts/qc_root_stability.py`
- [ ] `scripts/coordinate_mapper.py`（五套坐标系一致性断言）
- [ ] `scripts/select_sequences.py`
- [ ] `scripts/merge_alignments.py`
- [ ] `scripts/compare_trees.py`
- [ ] `scripts/aggregate_gap_blocks.py`
- [ ] `scripts/split_msa_by_module.py`
- [ ] `scripts/compute_meff.py`
- [ ] `scripts/downsample_to_meff.py`
- [ ] `scripts/run_plmc_batch.py`
- [ ] `scripts/dca_compare_meff_matched.py`
- [ ] `scripts/dca_significance.py`
- [ ] `scripts/extract_crossdomain_couplings.py`
- [ ] `scripts/qc_length.py`
- [ ] `scripts/minimal_trim.py`
- [ ] `scripts/extract_struct_subset.py`
- [ ] `scripts/drop_a2m_insertions.py`（fallback 用）
- [ ] `scripts/subset_msa_by_ids.py`

---

## 最终交付物

- [ ] `results/03_msa_core/qc_core_alignment.md`
- [ ] `results/04_phylogeny_asr/QC3_root_stability.md`
- [ ] `results/05_struct_md/qc_dynamics.md`
- [ ] `results/06_icdc/icdc_core_network.graphml`
- [ ] `results/06_icdc/icdc_module_network.graphml`
- [ ] `results/06_icdc/icdc_crosslayer_paths.tsv`
- [ ] `results/06_icdc/icdc_crossdomain_couplings.tsv`（若跑联合 DCA）
