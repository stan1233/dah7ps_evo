# DAH7PS 变构起源的演化动力学重建：V5.0 标准作业程序（SOP rev5）

**A FoldMason-Centred, Core–Module Decoupling Framework for Reconstructing the Evolutionary Origin of Allostery in DAH7PS**

---

## 摘要

DAH7PS（3-deoxy-D-arabino-heptulosonate 7-phosphate synthase）在保守的 (β/α)₈ TIM-barrel 催化核心之上，于不同亚型中独立演化出多样的变构调控元件——Type Iα 的 β 发夹插入、Type Iβ 的 ACT/铁氧还蛋白样调控域、Type II 的 N 端延伸与 α2β3 插片。V3.1 管线的"全长单一 MSA 驱动一切"策略在跨亚型远缘同源（序列一致性进入暮光区 <20%）时已产生严重比对膨胀（10–17:1），系统性污染下游 ASR、DCA 与动力学推断。

V5.0 在 V4.1 的**核心–模块分治（core–module decoupling）**范式基础上，根据 Phase 1–3.8 的实际执行经验做出关键范围调整：

1. **Phase 1–3.8（已完成，不变）**：以 FoldMason 结构字母表（3Di+AA）与逐列 LDDT 物理置信度为中枢，构建了 9,393 条序列 × 521 列的去膨胀核心 MSA、双版本修剪比对、5 类模块注释（strict/relaxed 双矩阵）。Phase 3.9（亚型内全长缝合比对）待执行——它是 Phase 4.4 嵌套 ASR 的前置依赖。
2. **Phase 4（系统发育 + ASR，保持完整）**：核心树定根、AA/3Di 交叉验证、分层 ASR、模块获得/丢失离散性状重建——这些是论文演化叙事的时间轴骨架，不可缩减。
3. **Phase 5（大幅缩减 → 验证性结构/MD）**：从"8 节点 × 4 条件全因子 MD"缩减为"2–4 个关键祖先节点的 Apo 结构预测 + 有限验证性 MD"，目标从"建立完整动力学通讯网络"降级为"为演化推断提供结构可行性验证"。
4. **Phase 6（DCA 聚焦核心层）**：仅在信息量充足的核心层（Meff/L ≥ 3.0）执行 DCA；模块 DCA 与联合跨域 DCA 因深度不足不作为主线证据。
5. **ICDC 重新定位**：从 V4.1 的"第六步主线结论"降级为 Discussion 中的"跨证据一致性展望（cross-evidence convergence outlook）"——点到即止地展示核心 DCA、祖先结构、有限 MD 之间的一致性信号，但不宣称建立了闭环证据链。

**论文主线叙事的证据支柱变为**：结构感知核心 MSA（Phase 3）→ 定根系统发育 + 模块离散性状 ASR（Phase 4）→ 核心层共进化信号（Phase 6 core DCA）→ 关键祖先节点的结构验证（Phase 5 精简版）。ICDC 融合仅作为未来方向展望。

---

## V5.0 修订要点（相对 V4.1）

| 编号 | 改动 | 理由 |
|------|------|------|
| 1 | Phase 5 祖先节点从 ≤8 缩减为 **2–4 个**，全部默认 **Apo-only** | 资源聚焦；Holo 预测的"现代配体幻觉"风险过高，收益不确定 |
| 2 | Phase 5 MD 从"2×2 因子实验"缩减为"Apo native-oligomer ≥200 ns × 2 rep"验证性轨迹；**新增 5.0 装配体判定** | 目标从"建立动力学网络"降为"验证装配体稳定性与口袋拓扑"；装配体状态不默认 |
| 3 | Phase 6 DCA **仅保留核心层**；模块 DCA、联合跨域 DCA、Meff 匹配比较 → 全部移入"可选探索" | ACT strict 仅 47 序列（Meff/L ≈ 0.2–0.3），远低于 3.0 门控 |
| 4 | ICDC 从"第六步主线产出"降级为 **Discussion 展望** | 无完整 MD 动力学网络，无法做 DCA×DCCM 正式融合 |
| 5 | 新增 **Phase 7：论文写作蓝图**，锁定主线叙事与图表清单 | 防止下游执行时叙事漂移 |
| 6 | Phase 1–3.8 标记为 ✅ 已完成，Phase 3.9 保留待执行 | 不可修改已完成步骤；3.9 是 4.4 的前置依赖 |

---

## 科学问题与工作假设

**核心问题：** DAH7PS 的变构调控元件如何在保守 TIM-barrel 核心之上被募集、稳定并与多聚体界面及动力学网络耦合，从而产生可选择的变构表型？

**工作假设：**

1. 跨亚型可可靠比较的演化"时间轴"主要由 TIM-barrel 核心提供；变构元件的起源在演化上表现为模块获得/丢失事件与模块内序列—结构协变。
2. 变构功能的形成与维持依赖多聚体组装态与跨亚基通讯路径（Cross et al., 2013; Lang et al., 2016）。
3. 核心层序列共进化（DCA）可独立于 MD 为变构相关位点对提供演化约束证据；完整的 DCA×动力学融合留待深度充足时实现。

---

## 方法学取舍声明

1. **FoldMason 定位：** 中枢工具，非唯一判据。保留多证据链交叉验证（催化位点几何一致性、AA 树 vs 3Di 树拓扑、外群定根稳定性）。
2. **MAFFT 角色：** 不适合跨亚型无约束全长裸跑，但在亚型内/模块子集高精度 MSA（E-INS-i）以及受约束增量映射（`--add --keeplength`）场景仍具方法学价值。
3. **位点保护：** 不对特定残基做硬编码保护（避免确认偏差）。以 FoldMason-LDDT 高置信列为一级保护标准。
4. **参数管理：** 所有阈值外部化至 `meta/params.json`，不硬编码于脚本中。
5. **MD 定位（V5.0 新增）：** MD 在本研究中仅作为"祖先结构可行性验证"，而非"建立完整的动力学通讯网络"。后者需要更大规模的计算投入，是本研究的自然延伸方向。
6. **DCA 定位（V5.0 新增）：** 仅核心层 DCA 具有充足信息量（Meff/L > 5）进入主线证据链。模块层和跨域联合 DCA 因有效序列深度不足，不作为论文主线结论的定量证据。

---

# ✅ 第零步：计算环境与项目规范 [已完成]

> **状态：2026-02-23 完成。** 以下内容保留供审计，不可修改。

## 0.1 可复现环境

**层 1：conda 核心环境**

```bash
mamba create -n dah7ps_v4 \
  python=3.11 hmmer mafft iqtree mmseqs2 foldmason seqkit cd-hit \
  -y
conda activate dah7ps_v4
pip install clipkit pastml
```

**层 2：HPC / 容器依赖**

```
工具              最低版本     安装方式
──────────────────────────────────────────────
GROMACS           ≥2023.x     HPC module 或 Singularity/Apptainer 容器
AlphaFold3        参见 DeepMind AF3 仓库
ESMFold           参见 ESM GitHub
PyMOL / ChimeraX  可选         结构可视化
```

## 0.2 软件版本锁定

已记录至 `results/meta/software_versions.tsv`（13 行，含 python/hmmer/mafft/iqtree/mmseqs/foldmason/seqkit/cd-hit/clipkit/pastml 等）。

3Di 模型文件（Q.3Di.AF / Q.3Di.LLM）已下载至 `meta/models/`，sha256 校验已写入 `results/meta/model_files.tsv`。

## 0.3 目录规范

```
data/
  seeds/                # 手工整理的种子序列与结构
  structures/
    panel_dah7ps/       # ⚠ 仅 DAH7PS 三亚型代表（核心列界定面板）
    panel_kdops/        # ⚠ 仅 KDOPS 外群结构（定根面板）
  db/                   # UniRef/UniProt 数据库索引
results/
  01_mining/            # HMM 搜索 + KDOPS 过滤
  02_qc/                # 长度过滤 + 去冗余
  03_msa_core/          # 核心层 MSA（FoldMason 骨架 + hmmalign 映射）
  03_msa_modules/       # 模块层 MSA（分模块独立比对）
  03_msa_full/          # 亚型内全长缝合 MSA
  04_phylogeny_asr/     # 系统发育 + 分层 ASR
  05_struct_valid/      # 结构预测 + 验证性 MD（V5.0 更名）
  06_dca/               # 核心层 DCA（V5.0 更名，聚焦核心）
  meta/                 # 参数文件、版本锁定、QC 报告
scripts/
```

## 0.4 参数文件

`meta/params.json` 已初始化，包含 mining, qc, core_definition, msa, phylogeny, dca, af3, asr 全部参数块。

---

# ✅ 第一步：全长序列挖掘与 KDOPS 免疫过滤 [已完成]

> **状态：2026-02-24 完成。** QC1 通过。

## 实际执行摘要

**种子扩充（V3.1 → V4.1）：** Ia 从 3 条扩至 13 条，Ib 从 3 条扩至 6 条，II 从 1 条扩至 14 条。扩充过程中发现并纠正 6 个错误 UniProt ID。KDOPS 外群 12 条，覆盖 9 个分类群。

**HMM 搜索结果：** Ia=10,071 / Ib(raw)=15,608 / II=18,529。Type II 命中数翻倍验证种子扩充有效。

**KDOPS 反向过滤（Ib）：** 双向竞争打分移除 7,739 条（49.6%），保留 7,869 条。得分差呈强双峰分布，边界清晰。4 条边缘序列（delta 13.7–16.6）隔离进 borderline 集。

**Gate 检查（计划外追加）：**
- Gate A：发现 Ia∩II = 8,561（85%），全部以 ≥20 bits 优势归属 Ia。Best-hit 归属后三亚型互斥。
- Gate B：Ib 边缘 4 条已隔离。

**Phase 1 最终互斥输入集合：** Ia=10,071 / Ib=7,869 / II=9,968（总计 27,908）

### QC1 产出

`results/01_mining/qc_mining_report.md` ✓

---

# ✅ 第二步：数据驱动 QC 与去冗余 [已完成]

> **状态：2026-02-25 完成。** QC2 通过。

## 实际执行摘要

**Phase 2.1 三箱过滤：** 引入 HMM 覆盖度 ≥ 0.70 + 长度窗口联合过滤。Type II 因 α2β3 内插片导致 hmmsearch 碎片化 hits（FRAG 高达 34.8%），通过 multi-domain stitching 修复（FRAG 降至 13.8%，rescued 2,093 条）。

| 亚型 | PASS_CANONICAL | PASS_LONG | FRAG |
|------|---------------|-----------|------|
| Ia | 9,022 (89.6%) | 182 (1.8%) | 867 (8.6%) |
| Ib | 6,229 (79.2%) | 172 (2.2%) | 1,468 (18.7%) |
| II | 8,457 (84.8%) | 140 (1.4%) | 1,371 (13.8%) |

**Phase 2.2 CD-HIT 80%：** nr80 产出 Ia=3,521 / Ib=3,073 / II=3,079（总计 9,673）

**Phase 2.3 CD-HIT 60% 种子：** seeds60 产出 Ia=581 / Ib=648 / II=649（总计 1,878）

**Phase 2.4 Stepping Stones：** MMseqs2 40% identity 聚类 → 258 个代表（Ia=104, Ib=46, II=108）。零跨亚型混簇。258 条保留为 Coverage Backbone；结构面板从中按 PDB > AFDB > ESMFold 优先级筛选。

### QC2 产出

`results/02_qc/qc_length_report.md` ✓

---

# ✅ 第三步：结构感知 MSA——核心/模块分治，FoldMason 为中枢 [Phase 3.1–3.8 已完成]

> **状态：Phase 3.1–3.8 于 2026-03-03 完成。Phase 3.9（缝合比对）待执行。**

## Phase 3.1 实际执行：结构面板（35 结构）

Selection Contract：目标 N=30 AFDB + PDB 锚点外置。实际产出 35 结构（30 AFDB + 5 PDB：1KFL/1RZM/3NV8/5CKV/2B7O）。配额 Ia=12 / Ib=5 / II=13。

## Phase 3.2 实际执行：FoldMason 结构骨架

`foldmason easy-msa` 产出 46 序列（30 AFDB 单链 + 16 PDB 多链拆分）× 1966 列骨架比对。同时产出 3Di 比对和引导树。

## Phase 3.3 实际执行：refinemsa

`refinemsa` 因 createdb 与 skeleton 序列数不匹配而 segfault。替代方案：直接使用 `msa2lddt` + `msa2lddtreport` 评估原始骨架质量（msaLDDT = 0.2391）。

## Phase 3.4 实际执行：逐列 LDDT 核心列界定

Auto_inflection knee LDDT 阈值 = 0.1814。基础核心列 302 → ±20 padding 后最终核心列 **521** 列（目标 400–600 ✓）。

## Phase 3.6 实际执行：全量去膨胀映射

全程走 Stockholm → `esl-alimask --rf-is-mask` → AFA 路径（未使用 `--outformat afa` 捷径）。

核心域提取含 hit stitching（1,069 条 = 11.4% 需要合并多段 hit）。pad=20 修复端部 gap 死区（改善有限，确认为真实生物学变异）。

**最终核心 MSA：** `core_global_matchonly.afa` — **9,393 seqs × 521 cols** ✓

## Phase 3.7 实际执行：双版本修剪

| 比对 | 列数 | 用途 |
|------|------|------|
| `core_tree.afa` | 436 | 系统发育建树 |
| `core_asr.afa` | 472 | ASR / DCA |

FoldMason msa2lddt 结构复核：30 条面板子集平均 LDDT = 0.2638（> inflection 0.1814 ✓）。

## Phase 3.8 实际执行：模块注释

5 类模块（N_ext / α2β3_insert / ACT_domain / CM_domain / C_tail），strict + relaxed 双矩阵（9,393 行）。Boundary confidence: high=79.1%, medium=20.8%, low=0.0%。

模块 MSA 已构建：ACT(47 seqs × 142 cols), CM(408 × 266), α2β3(172 × 582), N_ext(3,130 × 8,153), C_tail(360 × 3,617)。

### QC2/QC2b 产出

- `results/03_msa_core/qc_core_alignment.md` ✓
- `results/03_msa_modules/boundary_robustness.md` ✓

---

# Phase 3.9：轮廓锚定缝合比对（Profile-anchored Stitching MSA） [待执行]

> **目的：** 为 Phase 4.4 嵌套 ASR 与可选的联合 DCA 生成**高质量、非膨胀**的"亚型内全长比对"（例如 `msa_full_Ib_v4.afa`）。
>
> **依赖：** Phase 3.6（`core_asr.afa`）+ Phase 3.8（模块 MSA + `core_domain_coords.tsv`）

**原则：**

- **核心列绝不允许重新对齐**：核心段必须直接继承 `core_asr.afa`（或 `core_global_matchonly.afa` 子集），列坐标保持不变。
- **模块列绝不允许重新对齐**：模块段继承 Phase 3.8 的模块 MSA。
- **只有 linker 区域允许自由对齐**：linker 使用 MAFFT E-INS-i 在亚型内自由对齐。
- 仅对**同一结构架构**的序列构建 full-length MSA（例如 Type Iβ-ACT），避免强行缝到不匹配架构上。

### 3.9.1 定义亚型与架构子集

```bash
# 以 Type Iβ + ACT 为例
python scripts/select_sequences.py \
  --presence_table results/03_msa_modules/module_presence_absence_strict.tsv \
  --require_core 1 \
  --require_module ACT \
  --output results/03_msa_full/Ib_ACT.ids
```

### 3.9.2 提取 linker 片段

```bash
mkdir -p results/03_msa_full

python scripts/extract_linkers.py \
  --full_length_fasta results/02_qc/nr80_Ib.fasta \
  --seq_ids results/03_msa_full/Ib_ACT.ids \
  --core_coords results/03_msa_core/core_domain_coords.tsv \
  --module_coords results/03_msa_modules/ACT_domain_domain_coords.tsv \
  --output results/03_msa_full/Ib_ACT_linkers.fasta
```

### 3.9.3 亚型内 linker 自由对齐

```bash
mafft --genafpair --maxiterate 1000 --ep 0 --thread 20 \
  results/03_msa_full/Ib_ACT_linkers.fasta \
  > results/03_msa_full/Ib_ACT_linkers_einsi.afa
```

### 3.9.4 拼接 core + linker + module → full-length MSA

```bash
python scripts/stitch_full_length_msa.py \
  --seq_ids results/03_msa_full/Ib_ACT.ids \
  --core_msa results/03_msa_core/core_asr.afa \
  --linker_msa results/03_msa_full/Ib_ACT_linkers_einsi.afa \
  --module_name ACT \
  --module_msa results/03_msa_modules/ACT_msa.afa \
  --output results/03_msa_full/msa_full_Ib_v4.afa \
  --emit_column_map results/03_msa_full/msa_full_Ib_column_map.tsv \
  --assert_core_columns_unchanged
```

**产出：**

- `results/03_msa_full/msa_full_Ib_v4.afa`：Phase 4.4 嵌套 ASR 的输入
- `results/03_msa_full/msa_full_Ib_column_map.tsv`：列坐标 →（core/linker/module）分段映射表

**QC2b（强制断言）：**

1. full-length MSA 的 core 段列数与 `core_asr.afa` 完全一致
2. 任取 3 条代表序列，core 段字符串逐列完全相等
3. linker 段膨胀不会"泄漏"进 core 或 module 段（通过 column_map 检查）

---

# 第四步：系统发育与分层 ASR [待执行]

> **V5.0 注：本步保持完整，不缩减。这是论文演化叙事的时间轴骨架。**

## 4.1 核心树推断与外群定根（含 LBA 抗性检验）

```bash
# 1) KDOPS 外群核心域对齐后并入（仍走 Stockholm → esl-alimask → AFA）
hmmalign --trim --mapali results/03_msa_core/skeleton_core_aa.fa \
  -o results/04_phylogeny_asr/kdops_core_raw.sto \
  results/03_msa_core/core_global.hmm data/seeds/kdops_outgroup_core.fasta
esl-alimask --rf-is-mask results/04_phylogeny_asr/kdops_core_raw.sto \
  > results/04_phylogeny_asr/kdops_core_matchonly.sto
esl-reformat afa results/04_phylogeny_asr/kdops_core_matchonly.sto \
  | seqkit seq --upper-case > results/04_phylogeny_asr/kdops_core_aligned.afa

# 2) 合并核心 + 外群
python scripts/merge_alignments.py \
  --core results/03_msa_core/core_tree.afa \
  --outgroup results/04_phylogeny_asr/kdops_core_aligned.afa \
  --output results/04_phylogeny_asr/core_with_outgroup.afa

# 3) Baseline 全局树
iqtree -s results/04_phylogeny_asr/core_with_outgroup.afa \
  -m MFP -B 1000 -T AUTO \
  -o "KDOPS_P0A715,KDOPS_Q9ZFK4,KDOPS_O66496" \
  --prefix results/04_phylogeny_asr/CoreTree_rooted_MFP

# 4) LBA 抗性：位点异质性模型（至少跑一个）
iqtree -s results/04_phylogeny_asr/core_with_outgroup.afa \
  -m LG+C20+F+G -B 1000 -T AUTO \
  -o "KDOPS_P0A715,KDOPS_Q9ZFK4,KDOPS_O66496" \
  --prefix results/04_phylogeny_asr/CoreTree_rooted_LGC20

# 可选：交叉验证
iqtree -s results/04_phylogeny_asr/core_with_outgroup.afa \
  -m EX_EHO+F+G -B 1000 -T AUTO \
  -o "KDOPS_P0A715,KDOPS_Q9ZFK4,KDOPS_O66496" \
  --prefix results/04_phylogeny_asr/CoreTree_rooted_EXEHO
```

**根稳定性检验（QC3，强制）：**

- **外群子集敏感性**：更换 KDOPS 代表集（至少 2 套）重新定根，根位置必须一致。
- **模型敏感性**：MFP 与 LG+C20 的三亚型分化极性必须一致；不一致则声明"根不鲁棒"并在下游做敏感性分析。
- **中点定根一致性**：对去外群树做 midpoint rooting 交叉检查。

```bash
python scripts/qc_root_stability.py \
  --tree_mfp results/04_phylogeny_asr/CoreTree_rooted_MFP.treefile \
  --tree_c20 results/04_phylogeny_asr/CoreTree_rooted_LGC20.treefile \
  --outgroup_prefix KDOPS_ \
  --output results/04_phylogeny_asr/QC3_root_stability.md
```

## 4.2 结构系统发育交叉验证——AA 树 vs 3Di 树

```bash
# AA 树（标准模型选择）
iqtree -s results/03_msa_core/skeleton_refined_aa.fa \
  -m MFP -B 1000 -T AUTO \
  --prefix results/04_phylogeny_asr/SkeletonTree_AA

# 3Di 树（3Di 专用模型）
iqtree -s results/03_msa_core/skeleton_core_3di.fa \
  -m MFP -mset Q.3Di.AF,Q.3Di.LLM,GTR20 \
  -mdef meta/models/Q_3Di_models.nex \
  -B 1000 -T AUTO \
  --prefix results/04_phylogeny_asr/SkeletonTree_3Di

# 拓扑对比
python scripts/compare_trees.py \
  --tree1 results/04_phylogeny_asr/SkeletonTree_AA.treefile \
  --tree2 results/04_phylogeny_asr/SkeletonTree_3Di.treefile \
  --output results/04_phylogeny_asr/tree_comparison.md
```

**解释边界：** 3Di 树仅在骨架代表（46 序列）上构建，角色严格限定为验证 AA 树的关键深分支（三亚型分化节点、KDOPS 分离位置）。

## 4.3 核心氨基酸 ASR

> **⚠ 树–比对 tip 集一致性（硬约束）：** Phase 4.1 的定根树 `CoreTree_rooted_MFP.treefile` 包含 KDOPS 外群 tips，而 `core_asr.afa` 仅含 ingroup。IQ-TREE `-te` 要求树与比对的 tip 集严格一致，不匹配会报错或产生不可控的隐式行为。
>
> **选定方案：Prune 外群，保留根位置。** 从定根树中剪掉所有 `KDOPS_*` tips，产出仅含 ingroup 的有根树 `CoreTree_rooted_ingroup.treefile`。根位置（即 ingroup MRCA 上方的根节点）在 pruning 后自然保持。这比"把 KDOPS 加入 core_asr.afa"更干净，因为外群序列在 ASR 比对中只会引入噪声（极远缘序列的 gap 和替换污染祖先重建）。

```bash
# Prune KDOPS 外群，保留根位置
python scripts/prune_tree.py \
  --input results/04_phylogeny_asr/CoreTree_rooted_MFP.treefile \
  --remove_prefix KDOPS_ \
  --output results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
  --assert_rooted

# 断言 tip 集一致性
python scripts/assert_tip_match.py \
  --tree results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
  --msa results/03_msa_core/core_asr.afa \
  --assert_identical
# 断言失败 → 终止，回溯定位不匹配的 tip IDs

# ASR：使用 pruned ingroup tree
iqtree -s results/03_msa_core/core_asr.afa \
  -te results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
  -m MFP -asr \
  --prefix results/04_phylogeny_asr/ASR_core
```

> **同理适用于 LG+C20 定根树：** 若后续需要在 site-heterogeneous tree 上做 ASR 敏感性分析，同样先 prune 再跑。

## 4.4 嵌套 ASR（亚型内全长比对上的局部 ASR）

> **输入约束：** `msa_full_Ib_v4.afa` 必须来自 Phase 3.9 缝合比对。

```bash
iqtree -s results/03_msa_full/msa_full_Ib_v4.afa \
  -g results/04_phylogeny_asr/Subtree_Ib_topology.nwk \
  -m MFP -asr \
  --prefix results/04_phylogeny_asr/ASR_Ib_local
```

## 4.5 Gap 祖先态重建——分层策略

**方案 A（推荐，用于模块级）：** 模块存在/缺失 → PastML 离散性状 ASR（Phase 4.6）

**方案 B（核心列级 gap）：** 连续 gap 列聚合为 indel 区段 → PastML

```bash
python scripts/aggregate_gap_blocks.py \
  --input results/03_msa_core/core_asr.afa \
  --min_block_size 3 \
  --output results/04_phylogeny_asr/gap_blocks.tsv

pastml --tree results/04_phylogeny_asr/CoreTree_rooted_MFP.treefile \
  --data results/04_phylogeny_asr/gap_blocks.tsv \
  --out_data results/04_phylogeny_asr/gap_block_asr/
```

**方案 C（可选进阶）：** indelMaP / ARPIP 等 indel-aware ASR，用于需要精确 indel 事件重建的场景。

## 4.6 模块获得/丢失的离散性状 ASR

```bash
for version in strict relaxed; do
  pastml --tree results/04_phylogeny_asr/CoreTree_rooted_MFP.treefile \
    --data results/03_msa_modules/module_presence_absence_${version}.tsv \
    --out_data results/04_phylogeny_asr/module_origins_${version}/ \
    --html_compressed results/04_phylogeny_asr/module_origins_${version}_viz.html
done
```

**敏感性分析（强制）：** strict 与 relaxed 两版的模块起源节点是否一致 → 写入 QC3。

## 4.7 AltAll 系综采样

使用 maxaltall R 包处理 IQ-TREE `.state` 文件。阈值：PP₁ = 0.80, PP₂ = 0.20。

## 4.8 模块层局部 ASR

```bash
for module in ACT N_ext alpha2beta3; do
  iqtree -s results/03_msa_modules/${module}_msa.afa \
    -g results/04_phylogeny_asr/Subtree_${module}_topology.nwk \
    -m MFP -asr \
    --prefix results/04_phylogeny_asr/ASR_${module}
done
```

### QC3（强制产出）

`results/04_phylogeny_asr/qc_phylogeny_asr.md`：

- 关键节点 bootstrap 支持
- 外群定根稳定性
- AA 树 vs 3Di 树关键节点一致性
- 祖先序列平均 PP 分布
- 模块获得节点在不同定根方案/strict-relaxed 下是否稳定

---

# 第五步：关键祖先节点的结构验证 [待执行，V5.0 精简版]

> **V5.0 范围调整：** 从 V4.1 的"8 节点 × 2×2 因子 × 3 rep 全因子 MD"大幅缩减为"2–4 个关键祖先节点 × Apo × 有限验证性 MD"。目标从"建立动力学通讯网络"降级为"验证祖先结构在其功能装配体下的可行性与口袋拓扑演化"。
>
> **⚠ V5.0 硬约束：不默认四聚体。** 装配体状态（dimer / tetramer / 其他）本身可能是变构演化的性状之一（文献中 Type II 存在 dimeric 实例）。每个祖先节点的功能装配体必须通过 Phase 5.0 显式判定，不可写死。

## 5.0 装配体判定（Assembly Adjudication）⚠ [CHECK-08]

> **前置于所有结构预测与 MD。** 对每个目标 clade / 祖先节点，先确定"最合理的功能装配体"是 dimer 还是 tetramer（或其他），再决定送入 AF3 的拷贝数。

**判定流程：**

1. **文献/结构注释扫描：** 检查该节点所在 clade 的现存代表在 PDB / UniProt / PDBe PISA 中的生物学装配体注释。统计 clade 内 dimer vs tetramer 的比例。
2. **AF3 多装配体平行预测：** 对每个祖先序列同时提交 2-copy（dimer）和 4-copy（tetramer）两套 AF3 任务（均为 Apo）。
3. **界面质量评分：** 比较两种装配体的 ipTM、跨链 PAE、PISA 界面面积与 ΔG_dissoc。
4. **判定规则：**
   - 若 clade 内文献/注释一致指向同一装配体 → 采用文献共识（但仍做另一种装配体的 AF3 对照，不做 MD）
   - 若文献不一致或注释缺失 → 以 AF3 ipTM + PISA ΔG_dissoc 较优者为"native oligomer hypothesis"
   - 在任何情况下，判定结果和证据链写入 `results/05_struct_valid/assembly_adjudication.tsv`

```bash
python scripts/adjudicate_assembly.py \
  --ancestor_seqs results/05_struct_valid/ancestor_candidates.fasta \
  --clade_pdb_annotations results/05_struct_valid/clade_assembly_survey.tsv \
  --af3_dimer_dir results/05_struct_valid/af3_dimer/ \
  --af3_tetramer_dir results/05_struct_valid/af3_tetramer/ \
  --output results/05_struct_valid/assembly_adjudication.tsv
```

**产出：** 每个祖先节点一行：`node_id | native_oligomer | evidence_type | dimer_ipTM | tetramer_ipTM | dimer_PISA_dG | tetramer_PISA_dG | decision_rationale`

## 5.1 候选祖先节点选择（2–4 个）

从 Phase 4.6 的模块获得事件中选择**最关键的 2–4 个节点**：

| 节点类型 | 数量 | 选择标准 |
|---------|------|---------|
| Pre-gain（模块获得前） | 1–2 | 核心树上模块性状 0→1 的前一节点，bootstrap ≥ 70 |
| Post-gain（模块获得后） | 1–2 | 模块性状为 1 的最早节点，bootstrap ≥ 70 |

**选择优先级：**
1. 首选 Type Iβ 的 ACT 获得事件前后节点（文献最丰富、功能验证最充分的变构模块）
2. 若 ACT 获得节点 bootstrap < 70 或 ASR PP 太低，换选 CM 域或 N_ext

> **⚠ [CHECK-02] 祖先序列来源：** 送进 AF3 的全长祖先序列必须来自 **Phase 4.4 的亚型内全长嵌套 ASR**。不可将核心 ASR 与模块 ASR 片段硬拼接。

## 5.2 Apo 结构预测（装配体由 5.0 判定）

**仅做 Apo 预测（V5.0 硬约束）：** 不做 Holo 预测，消除"现代配体幻觉"风险。

```bash
# AlphaFold3：拷贝数由 Phase 5.0 判定结果决定（N_copies 读自 assembly_adjudication.tsv）
# 对每个候选祖先同时使用 ESMFold 交叉验证（单体 pLDDT）

# FoldMason 标记预测可信度
foldmason easy-msa ancestor_af3.pdb known_pdbs/*.pdb \
  results/05_struct_valid/struct_comparison tmp_struct --report-mode 1
```

## 5.2b 结构 QC 门控（强制）⚠ [CHECK-05]

```
门控标准：
1. AF3 ipTM ≥ 0.6 + 跨链 PAE < 15 Å
   → 不通过：换用 ML-best 序列重新预测
2. GROMACS 能量最小化 + 分阶平衡（NVT→NPT→逐步撤约束）
3. 10 ns 短轨迹筛选：RMSD 稳定 + 界面面积维持
   → 不通过：标记"结构不可信"，不进入正式 MD
```

## 5.3 验证性 MD（Apo，装配体由 5.0 判定）

> **V5.0 设计：** 不做全因子实验。目标是验证"祖先在其 native oligomer（由 5.0 判定）的 Apo 条件下是否动力学稳定"，并初步比较 Pre-gain 与 Post-gain 祖先的口袋拓扑差异。

```
实验矩阵（每个祖先节点，仅 Apo）：
  native-oligomer-apo:  ≥200 ns × 2 replicas（主分析）
  （可选）alternative-oligomer-apo: 10–20 ns × 1 replica
    → 仅用于排除性测试（证明"该装配体不稳定/不合理"），不作为机制主证据

资源预算：
  祖先节点数: ≤ 4
  AF3 预测: ≤ 8 结构（每个节点 dimer + tetramer 两套，全部 Apo）
  MD 主分析: native-oligomer × ≥200 ns × 2 rep × ≤4 节点
  MD 排除测试（可选）: alt-oligomer × 10–20 ns × 1 rep × ≤4 节点
  MD 总 GPU·时: ≤ 25 GPU-days
```

## 5.4 有限动力学读出

- RMSD/RMSF 稳定性验证
- 关键口袋/界面面积是否维持
- Pre-gain vs Post-gain 的口袋拓扑定性比较（Fpocket）
- Alternative oligomer 的界面崩溃/稳定性作为排除证据
- **不做完整 DCCM/网络流分析**（留待未来全因子 MD）

### QC4（强制产出）

`results/05_struct_valid/qc_struct_validation.md`：

- 装配体判定证据与结论（引用 `assembly_adjudication.tsv`）
- AF3 置信度（ipTM / PAE，native vs alternative oligomer 对比）
- 轨迹收敛性（RMSD 平台期）
- 界面面积稳定性（native oligomer 主分析 + alternative oligomer 排除测试）
- Pre-gain vs Post-gain 口袋拓扑定性比较
- 低置信区段标注
- **显式声明**：若 native oligomer 判定在某节点存在不确定性，在论文中报告为"oligomeric state ambiguous at this node"

---

# 第六步：核心层共进化分析（DCA） [待执行，V5.0 聚焦版]

> **V5.0 范围调整：** 仅保留核心层 DCA（Meff/L > 5，信息量充足）。模块 DCA、联合跨域 DCA、Meff 匹配比较全部移入"可选探索"。

## 6.1 核心层 DCA 输入准备

```bash
python scripts/prepare_dca_input.py \
  --input results/03_msa_core/core_asr.afa \
  --gap_col_max 0.50 \
  --gap_row_max 0.50 \
  --meff_min 3.0 \
  --output results/06_dca/core_dca.afa \
  --stats results/06_dca/core_dca_stats.tsv
```

**预期 Meff/L：** 基于 9,393 条 × 472 列的核心 MSA，预期 Meff/L 远超 5.0 门控。

## 6.2 核心层 DCA 执行

```bash
plmc -o results/06_dca/core.params \
  -c results/06_dca/core_couplings.txt \
  -le 16.0 -lh 0.01 -m 100 -t 0.8 \
  -f <focus_seq_id> -g \
  results/06_dca/core_dca.afa
```

## 6.3 显著性评估与接触验证

```bash
python scripts/dca_significance.py \
  --couplings results/06_dca/core_couplings.txt \
  --method top_L \
  --contact_pdb data/structures/panel_dah7ps/PDB-1KFL.cif \
  --contact_cutoff 8.0 \
  --min_separation 5 \
  --output results/06_dca/core_significant_couplings.tsv
```

**分析重点：**

1. **核心层 top-L 耦联的接触精度**：用实验结构（1KFL/1RZM/3NV8）验证 DCA 预测的残基对是否在结构中形成物理接触（Cβ–Cβ < 8 Å）。
2. **功能位点富集**：检查 top 耦联对是否富集在催化位点、亚基界面、已知变构铰链等功能区域。
3. **跨亚型保守 vs 亚型特异耦联**：将核心 DCA 耦联按亚型归属分层，检查哪些耦联在所有亚型中保守（→ 核心架构约束），哪些仅在特定亚型中显著（→ 可能与变构模块整合相关）。

### QC5（强制产出）

`results/06_dca/qc_core_dca.md`：

- Meff / L / Meff/L
- Top-L 接触精度（PPV）
- 催化位点 / 界面 / 铰链区域的耦联富集度
- 与现有文献已知的变构通讯路径的一致性

## 6.4 可选探索（不进入论文主线结论）

以下分析因有效序列深度不足（Meff/L < 3.0）或缺乏完整 MD 网络而不作为论文主线定量证据，但可作为补充材料或 Discussion 中的定性支持：

**6.4.1 模块层 DCA（探索性）：**

> ACT strict 仅 47 序列（L=142），Meff/L 预计 ≈ 0.2–0.3，远低于 3.0 门控。如执行，结果仅写入补充材料，明确标注"探索性，深度不足"。

**6.4.2 联合跨域 DCA（探索性）：**

> 同一架构亚型的 core+linker+module 联合比对 DCA，用于探测 core↔module 跨域耦联。受限于模块携带子集的序列数。

**6.4.3 "模块获得前后耦联变化"比较（探索性）：**

> `core_with_ACT` vs `core_without_ACT` 的 Meff 匹配下采样比较。需要充足的两组 Meff 才有意义。

---

# 第七步：论文写作蓝图与跨证据一致性展望 [待执行]

> **V5.0 新增。** 在 V4.1 中 ICDC 是"第六步主线产出"；在 V5.0 中降级为 Discussion 展望。本步定义论文的主线叙事、图表清单，以及 ICDC 在 Discussion 中的定位。

## 7.1 论文主线叙事

**叙事主轴：** "DAH7PS 的变构调控模块在保守 TIM-barrel 核心之上被独立募集——它们何时出现、如何在序列-结构层面被整合。"

**证据链条（按论文逻辑顺序）：**

1. **结构感知的核心 MSA 消除了暮光区膨胀**（Phase 3）→ 建立了跨亚型可靠的同源列坐标系
2. **定根系统发育确立时间轴**（Phase 4.1–4.2）→ 三亚型分化极性、KDOPS 外群验证
3. **模块离散性状 ASR 锁定获得/丢失时间节点**（Phase 4.6）→ 变构元件起源的演化时间线
4. **核心层 DCA 揭示保守架构约束与变构相关耦联**（Phase 6）→ 哪些位点对在核心中长期共进化，哪些与模块整合相关
5. **关键祖先节点的结构验证**（Phase 5）→ 装配体判定、祖先 native-oligomer 可行性、Pre-gain vs Post-gain 口袋拓扑变化

## 7.2 核心图表清单（预定）

| 图号 | 内容 | 数据来源 |
|------|------|---------|
| Fig 1 | DAH7PS 家族架构概览（三亚型 + 模块分布） | Phase 3.8 模块注释 |
| Fig 2 | FoldMason 骨架 vs V3.1 旧比对的 LDDT 对比 | Phase 3.4 / QC2 |
| Fig 3 | 定根核心树 + 模块获得/丢失时间线 | Phase 4.1 + 4.6 |
| Fig 4 | AA 树 vs 3Di 树关键节点一致性 | Phase 4.2 |
| Fig 5 | 核心层 DCA 耦联网络（投影到 1KFL 结构） | Phase 6 |
| Fig 6 | 关键祖先 Apo native-oligomer 结构对比（Pre-gain vs Post-gain） | Phase 5 |
| Fig S1–Sn | 补充：Gap ASR、模块 MSA 详情、DCA 接触精度等 | 各 Phase |

## 7.3 ICDC 在 Discussion 中的定位

> **写法指南：** Discussion 倒数第二段设置一个"跨证据一致性展望（Cross-Evidence Convergence Outlook）"小节。

**内容要点：**

1. 指出核心 DCA 耦联 × 祖先结构口袋拓扑 × 有限 MD 稳定性之间已观察到的定性一致性信号（如：top DCA 耦联对在 Pre-gain 祖先结构中距离较远，在 Post-gain 祖先结构中缩短 → 提示模块整合伴随着界面耦联强化）。
2. 明确声明：这些一致性观察是定性的，尚未达到正式 ICDC 融合所需的完整 MD 动力学网络支持。
3. 提出未来方向：全因子 MD（Apo/Holo × 单体/多聚体 × 多节点）→ 完整 DCCM/网络流 → 正式 DCA×DCCM ICDC 融合。

---

# 工程风险 Checklist（V5.0 精简版）

> V5.0 保留与当前范围相关的 CHECK 项，去除仅与 Holo-MD 或 ICDC 融合相关的项。

### [CHECK-01] HMM 域提取的"边缘破碎"效应

| | |
|---|---|
| **关联 Phase** | 3.6（`extract_core_domains.py`）— ✅ 已执行并验证 |
| **状态** | pad=20 已实施。端部 gap 死区确认为真实生物学变异，非提取截断。 |

### [CHECK-02] 祖先序列的"弗兰肯斯坦缝合"

| | |
|---|---|
| **关联 Phase** | 5.1（候选祖先集合定义）、5.2（结构预测） |
| **缓解** | 全长祖先序列必须来自 Phase 4.4 亚型内全长嵌套 ASR。不可手工硬拼接。 |
| **失败处理** | 若 linker PP < 0.3，报告"linker 不可信"并降权。 |

### [CHECK-03] DCA 的 Meff/L 门控

| | |
|---|---|
| **关联 Phase** | 6.1（DCA 输入准备） |
| **V5.0 调整** | 仅核心层 DCA 需要通过门控（预期 Meff/L >> 5）。模块 DCA 已移入可选探索，不影响主线。 |
| **参数** | `meta/params.json → dca.meff_min_main: 3.0; dca.meff_ideal: 5.0` |

### [CHECK-05] AlphaFold3 祖先多聚体的"立体位阻爆炸"

| | |
|---|---|
| **关联 Phase** | 5.2b（MD 前结构 QC 门控） |
| **缓解** | 三级门控不变：ipTM ≥ 0.6 → 能量最小化 + 分阶平衡 → 10 ns 短轨迹。 |
| **失败处理** | 不通过 → 换 ML-best 序列或标记"结构不可信"。 |

### [CHECK-06] Type II 内插片导致的 HMM Hit 断裂

| | |
|---|---|
| **关联 Phase** | 3.6（`extract_core_domains.py`）— ✅ 已执行并验证 |
| **状态** | hit stitching 已实施（1,069 条 = 11.4%），rescued 2,093 条 Type II 在 Phase 2.1。 |

### [CHECK-08] 装配体状态不可默认——Assembly Adjudication（V5.0 新增）

| | |
|---|---|
| **关联 Phase** | 5.0（装配体判定）、5.2（结构预测）、5.3（MD） |
| **风险** | DAH7PS 家族内装配体状态并非统一的四聚体：Type II 存在 dimeric 实例，部分深分支祖先的功能装配体未知。若将 tetramer 写死送入 AF3 和 MD，可能模拟一个不属于该节点的装配体，产生"四聚体界面崩溃"的假阳性或"虚假界面稳定"的假阴性。 |
| **缓解** | Phase 5.0 强制做 Assembly Adjudication：文献扫描 + dimer/tetramer 平行 AF3 预测 + PISA 界面评分，产出 `assembly_adjudication.tsv`。所有下游结构/MD 步骤读取该表决定拷贝数。 |
| **失败处理** | 若两种装配体的 AF3 指标无法区分（ipTM 差异 < 0.05 且 PISA ΔG 方向不一致）→ 标记为"oligomeric state ambiguous"，在论文中显式声明不确定性，对两种装配体均做短 MD 排除测试但不将任一作为主线机制证据。 |
| **参数记录** | `meta/params.json → assembly.ipTM_diff_threshold: 0.05` |

---

# 里程碑与交付物（V5.0 修订版）

| 里程碑 | 交付物 | 验收标准 | 状态 |
|--------|--------|---------|------|
| M1 | 结构骨架 + 核心列定义 | LDDT 报告 + 核心列 521（目标 400–600） | ✅ 完成 |
| M2 | 全量 core_global MSA | 9,393 seqs × 521 cols + 结构复核 LDDT=0.2638 | ✅ 完成 |
| M2b | 双版本修剪 + 模块注释 | core_tree(436) + core_asr(472) + 5 模块 MSA | ✅ 完成 |
| M2c | 亚型内全长缝合比对 | Phase 3.9 core 段列数一致性断言通过 | ⬜ 待执行 |
| M3 | 外群定根核心树 | 根位置对外群抽样鲁棒 + AA/3Di 一致 | ⬜ 待执行 |
| M4 | 模块性状 ASR + Gap ASR | 候选起源节点可复现 + strict/relaxed 一致 | ⬜ 待执行 |
| M5 | 核心层 DCA | Meff/L ≥ 5 + top-L PPV 合理 + 功能位点富集 | ⬜ 待执行 |
| M6 | 关键祖先 Apo 结构验证 | 2–4 节点 AF3 ipTM ≥ 0.6 + MD ≥200 ns 稳定 | ⬜ 待执行 |
| M7 | 论文初稿 | 主线叙事闭环 + Fig 1–6 齐备 | ⬜ 待执行 |

---

# 参考文献

1. Cross PJ, Dobson RCJ, Patchett ML, Parker EJ. The diversity of allosteric controls at the gateway to aromatic amino acid biosynthesis. *Biochem J*. 2013;451(2):127–136.
2. Eick GN, Bridgham JT, Anderson DP, Harms MJ, Thornton JW. Robustness of reconstructed ancestral protein functions to statistical uncertainty. *Mol Biol Evol*. 2017;34(2):247–261.
3. Gilchrist CLM, Mirdita M, Steinegger M. Multiple protein structure alignment at scale with FoldMason. *Science*. 2026; doi:10.1126/science.ads6733.
4. Hill MD, Shafiei R, Muñoz-Rojas T, Harms MJ. Topiary: pruning the manual labor from ancestral sequence reconstruction. *Protein Sci*. 2023;32(1):e4519.
5. Katoh K, Kuma K, Toh H, Miyata T. MAFFT version 5: improvement in accuracy of multiple sequence alignment. *Nucleic Acids Res*. 2005;33(2):511–518.
6. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7. *Mol Biol Evol*. 2013;30(4):772–780.
7. Lang EJM, Cross PJ, Sherwood OC, Parker EJ. Quaternary structure is an essential component that contributes to the sophisticated allosteric regulation mechanism in a key enzyme from *Mycobacterium tuberculosis*. *PLoS ONE*. 2016;12(1):e0169601.
8. Mariani V, Biasini M, Barbato A, Schwede T. lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests. *Bioinformatics*. 2013;29(21):2722–2728.
9. Minh BQ, Schmidt HA, Chernomor O, et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol*. 2020;37(5):1530–1534.
10. Muñiz-Trejo R, Park Y, Thornton JW. Robustness of ancestral sequence reconstruction to among-site and among-lineage evolutionary heterogeneity. *Mol Biol Evol*. 2025;42(4):msaf084.
11. Steinegger M, Söding J. Clustering huge protein sequence sets in linear time. *Nat Commun*. 2018;9:2542.
12. Steenwyk JL, Buida TJ III, Li Y, Shen XX, Rokas A. ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. *PLoS Biol*. 2020;18(12):e3001007.
13. van Kempen M, Kim SS, Tumescheit C, et al. Fast and accurate protein structure search with Foldseek. *Nat Biotechnol*. 2024;42:243–246.
14. Webby CJ, Baker HM, Lott JS, Baker EN, Parker EJ. Synergistic allostery, a sophisticated regulatory network for the control of aromatic amino acid biosynthesis in *Mycobacterium tuberculosis*. *J Biol Chem*. 2010;285(40):30567–30576.
15. Ishikawa SA, Zhukova A, Iwasaki W, Gascuel O. A fast likelihood method to reconstruct and visualize ancestral scenarios. *Mol Biol Evol*. 2019;36(9):2069–2085.
16. Jiao W, Lang EJM, Bai Y, Fan Y, Parker EJ. Diverse allosteric componentry and mechanisms control entry into aromatic amino acid biosynthesis. *Curr Opin Struct Biol*. 2020;65:159–167.
17. Seffernick JT, Pardo-Avila F, Shen H, et al. A general substitution matrix for structural phylogenetics. *Mol Biol Evol*. 2025;42(6):msaf124.
18. Ekeberg M, Lövkvist C, Lan Y, Weigt M, Aurell E. Improved contact prediction in proteins: using pseudolikelihoods to infer Potts models. *Phys Rev E*. 2013;87(1):012707.
19. Hopf TA, Green AG, Schuber B, et al. The EVcouplings Python framework for coevolutionary sequence analysis. *Bioinformatics*. 2019;35(9):1582–1584.
