# DAH7PS 变构起源的演化动力学重建：V5.1 标准作业程序（SOP rev6）

**A FoldMason-Centred, Core–Module Decoupling Framework for Reconstructing the Evolutionary Origin of Allostery in DAH7PS**  
**版本说明：基于 Phase 4.3 阶段审阅后修订，用于替换当前 `PLAN.md`。**

---

## 摘要

DAH7PS（3-deoxy-D-arabino-heptulosonate 7-phosphate synthase）在保守的 `(β/α)8` TIM-barrel 催化核心之上，于不同亚型中独立演化出多样的变构调控元件——Type Iα 的 β 发夹插入、Type Iβ 的 ACT/铁氧还蛋白样调控域、Type II 的 N 端延伸与 α2β3 插片。V3.1 以来的最大方法学风险一直是：在跨亚型远缘同源（序列一致性进入暮光区）条件下，用全长单一 MSA 同时驱动系统发育、ASR、DCA 和结构解释，会放大比对膨胀与错误同源关系，进而污染整条证据链。

V5.0 已经把主线收缩为“结构感知核心 MSA → rooted phylogeny + 模块 trait ASR → core DCA → 少量关键祖先的结构验证”，方向正确；但截至 Phase 4.3，项目的主要瓶颈已从 **tip set consistency** 转移为 **root robustness**。因此，V5.1 的核心调整不是推翻现有路线，而是：

1. **确认 Phase 3.9 已完成**，并把其纳入正式计划状态，而不是继续作为待执行项。  
2. **把 QC3 从普通 QC 升级为真正的 gate**：在 root robustness 评估完成前，可以继续跑工作树上的 ASR 与 DCA 准备，但不能锁定最深层历史叙事，也不能锁定 Phase 5 的主线祖先节点。  
3. **引入多 root scenario 框架**：KDOPS+MFP、KDOPS+site-heterogeneous model、outgroup-free midpoint/MAD（可选再加 reduced-set nonreversible/rootstrap）。  
4. **区分“working tree”与“narrative tree”**：当前 pruned ingroup rooted tree 可以继续用于计算，但只有通过 QC3 的 root-robust 事件才进入主文强结论。  
5. **收紧 Phase 5 节点选择标准**：不再用“bootstrap ≥ 70”作为主线节点准入；改为支持度、root stability、strict/relaxed stability、full-length ASR 可解释性四重共同门控。  
6. **新增文档与指标单一真源（single source of truth）**，避免 README / TASKS / PLAN / log 之间数字和状态漂移。  
7. **论文写作改为 root-robust narrative 优先**：主文强调“保守核心上的重复模块募集”，root-dependent 的最深层极性放入补充或 Discussion 条件性阐释。

V5.1 的目标是：**不停工、不断线、不硬拗。** 继续推进 Phase 4.4、4.6、6.1 等低耦合步骤，同时把最容易被审稿人质疑的 deep root 问题显式纳入分析设计。

---

## V5.1 相对 V5.0 的关键修改

| 编号 | 修改 | V5.1 处理 |
|---|---|---|
| 1 | Phase 3.9 在 PLAN 中仍写“待执行” | 改为 **已完成**，并列出产物与 QC2b 断言 |
| 2 | QC3 只是一般性 QC | 升级为 **强制 gate**，决定是否允许锁定历史叙事与 Phase 5 节点 |
| 3 | 根只依赖 KDOPS outgroup + 两套模型 | 扩展为 **多 root scenario**：KDOPS+MFP、KDOPS+LG+C20、midpoint/MAD，条件允许再加 nonreversible/rootstrap |
| 4 | 4.2 AA vs 3Di tree 更像补充检查 | 提升为 **高优先级正交证据**，并纳入 QC3 判定 |
| 5 | 4.3 解决 tip mismatch 后默认可直接进入后续强解释 | 改为：**可继续计算，不可直接写死最深层历史解释** |
| 6 | Phase 5 节点选择阈值偏松（bootstrap ≥ 70） | 改为：**UFBoot ≥ 95 为主阈值；若可补跑则 SH-aLRT ≥ 80 为优先 corroboration** |
| 7 | 文档中的状态与数字可能漂移 | 新增 `results/meta/metrics_manifest.tsv` 与 `results/meta/progress_snapshot.md` |
| 8 | 主文叙事默认单一时间轴 | 改为 **root-robust 主叙事 + root-sensitive 条件性解释** |

---

## 当前状态快照（以 Phase 4.3 审阅时点为准）

### 已完成
- Phase 0：环境、参数文件、模型文件记录、目录规范。
- Phase 1：全长挖掘、KDOPS 反向过滤、扩充种子重建 HMM。
- Phase 2：长度/覆盖度三箱过滤、Type II hit stitching、NR80/seed60 去冗余。
- Phase 3.1–3.8：结构面板、FoldMason 骨架、LDDT 核心列、全量核心映射、双版本 trim、模块注释。
- **Phase 3.9：Profile-anchored stitching 已完成**：
  - `results/03_msa_full/Ib_ACT.ids`
  - `results/03_msa_full/Ib_ACT_linkers.fasta`
  - `results/03_msa_full/Ib_ACT_linkers_einsi.afa`
  - `results/03_msa_full/msa_full_Ib_v4.afa`
  - `results/03_msa_full/msa_full_Ib_column_map.tsv`
  - QC2b：core 段 472 列与 `core_asr.afa` 完全一致。

### 已进入进行中
- Phase 4.1：
  - KDOPS 外群并入核心比对完成。
  - MFP 树完成（集群，Q.PFAM+F+R10，LogL=-2835499.68）。
  - **LG+C20 抗 LBA 树未完成**（本地 OOM killed，需在集群补跑）。
- Phase 4.3：
  - KDOPS 已从 rooted tree 剪除，得到 ingroup rooted working tree。
  - `assert_tip_match.py` 已通过，`core_asr.afa` 与 pruned ingroup tree 的 tips 严格一致。

### 当前真正瓶颈
- **不是 core MSA 崩坏，也不是 tip mismatch。**
- **而是 rooted history 的稳健性：KDOPS 在 ML 树中的行为提示根位置仍需敏感性分析。**

### 立刻可以并行推进的工作
- 4.4 嵌套 full-length ASR（优先 Iβ-ACT）
- 4.6 strict / relaxed 模块 trait ASR 的数据准备
- 6.1 core DCA 输入准备
- 5.0 Assembly adjudication 的文献/注释扫描准备

### 暂时不能锁死的工作
- 最深层模块获得/丢失的唯一历史箭头
- Phase 5 的主线祖先节点最终名单
- 主文 Fig.3 的最终 rooted history 叙事

---

## 科学问题与工作假设

### 核心问题
DAH7PS 的变构调控元件如何在保守 TIM-barrel 核心之上被重复募集、稳定并与多聚体界面及构象动力学耦合，从而产生可选择的变构表型？

### 工作假设
1. 跨亚型可可靠比较的演化时间轴主要由 **TIM-barrel core** 提供；模块起源体现为 gain/loss 事件与局部序列—结构整合。
2. DAH7PS 家族的 allostery 是 **模块化、可重复招募** 的，而不是只在单一路径上发生一次后保守继承。
3. 核心层 DCA 可以独立提供“长期保守耦联约束”的证据；它不依赖完整 MD 网络才有价值。
4. 根的不确定性主要影响 **历史极性与 earliest event 的解释**，不必然推翻 subtype-local 或 root-robust 的结论。
5. 祖先结构预测与有限 MD 在本研究中只承担 **可行性验证** 角色，不承担“闭环证明动力学通信网络”的任务。

---

## 方法学立场与解释边界

1. **FoldMason 是中枢，不是唯一真理。**  
   核心列定义依然以结构感知 MSA 和 LDDT 为中心，但系统发育解释必须额外通过 root robustness 与 3Di 交叉验证。

2. **Outgroup rooting 不是唯一根。**  
   KDOPS 在数据挖掘阶段仍然成功，但在 Phase 4 中不得自动享有“唯一正确根”的地位。V5.1 必须显式比较多个 rooting scenario。

3. **当前最佳 rooted tree 只是 working tree。**  
   在 QC3 通过前，它适合做 ASR 计算和后续准备，但不适合直接作为论文中 deepest origin 的唯一时间轴。

4. **主线证据只保留信息量足够的层。**  
   Core DCA 进入主线；module/joint DCA 继续作为探索性分析，不抢主线。

5. **结构验证只回答“可不可能”，不回答“全部怎么发生”。**  
   AF3 + 有限 MD 的任务是验证祖先结构与装配体假设是否可行，而不是完整重建 allosteric communication network。

6. **root-robust claim 优先于全局统一故事。**  
   如果最深层极性在不同根方案下变化，则主文只写对多个根方案都成立的结论。

---

## 0. 文档同步与单一真源（立即执行，不阻塞）

> 新增的跨阶段要求。目标：以后所有 Methods、Figure legend、README、TASKS、论文草稿中的数字，只引用同一份 manifest。

### 0.1 新增交付物
- `results/meta/metrics_manifest.tsv`
- `results/meta/progress_snapshot.md`
- `results/04_phylogeny_asr/root_scenarios.tsv`
- `results/04_phylogeny_asr/node_selection_registry.tsv`

### 0.2 `metrics_manifest.tsv` 建议字段
| field | 含义 |
|---|---|
| metric_id | 唯一指标名 |
| value | 数值或字符串 |
| unit | 单位 |
| source_file | 来自哪个结果文件 |
| phase | 所属阶段 |
| status | final / provisional / exploratory |
| updated_on | 日期 |
| note | 备注 |

### 0.3 立即补记的关键指标
- PASS / NR80 / seed60 数量
- core MSA 行列数
- tree trim / ASR trim 列数
- module counts
- `msa_full_Ib_v4.afa` 行列数与 column map
- rooted working tree tip 数
- QC3 root scenario 列表
- DCA 的 Meff / L / Meff/L
- Phase 5 候选节点 registry

### 0.4 规则
- PLAN / TASKS / README / 论文草稿中的数字，一律以 `metrics_manifest.tsv` 为准。
- 若任一文档与 manifest 冲突，先更新 manifest，再回填其他文档。
- `progress_snapshot.md` 每次推进阶段时同步更新，不再只靠 `TASKS.md` 口头描述。

---

# ✅ 第零步：计算环境与项目规范 [已完成]

> 状态：已完成。以下内容保留供审计。  
> V5.1 只新增一个补记要求：把 root-scenario 与 metrics manifest 纳入 `results/meta/` 管理。

### 0.1 环境与版本锁定
- conda / mamba 环境：`dah7ps_v4`
- 已安装：HMMER、MAFFT、IQ-TREE 2、MMseqs2、FoldMason、SeqKit、CD-HIT、ClipKIT、PastML、plmc
- `results/meta/software_versions.tsv` 已存在
- `results/meta/model_files.tsv` 已存在
- `meta/params.json` 已初始化

### 0.2 目录规范（补充）
```text
results/
  03_msa_core/
  03_msa_modules/
  03_msa_full/
  04_phylogeny_asr/
  05_struct_valid/
  06_dca/
  meta/
```

### 0.3 V5.1 新增参数块（写入 `meta/params.json`）
```json
{
  "rooting": {
    "primary_scenarios": ["MFP_KDOPS", "LGC20_KDOPS", "MIDPOINT_INGROUP", "MAD_INGROUP"],
    "support_main_min_ufboot": 95,
    "support_preferred_min_shalrt": 80,
    "trait_stability_min_scenarios": 2
  },
  "manuscript": {
    "claim_tiers": ["root_robust", "root_sensitive", "exploratory"]
  },
  "metrics": {
    "manifest_path": "results/meta/metrics_manifest.tsv"
  }
}
```

---

# ✅ 第一步：全长序列挖掘与 KDOPS 免疫过滤 [已完成]

> 状态：已完成。QC1 通过。  
> V5.1 不改 Phase 1 的结论，也不要求回滚重做。

### 已确认的要点
- 扩充种子后重建 HMM 成功提升 Type II 召回。
- KDOPS 反向过滤在 Phase 1 的角色仍然有效。
- Ia/II 的 overlap 已通过 best-hit competition 清理为互斥集合。
- KDOPS borderline 序列已隔离，不需要回流主集合。

### V5.1 解释边界
- **KDOPS 在 Phase 1 中是成功的负选择工具。**
- **KDOPS 在 Phase 4 中不再自动等价于唯一可信的深根。**

---

# ✅ 第二步：数据驱动 QC 与去冗余 [已完成]

> 状态：已完成。QC2 通过。

### 已完成且继续沿用
- 2.1 长度 + HMM 覆盖度三箱过滤
- Type II multi-domain stitching
- 2.2 NR80 去冗余
- 2.3 seed60
- 2.4 stepping-stone 结构候选筛选

### V5.1 不改动
- 所有 Phase 2 产出继续作为下游标准输入
- 不因当前 root 问题返工序列挖掘或去冗余

---

# ✅ 第三步：结构感知 MSA——核心/模块分治 [已完成]

> 状态：Phase 3.1–3.9 全部完成。  
> 这是当前项目最稳固的一层，也是 V5.1 继续向前推进的基础。

## 3.1–3.8 已完成摘要
- 35 结构面板
- FoldMason 骨架：46 序列
- core columns：521
- `core_global_matchonly.afa`：9,393 × 521
- `core_tree.afa`：436 列
- `core_asr.afa`：472 列
- 5 类模块注释 strict / relaxed 双矩阵
- 模块 MSA：ACT / CM / α2β3 / N_ext / C_tail

## 3.9 Profile-anchored stitching [已完成]
### 已完成产物
- `results/03_msa_full/Ib_ACT.ids`
- `results/03_msa_full/Ib_ACT_linkers.fasta`
- `results/03_msa_full/Ib_ACT_linkers_einsi.afa`
- `results/03_msa_full/msa_full_Ib_v4.afa`
- `results/03_msa_full/msa_full_Ib_column_map.tsv`

### QC2b（已通过）
1. full-length MSA 的 core 段列数与 `core_asr.afa` 完全一致  
2. 代表序列的 core 段逐列完全相等  
3. linker 膨胀未泄漏至 core / module 段  

### V5.1 新解释
- `msa_full_Ib_v4.afa` 已经足够支持 4.4 Iβ-ACT nested ASR。
- 因此，4.4 不再需要等待 3.9。

---

# 第四步：系统发育与分层 ASR [进行中，V5.1 重点修订]

> **这是当前项目的核心阶段。**  
> V5.1 的原则是：**计算可以继续，叙事必须分层。**

---

## 4.0 Phase 4 的三条总规则

### 规则 A：working tree 与 narrative tree 分离
- `CoreTree_rooted_ingroup.treefile`（当前默认 working tree）可以用于继续计算。
- 但只有经过 QC3 判定为 **root-robust** 的事件与节点，才能进入主文主结论。

### 规则 B：多个 root scenario 并行管理
V5.1 将所有 rooted 分析统一纳入 scenario 框架，至少包含：

| Scenario ID | 描述 | 角色 |
|---|---|---|
| S1 | KDOPS + MFP rooted tree → prune outgroup | 当前 working scenario |
| S2 | KDOPS + LG+C20 rooted tree → prune outgroup | site-heterogeneous 对照 |
| S3 | ingroup-only midpoint root | outgroup-free 最低成本对照 |
| S4 | ingroup-only MAD root | 推荐的 outgroup-free robustness 对照 |
| S5 | reduced representative set nonreversible/rootstrap | 可选高价值 corroboration |

### 规则 C：哪些任务现在就能做，哪些必须等 QC3
| 任务 | 现在就能做 | 必须等 QC3 |
|---|---:|---:|
| 4.3 core ASR 计算 | ✅ | |
| 4.4 nested ASR | ✅ | |
| 4.6 trait table / PastML 运行 | ✅ | |
| 6.1 DCA 输入准备 | ✅ | |
| Phase 5 节点最终锁定 | | ✅ |
| deepest module origin 的唯一历史叙事 | | ✅ |
| 主文 rooted history 图定稿 | | ✅ |

---

## 4.1 核心树推断与 root scenario 管理（V5.1 强化版）

### 4.1.1 现有成果沿用
- `CoreTree_rooted_MFP.treefile`
- `CoreTree_rooted_LGC20.treefile`
- `CoreTree_rooted_ingroup.treefile`（作为 S1 alias）

### 4.1.2 V5.1 新要求：支持度与 scenario 统一记录
若当前 full-tree 结果尚未同时包含 UFBoot 和 SH-aLRT，则补跑支持度版本：

```bash
# S1: KDOPS + MFP
iqtree -s results/04_phylogeny_asr/core_with_outgroup.afa \
  -m Q.PFAM+F+R10 -alrt 1000 -B 1000 -T AUTO \
  -o "KDOPS_P0A715,KDOPS_Q9ZFK4,KDOPS_O66496" \
  --prefix results/04_phylogeny_asr/CoreTree_rooted_MFP_support

# S2: KDOPS + LG+C20
iqtree -s results/04_phylogeny_asr/core_with_outgroup.afa \
  -m LG+C20+F+G -alrt 1000 -B 1000 -T AUTO \
  -o "KDOPS_P0A715,KDOPS_Q9ZFK4,KDOPS_O66496" \
  --prefix results/04_phylogeny_asr/CoreTree_rooted_LGC20_support
```

> 若 full-tree 补跑成本过高，至少在后续用于 Phase 5 的 focal subtree 上补做 `-alrt 1000 -B 1000`。

### 4.1.3 scenario 产物规范
在 `results/04_phylogeny_asr/` 下维护：

```text
CoreTree_rooted_MFP.treefile
CoreTree_rooted_LGC20.treefile
CoreTree_rooted_MFP_ingroup.treefile
CoreTree_rooted_LGC20_ingroup.treefile
CoreTree_rooted_midpoint_ingroup.treefile
CoreTree_rooted_MAD_ingroup.treefile
root_scenarios.tsv
```

`root_scenarios.tsv` 字段建议：
- `scenario_id`
- `source_alignment`
- `tree_model`
- `rooting_method`
- `outgroup_set`
- `tip_count`
- `support_status`
- `status` (`working` / `accepted_for_sensitivity` / `rejected`)
- `note`

### 4.1.4 outgroup-free rooting（V5.1 新增）
对 ingroup-only tree 至少构建：
- S3：midpoint rooted tree
- S4：MAD rooted tree

> MAD 的价值在于：它不依赖外群，且对分支长度信息更敏感，适合作为当前 KDOPS 根不稳时的对照根。  
> 如果实现成本高，可先完成 midpoint，再补 MAD；但在论文定稿前，MAD 应尽量完成。

### 4.1.5 可选高价值增强：reduced-set nonreversible/rootstrap
在 reduced representative set（例如 46 skeleton 或 183/258 representative set）上，执行 nonreversible rooting / rootstrap 作为额外 corroboration：

- 目的：不是替代全树，而是检查 deep split 方向是否与 S1/S2/S3/S4 相冲突。
- 若该 reduced-set 结果与 S1/S2/S3/S4 中至少一个 outgroup-free scenario 一致，则增强对 root-robust clade 解释的信心。

### 4.1.6 禁止事项
- 不再以“KDOPS 是外群，所以它给出的根一定对”作为默认前提。
- 不因 KDOPS 在 ML 树中的异常行为就回滚 Phase 1–3。

---

## 4.2 结构系统发育交叉验证：AA 树 vs 3Di 树 [高优先级，纳入 QC3]

> V5.1 将 4.2 从“有益补充”提升为“正交证据”。

```bash
# AA tree
iqtree -s results/03_msa_core/skeleton_refined_aa.fa \
  -m MFP -alrt 1000 -B 1000 -T AUTO \
  --prefix results/04_phylogeny_asr/SkeletonTree_AA

# 3Di tree
iqtree -s results/03_msa_core/skeleton_core_3di.fa \
  -m MFP -mset Q.3Di.AF,Q.3Di.LLM,GTR20 \
  -mdef meta/models/Q_3Di_models.nex \
  -alrt 1000 -B 1000 -T AUTO \
  --prefix results/04_phylogeny_asr/SkeletonTree_3Di
```

### 4.2 判读规则
- 4.2 的任务不是证明每一条枝都一致。
- 只问三个问题：
  1. 三大亚型的深分化方向是否一致？
  2. KDOPS / outgroup 分离位置是否与 AA 树大体兼容？
  3. 与模块 gain/loss 解释相关的 focal deep split 是否发生冲突？

### 4.2 输出
- `results/04_phylogeny_asr/tree_comparison.md`
- `results/04_phylogeny_asr/tree_comparison.tsv`

---

## 4.3 核心氨基酸 ASR（可继续，解释分层）

> **4.3 已解决 tip-set mismatch；V5.1 的变化在于解释边界，而不是停止执行。**

### 4.3.1 当前 working tree
- 默认 working tree：`CoreTree_rooted_ingroup.treefile`
- 建议显式保留 alias：
  - `CoreTree_rooted_ingroup.treefile` → 当前默认 scenario
  - `CoreTree_rooted_MFP_ingroup.treefile` → S1 正式文件名

### 4.3.2 运行策略
**主运行：** 在当前 working scenario 上先完成 core ASR。  
**敏感性运行：** 至少在一个 alternative root scenario（优先 S2 或 S4）上补跑 focal-node ASR 或整树 ASR。

```bash
# primary
iqtree -s results/03_msa_core/core_asr.afa \
  -te results/04_phylogeny_asr/CoreTree_rooted_MFP_ingroup.treefile \
  -m MFP -asr \
  --prefix results/04_phylogeny_asr/ASR_core_S1

# sensitivity
iqtree -s results/03_msa_core/core_asr.afa \
  -te results/04_phylogeny_asr/CoreTree_rooted_LGC20_ingroup.treefile \
  -m MFP -asr \
  --prefix results/04_phylogeny_asr/ASR_core_S2
```

### 4.3.3 解释规则
- 在 QC3 之前，`ASR_core_S1` 可用于：
  - 下游结构准备
  - DCA 映射准备
  - focal residue 的初步检查
- 在 QC3 之后，只有下列节点可进入 Phase 5 主线：
  - 节点解释在 ≥2 root scenarios 下不变
  - 支持度满足 Phase 5.1 门控
  - 对应的 gain/loss 事件在 strict/relaxed 注释下稳定

### 4.3.4 产物
- `ASR_core_S1.*`
- `ASR_core_S2.*`
- `asr_node_summary.tsv`

---

## 4.4 嵌套 ASR（亚型内 full-length local ASR）[立即并行推进]

> 这是 V5.1 推荐的第一优先并行任务之一。  
> 原因：它依赖 3.9，而 3.9 已经完成；同时它对最深层 root 的依赖比 family-wide story 小得多。

### 4.4.1 主 exemplar
- **Iβ-ACT architecture** 作为主 exemplar
- 使用 `msa_full_Ib_v4.afa`

```bash
iqtree -s results/03_msa_full/msa_full_Ib_v4.afa \
  -g results/04_phylogeny_asr/Subtree_Ib_topology.nwk \
  -m MFP -asr \
  --prefix results/04_phylogeny_asr/ASR_Ib_local
```

### 4.4.2 V5.1 的新要求
- 必须对 focal Pre-gain / Post-gain 边界节点导出：
  - ML 序列
  - AltAll 序列
  - 每个位点 PP 概要
- `msa_full_Ib_column_map.tsv` 必须参与后续结构注释和论文图注，不得只当中间文件。

### 4.4.3 解释定位
- 4.4 是本项目最可能形成“高可信 exemplar 机制故事”的部分。
- 若 deepest family root 仍不稳，Iβ-ACT nested ASR 仍然可以独立支撑一段强结果。

---

## 4.5 Gap 祖先态重建（非阻塞，继续保留）

### 方案分层
- **方案 A（主线优先）**：模块 presence/absence 走 PastML（见 4.6）
- **方案 B（次优先）**：核心连续 gap block 聚合后走 PastML
- **方案 C（探索性）**：局部 linker/insert 特征只在 focal clade 中解释

### V5.1 修改
- 所有 gap ASR 默认使用 **pruned ingroup rooted scenario trees**，不直接用含 KDOPS 的原始 rooted tree。
- 4.5 不阻塞 4.4 / 6.1，也不阻塞 QC3。

---

## 4.6 模块获得/丢失的离散性状 ASR（PastML）[V5.1 关键加强]

> 这是最直接决定“模块何时出现/消失”的步骤。  
> V5.1 规定：**必须在 strict / relaxed × multiple root scenarios 的组合上运行。**

### 4.6.1 运行矩阵
至少运行以下四组：
- strict × S1
- strict × S2（或 S4）
- relaxed × S1
- relaxed × S2（或 S4）

```bash
for version in strict relaxed; do
  for scenario in S1 S2; do
    pastml \
      --tree results/04_phylogeny_asr/CoreTree_rooted_${scenario}_ingroup.treefile \
      --data results/03_msa_modules/module_presence_absence_${version}.tsv \
      --out_data results/04_phylogeny_asr/module_origins_${version}_${scenario}/ \
      --html_compressed results/04_phylogeny_asr/module_origins_${version}_${scenario}.html
  done
done
```

> 实际文件名可按你的脚本习惯调整；关键是“版本 × 场景”必须显式保留，不可覆盖写。

### 4.6.2 V5.1 新增稳定性分级
对于每一种模块事件，输出 `module_origin_stability.tsv`，并分成：

| stability_class | 定义 |
|---|---|
| robust | strict/relaxed 一致，且在 ≥2 root scenarios 下指向同一 biological event |
| semi-robust | strict/relaxed 一致，但 root scenarios 之间 node 略漂移，biological interpretation 仍同向 |
| root-sensitive | 严重依赖根位置，event polarity 改变 |
| annotation-sensitive | strict / relaxed 差异过大 |
| unresolved | 支持或信息不足 |

### 4.6.3 只有哪些事件能进入主线
仅 `robust` 和部分 `semi-robust` 事件可进入主文。  
`root-sensitive` 事件只能进入补充或 Discussion，并明确加条件语气。

---

## 4.7 AltAll 系综采样（保留，但对准 shortlist 节点）

### V5.1 调整
- 不再对大量节点“平均用力”
- 只对下列节点做 AltAll：
  1. Phase 5 的候选节点
  2. 模块 gain/loss 边界节点
  3. 4.4 exemplar 的关键祖先节点

### 默认阈值
- `PP1 = 0.80`
- `PP2 = 0.20`

### 输出
- `altall_sequences.tsv`
- `focal_nodes_altall_summary.md`

---

## 4.8 模块层局部 ASR（探索性保留）

> 保留，但不再作为主线 gating。

### 定位
- 用于解释特定模块内部保守位点与已知功能位点
- 不用于单独支持 deep origin claim
- 若模块深度过低，结果只能作为定性补充

---

## QC3：系统发育 / ASR 统一 gate（V5.1 强制执行）

> **这是 V5.1 的核心。**  
> QC3 的任务不是问“哪棵树看起来最好”，而是问：  
> **“哪些结论在合理的 root 变动下仍然成立？”**

### QC3 必须包含的内容
1. root scenario 列表与对应模型 / 方法  
2. UFBoot（以及可获得时的 SH-aLRT）  
3. KDOPS subset sensitivity  
4. midpoint / MAD 对照  
5. AA tree vs 3Di tree 对照  
6. strict vs relaxed 模块事件稳定性  
7. focal node 的解释是否在多个场景下保持同向  
8. 对 Phase 5 的推荐：Go / Conditional Go / Hold

### QC3 输出
- `results/04_phylogeny_asr/qc_phylogeny_asr.md`
- `results/04_phylogeny_asr/qc_phylogeny_asr.tsv`
- `results/04_phylogeny_asr/root_scenarios.tsv`
- `results/04_phylogeny_asr/module_origin_stability.tsv`

### QC3 判定标准
#### QC3-GREEN（可进入主线）
- 至少一个 focal event 在 ≥2 root scenarios 下保持同一 biological interpretation
- 主线节点支持度满足 Phase 5.1 门控
- AA / 3Di 不发生致命冲突
- strict / relaxed 不发生 polarity 逆转

#### QC3-YELLOW（可推进，但需降调）
- deepest root 仍不稳
- 但 focal exemplar（如 Iβ-ACT）在 subtype-local 或局部根方案下稳定
- 允许继续 4.4、6.1、5.0；Phase 5 只选 root-robust 节点

#### QC3-RED（暂停历史锁定）
- 关键模块事件在不同根方案下方向翻转
- 3Di 与 AA 对 deep split 出现严重矛盾
- 则冻结：
  - Phase 5 主线节点锁定
  - 主文 rooted history 定稿
- 但仍可继续：
  - 4.4 nested ASR
  - 6.1–6.3 core DCA
  - 5.0 assembly adjudication

---

# 第五步：关键祖先节点的结构验证 [待执行，V5.1 更严格门控]

> V5.1 继续坚持：**2–4 节点、Apo-only、native-oligomer 优先。**  
> 但节点进入 Phase 5 的条件被显著收紧。

---

## 5.0 装配体判定（Assembly Adjudication）[保留为硬门控]

### 步骤
1. 文献 / PDB / UniProt / PDBe PISA 扫描
2. dimer / tetramer 平行 AF3（Apo）
3. PISA 界面评分
4. 写入 `assembly_adjudication.tsv`

### 输出
- `results/05_struct_valid/assembly_adjudication.tsv`

### 解释规则
- 不默认四聚体
- 若装配体状态不清晰，可做短 MD 排除，但不强写机制结论

---

## 5.1 候选祖先节点选择（V5.1 关键修改）

> **废止“bootstrap ≥ 70 即可进入主线”的旧标准。**

### 5.1.1 主线节点必须同时满足
1. **事件稳健性**  
   对应模块 gain/loss 事件属于 `robust` 或强 `semi-robust`
2. **支持度**  
   - 主阈值：UFBoot ≥ 95  
   - 优先 corroboration：SH-aLRT ≥ 80（若 full-tree 未补跑，则在 focal subtree 上补做）
3. **root stability**  
   同一 biological interpretation 在 ≥2 root scenarios 下保持一致
4. **annotation stability**  
   strict / relaxed 不改变事件方向
5. **ASR 可解释性**  
   focal region 的 PP 分布可接受；若 linker 极低置信，仅做降权，不可硬写
6. **assembly 不模糊**  
   装配体判定不能完全 unresolved

### 5.1.2 节点分级
| class | 进入范围 |
|---|---|
| Tier-1 | 主文主图，可进入 AF3 + validation MD 主线 |
| Tier-2 | 补充或 reserve node，可做结构预测，MD 视资源决定 |
| Reject | 不进入 Phase 5 主线，仅保留在 trait ASR 表中 |

### 5.1.3 数量控制
- 默认 **2 个 Tier-1 主线节点**
- 最多再加 **1–2 个 Tier-2 reserve 节点**
- 若 QC3 只给出 1 个 truly robust event，则宁可只做 1 个 exemplar，也不要为凑数选弱节点

### 5.1.4 强制记录
- `results/04_phylogeny_asr/node_selection_registry.tsv`

字段建议：
- `node_id`
- `event_type`
- `scenario_support`
- `strict_relaxed_stability`
- `ufboot`
- `shalrt`
- `assembly_status`
- `tier`
- `decision_note`

---

## 5.2 Apo 结构预测（沿用 V5.0，但绑定 tier）

### Tier-1 必做
- AF3：按 `assembly_adjudication.tsv` 指定拷贝数
- ESMFold：单体或单链交叉检查
- 导出：
  - pLDDT / ipTM
  - PAE
  - 接口面积
  - 催化口袋几何

### Tier-2 可选
- 仅在资源允许时执行 AF3
- 不默认进入 MD

---

## 5.2b 结构 QC 门控

### 三级门控
1. AF3 / 结构预测质量  
   - ipTM ≥ 0.6（多聚体）
   - 跨链 PAE 合理
2. 能量最小化 + 分级平衡  
3. 10 ns short sanity MD（仅作筛选）

### 失败处理
- 若模型不通过：
  - 换 ML-best / AltAll 代表序列重跑
  - 或将该节点降级为 Tier-2 / Reject
- 不允许为了“保留故事完整性”硬保留低可信模型

---

## 5.3 验证性 MD（Apo-only）

### 主分析矩阵
- Tier-1 节点：
  - `native-oligomer-apo`：≥200 ns × 2 replicas

### 排除性测试（可选）
- `alternative-oligomer-apo`：10–20 ns × 1
- 只用于证明不合理装配体，不用于主机制叙事

### V5.1 强调
- MD 是 **结构可行性验证**
- 不是跨证据闭环证明
- 不再扩张为全因子设计

---

## 5.4 有限动力学读出

### 主线读出
- RMSD / RMSF
- 界面面积与稳定性
- 关键口袋体积 / Fpocket
- 活性位点通道几何

### 可写进主文的前提
- 只比较 **root-robust、tier-1** 节点
- 结论写成“与模块整合相容的结构变化”  
  而不是“已完整证明 allosteric communication path”

### QC4 产物
- `results/05_struct_valid/qc_struct_validation.md`

---

# 第六步：核心层共进化分析（DCA） [待执行，V5.1 继续聚焦]

> V5.1 明确：**不要把 module/joint DCA 拉回主线。**

## 6.1 核心层 DCA 输入准备 [可立即并行]
```bash
python scripts/prepare_dca_input.py \
  --input results/03_msa_core/core_asr.afa \
  --gap_col_max 0.50 \
  --gap_row_max 0.50 \
  --meff_min 3.0 \
  --output results/06_dca/core_dca.afa \
  --stats results/06_dca/core_dca_stats.tsv
```

### V5.1 说明
- 6.1 可以在 QC3 完成前就启动
- 因为 DCA 对 deepest root 的依赖很弱
- 但若要把耦联变化与具体 gain event 绑定，仍需等 QC3 后再写强结论

## 6.2 核心层 DCA 执行
```bash
plmc -o results/06_dca/core.params \
  -c results/06_dca/core_couplings.txt \
  -le 16.0 -lh 0.01 -m 100 -t 0.8 \
  -f <focus_seq_id> -g \
  results/06_dca/core_dca.afa
```

## 6.3 显著性评估与接触验证
- Top-L / Top-2L 接触精度
- 映射到代表性结构（1KFL / 1RZM / 3NV8 等）
- 催化位点、界面位点、铰链位点富集
- root-independent 的保守耦联优先解释

### QC5
- `results/06_dca/qc_core_dca.md`

## 6.4 可选探索（不入主线）
- module DCA
- joint core+module DCA
- gain 前后耦联比较

### V5.1 解释边界
- 这些结果最多作为补充或 Discussion 的定性线索
- 不得反向主导主文叙事

---

# 第七步：论文写作蓝图 [待执行，V5.1 重写主叙事]

> **V5.1 的写作核心：从“唯一深根故事”改为“root-robust 主叙事”。**

## 7.1 主文主轴（建议标题方向）
**Conserved catalytic core, recurrent module recruitment, and diversification of allosteric control in DAH7PS**

### 主叙事
1. 结构感知核心 MSA 解决了跨亚型 full-length MSA 失真问题
2. DAH7PS 家族共享保守 core，但反复招募不同 allosteric modules
3. 模块 gain/loss 的一部分事件对 root 方案稳健，可形成可信的历史解释
4. Iβ-ACT 的 nested full-length ASR 提供具体 exemplar
5. core DCA 给出长期保守耦联层
6. 关键祖先结构验证说明这些推断在结构上可行

## 7.2 claim tiers（V5.1 新增）
| tier | 写作位置 | 允许的结论强度 |
|---|---|---|
| root_robust | 主文 Results | 直接陈述 |
| root_sensitive | Supplement / Discussion | 条件性表述 |
| exploratory | Supplement / Outlook | 明确降调 |

### 写作规则
- 如果一个事件在不同根方案下 polarity 翻转，则它不属于 `root_robust`
- 如果一个事件 node 漂移但 biological interpretation 不变，可列为 `semi-robust`，主文需加限定词
- 如果事件只在单一 scenario 下成立，只能放 Supplement 或 Discussion

## 7.3 核心图表清单（V5.1 建议版）
| 图号 | 内容 | 备注 |
|---|---|---|
| Fig 1 | DAH7PS 家族架构与模块分布 | 基于 3.8 |
| Fig 2 | FoldMason 骨架与核心列 QC | 基于 3.1–3.7 |
| Fig 3 | root scenario 比较 + stable / unstable module events | V5.1 新增重点图 |
| Fig 4 | Iβ-ACT nested full-length ASR exemplar | 4.4 |
| Fig 5 | core DCA 网络与代表结构映射 | 6.1–6.3 |
| Fig 6 | Tier-1 ancestral structures / validation MD summary | 5.0–5.4 |
| Fig S1–Sn | gap ASR、AltAll、补充根方案、exploratory DCA | Supplement |

## 7.4 ICDC 在 Discussion 中的定位
- 保留为 **Cross-Evidence Convergence Outlook**
- 只写“定性一致性”
- 不再宣称已经建立 DCA × MD 的正式闭环

---

# 工程风险 Checklist（V5.1）

## [CHECK-01] Core 域提取边缘破碎
- 状态：已通过
- 不改

## [CHECK-02] 祖先全长“弗兰肯斯坦拼接”
- 缓解：全长祖先必须来自 4.4 nested ASR
- 不允许手工拼 core + module

## [CHECK-03] DCA 的 Meff/L 门控
- 主线仅 core DCA
- module/joint DCA 降为探索

## [CHECK-04] Deep root 过度解释（V5.1 新增重点）
- 风险：把单一 outgroup rooted tree 当成唯一历史
- 缓解：QC3 gate + claim tiers
- 失败处理：冻结 Phase 5 节点锁定与 Fig 3 定稿

## [CHECK-05] 结构预测位阻爆炸 / 装配体错误
- 关联：5.0 / 5.2 / 5.2b
- 缓解：assembly adjudication + QC 门控

## [CHECK-06] 文档与数值漂移（V5.1 新增）
- 风险：PLAN / TASKS / README / 论文草稿数字不一致
- 缓解：`metrics_manifest.tsv`
- 失败处理：先修 manifest，再同步全部文档

## [CHECK-07] Phase 4 脚本债拖慢收尾
- 高优先脚本见下节
- 任何新结果必须先能复现，再写进主文

---

# 脚本开发优先级（V5.1）

## 已有或已被使用
- `extract_core_domains.py`
- `define_core_columns.py`
- `annotate_modules.py`
- `extract_module_seqs.py`
- `minimal_trim.py`
- `extract_struct_subset.py`
- `extract_linkers.py`
- `stitch_full_length_msa.py`
- `select_sequences.py`
- `prune_tree.py`
- `assert_tip_match.py`

## 高优先补齐
1. `qc_root_stability.py`
2. `compare_trees.py`
3. `adjudicate_assembly.py`
4. `prepare_dca_input.py`
5. `compute_meff.py`
6. `dca_significance.py`
7. `coordinate_mapper.py`
8. `build_metrics_manifest.py`（V5.1 新增，若无则先手工维护 manifest）
9. `summarize_module_stability.py`（推荐新增）

---

# 里程碑与验收标准（V5.1 修订版）

| 里程碑 | 交付物 | 验收标准 | 状态 |
|---|---|---|---|
| M1 | 结构骨架 + 核心列定义 | core_len 落在可解释范围，结构 QC 通过 | ✅ |
| M2 | 全量 core MSA | `9,393 × 521` match-only core MSA | ✅ |
| M2b | 双版本 trim + 模块注释 | `core_tree.afa` / `core_asr.afa` / module matrices | ✅ |
| M2c | full-length stitched MSA | `msa_full_Ib_v4.afa` + QC2b 断言 | ✅ |
| M3a | rooted working trees | S1 working tree + tip identity PASS | ✅ |
| M3b | QC3 root robustness gate | 至少 2 个 scenario 完成并给出 Go / Conditional Go / Hold | ⬜ |
| M4 | trait ASR stability matrix | `module_origin_stability.tsv` | ⬜ |
| M5 | core DCA | Meff/L 合格 + contact validation | ⬜ |
| M6 | Tier-1 ancestral structure validation | 2 个主线节点结构与 MD QC 通过 | ⬜ |
| M7 | 论文初稿 | claim tiers 落地，Fig 1–6 成形 | ⬜ |

---

# 立即执行顺序（V5.1 推荐）

## A. 不停工并行推进
1. 完成 4.2 AA vs 3Di
2. 完成 4.4 Iβ-ACT nested ASR
3. 跑 4.6 strict / relaxed × S1/S2（或 S4）
4. 启动 6.1 core DCA 输入准备
5. 开始 5.0 assembly adjudication 文献/数据库扫描

## B. 同步完成 QC3 gate
1. 整理 S1 / S2 已有 rooted trees
2. 补建 S3 midpoint
3. 尽量补建 S4 MAD
4. 汇总 `module_origin_stability.tsv`
5. 输出 `qc_phylogeny_asr.md`
6. 给出 Go / Conditional Go / Hold

## C. QC3 之后再做
1. 锁定 Tier-1 / Tier-2 节点
2. 进入 5.2–5.4
3. 定稿 Fig 3 rooted history 叙事
4. 将 root-sensitive 结果移入 Supplement / Discussion

---

# 参考文献（精简版）

1. Cross PJ, Dobson RCJ, Patchett ML, Parker EJ. The diversity of allosteric controls at the gateway to aromatic amino acid biosynthesis. *Biochem J*. 2013;451(2):127–136.  
2. Lang EJM, Cross PJ, Sherwood OC, Parker EJ. Quaternary structure is an essential component that contributes to the sophisticated allosteric regulation mechanism in a key enzyme from *Mycobacterium tuberculosis*. *PLoS ONE*. 2016;12(1):e0169601.  
3. Jiao W, Lang EJM, Bai Y, Fan Y, Parker EJ. Diverse allosteric componentry and mechanisms control entry into aromatic amino acid biosynthesis. *Curr Opin Struct Biol*. 2020;65:159–167.  
4. Gilchrist CLM, Mirdita M, Steinegger M. Multiple protein structure alignment at scale with FoldMason. *Science*. 2026.  
5. Steenwyk JL, Buida TJ III, Li Y, Shen XX, Rokas A. ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. *PLoS Biol*. 2020;18(12):e3001007.  
6. Minh BQ, Schmidt HA, Chernomor O, et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol*. 2020;37(5):1530–1534.  
7. Tria FDK, Landan G, Dagan T. Phylogenetic rooting using minimal ancestor deviation. *Nat Ecol Evol*. 2017;1:0193.  
8. Naser-Khdour S, Minh BQ, Lanfear R. Assessing confidence in root placement on phylogenies: an empirical study using nonreversible models. *Syst Biol*. 2022;71(4):959–972.  
9. Eick GN, Bridgham JT, Anderson DP, Harms MJ, Thornton JW. Robustness of reconstructed ancestral protein functions to statistical uncertainty. *Mol Biol Evol*. 2017;34(2):247–261.  
10. Muñiz-Trejo R, Park Y, Thornton JW. Robustness of ancestral sequence reconstruction to among-site and among-lineage evolutionary heterogeneity. *Mol Biol Evol*. 2025;42(4):msaf084.  
11. Ishikawa SA, Zhukova A, Iwasaki W, Gascuel O. A fast likelihood method to reconstruct and visualize ancestral scenarios. *Mol Biol Evol*. 2019;36(9):2069–2085.  
12. Seffernick JT, Pardo-Avila F, Shen H, et al. A general substitution matrix for structural phylogenetics. *Mol Biol Evol*. 2025;42(6):msaf124.  
13. Ekeberg M, Lövkvist C, Lan Y, Weigt M, Aurell E. Improved contact prediction in proteins: using pseudolikelihoods to infer Potts models. *Phys Rev E*. 2013;87(1):012707.  

---

## 结论性说明

V5.1 的核心不是“再找一棵更好看的树”，而是把项目从 **单一 rooted story** 升级为 **可审稿的 root-robust inference framework**。  
它允许你继续推进大部分工作，同时避免在最脆弱的 deep-root 问题上过度下注。对当前阶段来说，这是最稳、也最容易形成论文的路径。
