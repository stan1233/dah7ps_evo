# DAH7PS 变构起源演化项目计划：V6 审计修订版（SOP rev7）

**A FoldMason-Centred, Core–Module Decoupling Framework with Audit-Driven Root Control**  
**修订日期：2026-04-08**  
**适用范围：替换当前仓库中的 `PLAN.md`，用于审计后路线重排**

---

## 摘要

本项目的科学问题没有改变：我们仍然要重建 DAH7PS 在保守 TIM-barrel 核心之上演化变构调控的历史。  
真正发生变化的是**优先级**：

- 旧瓶颈：全长远缘同源下的 MSA 崩坏与 tip-set mismatch  
- 新瓶颈：**root scenario 的 provenance、可追溯性与敏感性**

审计表明，Phase 0–3 已经形成足够稳固的基础；因此本次修订**不是返工**，而是把路线改成：

> **先修 root scenario 与 cross-scenario ASR，再推进模块历史和结构验证。**

本计划的核心原则是：

1. **不推翻现有主线**：sequence mining、FoldMason 核心、core MSA、模块注释、core DCA、少量结构验证仍然成立。  
2. **冻结最深层叙事**：在 QC3 修复前，允许继续计算，不允许把 deepest history 写成唯一故事。  
3. **以 provenance 为先**：有结果不等于有证据链；没有 treefile / summary / log / md5 一致性，就不能进入 narrative gate。  
4. **以 scenario 比较替代单根叙事**：S1 / S2 / S3 / S4 修复后共同决定可写结论，必要时上 S5。  
5. **Phase 5 节点要晚于 QC3**：AF3 / MD 只服务于 root-robust 节点，不再反向“证明”深根。

---

## 一、经过审计后确认的项目状态

### 1. 稳定层（不返工）

以下层级保持不变，继续作为主线输入：

- Phase 0：环境、版本、参数、目录规范
- Phase 1：全库挖掘、KDOPS 负筛、互斥 subtype 集合
- Phase 2：长度/覆盖度过滤、Type II stitching、NR80 / seeds60 / stepping-stone
- Phase 3：FoldMason 骨架、LDDT core columns、全量核心映射、双版本 trim、模块 strict/relaxed 矩阵、Ib-ACT full-length stitching
- S1 rooted tree 和 S1 core ASR

### 2. 当前阻塞层（必须优先修）

#### 2.1 S4 provenance 冲突
当前所谓 S4（legacy filename 含 `MAD`）存在以下问题：

- 不同脚本生成过不同 root split
- treefile、summary、log 之间记录不一致
- top-500 与 full-search 输出没有被显式区分
- 在这种状态下，S4 不能被当作正式 root evidence 使用

#### 2.2 S2 ingroup prune + ASR 缺失
虽然 S2 rooted tree 已存在，但缺少：

- KDOPS prune
- tree / alignment tip-set consistency 断言
- S2 core ASR state file
- 与 S1 的 residue-level / node-level 对照

这使得当前无法区分：

- model-sensitive 变化
- root-sensitive 变化
- 仅仅是 documentation / provenance artifact

#### 2.3 模块边界敏感性未文档化
strict / relaxed 矩阵已经存在，但对于以下模块：

- `N_ext`
- `C_tail`

两种定义差异很大，说明当前模块历史不能直接解释为“真实 gain/loss 时间线”，必须先补：

- 阈值规范
- 边界敏感性表
- 必要的抽样人工复核

#### 2.4 DCA 尚未进入审计
DCA 没有任何产物，但输入已经齐备。  
因此它不是阻塞项，而是**可并行启动的审计分支**。

---

## 二、修订后的总路线

### 路线总纲

**稳定基础层（已完成）**  
→ **Phase 4.0 审计热修复：root scenario provenance**  
→ **Phase 4.1 S2 ingroup ASR 与 cross-scenario sensitivity**  
→ **Phase 4.2 模块 trait ASR 与 annotation sensitivity**  
→ **Phase 4.3 DCA 审计并行**  
→ **Phase 5 仅对通过 gate 的节点做结构验证**  
→ **Phase 6 论文叙事按 claim tiers 分层**

---

## 三、执行原则与解释边界

### 3.1 working tree ≠ narrative tree
- 现有 tree 可以继续用来跑计算
- 但不能因为“能算”就等于“能写结论”

### 3.2 provenance 缺陷优先于统计结果
- 支持度高但 provenance 混乱的结果，不可用于 gate
- provenance 清晰但结论敏感的结果，可以进入 `root_sensitive`

### 3.3 S4 在修复前降级命名
- 旧文件名保留用于兼容
- 新文档中统一写作 **S4 provisional ingroup reroot scenario**
- 不再称作正式 MAD，直到方法学标签和输出链条统一

### 3.4 模块历史先看稳定性，再谈极性
- 没有 `module_origin_stability.tsv`
- 就没有资格写“某模块先于另一模块获得”

### 3.5 结构验证只服务于 root-robust 节点
- AF3 / MD 不是用来替代 phylogeny 的
- 不能拿结构模型反推哪一个 deepest root 更真

---

## 四、分阶段工作包

# Phase 0–3：稳定基础层（已完成，保留不变）

### 已完成的关键产物
- `results/01_mining/`
- `results/02_qc/`
- `results/03_msa_core/`
- `results/03_msa_modules/`
- `results/03_msa_full/`
- `results/meta/metrics_manifest.tsv`
- `results/meta/progress_snapshot.md`

### 本阶段不再追加返工任务
- 不返工 sequence mining
- 不返工 FoldMason skeleton
- 不返工 core MSA
- 不返工 full-length stitching

---

# Phase 4.0：审计热修复（新增，最高优先级）

## 4.0.1 冻结并整理现有 S4 产物

### 目标
把当前混乱的 S4 输出链条变成**可追溯、可复核、不可覆盖**的状态。

### 动作
1. 冻结现有所有 S4 相关 treefile、summary、log
2. 用显式文件名区分：
   - `..._top500...`
   - `..._fullsearch...`
   - 如仍保留 legacy 文件，必须在 `root_scenarios.tsv` 注明它实际指向哪一个
3. 为每个文件记录：
   - source script
   - 运行参数
   - search space
   - root split
   - rho / objective value
   - md5
   - 生成时间
4. 重建 `root_scenarios.tsv`
5. 补写 `log.md` 中缺失的 S3 / S4 条目

### Done 条件
- 任一 scenario 都不允许“同名不同内容”
- treefile / summary / log / manifest 四者一致
- legacy filename 的实际身份可被唯一追溯

## 4.0.2 重新定义 S4 的方法学身份

### 原则
在正式确认实现前，S4 统一视为：

> **ingroup reroot proxy**

### 允许的两种结局
- **A.** 若验证后为 formal MAD：恢复 MAD 命名  
- **B.** 若本质是 clock-balance / root-to-tip variance proxy：保留 proxy 命名，不再冒充 formal MAD

### Done 条件
- 方法学标签与脚本实现一致
- 文档不再混用“proxy”与“formal MAD”

---

# Phase 4.1：S2 ingroup prune + model-sensitive core ASR（最高优先级）

## 4.1.1 S2 KDOPS prune
### 输入
- `results/04_phylogeny_asr/CoreTree_rooted_LGC20.treefile`

### 输出
- `results/04_phylogeny_asr/CoreTree_rooted_LGC20_ingroup.treefile`

## 4.1.2 tree / alignment tip-set 断言
### 输入
- `CoreTree_rooted_LGC20_ingroup.treefile`
- `results/03_msa_core/core_asr.afa`

### 输出
- `assert_tip_match` 通过记录
- 断言结果写入 `log.md`

## 4.1.3 S2 core ASR
### 输出
- `results/04_phylogeny_asr/ASR_core_S2.state`
- `results/04_phylogeny_asr/ASR_core_S2.iqtree`
- `results/04_phylogeny_asr/ASR_core_S2.treefile`
- `results/04_phylogeny_asr/ASR_core_S2.log`

## 4.1.4 S1 vs S2 对比摘要
### 输出
- `results/04_phylogeny_asr/asr_sensitivity_summary.tsv`
- `results/04_phylogeny_asr/asr_node_summary.tsv`

### 至少包含
- 节点映射方式说明
- 每个 focal node 的 residue disagreement rate
- 高 PP 位点是否稳定
- 哪些节点可进入后续 module / structure 评估

### Done 条件
- 至少有一套 alternative model ASR（S2）
- 能够把“model-sensitive”与“仅 root-sensitive”粗分出来

---

# Phase 4.2：cross-scenario root / trait sensitivity（审计主线）

## 4.2.1 root scenario matrix
### 场景集合
- S1：KDOPS + MFP
- S2：KDOPS + LG+C20
- S3：ingroup midpoint
- S4：repaired provisional reroot
- S5：仅在必要时触发

### 输出
- `results/04_phylogeny_asr/root_scenarios.tsv`
- `results/04_phylogeny_asr/qc_root_stability.md`
- `results/04_phylogeny_asr/root_sensitivity_matrix.tsv`

### 必须回答的问题
1. 哪些结论在 S1/S2/S3/S4 下方向一致？
2. 哪些结论只依赖 S4 或只依赖 outgroup rooting？
3. 哪些冲突来自 provenance defect，而不是 biological signal？

## 4.2.2 S5 触发条件
满足任一条件时启动：

- S1 / S2 / S3 / repaired S4 的深节点方向持续冲突
- outgroup 与 outgroup-free 场景无法给出最小稳定共识
- Phase 5 节点选择依然高度依赖单一根位

### S5 目标
- reduced representative set
- taxon-balanced
- nonreversible rooting / rootstrap
- 给出“可写 / 不可写”的最终边界，而不是执着于唯一根

---

# Phase 4.3：模块 trait ASR 与 annotation sensitivity（高优先级）

## 4.3.1 补模块定义文档
### 输出
- `results/03_msa_modules/module_definition_spec.md`
- `results/03_msa_modules/module_boundary_sensitivity.tsv`

### 至少包括
- strict / relaxed 的长度阈值、坐标阈值、HMM 规则
- 每个模块在 strict vs relaxed 下的数量差异
- 对 `N_ext` 与 `C_tail` 的重点解释
- 需要时抽样人工复核实例

## 4.3.2 PastML trait ASR
### 最低运行组合
- strict × S1
- strict × S2
- relaxed × S1
- relaxed × S2

### S4 修复后可追加
- strict × S4
- relaxed × S4

### 输出
- `results/04_phylogeny_asr/module_origin_stability.tsv`
- `results/04_phylogeny_asr/module_trait_asr_summary.md`

### 稳定性标签
- `root_robust`
- `root_sensitive`
- `annotation_sensitive`
- `unresolved`

### Done 条件
- 不再用启发式语言判断模块历史
- 所有模块 gain/loss 结论都有稳定性标签

---

# Phase 4.4：DCA 审计并行（可立刻启动）

## 4.4.1 输入准备
### 输入
- `results/03_msa_core/core_asr.afa`

### 输出
- `results/06_dca/core_dca.afa`
- `results/06_dca/core_dca_stats.tsv`

### 必须记录
- 序列数
- 列数
- gap 过滤策略
- Meff
- L
- Meff/L

## 4.4.2 门控
- 若 `Meff/L >= 3.0`：允许进入主线 plmc
- 若 `< 3.0`：停止在审计层，不进入主线解释

## 4.4.3 主线输出
- `results/06_dca/core_couplings.txt`
- `results/06_dca/qc_core_dca.md`

### Done 条件
- DCA 有正式 gate，不再停留在 README 口头承诺

---

# Phase 5：结构验证（QC3 后置，当前仅做准备）

## 当前允许的动作
- 文献与 PDB 装配体注释整理
- PISA 先验信息表
- 结构模板盘点
- 目录与脚本准备

## 当前不允许的动作
- 正式提名 Tier-1 / Tier-2 节点
- 以未通过 gate 的 deepest node 进入 AF3 / MD 主线
- 因为“节点看起来有趣”就越过 QC3

## 进入 Phase 5 的门槛
某节点必须同时满足：

1. 支持度门槛  
2. root stability 门槛  
3. strict / relaxed annotation stability  
4. ASR interpretability（关键位点后验概率可解释）

### 输出
- `results/05_struct_md/assembly_adjudication.tsv`
- `results/04_phylogeny_asr/node_selection_registry.tsv`

---

# Phase 6：论文叙事与图表

## 主文只写三类内容
1. 稳定基础层：核心 MSA、FoldMason、家族框架  
2. `root_robust` 事件：跨场景方向一致  
3. DCA / 结构验证中与上述事件一致的支持证据  

## 补充材料与 Discussion 承担
- `root_sensitive` 事件
- `annotation_sensitive` 事件
- S5 仍无法裁决的 deepest root 争议
- 探索性 module / structure 观察

---

## 五、门控体系（Gate System）

### Gate P：provenance gate
**目标**：scenario 身份清晰  
**通过条件**：
- 有且仅有一个官方输出链条
- treefile / summary / log / md5 一致
- S3 / S4 写入 `log.md`

### Gate M：model-sensitivity gate
**目标**：S2 ASR 到位  
**通过条件**：
- S2 prune + assert_tip_match + ASR 完成
- S1 vs S2 差异摘要完成

### Gate T：trait-stability gate
**目标**：模块历史可分层  
**通过条件**：
- strict / relaxed × S1 / S2 完整
- `module_origin_stability.tsv` 生成

### Gate N：node-selection gate
**目标**：Phase 5 候选节点真实可选  
**通过条件**：
- 支持度 + root stability + annotation stability + interpretability 同时通过

### Gate R：narrative gate
**目标**：主文允许写 deepest-history 相关内容  
**通过条件**：
- 结论至少达到 `root_robust`
- provenance 无缺陷

---

## 六、文件与命名规范

### 6.1 scenario 文件不可覆盖
任何 rooting scenario 一旦输出，文件名必须包含其身份信息，例如：

- model / reroot method
- top500 / fullsearch
- ingroup / outgroup
- 日期或 run id（可选）

### 6.2 必须有 summary + log
每个 scenario 至少要有：

- `.treefile`
- `.summary.txt` 或 `.summary.tsv`
- `.log`
- `root_scenarios.tsv` 记录行

### 6.3 manifest 单一真源
- 数字以 `results/meta/metrics_manifest.tsv` 为准
- 状态以 `results/meta/progress_snapshot.md` 为准
- README / TASKS / PLAN 只能回显，不能自创数字

---

## 七、当前优先级排序

### P0（阻塞性）
1. S4 provenance repair
2. S2 prune + ASR
3. `root_scenarios.tsv` 与 `log.md` 修复

### P1（高优先）
4. S1 vs S2 ASR 差异摘要
5. 模块定义文档 + boundary sensitivity
6. strict / relaxed × scenario trait ASR

### P2（并行）
7. DCA 输入准备 + Meff/L 审计
8. Assembly adjudication 文献整理

### P3（后置）
9. Phase 5 节点选择
10. AF3 / MD
11. 主文 deepest-history 叙事

---

## 八、非目标（本版明确不做）

- 不重跑 Phase 1–3
- 不用结构模型替代 phylogeny
- 不在 QC3 修复前给 deepest history 定唯一箭头
- 不让 placeholder node registry 伪装成“已筛选节点集”

---

## 九、成功标准

本轮修订成功，不是指“已经找到唯一根”，而是指：

1. **所有 scenario 都可追溯**
2. **S2 ASR 能与 S1 做正式比较**
3. **模块历史有稳定性标签**
4. **DCA 有实际 gate 和实际输出**
5. **Phase 5 的启动条件被真正执行**
6. **主文只写 root-robust 结论**

只要做到这一点，项目就从“可能因为深根而失控”转为“即使深根不唯一，也能稳健产出可发表结果”。

---

## 附：建议的最小执行顺序

```text
S4 provenance repair
  -> S2 prune + ASR
  -> S1 vs S2 ASR sensitivity summary
  -> module_definition_spec + boundary sensitivity
  -> strict/relaxed × S1/S2 trait ASR
  -> DCA Meff/L audit
  -> (if needed) S5 reduced-set rooting
  -> node selection registry
  -> assembly adjudication
  -> AF3 / MD on gated nodes
```
