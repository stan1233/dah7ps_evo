# TASKS.md — DAH7PS V6 审计修订版任务清单

> 本文件基于 2026-04-08 诊断审计重排优先级。  
> 所有数字以 `results/meta/metrics_manifest.tsv` 为准；所有状态以 `results/meta/progress_snapshot.md` 为准。  
> 状态标记：`[x]` 完成 / `[/]` 进行中 / `[ ]` 待做 / `[-]` 暂停或被 gate 阻塞

---

## 0. 全局规则（新增）

- [x] 将 Phase 0–3 视为稳定基础层，不因当前 root 问题返工
- [x] 将“working tree ≠ narrative tree”写入 README / PLAN / TASKS
- [x] 将 Phase 5 改为 **QC3 后置**
- [ ] 所有 scenario 强制记录 `treefile + summary + log + source script + md5`
- [ ] 将 S3 / S4 缺失的执行记录补入 `log.md`
- [ ] 统一 Phase 5 目录名为 `results/05_struct_md/`
- [ ] 清理旧文档中的过期状态与路径漂移

---

## 1. 已稳定完成的基础层（不返工）

### Phase 0–3 stable baseline

- [x] 环境、软件版本、参数文件、目录规范
- [x] UniRef90 序列挖掘与 subtype HMM 流程
- [x] KDOPS 负筛与 Ib 清洗
- [x] Ia / II overlap 竞争归属
- [x] Phase 2 长度 / 覆盖度 QC
- [x] Type II stitching rescue
- [x] NR80 / seeds60 / stepping-stone
- [x] FoldMason 结构面板与 skeleton
- [x] LDDT-based core column definition
- [x] `core_global_matchonly.afa`、`core_tree.afa`、`core_asr.afa`
- [x] strict / relaxed 模块矩阵
- [x] Ib-ACT full-length profile-anchored stitching
- [x] S1 rooted tree
- [x] S1 core ASR

> 说明：以上层级默认继续沿用，不作为当前阻塞项。

---

## 2. P0 Critical Path（当前真正阻塞主线的任务）

### 2.1 S4 provenance 修复（最高优先级）

- [/] 盘点现有 S4 相关文件（treefile / summary / log / bak / alternative output）
- [ ] 冻结现有 legacy 文件，不再允许覆盖写入
- [ ] 给 S4 各版本显式命名：
  - [ ] `..._top500...`
  - [ ] `..._fullsearch...`
  - [ ] 若存在 legacy 官方文件，写清它当前实际对应哪一版
- [ ] 为每个 S4 输出记录 md5
- [ ] 在 `root_scenarios.tsv` 中加入：
  - [ ] scenario_id
  - [ ] label
  - [ ] source_script
  - [ ] search_space
  - [ ] objective / rho
  - [ ] root split
  - [ ] output_treefile
  - [ ] md5
  - [ ] notes
- [ ] 补写 `log.md` 的 S3 / S4 条目
- [ ] 重新生成或修复：
  - [ ] `CoreTree_rooted_MAD_ingroup_summary.txt`（若继续保留）
  - [ ] `qc_root_stability.md`
- [ ] 在文档中把 S4 暂时改写为 **provisional ingroup reroot scenario**
- [ ] 判断当前 S4 的方法学身份：
  - [ ] formal MAD
  - [ ] balance / clock-like reroot proxy
  - [ ] 未定（则继续用 provisional）

**Done 条件：**
- [ ] 同名文件不再对应不同内容
- [ ] treefile / summary / log / manifest 一致
- [ ] S4 的名字和方法学身份一致

---

### 2.2 S2 KDOPS prune + assert_tip_match + ASR（最高优先级）

- [x] `CoreTree_rooted_LGC20.treefile` 已存在
- [ ] KDOPS prune → `CoreTree_rooted_LGC20_ingroup.treefile`
- [ ] `assert_tip_match.py` 通过并写入日志
- [ ] 提交 S2 core ASR
- [ ] 获得以下文件：
  - [ ] `ASR_core_S2.state`
  - [ ] `ASR_core_S2.iqtree`
  - [ ] `ASR_core_S2.treefile`
  - [ ] `ASR_core_S2.log`
- [ ] 更新 `progress_snapshot.md`
- [ ] 将 S2 状态从“tree-only”升级为“ASR-ready”

**审计报告中建议的命令序列：**
```bash
conda run -n dah7ps_v4 python scripts/prune_tree.py \
  --input results/04_phylogeny_asr/CoreTree_rooted_LGC20.treefile \
  --remove KDOPS \
  --output results/04_phylogeny_asr/CoreTree_rooted_LGC20_ingroup.treefile

conda run -n dah7ps_v4 python scripts/assert_tip_match.py \
  --tree results/04_phylogeny_asr/CoreTree_rooted_LGC20_ingroup.treefile \
  --alignment results/03_msa_core/core_asr.afa

iqtree -s results/03_msa_core/core_asr.afa \
  -te results/04_phylogeny_asr/CoreTree_rooted_LGC20_ingroup.treefile \
  -m LG+C20+F+G --asr -T 20 \
  --prefix results/04_phylogeny_asr/ASR_core_S2
```

**Done 条件：**
- [ ] S2 至少能与 S1 做正式 ASR 对比
- [ ] 不再只有 S1 一个场景能谈 residue-level 祖先态

---

### 2.3 Cross-scenario ASR sensitivity 摘要

- [ ] 建立 S1 vs S2 节点对应策略
- [ ] 生成 `asr_node_summary.tsv`
- [ ] 生成 `asr_sensitivity_summary.tsv`
- [ ] 统计每个 focal node 的：
  - [ ] residue disagreement rate
  - [ ] 高 PP 位点是否稳定
  - [ ] 结构解释关键位点是否一致
- [ ] 将节点分为：
  - [ ] 可进入 Phase 5 候选池
  - [ ] 仅 root-sensitive / model-sensitive
  - [ ] 暂不解释

**Done 条件：**
- [ ] 可以回答“某节点不稳定，到底是根问题还是模型问题”

---

## 3. P1 高优先任务（紧随 P0 之后）

### 3.1 模块定义文档化与边界敏感性

- [ ] 产出 `results/03_msa_modules/module_definition_spec.md`
- [ ] 产出 `results/03_msa_modules/module_boundary_sensitivity.tsv`
- [ ] 为每个模块记录：
  - [ ] strict 规则
  - [ ] relaxed 规则
  - [ ] HMM 判定条件
  - [ ] 坐标判定条件
  - [ ] 数量差异
- [ ] 重点解释：
  - [ ] `N_ext`
  - [ ] `C_tail`
- [ ] 对边界变化最大的模块抽样人工复核
- [ ] 将结论写入 `module_trait_asr_summary.md` 的方法部分

**Done 条件：**
- [ ] strict / relaxed 不再只是两个黑箱矩阵
- [ ] 后续 trait ASR 能够解释 `annotation_sensitive`

---

### 3.2 模块 trait ASR（PastML）

- [ ] 运行 `strict × S1`
- [ ] 运行 `strict × S2`
- [ ] 运行 `relaxed × S1`
- [ ] 运行 `relaxed × S2`
- [ ] 若 S4 修复完成：
  - [ ] `strict × S4`
  - [ ] `relaxed × S4`
- [ ] 生成 `module_origin_stability.tsv`
- [ ] 将事件标注为：
  - [ ] `root_robust`
  - [ ] `root_sensitive`
  - [ ] `annotation_sensitive`
  - [ ] `unresolved`

**Done 条件：**
- [ ] 模块 gain/loss 不再靠启发式描述
- [ ] 最深层模块极性若不稳，必须显式标出来

---

### 3.3 root scenario matrix 与 QC3 修复

- [ ] 修复 `root_scenarios.tsv`
- [ ] 更新 `qc_root_stability.md`
- [ ] 生成 `root_sensitivity_matrix.tsv`
- [ ] 在 QC3 中单独列出：
  - [ ] provenance pass / fail
  - [ ] S1 vs S2
  - [ ] S3 vs repaired S4
  - [ ] strict vs relaxed 模块稳定性
  - [ ] 是否需要 S5
- [ ] 给出新的 QC3 结论：
  - [ ] GO
  - [ ] Conditional GO
  - [ ] HOLD

**Done 条件：**
- [ ] QC3 不再依赖 stale summary
- [ ] QC3 结论能够真正支配 Phase 5 和 manuscript narrative

---

## 4. P2 并行分支（不阻塞 P0/P1，但可以立即开始）

### 4.1 DCA 审计

- [ ] 准备 `results/06_dca/core_dca.afa`
- [ ] 计算 `Meff`
- [ ] 计算 `L`
- [ ] 计算 `Meff/L`
- [ ] 生成 `results/06_dca/core_dca_stats.tsv`
- [ ] 判断是否满足主线门槛 `Meff/L >= 3.0`
- [ ] 若通过：
  - [ ] 运行 plmc
  - [ ] 生成 `core_couplings.txt`
  - [ ] 生成 `qc_core_dca.md`
- [ ] 若不通过：
  - [ ] 在 `qc_core_dca.md` 中终止并说明原因

**Done 条件：**
- [ ] DCA 从“README 里的计划”变成“有实际 gate 的分析分支”

---

### 4.2 Assembly adjudication 准备（Phase 5 前置准备）

- [ ] 建立 `results/05_struct_md/`
- [ ] 收集 extant structure / literature 中的装配体证据
- [ ] 建立 `assembly_adjudication.tsv`
- [ ] 记录：
  - [ ] PDB biological assembly
  - [ ] 文献描述的装配状态
  - [ ] PISA 先验信息
  - [ ] 可能的 native oligomer 候选
- [ ] 仅做准备，不提名祖先节点

**Done 条件：**
- [ ] Phase 5 启动前已有装配体证据底稿
- [ ] 不会出现“先跑 AF3，后补装配体依据”的倒序问题

---

## 5. 条件触发分支

### 5.1 S5 reduced-set rooting（仅在必要时启动）

- [ ] 设计 reduced representative set
- [ ] 保证 subtype-balance 与 taxon-balance
- [ ] 运行 nonreversible rooting / rootstrap
- [ ] 将结果写入 `root_sensitivity_matrix.tsv`
- [ ] 更新 QC3 结论

**触发条件：**
- [ ] S1 / S2 / S3 / repaired S4 仍持续冲突
- [ ] 或者 Phase 5 节点选择仍高度依赖单一 root

---

## 6. 当前被 gate 阻塞的任务

### 6.1 Phase 5 节点选择

- [-] 正式提名 Tier-1 节点
- [-] 正式提名 Tier-2 节点
- [-] 围绕 deepest node 启动 AF3 主线
- [-] 围绕未 gated 节点启动 MD 主线

**解除阻塞条件：**
- [ ] provenance gate 通过
- [ ] S2 ASR 完成
- [ ] module stability 表完成
- [ ] node_selection_registry.tsv 有真实评估而非 placeholder

---

### 6.2 主文 deepest-history 叙事

- [-] 锁定单一 deepest root 故事
- [-] 在主文中给出唯一模块起源顺序
- [-] 把 root-sensitive 结论写成硬结论

**解除阻塞条件：**
- [ ] 至少达到 `root_robust`
- [ ] 或者明确降级到 Supplement / Discussion

---

## 7. 文档与可复现性维护

- [ ] README 更新为 2026-04-08 审计状态
- [ ] PLAN 更新为 V6 审计修订版
- [ ] TASKS 更新为审计版 critical path
- [ ] `metrics_manifest.tsv` 中补齐：
  - [ ] DCA 行
  - [ ] S2 ASR 行
  - [ ] S4 provenance 行
- [ ] `progress_snapshot.md` 与 README / PLAN / TASKS 一致
- [ ] `root_scenarios.tsv` 与 QC3 / summary / log 一致

**Done 条件：**
- [ ] 文档不再出现“README 还停在旧状态，但 snapshot 已更新”的漂移

---

## 8. 当前建议的执行顺序

```text
[1] S4 provenance repair
[2] S2 prune + assert_tip_match + ASR
[3] S1 vs S2 ASR sensitivity summary
[4] module_definition_spec + boundary sensitivity
[5] strict/relaxed × S1/S2 trait ASR
[6] QC3 repair and new verdict
[7] DCA Meff/L audit (parallel)
[8] assembly_adjudication preparation (parallel)
[9] if needed -> S5 reduced-set rooting
[10] node_selection_registry
[11] AF3 / MD on gated nodes
[12] manuscript claim-tier lock
```

---

## 9. 本轮任务完成的标志

- [ ] S4 不再是 provenance 黑箱
- [ ] S2 不再只是 tree-only scenario
- [ ] 模块历史不再只有启发式语言
- [ ] DCA 不再是空白
- [ ] Phase 5 的启动时点被真实 gate 控制
- [ ] 主文不再提前写死 deepest-history narrative
