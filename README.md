# :dna: DAH7PS Allosteric Evolution & Ancestral Sequence Reconstruction

> **Audit-adjusted project status (V6, 2026-04-08)**
>
> 当前工作模式：**计算推进 = GO；深根叙事锁定 = HOLD；Phase 5 节点锁定 = HOLD**  
> 当前版本基于 2026-04-08 诊断审计结果重写，用来替换仓库内过期的 README / PLAN / TASKS 状态描述。

---

## 项目目标

本项目旨在回答一个核心问题：

**DAH7PS 如何在保守的 TIM-barrel 催化核心之上，反复募集不同的变构模块，并最终形成稳定、可选择的调控表型？**

项目主线仍然保持不变：

**结构感知核心 MSA → 多场景系统发育 / ASR → 模块 trait ASR → 核心层 DCA → 少量关键祖先的结构验证**

需要调整的不是科学问题本身，而是**路线优先级**和**解释边界**。

---

## 2026-04-08 审计后的核心结论

### 结论一句话版

- **Phase 0–3 是稳的，不需要返工。**
- **真正的瓶颈已经从 MSA / tip-set consistency 转移为 root scenario 的可追溯性与敏感性。**
- **当前不能继续按“深根已基本锁定”的假设推进。**

### 当前项目状态

| 层级 | 状态 | 当前判断 |
|---|---|---|
| Phase 0–3（挖掘、QC、去冗余、核心 MSA、FoldMason、模块矩阵） | ✅ 稳定 | 作为后续全部分析的基础继续沿用 |
| S1：KDOPS + MFP rooted tree | ✅ 可用 | 已完成，可继续作为 working scenario |
| S1：核心 ASR | ✅ 可用 | `ASR_core.state` 已存在，可作为当前唯一完整 ASR 基线 |
| S2：KDOPS + LG+C20 tree | ⚠️ 部分完成 | 树已存在，但 **未 prune KDOPS、未 assert tip match、未跑 ASR** |
| S3：ingroup midpoint | ✅ 可用 | 可作为 outgroup-free 对照场景 |
| S4：ingroup reroot（原称 MAD） | 🔴 阻塞 | **treefile / summary / log 之间存在 provenance 冲突**，在修复前不能作为正式 narrative gate |
| 模块 strict / relaxed 矩阵 | ⚠️ 可用但需审计 | N_ext 和 C_tail 对边界定义高度敏感，必须先补阈值文档与稳定性分析 |
| DCA | ❌ 未开始 | 可以并行启动，但必须先完成 Meff/L 审计 |
| Phase 5 结构验证 | ❌ 未开始 | 当前只允许做 assembly adjudication 准备，不允许锁定节点名单 |

---

## 本次路线调整的核心变化

### 1. 深根叙事立即冻结

从现在开始，所有依赖 deepest root 极性的结论都必须降级处理：

- **主文只接受 `root_robust` 结论**
- **`root_sensitive` 结论只能进入补充材料或 Discussion**
- **存在 provenance 缺陷的场景不能参与 narrative gating**

### 2. S4 不再直接被当作“formal MAD”

在当前实现和结果登记修复之前：

- `CoreTree_rooted_MAD_ingroup.treefile` 只能被视为 **S4 provisional ingroup reroot scenario**
- 不再把它直接描述成正式 MAD 证据
- 后续如果完成统一脚本、统一日志、统一 summary，并确认方法学标签无误，才允许恢复正式命名

### 3. S2 prune + ASR 上升为最高优先级

当前最关键的不是再加更多探索图，而是先回答：

> **祖先位点差异到底是 model-sensitive，还是 root-sensitive？**

这只有在 S2 完成以下步骤后才能判断：

1. KDOPS prune  
2. tree / alignment tip-set 一致性断言  
3. S2 场景的核心 ASR  
4. 与 S1 做节点和位点级比较  

### 4. 模块历史解释必须延后到 strict/relaxed × scenario 稳定性矩阵之后

模块 gain/loss 的时间顺序、早晚关系和极性，必须以以下结果为前提：

- `strict × S1`
- `strict × S2`
- `relaxed × S1`
- `relaxed × S2`
- 若 S4 修复完成，再加入 `× S4`

没有这一步，模块故事只能叫**启发式观察**，不能叫正式结果。

### 5. DCA 可以并行，但只能先做审计分支

DCA 现在可以启动，因为它与 root 的耦合较弱；但必须按以下顺序推进：

1. 准备 `core_dca.afa`
2. 计算 Meff、L、Meff/L
3. 只有在 **Meff/L ≥ 3.0** 时才允许进入 plmc 主分析
4. 若门槛不达标，则直接停在审计报告，不做主线解释

### 6. Phase 5 只允许做“准备工作”，不允许锁定候选祖先

在 QC3 修复完成前，Phase 5 只能做：

- 文献装配体证据整理
- PDB / PISA 先验注释
- 结构模板盘点
- 需要的脚本和目录准备

**不允许做的事：**

- 正式提名 Tier-1 / Tier-2 节点
- 围绕某个 deepest ancestor 展开 AF3 / MD 主线投入
- 把当前 placeholder registry 误当作已筛选名单

### 7. 若 S1 / S2 / S3 / 修复后的 S4 仍不稳定，启动 S5

触发条件：

- S4 provenance 修复完成后，root 结论仍然大幅不一致
- 或者 S1 / S2 的深节点方向持续冲突，无法形成最小稳健共识

则新增：

- **reduced representative set**
- **nonreversible rooting / rootstrap**
- 目标不是强行“找唯一根”，而是判定哪些结论永远不该写成硬结论

---

## 审计后保留不变的主策略

| 组件 | 保留策略 | 当前说明 |
|---|---|---|
| Sequence mining | UniRef90 + subtype HMM + KDOPS negative selection | 有效，不返工 |
| Structural skeleton | FoldMason + AFDB + PDB anchor | 有效，不返工 |
| Core MSA | hmmalign + Stockholm + strip inserts | 有效，不返工 |
| Module annotation | 坐标规则 + HMM 域命中 | 保留，但要补边界敏感性文档 |
| Tree inference | 多 scenario 并行 | 保留，但 S4 必须先修 provenance |
| Core ASR | ingroup pruned trees | 保留，当前重点转向补 S2 |
| Module trait ASR | strict / relaxed × scenarios | 保留，优先级提高 |
| DCA | core-only mainline | 保留，但先做 Meff/L gate |
| Structural validation | 少量节点、Apo-only、native-oligomer | 保留，但等待 QC3 |
| Manuscript logic | root-robust main text | 保留，并且现在必须严格执行 |

---

## 当前审计通过的数据快照

> 所有数值应继续以 `results/meta/metrics_manifest.tsv` 为单一真源；README 只做摘要展示。

| 指标 | 当前值 |
|---|---|
| PASS sequences | 24,202 |
| NR80 | 9,673 |
| Core mapped sequences | 9,393 |
| Core match-only MSA | 9,393 × 521 |
| Tree alignment | 9,393 × 436 |
| ASR / DCA alignment | 9,393 × 472 |
| Structure panel | 35 structures（30 AFDB + 5 PDB） |
| S1 ASR state file | 774 MB |
| S2 rooted tree | 已存在，但未做 ingroup prune / ASR |
| S3 midpoint tree | 已存在 |
| S4 provisional reroot | 已存在，但 provenance 冲突未修复 |
| Phase 5 registry | 仅 placeholder，未正式评估节点 |
| DCA outputs | 无 |

---

## 场景命名与允许用途

| Scenario | 定义 | 当前状态 | 允许用途 |
|---|---|---|---|
| S1 | KDOPS outgroup + MFP rooted tree | ✅ 完整 | working tree、基线 ASR、对照 |
| S2 | KDOPS outgroup + LG+C20 rooted tree | ⚠️ 树完整，ASR 缺失 | 仅在补齐 prune + ASR 后用于敏感性比较 |
| S3 | ingroup midpoint reroot | ✅ 完整 | outgroup-free 对照 |
| S4 | ingroup reroot proxy（legacy filename 仍含 `MAD`） | 🔴 provenance blocked | 修复前不得用于正式 gate |
| S5 | reduced representative set + nonreversible/rootstrap | ⬜ 条件触发 | 仅在 S1–S4 仍无法收敛时启动 |

---

## 解释层级（Claim Tiers）

| 层级 | 进入哪里 | 判定标准 |
|---|---|---|
| `root_robust` | 主文 Results | 多 scenario、strict/relaxed 或关键敏感性检查下方向一致 |
| `root_sensitive` | Supplement / Discussion | 结论依赖 rooting 选择，但 provenance 无缺陷 |
| `annotation_sensitive` | Supplement | 结论主要依赖 strict/relaxed 边界定义 |
| `exploratory` | Outlook / 附录 | 数据不足但保留为探索线索 |
| `blocked` | 不写入结果 | provenance 不清、日志不全或 gate 未完成 |

---

## 近期必须产出的文件

### Root / ASR 线

- `results/04_phylogeny_asr/root_scenarios.tsv`（修复并补完整 provenance）
- `results/04_phylogeny_asr/CoreTree_rooted_LGC20_ingroup.treefile`
- `results/04_phylogeny_asr/ASR_core_S2.state`
- `results/04_phylogeny_asr/asr_sensitivity_summary.tsv`
- `results/04_phylogeny_asr/module_origin_stability.tsv`

### Module annotation 线

- `results/03_msa_modules/module_definition_spec.md`
- `results/03_msa_modules/module_boundary_sensitivity.tsv`

### DCA 线

- `results/06_dca/core_dca.afa`
- `results/06_dca/core_dca_stats.tsv`
- `results/06_dca/qc_core_dca.md`

### Phase 5 准备线

- `results/05_struct_md/assembly_adjudication.tsv`
- `results/04_phylogeny_asr/node_selection_registry.tsv`（仍可为空，但状态字段必须真实）

---

## 当前不允许做的事

- 不允许继续把 S4 当作正式 MAD 证据写进结论
- 不允许在没有 S2 ASR 的情况下讨论 deepest ancestral residue 的稳健性
- 不允许在没有 `module_origin_stability.tsv` 的情况下写模块 gain/loss 时间顺序
- 不允许在 QC3 未修复前提名 Tier-1 / Tier-2 祖先节点
- 不允许让 README / TASKS / PLAN / summary / log 彼此漂移

---

## 项目目录（建议以此为准）

```text
results/
  01_mining/
  02_qc/
  03_msa_core/
  03_msa_modules/
  03_msa_full/
  04_phylogeny_asr/
  05_struct_md/
  06_dca/
  meta/
scripts/
meta/
data/
workflow/
```

> `results/05_struct_md/` 作为 Phase 5 的规范目录名。旧文档中的 `05_struct_valid/` 不再作为主命名使用。

---

## 环境

```bash
conda activate dah7ps_v4
# HMMER, MAFFT, IQ-TREE, MMseqs2, FoldMason, SeqKit, CD-HIT, ClipKIT, PastML, plmc, GROMACS
```

---

## 使用原则

1. **所有数字以 `results/meta/metrics_manifest.tsv` 为准。**
2. **所有 scenario 必须有 treefile + summary + log + source script + md5。**
3. **所有深历史叙事都必须先经过 QC3。**
4. **working tree 可以继续算，但不能提前写死历史故事。**

---

## License

Academic research use.
