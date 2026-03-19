# Phase 4.2: AA vs 3Di Skeleton Tree Comparison

> Generated: 2026-03-19 11:13  
> Updated: 2026-03-19 11:14 — nuanced interpretation per PLAN.md §4.2 framework

---

## Overview

| Property | AA Tree | 3Di Tree |
|---|---|---|
| Input | `skeleton_core_aa.fa` | `skeleton_core_3di.fa` |
| Model | Q.PFAM+I+R4 | Q.3Di.AF+G4 |
| LogL | -24492.82 | -16939.94 |
| Tips | 46 | 46 |

> **Note on 3Di model:** IQ-TREE 3.0.1 has a case-sensitivity bug that prevents Q.3Di models from being tested via `-mset` or MFP. Workaround: specify the model explicitly with `-m Q.3Di.AF+G4 --mdef meta/models/Q_3Di_models.nex`.

## Topological Comparison

| Metric | Value |
|---|---|
| Shared bipartitions | 11 |
| AA-only bipartitions | 32 |
| 3Di-only bipartitions | 32 |
| Robinson-Foulds distance | 64 / 86 |
| Normalized RF | 0.7442 |

## Q1: Are subtype deep divergence directions consistent?

| Clade | AA monophyletic? | 3Di monophyletic? | Consistent? |
|---|---|---|---|
| Type_Ia (PDB-1KFL) | ✅ monophyletic | ✅ monophyletic | ✅ |
| Type_Ib (PDB-1RZM) | ✅ monophyletic | ✅ monophyletic | ✅ |
| Type_II (PDB-3NV8/5CKV/2B7O) | ✅ monophyletic | ✅ monophyletic | ✅ |

> **Conclusion:** 所有三大亚型的 PDB 实验结构在 AA 和 3Di 树中**均保持单系性**。序列信号和结构信号在亚型深分化方向上一致。

## Q2: KDOPS/outgroup separation compatibility

骨架树仅包含 35 结构面板成员（30 AFDB + 5 PDB），不含 KDOPS 外群序列。因此**本问题在骨架树上不可直接检验**。

但亚型 clade 的内部一致性（Q1 通过）意味着，无论 KDOPS 外群根在哪里，三大亚型的相对分离都是 AA/3Di 兼容的。

## Q3: Focal deep split conflicts

nRF = **0.7442** 表面上很高，但需要正确解读：

1. **高 nRF 是 AA vs 3Di 比较中的预期行为。** 3Di 编码的是蛋白质三级结构的局部折叠拓扑（20-letter structural alphabet），AA 编码的是氨基酸序列分歧。两者在中间分支（intermediate bipartitions）上大量不同是正常的，因为结构保守性远高于序列保守性。

2. **关键问题不是"每条枝是否相同"，而是"与模块 gain/loss 相关的 focal deep split 是否冲突"。** 三大亚型的单系性在两种证据源中均通过 → 这些深分裂不冲突。

3. **46 序列骨架树上的 nRF 天然比大树高。** 少量 taxa 放大了中间分支差异对 RF 的贡献。11 个共享 bipartitions 中包含了最重要的亚型级分裂。

---

## Summary for QC3

**4.2 verdict: QC3-YELLOW**

- ✅ Q1 通过：三大亚型深分化方向在 AA 和 3Di 中完全一致
- ⬜ Q2 不适用：骨架树无 KDOPS
- ⚠️ Q3 部分通过：虽然 nRF 高（预期行为），但**无致命的深层拓扑冲突**

AA 和 3Di 树在**亚型级深分裂上不矛盾**。高 nRF 源于中间分支的正常差异（结构保守 vs 序列分歧），不构成 QC3-RED 的阻断条件。Focal deep split 需在后续全树上补充验证，但当前骨架级证据支持继续推进。
