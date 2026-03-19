# DAH7PS V5.1 Progress Snapshot

> 最后更新：2026-03-18  
> 所有数字以 `results/meta/metrics_manifest.tsv` 为准。

---

## Phase 状态总览

| Phase | 状态 | 关键产物 |
|---|---|---|
| 0 环境与可复现性 | ✅ 完成 | `params.json`, `software_versions.tsv`, `model_files.tsv` |
| 1 数据挖掘 | ✅ 完成 | PASS 24,202 seqs (Ia=9,204 / Ib=6,401 / II=8,597) |
| 2 QC 与去冗余 | ✅ 完成 | NR80=9,673 / seeds60=1,878 / stepping-stone=258 |
| 3.1–3.8 结构感知核心 MSA | ✅ 完成 | `core_global_matchonly.afa` (9,393×521), `core_tree.afa` (436 cols), `core_asr.afa` (472 cols), 5 类模块 strict/relaxed 双矩阵 |
| 3.9 Full-length stitching | ✅ 完成 | `msa_full_Ib_v4.afa` (47×1040), QC2b 断言通过 |
| 4.1 Rooted tree | 🏃 进行中 | MFP ✅ (Q.PFAM+F+R10, LogL=-2835379.61); **LG+C20 正在本地运行** |
| 4.2 AA vs 3Di tree | ✅ 完成 | AA: Q.PFAM+I+R4; 3Di: Q.3Di.AF+G4; nRF=0.74 但亚型单系性一致 → **QC3-YELLOW** |
| 4.3 Core ASR (S1) | ✅ 完成 | `ASR_core.*` (LogL=-2904285.44) |
| 4.3 Core ASR (S2+) | ⬜ 待做 | 需等 LG+C20 树完成 |
| 4.4 Nested ASR (Iβ-ACT) | ⬜ 待做 | 可立即启动 |
| 4.6 Module trait ASR | ⬜ 待做 | 可立即启动 (S1) |
| QC3 Root robustness gate | ⬜ 待做 | 需 S1+S2+S3+S4 全部完成 |
| 5 结构验证 | ⬜ 待做 | 需先通过 QC3 |
| 6 Core DCA | ⬜ 待做 | 6.1 可立即启动 |
| 7 论文蓝图 | ⬜ 待做 | — |

## 当前瓶颈

- **LG+C20 树（S2）正在本地计算中**，候选树搜索阶段，预计数天。
- Root robustness（QC3）无法在 S2 完成前锁定。
- **Working tree (S1) 可继续支撑 ASR、trait 准备、DCA 输入准备。**

## 当前可并行任务

1. S3 midpoint / S4 MAD rooting（仅需 ingroup tree）
2. 4.2 AA vs 3Di 比较
3. 4.4 Iβ-ACT nested ASR
4. 4.6 trait ASR 数据准备（S1）
5. 6.1 core DCA 输入准备 + Meff/L
6. 5.0 assembly adjudication 文献准备
