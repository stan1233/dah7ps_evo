# DAH7PS V5.1 Progress Snapshot

> 最后更新：2026-04-08  
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
| 4.1 S1 MFP rooted tree | ✅ 完成 | Q.PFAM+F+R10, LogL=-2835379.61, UFBoot=1000 |
| 4.1 S2 LG+C20 rooted tree | ✅ **完成（2026-04-05）** | `CoreTree_rooted_LGC20.treefile`, LogL=-2820118.09 |
| 4.1 S3 Midpoint ingroup | ✅ 完成 | `CoreTree_rooted_midpoint_ingroup.treefile` |
| 4.1 S4 MAD ingroup | ✅ 完成 | `CoreTree_rooted_MAD_ingroup.treefile`, rho=0.163617 |
| 4.2 AA vs 3Di tree | ✅ 完成 | nRF=0.74, 三大亚型单系一致 → QC3-YELLOW |
| 4.3 Core ASR (S1) | ✅ 完成 | `ASR_core.*` (LogL=-2904285.44) |
| 4.3 Core ASR (S2–S4) | ⬜ 待做 | 需先 prune LGC20 tree，补跑 S2/S4 |
| **QC3 Root robustness gate** | ✅ **完成（2026-04-08）** | **YELLOW — root_sensitive** |
| 4.4 Nested ASR (Iβ-ACT) | ⬜ 待做 | 可立即启动 |
| 4.6 Module trait ASR | ⬜ 待做 | 可立即启动 (S1 + S2) |
| 5 结构验证 | ⬜ 待做 | QC3 已通过 YELLOW；可有条件开始 5.0 |
| 6 Core DCA | ⬜ 待做 | 6.1 可立即启动 |
| 7 论文蓝图 | ⬜ 待做 | — |

---

## QC3 关键发现（2026-04-08）

| 比较对 | RF | nRF | 解读 |
|---|---|---|---|
| S1 MFP vs S2 LGC20 | 7748 | 0.413 | **两棵 outgroup-rooted 树差异显著** |
| S1 MFP vs S3 Midpoint | 54 | 0.003 | 几乎完全一致 |
| S1 MFP vs S4 MAD | 52 | 0.003 | 几乎完全一致 |
| S3 Midpoint vs S4 MAD | 2 | 0.000 | 拓扑等同 |

**根位一致性判断（root_partition_ratio）：**
- S1 (MFP): 1.0 — KDOPS 几乎全部在同一根子支（KDOPS 多系）
- S2 (LGC20): 1.0 — 同上
- S3 (Midpoint): 0.326
- S4 (MAD): 0.300

> **核心解读**：S1 与 S2 两棵 outgroup-rooted 树之间 nRF=0.413，表明换用抗 LBA 模型（LG+C20）后拓扑发生了中等程度的变化，但 S1/S2 与 S3/S4 的根位模式存在差异（ratio=1.0 vs ~0.3）。这反映了 KDOPS 多系性导致的 outgroup rooting 不稳定，是预期的结果，也正是 V5.1 多场景管理的核心原因。

## 模块事件 Root-Sensitivity 标签

| 模块 | Strict 比例 | Root-Sensitivity |
|---|---|---|
| N_ext | 33.3% | ⚠️ root_sensitive |
| alpha2beta3_insert | 1.8% | ✅ root_robust（极稀有，解释不依赖深根） |
| ACT_domain | 0.5% | ✅ root_robust（极稀有）|
| CM_domain | 4.3% | ⚠️ root_sensitive |
| C_tail | 3.8% | ⚠️ root_sensitive |

> **注意**：以上标签是基于流行率的启发式估计。精确的 root-sensitivity 判定必须等 Phase 4.6 trait ASR 在 S1 + S2 树上运行完毕后才能锁定。

## 当前瓶颈

- **QC3 已通过（YELLOW）**，Phase 5 可有条件推进。
- S2 (LGC20) ingroup tree 尚未 prune + assert_tip_match，需尽快完成以支持 ASR 补跑。
- Phase 4.6 trait ASR 是将 module event root-sensitivity 从"启发式"升级为"定量确认"的关键步骤。

## 当前可并行任务

1. **S2 ingroup tree 准备**（prune KDOPS → assert_tip_match → 归档）
2. **4.4 Iβ-ACT nested ASR**（输入 `msa_full_Ib_v4.afa` 已就绪）
3. **4.6 trait ASR 数据准备**（strict/relaxed × S1 + S2）
4. **6.1 core DCA 输入准备** + Meff/L 审计
5. **5.0 assembly adjudication 文献准备**（QC3 YELLOW，条件性开放）
