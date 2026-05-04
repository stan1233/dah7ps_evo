# DAH7PS V6.1 Progress Snapshot

> 最后更新：2026-04-29  
> 所有数字以 `results/meta/metrics_manifest.tsv` 为准。  
> 所有 provenance 以 `results/04_phylogeny_asr/artifact_manifest.tsv` 与 `results/04_phylogeny_asr/root_scenarios.tsv` 为准。

---

## Phase 状态总览

| Phase | 状态 | 关键产物 / 现状 |
|---|---|---|
| 0 环境与可复现性 | ✅ 完成 | `params.json`, `software_versions.tsv`, `model_files.tsv` |
| 1 数据挖掘 | ✅ 完成 | PASS 24,202 seqs |
| 2 QC 与去冗余 | ✅ 完成 | NR80=9,673 / seeds60=1,878 / stepping-stone=258 |
| 3.1–3.8 结构感知核心 MSA | ✅ 完成 | `core_global_matchonly.afa`, `core_tree.afa`, `core_asr.afa` |
| 3.9 Full-length stitching | ✅ 完成 | `msa_full_Ib_v4.afa`, column drift pass |
| 4.0 provenance repair | [/] 进行中 | `artifact_manifest.tsv` 已建立；S4 正式拆为 `S4a_top500` / `S4b_fullsearch` |
| 4.1 S1 rooted tree + ASR | ✅ 完成 | `CoreTree_rooted_MFP.*`, `ASR_core.*` |
| 4.1 S2 rooted tree | ✅ 完成 | `CoreTree_rooted_LGC20.*` |
| 4.1 S2 prune + assert_tip_match + ASR | ⛔ HOLD / DEFERRED | 按 2026-04-29 interim 决策暂不重跑 S2；旧 S2 post-hoc noO66496 仅用于 visualization/provenance consistency |
| 4.1 S3 midpoint | ✅ 完成 | `CoreTree_rooted_midpoint_ingroup.treefile`, repaired summary/log |
| 4.1 S4 reroot proxy | [/] repaired | `S4a` / `S4b` 显式分拆；两者同 rho 不同 root identity |
| 4.2 AA vs 3Di tree | ✅ 完成 / metadata repaired | `2B7O/3NV8/5CKV` 已标记为同一 O53512 anchor group；AA/3Di 作为 QC3 正交证据保留 |
| 4.x cross-scenario ASR sensitivity | ⬜ 待做 | 脚本接口已补；待 S2 ASR 后运行 |
| 4.x orthogonal trait encoding | [/] 进行中 | `module_feature_registry.tsv`, `module_feature_matrix.tsv`, `panel35_feature_calibration.tsv` 已生成；`pdb_anchor_registry.tsv` 已记录 O53512 state variants |
| QC3 multi-dimensional gate | ⛔ HOLD | provenance=PASS; topology/model=HOLD; root tie=HOLD; annotation=HOLD |
| 5 结构验证 | ⛔ HOLD | S2 ASR 未完成前不得锁定节点，不得推进主线 Phase 5 |
| 6 Core DCA | ⬜ 可并行准备 | 允许准备输入与 Meff/L 审计；不受 deep-root narrative 影响 |
| 7 论文蓝图 | ⛔ HOLD | 主文只允许 root-robust 结论；当前 deepest root 相关叙事冻结 |

---

## 当前关键变化（V6.1）

1. `provenance repair` 已被插入为正式小阶段，位置在 S1–S4 结果解释之前。
2. 原单一 `S4` 已拆分为：
   - `S4A_TOP500_PROXY`
   - `S4B_FULLSEARCH_PROXY`
3. 旧版“QC3 已完成 / YELLOW / Phase 5 可条件开始”口径已废止。
4. trait ASR 主编码从 legacy 5-module 黑箱转向 orthogonal feature registry。
5. `C_tail` 已被降级为 residual secondary character，不得直接充当主要 evolutionary character。

---

## QC3 现状态（2026-04-10）

| 维度 | 状态 | 当前判定 |
|---|---|---|
| provenance | PASS | S1/S2/S3/S4a/S4b 三联件与 md5/commit/command 已登记 |
| topology / model sensitivity | HOLD | S1 vs S2 nRF=0.413，模型敏感性仍显著 |
| root tie / identity | HOLD | S4a 与 S4b rho 同为 0.163617，但根身份不同 |
| annotation sensitivity | HOLD | orthogonal feature schema 已建立，但 35-panel 尚未人工校准 |

**结论：QC3 当前为 `HOLD`，不是 `YELLOW GO`。**

---

## 当前绝对优先级

1. 等待正在运行的 `S1_KDOPS11_noO66496` 正式 MFP rerun 完成，再与 legacy S1 ingroup tree 比较。
2. 旧 S2 仅保留 `legacy12_posthoc_noO66496`，不得视为 KDOPS11 S2 rerun。
3. 35-panel feature calibration + annotation sensitivity，尤其是 O53512 PDB anchor group。
4. 若未来要解除 QC3 或推进 Phase 5，仍需正式解决 S2/KDOPS11 model-sensitivity 缺口。

---

## 当前禁止事项

- 不得继续写 deepest-root narrative。
- 不得定稿 trait ASR。
- 不得锁定 Phase 5 节点。
- 不得把 `S4b` 或任何 single-root scenario 当作唯一真相。
- 不得把 `C_tail` residual bin 当作主 evolutionary character。

---

## 当前可并行任务

1. `Ibeta-ACT` nested full-length ASR 的输入准备
2. core DCA 输入准备与 `Meff/L` 审计
3. assembly adjudication 的文献与模板整理
4. `cross_scenario_asr_sensitivity.py` 的 node-map 准备
