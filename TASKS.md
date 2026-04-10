# TASKS.md — DAH7PS V6.1 Live Checklist

> 状态标记：`[x]` 完成 / `[/]` 进行中 / `[ ]` 待做 / `[-]` 被 gate 阻塞

---

## 0. 已执行的 V6.1 调整

- [x] 插入 `provenance repair` 小阶段
- [x] 建立 `results/04_phylogeny_asr/artifact_manifest.tsv`
- [x] 给 S1–S4 正式/暂定产物登记 `md5 + script + commit + command + input md5 + generated_at + formal/provisional`
- [x] 将旧单一 `S4` 拆为 `S4A_TOP500_PROXY` 与 `S4B_FULLSEARCH_PROXY`
- [x] 重写 `root_scenarios.tsv`
- [x] 重新生成并同步 S4a/S4b summary/log triads
- [x] 将 `S2 prune + assert_tip_match + ASR` 提升为绝对第一优先级
- [x] 新增 `scripts/cross_scenario_asr_sensitivity.py`
- [x] 新增 orthogonal trait encoding schema
- [x] 生成 `module_feature_registry.tsv`
- [x] 生成 `panel35_feature_calibration.tsv`
- [x] 重写 QC3 为四维 gate
- [x] 将 `progress_snapshot.md` 改为 V6.1 口径

---

## 1. P0 Blocking Path

### 1.1 S2 prune + assert_tip_match + ASR

- [ ] prune `CoreTree_rooted_LGC20.treefile`
- [ ] 对 `CoreTree_rooted_LGC20_ingroup.treefile` 运行 `assert_tip_match.py`
- [ ] 补跑 `ASR_core_S2.state / .iqtree / .treefile / .log`
- [ ] 将 S2 从 `tree-only` 升级到 `ASR-ready`
- [ ] 把结果写回 `artifact_manifest.tsv`
- [ ] 更新 `metrics_manifest.tsv`
- [ ] 更新 `log.md`

### 1.2 Cross-scenario ASR sensitivity

- [ ] 建立 node mapping 或确认共用 node labels
- [ ] 跑 `S1/S2/S4A/S4B` site-level 比较
- [ ] 输出 MAP consistency
- [ ] 输出 posterior delta
- [ ] 输出 information bits / information range
- [ ] 生成 node-level summary
- [ ] 明确区分 `model-sensitive` / `root-sensitive` / `stable`

---

## 2. P1 High Priority

### 2.1 Orthogonal trait encoding + calibration

- [x] 从 legacy module matrix 拆出 orthogonal features
- [x] 将 `C_tail` 降级为 `c_residual` secondary character
- [x] 建立 `panel35_feature_calibration.tsv`
- [ ] 完成 35-panel 人工校准
- [ ] 形成 `annotation_sensitivity` 报告
- [ ] trait ASR 改为基于 orthogonal features，而不是 residual class

### 2.2 QC3 closure

- [x] provenance dimension
- [x] topology/model sensitivity dimension
- [x] root tie/identity dimension
- [x] annotation sensitivity dimension
- [ ] 在 S2 ASR 完成后重跑 QC3
- [ ] 只在四维 gate 都过线后解除 HOLD

---

## 3. Conditional Trigger

### 3.1 S5 nonreversible/rootstrap

- [ ] 在 repaired S4 + S2 ASR 后重新判断最深层结论是否仍不稳
- [ ] 若仍不稳，设计 reduced representative set
- [ ] 运行 nonreversible rooting
- [ ] 运行 rootstrap
- [ ] 结论按 `main-text eligible` / `conditional only` / `reject` 分层

---

## 4. Parallel Safe Work

- [ ] core DCA 输入准备
- [ ] `Meff/L` 审计
- [ ] `Ibeta-ACT` nested full-length ASR 准备
- [ ] assembly adjudication 文献整理

---

## 5. Explicitly Frozen

- [-] deep-root narrative
- [-] trait ASR final wording
- [-] Phase 5 node locking
- [-] 用 `C_tail` 直接写主 evolutionary character
- [-] 用单一 `S4` 代表全部 reroot evidence
