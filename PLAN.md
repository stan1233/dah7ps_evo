# DAH7PS 变构起源演化计划：V6.1

> 修订日期：2026-04-10  
> 本版目标不是“硬找唯一真根”，而是把哪些结论可以进主文、哪些只能保留为条件性结论明确分层。

---

## 一、V6.1 的核心变化

### 1. 插入 provenance repair 小阶段

在解释任何 S1–S4 root scenario 之前，先要求：

- `artifact_manifest.tsv`
- 每个 scenario 的 `treefile + summary + log`
- md5 / script / commit / command / input md5 / generated_at / formal-or-provisional

### 2. S2 ASR 上升为绝对第一优先级

在以下步骤完成前，禁止锁定深根叙事：

1. prune `S2`
2. `assert_tip_match`
3. `ASR_core_S2`

### 3. 新增 cross-scenario ASR sensitivity

不再只比较一条祖先序列。  
必须输出：

- 每个节点
- 每个位点
- MAP 一致性
- posterior 差值
- 信息量

比较对象至少包括：

- `S1`
- `S2`
- `S4a`
- `S4b`

### 4. 模块 trait 改为 orthogonal feature 编码

trait ASR 不应再围绕 legacy 5-module 黑箱直接展开。  
尤其：

- `C_tail` 只保留为 residual secondary character
- 不得再把它当主要 evolutionary character
- 必须用 35 个结构面板做校准

### 5. QC3 改为四维 gate

QC3 至少包含：

1. provenance
2. topology / model sensitivity
3. root tie / identity
4. annotation sensitivity

### 6. S5 是条件触发，不是默认下一步

只有在补完 `S2` 与 repaired `S4` 后，最深层结论仍不稳时，才触发：

- nonreversible rooting
- rootstrap

---

## 二、当前状态

- Phase 0–3.9：完成
- S1 rooted tree + ASR：完成
- S2 rooted tree：完成
- S2 prune + assert_tip_match + ASR：未完成
- S3 midpoint：完成，provenance 已修
- S4：已拆分为 `S4a_top500` 与 `S4b_fullsearch`
- QC3：`HOLD`
- Phase 5：`HOLD`

---

## 三、当前执行顺序

### P0

1. `provenance repair`
2. `S2 prune + assert_tip_match + ASR`
3. `cross-scenario ASR sensitivity`

### P1

4. orthogonal feature calibration
5. trait ASR 改为 feature-based
6. 四维 QC3 重评估

### P2

7. core DCA 准备与 `Meff/L`
8. nested ASR 准备
9. assembly adjudication 准备

### P3

10. 若仍不稳，触发 `S5`
11. 之后才讨论 Phase 5 节点

---

## 四、当前正式 root scenario 集

| Scenario | 说明 | 状态 |
|---|---|---|
| `S1_MFP_KDOPS` | KDOPS outgroup + MFP | formal |
| `S2_LGC20_KDOPS` | KDOPS outgroup + LG+C20 | formal |
| `S3_MIDPOINT_INGROUP` | midpoint reroot | formal |
| `S4A_TOP500_PROXY` | top-500 reroot proxy | provisional |
| `S4B_FULLSEARCH_PROXY` | full-search reroot proxy | provisional |
| `S5_NONREV_ROOTSTRAP_REDUCED` | reduced nonreversible/rootstrap | conditional |

**关键说明：**

- `S4a` 与 `S4b` 当前 `rho` 相同
- 但根身份不同
- 因此必须并列存在

---

## 五、Gate System

### Gate P: Provenance

通过条件：

- `artifact_manifest.tsv` 完整
- S1–S4 当前正式/暂定场景有 triad
- legacy 文件的身份已写清

### Gate M: Model / ASR sensitivity

通过条件：

- `S2` 已完成 prune + assert + ASR
- `S1/S2/S4a/S4b` 已完成 site-level 与 node-level sensitivity 输出

### Gate A: Annotation sensitivity

通过条件：

- orthogonal features 取代 legacy 5-module 直推
- 35-panel calibration 有明确人工复核
- `C_tail` residual 不再承担主结论

### Gate R: Root identity

通过条件：

- repaired S4 的 root tie 已被显式处理
- 不再把同 rho 不同 identity 强行归并

### Gate N: Narrative

通过条件：

- 只有 `root_robust` 结论进主文
- `root_sensitive` / `annotation_sensitive` 仅保留为条件性结论

---

## 六、Phase 5 开放条件

以下全部满足前，Phase 5 保持关闭：

1. `S2 ASR` 已完成
2. cross-scenario ASR sensitivity 已完成
3. QC3 不再是 `HOLD`
4. 节点在多场景下至少具备最小一致性
5. assembly adjudication 有明确 oligomer decision

---

## 七、必备输出

### Root / Provenance

- `results/04_phylogeny_asr/artifact_manifest.tsv`
- `results/04_phylogeny_asr/root_scenarios.tsv`
- `results/04_phylogeny_asr/QC3_root_stability.md`

### ASR

- `results/04_phylogeny_asr/ASR_core_S2.state`
- `results/04_phylogeny_asr/asr_cross_scenario_site.tsv`
- `results/04_phylogeny_asr/asr_cross_scenario_node.tsv`

### Trait encoding

- `results/03_msa_modules/module_feature_registry.tsv`
- `results/03_msa_modules/module_feature_matrix.tsv`
- `results/03_msa_modules/panel35_feature_calibration.tsv`

---

## 八、成功标准

本轮修订成功的定义是：

1. provenance 不再漂移
2. `S2` 从 tree-only 变成 ASR-ready
3. S4 tie 被显式分拆而非隐藏
4. trait ASR 建立在 orthogonal features 上
5. QC3 变成真正 gate
6. 若 deepest root 仍不稳，项目仍能给出分层、可写、可审稿的结论
