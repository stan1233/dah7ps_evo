# AGENTS.md — DAH7PS V6.1 执行契约

> **唯一真相源**：`PLAN.md`  
> 本文件只负责把 V6.1 计划翻译成执行顺序、门控规则和禁止事项。

---

## 0. 文档层级

1. `PLAN.md`
2. `results/meta/metrics_manifest.tsv`
3. `results/meta/progress_snapshot.md`
4. `AGENTS.md`
5. `TASKS.md`
6. `log.md`

---

## 1. 当前状态

- Phase 0–3.9：完成
- 已新增 `provenance repair`
- S1 rooted tree + ASR：完成
- S2 rooted tree：完成
- S2 prune + assert_tip_match + ASR：未完成
- `S4` 已正式拆成 `S4a` / `S4b`
- QC3：`HOLD`
- Phase 5：`HOLD`

---

## 2. 总原则

### 2.1 S2 ASR 之前冻结三件事

在 `S2 prune + assert_tip_match + ASR` 完成前：

- 不写 deep-root narrative
- 不定稿 trait ASR
- 不锁定 Phase 5 节点

### 2.2 provenance 先于解释

若没有：

- artifact manifest
- explicit triad
- md5 / commit / command / input md5

则该 scenario 不得进入 narrative gate。

### 2.3 S4 tie 不得合并

如果 `S4a` 与 `S4b`：

- rho 相同
- root identity 不同

则必须保留为两个 scenario。

### 2.4 模块 trait 不再以 residual bin 为主

- trait ASR 应基于 orthogonal features
- `C_tail` 仅允许保留为 `secondary / conditional`

### 2.5 QC3 不是单指标

QC3 至少要分成：

1. provenance
2. topology/model sensitivity
3. root tie/identity
4. annotation sensitivity

### 2.6 S5 仅条件触发

只在 repaired S4 + S2 ASR 完成后仍不稳时才启动。

---

## 3. 当前执行顺序

### P0

1. `artifact_manifest.tsv`
2. `root_scenarios.tsv`
3. `S2 prune + assert_tip_match + ASR`
4. `cross_scenario_asr_sensitivity.py`

### P1

5. orthogonal feature calibration
6. QC3 重评估

### P2

7. DCA 准备
8. nested ASR 准备
9. assembly adjudication 准备

### P3

10. 必要时触发 `S5`
11. 之后才讨论 Phase 5 节点

---

## 4. 硬约束

1. 仍然禁止用 `hmmalign --outformat afa` 直接导出核心 AFA。
2. 全长祖先序列仍然只能来自正式 nested ASR。
3. `assert_tip_match.py` 仍然是 ASR 前硬门槛。
4. 默认 20 线程，不使用 `-1`。
5. 每次关键运行必须即时写入 `log.md`。
6. 正式结果不得静默覆盖；需要 `.bak` 或新显式文件名。
7. 需要绘制进化树及相关树图时，统一使用 `dah7ps_ggtree` 环境中的 R / `ggtree` 工作流。
8. 树图绘制脚本统一保存在 `scripts/`，不得把正式绘图脚本放在 `figures/` 或其他临时目录。
9. 正式绘图输出统一保存在 `figures/`，且每次正式出图必须同时生成 `PDF` 和 `PNG` 各一份。

---

## 5. 当前必备文件

### Root / Provenance

- `results/04_phylogeny_asr/artifact_manifest.tsv`
- `results/04_phylogeny_asr/root_scenarios.tsv`
- `results/04_phylogeny_asr/QC3_root_stability.md`

### Trait encoding

- `results/03_msa_modules/module_feature_registry.tsv`
- `results/03_msa_modules/module_feature_matrix.tsv`
- `results/03_msa_modules/panel35_feature_calibration.tsv`

### Node gate

- `results/04_phylogeny_asr/node_selection_registry.tsv`

---

## 6. 必须具备的脚本接口

- `prune_tree.py`
- `assert_tip_match.py`
- `build_artifact_manifest.py`
- `cross_scenario_asr_sensitivity.py`
- `recode_module_features.py`
- `qc_root_stability.py`

---

## 7. Phase 5 门控

任何候选节点在满足以下全部条件前，一律 `hold`：

1. `S2 ASR` 已完成
2. cross-scenario ASR sensitivity 已完成
3. QC3 不再是 `HOLD`
4. annotation sensitivity 可解释
5. assembly adjudication 明确

---

## 8. 最终提醒

V6.1 的成功标准不是找到唯一真根，而是：

- 让 provenance 可追溯
- 让不确定性可量化
- 让主文只承诺可支撑的结论
