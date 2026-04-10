# DAH7PS Allosteric Evolution Project

> **Status: V6.1 workflow realigned on 2026-04-10**  
> 当前口径：`计算准备 = GO`，`deep-root narrative = HOLD`，`trait ASR 定稿 = HOLD`，`Phase 5 节点锁定 = HOLD`

---

## 当前项目结论

本项目没有回到 Phase 0–3 重做。  
当前真正的调整是把 Phase 4 前半段改成：

1. 先做 `provenance repair`
2. 把 `S2 prune + assert_tip_match + ASR` 提到绝对第一优先级
3. 新增 `cross-scenario ASR sensitivity`
4. 把模块 trait 编码改成 `orthogonal features + 35-panel calibration`
5. 把 QC3 改成四维 gate，而不是单一 `root_partition_ratio`
6. 只有在 repaired S4 与 S2 ASR 之后仍不稳时，才触发 `S5 nonreversible/rootstrap`

---

## 当前主状态

| 层级 | 状态 | 说明 |
|---|---|---|
| Phase 0–3.9 | ✅ 完成 | 不返工 |
| S1 rooted tree + ASR | ✅ 完成 | 当前唯一完整 ASR 基线 |
| S2 rooted tree | ✅ 完成 | 但 `prune + assert_tip_match + ASR` 未完成 |
| S3 midpoint | ✅ 完成 | provenance repair 已补 summary/log |
| S4 | [/] repaired | 已正式拆成 `S4a_top500` / `S4b_fullsearch` |
| artifact manifest | ✅ 已建立 | `results/04_phylogeny_asr/artifact_manifest.tsv` |
| cross-scenario ASR sensitivity | ⬜ 待运行 | 脚本已补，等 S2 ASR |
| orthogonal trait encoding | [/] 进行中 | registry + matrix + panel35 calibration 已建 |
| QC3 | ⛔ HOLD | 多维 gate 未关闭 |
| Phase 5 | ⛔ HOLD | 禁止锁节点 |

---

## Root Scenario 现定义

| Scenario | 身份 | 当前地位 |
|---|---|---|
| `S1_MFP_KDOPS` | KDOPS outgroup + MFP | formal |
| `S2_LGC20_KDOPS` | KDOPS outgroup + LG+C20 | formal |
| `S3_MIDPOINT_INGROUP` | midpoint ingroup reroot | formal |
| `S4A_TOP500_PROXY` | ingroup reroot proxy, top-500 | provisional |
| `S4B_FULLSEARCH_PROXY` | ingroup reroot proxy, full-search | provisional |
| `S5_NONREV_ROOTSTRAP_REDUCED` | reduced nonreversible/rootstrap | conditional trigger |

**关键事实：**

- `S4a` 与 `S4b` 的 `rho` 同为 `0.163617`
- 但两者不是同一个 root identity
- 因此不能再强行保留一个单一 `S4`

---

## 当前新增产物

### Provenance

- `results/04_phylogeny_asr/artifact_manifest.tsv`
- `results/04_phylogeny_asr/root_scenarios.tsv`
- `results/04_phylogeny_asr/CoreTree_rooted_S3_midpoint.summary.txt`
- `results/04_phylogeny_asr/CoreTree_rooted_S4a_top500.summary.txt`
- `results/04_phylogeny_asr/CoreTree_rooted_S4b_fullsearch.summary.txt`

### ASR / QC

- `scripts/cross_scenario_asr_sensitivity.py`
- `scripts/qc_root_stability.py` (V6.1 multi-dimensional gate)
- `results/04_phylogeny_asr/QC3_root_stability.md`

### Trait Encoding

- `scripts/recode_module_features.py`
- `results/03_msa_modules/module_feature_registry.tsv`
- `results/03_msa_modules/module_feature_matrix.tsv`
- `results/03_msa_modules/panel35_feature_calibration.tsv`

---

## 现在绝对先做什么

```bash
conda run -n dah7ps_v4 python scripts/prune_tree.py \
  --input results/04_phylogeny_asr/CoreTree_rooted_LGC20.treefile \
  --remove_prefix KDOPS_ \
  --output results/04_phylogeny_asr/CoreTree_rooted_LGC20_ingroup.treefile \
  --assert_rooted

conda run -n dah7ps_v4 python scripts/assert_tip_match.py \
  --tree results/04_phylogeny_asr/CoreTree_rooted_LGC20_ingroup.treefile \
  --msa results/03_msa_core/core_asr.afa \
  --assert_identical

iqtree -s results/03_msa_core/core_asr.afa \
  -te results/04_phylogeny_asr/CoreTree_rooted_LGC20_ingroup.treefile \
  -m LG+C20+F+G --asr -T 20 \
  --prefix results/04_phylogeny_asr/ASR_core_S2
```

随后再跑：

```bash
python scripts/cross_scenario_asr_sensitivity.py \
  --state S1=results/04_phylogeny_asr/ASR_core.state \
  --state S2=results/04_phylogeny_asr/ASR_core_S2.state \
  --state S4A=results/04_phylogeny_asr/ASR_core_S4a.state \
  --state S4B=results/04_phylogeny_asr/ASR_core_S4b.state \
  --out_prefix results/04_phylogeny_asr/asr_cross_scenario
```

---

## 当前禁止事项

- 不得再写 “QC3 已通过，可以开始 Phase 5”
- 不得在 S2 ASR 前定稿 deepest-root narrative
- 不得在 `C_tail` residual bin 上直接讲主 evolutionary story
- 不得锁定 Phase 5 节点
- 不得把 `S4b` 当唯一真根

---

## Manuscript 边界

主文只允许：

1. 稳定基础层
2. `root_robust` 结论
3. 与之同向的 DCA / 结构支持

补充或 Discussion：

1. `root_sensitive`
2. `annotation_sensitive`
3. `S4 tie`
4. 若触发 `S5` 后仍 unresolved 的 deepest-root 争议
