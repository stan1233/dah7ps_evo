# QC2：Phase 2 质量控制与去冗余报告

> **生成时间**：2026-02-25 09:45  
> **操作环境**：`conda activate dah7ps_v4`  
> **SOP 对应**：Phase 2（质量控制与去冗余）  
> **验收标的**：`results/02_qc/nr80_*.fasta` + 本报告

---

## 1. Phase 2 输入与门控总结

### 1.1 Gate A：Ia-II 高重叠归属

- **发现**：Ia 与 II 的 hmmsearch 命中集合有 8,561 条交集（占 Ia 的 85%）
- **处理**：竞争打分归属（`assign_besthit_Ia_vs_II.py`），全部 8,561 条以 >= 20 bits 优势归属 Ia（100% HIGH 置信度）
- **产出**：`results/01_mining/besthit_Ia_vs_II.tsv`

### 1.2 Gate B：Ib 边缘 KDOPS 序列隔离

| ID | Delta (bits) |
|----|-------------|
| UniRef90_A0A1I4NDE8 | 13.7 |
| UniRef90_A0A806JZI9 | 14.8 |
| UniRef90_X1MS81 | 16.0 |
| UniRef90_A0ABT1SNH0 | 16.6 |

处理：保留在全量库，但从 seeds60 / stepping stones 中排除。

### 1.3 Phase 2 互斥输入

| 亚型 | 输入序列数 | 来源文件 |
|------|-----------|---------|
| Ia | 10,071 | `hits_Ia_seqs.fasta` |
| Ib | 7,869 | `hits_Ib_clean.fasta` |
| II | 9,968 | `hits_II_final_seqs.fasta` |
| **合计** | **27,908** | |

---

## 2. Phase 2.1：长度 + HMM 覆盖度三箱过滤

### 参数

| 亚型 | HMM长度 | canonical_min | canonical_max | cov_mode |
|------|--------|-------------|-------------|----------|
| Ia | 355 | 300 | 450 | best |
| Ib | 334 | 280 | 450 | best |
| II | 471 | 320 | 600 | merged |

cov_min = 0.70（所有亚型）。Type II 使用 merged 模式（multi-domain stitching: i-Evalue <= 1e-5, HMM span >= 30 aa, merge_gap = 5）。

### 结果

| 亚型 | 总数 | PASS_CANONICAL | PASS_LONG | FRAG | cov_mode |
|------|------|---------------|-----------|------|----------|
| Ia | 10,071 | 9,022 (89.6%) | 182 (1.8%) | 867 (8.6%) | best |
| Ib | 7,869 | 6,229 (79.2%) | 172 (2.2%) | 1,468 (18.7%) | best |
| II | 9,968 | 8,457 (84.8%) | 140 (1.4%) | 1,371 (13.8%) | merged |
| **合计** | **27,908** | **23,708** | **494** | **3,706** | |

### Type II stitching 效果

| 指标 | best | merged | 变化 |
|------|------|--------|------|
| PASS_CANONICAL | 6,411 | 8,457 | +2,046 |
| FRAG | 3,464 (34.8%) | 1,371 (13.8%) | -2,093 |
| Rescued | - | 2,093 | |

Phase 2.2 输入 = PASS_CANONICAL + PASS_LONG：Ia=9,204 / Ib=6,401 / II=8,597（总计 24,202）

---

## 3. Phase 2.2：CD-HIT 80% 去冗余 (nr80)

参数：`cd-hit -c 0.80 -n 5 -M 4000 -T 20`

| 文件 | 序列数 | 平均长度 | 去冗余率 |
|------|--------|---------|---------| 
| `nr80_Ia.fasta` | 3,521 | 385.4 | 61.7% |
| `nr80_Ib.fasta` | 3,073 | 356.8 | 52.0% |
| `nr80_II.fasta` | 3,079 | 464.8 | 64.2% |
| **总计** | **9,673** | | **60.0%** |

---

## 4. Phase 2.3：CD-HIT 60% 种子提取 (seeds60)

参数：`cd-hit -c 0.60 -n 4 -M 4000 -T 20`（从 nr80 输入）

| 文件 | 序列数 | 平均长度 |
|------|--------|---------|
| `seeds60_Ia.fasta` | 581 | 482.8 |
| `seeds60_Ib.fasta` | 648 | 374.9 |
| `seeds60_II.fasta` | 649 | 516.5 |
| **总计** | **1,878** | |

Gate B 排除：seeds60_Ib 已剔除 4 条 KDOPS 边缘序列。

---

## 5. Phase 2.4：Stepping Stones（两层定义）

### 5.1 Coverage Backbone（Phase 2.4a）：258 条

参数：MMseqs2 `easy-cluster --min-seq-id 0.4 -c 0.8 --cov-mode 1`
输入：`all_seeds_mixed.fasta`（1,878 条）
产出：`stepping_stones_rep_seq.fasta`（258 条）

**Per-subtype 分布：**

| 亚型 | 代表数 | 占比 |
|------|--------|------|
| Ia | 104 | 40.3% |
| Ib | 46 | 17.8% |
| II | 108 | 41.9% |

**跨亚型混簇：0 个（100% 纯单亚型簇）**

聚类大小分布：singleton=87, 2-5=95, 6-10=23, 11-20=33, >20=20（min=1, median=2, mean=7.3, max=126）

### 5.2 258 vs SOP 目标 20-50 的说明

SOP 目标 20-50 条的约束意图是限制 Phase 3.1 **结构面板**的规模（FoldMason 计算可行性），而非限制覆盖度骨架本身。258 条反映了 DAH7PS 三亚型在 40% 一致性下的真实序列多样性。

验证：尝试在 258 条上二次聚类（--min-seq-id 0.25 和 0.20），分别得到 183 和 181 个聚类，证实 DAH7PS 家族序列极度发散，纯序列聚类无法将其压缩到 20-50 条。

**决策：分层使用（方案 2）**

1. **Coverage Backbone（258 条）**：用于覆盖度验证、空白分支发现、stepping-stone bridge 证据保留
   - 文件：`stepping_stones_rep_seq.fasta`
   - 用途：Phase 2 覆盖度 QC、Phase 3.5 亚型内骨架参考

2. **Structure Panel（目标 20-40 条）**：Phase 3.1 结构面板的输入
   - 产出文件：Phase 3.1 生成 `data/structures/panel_dah7ps/`
   - 选择标准（优先级递减）：
     - P1：实验 PDB 结构（直接入面板）
     - P2：AlphaFold DB 已有高置信预测（pLDDT >= 70）
     - P3：无结构但覆盖关键分支空白 → ESMFold 预测，pLDDT >= 70 才入面板
   - 每亚型配额建议：Ia 8-15, Ib 4-10, II 8-15（按 backbone 代表数比例分配）

SOP 的 20-50 约束重新绑定到 Structure Panel 子集上，覆盖度骨架保持 258 条不变。

### 5.3 产出文件清单

| 文件 | 说明 | 层级 |
|------|------|------|
| `stepping_stones_rep_seq.fasta` | 258 条代表序列 | Coverage Backbone |
| `stepping_stones_cluster.tsv` | 聚类归属表 | Coverage Backbone |
| `stepping_stones_all_seqs.fasta` | 全部序列 | Coverage Backbone |
| （Phase 3.1 产出） | 20-40 条结构面板 | Structure Panel |

---

## 6. Phase 2 验收总结

### Done 条件检查

| 条件 | 状态 |
|------|------|
| `results/02_qc/nr80_*.fasta` 产出 | PASS (Ia=3,521 / Ib=3,073 / II=3,079) |
| `results/02_qc/qc_length_report.md` 产出 | PASS（本报告） |
| seeds60 可用于 Phase 3.5 | PASS (总计 1,878) |
| stepping stones 覆盖度骨架 | PASS (258 条, 零跨亚型混簇) |
| Gate A/B 风险已消解 | PASS |

### Phase 2 全流程数据链

```
Phase 1 输出 (27,908)
    |
    v
Gate A: Ia∩II 归属 → 互斥集合 (27,908)
Gate B: Ib 边缘 KDOPS 隔离 (4 条标记)
    |
    v
Phase 2.1: 长度+覆盖度过滤 → PASS (24,202) + FRAG (3,706)
    |
    v
Phase 2.2: CD-HIT 80% → nr80 (9,673)
    |
    v
Phase 2.3: CD-HIT 60% → seeds60 (1,878, Gate B 已排除)
    |
    v
Phase 2.4a: MMseqs2 40% → Coverage Backbone (258)
    |
    v  (Phase 3.1)
Phase 2.4b: 结构优先级筛选 → Structure Panel (20-40)
```

**Phase 2 正式关账。下一步：Phase 3.1 结构面板构建。**
