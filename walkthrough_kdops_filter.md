# KDOPS 反向过滤 Type Iβ 序列 — 总结

## 背景

DAH7PS（3-deoxy-7-phosphoheptulonate synthase）与 KDOPS（KDO8P synthase）同属 TIM barrel 超家族，序列相似度较高。用 DAH7PS Type Iβ 的 HMM 从 UniRef90 中搜库时，不可避免地会"钓"进一批 KDOPS 序列。需要构建 KDOPS 专属 HMM，对捞出的 Type Iβ 序列做**反向筛选（Negative Selection）**，将 KDOPS 污染序列剔除。

## 流程总览

```
kdo8ps_uniprot.fasta (334 seqs, UniProt 下载)
        │
        ▼
   CD-HIT 90% ──→ kdo8ps_nr90.fasta (119 seqs)
        │
        ▼
  MAFFT L-INS-i ──→ kdo8ps_aligned.afa
        │
        ▼
    hmmbuild ──→ kdo8ps.hmm (280 match states)
        │
        ▼
  ┌─────┴─────┐
  │            │
  ▼            ▼
hmmsearch    hmmsearch
kdo8ps.hmm   model_Ib.hmm
  vs           vs
raw_full_Ib  raw_full_Ib
  │            │
  └─────┬──────┘
        ▼
  比较 full-sequence score
  KDOPS score ≥ DAH7PS score → 剔除
        │
        ▼
  seqkit grep -v ──→ raw_full_Ib_clean.fasta (15,731 seqs)
```

## 具体步骤

### 1. CD-HIT 去冗余

```bash
cd-hit -i kdo8ps_uniprot.fasta -o kdo8ps_nr90.fasta -c 0.9 -n 5 -M 4000 -T 8
```

将 UniProt 下载的 334 条 KDOPS 序列（含 KDSA + KDSC）在 90% 相似度下聚类，得到 **119 条**代表序列。

> 注：KDSC 是磷酸酶，与 DAH7PS 结构无关，混入不影响反向过滤结果。

### 2. MAFFT 多序列比对

```bash
mafft --auto kdo8ps_nr90.fasta > kdo8ps_aligned.afa
```

MAFFT `--auto` 自动选择了 **L-INS-i**（最精确模式），输出比对文件。

### 3. hmmbuild 构建 HMM

```bash
hmmbuild kdo8ps.hmm kdo8ps_aligned.afa
```

得到 KDOPS HMM profile，**280 个 match states**。

### 4. 双向 hmmsearch 扫描

```bash
# KDOPS 模型扫 Type Iβ 序列
hmmsearch --cpu 8 --domtblout domhits_Ib_vs_kdops.tbl kdo8ps.hmm raw_full_Ib.fasta > hmmsearch_kdops_vs_Ib.log 2>&1

# DAH7PS Ib 模型扫 Type Iβ 序列（重新跑一遍，确保打分环境一致）
hmmsearch --cpu 8 --domtblout domhits_Ib_vs_dah7ps.tbl model_Ib.hmm raw_full_Ib.fasta > hmmsearch_dah7ps_vs_Ib.log 2>&1
```

> 虽然之前已有 `domhits_Ib.tbl`（搜 UniRef90 的结果），但为了打分环境完全一致（数据库大小影响 E-value），对同一组 16,211 条序列分别用两个 HMM 重新跑了一次。

### 5. 分数比较与过滤

对每条序列取 **full-sequence score**（domtblout 第 8 列），比较两个模型的打分：

| 统计项 | 数量 |
|--------|------|
| Type Iβ 原始序列 | 16,211 |
| 命中 KDOPS 模型 | 12,074 |
| 命中 DAH7PS Ib 模型 | 16,211 |
| KDOPS score ≥ DAH7PS score | **480** |
| 清洗后剩余 | **15,731** |

过滤脚本 `filter_kdops.py` 完整代码如下：

```python
#!/usr/bin/env python3
"""
Compare hmmsearch scores: KDOPS HMM vs DAH7PS Ib HMM on Type Iβ sequences.
Remove sequences where KDOPS score >= DAH7PS Ib score (or only hit KDOPS).
"""
import sys

def parse_domtbl(path):
    """Parse domtblout, return dict: seqid -> best full-sequence score (column 7, 0-indexed)."""
    best = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) < 22:
                continue
            seqid = cols[0]
            score = float(cols[7])  # full sequence score
            if seqid not in best or score > best[seqid]:
                best[seqid] = score
    return best

# Parse both result tables
kdops_scores = parse_domtbl('domhits_Ib_vs_kdops.tbl')
dah7ps_scores = parse_domtbl('domhits_Ib_vs_dah7ps.tbl')

# Find sequences to remove: KDOPS score >= DAH7PS score, or only in KDOPS
remove_ids = set()
for seqid, kscore in kdops_scores.items():
    dscore = dah7ps_scores.get(seqid, -999999)
    if kscore >= dscore:
        remove_ids.add(seqid)

# Report
all_seqs = set(dah7ps_scores.keys()) | set(kdops_scores.keys())
only_kdops = set(kdops_scores.keys()) - set(dah7ps_scores.keys())
both = set(kdops_scores.keys()) & set(dah7ps_scores.keys())
kdops_wins = {s for s in both if kdops_scores[s] >= dah7ps_scores[s]}

print(f"=== KDOPS Negative Selection Report ===")
print(f"Total unique sequences with hits:  {len(all_seqs)}")
print(f"Sequences with KDOPS hits:         {len(kdops_scores)}")
print(f"Sequences with DAH7PS Ib hits:     {len(dah7ps_scores)}")
print(f"Only hit KDOPS (no DAH7PS hit):    {len(only_kdops)}")
print(f"Hit both, KDOPS score >= DAH7PS:   {len(kdops_wins)}")
print(f"Total to remove:                   {len(remove_ids)}")
print(f"Remaining after filtering:         {len(all_seqs) - len(remove_ids)}")

# Write IDs to remove
with open('kdops_contaminants.txt', 'w') as f:
    for seqid in sorted(remove_ids):
        f.write(seqid + '\n')

print(f"\nContaminant IDs written to: kdops_contaminants.txt")

# Also write the list of IDs to keep (from raw_full_Ib.fasta)
# We need all IDs from the original FASTA, then subtract contaminants
print(f"\nNow use seqkit to filter:")
print(f"  seqkit grep -v -f kdops_contaminants.txt raw_full_Ib.fasta > raw_full_Ib_clean.fasta")
```

### 6. seqkit 输出清洁序列

```bash
seqkit grep -v -f kdops_contaminants.txt raw_full_Ib.fasta > raw_full_Ib_clean.fasta
```

## 产出文件

| 文件 | 说明 |
|------|------|
| `raw_full_Ib_clean.fasta` | ✅ **最终产出**：去除 KDOPS 后的 Type Iβ 序列（15,731 条） |
| `kdops_contaminants.txt` | 被剔除的 480 条 KDOPS 序列 ID |
| `kdo8ps.hmm` | KDOPS HMM profile |
| `kdo8ps_aligned.afa` | KDOPS 多序列比对 |
| `kdo8ps_nr90.fasta` | CD-HIT 去冗余后的 KDOPS 序列 |
| `filter_kdops.py` | 分数比较脚本（代码见上方） |
| `domhits_Ib_vs_kdops.tbl` | KDOPS HMM vs Type Iβ hmmsearch 结果 |
| `domhits_Ib_vs_dah7ps.tbl` | DAH7PS Ib HMM vs Type Iβ hmmsearch 结果 |
