# 🧪 DAH7PS ASR 项目实验记录 (Experiment Log)

> 本文件记录所有操作的精确命令、参数、输出文件和结果，确保完全可重现。
> 工作目录：`/home/tynan/0218/`
> 虚拟环境：`conda activate dah7ps_evo`
> 系统：Ubuntu Linux, 28 核 CPU

---

## 2026-02-18 项目启动日

---

### 09:35 — 步骤 0：环境确认

环境已提前搭建完成。确认状态：

```bash
(base) tynan@OUC-Desktop:~/0218$ conda env list
# conda environments:
base                 *   /home/tynan/miniforge3
dah7ps_evo               /home/tynan/miniforge3/envs/dah7ps_evo
```

`dah7ps_evo` 环境内已安装：`hmmer`, `mafft`, `iqtree`, `trimal`, `cd-hit`, `seqkit`, `matplotlib`

---

### 09:35 — 步骤 1.1：HMM 模型构建（已提前完成）

**精确命令：**

```bash
conda activate dah7ps_evo

# 种子序列 L-INS-i 高精度比对
mafft --maxiterate 1000 --localpair seeds_Ia.fasta > aligned_seeds_Ia.fasta
mafft --maxiterate 1000 --localpair seeds_Ib.fasta > aligned_seeds_Ib.fasta
mafft --maxiterate 1000 --localpair seeds_II.fasta > aligned_seeds_II.fasta

# 构建 HMM 模型
hmmbuild model_Ia.hmm aligned_seeds_Ia.fasta
hmmbuild model_Ib.hmm aligned_seeds_Ib.fasta
hmmbuild model_II.hmm aligned_seeds_II.fasta
```

**输出文件：**

| 文件 | 大小 | 说明 |
|------|------|------|
| `aligned_seeds_Ia.fasta` | 1.5 KB | Iα 种子比对 |
| `aligned_seeds_Ib.fasta` | 1.8 KB | Iβ 种子比对 |
| `aligned_seeds_II.fasta` | 626 B | II 种子比对 |
| `model_Ia.hmm` | 163 KB | Iα HMM 模型 |
| `model_Ib.hmm` | 176 KB | Iβ HMM 模型 |
| `model_II.hmm` | 216 KB | II HMM 模型 |

---

### 09:51 — 步骤 1.2：hmmsearch 扫描 UniRef90

**精确命令：**

```bash
# 通过 nohup 后台串行执行 Ia → Ib → II
nohup bash -c 'eval "$(conda shell.bash hook)" && conda activate dah7ps_evo && \
  hmmsearch --cpu 20 --domtblout domhits_Ia.tbl model_Ia.hmm uniref90.fasta.gz > hmmsearch_Ia.log 2>&1 && \
  hmmsearch --cpu 20 --domtblout domhits_Ib.tbl model_Ib.hmm uniref90.fasta.gz > hmmsearch_Ib.log 2>&1 && \
  hmmsearch --cpu 20 --domtblout domhits_II.tbl model_II.hmm uniref90.fasta.gz > hmmsearch_II.log 2>&1 && \
  echo "ALL_DONE"' > nohup_hmmsearch.log 2>&1 &
```

**⏱ 耗时：** ~1h13min（09:51 → ~11:04）

**输出文件：**

| 文件 | 大小 | 行数 |
|------|------|------|
| `domhits_Ia.tbl` | 3.0 MB | 10,553 |
| `domhits_Ib.tbl` | 7.8 MB | 27,053 |
| `domhits_II.tbl` | 3.6 MB | 12,528 |
| `hmmsearch_Ia.log` | 24.9 MB | — |
| `hmmsearch_Ib.log` | 53.8 MB | — |
| `hmmsearch_II.log` | 28.5 MB | — |

---

### 11:04 — 步骤 1.2 续：序列 ID 提取 + SeqKit 全长提取

**精确命令：**

```bash
conda activate dah7ps_evo

# 提取高置信度序列 ID（domain conditional E-value < 1e-10，即 $12 列）
awk '!/^#/ && $12 < 1e-10 {print $1}' domhits_Ia.tbl | sort | uniq > ids_Ia.txt
awk '!/^#/ && $12 < 1e-10 {print $1}' domhits_Ib.tbl | sort | uniq > ids_Ib.txt
awk '!/^#/ && $12 < 1e-10 {print $1}' domhits_II.tbl | sort | uniq > ids_II.txt

# SeqKit 从 UniRef90 提取全长序列
seqkit grep -f ids_Ia.txt uniref90.fasta.gz > raw_full_Ia.fasta
seqkit grep -f ids_Ib.txt uniref90.fasta.gz > raw_full_Ib.fasta
seqkit grep -f ids_II.txt uniref90.fasta.gz > raw_full_II.fasta
```

**结果：**

| 文件 | 序列 ID 数 | 全长文件大小 |
|------|-----------|-------------|
| `ids_Ia.txt` → `raw_full_Ia.fasta` | 10,102 | 5.0 MB |
| `ids_Ib.txt` → `raw_full_Ib.fasta` | 16,211 | 7.9 MB |
| `ids_II.txt` → `raw_full_II.fasta` | 9,862 | 5.6 MB |

---

### 下午 — 漏洞发现 & KDOPS 反向过滤（步骤 1.3）

**漏洞 1（Jiao 2020）：** Type Iα 和 II 的变构元件是 TIM 桶内部插片，物理切除会断裂多肽链。
**漏洞 2（Yokoyama 2025）：** KDOPS 是 Type Iβ 的姊妹枝，hmmsearch 必然捞出假阳性。

**KDOPS 反向过滤精确命令：**

```bash
conda activate dah7ps_evo

# 1. KDOPS 参考序列准备（kdo8ps_uniprot.fasta 从 UniProt 手动下载）
cd-hit -i kdo8ps_uniprot.fasta -o kdo8ps_nr90.fasta -c 0.9 -n 5 -M 4000 -T 20

# 2. 构建 KDOPS HMM 诱饵模型
mafft --auto kdo8ps_nr90.fasta > kdo8ps_aligned.afa
hmmbuild kdo8ps.hmm kdo8ps_aligned.afa

# 3. 双向竞争打分
hmmsearch --cpu 20 --domtblout domhits_Ib_vs_kdops.tbl kdo8ps.hmm raw_full_Ib.fasta > hmmsearch_kdops_vs_Ib.log 2>&1
hmmsearch --cpu 20 --domtblout domhits_Ib_vs_dah7ps.tbl model_Ib.hmm raw_full_Ib.fasta > hmmsearch_dah7ps_vs_Ib.log 2>&1

# 4. Python 过滤脚本：比较双向打分，剔除 KDOPS 打分更高者
python filter_kdops.py
# → 输出 kdops_contaminants.txt（被识别为 KDOPS 的污染 ID）

# 5. 反向提取干净序列
seqkit grep -v -f kdops_contaminants.txt raw_full_Ib.fasta > raw_full_Ib_clean.fasta
```

**输出文件：**

| 文件 | 大小 | 说明 |
|------|------|------|
| `kdo8ps_uniprot.fasta` | 143 KB | KDOPS UniProt 参考序列 |
| `kdo8ps_nr90.fasta` | 51 KB | 90% 去冗余后 |
| `kdo8ps_aligned.afa` | 64 KB | MAFFT 比对 |
| `kdo8ps.hmm` | 131 KB | KDOPS HMM 模型 |
| `domhits_Ib_vs_kdops.tbl` | 3.6 MB | Iβ vs KDOPS 打分 |
| `domhits_Ib_vs_dah7ps.tbl` | 4.7 MB | Iβ vs DAH7PS 打分 |
| `kdops_contaminants.txt` | 10 KB | 污染 ID 列表 |
| `raw_full_Ib_clean.fasta` | 7.7 MB | 净化后 Iβ 全长序列 |

---

### 下午 — 长度直方图分析

**精确命令：**

```bash
python plot_length_hist.py
# → 输出 length_distribution.png
```

**直方图数据驱动结论：**

| Type | 主峰 | 过滤区间 | 科学依据 |
|------|------|----------|----------|
| Iα | ~350 aa | 300–450 aa | 极度保守 |
| Iβ (clean) | ~370 aa | 300–480 aa | 尾巴 400–450 为 ACT/CM 附件 |
| II | ~450 aa | 380–650 aa | 肩峰 480–600 为植物叶绿体转运肽 |

---

### 16:30 — 步骤 QC1：长度过滤 + CD-HIT 80% 去冗余

**精确命令：**

```bash
conda activate dah7ps_evo

# 长度过滤
seqkit seq -m 300 -M 450 raw_full_Ia.fasta > qc_len_Ia.fasta
seqkit seq -m 300 -M 480 raw_full_Ib_clean.fasta > qc_len_Ib.fasta
seqkit seq -m 380 -M 650 raw_full_II.fasta > qc_len_II.fasta

# CD-HIT 80% 去冗余
cd-hit -i qc_len_Ia.fasta -o caseA_full_Ia.fasta -c 0.8 -n 5 -M 4000 -T 20
cd-hit -i qc_len_Ib.fasta -o caseA_full_Ib.fasta -c 0.8 -n 5 -M 4000 -T 20
cd-hit -i qc_len_II.fasta -o caseA_full_II.fasta -c 0.8 -n 5 -M 4000 -T 20
```

**seqkit stats 结果：**

```
=== 长度过滤后 ===
file             num_seqs  min_len  avg_len  max_len
qc_len_Ia.fasta     9,071      300    361.9      450
qc_len_Ib.fasta    13,879      300    356.8      479
qc_len_II.fasta     8,446      380    454.7      650

=== CD-HIT 80% 后 ===
file                 num_seqs  min_len  avg_len  max_len
caseA_full_Ia.fasta     3,473      300    363.3      450
caseA_full_Ib.fasta     5,728      300    356.8      479
caseA_full_II.fasta     3,064      380    451.1      650
```

**输出文件：**

| 文件 | 序列数 | 大小 |
|------|--------|------|
| `caseA_full_Ia.fasta` | 3,473 | 1.7 MB |
| `caseA_full_Ib.fasta` | 5,728 | 2.7 MB |
| `caseA_full_II.fasta` | 3,064 | 1.8 MB |

---

### 17:08 — V3.0 范式跃迁

识别新漏洞并升级：
- trimAl -gappyout 丢失变构铰链 Gap → **ClipKIT**
- ML 单一祖先过度稳定化 → **AltAll 系综**
- AlphaFold 静态快照 → **MD ≥500ns×3**
- 纯 DCA → **ICDC (DCA × DCCM)**

**ClipKIT 安装命令：**

```bash
conda activate dah7ps_evo
pip install clipkit
# → 安装 clipkit-2.10.2, biopython-1.86, numpy-1.26.4
```

---

## 2026-02-19 MSA 执行日

---

### 04:15 — V3.1 MSA 策略专家评审

与专业人士讨论确定：
- **情况 A**：Seed & Add（CD-HIT 60% 种子 → L-INS-i → MAFFT --add）
- **情况 B**：PROMALS3D 混合骨架（PDB + 进化踏脚石 → MAFFT --add）
- **修剪**：ClipKIT `kpi-smart-gap`

---

### 04:20 — Phase 1 Step 1: CD-HIT 60% 提取种子

**精确命令：**

```bash
conda activate dah7ps_evo

cd-hit -i caseA_full_Ia.fasta -o seeds60_Ia.fasta -c 0.6 -n 4 -M 4000 -T 20
cd-hit -i caseA_full_Ib.fasta -o seeds60_Ib.fasta -c 0.6 -n 4 -M 4000 -T 20
cd-hit -i caseA_full_II.fasta -o seeds60_II.fasta -c 0.6 -n 4 -M 4000 -T 20
```

**seqkit stats 结果：**

```
file              num_seqs  min_len  avg_len  max_len
seeds60_Ia.fasta       537      300    370.2      450
seeds60_Ib.fasta       963      300    363.3      479
seeds60_II.fasta       652      380    458.0      650
```

**输出文件：**

| 文件 | 种子数 | 大小 |
|------|--------|------|
| `seeds60_Ia.fasta` | 537 | 267 KB |
| `seeds60_Ib.fasta` | 963 | 472 KB |
| `seeds60_II.fasta` | 652 | 380 KB |

---

### 05:40 — Phase 1 Step 2: L-INS-i 种子骨架比对

**精确命令：**

```bash
# nohup 后台串行执行（注意：此次使用 --thread -1，后续应改为 --thread 20）
nohup bash -c 'eval "$(conda shell.bash hook)" && conda activate dah7ps_evo && \
  mafft --localpair --maxiterate 1000 --thread -1 seeds60_Ia.fasta > aligned_seeds60_Ia.afa 2>mafft_seeds60_Ia.log && \
  mafft --localpair --maxiterate 1000 --thread -1 seeds60_Ib.fasta > aligned_seeds60_Ib.afa 2>mafft_seeds60_Ib.log && \
  mafft --localpair --maxiterate 1000 --thread -1 seeds60_II.fasta > aligned_seeds60_II.afa 2>mafft_seeds60_II.log && \
  echo "ALL_LINSI_DONE"' > nohup_linsi.log 2>&1 &
```

> ⚠️ **注意：** `--thread -1` 实际只使用了 8 核（`-C 8`），后续命令已统一改为 `--thread 20`

**已完成结果（Iα）：**

```
[Thu Feb 19 05:40:16 CST 2026] Starting Ia L-INS-i...
[Thu Feb 19 05:43:02 CST 2026] Ia DONE        # ⏱ ~3 分钟
```

| 文件 | 序列数 | 比对长度 | 大小 |
|------|--------|---------|------|
| `aligned_seeds60_Ia.afa` | 537 | 1,920 aa | 1.1 MB |

**完成结果（全部）：**

```
[Thu Feb 19 05:40:16 CST 2026] Starting Ia L-INS-i...
[Thu Feb 19 05:43:02 CST 2026] Ia DONE        # ⏱ ~3 分钟
[Thu Feb 19 05:43:02 CST 2026] Starting Ib L-INS-i...
[Thu Feb 19 06:18:59 CST 2026] Ib DONE        # ⏱ ~36 分钟
[Thu Feb 19 06:18:59 CST 2026] Starting II L-INS-i...
[Thu Feb 19 06:33:12 CST 2026] II DONE        # ⏱ ~14 分钟
ALL_LINSI_DONE
```

| 文件 | 序列数 | 比对长度 | 大小 |
|------|--------|---------|------|
| `aligned_seeds60_Ia.afa` | 537 | 1,920 aa | 1.1 MB |
| `aligned_seeds60_Ib.afa` | 963 | 3,239 aa | 3.2 MB |
| `aligned_seeds60_II.afa` | 652 | 4,493 aa | 3.0 MB |

---

### 09:00 — Phase 1 Step 3: 增量映射（MAFFT --add）

**精确命令：**

```bash
conda activate dah7ps_evo

# 提取种子 ID
seqkit seq -n seeds60_Ia.fasta | awk '{print $1}' > seed_ids_Ia.txt
seqkit seq -n seeds60_Ib.fasta | awk '{print $1}' > seed_ids_Ib.txt
seqkit seq -n seeds60_II.fasta | awk '{print $1}' > seed_ids_II.txt

# 提取非种子的剩余序列
seqkit grep -v -f seed_ids_Ia.txt caseA_full_Ia.fasta > remaining_Ia.fasta
seqkit grep -v -f seed_ids_Ib.txt caseA_full_Ib.fasta > remaining_Ib.fasta
seqkit grep -v -f seed_ids_II.txt caseA_full_II.fasta > remaining_II.fasta
```

**剩余序列数：** Ia=2,936 / Ib=4,765 / II=2,412

```bash
# 增量映射到骨架（--thread 20）
mafft --add remaining_Ia.fasta --thread 20 aligned_seeds60_Ia.afa > msa_full_Ia.afa 2>mafft_add_Ia.log
mafft --add remaining_Ib.fasta --thread 20 aligned_seeds60_Ib.afa > msa_full_Ib.afa 2>mafft_add_Ib.log
mafft --add remaining_II.fasta --thread 20 aligned_seeds60_II.afa > msa_full_II.afa 2>mafft_add_II.log
```

**时间戳：**

```
[Thu Feb 19 09:00:17 CST 2026] Starting Ia --add mapping...
[Thu Feb 19 09:00:30 CST 2026] Ia DONE        # ⏱ ~13 秒
[Thu Feb 19 09:00:30 CST 2026] Starting Ib --add mapping...
[Thu Feb 19 09:01:24 CST 2026] Ib DONE        # ⏱ ~54 秒
[Thu Feb 19 09:01:24 CST 2026] Starting II --add mapping...
[Thu Feb 19 09:01:47 CST 2026] II DONE        # ⏱ ~23 秒
ALL_ADD_DONE                                   # 总计 < 2 分钟
```

**✅ Phase 1 最终产出（情况 A 全长 MSA）：**

| 文件 | 序列数 | 比对长度 | 大小 | 用途 |
|------|--------|---------|------|------|
| `msa_full_Ia.afa` | 3,473 | 3,544 aa | 13 MB | DCA/EVcouplings |
| `msa_full_Ib.afa` | 5,728 | 5,185 aa | 30 MB | DCA/EVcouplings |
| `msa_full_II.afa` | 3,064 | 7,689 aa | 24 MB | DCA/EVcouplings |

---

### 待记录

- [ ] Phase 2: 混合骨架 + PROMALS3D + 万级映射
- [ ] Phase 3: ClipKIT kpi-smart-gap
- [ ] QC2: Jalview 催化残基核验
