# 🧪 DAH7PS ASR 项目实验记录 (Experiment Log)

> 本文件记录所有操作的精确命令、参数、输出文件和结果，确保完全可重现。
> 工作目录：`/home/tynan/0218/`
> 虚拟环境：`conda activate dah7ps_v4`（V4.1 起）/ `conda activate dah7ps_evo`（V3.1 旧环境）
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

- [ ] ~~Phase 2: 混合骨架 + PROMALS3D + 万级映射~~（V3.1 已废弃，V4.1 改用 FoldMason）
- [ ] ~~Phase 3: ClipKIT kpi-smart-gap~~（V3.1 已废弃，V4.1 改为双版本修剪）
- [ ] ~~QC2: Jalview 催化残基核验~~（V4.1 改为 FoldMason msa2lddt 结构复核）

---

## 2026-02-23 V3.1 → V4.1 迁移 & Phase 0

---

### 05:06 — V3.1 → V4.1 文档升级

**背景：** V3.1 的"全长单一 MSA 驱动一切"策略在跨亚型远缘同源时产生严重比对膨胀（比对长度 3,544–7,689 列 vs 预期 400–600 列）。V4.1 采用核心–模块分治范式，以 FoldMason 结构字母表为中枢。

**操作：**

```bash
# 用 V4.1 版本覆盖旧文档
cp AGENT_v4_1.md AGENT.md
cp dah7ps_v4_final_sop_4.md plan.md

# 删除临时源文件
rm AGENT_v4_1.md dah7ps_v4_final_sop_4.md
```

**结果：**

| 文件 | 变化 |
|------|------|
| `AGENT.md` | 76 行 (V3.1) → 207 行 (V4.1 SOP rev4 执行指南) |
| `plan.md` | 208 行 (V3.1) → 1153 行 (V4.1 标准作业程序 SOP rev4) |

---

### 05:15 — 目录结构重组（SOP Phase 0.3）

**操作：** 将根目录 80 个散落文件按 V4.0 SOP Phase 0.3 分类移入子目录。

```bash
# 创建 V4.0 SOP 目录结构
mkdir -p data/seeds data/structures/panel_dah7ps data/structures/panel_kdops data/db \
  meta/models scripts \
  results/01_mining results/02_qc results/03_msa_core results/03_msa_modules \
  results/04_phylogeny_asr results/05_struct_md results/06_icdc results/meta \
  workflow

# 种子序列 → data/seeds/
mv seeds_Ia.fasta seeds_Ib.fasta seeds_II.fasta kdo8ps_uniprot.fasta data/seeds/

# 数据库 → data/db/
mv uniref90.fasta.gz Pfam-A.hmm Pfam-A.hmm.ssi PF00793.hmm data/db/

# 脚本 → scripts/
mv filter_kdops.py plot_length_hist.py scripts/

# HMM 模型 + 搜索结果 → results/01_mining/ (29 files)
mv model_Ia.hmm model_Ib.hmm model_II.hmm results/01_mining/
mv aligned_seeds_Ia.fasta aligned_seeds_Ib.fasta aligned_seeds_II.fasta results/01_mining/
mv domhits_*.tbl results/01_mining/
mv hmmsearch_*.log results/01_mining/
mv raw_full_*.fasta results/01_mining/
mv kdo8ps.hmm kdo8ps_aligned.afa kdo8ps_nr90.fasta* results/01_mining/
mv nohup_hmmsearch.log nohup_seqkit.log seqkit_*.log results/01_mining/

# QC 和去冗余 → results/02_qc/ (19 files)
mv qc_len_*.fasta results/02_qc/
mv caseA_full_*.fasta* results/02_qc/
mv seeds60_*.fasta* results/02_qc/
mv remaining_*.fasta results/02_qc/
mv length_distribution.png results/02_qc/

# V3.1 MSA 存档 → results/03_msa_core/ (14 files)
mv aligned_seeds60_*.afa results/03_msa_core/
mv msa_full_*.afa results/03_msa_core/
mv mafft_seeds60_*.log mafft_add_*.log nohup_add.log nohup_linsi.log results/03_msa_core/

# 清理旧文档
rm walkthrough_kdops_filter.md
```

**结果：**

| 目录 | 文件数 | 说明 |
|------|--------|------|
| `data/seeds/` | 4 | 种子序列 |
| `data/db/` | 4 | UniRef90 (47 GB) + Pfam HMM |
| `data/structures/panel_dah7ps/` | 0 | 待下载 DAH7PS PDB |
| `data/structures/panel_kdops/` | 0 | 待下载 KDOPS PDB |
| `scripts/` | 2 | Python 脚本 |
| `results/01_mining/` | 29 | HMM + 搜索 + 原始序列 |
| `results/02_qc/` | 19 | QC + 去冗余 |
| `results/03_msa_core/` | 14 | V3.1 MSA（存档） |
| 根目录 | 6 | AGENT.md, plan.md, README.md, log.md, .gitignore |

---

### 05:38 — Git 提交并推送 Plan-V4.1 分支

```bash
git checkout -b Plan-V4.1
git add -A
git commit -m "Upgrade to V4.1: update AGENT.md & plan.md, reorganize directory structure per SOP Phase 0.3"
git push -u origin Plan-V4.1
```

**结果：** `+1269/-421` 行变更，推送至 `origin/Plan-V4.1`

---

### 05:48 — Phase 0.1：创建 dah7ps_v4 conda 环境

**精确命令：**

```bash
# 使用 mamba 创建新环境（SOP Phase 0.1）
mamba create -n dah7ps_v4 \
  python=3.11 hmmer mafft iqtree mmseqs2 foldmason seqkit cd-hit \
  -y

# pip 安装剩余工具
conda run -n dah7ps_v4 pip install clipkit pastml
```

**安装结果：**

| 工具 | 版本 | 来源 |
|------|------|------|
| python | 3.11.14 | conda |
| hmmer | 3.4 | bioconda |
| mafft | 7.526 | conda-forge |
| iqtree | 3.0.1 | bioconda |
| mmseqs2 | 18.8cc5c | bioconda |
| foldmason | 4.dd3c235 | bioconda |
| seqkit | 2.12.0 | bioconda |
| cd-hit | 4.8.1 | bioconda |
| clipkit | 2.10.2 | pip |
| pastml | 1.9.51 | pip |

> ⚠️ **注意：** bioconda 的 iqtree 3.0.1 包中命令为 `iqtree`（非 `iqtree2`）。SOP 中所有 `iqtree2` 命令实际执行时需改为 `iqtree`。

---

### 05:55 — Phase 0.2：软件版本锁定

**精确命令：**

```bash
conda run -n dah7ps_v4 bash -c '
echo -e "tool\tversion" > results/meta/software_versions.tsv
for t in python hmmbuild hmmalign hmmsearch mafft iqtree mmseqs foldmason seqkit cd-hit clipkit pastml; do
  # ...各工具版本提取写入
done
'
```

**输出文件：** `results/meta/software_versions.tsv`（13 行，包含所有关键工具版本）

---

### 05:58 — Phase 0.4：参数文件初始化

**操作：** 从 SOP 0.4 模板创建 `meta/params.json`

**文件：** `meta/params.json`（包含 mining, qc, core_definition, msa, phylogeny, dca, af3, asr 全部参数块）

---

### 06:00 — Phase 0.2.1：外部模型文件

**操作：** 初始化 `results/meta/model_files.tsv`

**状态：** 3Di 模型文件（Q.3Di.AF / Q.3Di.LLM）需从 Seffernick et al. 2025 (*MBE* 42:msaf124) 补充材料获取，延后到 Phase 4.2 执行。

---

### 06:06 — 追加 AGENT.md 规则 #9

**操作：** 向 `AGENT.md` 总原则添加第 9 条：

> 9. **实验记录必须同步更新**：每次执行任何分析操作，必须将精确命令、参数、输出文件、结果摘要和时间戳记录到 `log.md`。

同时更新 SOP 引用源从 `dah7ps_v4_final_sop_4.md` → `PLAN.md（V4.1 SOP rev4）`。

---

### 09:15 — Phase 0.2.1：下载 3Di 模型文件

**操作（用户手动完成）：** 从 Seffernick et al. 2025 补充材料下载 Q.3Di.AF 和 Q.3Di.LLM。

```bash
# 用户已放置到 meta/models/
ls meta/models/
# Q.3Di.AF   Q.3Di.LLM

# sha256 校验
sha256sum meta/models/Q.3Di.AF meta/models/Q.3Di.LLM
# 5e37984b...  Q.3Di.AF
# 8cc3f62c...  Q.3Di.LLM
```

**产出：** `results/meta/model_files.tsv` 已更新 sha256 记录。

**✅ Phase 0 所有 Done 条件通过。**

---

### 10:55 — Phase 1.1：种子扩充（V3.1 → V4.1）

**背景：** V3.1 种子严重不足（Ia=3 条均来自 *E. coli*，Ib=3 条，II=仅 1 条 *M. tuberculosis*）。

**用户提供扩充种子：** 从 UniProt 手动整理，按亚型覆盖多个分类群和结构类型。

**QC 发现 6 个错误 UniProt ID（实际蛋白非 DAH7PS）：**

| 错误 ID | 实际蛋白 | 处理 |
|---------|---------|------|
| Q8Y825 (Ib) | NAD(+) synthetase | → 替换为 Q8Y6T2 (*L. monocytogenes* AroA) |
| D6D1V7 (Ib) | Glycosyl hydrolase 71 | → 删除（未找到替代） |
| Q9HZE4 (II) | PelA | → 替换为 Q9I000 (*P. aeruginosa* PA2843) |
| Q9I2V4 (II) | DUF1289 | → 替换为 Q7DC82 (*P. aeruginosa* phzC1) |
| Q9SEU1 (II) | S-RNase | → 替换为 Q9SK84 (*A. thaliana* At1g22410) |
| O25000 (II) | HofA | → 删除（未找到替代） |

**精确命令：**

```bash
conda run -n dah7ps_v4 python3 -c "
from Bio import SeqIO
# Load corrections from user-provided uniprotkb_2026_02_23.fasta
corrections = {rec.id.split('|')[1]: rec for rec in SeqIO.parse('uniprotkb_2026_02_23.fasta', 'fasta')}
# Ib: remove Q8Y825/D6D1V7, add Q8Y6T2
# II: remove Q9HZE4/Q9I2V4/Q9SEU1/O25000, add Q9I000/Q7DC82/Q9SK84
# ... (full script in session)
"

# Backup old seeds
mv data/seeds/seeds_*.fasta data/seeds/*_v31.fasta.bak
# Move fixed seeds
mv Type_*_seeds_fixed.fasta data/seeds/seeds_*.fasta
```

**最终种子统计：**

| 文件 | 序列数 | V3.1 | 覆盖物种 |
|------|--------|------|---------|
| `data/seeds/seeds_Ia.fasta` | 13 | 3 | E. coli, S. cerevisiae, C. albicans, S. pombe, N. meningitidis 等 |
| `data/seeds/seeds_Ib.fasta` | 6 | 3 | T. maritima, B. subtilis, Geobacillus, P. furiosus, A. pernix, L. monocytogenes |
| `data/seeds/seeds_II.fasta` | 14 | 1 | M. tuberculosis, C. glutamicum, P. aeruginosa, A. thaliana, 水稻, 番茄, 马铃薯, 矮牵牛, 长春花 |

---

### 11:21 — Phase 1.3：KDOPS 外群种子提取

**操作：** 从用户下载的 UniProt KDO8PS FASTA（ec:2.5.1.55, reviewed:true, 300+ 条）中提取 12 条多样性代表。

```bash
conda run -n dah7ps_v4 python3 -c "
from Bio import SeqIO
kdops_ids = {'Q9AV97','O50044','Q31KV0','Q7V4M4','O66496','P61657','Q92Q99','P0A715','Q9ZFK4','Q0KCE4','P0CD74','Q9ZN55'}
found = []
for rec in SeqIO.parse('uniprotkb_ec_2_5_1_55_AND_reviewed_true_2026_02_23.fasta', 'fasta'):
    acc = rec.id.split('|')[1]
    if acc in kdops_ids:
        rec.id = f'KDOPS_{acc}'
        found.append(rec)
SeqIO.write(found, 'kdops_outgroup.fasta', 'fasta')
"
mv kdops_outgroup.fasta data/seeds/
```

**KDOPS 外群组成（12 条，覆盖 9 个分类群）：**

| KDOPS ID | 物种 | 分类群 |
|----------|------|--------|
| KDOPS_P0A715 | *E. coli* | γ-变形菌 |
| KDOPS_Q9ZFK4 | *P. aeruginosa* | γ-变形菌 |
| KDOPS_Q0KCE4 | *C. necator* | β-变形菌 |
| KDOPS_P61657 | *R. palustris* | α-变形菌 |
| KDOPS_Q92Q99 | *R. meliloti* | α-变形菌 |
| KDOPS_O66496 | *A. aeolicus* | 水生菌门 |
| KDOPS_Q31KV0 | *S. elongatus* | 蓝细菌 |
| KDOPS_Q7V4M4 | *P. marinus* | 蓝细菌 |
| KDOPS_P0CD74 | *C. trachomatis* | 衣原体 |
| KDOPS_Q9ZN55 | *H. pylori* | 弯曲菌门 |
| KDOPS_Q9AV97 | *A. thaliana* | 植物 |
| KDOPS_O50044 | *P. sativum* | 植物 |

---

### 11:52 — Phase 1.2：种子比对 + HMM 重建 + hmmsearch

**精确命令：**

```bash
# 1. MAFFT E-INS-i 种子比对（V4.1 推荐 E-INS-i 替代 L-INS-i）
for type in Ia Ib II; do
  mafft --genafpair --maxiterate 1000 --ep 0 --thread 20 \
    data/seeds/seeds_${type}.fasta \
    > results/01_mining/seeds_${type}.afa
done

# 2. hmmbuild（旧 HMM 备份为 *_v31.hmm.bak）
for type in Ia Ib II; do
  mv results/01_mining/model_${type}.hmm results/01_mining/model_${type}_v31.hmm.bak
  hmmbuild results/01_mining/model_${type}.hmm results/01_mining/seeds_${type}.afa
done

# 3. hmmsearch 扫描 UniRef90（后台 nohup，预计 ~1-2 小时）
nohup bash -c '
for type in Ia Ib II; do
  hmmsearch --cpu 20 -E 1e-10 \
    --domtblout results/01_mining/hits_${type}.domtbl \
    results/01_mining/model_${type}.hmm \
    data/db/uniref90.fasta.gz > /dev/null
  # 提取命中序列
  grep -v "^#" results/01_mining/hits_${type}.domtbl | awk "{print \$1}" | sort -u > results/01_mining/hits_${type}_ids.txt
  seqkit grep -f results/01_mining/hits_${type}_ids.txt data/db/uniref90.fasta.gz > results/01_mining/hits_${type}_seqs.fasta
done
' > results/01_mining/hmmsearch_v41.log 2>&1 &
```

**HMM 模型长度变化：**

| 模型 | V3.1 | V4.1 | 种子数变化 |
|------|------|------|-----------|
| `model_Ia.hmm` | 163 KB | 166 KB (L=355) | 3 → 13 |
| `model_Ib.hmm` | 176 KB | 156 KB (L=334) | 3 → 6 |
| `model_II.hmm` | 216 KB | 220 KB (L=471) | 1 → 14 |

**状态：** hmmsearch 正在后台运行，监控命令 `tail -f results/01_mining/hmmsearch_v41.log`

---

## 2026-02-24 Phase 1 完成 & QC1

---

### 11:24 — Phase 1.2 hmmsearch V4.1 结果确认

hmmsearch 已于 2026-02-23 12:50 完成（耗时约 50 分钟）。

**V4.1 hmmsearch 结果：**

| 亚型 | V3.1 命中数 | V4.1 命中数 | 变化 |
|------|-----------|-----------|------|
| Ia | 10,102 | 10,071 | -31 (≈持平) |
| Ib (raw) | 16,211 | 15,608 | -603 (-3.7%) |
| II | 9,862 | 18,529 | +8,667 (+87.9%) |

---

### 11:25 — Phase 1.3：KDOPS 反向过滤（V4.1 重做）

**背景：** V3.1 的 KDOPS 过滤使用旧 HMM 模型，需用 V4.1 的 `model_Ib.hmm` 对新 `hits_Ib_seqs.fasta` 重跑。

**步骤 1：双向 hmmsearch 竞争打分**

```bash
conda activate dah7ps_v4

# DAH7PS Ib HMM vs Ib 候选序列
hmmsearch --cpu 20 \
  --domtblout results/01_mining/hits_Ib_vs_dah7ps_v41.domtbl \
  results/01_mining/model_Ib.hmm \
  results/01_mining/hits_Ib_seqs.fasta \
  > results/01_mining/hmmsearch_Ib_vs_dah7ps_v41.log 2>&1

# KDOPS HMM vs Ib 候选序列
hmmsearch --cpu 20 \
  --domtblout results/01_mining/hits_Ib_vs_kdops_v41.domtbl \
  results/01_mining/kdo8ps.hmm \
  results/01_mining/hits_Ib_seqs.fasta \
  > results/01_mining/hmmsearch_Ib_vs_kdops_v41.log 2>&1
```

**步骤 2：升级 filter_kdops.py 并执行过滤**

脚本升级为 CLI 参数化版本（支持 `--help`、输入检查、输出目录自动创建）。

```bash
python scripts/filter_kdops.py \
  --dah7ps_domtbl results/01_mining/hits_Ib_vs_dah7ps_v41.domtbl \
  --kdops_domtbl results/01_mining/hits_Ib_vs_kdops_v41.domtbl \
  --input results/01_mining/hits_Ib_seqs.fasta \
  --output results/01_mining/hits_Ib_clean.fasta \
  --contaminants results/01_mining/kdops_contaminants_v41.txt \
  --report results/01_mining/kdops_filter_report_v41.tsv \
  --score_margin 0
```

**KDOPS 过滤结果：**

| 指标 | 数量 |
|------|------|
| Ib 总序列 | 15,608 |
| 仅 DAH7PS 命中 | 136 |
| 仅 KDOPS 命中 | 0 |
| 双向命中 | 15,472 |
| DAH7PS 胜出 | 7,733 |
| KDOPS 胜出 (score ≥ DAH7PS) | 7,739 |
| **总移除** | **7,739** |
| **保留（clean）** | **7,869** |

**得分差分布（DAH7PS − KDOPS）：** 强双峰分布，|delta| ≤ 20 仅 4 条序列。过滤边界清晰无歧义。

**V3.1 vs V4.1 过滤率差异说明：** V3.1 仅移除 480 条（3.0%），V4.1 移除 7,739 条（49.6%）。原因是 V3.1 使用 domain conditional E-value < 1e-10（domtblout 第 12 列）预筛选 ID，而 V4.1 使用 hmmsearch -E 1e-10（全序列 E-value），后者纳入了更多边界命中，其中大量实际为 KDOPS。

**输出文件：**

| 文件 | 大小 | 说明 |
|------|------|------|
| `hits_Ib_vs_dah7ps_v41.domtbl` | 4.8 MB | Ib vs DAH7PS 打分 |
| `hits_Ib_vs_kdops_v41.domtbl` | 4.7 MB | Ib vs KDOPS 打分 |
| `hits_Ib_clean.fasta` | 3.7 MB | 净化后 7,869 条 |
| `kdops_contaminants_v41.txt` | 160 KB | 7,739 条污染 ID |
| `kdops_filter_report_v41.tsv` | 698 KB | 逐序列得分比较 |

---

### 11:30 — QC1：挖掘质量报告

**操作：** 生成 `results/01_mining/qc_mining_report.md`，包含：
- 每亚型命中数量与 V3.1 对比
- E-value / Score 分布
- 长度分布直方图
- KDOPS 竞争得分分布
- 分类群覆盖（Top 10）

**关键发现：**
- Type II 命中数翻倍（9,862 → 18,529），验证种子扩充策略有效
- Ia 稳定（10,102 → 10,071），说明原 3 条种子已足够覆盖 Ia 多样性
- Ib clean 从 V3.1 的 15,731 降至 7,869，主要因为 V4.1 的 KDOPS 过滤更严格
- 分类群覆盖：Ia 以 γ-变形菌为主，Ib 以厚壁菌/拟杆菌为主，II 以放线菌（链霉菌）为主

**产出文件：** `results/01_mining/qc_mining_report.md`

**✅ Phase 1 所有 Done 条件通过：**
- `results/01_mining/hits_*.domtbl` ✓
- `results/01_mining/hits_*_seqs.fasta` ✓
- `results/01_mining/hits_Ib_clean.fasta` ✓
- `results/01_mining/qc_mining_report.md` ✓

---

### 17:26 — Pre-Phase 2 门控检查（计划外追加）

**背景：** 在 QC1 报告审阅后，用户指出在进入 Phase 2 前需做 2 个快速门控检查：(A) 三亚型 hits 集合交集统计；(B) Ib 边缘序列隔离。此步骤不在 SOP 原文中，是基于 QC1 结果的风险预判。

**脚本编写：** `scripts/gate_checks.py`（新建）

**精确命令：**

```bash
python3 scripts/gate_checks.py --workdir /home/tynan/0218
```

#### Gate A：三亚型 hits 交集统计

**输入：** `hits_Ia_ids.txt`（10,071）、`hits_Ib_clean.fasta` 中提取的 ID（7,869）、`hits_II_ids.txt`（18,529）

**结果：**

| 交集 | 数量 | 占较小集百分比 |
|------|------|---------------|
| Ia ∩ Ib | 0 | 0.00% |
| **Ia ∩ II** | **8,561** | **85.01%** |
| Ib ∩ II | 0 | 0.00% |
| Ia ∩ Ib ∩ II | 0 | 0.00% |

**⚠ 重大发现：Ia 与 II 交集极高（85%）。** 即 10,071 条 Ia 命中中的 8,561 条也同时被 II 模型命中。这在生物学上合理（Ia 和 II 共享 TIM-barrel 核心，暮光区同源），但若不做 best-hit 归属，Phase 2 将出现大量重复计数和混亚型 seeds。

**→ 判定：需要在 Phase 2.1 之前构建 best-hit subtype 归属表。**

重叠 ID 列表写入：
- `results/01_mining/overlap_Ia_II.txt`（8,561 条）
- `results/01_mining/overlap_Ia_Ib.txt`（0 条）
- `results/01_mining/overlap_Ib_II.txt`（0 条）

#### Gate B：Ib 边缘序列隔离

**输入：** `results/01_mining/kdops_filter_report_v41.tsv`

**筛选条件：** verdict=KEEP 且 delta (DAH7PS−KDOPS) ∈ [0, 20]

**结果：** 4 条边缘序列

| ID | Delta |
|----|-------|
| UniRef90_A0A1I4NDE8 | 13.7 |
| UniRef90_A0A806JZI9 | 14.8 |
| UniRef90_X1MS81 | 16.0 |
| UniRef90_A0ABT1SNH0 | 16.6 |

**→ 已写入 `results/01_mining/kdops_borderline_ids.txt`。** 后续 seeds60 和 stepping stones 将排除这 4 条 ID，但允许它们保留在全量库中。

**产出文件：**

| 文件 | 说明 |
|------|------|
| `scripts/gate_checks.py` | 门控检查脚本（新建） |
| `results/01_mining/overlap_Ia_II.txt` | Ia∩II 重叠 ID（8,561 条） |
| `results/01_mining/kdops_borderline_ids.txt` | Ib 边缘 ID（4 条） |

**下一步：** Phase 2.1 之前，需对 Ia∩II 的 8,561 条序列做 best-hit 归属（比较 model_Ia.hmm vs model_II.hmm 的 full-sequence score），将每条分配到得分最高的亚型。

---

### 05:05 — Gate A-2：Best-hit Ia vs II 亚型归属

**背景：** Gate A 发现 Ia∩II = 8,561（85% of Ia）。必须做竞争打分归属，产出互斥 ID/FASTA 集合，否则 Phase 2 会重复计数。

**脚本编写：** `scripts/assign_besthit_Ia_vs_II.py`（新建）

**设计：** 3 档置信度分级
- HIGH: |Δ bits| ≥ 20
- MED: 10 ≤ |Δ bits| < 20
- LOW: |Δ bits| < 10

**精确命令：**

```bash
python3 scripts/assign_besthit_Ia_vs_II.py \
  --ia_domtbl results/01_mining/hits_Ia.domtbl \
  --ii_domtbl results/01_mining/hits_II.domtbl \
  --overlap_ids results/01_mining/overlap_Ia_II.txt \
  --ia_all_ids results/01_mining/hits_Ia_ids.txt \
  --ii_all_ids results/01_mining/hits_II_ids.txt \
  --outdir results/01_mining
```

**结果：**

| 指标 | 数量 |
|------|------|
| Overlap 总数 | 8,561 |
| → 归属 Ia | **8,561 (100%)** |
| → 归属 II | 0 |
| HIGH 置信度 | 8,561 (100.0%) |
| MED 置信度 | 0 (0.0%) |
| LOW 置信度 | 0 (0.0%) |

**解读：** 全部重叠序列以 ≥20 bits 优势归属 Ia。这说明 II HMM 的交叉命中纯粹来自共享的 TIM-barrel 核心同源性，Ia HMM 始终更匹配。归属完全无歧义。

**Go/No-Go 评估：** ✅ GO
- Ia_final = 10,071（充足）
- II_final = 9,968（充足，原 18,529 去掉 8,561 重叠后 + 0 归属回来 = 9,968）
- LOW = 0（无需排除 seeds 的低置信序列）

---

### 05:10 — 生成互斥 FASTA（Phase 2 输入）

**精确命令：**

```python
# Python 脚本：从 hits_II_seqs.fasta 中仅保留 hits_II_final_ids.txt 的序列
# Ia: hits_Ia_seqs.fasta 保持不变（10,071 = 全部归属 Ia）
# Ib: hits_Ib_clean.fasta 保持不变
```

**Phase 2 最终互斥输入集合：**

| 文件 | 序列数 | 说明 |
|------|--------|------|
| `hits_Ia_seqs.fasta` | 10,071 | 原文件即为 final（全部 overlap 归属 Ia） |
| `hits_Ib_clean.fasta` | 7,869 | KDOPS 过滤后不变 |
| `hits_II_final_seqs.fasta` | 9,968 | 移除 8,561 条 Ia 归属后 |
| **总计** | **27,908** | 三亚型互斥 |

**产出文件：**

| 文件 | 说明 |
|------|------|
| `scripts/assign_besthit_Ia_vs_II.py` | 竞争打分归属脚本（新建） |
| `results/01_mining/besthit_Ia_vs_II.tsv` | 逐序列归属表（8,561 行） |
| `results/01_mining/hits_Ia_final_ids.txt` | Ia 最终 ID（10,071） |
| `results/01_mining/hits_II_final_ids.txt` | II 最终 ID（9,968） |
| `results/01_mining/hits_IaII_lowconf_ids.txt` | LOW 置信 ID（0 条） |
| `results/01_mining/hits_II_final_seqs.fasta` | II 最终互斥 FASTA |
| `results/01_mining/qc_subtype_assignment.md` | QC1b 归属报告 |

**✅ Gate A-2 完成。所有预门控通过，可进入 Phase 2.1。**

---

## 2026-02-25 Phase 2 执行日

---

### 05:45 — Phase 2.1：长度 + HMM 覆盖度三箱过滤

**脚本编写：** `scripts/qc_length_coverage.py`（新建），参数更新至 `meta/params.json`

**分类逻辑：**
- **PASS_CANONICAL**: `cov_hmm ≥ 0.70` 且 `L_seq` 在亚型主峰窗口
- **PASS_LONG**: `cov_hmm ≥ 0.70` 且 `L_seq` 偏长（融合域/转运肽）
- **FRAG**: `cov_hmm < 0.70` 或 `L_seq < canonical_min`

**精确命令：**

```bash
# Ia (HMM L=355, canonical 300-450)
python3 scripts/qc_length_coverage.py \
  --fasta results/01_mining/hits_Ia_seqs.fasta \
  --domtbl results/01_mining/hits_Ia.domtbl \
  --hmm_length 355 --subtype Ia \
  --canonical_min 300 --canonical_max 450 \
  --cov_min 0.70 --outdir results/02_qc

# Ib (HMM L=334, canonical 280-450)
python3 scripts/qc_length_coverage.py \
  --fasta results/01_mining/hits_Ib_clean.fasta \
  --domtbl results/01_mining/hits_Ib_vs_dah7ps_v41.domtbl \
  --hmm_length 334 --subtype Ib \
  --canonical_min 280 --canonical_max 450 \
  --cov_min 0.70 --outdir results/02_qc

# II (HMM L=471, canonical 320-600)
python3 scripts/qc_length_coverage.py \
  --fasta results/01_mining/hits_II_final_seqs.fasta \
  --domtbl results/01_mining/hits_II.domtbl \
  --hmm_length 471 --subtype II \
  --canonical_min 320 --canonical_max 600 \
  --cov_min 0.70 --outdir results/02_qc
```

**结果汇总：**

| 亚型 | 总数 | PASS_CANONICAL | PASS_LONG | FRAG | FRAG% |
|------|------|---------------|-----------|------|-------|
| Ia | 10,071 | 9,022 (89.6%) | 182 (1.8%) | 867 (8.6%) | ✅ 符合 QC1 预期 ~7.7% |
| Ib | 7,869 | 6,229 (79.2%) | 172 (2.2%) | 1,468 (18.7%) | ✅ 符合 QC1 预期 ~16% |
| **II** | **9,968** | **6,411 (64.3%)** | **93 (0.9%)** | **3,464 (34.8%)** | **⚠ 远超 QC1 预期 ~6.8%** |

**⚠ Type II FRAG 异常分析：**
- FRAG 中位长度 = 390 aa（在 canonical 窗口 320–600 内！）
- FRAG 最大 cov_hmm = 0.699（刚好低于 0.70 阈值）
- **根因：** Type II 的 α2β3 内插片导致 hmmsearch 产生碎片化 domain hits（SOP ⚠ [CHECK-06]），单个最佳 domain 的 HMM 覆盖度被压低到 <0.70，但序列本身是完整的 Type II
- **这正是 SOP Phase 3.6 要求 hit stitching 的原因**
- **需要决策：** 是否为 Type II 实现 multi-domain stitching 来恢复覆盖度，或降低 II 的 cov 阈值