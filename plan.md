# 🧬 DAH7PS 变构起源的演化动力学重建 SOP (V3.1)

**科学问题：** DAH7PS 酶的变构调节作用是如何在保守的 TIM barrel 基础上演化出来的？

**V3.1 核心策略：**

1. **防 KDOPS 污染**：双 HMM 竞争打分，剔除 Type Iβ 中 KDOPS 假阳性
2. **数据驱动 QC**：基于长度直方图设定过滤阈值
3. **Seed & Add MSA**：CD-HIT 60% 提取代表种子 → L-INS-i 精密骨架 → 增量映射
4. **PROMALS3D 混合骨架**：PDB 结构锚点 + 进化踏脚石 → 万级映射
5. **ClipKIT kpi-smart-gap**：保留进化信息 Gap，守住变构铰链
6. **AltAll 系综采样**：对抗 ML 过度稳定化偏差
7. **全原子 MD**：打破 AlphaFold 静态偏差，捕捉动力学预适应
8. **ICDC 映射**：DCA × DCCM 交汇 = 变构通讯电缆现形

**运行环境：** Ubuntu Linux / Conda `dah7ps_evo`

------

## 💻 第零步：计算环境搭建

### 0.1 安装 Miniforge

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
source ~/.bashrc
```

### 0.2 创建并激活虚拟环境

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

mamba create -n dah7ps_evo python=3.10 hmmer mafft iqtree trimal cd-hit seqkit matplotlib -y
conda activate dah7ps_evo
pip install clipkit
```

------

## 🔍 第一步：全长序列挖掘与 KDOPS 免疫过滤

### 1.1 构建 HMM 模型

```bash
mafft --maxiterate 1000 --localpair seeds_Ia.fasta > aligned_seeds_Ia.fasta
mafft --maxiterate 1000 --localpair seeds_Ib.fasta > aligned_seeds_Ib.fasta
mafft --maxiterate 1000 --localpair seeds_II.fasta > aligned_seeds_II.fasta

hmmbuild model_Ia.hmm aligned_seeds_Ia.fasta
hmmbuild model_Ib.hmm aligned_seeds_Ib.fasta
hmmbuild model_II.hmm aligned_seeds_II.fasta
```

### 1.2 撒网搜库与序列提取

```bash
hmmsearch --cpu 20 --domtblout domhits_Ia.tbl model_Ia.hmm uniref90.fasta.gz > /dev/null
hmmsearch --cpu 20 --domtblout domhits_Ib.tbl model_Ib.hmm uniref90.fasta.gz > /dev/null
hmmsearch --cpu 20 --domtblout domhits_II.tbl model_II.hmm uniref90.fasta.gz > /dev/null

awk '!/^#/ && $12 < 1e-10 {print $1}' domhits_Ia.tbl | sort | uniq > ids_Ia.txt
awk '!/^#/ && $12 < 1e-10 {print $1}' domhits_Ib.tbl | sort | uniq > ids_Ib.txt
awk '!/^#/ && $12 < 1e-10 {print $1}' domhits_II.tbl | sort | uniq > ids_II.txt

seqkit grep -f ids_Ia.txt uniref90.fasta.gz > raw_full_Ia.fasta
seqkit grep -f ids_Ib.txt uniref90.fasta.gz > raw_full_Ib.fasta
seqkit grep -f ids_II.txt uniref90.fasta.gz > raw_full_II.fasta
```

### 1.3 KDOPS 反向过滤（仅 Type Iβ）

```bash
cd-hit -i kdo8ps_uniprot.fasta -o kdo8ps_nr90.fasta -c 0.9 -n 5 -M 4000 -T 8
mafft --auto kdo8ps_nr90.fasta > kdo8ps_aligned.afa
hmmbuild kdo8ps.hmm kdo8ps_aligned.afa

hmmsearch --cpu 8 --domtblout domhits_Ib_vs_kdops.tbl kdo8ps.hmm raw_full_Ib.fasta > /dev/null
hmmsearch --cpu 8 --domtblout domhits_Ib_vs_dah7ps.tbl model_Ib.hmm raw_full_Ib.fasta > /dev/null

python filter_kdops.py
seqkit grep -v -f kdops_contaminants.txt raw_full_Ib.fasta > raw_full_Ib_clean.fasta
```

------

## 📊 第二步：数据驱动 QC1

### 2.1 长度过滤 + CD-HIT 去冗余

```bash
# 数据驱动阈值（基于直方图）
seqkit seq -m 300 -M 450 raw_full_Ia.fasta > qc_len_Ia.fasta       # Iα
seqkit seq -m 300 -M 480 raw_full_Ib_clean.fasta > qc_len_Ib.fasta  # Iβ
seqkit seq -m 380 -M 650 raw_full_II.fasta > qc_len_II.fasta        # II

# CD-HIT 80% 去冗余
cd-hit -i qc_len_Ia.fasta -o caseA_full_Ia.fasta -c 0.8 -n 5 -M 4000 -T 8
cd-hit -i qc_len_Ib.fasta -o caseA_full_Ib.fasta -c 0.8 -n 5 -M 4000 -T 8
cd-hit -i qc_len_II.fasta -o caseA_full_II.fasta -c 0.8 -n 5 -M 4000 -T 8
```

------

## 🧵 第三步：结构感知 MSA — Seed & Add + 混合骨架

### 🟢 Phase 1: 情况 A（亚型内 MSA）— Seed & Add 极限精度

*(以 Type Iβ 5,728 条为例，Iα 和 II 同理)*

**原理：** 将 O(N²) 暴力计算降维为 O(N)。先用 L-INS-i 死磕代表种子骨架，再增量映射。

```bash
# 1. CD-HIT 60% 提取进化均匀分布的种子
cd-hit -i caseA_full_Ib.fasta -o seeds60_Ib.fasta -c 0.6 -n 4 -M 4000 -T 8
# (预期 ~500-800 条种子)

# 2. L-INS-i 最高精度比对种子骨架
mafft --localpair --maxiterate 1000 --thread -1 seeds60_Ib.fasta > aligned_seeds60_Ib.afa

# 3. 提取剩余序列
seqkit seq -n seeds60_Ib.fasta | awk '{print $1}' > seed_ids_Ib.txt
seqkit grep -v -f seed_ids_Ib.txt caseA_full_Ib.fasta > remaining_Ib.fasta

# 4. 增量映射（保护骨架不变）
mafft --add remaining_Ib.fasta --thread -1 aligned_seeds60_Ib.afa > msa_full_Ib.afa
```

*(✅ 产出：`msa_full_Ia/Ib/II.afa` — EVcouplings/DCA 专用全长矩阵)*

### 🔴 Phase 2: 情况 B（跨三大亚型 1.2 万条）— 结构辅助混合骨架

**原理：** PDB 结构锚点提供三维刚性约束 + 进化踏脚石（Stepping stones）消除 Mapping Cliff 风险。

```bash
# 1. 合并三亚型种子，深度聚类为 ~300 条踏脚石
cat seeds60_Ia.fasta seeds60_Ib.fasta seeds60_II.fasta > all_seeds_mixed.fasta
cd-hit -i all_seeds_mixed.fasta -o stepping_stones.fasta -c 0.4 -n 2 -M 4000 -T 8

# 2. 与 PDB 种子合并
cat stepping_stones.fasta PDB_seeds.fasta > skeleton_raw_12k.fasta
```

**👉 手动干预：** 将 `skeleton_raw_12k.fasta` 提交 **PROMALS3D 网页端**，在高级选项中输入 PDB ID（1KFL, 1RZM, 3NV8 等）启用 3D 锚定。下载结果保存为 `promals3d_skeleton.afa`。

```bash
# 3. 提取剩余 1.1 万条
cat caseA_full_*.fasta > all_12k_mixed.fasta
seqkit seq -n skeleton_raw_12k.fasta | awk '{print $1}' > skel_ids_12k.txt
seqkit grep -v -f skel_ids_12k.txt all_12k_mixed.fasta > remaining_12k.fasta

# 4. 万级映射到 3D 骨架
mafft --add remaining_12k.fasta --thread -1 promals3d_skeleton.afa > global_alignment_raw.afa
```

### 🟣 Phase 3: ClipKIT 进化保护级修剪

```bash
# kpi-smart-gap：保留系统发育信息位点，哪怕充满 Gap 的变构 Loop 只要有进化特征就保留
clipkit global_alignment_raw.afa -m kpi-smart-gap -o msa_global_smart.afa
```

### 🛡️ 【QC 2：催化残基核验】

将 `msa_global_smart.afa` 导入 **Jalview**：K97、R165 等催化残基在跨三大 Type 的几千条序列中必须排在同一垂直列。

------

## 🌳 第四步：嵌套式 ASR + AltAll 采样

### 4.1 嵌套式 ASR

```bash
# 全局树
iqtree2 -s msa_global_smart.afa -m MFP -B 1000 -T AUTO --prefix Global_Tree

# 局部 ASR（含完整变构区的祖先）
iqtree2 -s msa_full_Ib.afa -te Subtree_Ib.treefile -m MFP -asr --prefix ASR_Ib
```

### 4.2 AltAll 序列生成

`generate_altall_ancestor.py`：解析 `.state` 文件

- **Best-ML**：每位点取 PP 最高氨基酸
- **AltAll**：0.2 < PP < 0.8 的模糊位点强制替换为第二高概率氨基酸

### 🛡️ 【QC 3】Bootstrap ≥ 95% / 变构面残基 PP > 0.85

------

## ⚛️ 第五步：MD 系综与动力学通讯

### 5.1 AlphaFold 3 → Anc_Pre / Anc_Post 初始 PDB
### 5.2 GROMACS/OpenMM：Apo vs Holo，≥500 ns × 3
### 5.3 RMSF（动力学预适应）/ 隐蔽口袋 / DCCM

------

## 🕸️ 第六步：ICDC — DCA × DCCM 交汇映射

进化偶联网络 × 物理动态网络 = 变构通讯电缆现形

### 🛡️ 【QC 4】pLDDT > 85 / MD RMSD 收敛