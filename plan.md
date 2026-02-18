# 🧬 DAH7PS 酶变构演化与祖先序列重建 (ASR) 项目 SOP

**科学问题：** DAH7PS 酶的变构调节作用是如何在保守的催化核心（TIM barrel）基础上演化出来的？变构信号是如何与 TIM 桶建立物理通讯的？

**核心策略：** “剥离附件，保留主干” + 分治多序列比对 + 祖先重现。

**当前输入文件：** `seeds_Ia.fasta`, `seeds_Ib.fasta`, `seeds_II.fasta`, `uniref90.fasta.gz`

**运行环境：** Ubuntu Linux

------

## 💻 第零步：Ubuntu 生信计算环境搭建

为了避免环境冲突和漫长的解析等待，我们使用 `Miniforge`（自带提速版的 Mamba）来管理环境。

### 0.1 安装 Miniforge

打开 Ubuntu 终端，执行以下命令：

Bash

```
# 下载安装脚本
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

# 运行安装（一路按 Enter，遇到 prompt 输入 yes）
bash Miniforge3-Linux-x86_64.sh

# 安装完成后，刷新环境变量
source ~/.bashrc
```

### 0.2 创建并激活项目专属虚拟环境

我们将项目所需的底层软件和序列处理神器 `seqkit` 一次性装好：

Bash

```
# 配置 Bioconda 通道
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# 创建名为 dah7ps_evo 的环境，并安装核心软件
mamba create -n dah7ps_evo python=3.10 hmmer mafft iqtree trimal cd-hit seqkit -y

# 激活环境（⚠️ 注意：以后每次开始做实验前，都要先运行这一句！）
conda activate dah7ps_evo
```

------

## 🔍 第一步：全长序列的系统挖掘与清洗

**目标：** 利用 3 个种子库，从庞大的 UniRef90 中精准“钓”出所有 DAH7PS 的全长同源序列。

### 1.1 构建 HMM 模型

Bash

```
# 1. 对三种 Type 的种子分别进行高精度多序列比对 (MSA)
mafft --maxiterate 1000 --localpair seeds_Ia.fasta > aligned_seeds_Ia.fasta
mafft --maxiterate 1000 --localpair seeds_Ib.fasta > aligned_seeds_Ib.fasta
mafft --maxiterate 1000 --localpair seeds_II.fasta > aligned_seeds_II.fasta

# 2. 将比对结果转化为隐马尔可夫模型 (.hmm)
hmmbuild model_Ia.hmm aligned_seeds_Ia.fasta
hmmbuild model_Ib.hmm aligned_seeds_Ib.fasta
hmmbuild model_II.hmm aligned_seeds_II.fasta
```

### 1.2 撒网搜库与序列提取

*注意：扫描大数据库计算量极大，建议保持电脑不休眠。此处以 Type Iα 为例演示流程（Ib 和 II 同理操作）。*

Bash

```
# 1. 扫描 UniRef90，输出包含结构域边界的详细表格 (--domtblout)
hmmsearch --cpu 8 --domtblout domhits_Ia.tbl model_Ia.hmm uniref90.fasta.gz > /dev/null

# 2. 提取高置信度（E-value < 1e-10）的序列 ID，并去重
awk '!/^#/ && $12 < 1e-10 {print $1}' domhits_Ia.tbl | sort | uniq > ids_Ia.txt

# 3. 使用高效的 SeqKit 瞬间提取全长序列
seqkit grep -f ids_Ia.txt uniref90.fasta.gz > raw_full_Ia.fasta
```

### 🛡️ 【QC 1：质量控制与降噪】（防垃圾输入）

提取的序列中必然混有测序碎片或异常长序列。

Bash

```
# 1. 长度过滤 (假设 Type Iα 正常长度约 350 aa，保留 300~400 之间的序列)
seqkit seq -m 300 -M 400 raw_full_Ia.fasta > length_qc_Ia.fasta

# 2. CD-HIT 去冗余：防止进化树被大肠杆菌等过度测序的物种“绑架”，聚类到 80% 相似度
cd-hit -i length_qc_Ia.fasta -o caseA_full_Ia.fasta -c 0.8 -n 5 -M 4000 -T 8

# ✅ 最终输出的 `caseA_full_Ia/Ib/II.fasta` 即为【情况 A】的高质量输入。
```

------

## ✂️ 第二步：“剥离附件”与分治法多序列比对

**目标：** 切除非同源的别构附件，仅保留纯 TIM 桶用于构建全局进化骨架树（情况 B）。

### 2.1 极速智能裁切 TIM 桶 (制备情况 B)

这里教你一个**高阶生信技巧**：第一步生成的 `domhits_Ia.tbl` 表格中，第 18 和 19 列精确记录了 TIM 桶在全长序列上的起止坐标！我们可以直接用命令将其精准切下，无需手写复杂的 Python 脚本：

Bash

```
# 1. 提取 ID 和 TIM 桶的起止坐标，制作 BED 格式文件 (注意 awk 中坐标减 1 适应 0-based 规则)
awk '!/^#/ && $12 < 1e-10 {print $1"\t"$18-1"\t"$19}' domhits_Ia.tbl > tim_boundaries_Ia.bed

# 2. 使用 SeqKit 根据坐标精准切除别构附件，只保留 TIM 桶
seqkit subseq --bed tim_boundaries_Ia.bed caseA_full_Ia.fasta > pure_TIM_Ia.fasta
# (同理对 Ib 和 II 剥离附件，得到 pure_TIM_Ib.fasta, pure_TIM_II.fasta)
```

### 2.2 分治法执行 MSA

Bash

```
# 【情况 A】全长序列独立比对：仅在同 Type 内部比对，用于后续共进化分析
mafft --maxiterate 1000 --localpair caseA_full_Ia.fasta > msa_full_Ia.fasta

# 【情况 B】纯 TIM 桶全局比对：合并三大亚型，进行跨 Type 极难比对
cat pure_TIM_*.fasta > all_pure_TIM.fasta
mafft --maxiterate 1000 --localpair all_pure_TIM.fasta > msa_global_TIM.fasta
```

### 🛡️ 【QC 2：比对剪裁与残基核验】（生死攸关的一步）

未修剪的进化比对包含大量无信息的 Gap（-），这会让建树模型崩溃。

Bash

```
# 1. 使用 trimAl 自动修剪全是 Gap 的劣质列
trimal -in msa_global_TIM.fasta -out trimmed_global_TIM.fasta -gappyout
```

- **2. 人工严酷质控：** 在你的个人电脑上安装 **Jalview** 软件。导入 `trimmed_global_TIM.fasta`。仔细寻找结合底物 PEP 的极度保守残基（例如大肠杆菌 AroG 序列中的 K97 或 R165）。
- **通过标准：这些关键催化残基必须在几千条序列中，严丝合缝地排在完美的一条垂直列上！** 如果错位了，立刻停止，回去重新调整 HMM 边界。

------

## 🌳 第三步：系统发育树构建与祖先序列重建 (ASR)

**目标：** 基于纯 TIM 桶恢复亿万年的垂直演化史，并计算古老祖先的“氨基酸字母”。

### 3.1 构建高精度系统发育树

使用目前学术界公认最强的建树软件 `IQ-TREE 2`。

Bash

```
# -s: 输入情况 B 修剪后的比对文件
# -m MFP: 自动测试并寻找最佳氨基酸演化模型 (机器自动挑，极其关键)
# -B 1000: 运行 1000 次 UltraFast Bootstrap 置信度检验
# -asr: 核心魔法参数！开启祖先序列推断
# -T AUTO: 自动调用所有 CPU 核心
iqtree2 -s trimmed_global_TIM.fasta -m MFP -B 1000 -asr -T AUTO --prefix Evolution_TIM
```

### 🛡️ 【QC 3：历史置信度评估】

运行结束后，你会得到 `.treefile` (进化树) 和 `.state` (祖先概率文件)。

1. **拓扑结构 QC：** 将 `.treefile` 拖入在线可视化网站 **iTOL** (itol.embl.de)。定位到三大亚型分道扬镳的“树杈”，其 Bootstrap 支持率必须 **$\ge 95\%$**。
2. **ASR 概率 QC：** 打开 `.state` 文件，找到“获得别构域前夕”的祖先节点。检查那些未来将变成“别构结合面”的氨基酸位点，其推断的**后验概率（Posterior Probability）最好 $> 0.85$**。如果概率过低，说明该位点历史信息已模糊，不能用于后续解释。

------

## 🏗️ 第四步：结构复原与共进化网络分析

**目标：** 在 3D 物理空间中揭示变构信号的通讯线缆是如何铺设的。

### 4.1 AlphaFold 3 祖先化石结构复原

1. 从上一步提取两个关键节点序列：**Anc_Pre** (刚获得别构域前夕的纯 TIM 桶祖先) 和 **Anc_Post** (刚接上别构域后的嵌合早期祖先)。
2. 将序列提交至 **AlphaFold 3** 官方网页端预测 3D 结构。
3. **机制洞察：** 在 PyMOL 中对齐两个结构。仔细观察在 Anc_Pre 时期的 TIM 桶表面，是否已经预先演化出了一片“疏水斑块”或“柔性 Loop”，随时准备迎接未来别构附件的降临？这就是进化中迷人的**“预适应 (Pre-adaptation)”** 现象。

### 4.2 直接偶联分析 (DCA) 揭示通讯网络

1. 取出**情况 A 的 MSA**（如 `msa_full_Ib.fasta`），提交至 **EVcouplings** 网页服务器进行共演化分析。
2. **核心关注：** 寻找跨越“变构附件”和“TIM 桶催化中心”的强烈协同突变残基对（Epistasis）。如果残基 A（在别构域）突变，残基 B（在催化中心）必定跟着突变以维持稳定，你就成功挖出了大自然设计的一条隐藏变构信号“传导高速公路”。

### 🛡️ 【QC 4：结构置信度核验】

- 检查 AlphaFold 的 **pLDDT** 分数。核心 TIM 桶必须是深蓝色（$> 85$ 分）。如果是面条状的红色（$< 50$ 分），说明序列重建不合理导致蛋白无法折叠。