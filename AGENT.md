# AGENT.md — DAH7PS 变构起源演化动力学重建 (V3.1)

## 项目概述

研究 DAH7PS 酶变构调节在 TIM barrel 上的演化起源。

**V3.1 策略：**
1. **KDOPS 反向过滤** — 双 HMM 竞争打分剔除假阳性
2. **Seed & Add MSA** — CD-HIT 60% 种子 → L-INS-i 骨架 → 增量映射（O(N²)→O(N)）
3. **PROMALS3D 混合骨架** — PDB 锚点 + 进化踏脚石，消除 Mapping Cliff
4. **ClipKIT kpi-smart-gap** — 保留变构铰链进化信息 Gap
5. **AltAll 系综** — 对抗 ML 过度稳定化偏差
6. **MD + ICDC** — AF3 → GROMACS ≥500ns → DCA×DCCM

**关键文献：** Jiao 2020 (内部插片) / Yokoyama 2025 (KDOPS) / Steenwyk 2020 (ClipKIT)

## 运行环境

- **系统：** Ubuntu Linux，28 核 CPU
- **Conda 环境：** `dah7ps_evo`（`/home/tynan/miniforge3/envs/dah7ps_evo`）
- **激活：** `conda activate dah7ps_evo`
- **CPU 规范：** 所有多线程工具统一使用 **20 核**
  - `hmmsearch --cpu 20`
  - `mafft --thread 20`
  - `cd-hit -T 20`
  - `iqtree2 -T 20`

## 工具链

| 工具 | 用途 | CPU 参数 |
|------|------|----------|
| HMMER | HMM 建模 & 序列搜索 | `--cpu 20` |
| MAFFT | Seed & Add MSA（L-INS-i + --add）| `--thread 20` |
| PROMALS3D | 结构引导骨架比对 | 网页端 |
| ClipKIT | 智能修剪（kpi-smart-gap）| 单线程 |
| IQ-TREE 2 | 嵌套式 ASR | `-T 20` |
| CD-HIT | 去冗余 & 种子提取 | `-T 20` |
| SeqKit | 序列处理 | 自动 |
| GROMACS/OpenMM | MD 模拟 | 待配置 |
| EVcouplings/plmDCA | DCA 共进化 | 待配置 |

## 目录结构

```
/home/tynan/0218/
├── seeds_Ia/Ib/II.fasta         # 原始种子
├── model_Ia/Ib/II.hmm           # DAH7PS HMM
├── uniref90.fasta.gz            # UniRef90 (~47 GB)
├── kdo8ps.hmm                   # KDOPS 诱饵 HMM
├── filter_kdops.py              # KDOPS 过滤脚本
├── kdops_contaminants.txt       # KDOPS 污染 ID
├── raw_full_Ib_clean.fasta      # 干净 Iβ 序列
├── caseA_full_*.fasta           # 情况 A: QC 后全长序列
├── seeds60_*.fasta              # 60% 聚类种子（Seed & Add 用）
├── aligned_seeds60_*.afa        # L-INS-i 种子骨架
├── msa_full_*.afa               # 情况 A MSA（DCA 用）
├── stepping_stones.fasta        # 跨亚型踏脚石（40% 聚类）
├── skeleton_raw_12k.fasta       # 踏脚石 + PDB 混合骨架
├── promals3d_skeleton.afa       # PROMALS3D 结构对齐骨架
├── global_alignment_raw.afa     # 情况 B 万级原始比对
├── msa_global_smart.afa         # ClipKIT 修剪后比对
├── plan.md / AGENT.md / log.md  # 文档
└── length_distribution.png      # 直方图
```

## 关键禁止事项

- ❌ 物理切断 Type Iα/II 序列（内部插片）
- ❌ 使用 `trimAl -gappyout`（丢失变构铰链）
- ❌ 对 3k+ 序列直接跑 `--localpair`（内存溢出风险）
- ❌ 跨亚型骨架仅用 PDB（Mapping Cliff 风险）
- ❌ 使用 `--thread -1`（MAFFT 自动检测不准，默认 8 核）
- ✅ 统一使用 `--thread 20` / `--cpu 20` / `-T 20`
- ✅ Type Iβ 必须用 KDOPS 过滤后的 `raw_full_Ib_clean.fasta`
- ✅ ASR 同时生成 Best-ML + AltAll 做鲁棒性对照
