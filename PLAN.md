# DAH7PS 变构起源的演化动力学重建：V4.1 标准作业程序（SOP rev4）

**A FoldMason-Centred, Core–Module Decoupling Framework for Reconstructing the Evolutionary Origin of Allostery in DAH7PS**

---

## 摘要

DAH7PS（3-deoxy-D-arabino-heptulosonate 7-phosphate synthase）在保守的 (β/α)₈ TIM-barrel 催化核心之上，于不同亚型中独立演化出多样的变构调控元件——Type Iα 的 β 发夹插入、Type Iβ 的 ACT/铁氧还蛋白样调控域、Type II 的 N 端延伸与 α2β3 插片。V3.1 管线的"全长单一 MSA 驱动一切"策略在跨亚型远缘同源（序列一致性进入暮光区 <20%）时已产生严重比对膨胀（10–17:1），系统性污染下游 ASR、DCA 与动力学推断。

V4.1 采用**核心–模块分治（core–module decoupling）**重构分析范式：以 FoldMason 结构字母表（3Di+AA）与逐列 LDDT 物理置信度作为同源性定义与比对质量控制的中枢度量（非唯一判据），结合催化位点几何一致性、AA/3Di 树拓扑交叉验证、外群定根稳定性等多条独立证据链；将跨亚型可可靠同源的 TIM-barrel 核心与变构调控模块分开建模；**但在需要直接检测跨域共进化信号时（core↔module），在同一架构亚型内额外构建 profile-anchored 拼接比对并运行联合 DCA（Phase 3.9 & 6.2.3）**；在核心树上进行定根系统发育与氨基酸/gap 双轨 ASR；将变构元件视为离散性状追踪获得/丢失时间节点；最终通过分层 DCA × 多聚体 MD 动力学通讯网络整合（ICDC）建立演化耦联—物理通讯的一致证据链。

---

## V4.1 修订要点（相对 V4.0）

本版本将你提出的“致命断层/高风险陷阱/工程排雷”建议**正式并入 SOP**，并把关键点落实为可执行的 Phase、脚本契约与 QC 断言。主要改动如下（按影响优先级）：

1. **修复 `hmmalign --outformat afa` 的插入列膨胀**：核心映射改为 *Stockholm → 剥离 Insert 列 → 再转 AFA*，保证核心 MSA 长度严格等于 HMM match-state 长度（Phase 3.6，新增 `esl-alimask --rf-is-mask` 流程与 fallback）。
2. **补齐“亚型内全长 MSA”的生成链路**：新增 *Profile-anchored Stitching MSA*（Phase 3.9），明确 `msa_full_Ib_v4.afa` 等全长比对来自“核心锚定 + 模块锚定 + linker 自由对齐 + 拼接”，避免 V3.1 式膨胀噩梦回归。
3. **新增亚型内跨域联合 DCA（Joint DCA）**：在不破坏“核心–模块分治”总体范式的前提下，对同一架构亚型（如 Type Iβ-ACT）构建“核心+linker+模块”的联合比对，直接提取 core↔module 的耦联证据（Phase 6.2.X）。
4. **DCA 比较引入 Meff 匹配下采样与 Z-score**：`with_ACT` vs `without_ACT` 的耦联变化比较升级为 *Meff-matched downsampling + bootstrap 稳定性 +（可选）IDR*，避免系统发育/有效深度混杂（Phase 6.2.5）。
5. **外群定根加入位点异质性模型抗 LBA**：全局根位置除 `-m MFP` 外，强制再跑至少一个 site-heterogeneous mixture（如 `LG+C20+F+G` 或 `EX_EHO+F+G`）并把“根一致性”写入 QC3（Phase 4.1）。
6. **祖先 AF3 结构预测改为 Apo-first 门控**：祖先节点默认先做 Apo 预测；仅当口袋/对接支持“可容纳”时才允许做 Holo 预测与 Holo-MD，降低“现代配体幻觉”风险（Phase 5.2–5.3，新增 CHECK-07）。
7. **工程补丁**：新增 CHECK-06（Type II 内插片导致 HMM hit 断裂的 stitching 逻辑）；3Di 模型调用改为 `-mdef` 显式加载 Nexus model file；模块 DCA 的 Meff/L 底线提升到 3.0（理想 5.0），不足则禁跑模块 DCA。

---

## 科学问题与工作假设

**核心问题：** DAH7PS 的变构调控元件如何在保守 TIM-barrel 核心之上被募集、稳定并与多聚体界面及动力学网络耦合，从而产生可选择的变构表型？

**工作假设：**

1. 跨亚型可可靠比较的演化"时间轴"主要由 TIM-barrel 核心提供；变构元件的起源在演化上表现为模块获得/丢失事件与模块内序列—结构协变。
2. 变构功能的形成与维持依赖多聚体组装态与跨亚基通讯路径（Cross et al., 2013; Lang et al., 2016）。
3. 序列共进化（DCA）与动力学相关（DCCM/网络流）在关键演化节点应呈现可对齐的跨证据一致性。

---

## 方法学取舍声明

1. **FoldMason 定位：** 中枢工具，非唯一判据。"完美对齐"并非可验证的普遍命题；V4.1 保留多证据链交叉验证（催化位点几何一致性、AA 树 vs 3Di 树拓扑、外群定根稳定性、对齐长度与 gap 谱稳定性）。
2. **MAFFT 角色：** 不适合跨亚型无约束裸跑全长比对，但在亚型内/模块子集高精度 MSA（E-INS-i）以及受约束增量映射（`--add --keeplength`）场景仍具方法学价值。V4.1 采用"FoldMason 结构骨架 + hmmalign profile 映射（首选）+ MAFFT 作为备选/模块内工具"的组合。
3. **位点保护：** 不对特定残基三联体进行硬编码保护（避免确认偏差）。以 FoldMason-LDDT 高置信列为一级保护标准，功能注释/结构证据为二级验证。
4. **参数管理：** 所有阈值（LDDT 截断、gap fraction、聚类一致性等）外部化至 `meta/params.json`，不硬编码于脚本中。

---

# 第零步：计算环境与项目规范

## 0.1 可复现环境

**层 1：conda 核心环境（所有成员必须安装）**

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict

mamba create -n dah7ps_v4 \
  python=3.11 hmmer mafft iqtree2 mmseqs2 foldmason seqkit cd-hit \
  -y
conda activate dah7ps_v4

pip install clipkit pastml topiary-asr --break-system-packages
```

**层 2：HPC / 容器依赖（不在 conda 内安装，但在 `meta/software_versions.tsv` 中强制记录）**

```
工具              最低版本     安装方式
──────────────────────────────────────────────
GROMACS           ≥2023.x     HPC module 或 Singularity/Apptainer 容器
AlphaFold3        参见 DeepMind AF3 仓库
ESMFold           参见 ESM GitHub
PyMOL / ChimeraX  可选         结构可视化
```

## 0.2 软件版本锁定（强制）

```bash
mkdir -p results/meta
echo -e "tool\tversion" > results/meta/software_versions.tsv
# 层 1：conda 工具
for tool in hmmbuild hmmalign hmmsearch mafft iqtree2 mmseqs foldmason seqkit cd-hit clipkit; do
  echo -e "${tool}\t$($tool --version 2>&1 | head -1)" >> results/meta/software_versions.tsv
done
# 层 2：HPC 工具（若可用）
gmx --version 2>&1 | grep "GROMACS version" >> results/meta/software_versions.tsv || true
```

### 0.2.1 外部模型/矩阵文件锁定（强制）

除可执行程序版本外，**凡是通过下载获得的模型文件/矩阵文件**（例如 3Di 的 `Q.3Di.AF`/`Q.3Di.LLM` Nexus 定义文件）必须：

1. 统一存放到 `meta/models/`；
2. 记录来源（URL/commit/tag）、下载日期；
3. 记录校验和（建议 sha256）到 `results/meta/model_files.tsv`，确保可复现。

```bash
mkdir -p meta/models results/meta
echo -e "file	sha256	source" > results/meta/model_files.tsv
# 示例（将 <URL> 替换为实际来源链接）
sha256sum meta/models/Q_3Di_models.nex | awk '{print $2"\t"$1"\t""<URL or commit>"}' >> results/meta/model_files.tsv
```


## 0.3 目录规范（强制）

```
data/
  seeds/                # 手工整理的种子序列与结构
  structures/
    panel_dah7ps/       # ⚠ 仅 DAH7PS 三亚型代表（核心列界定面板）
    panel_kdops/        # ⚠ 仅 KDOPS 外群结构（定根面板）
  db/                   # UniRef/UniProt 数据库索引
results/
  01_mining/            # HMM 搜索 + KDOPS 过滤
  02_qc/                # 长度过滤 + 去冗余
  03_msa_core/          # 核心层 MSA（FoldMason 骨架 + hmmalign 映射）
  03_msa_modules/       # 模块层 MSA（分模块独立比对）
  04_phylogeny_asr/     # 系统发育 + 分层 ASR
  05_struct_md/         # 结构预测 + 多聚体 MD
  06_icdc/              # 分层 DCA × DCCM 整合
  meta/                 # 参数文件、版本锁定、QC 报告
workflow/
  Snakefile             # 建议使用 Snakemake 管理（非强制）
scripts/
```

> **路径一致性原则（关键）：** 结构面板采用两阶段策略后，必须严格区分：
> - `data/structures/panel_dah7ps/` — 仅含 DAH7PS 三亚型代表，用于 Phase 3.1-3.4 核心列界定
> - `data/structures/panel_kdops/` — 仅含 KDOPS 外群，在 Phase 4.1 通过 hmmalign 并入用于定根
> - 所有 FoldMason、结构复核、DCA contact_pdb、ICDC anchor_pdb 均应引用 `panel_dah7ps/1KFL.pdb`
> - 任何混用 `panel/` 或不带后缀的路径都会导致科学结果污染

## 0.4 参数文件模板

```json
// meta/params.json  （示例模板：所有阈值/模型选择必须在此文件中可追溯）
{
  "mining": {
    "evalue_threshold": 1e-10,
    "kdops_score_margin": 0
  },
  "qc": {
    "length_lower_percentile": 5,
    "length_upper_percentile": 95,
    "cdhit_identity_phase1": 0.80,
    "cdhit_identity_seeds": 0.60
  },
  "core_definition": {
    "pad_residues": 20,
    "lddt_min": "auto_inflection",
    "gap_fraction_max": 0.30,
    "note": "lddt_min 由 FoldMason LDDT 分布拐点确定，记录于 QC2"
  },
  "msa": {
    "hmmalign": {
      "strip_insertions": true,
      "strip_method": "esl-alimask --rf-is-mask",
      "fallback_method": "a2m_drop_lowercase"
    },
    "full_length_stitch": {
      "linker_align_method": "mafft_einsi",
      "architecture_scope": "subtype_specific_homogeneous"
    }
  },
  "phylogeny": {
    "bootstrap_replicates": 1000,
    "outgroup_prefix": "KDOPS_",
    "root_model_baseline": "MFP",
    "root_model_site_hetero": "LG+C20+F+G",
    "root_model_site_hetero_alt": "EX_EHO+F+G"
  },
  "dca": {
    "meff_min_main": 3.0,
    "meff_ideal": 5.0,
    "reweight_threshold": 0.8,
    "regularization_le": 16.0,
    "regularization_lh": 0.01,
    "comparison": {
      "downsample_reps": 20,
      "target_meff": "min(group_meff)"
    }
  },
  "af3": {
    "apo_first": true,
    "holo_requires_pocket": true
  },
  "asr": {
    "altall_pp1": 0.80,
    "altall_pp2": 0.20
  }
}
```


---

# 第一步：全长序列挖掘与 KDOPS 免疫过滤

## 1.1 种子策略（关键修正：Type II 种子多样性）

V3.1 的 `seeds_II.fasta` 仅含 1 条来自 *Mycobacterium tuberculosis* 的序列。以单一种子构建的 HMM 模型对 Type II 多样性覆盖不足，可能系统性遗漏来自放线菌、植物叶绿体等分支的变体。

**V4.1 要求：** 每一大类（Iα / Iβ / II）种子集覆盖主要分支多样性，建议 5–20 条。从 UniProt/InterPro（PF00793 家族）中按系统分类均匀取样。KDOPS 外群同样需多样性覆盖（10–20 条），用于后续定根与稳定性检验。

## 1.2 HMM 构建与搜库

```bash
# 种子对齐：小规模，精度优先
mafft --localpair --maxiterate 1000 data/seeds/seeds_Ia.fasta > results/01_mining/seeds_Ia.afa
mafft --localpair --maxiterate 1000 data/seeds/seeds_Ib.fasta > results/01_mining/seeds_Ib.afa
mafft --localpair --maxiterate 1000 data/seeds/seeds_II.fasta > results/01_mining/seeds_II.afa

# 构建 HMM
hmmbuild results/01_mining/model_Ia.hmm results/01_mining/seeds_Ia.afa
hmmbuild results/01_mining/model_Ib.hmm results/01_mining/seeds_Ib.afa
hmmbuild results/01_mining/model_II.hmm results/01_mining/seeds_II.afa

# 搜索 UniRef90
for type in Ia Ib II; do
  hmmsearch --cpu 20 -E 1e-10 \
    --domtblout results/01_mining/hits_${type}.domtbl \
    results/01_mining/model_${type}.hmm data/db/uniref90.fasta > /dev/null
done

# 从 domtblout 提取命中序列 ID → 从 UniRef90 中取出对应全长序列
for type in Ia Ib II; do
  # 提取去重的命中 ID（domtblout 第一列为 target name）
  grep -v "^#" results/01_mining/hits_${type}.domtbl \
    | awk '{print $1}' | sort -u > results/01_mining/hits_${type}_ids.txt

  # 从数据库中提取全长序列（seqkit 或 esl-sfetch）
  seqkit grep -f results/01_mining/hits_${type}_ids.txt \
    data/db/uniref90.fasta \
    > results/01_mining/hits_${type}_seqs.fasta
done
```

## 1.3 KDOPS 免疫过滤（双用途：过滤 + 外群回收）

对 Iβ 候选集使用 V3.1 已验证的 dual-HMM 竞争得分策略剔除 KDOPS 假阳性。同时，**独立构建高置信 KDOPS 外群小集合**，仅在第四步全局树定根时加入，downstream 分析（DCA、MD）一律剔除。

```bash
# 沿用 V3.1 的 filter_kdops.py 对 Iβ 做过滤
python scripts/filter_kdops.py \
  --dah7ps_hmm results/01_mining/model_Ib.hmm \
  --kdops_hmm data/seeds/model_kdops.hmm \
  --input results/01_mining/hits_Ib_seqs.fasta \
  --output results/01_mining/hits_Ib_clean.fasta

# 独立获取 KDOPS 外群（10–20 条，覆盖主要分支）
# → data/seeds/kdops_outgroup.fasta
```

### QC1（强制产出）

`results/01_mining/qc_mining_report.md`：
- 每亚型命中数量、e-value 分布、长度分布直方图
- KDOPS 竞争得分分布（Iβ）
- 分类群覆盖（门/纲级别）

---

# 第二步：数据驱动 QC 与去冗余

## 2.1 长度过滤（数据驱动阈值）

以每类长度直方图确定下限与上限，阈值由脚本记录至 `meta/params.json`，不硬编码。

**片段序列处理：** 明显缺失核心片段的序列进入"缺失集（fragments bin）"——不进入核心树推断，但可用于模块存在率统计（避免系统偏倚）。

```bash
python scripts/qc_length.py \
  --input results/01_mining/hits_Ib_clean.fasta \
  --output_pass results/02_qc/qc1_Ib.fasta \
  --output_fragments results/02_qc/fragments_Ib.fasta \
  --params meta/params.json
```

## 2.2 去冗余（两级）

```bash
# 第一级：80% 去冗余，降低计算成本
cd-hit -i results/02_qc/qc1_Ib.fasta -o results/02_qc/nr80_Ib.fasta -c 0.80 -T 20

# 第二级：60% 种子代表（用于种子骨架构建）
cd-hit -i results/02_qc/nr80_Ib.fasta -o results/02_qc/seeds60_Ib.fasta -c 0.60 -T 20
```

## 2.3 跨亚型 stepping stones（MMseqs2 替代 CD-HIT）

CD-HIT 在 <60% 一致性时敏感度急剧下降（Steinegger & Söding, 2018）。跨亚型代表选取使用 MMseqs2：

```bash
# 合并三亚型种子代表
cat results/02_qc/seeds60_Ia.fasta \
    results/02_qc/seeds60_Ib.fasta \
    results/02_qc/seeds60_II.fasta > results/02_qc/all_seeds_mixed.fasta

mmseqs easy-cluster results/02_qc/all_seeds_mixed.fasta \
  results/02_qc/stepping_stones tmp_mmseqs \
  --min-seq-id 0.4 -c 0.8 --cov-mode 1
```

**定量验收标准（强制）：**

1. **数量范围：** 每亚型 5–15 个代表，跨亚型总计 20–50 个（结构面板的计算可行范围）。若聚类结果超出此范围，调整 `--min-seq-id`（降低以减少代表数，升高以增加代表数）并记录至 `meta/params.json`。
2. **覆盖度验证：** 将 stepping stones 映射到亚型内初始树（Phase 1 种子建树即可）上，检查每个主要分支（bootstrap ≥ 70 的内部节点）是否至少有 1 个代表。若某主要分支无覆盖，手动补充该分支的代表序列。
3. **结构面板衔接规则：** stepping stones 中有实验结构的直接纳入 Phase 3.1 面板；无实验结构的用 ESMFold/ColabFold 预测后纳入，但需 pLDDT ≥ 70（全局平均），pLDDT < 70 的预测结构不进入面板。
4. **参数记录：** 最终代表数、每亚型分布、手动补充记录写入 `results/02_qc/stepping_stones_report.md`。

---

# 第三步：结构感知 MSA——核心/模块分治，FoldMason 为中枢

> 本步产出分为两套：
> **(A) core MSA**：跨亚型可靠同源的 TIM-barrel 核心 → 核心树、核心 ASR、核心 DCA
> **(B) module MSAs**：按模块类型与携带子集分别构建 → 模块内 ASR、模块内 DCA

## Phase 3.1：结构面板构建（两阶段策略）

> **设计原则：** "核心列界定"与"外群定根"的目标不同，需要在面板构建层面被显式解耦。若将 KDOPS 纳入核心列界定面板，KDOPS 与 DAH7PS 在功能性 loop 和界面片段上的结构差异会系统性拉低这些位置的 LDDT，导致核心列集合 C 被过度压缩为"DAH7PS + KDOPS 共有的最保守骨架"，丢失与变构耦联相关但仅在 DAH7PS 内保守的关键区段。

**阶段 A（核心列界定面板）：仅 DAH7PS 三亚型代表，不含 KDOPS。**

来源优先级：

1. 实验结构（PDB）：1KFL（Iα）、1RZM（Iβ）、3NV8/5CKV（II）等
2. AlphaFold Database 中已有的高可信预测结构（pLDDT ≥ 70）
3. 对缺口分支的 stepping stone 代表进行 ESMFold/ColabFold 补充预测

目标：20–40 个 DAH7PS 结构，覆盖三大亚型的主要分支。

```bash
mkdir -p data/structures/panel_dah7ps
# 下载/预测 DAH7PS 结构到此目录（PDB/mmCIF 格式，可 gzip）
# KDOPS 结构单独存放：
mkdir -p data/structures/panel_kdops
```

**阶段 B（外群定根面板）：** KDOPS 结构独立存放，仅在 Phase 4.1 通过 hmmalign 映射到已有核心 profile 后并入全局树用于定根。KDOPS 不参与核心列的定义过程。

## Phase 3.2：FoldMason 结构多重比对 → 结构骨架

```bash
# FoldMason easy-msa：仅使用 DAH7PS 面板（不含 KDOPS）
foldmason easy-msa data/structures/panel_dah7ps/* \
  results/03_msa_core/skeleton tmp_foldmason \
  --report-mode 1

# 输出（自动生成）：
# results/03_msa_core/skeleton_aa.fa    — AA 序列骨架比对
# results/03_msa_core/skeleton_3di.fa   — 3Di 字母表比对
# results/03_msa_core/skeleton.nw       — 引导树
# results/03_msa_core/skeleton.html     — 交互式逐列 LDDT 可视化
```

## Phase 3.3：FoldMason refinemsa 迭代优化

对骨架执行 LDDT 最大化 refinement（1000 轮迭代）。保留 refine 前版本作为对照。

```bash
# 创建结构数据库（refine 需要）
foldmason createdb data/structures/panel_dah7ps/* results/03_msa_core/panelDb

# 迭代优化
foldmason refinemsa results/03_msa_core/panelDb \
  results/03_msa_core/skeleton_aa.fa \
  results/03_msa_core/skeleton_refined_aa.fa \
  --refine-iters 1000

# 对比 refine 前后 LDDT
foldmason msa2lddt results/03_msa_core/panelDb results/03_msa_core/skeleton_aa.fa
foldmason msa2lddt results/03_msa_core/panelDb results/03_msa_core/skeleton_refined_aa.fa

# 生成 refine 后的 HTML 报告
foldmason msa2lddtreport results/03_msa_core/panelDb \
  results/03_msa_core/skeleton_refined_aa.fa \
  results/03_msa_core/skeleton_refined.html
```

**决策规则：** 若 refine 后关键催化位点列的对齐无实质变化（说明原始骨架已足够好），以原始骨架为准；若 refine 纠正了可视的催化位点错配，以 refined 版本为准。

## Phase 3.4：基于逐列 LDDT 的核心列界定

**原则：** 用结构一致性而非经验阈值定义"跨亚型可靠同源"列。

```python
# scripts/define_core_columns.py（伪代码）
# 1. 从 msa2lddt 或 HTML 解析逐列 LDDT 与 gap fraction
# 2. 绘制 LDDT 分布直方图——寻找双峰/拐点
# 3. 高 LDDT 列 = 结构保守核心；低 LDDT 列 = 可变插入/loop
# 4. 核心列集合 C: LDDT(col) ≥ inflection_point 且 gap_fraction(col) ≤ 0.30
# 5. 阈值写入 meta/params.json
# 6. 输出 skeleton_core_aa.fa（仅保留列集合 C）
```

## Phase 3.5：亚型内种子骨架重建（E-INS-i）

V3.1 使用 L-INS-i（`--localpair`）进行种子骨架对齐。DAH7PS 包含多个保守 motif（八条 β 链、八条 α 螺旋）被可变长度 loop 和调控域插入分隔，符合 E-INS-i 的设计场景——广义仿射 gap 代价模型允许非同源区域保持未对齐（Katoh et al., 2005）。

```bash
# 亚型内种子骨架（用于亚型内 hmmbuild）
for type in Ia Ib II; do
  mafft --genafpair --maxiterate 1000 --ep 0 --thread 20 \
    results/02_qc/seeds60_${type}.fasta \
    > results/03_msa_core/seeds60_${type}_einsi.afa
done
```

## Phase 3.6：全量序列去膨胀映射——hmmalign（首选，V4.1 修复版）

核心思想：先用核心骨架构建 profile HMM，再将每条序列的核心域片段映射到 profile 的 **match states**，从根本上消除“全长裸跑 MSA”带来的跨亚型膨胀。

> **V4.1 关键修正（必须遵守）：** `hmmalign` 若直接输出 `afa`，会把某些序列的长插入（Insert states）展开为显式对齐列，从而在“其他所有序列”中注入成百上千 gap，导致核心比对再次膨胀。  
> 因此核心映射必须采用：**Stockholm 输出 → 剥离 Insert 列 → 再转 AFA**，最终得到严格等于 HMM 骨架长度的 `core_global_matchonly.afa`。

```bash
# 1) 用 FoldMason 核心骨架构建跨亚型核心 HMM
hmmbuild results/03_msa_core/core_global.hmm results/03_msa_core/skeleton_core_aa.fa

# 2) 定位核心域区间并提取核心域片段
# ⚠ [CHECK-01] pad 20 aa：防止边界锯齿切碎 β1/α8
# ⚠ [CHECK-06] Type II α2β3 内插片会导致 hmmsearch 产生碎片化 hits：必须启用 hit stitching
python scripts/extract_core_domains.py --domtblout results/01_mining/hits_*.domtbl --sequences results/02_qc/nr80_*.fasta --pad 20 --hit_stitching auto --output results/03_msa_core/all_core_only.fasta --coords_out results/03_msa_core/core_domain_coords.tsv

# 3) hmmalign 映射：先输出 Stockholm（保留 Match/Insert 状态信息）
hmmalign --trim --mapali results/03_msa_core/skeleton_core_aa.fa -o results/03_msa_core/core_global_raw.sto results/03_msa_core/core_global.hmm results/03_msa_core/all_core_only.fasta

# 4) 剥离所有 Insert 列，只保留 HMM match-state 骨架长度（依赖 '#=GC RF' 注释）
esl-alimask --rf-is-mask results/03_msa_core/core_global_raw.sto > results/03_msa_core/core_global_matchonly.sto

# 5) 转换为 aligned FASTA (afa) 并统一大写
esl-reformat afa results/03_msa_core/core_global_matchonly.sto | seqkit seq --upper-case > results/03_msa_core/core_global_matchonly.afa
```

**Fallback（当 Stockholm 不含 `#=GC RF` 时）：** 改用 `a2m` 并显式删除插入态字符（插入态通常为小写）。

```bash
hmmalign --outformat a2m --trim --mapali results/03_msa_core/skeleton_core_aa.fa results/03_msa_core/core_global.hmm results/03_msa_core/all_core_only.fasta > results/03_msa_core/core_global_raw.a2m
python scripts/drop_a2m_insertions.py --input results/03_msa_core/core_global_raw.a2m --output results/03_msa_core/core_global_matchonly.a2m
esl-reformat afa results/03_msa_core/core_global_matchonly.a2m | seqkit seq --upper-case > results/03_msa_core/core_global_matchonly.afa
```

**备选方案（MAFFT --keeplength）：** 当 hmmalign 核心域定位失败率偏高时使用（列数由骨架锁死，不会膨胀）。

```bash
mafft --add results/03_msa_core/all_core_only.fasta --keeplength --thread 20 results/03_msa_core/skeleton_core_aa.fa > results/03_msa_core/core_global_matchonly_mafft.afa
```

esults/03_msa_core/core_global_matchonly_mafft.afa
```

## Phase 3.7：修剪（双版本）

系统发育建树与 ASR/DCA 的修剪目标不同，维护两套比对：

```bash
# 树推断用（kpic-smart-gap：保留信息位点 + 恒定位点 + 智能 gap 过滤）
clipkit results/03_msa_core/core_global_matchonly.afa \
  -m kpic-smart-gap \
  -o results/03_msa_core/core_tree.afa \
  --complementary

# ASR/DCA 用（最小修剪——仅去极端 gap 列，保留更多位点供模型估计）
python scripts/minimal_trim.py \
  --input results/03_msa_core/core_global_matchonly.afa \
  --gap_col_threshold 0.95 \
  --output results/03_msa_core/core_asr.afa
```

### 结构复核（用 FoldMason msa2lddt）

```bash
# 从修剪后比对中提取有结构的代表子集
python scripts/extract_struct_subset.py \
  --msa results/03_msa_core/core_tree.afa \
  --structures data/structures/panel_dah7ps/ \
  --output results/03_msa_core/core_tree_struct_subset.fa

foldmason msa2lddt results/03_msa_core/panelDb \
  results/03_msa_core/core_tree_struct_subset.fa

# 仅评估无 gap 列的 LDDT（--pair-threshold 1.0）
foldmason msa2lddt results/03_msa_core/panelDb \
  results/03_msa_core/core_tree_struct_subset.fa \
  --pair-threshold 1.0
```

### QC2（强制产出）

`results/03_msa_core/qc_core_alignment.md`：
- 比对长度（预期 400–600 列范围，vs V3.1 的 3,544–7,689）
- gap 谱分布
- 平均 LDDT（全列 + 仅无 gap 列）
- V3.1 旧比对 vs V4.1 新比对的 LDDT 定量对比
- 核心保守位点（基于结构/功能注释）在跨亚型是否保持同列一致

## Phase 3.8：模块定义与模块数据集构建

将变构元件从全局坐标系中剥离，构建三类产物：

```bash
# 模块注释（多证据融合）
python scripts/annotate_modules.py \
  --domtblout results/01_mining/hits_*.domtbl \
  --lddt_cols results/03_msa_core/lddt_per_column.tsv \
  --output results/03_msa_modules/module_presence_absence.tsv
  # 输出格式：seq_id | ACT_domain | N_ext | alpha2beta3_insert | CM_domain | transit_peptide

# 模块边界来源：
# - HMMER domain 命中坐标（env_from/env_to）
# - FoldMason LDDT 低置信区段
# - 结构注释/二级结构一致性
```

**模块边界稳健性 QC（强制）：**

> 模块获得/丢失事件将构成科学叙事的时间轴。任何系统性的模块边界漂移都可能被误读为"演化事件"，因此需显式记录边界稳健性。

1. **三证据一致性检查：** 对每条序列，独立用三种证据（HMMER 域坐标、LDDT 低置信区段、Pfam/InterPro 注释）分别判定模块边界，计算三种判据的一致率。一致率 < 80% 的序列标记为"边界模糊（boundary_ambiguous）"。

2. **敏感性分析（双版本模块矩阵）：** 分别构建"严格版"（三种证据一致才标记为模块存在）和"宽松版"（任一证据即标记为模块存在）的 `module_presence_absence.tsv`，后续 Phase 4.6 的 PastML 均在两版上运行。若模块起源节点在两版之间不稳定（推断出不同获得节点），该模块起源时间标记为"不确定"。

3. **构象态控制：** 结构面板中的实验结构标注 apo/holo/ligand-bound 状态。若模块边界的 LDDT 与构象态系统相关（如所有 holo 态结构的某段 loop LDDT 显著低于 apo 态），则该区段的模块归属应排除构象态效应后重新判定。

```bash
python scripts/annotate_modules.py \
  --domtblout results/01_mining/hits_*.domtbl \
  --lddt_cols results/03_msa_core/lddt_per_column.tsv \
  --pfam_annotations data/db/pfam_annotations.tsv \
  --output results/03_msa_modules/module_presence_absence_strict.tsv \
  --output_relaxed results/03_msa_modules/module_presence_absence_relaxed.tsv \
  --boundary_report results/03_msa_modules/boundary_robustness.md
  # 输出格式：seq_id | ACT_domain | N_ext | alpha2beta3_insert | CM_domain | transit_peptide | boundary_confidence
```

对每个模块，提取携带子集的模块片段：

```bash
for module in ACT N_ext alpha2beta3 CM; do
  python scripts/extract_module_seqs.py \
    --module_name ${module} \
    --presence results/03_msa_modules/module_presence_absence.tsv \
    --sequences results/02_qc/nr80_*.fasta \
    --output results/03_msa_modules/${module}_seqs.fasta \
    --coords_out results/03_msa_modules/${module}_domain_coords.tsv

  # 模块内部 MSA（MAFFT E-INS-i——适合含大型间隔的模块）
  mafft --genafpair --maxiterate 1000 --ep 0 --thread 20 \
    results/03_msa_modules/${module}_seqs.fasta \
    > results/03_msa_modules/${module}_msa.afa
done
```

**若模块结构可得，追加 FoldMason 模块层对齐：**

```bash
# 例如：ACT 域结构面板
foldmason easy-msa data/structures/ACT_panels/* \
  results/03_msa_modules/ACT_skeleton tmp_act --report-mode 1
```

## Phase 3.9：轮廓锚定缝合比对（Profile-anchored Stitching MSA）

**目的：** 为嵌套 ASR（Phase 4.4）与亚型内跨域联合 DCA（Phase 6.2.X）生成**高质量、非膨胀**的“亚型内全长比对”（例如 `msa_full_Ib_v4.afa`），补齐 V4.0 文档中隐含但未定义的关键输入链路。

**原则：**
- **核心列绝不允许重新对齐**：核心段必须直接继承 `core_global_matchonly.afa`（或其子集），列坐标保持不变。
- **模块列绝不允许重新对齐**：模块段继承 Phase 3.8 的模块 MSA（例如 `ACT_msa.afa`）。
- **只有 linker 区域允许自由对齐**：linker 使用 MAFFT E-INS-i（或同等长插入友好算法）在亚型内进行自由对齐。
- 仅对**同一结构架构**的序列构建 full-length MSA（例如 Type Iβ-ACT），避免把“无模块序列”强行缝到“有模块架构”的 linker 上引入幻觉 gap。

### 3.9.1 定义亚型与架构子集

```bash
# 以 Type Iβ + ACT 为例：从模块注释表中筛出同时具备核心与 ACT 的序列
python scripts/select_sequences.py   --presence_table results/03_msa_modules/module_presence_absence_strict.tsv   --require_core 1   --require_module ACT   --output results/03_msa_full/Ib_ACT.ids
```

### 3.9.2 提取 linker 片段（核心尾 ↔ 模块头）

linker 的定义必须可追溯：来源于 **同一条全长序列** 中“核心域 envelope 末端”到“模块域 envelope 起点”的区间。

```bash
mkdir -p results/03_msa_full

# 需要在 Phase 3.6/3.8 的提取脚本中同时输出坐标表：
# results/03_msa_core/core_domain_coords.tsv   （seq_id, core_from, core_to）
# results/03_msa_modules/ACT_domain_coords.tsv （seq_id, act_from,  act_to）

python scripts/extract_linkers.py   --full_length_fasta results/02_qc/nr80_Ib.fasta   --seq_ids results/03_msa_full/Ib_ACT.ids   --core_coords results/03_msa_core/core_domain_coords.tsv   --module_coords results/03_msa_modules/ACT_domain_coords.tsv   --output results/03_msa_full/Ib_ACT_linkers.fasta
```

### 3.9.3 亚型内 linker 自由对齐（允许膨胀，但只发生在 linker）

```bash
mafft --genafpair --maxiterate 1000 --ep 0 --thread 20   results/03_msa_full/Ib_ACT_linkers.fasta   > results/03_msa_full/Ib_ACT_linkers_einsi.afa
```

### 3.9.4 拼接 core + linker + module → full-length MSA

```bash
python scripts/stitch_full_length_msa.py   --seq_ids results/03_msa_full/Ib_ACT.ids   --core_msa results/03_msa_core/core_asr.afa   --linker_msa results/03_msa_full/Ib_ACT_linkers_einsi.afa   --module_name ACT   --module_msa results/03_msa_modules/ACT_msa.afa   --output results/03_msa_core/msa_full_Ib_v4.afa   --emit_column_map results/03_msa_full/msa_full_Ib_column_map.tsv   --assert_core_columns_unchanged
```

**产出：**
- `results/03_msa_core/msa_full_Ib_v4.afa`：用于 Phase 4.4 的嵌套 ASR、Phase 6.2.X 的联合 DCA；
- `results/03_msa_full/msa_full_Ib_column_map.tsv`：full-length 列坐标 →（core/linker/module）分段映射表（后续跨域耦联与结构映射必需）。

**QC2b（强制）：**  
`stitch_full_length_msa.py` 必须断言：
1. full-length MSA 的 core 段列数与 `core_asr.afa` 完全一致；
2. 任取 3 条代表序列，core 段字符串逐列完全相等（允许 gap 不同，但列位置一致）；
3. linker 段膨胀不会“泄漏”进 core 或 module 段（通过 column_map 检查）。

---

---

# 第四步：系统发育与分层 ASR

## 4.1 核心树推断与外群定根（含 LBA 抗性检验）

```bash
# 1) 将 KDOPS 外群核心域对齐后并入核心比对
#    ⚠ 仍需遵循 Phase 3.6 的“先 sto → 剥离 Insert → 再 afa”，避免 outgroup 引入插入列膨胀
hmmalign --trim --mapali results/03_msa_core/skeleton_core_aa.fa -o results/04_phylogeny_asr/kdops_core_raw.sto results/03_msa_core/core_global.hmm data/seeds/kdops_outgroup_core.fasta
esl-alimask --rf-is-mask results/04_phylogeny_asr/kdops_core_raw.sto > results/04_phylogeny_asr/kdops_core_matchonly.sto
esl-reformat afa results/04_phylogeny_asr/kdops_core_matchonly.sto | seqkit seq --upper-case > results/04_phylogeny_asr/kdops_core_aligned.afa

# 2) 合并（核心树用 core_tree.afa；外群仅用于定根）
python scripts/merge_alignments.py --core results/03_msa_core/core_tree.afa --outgroup results/04_phylogeny_asr/kdops_core_aligned.afa --output results/04_phylogeny_asr/core_with_outgroup.afa

# 3) 全局树（baseline：ModelFinder）
iqtree2 -s results/04_phylogeny_asr/core_with_outgroup.afa -m MFP -B 1000 -T AUTO -o "KDOPS_1,KDOPS_2,KDOPS_3" --prefix results/04_phylogeny_asr/CoreTree_rooted_MFP

# 4) 全局树（LBA 抗性：位点异质性/混合模型；至少跑其中之一）
iqtree2 -s results/04_phylogeny_asr/core_with_outgroup.afa -m LG+C20+F+G -B 1000 -T AUTO -o "KDOPS_1,KDOPS_2,KDOPS_3" --prefix results/04_phylogeny_asr/CoreTree_rooted_LGC20

# 可选：再跑一个不同混合模型作交叉验证（更严格）
iqtree2 -s results/04_phylogeny_asr/core_with_outgroup.afa -m EX_EHO+F+G -B 1000 -T AUTO -o "KDOPS_1,KDOPS_2,KDOPS_3" --prefix results/04_phylogeny_asr/CoreTree_rooted_EXEHO
```

**根稳定性检验（QC3，强制）：**

- **外群子集敏感性**：更换 KDOPS 代表集（至少 2 套）重新定根，根位置必须一致或差异可解释（例如某 KDOPS 序列明显异常长枝）。
- **模型敏感性（LBA 诊断）**：`CoreTree_rooted_MFP` 与 `CoreTree_rooted_LGC20/EXEHO` 的三大亚型分化极性（根位置）必须一致；若不一致，必须在结果中声明“根不鲁棒”，并在下游 PastML/ASR 里对该不确定性做敏感性分析（例如在两种根假设下分别重建一次）。
- **中点定根一致性**：对去外群树做 midpoint rooting，检查其与外群定根是否冲突；冲突通常提示 LBA 或外群选择问题。

```bash
python scripts/qc_root_stability.py --tree_mfp results/04_phylogeny_asr/CoreTree_rooted_MFP.treefile --tree_c20 results/04_phylogeny_asr/CoreTree_rooted_LGC20.treefile --tree_exeho results/04_phylogeny_asr/CoreTree_rooted_EXEHO.treefile --outgroup_prefix KDOPS_ --output results/04_phylogeny_asr/QC3_root_stability.md
```

## 4.2 结构系统发育交叉验证——AA 树 vs 3Di 树

FoldMason 的 3Di 比对（`skeleton_3di.fa`）编码三级相互作用模式，在暮光区以下仍能捕捉结构保守性（van Kempen et al., 2024）。当 AA 序列信号衰减时，3Di 提供独立于序列的结构演化证据。

> **模型选择注意事项：** 3Di 虽为 20 状态字母表，但其替换过程的物理含义完全不同于氨基酸替换——3Di 编码三级结构相互作用模式，而非化学侧链性质。直接使用 IQ-TREE 默认 ModelFinder（`-m MFP`）会在氨基酸替换矩阵集合中搜索，存在系统性模型错配风险。Seffernick et al. (2025, *MBE* 42:msaf124) 已专门为 3Di 推导了替换矩阵 **Q.3Di.AF**（基于 AlphaFold 结构）和 **Q.3Di.LLM**（基于 ProtT5 翻译），benchmark 表明在绝大多数蛋白家族中显著优于 GTR20 和所有标准氨基酸模型。

```bash
# 骨架代表的 AA 树（标准氨基酸模型选择）
iqtree2 -s results/03_msa_core/skeleton_refined_aa.fa -m MFP -B 1000 -T AUTO --prefix results/04_phylogeny_asr/SkeletonTree_AA

# 骨架代表的 3Di 树（指定 3Di 专用模型集）
# ⚠ IQ-TREE 2 标准发行版通常不内置 Q.3Di.AF/Q.3Di.LLM，需要用 -mdef 显式加载 Nexus model file
# 将模型文件放在 meta/models/Q_3Di_models.nex，并在 results/meta/model_files.tsv 记录 sha256 与来源
iqtree2 -s results/03_msa_core/skeleton_3di.fa -m MFP -mset Q.3Di.AF,Q.3Di.LLM,GTR20 -mdef meta/models/Q_3Di_models.nex -B 1000 -T AUTO --prefix results/04_phylogeny_asr/SkeletonTree_3Di

# 若仍不支持（版本过旧或模型文件不兼容），使用 GTR20 作为备选：
# iqtree2 -s results/03_msa_core/skeleton_3di.fa -m GTR20 -B 1000 -T AUTO --prefix results/04_phylogeny_asr/SkeletonTree_3Di_GTR20

# 拓扑对比（Robinson–Foulds 距离、quartet 距离、关键节点一致性）
python scripts/compare_trees.py --tree1 results/04_phylogeny_asr/SkeletonTree_AA.treefile --tree2 results/04_phylogeny_asr/SkeletonTree_3Di.treefile --output results/04_phylogeny_asr/tree_comparison.md
```

**解释边界（强制声明）：** 3Di 树仅在骨架代表（有结构的 20–40 序列）上构建，不对全量 12k 序列构建。其角色严格限定为验证 AA 树的关键深分支（三亚型分化节点、KDOPS 外群分离位置）是否鲁棒，而非替代全数据 AA 树。若 AA 树与 3Di 树在某深分支不一致，应解释为"该节点在两种证据下不鲁棒"，而非简单取其中一棵为准。

## 4.3 核心氨基酸 ASR

```bash
iqtree2 -s results/03_msa_core/core_asr.afa \
  -te results/04_phylogeny_asr/CoreTree_rooted.treefile \
  -m MFP -asr \
  --prefix results/04_phylogeny_asr/ASR_core
```

此处 `-te` 合理，因为 `core_asr.afa` 与 `CoreTree_rooted` 基于同一核心比对构建。

## 4.4 嵌套 ASR（局部子树 + 亚型内全长比对）

当需要在亚型内全长比对上进行局部 ASR 以获取完整变构区域信息时：

> **输入约束（强制）：** `results/03_msa_core/msa_full_Ib_v4.afa` 必须来自 **Phase 3.9 轮廓锚定缝合比对**（核心列继承 `core_asr.afa`，模块列继承模块 MSA，仅 linker 自由对齐）。不得用全长序列直接 MAFFT 裸跑生成，否则比对膨胀与 linker 幻觉将回归。

```bash
# 关键修正：用 -g 仅约束拓扑，重新估计枝长和模型参数
iqtree2 -s results/03_msa_core/msa_full_Ib_v4.afa \
  -g results/04_phylogeny_asr/Subtree_Ib_topology.nwk \
  -m MFP -asr \
  --prefix results/04_phylogeny_asr/ASR_Ib_local

# 在子树中保留 3–5 条来自其他亚型的外群序列以锚定根重建
```

## 4.5 Gap 祖先态重建——三方案分层策略

IQ-TREE 将 gap 视为缺失数据，无法重建 indel 事件。这对"变构元件何时插入 TIM barrel"的核心科学问题构成致命盲区。V4.1 采用分层策略，在不同粒度上使用不同方法。

> **方法学澄清：** V4.0 早期版本曾将 Topiary 描述为"集成 DOWNPASS gap 重建 + AltAll"并使用虚构的 `topiary-asr` 命令行接口。经查证，Topiary 的实际 CLI 入口为 `topiary-seed-to-alignment` → `topiary-alignment-to-ancestors` → `topiary-bootstrap-reconcile` 的多步管线。Topiary 确实使用 DOWNPASS（via PastML）为祖先序列分配 gap 态，并生成 gapped 的 ML + altAll 祖先序列（Orlandi et al., 2023, *Protein Sci* 32:e4519）。但 Topiary 的 gap 处理本质上仍是逐列的二态简约法（gap/non-gap），**不是真正的 indel-aware 模型**（不建模插入/缺失事件的长度分布和位置概率）。因此，Topiary 适合回答"此列在该祖先是否存在"，但不适合回答"此 indel 事件何时发生、长度如何变化"。

**方案 A（粗粒度，推荐用于模块级问题）：模块存在/缺失 → PastML 离散性状 ASR**

这是 Phase 4.6 的内容，天然适合"模块何时获得"的科学问题——将整个 ACT 域、N 端延伸等视为单一二态性状（有/无），在核心树上用 PastML 的 MPPA 或 MAP 方法重建。维度低（每个模块 1 列），无多重比较问题。

**方案 B（中粒度，核心列级 gap）：gap 二值化 → 分区段聚合 → PastML**

对核心 MSA 的逐列 gap 直接做 PastML 会面临维度爆炸（数百列 × 数千 tip）。正确做法是先将连续 gap 列聚合为"indel 区段（gap block）"，每个区段作为一个二态性状。

```bash
# Step 1：将连续 gap 列聚合为 indel 区段
python scripts/aggregate_gap_blocks.py \
  --input results/03_msa_core/core_asr.afa \
  --min_block_size 3 \
  --output results/04_phylogeny_asr/gap_blocks.tsv
  # 输出格式：seq_id | block_1(cols 45-52) | block_2(cols 120-135) | ...
  # 值：1 = 该区段全为 gap, 0 = 该区段有非 gap 残基

# Step 2：在核心树上对每个区段做离散性状 ASR
pastml --tree results/04_phylogeny_asr/CoreTree_rooted.treefile \
  --data results/04_phylogeny_asr/gap_blocks.tsv \
  --out_data results/04_phylogeny_asr/gap_block_asr/
```

**方案 C（细粒度，真正的 indel-aware ASR，可选进阶）：indelMaP 或 ARPIP**

若需要为关键祖先节点重建精确的 indel 事件（插入 vs 缺失区分、indel 长度分布），应使用具有显式 indel 模型的方法：

- **indelMaP**（Jia et al., 2024, *MBE* 41:msae109）：indel-aware 简约法，支持长 indel，Rust 实现高效。
- **ARPIP**（Jowkar et al., 2023, *Syst Biol*）：基于 Poisson Indel Process 的似然法 ASR，可建模单残基 indel 事件。
- **Historian**（Holmes, 2017）：时间依赖进化模型的 indel-aware MSA+ASR。

```bash
# indelMaP 示例（需单独安装）
indelmap \
  --msa results/03_msa_core/core_asr.afa \
  --tree results/04_phylogeny_asr/CoreTree_rooted.treefile \
  --output results/04_phylogeny_asr/indelmap_ancestors/
```

**方案选择决策树：**
- "模块何时获得" → 方案 A（Phase 4.6）
- "核心 TIM barrel 的哪些 loop 在哪个节点发生长度变化" → 方案 B
- "精确区分特定 indel 是插入还是缺失、重建 indel 历史" → 方案 C

## 4.6 模块获得/丢失的离散性状 ASR

```bash
# 使用严格版和宽松版模块矩阵分别运行（敏感性分析，见 Phase 3.8）
for version in strict relaxed; do
  pastml --tree results/04_phylogeny_asr/CoreTree_rooted.treefile \
    --data results/03_msa_modules/module_presence_absence_${version}.tsv \
    --out_data results/04_phylogeny_asr/module_origins_${version}/ \
    --html_compressed results/04_phylogeny_asr/module_origins_${version}_viz.html
done
# 比较两版结果：起源节点是否一致 → 写入 QC3
```

输出：每一模块的"获得/丢失"最可能节点集合。

## 4.7 AltAll 系综采样

```bash
# AltAll 系综采样：使用 maxaltall R 包处理 IQ-TREE .state 文件
# 阈值：PP₁ = 0.80, PP₂ = 0.20 (Eick et al., 2017; Muñiz-Trejo et al., 2025)
```

## 4.8 模块层局部 ASR

对每个模块：提取携带子集 → 模块片段 MSA → 核心树对应子树拓扑 → 局部 ASR。

```bash
for module in ACT N_ext alpha2beta3; do
  iqtree2 -s results/03_msa_modules/${module}_msa.afa \
    -g results/04_phylogeny_asr/Subtree_${module}_topology.nwk \
    -m MFP -asr \
    --prefix results/04_phylogeny_asr/ASR_${module}
done
```

### QC3（强制产出）

`results/04_phylogeny_asr/qc_phylogeny_asr.md`：
- 关键节点（三亚型分化、模块获得前/后）的 bootstrap 支持
- 外群定根稳定性（不同 KDOPS 子集）
- AA 树 vs 3Di 树关键节点一致性
- 祖先序列的平均 PP 分布
- 模块获得节点在不同定根方案下是否稳定

---

# 第五步：结构重建与多聚体动力学实验

## 5.1 候选祖先集合定义

至少包含三类祖先：

1. **Pre-gain（模块获得前）：** 核心树上模块性状从 0→1 的前一节点附近
2. **Post-gain（模块获得后）：** 模块性状为 1 的最早节点附近
3. **对照祖先：** 同深度但无模块转换事件的节点，排除时间深度混淆效应

> **⚠ [CHECK-02] 祖先序列来源：** 核心层 ASR + 模块离散性状 ASR 的作用是**锁定哪些节点值得深入分析**（起源候选集）。真正送进 AlphaFold3 的全长祖先序列必须来自 **Phase 4.4 的亚型内全长嵌套 ASR**——该比对保留了 linker 区域，重建出的祖先是一条连续、演化自洽的全长多肽。**绝不可**将核心 ASR 与模块 ASR 的片段手工硬拼接，否则 linker 处序列无演化依据，AF3 预测将在连接处严重解折叠。

## 5.2 结构预测（多聚体）

DAH7PS 各亚型的生理组装态均为同源四聚体（Cross et al., 2013; Lang et al., 2016）。变构位点位于亚基界面（Webby et al., 2010）。

```bash
# AlphaFold3 预测：输入 4 拷贝祖先序列以预测四聚体
# 对关键祖先同时使用 ESMFold 交叉验证

# 用 FoldMason 标记预测结构的可信度
foldmason easy-msa ancestor_af3.pdb known_pdbs/*.pdb \
  results/05_struct_md/struct_comparison tmp_struct --report-mode 1
# 低 LDDT 区域 = 预测不可靠 → MD 分析时降权
```
## 5.2a 祖先配体态预测的 Apo-first 门控（强制）⚠ [CHECK-07]

> **为什么必须门控：** AF3 的复合物预测在“给定配体”的条件下可能发生强烈的诱导契合偏置。对于尚未形成成熟别构口袋的古老祖先，直接喂入 Trp/Phe 等现代效应物，容易得到“被模型硬捏出来”的假口袋与假构象，进而污染 MD。

**门控策略（Apo-first + 物理验证）：**

1. **先做 Apo**：每个祖先节点先预测 apo 四聚体结构（不提供任何配体）。
2. **再评估口袋**：用口袋检测（Fpocket 等）在 apo 结构上寻找潜在结合位点，并与现代已知效应物口袋位置作几何对应（同源残基/界面位置）。
3. **再做物理对接**：仅在 apo 口袋可达的情况下，对 Trp/Phe 做分子对接（AutoDock Vina/Glide 等），检查是否存在合理的结合姿势与能量。
4. **通过才允许 Holo**：只有在“口袋存在 + 对接合理 + 位点同源合理”的情况下，才进入 holo AF3 预测与 holo-MD；否则该祖先节点 **不做 holo 条件**，并在论文中声明“该节点尚不具备成熟配体容纳能力”。

示例命令（可替换为本地既有流程）：

```bash
# 口袋检测（Fpocket）
fpocket -f results/05_struct_md/af3/apo/<node>/model.pdb -o results/05_struct_md/pocket/<node>/

# 对接（AutoDock Vina；示例）
vina --receptor results/05_struct_md/af3/apo/<node>/model.pdb      --ligand data/ligands/trp.pdbqt      --center_x <x> --center_y <y> --center_z <z>      --size_x 22 --size_y 22 --size_z 22      --out results/05_struct_md/docking/<node>/trp_best.pdbqt      --log results/05_struct_md/docking/<node>/trp_vina.log
```

**产出（强制）：**
- `results/05_struct_md/pocket_gating/<node>_gating.md`：记录 apo 口袋、对接结果、是否允许 holo 的结论与依据；
- 对未通过门控的节点：在后续 5.3 的实验矩阵中跳过 holo 条件（只跑 apo）。



## 5.2b MD 前结构 QC 门控（强制）⚠ [CHECK-05]

AltAll 系综的位点独立性假设可能将互相排斥的大侧链塞进同一亚基界面，导致 AF3 预测出高度扭曲的四聚体。在提交正式 MD 之前必须通过以下门控：

```
1. AF3 置信度筛选：
   - ipTM ≥ 0.6（界面置信度）
   - 跨链 PAE < 15 Å（界面残基对）
   - 不通过 → 重新检查 AltAll 替换位点是否位于界面，必要时换用 ML-best 序列

2. 能量最小化与分阶平衡（GROMACS）：
   - 严苛能量最小化（steep + cg，Fmax < 500 kJ/mol/nm）
   - NVT 平衡（100 ps，位置限制 fc=1000 kJ/mol/nm²）
   - NPT 平衡（200 ps，位置限制 fc=1000）
   - NPT 平衡（200 ps，位置限制 fc=100）
   - NPT 平衡（200 ps，无位置限制）

3. 短轨迹筛选（10 ns）：
   - 四聚体 RMSD 是否稳定（无爆炸性漂移）
   - 界面面积是否维持（SASA 变化 < 20%）
   - 通过 → 升级为 ≥500 ns 正式轨迹
   - 不通过 → 终止该祖先的 MD，标记为"结构不可信"
```

## 5.3 多聚体 MD：2×2 因子实验

```
因子 1（组装态）：单体 vs 四聚体
因子 2（配体态）：apo vs effector（Trp+Phe 或对应亚型的调控配体）

实验矩阵（每个祖先节点；**holo 条件需先通过 ⚠ [CHECK-07] 门控**）：
  四聚体-apo:       ≥500 ns × 3 replicas（必须，主分析）
  四聚体-effector:  ≥500 ns × 3 replicas（仅当门控允许 holo）
  单体-apo:          200 ns × 2 replicas（对照）
  单体-effector:     200 ns × 2 replicas（可选；同样需门控允许 holo）

若某祖先节点未通过 [CHECK-07]（口袋/对接不支持配体容纳），则该节点只跑 apo（四聚体-apo + 单体-apo），并在结果中明确标注“无 holo 条件”。
```

单体对照的科学价值：证明"变构信号依赖界面"——若单体 DCCM 中活性位点-变构位点通讯消失而四聚体中存在，则界面对变构通讯至关重要。

## 5.4 动力学读出

- RMSD/RMSF、二级结构稳定性
- DCCM（动态互相关矩阵）
- 残基相互作用网络、信息流/最短路径（核心→模块、跨亚基界面）
- **降权策略：** 结构预测低置信区段（FoldMason LDDT 标注）不主导网络分析结论

### QC4（强制产出）

`results/05_struct_md/qc_dynamics.md`：
- 轨迹收敛性（RMSD 平台期、多重复一致性）
- 关键口袋/界面面积稳定
- 低置信区段标注与处理方式记录

---

# 第六步：ICDC——分层 DCA × 动力学通讯网络整合

## 6.1 DCA 输入准备（核心/模块分离 + 亚型内联合）

> **核心规则：** DCA 的可解释性高度依赖有效信息量。此 SOP 将 **Meff/L ≥ 3.0** 设为“可进入主分析/论文证据链”的最低门槛（理想 ≥ 5.0）。低于门槛的 DCA 结果一律视为探索性噪声，不得写入 ICDC 主网络。

```bash
mkdir -p results/06_icdc

# 6.1.1 核心层 DCA 输入：去极端 gap 列 + 去极端 gap 行 + 序列降权
python scripts/prepare_dca_input.py   --input results/03_msa_core/core_asr.afa   --gap_col_max 0.50   --gap_row_max 0.50   --meff_min 3.0   --output results/06_icdc/core_dca.afa   --stats results/06_icdc/core_dca_stats.tsv
  # 核心通常应达到 Meff/L > 5；若 <3 则停止，回溯到 MSA/QC2

# 6.1.2 模块层 DCA 输入（仅在携带模块的子集内做）
# ⚠ [CHECK-03] Meff/L < 3.0 → 自动跳过该模块 DCA（禁止把噪声耦联写入网络）
for module in ACT N_ext; do
  python scripts/prepare_dca_input.py     --input results/03_msa_modules/${module}_msa.afa     --gap_col_max 0.60     --gap_row_max 0.50     --meff_min 3.0     --skip_if_insufficient     --output results/06_icdc/${module}_dca.afa     --stats results/06_icdc/${module}_dca_stats.tsv
done

# 6.1.3 亚型内联合 DCA 输入（用于跨域 core↔module 耦联；例如 Type Iβ-ACT）
# 联合 DCA 仅对“同一架构且同时具备 core+module”的序列子集运行（避免缺失模块造成无意义的 gap）
python scripts/subset_msa_by_ids.py   --msa results/03_msa_core/msa_full_Ib_v4.afa   --ids results/03_msa_full/Ib_ACT.ids   --output results/06_icdc/Ib_ACT_full_for_dca.afa

python scripts/prepare_dca_input.py   --input results/06_icdc/Ib_ACT_full_for_dca.afa   --gap_col_max 0.60   --gap_row_max 0.50   --meff_min 3.0   --skip_if_insufficient   --output results/06_icdc/Ib_ACT_joint_dca.afa   --stats results/06_icdc/Ib_ACT_joint_dca_stats.tsv
```

## 6.2 DCA 执行（完整参数化规范）

### 6.2.1 核心层 DCA

使用 plmc（EVcouplings 后端），参数遵循 Ekeberg et al. (2013/2014) 与 Hopf et al. (2017) 的推荐值：

```bash
plmc -o results/06_icdc/core.params   -c results/06_icdc/core_couplings.txt   -le 16.0 -lh 0.01 -m 100 -t 0.8   -f <focus_seq_id> -g   results/06_icdc/core_dca.afa
```

### 6.2.2 模块层 DCA（仅在 Phase 6.1 通过 Meff/L 门槛的模块上执行）

```bash
for module in ACT N_ext; do
  if [ -f results/06_icdc/${module}_dca.afa ]; then
    plmc -o results/06_icdc/${module}.params       -c results/06_icdc/${module}_couplings.txt       -le 16.0 -lh 0.01 -m 100 -t 0.8       -f <module_focus_seq> -g       results/06_icdc/${module}_dca.afa
  fi
done
```

### 6.2.3 亚型内跨域/联合 DCA（Subtype-specific Joint DCA）

> **核心修补点：** 若核心与模块被拆成两个独立 MSA 分别跑 DCA，将永远无法得到 core 残基 *i* 与模块残基 *j* 的直接耦联。  
> 联合 DCA 在**同一拼接比对矩阵**上估计 MRF，可直接读出跨越“core–linker–module”边界的耦联对，是证明模块被物理整合进变构网络的关键演化证据。

```bash
# 以 Type Iβ-ACT 为例：联合 DCA
if [ -f results/06_icdc/Ib_ACT_joint_dca.afa ]; then
  plmc -o results/06_icdc/Ib_ACT_joint.params     -c results/06_icdc/Ib_ACT_joint_couplings.txt     -le 16.0 -lh 0.01 -m 100 -t 0.8     -f <Ib_ACT_focus_seq> -g     results/06_icdc/Ib_ACT_joint_dca.afa

  # 提取“跨域耦联”（i 在 core 段，j 在 ACT 段；linker 可选）
  python scripts/extract_crossdomain_couplings.py     --couplings results/06_icdc/Ib_ACT_joint_couplings.txt     --column_map results/03_msa_full/msa_full_Ib_column_map.tsv     --segment_a core     --segment_b ACT     --top_k 200     --output results/06_icdc/Ib_ACT_crossdomain_top200.tsv
fi
```

### 6.2.4 显著性评估与接触判据

```bash
# 接触定义：Cβ–Cβ 距离 < 8 Å（glycine 用 Cα）
# 短程排除：|i - j| > 5（排除序列近邻的非特异信号）
python scripts/dca_significance.py   --couplings results/06_icdc/core_couplings.txt   --method top_L   --contact_pdb data/structures/panel_dah7ps/1KFL.pdb   --contact_cutoff 8.0   --min_separation 5   --output results/06_icdc/core_significant_couplings.tsv
```

### 6.2.5 “模块获得前后耦联变化”的比较策略（Meff 匹配下采样 + Z-score）

> 仅用 `core_with_ACT` vs `core_without_ACT` 的原始 Rank 差异会遭遇系统发育背景与 Meff 混杂（Simpson’s paradox）。  
> **V4.1 强制：** 比较必须在 **Meff 匹配** 的信息量下进行，并输出效应量（Z-score）与稳定性（bootstrap/重复下采样 p）。

```bash
# 1) 拆分两组
python scripts/split_msa_by_module.py   --msa results/06_icdc/core_dca.afa   --module_table results/03_msa_modules/module_presence_absence_strict.tsv   --module ACT   --output_with results/06_icdc/core_with_ACT.afa   --output_without results/06_icdc/core_without_ACT.afa

# 2) 计算两组 Meff，并设定 target = min(Meff_with, Meff_without)
python scripts/compute_meff.py   --msa results/06_icdc/core_with_ACT.afa   --reweight_threshold 0.8   --out results/06_icdc/meff_core_with_ACT.tsv
python scripts/compute_meff.py   --msa results/06_icdc/core_without_ACT.afa   --reweight_threshold 0.8   --out results/06_icdc/meff_core_without_ACT.tsv

# 3) 对 Meff 较大的那一组做多次下采样，匹配到 target Meff（建议 20 次）
python scripts/downsample_to_meff.py   --msa_high results/06_icdc/core_with_ACT.afa   --msa_low  results/06_icdc/core_without_ACT.afa   --target_strategy min_group   --reps 20   --out_dir results/06_icdc/downsampled_meff_matched/

# 4) 对每个下采样 replicate 运行 DCA，并汇总得到每个位点对的耦联分布
python scripts/run_plmc_batch.py   --msa_dir results/06_icdc/downsampled_meff_matched/   --le 16.0 --lh 0.01 --m 100 --t 0.8 --g   --focus <focus_seq_id>   --out_dir results/06_icdc/downsampled_meff_matched/plmc_runs/

python scripts/dca_compare_meff_matched.py   --runs_dir results/06_icdc/downsampled_meff_matched/plmc_runs/   --group_labels with_ACT without_ACT   --report_zscore   --report_rank_shift   --output results/06_icdc/coupling_shift_ACT_meff_matched.tsv
```

**输出解释（强制）：**
- `ΔRank`：只能作为“排序变化”的直观指标；
- `Z-score(ΔCoupling)`：在等效信息量下的效应量（主指标）；
- `p_bootstrap / stability`：在重复下采样下的可重复性（过滤虚假耦联）。

## 6.3 ICDC 融合

> **⚠ [CHECK-04] 坐标系映射验证（强制前置）：** ICDC 融合涉及五套坐标系（全长序列编号、核心 MSA 列号、模块 MSA 列号、PDB 残基编号、GROMACS 拓扑编号）。在执行融合之前，必须运行 `coordinate_mapper.py`，以 1KFL（E. coli AroG）为绝对锚点，断言其催化残基（K97、R165 等）在所有五套坐标系中的映射完全等价。任何断言失败则终止管线。

```bash
python scripts/coordinate_mapper.py \
  --anchor_pdb data/structures/panel_dah7ps/1KFL.pdb \
  --core_msa results/03_msa_core/core_asr.afa \
  --module_msa results/03_msa_modules/ACT_msa.afa \
  --gromacs_top results/05_struct_md/topol.top \
  --known_sites "K97,R165,H268,C61,E302" \
  --output results/06_icdc/coordinate_map.tsv \
  --assert_consistent
# 断言失败 → 终止；通过 → 生成统一坐标映射表供下游使用
```

ICDC 边集定义为满足以下条件的残基对/路径：
1. 具有显著 DCA 耦联（核心层 / 模块层 / **亚型内联合 DCA 的跨域 core↔module 耦联**）
2. 在 MD 中呈现稳定的动力学相关或信息流贡献
3. 在模块获得前/后发生系统性变化（在 Meff 匹配比较下仍成立；效应大小与方向可解释）

```bash
# 输出
# results/06_icdc/icdc_core_network.graphml
# results/06_icdc/icdc_module_network.graphml
# results/06_icdc/icdc_crossdomain_couplings.tsv（联合 DCA：core↔module 耦联）
# results/06_icdc/icdc_crosslayer_paths.tsv（核心↔模块↔界面）
```

**科学叙事整合：** 以"模块获得事件"为时间轴，在核心树上投影网络指标变化，而非把全长硬映射到同一 MSA 坐标。

---

# 第七步：工程风险 Checklist

> 核心/模块分治范式在理论上清晰，但在实操中存在七处隐蔽的工程暗礁。以下 Checklist 项已内嵌于对应 Phase 的具体命令中（以 `⚠ [CHECK-xx]` 标注），此处集中列出供执行时核查。

### [CHECK-01] HMM 域提取的"边缘破碎"效应

| | |
|---|---|
| **关联 Phase** | 3.6（`extract_core_domains.py`） |
| **风险** | `hmmsearch` 的 `env_from/env_to` 在远缘序列上边界呈锯齿状。硬切可能截断 TIM-barrel 第一条 β1 链或最后一条 α8 螺旋，导致 ASR 在这些位置重建为 gap 而非真实氨基酸。 |
| **缓解** | 向 `env_from` 和 `env_to` 两端各外扩 20 aa 缓冲带（`--pad 20`）。多余的非同源残基由下游 `hmmalign --trim` 自动软修剪（标记为小写/插入态），不进入 match columns。 |
| **失败处理** | 若缓冲带后仍有 >5% 序列的 β1/α8 被截断（可通过检查 hmmalign 输出中首尾 match column 的 gap 率判断），则回溯检查 HMM 模型质量或补充种子多样性。 |
| **参数记录** | `meta/params.json → core_definition.pad_residues: 20` |

### [CHECK-02] 祖先序列的"弗兰肯斯坦缝合"

| | |
|---|---|
| **关联 Phase** | 5.1（候选祖先集合定义）、5.2（结构预测） |
| **风险** | 核心 ASR 与模块 ASR 分别在不同比对/列坐标系上执行。手工硬拼接两者的输出序列时，linker 区域无任何 ASR 覆盖——该片段没有演化依据，AF3 预测将在连接处严重解折叠。 |
| **缓解** | 核心层 ASR + 模块离散性状 ASR 仅用于**锁定起源候选节点**。真正送进 AF3 的全长祖先序列必须来自 **Phase 4.4 的亚型内全长嵌套 ASR**，该比对保留了 linker，重建出的祖先是一条连续、演化自洽的全长多肽。 |
| **失败处理** | 若 Phase 4.4 的全长嵌套 ASR 在目标节点的 linker 区域平均 PP 极低（<0.3），说明 linker 演化信号不足；此时应报告该节点的全长祖先"linker 不可信"，在 MD 分析中对 linker 区域降权或截断分析。 |

### [CHECK-03] DCA 的“有效序列枯竭”与欠定陷阱（Meff/L 门控）

| | |
|---|---|
| **关联 Phase** | 6.1（DCA 输入准备）、6.2（DCA 执行） |
| **风险** | 单独切出的模块（如 ACT 域）或“亚型内联合比对”（core+linker+module），在去冗余与高 gap 过滤后，有效序列深度 **M_eff** 可能显著下降。若 **M_eff/L < 3.0**，plmc 的伪似然推断高度欠定，耦联矩阵会被正则化先验与采样噪声主导，假阳性泛滥。 |
| **缓解** | 在 `prepare_dca_input.py` 中强制打印并写入 `*_dca_stats.tsv`：`L_eff`、`M_eff`、`M_eff/L`。当 **M_eff/L < 3.0** 时必须自动跳过该 DCA（`--meff_min 3.0 --skip_if_insufficient`）。理想门槛为 **≥ 5.0**。 |
| **失败处理** | **禁止**让低深度 DCA 进入 ICDC 主网络：此类模块/联合 DCA 仅保留 MD/DCCM/结构证据，并在 ICDC 报告中显式标注“无可靠 DCA 支撑”。如必须给出演化证据，可改用更低维的离散性状 ASR（PastML）或扩大采样补深度。 |
| **参数记录** | `meta/params.json → dca.meff_min_main: 3.0; dca.meff_ideal: 5.0` |

### [CHECK-04] 五套坐标系的"地狱级映射"

| | |
|---|---|
| **关联 Phase** | 6.3（ICDC 融合） |
| **风险** | ICDC 融合同时涉及五套完全不同的坐标系：全长序列编号、核心 MSA 列号、模块 MSA 列号、PDB 物理残基编号、GROMACS 拓扑编号。任意一对映射错位 1 个残基，DCA 耦联与 DCCM 动力学图谱将永远无法缝合。 |
| **缓解** | Phase 6 执行前，必须运行 `coordinate_mapper.py`。以 1KFL（E. coli AroG）为绝对锚点，断言其催化残基（K97、R165 等已知位点）在所有五套坐标系中的映射完全等价。 |
| **失败处理** | 任何断言失败 → 立即终止管线，回溯定位映射断裂点。常见原因：(i) PDB 残基跳号/insertion code 未处理；(ii) GROMACS `pdb2gmx` 丢弃了非标准残基导致编号偏移；(iii) ClipKIT 修剪改变了 MSA 列号但未更新映射表。 |

### [CHECK-05] AlphaFold3 祖先多聚体的"立体位阻爆炸"

| | |
|---|---|
| **关联 Phase** | 5.2b（MD 前结构 QC 门控） |
| **风险** | ASR 的位点独立性假设可能将两个互相排斥的大侧链强行塞进同一亚基界面（尤其在 AltAll 替换中）。AF3 预测出的四聚体内部张力极大，MD 在 0.1 ns 内即可"爆炸"。 |
| **缓解** | 三级门控：(1) AF3 ipTM ≥ 0.6 + 跨链 PAE < 15 Å；(2) GROMACS 严苛能量最小化 + 分阶位置限制平衡（NVT→NPT→逐步撤约束）；(3) 10 ns 短轨迹筛选——四聚体 RMSD 稳定 + 界面面积维持。 |
| **失败处理** | 门控不通过 → 不执行正式 MD。首先检查 AltAll 替换位点是否集中在界面；若是，换用 ML-best 序列重新预测。若 ML-best 也不通过，标记该节点为"结构不可信"并在论文中报告。 |

### [CHECK-06] 内部插入导致的 HMM Hit 断裂：必须做 Hit Stitching

| | |
|---|---|
| **关联 Phase** | 3.6（`extract_core_domains.py`） |
| **风险** | 对于 Type II 的 `α2β3` 内插片（硬插在 TIM-barrel 内部），`hmmsearch` 常把核心域识别为**两个碎片化 hits**（TIM 前半段 + TIM 后半段）。若脚本只取 E-value 最优的单一 hit，将直接截掉半个催化核心，导致后续 MSA/ASR/DCA 全线失真。 |
| **缓解** | `extract_core_domains.py` 必须实现 **hit stitching**：对同一序列上、同一 HMM 的多段 hits，若满足“顺序合法、方向一致、间隔区段为潜在插入、合并后覆盖率显著提高”，则合并为单一 envelope（`from=min(env_from)`，`to=max(env_to)`）。同时输出合并日志（合并前后覆盖率、gap 长度、evalue）。 |
| **失败处理** | 若出现“大量序列需要 stitching”且合并间隔极不一致，提示核心 HMM 可能过窄或种子多样性不足：回溯 Phase 1 的种子集，补充 Type II 多样性并重建 HMM。 |
| **参数记录** | `meta/params.json → core_definition.hit_stitching: "auto"` |

### [CHECK-07] AlphaFold3 祖先预测的“现代配体幻觉”：Apo-first + 物理门控

| | |
|---|---|
| **关联 Phase** | 5.2（结构预测）、5.3（MD 因子实验） |
| **风险** | 若对“尚未形成成熟变构口袋”的祖先节点强行喂入现代效应物（Trp/Phe），AF3 可能通过诱导契合偏置捏造结合位点（侧链非物理性重排），从而把 MD 起始构象推入非自然态，污染下游动力学差异。 |
| **缓解** | **Apo-first**：祖先节点默认只做 Apo 预测。随后做口袋检测（Fpocket 等）+ 物理对接（AutoDock Vina/Glide 等）。仅当“口袋可达 + 对接得分合理 + 位点与现代结构同源/合理”时，才允许进入 Holo 预测与 Holo-MD；否则跳过该节点的 Holo 条件，并在结果中明确宣告“该节点不具备成熟配体容纳能力”。 |
| **失败处理** | 若 Apo 预测本身不稳定（低 ipTM/高 PAE/严重 clash），先回退到 ML-best 序列或减少 AltAll 替换；仍不稳定则标记为“结构不可信节点”，不进入 MD 主分析。 |
| **参数记录** | `meta/params.json → af3.apo_first: true; af3.holo_requires_pocket: true` |



---

# 里程碑与交付物

| 里程碑 | 交付物 | 验收标准 |
|--------|--------|---------|
| M1 | 结构骨架 + 核心列定义 | LDDT 报告 + 比对长度合理 |
| M2 | 全量 core_global MSA | 膨胀比 <3:1 + 结构复核通过 |
| M3 | 外群定根核心树 | 根位置对外群抽样鲁棒 |
| M4 | 模块性状 ASR | 候选起源节点可复现 |
| M5 | 关键祖先 2×2 MD | 动力学差异与模块起源节点一致 |
| M6 | 分层 ICDC | 演化耦联与动力学通讯闭环证据 |

---

# 参考文献

1. Cross PJ, Dobson RCJ, Patchett ML, Parker EJ. The diversity of allosteric controls at the gateway to aromatic amino acid biosynthesis. *Biochem J*. 2013;451(2):127–136.
2. Eick GN, Bridgham JT, Anderson DP, Harms MJ, Thornton JW. Robustness of reconstructed ancestral protein functions to statistical uncertainty. *Mol Biol Evol*. 2017;34(2):247–261.
3. Gilchrist CLM, Mirdita M, Steinegger M. Multiple protein structure alignment at scale with FoldMason. *Science*. 2026; doi:10.1126/science.ads6733.
4. Hill MD, Shafiei R, Muñoz-Rojas T, Harms MJ. Topiary: pruning the manual labor from ancestral sequence reconstruction. *Protein Sci*. 2023;32(1):e4519.
5. Katoh K, Kuma K, Toh H, Miyata T. MAFFT version 5: improvement in accuracy of multiple sequence alignment. *Nucleic Acids Res*. 2005;33(2):511–518.
6. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7. *Mol Biol Evol*. 2013;30(4):772–780.
7. Lang EJM, Cross PJ, Sherwood OC, Parker EJ. Quaternary structure is an essential component that contributes to the sophisticated allosteric regulation mechanism in a key enzyme from *Mycobacterium tuberculosis*. *PLoS ONE*. 2016;12(1):e0169601.
8. Mariani V, Biasini M, Barbato A, Schwede T. lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests. *Bioinformatics*. 2013;29(21):2722–2728.
9. Minh BQ, Schmidt HA, Chernomor O, et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol*. 2020;37(5):1530–1534.
10. Muñiz-Trejo R, Park Y, Thornton JW. Robustness of ancestral sequence reconstruction to among-site and among-lineage evolutionary heterogeneity. *Mol Biol Evol*. 2025;42(4):msaf084.
11. Steinegger M, Söding J. Clustering huge protein sequence sets in linear time. *Nat Commun*. 2018;9:2542.
12. Steenwyk JL, Buida TJ III, Li Y, Shen XX, Rokas A. ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. *PLoS Biol*. 2020;18(12):e3001007.
13. van Kempen M, Kim SS, Tumescheit C, et al. Fast and accurate protein structure search with Foldseek. *Nat Biotechnol*. 2024;42:243–246.
14. Webby CJ, Baker HM, Lott JS, Baker EN, Parker EJ. Synergistic allostery, a sophisticated regulatory network for the control of aromatic amino acid biosynthesis in *Mycobacterium tuberculosis*. *J Biol Chem*. 2010;285(40):30567–30576.
15. Ishikawa SA, Zhukova A, Iwasaki W, Gascuel O. A fast likelihood method to reconstruct and visualize ancestral scenarios. *Mol Biol Evol*. 2019;36(9):2069–2085.
16. Jiao W, Lang EJM, Bai Y, Fan Y, Parker EJ. Diverse allosteric componentry and mechanisms control entry into aromatic amino acid biosynthesis. *Curr Opin Struct Biol*. 2020;65:159–167.
17. Seffernick JT, Pardo-Avila F, Shen H, et al. A general substitution matrix for structural phylogenetics. *Mol Biol Evol*. 2025;42(6):msaf124.
18. Ekeberg M, Lövkvist C, Lan Y, Weigt M, Aurell E. Improved contact prediction in proteins: using pseudolikelihoods to infer Potts models. *Phys Rev E*. 2013;87(1):012707.
19. Hopf TA, Green AG, Schuber B, et al. The EVcouplings Python framework for coevolutionary sequence analysis. *Bioinformatics*. 2019;35(9):1582–1584.
