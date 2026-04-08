# AGENTS.md — DAH7PS 变构起源演化动力学重建（V6 SOP rev7 执行指南）

> **唯一真相源（Single Source of Truth）**：`PLAN.md`（V6 SOP rev7）  
> **本文件的作用**：把 `PLAN.md` 翻译成可执行的工程契约、阶段门控、验收标准与禁止事项，确保自动化执行时不偏航、不偷换问题、不提前下结论。

---

## 0. 文档层级与当前状态

### 0.1 文档层级（必须遵守）

1. **`PLAN.md`**：研究策略、阶段定义、科学边界、门控标准的最高依据。
2. **`results/meta/metrics_manifest.tsv`**：所有数字、行列数、样本量、阈值状态的唯一指标真源。
3. **`results/meta/progress_snapshot.md`**：当前推进状态的简明快照。
4. **`AGENTS.md`**：执行契约与验收清单。
5. **`TASKS.md`**：动态任务清单与完成状态。
6. **`log.md`**：不可追认修改的执行记录；每次运行必须即时记录。

### 0.2 当前状态快照（写给执行代理）

截至本版：

- Phase 0–3.8：**已完成**。
- **Phase 3.9：已完成**，不得再写成待执行。
- Phase 4.1：MFP rooted tree 已完成（集群）；**LG+C20 树未完成**（本地 OOM，需在集群补跑）。
- Phase 4.3：**tip-set consistency 已解决**；当前 pruned ingroup tree 与 `core_asr.afa` 已严格匹配。
- 当前主风险：**root robustness**，不是 core MSA 崩坏，也不是 tip mismatch。
- 在 QC3 完成前：
  - 可以继续跑 working-tree 基础上的 ASR、trait 准备、DCA 输入准备；
  - **不得**把 deepest root 当作唯一真相；
  - **不得**锁定 Phase 5 主线祖先节点。

### 0.3 执行代理必须理解的核心变化

V6 的关键不是“返工前半程”，而是：

- 用 **多 root scenario** 管理不确定性；
- 把 **QC3** 提升为真正 gate；
- 把 **4.2 AA vs 3Di tree** 提升为高优先级正交证据；
- 把 **Phase 5 节点选择标准** 收紧；
- 把数字与状态集中到 `metrics_manifest.tsv`，避免 README / PLAN / TASKS / log 漂移。

---

## 1. 总原则（硬约束）

1. **不得跳过 Phase gate。**  
   除非 `PLAN.md` 明确写为 optional / exploratory，否则不得任意重排会改变解释方向的步骤。

2. **QC 失败即停止。**  
   任意断言失败、tip mismatch、column drift、Meff/L 不达标、装配体判定缺证据、root scenario 未完成时，禁止带病推进到会放大错误解释的下一阶段。

3. **工作树不等于叙事树。**  
   当前 rooted ingroup tree 可作为 **working tree** 执行 ASR 与下游准备；但在 QC3 通过前，禁止把它当作论文中 deepest history 的唯一依据。

4. **所有数值必须写入单一真源。**  
   任何脚本运行后，如产生新的最终数字、样本量、行列数、阈值状态，必须同步更新：
   - `results/meta/metrics_manifest.tsv`
   - `results/meta/progress_snapshot.md`
   - `log.md`

5. **禁止覆盖结果。**  
   正式运行不得直接覆盖旧产物。必须使用：
   - 时间戳子目录，或
   - 旧文件改名 `.bak` 后再写新文件。

6. **参数、阈值、模型必须可追溯。**  
   统一写入 `meta/params.json`；随机性流程必须记录 seed 与选中 ID 列表。

7. **核心 MSA 严禁 insert 膨胀。**  
   只能走：
   `Stockholm -> RF mask -> esl-alimask --rf-is-mask -> AFA`  
   禁止用 `hmmalign --outformat afa` 直接输出核心 AFA。

8. **全长祖先序列禁止手工拼接。**  
   发给 AF3 的全长祖先序列只能来自 **Phase 4.4 subtype-local nested ASR** 正式输出；禁止人工把 core ASR 片段和模块片段拼成“祖先全长序列”。

9. **DCA 门槛是硬门槛。**  
   - 主线仅允许 **core DCA**。
   - 主线最低要求：`Meff/L >= 3.0`；理想目标：`>= 5.0`。
   - ACT / module / joint DCA 因深度不足仅属 exploratory，不得上升为论文主证据。

10. **装配体状态不得默认。**  
    不得把任何祖先节点自动视为 tetramer。每个候选祖先必须先经过：
    - 文献与注释扫描
    - dimer / tetramer 平行 AF3
    - PISA / 界面证据综合判定  
      产出 `assembly_adjudication.tsv` 后，后续 AF3 / MD 才可读取其拷贝数。

11. **树–比对 tip 集必须严格一致。**  
    IQ-TREE `-te` 使用的树与比对文件 tip 集必须完全相同。外群 KDOPS 只出现在 rooted tree 阶段，不进入 `core_asr.afa`。ASR 前必须：
    - prune KDOPS outgroup
    - 生成 `CoreTree_rooted_ingroup.treefile`
    - 运行 `assert_tip_match.py`

12. **Phase 5 全部 Apo-only。**  
    当前项目不做 Holo 祖先结构预测，不用现代配体去“诱导”祖先态。

13. **ICDC 不作为主线定量结论。**  
    仅允许作为 Discussion / Outlook 的跨证据一致性展望。

14. **CPU 与线程统一。**  
    默认 20 线程：`--thread 20` / `--cpu 20` / `-T 20` / `--threads 20` / `max_workers=20`。禁止 `--thread -1`。

15. **执行记录必须即时写入 `log.md`。**  
    禁止“跑完再补”。未记入日志，视为未执行。

---

## 2. 目录与文件契约

### 2.1 必须存在的顶层目录

```text
meta/
scripts/
results/
data/
```

### 2.2 `results/` 子目录职责

```text
results/
  01_mining/            # Phase 1 原始挖掘与命中表
  02_qc/                # 去冗余、长度/覆盖度 QC、中间审计
  03_msa_core/          # core-only MSA 与 QC
  03_msa_modules/       # 模块注释矩阵、模块序列与模块 MSA
  03_msa_full/          # subtype-local full-length stitched MSA
  04_phylogeny_asr/     # rooted trees, ASR, trait ASR, QC3
  05_struct_valid/      # assembly adjudication, AF3, ESMFold, MD, QC
  06_dca/               # core DCA 正式产物 + exploratory 附录
  meta/                 # manifest、snapshot、软件版本、模型文件
```

### 2.3 禁止事项

- 禁止把 `msa_full_*.afa` 放进 `03_msa_core/`。
- 禁止把 exploratory DCA 与正式 core DCA 混在同一目录下不做标识。
- 禁止在未标注来源与版本的情况下替换核心树文件。
- 禁止使用旧 V3.1 `mafft --add` 全长大膨胀比对作为任何 Phase 4–6 输入。

### 2.4 最小必备元文件

- `PLAN.md`
- `AGENTS.md`
- `TASKS.md`
- `log.md`
- `meta/params.json`
- `results/meta/software_versions.tsv`
- `results/meta/model_files.tsv`
- `results/meta/metrics_manifest.tsv`
- `results/meta/progress_snapshot.md`

---

## 3. 当前阶段优先级（执行顺序）

> 当前不是从头开始，而是从 **Phase 4.3 之后** 继续推进。

### P0 — 立即补同步（不阻塞计算）

1. 更新 `metrics_manifest.tsv`
2. 更新 `progress_snapshot.md`
3. 同步 `README.md` / `TASKS.md` / `CLAUDE.md` / `AGENTS.md`
4. 把 `root_scenarios.tsv` 与 `node_selection_registry.tsv` 初始化

### P1 — QC3 root robustness gate（最高优先级）

1. MFP rooted tree 归档与指标提取
2. LG+C20 rooted tree 归档与指标提取
3. 4.2 AA vs 3Di tree 对比
4. midpoint rooting
5. MAD rooting
6. 可选：reduced representative set 上 nonreversible / rootstrap
7. 形成 `QC3_root_stability.md`
8. 形成 `root_scenarios.tsv`
9. 给每个候选 gain/loss 事件标记：
   - root_robust
   - root_sensitive
   - unsupported

### P2 — 不依赖深根唯一性的并行工作

1. 4.4 Iβ-ACT nested full-length ASR
2. 4.6 strict / relaxed trait matrix 清洗与输入准备
3. 6.1 core DCA 输入准备与 Meff/L 审计
4. 5.0 assembly adjudication 的文献与模板准备

### P3 — 只有在 QC3 后才能锁定的工作

1. 4.6 模块历史箭头最终叙事
2. 5.1 Phase 5 主线祖先节点最终名单
3. 主文 rooted history 图与核心叙事措辞

---

## 4. 分阶段执行与验收标准

### Phase 0：环境与可复现性【已完成】

**Done 条件**

- `meta/params.json`
- `results/meta/software_versions.tsv`
- `results/meta/model_files.tsv`

**V6 追加要求**

- `results/meta/metrics_manifest.tsv`
- `results/meta/progress_snapshot.md`

---

### Phase 1：数据挖掘【已完成】

**Done 条件**

- `results/01_mining/*.fasta`
- `results/01_mining/hits_*.domtbl`
- QC1 报告存在

**解释边界**

- KDOPS 在本阶段是有效的负选择与边界过滤工具；
- 不代表它自动是唯一可信的 deep-root 参考。

---

### Phase 2：质量控制与去冗余【已完成】

**Done 条件**

- `results/02_qc/nr80_*.fasta`
- `qc_length_report.md`
- seeds60 与 stepping-stone 集合文件

**硬约束**

- 任何后续代表性抽样都必须记录 seed 与 ID 列表。

---

### Phase 3.1–3.8：结构感知核心 MSA + 模块注释【已完成】

**Done 条件**

- `panel_candidates.tsv`
- `panel_manifest.tsv`
- `skeleton_core_aa.fa`
- `core_columns.mask`
- `core_global_matchonly.afa`
- `core_tree.afa`
- `core_asr.afa`
- `module_presence_absence_strict.tsv`
- `module_presence_absence_relaxed.tsv`
- QC2 报告

**硬约束**

- 核心列只能来自结构感知 + LDDT + RF-mask 路径。
- 旧 V3.1 膨胀全长 MSA 只可归档，不可再用。

---

### Phase 3.9：Profile-anchored full-length stitching【已完成】

**正式状态：已完成，不得写成待执行。**

**正式产物**

- `results/03_msa_full/msa_full_Ib_v4.afa`
- `results/03_msa_full/msa_full_Ib_column_map.tsv`
- 相关 linker / ids / 中间文件

**Done 条件**

- core 段列数与 `core_asr.afa` 完全一致
- column drift 断言通过
- 结果放在 `03_msa_full/`，不得放错目录

**解释边界**

- 当前正式 full-length stitching 仅对 subtype-local 架构执行；
- 不恢复跨亚型单一全长 MSA 主线。

---

### Phase 4：系统发育与分层 ASR【进行中】

#### Phase 4.1：rooted tree 建立【已进行】

**当前已知产物**

- `CoreTree_rooted_MFP.treefile`
- `CoreTree_rooted_LGC20.treefile`

**必须补充的伴随文件**

- 模型、对数似然、支持度摘要写入 `metrics_manifest.tsv`
- root scenario 登记到 `root_scenarios.tsv`

#### Phase 4.2：AA vs 3Di tree comparison【高优先级，必须做】

**地位提升说明**

- 这不是补充图性质的小检查；
- 当前 deep root 不稳时，它是 QC3 的正交证据之一。

**Done 条件**

- 至少在结构骨架代表集上构建 AA tree 与 3Di tree
- 输出拓扑差异、关键 split 一致性、root-adjacent 信号比较
- 写入 `QC3_root_stability.md`

#### Phase 4.3：树–比对 tip 集一致性【已完成】

**正式状态：已通过。**

**Done 条件**

- rooted tree 中 KDOPS 外群已 prune
- 生成 `CoreTree_rooted_ingroup.treefile`
- `assert_tip_match.py` 对 `core_asr.afa` 通过

**后续约束**

- 这一步完成只意味着 ASR 计算可以继续；
- **不意味着** deepest-root 历史解释可以直接锁定。

#### QC3：root robustness gate【V6 核心新增】

**目标**
判断哪些历史事件是 root-robust，哪些是 root-sensitive。

**必须纳入的 root scenarios**

1. `MFP_KDOPS`
2. `LGC20_KDOPS`
3. `MIDPOINT_INGROUP`
4. `MAD_INGROUP`

**可选增强**

5. `NONREV_REDUCED`
6. `ROOTSTRAP_REDUCED`

**必须输出**

- `results/04_phylogeny_asr/QC3_root_stability.md`
- `results/04_phylogeny_asr/root_scenarios.tsv`
- `results/04_phylogeny_asr/rooted_working_tree.treefile`
- `results/04_phylogeny_asr/node_selection_registry.tsv`

**QC3 最低验收标准**

- root scenario 至少 4 套
- 对每个候选模块 gain/loss 事件给出状态：
  - `root_robust`
  - `root_sensitive`
  - `ambiguous`
- 对每个候选 Phase 5 节点给出状态：
  - `eligible`
  - `hold`
  - `reject`

#### Phase 4.4：nested full-length ASR【可并行推进】

**当前优先目标**

- `Iβ-ACT` subtype-local nested ASR

**Done 条件**

- 局部子树定义清晰
- 输入 full-length MSA 明确来自 `results/03_msa_full/`
- 祖先全长序列输出可追溯
- 与 `node_selection_registry.tsv` 关联

**禁止事项**

- 不得用手工拼接序列代替 nested ASR 输出

#### Phase 4.6：模块 trait ASR【可并行准备，最终解释需等 QC3】

**任务**

- strict / relaxed 双版本 presence-absence 表
- PastML 或等价离散性状重建
- 事件表与 root scenarios 交叉标注

**Done 条件**

- `module_gain_loss_events.tsv`
- 每个事件标明 `root_robust` 或 `root_sensitive`

---

### Phase 5：关键祖先节点结构验证【未最终开放，需通过 QC3】

> 范围仍为：2–4 个关键节点 × Apo-only × native-oligomer × 有限验证性 MD。

#### Phase 5.0：Assembly adjudication【必须先做】

**任务**

- 文献与结构注释扫描
- dimer / tetramer 平行 AF3
- PISA 或界面证据评分
- 形成 `assembly_adjudication.tsv`

**Done 条件**

- 每个候选节点都有明确的 oligomer decision 与证据链

#### Phase 5.1：候选节点选择【V6 收紧版】

**只有符合以下全部条件才可进入主线：**

1. 事件或节点在至少 **2 个 root scenarios** 下解释一致
2. 支持度达到主阈值：**UFBoot >= 95**
3. 若可提供 SH-aLRT，则优先要求 **SH-aLRT >= 80**
4. strict / relaxed trait 结论不互相打架
5. 若需全长结构验证，则存在可解释的 nested full-length ASR 结果
6. 能在 `assembly_adjudication.tsv` 中得到明确装配体判定

**拒绝进入主线的情形**

- 只在单一根方案下成立
- 支持度偏低但被“故事性”诱导想保留
- full-length ancestor 来源不规范
- 需要靠 Holo 才显得合理

#### Phase 5.2–5.4【在 5.1 通过后才能正式启动】

- AF3 / ESMFold
- 结构 QC
- 有限 MD
- 口袋与界面稳定性读出

**Done 条件**

- `qc_struct_validation.md`
- 与 `node_selection_registry.tsv` 中节点状态一致

---

### Phase 6：Core DCA【可并行推进】

**主线仅限 core DCA。**

#### Phase 6.1：输入准备

**Done 条件**

- `core_dca.afa`
- `core_dca_stats.tsv`
- `Meff/L` 已登记到 manifest

#### Phase 6.2–6.3：正式分析

**Done 条件**

- `core_significant_couplings.tsv`
- `qc_core_dca.md`

**主线门槛**

- `Meff/L >= 3.0` 为最低允许
- `>= 5.0` 为理想写作门槛

**禁止事项**

- 模块 DCA 不得冒充主线结果
- joint DCA 不得在深度不足时强行讲跨域通信故事

---

### Phase 7：论文蓝图与叙事门控【必须显式分层】

**主文只能承诺 root-robust 的结论。**

#### Claim tiers

1. `root_robust`：可入主文主结果
2. `root_sensitive`：可入补充或 Discussion 条件性解释
3. `exploratory`：只可作为展望、方法探索或附录

**Done 条件**

- 每一张主文图都能映射到 claim tier
- 最深层极性不稳的内容不得出现在标题级结论中

---

## 5. 关键脚本接口契约

所有脚本必须：

- 支持 `--help`
- 检查输入存在性
- 自动创建输出目录
- 失败时非 0 退出码
- 写入关键参数与摘要到日志或 companion report

### 5.1 当前必须具备或补齐的脚本

- `extract_core_domains.py`
- `annotate_modules.py`
- `stitch_full_length_msa.py`
- `prune_tree.py`
- `assert_tip_match.py`
- `compare_trees.py`
- `qc_root_stability.py`
- `prepare_dca_input.py`
- `compute_meff.py`
- `dca_significance.py`
- `adjudicate_assembly.py`
- `coordinate_mapper.py`

### 5.2 脚本优先级

**最高优先级**

- `qc_root_stability.py`
- `compare_trees.py`
- `prepare_dca_input.py`
- `adjudicate_assembly.py`

**次优先级**

- `compute_meff.py`
- `dca_significance.py`
- `coordinate_mapper.py`

---

## 6. 常见失败与处理规则

| 症状                                     | 真实问题                                        | 正确处理                                      |
| ---------------------------------------- | ----------------------------------------------- | --------------------------------------------- |
| `core` 比对列数突然 >1000                | insert 未剥离 / 用错 `hmmalign --outformat afa` | 回到 Phase 3.6，重走 Stockholm -> RF-mask     |
| rooted tree 与 `core_asr.afa` tip 不匹配 | KDOPS 外群未 prune                              | 先 `prune_tree.py`，再 `assert_tip_match.py`  |
| KDOPS 在 rooted tree 中不单系            | 远外群 / LBA / root 不稳                        | 进入 QC3，多 scenario 处理，不强拗 KDOPS 单系 |
| MFP 与 LG+C20 根位点不同                 | 根不稳                                          | 主文只保留 root-robust 结论                   |
| 想直接把某祖先送进 AF3                   | 节点门控没过                                    | 先过 QC3 与 5.1                               |
| module DCA 看起来有热点                  | 低深度伪信号                                    | 只能归 exploratory                            |
| 结构验证看起来需要配体才成立             | 想靠 Holo 解释祖先态                            | 当前项目不允许，退回说明边界                  |
| README / TASKS / PLAN 数字不一致         | 文档漂移                                        | 以 `metrics_manifest.tsv` 为准，再统一回填    |

---

## 7. 执行完成的最低交付标准

以下项目缺任一项，都不能宣称 V6 主线完成：

1. `metrics_manifest.tsv` 与 `progress_snapshot.md` 已建立并持续更新
2. Phase 3.9 已在所有文档中正式改为已完成
3. QC3 root robustness gate 已完成
4. 4.2 AA vs 3Di 对比已纳入 QC3
5. `node_selection_registry.tsv` 已建立
6. Phase 5 节点仅从 `eligible` 条目中选择
7. core DCA 达到最低门槛并形成正式 QC
8. 论文主线按 `root_robust / root_sensitive / exploratory` 分层完成

---

## 8. 最终提醒（写给任何自动化代理）

这个项目当前最重要的不是“再多跑一点”，而是：

- 不要把 working tree 当作唯一历史真相；
- 不要为了完整故事牺牲证据边界；
- 不要让文档状态继续漂移；
- 不要让 exploratory 分析篡位成主线结论。

**V6 的成功标准不是“讲出最宏大的故事”，而是“只讲能够被多证据、可复现、可审稿地支撑的故事”。**
