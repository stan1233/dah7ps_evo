# :dna: DAH7PS Allosteric Evolution & Ancestral Sequence Reconstruction

> **Reconstructing the evolutionary dynamics of allosteric regulation in DAH7PS enzymes**

## Scientific Question

How did allosteric regulation evolve on top of the conserved TIM barrel catalytic core in DAH7PS? How did allosteric signals establish physical communication with the barrel?

## V5.1 Core Strategy

| Component | Approach |
|-----------|----------|
| **Sequence Mining** | HMM search of UniRef90 + KDOPS negative selection (dual-HMM) |
| **Structural Skeleton** | FoldMason easy-msa (AFDB panel + PDB anchors) → per-column LDDT core columns |
| **Core MSA** | hmmalign profile mapping → Stockholm → esl-alimask strip inserts → match-only AFA |
| **Module Annotation** | Coordinate-based (N_ext, α2β3, C_tail) + HMM-based (ACT, CM) → strict/relaxed matrices |
| **Smart Trimming** | Dual-version: ClipKIT kpic-smart-gap (tree) + minimal gap trim (ASR/DCA) |
| **Phylogeny & ASR** | Multi root scenario framework (KDOPS+MFP, KDOPS+LG+C20, midpoint, MAD) + pruned ingroup ASR + nested full-length ASR + PastML module trait ASR (strict/relaxed × scenarios) |
| **Root Robustness (QC3)** | Mandatory gate: Go / Conditional Go / Hold — controls narrative lock-in and Phase 5 node selection |
| **Structural Validation** | Assembly adjudication → AF3 Apo native-oligomer → validation MD ≥200 ns × 2 rep (Apo-only); Tier-1/Tier-2 node gating |
| **Core DCA** | plmc on core MSA (Meff/L >> 5) → contact validation + functional site enrichment |
| **Integration** | ICDC demoted to Discussion outlook: qualitative convergence of core DCA × ancestral structure × limited MD |
| **Manuscript Narrative** | Root-robust claims in main text; root-sensitive claims in Supplement/Discussion with conditional language |

## V5.1 Key Changes (vs V5.0)

- **QC3 upgraded to mandatory gate**: root robustness must be evaluated before locking deep history or Phase 5 nodes
- **Multi root scenario framework (new)**: S1 KDOPS+MFP, S2 KDOPS+LG+C20, S3 midpoint, S4 MAD, S5 optional nonreversible/rootstrap
- **Working tree vs narrative tree distinction**: current pruned ingroup tree used for computation; only root-robust events enter main text
- **Phase 5 node selection tightened**: UFBoot ≥ 95 (main), SH-aLRT ≥ 80 (preferred corroboration), root stability, annotation stability, ASR interpretability — replaces old bootstrap ≥ 70
- **Module trait ASR must run strict/relaxed × multiple root scenarios**: stability classification (robust / semi-robust / root-sensitive / annotation-sensitive / unresolved)
- **Claim tiers for manuscript (new)**: `root_robust` (main Results), `root_sensitive` (Supplement/Discussion), `exploratory` (Supplement/Outlook)
- **Metrics single source of truth (new)**: `results/meta/metrics_manifest.tsv` + `results/meta/progress_snapshot.md` — all documents reference manifest
- **Phase 3.9 confirmed complete**: profile-anchored stitching done, QC2b passed
- **4.2 AA vs 3Di tree elevated**: high-priority orthogonal evidence, feeds into QC3
- **CHECK-04 (new)**: deep root over-interpretation risk, mitigated by QC3 gate + claim tiers
- **CHECK-06 (new)**: document/numeric drift risk, mitigated by metrics manifest

## Current Progress (2026-03-19)

> All numbers from `results/meta/metrics_manifest.tsv` (single source of truth).

| Phase | Status | Key Milestone |
|---|---|---|
| 0 环境与可复现性 | ✅ | `params.json`, `software_versions.tsv`, `metrics_manifest.tsv` |
| 1 数据挖掘 | ✅ | PASS 24,202 seqs (Ia=9,204 / Ib=6,401 / II=8,597) |
| 2 QC 与去冗余 | ✅ | NR80=9,673 / seeds60=1,878 / stepping-stone=258 |
| 3.1–3.8 核心 MSA + 模块 | ✅ | `core_asr.afa` (9,393×472), 5 类模块 strict/relaxed |
| 3.9 Full-length stitching | ✅ | `msa_full_Ib_v4.afa` (47×1040), QC2b 通过 |
| 4.1 Rooted tree | 🏃 | S1 MFP ✅; S2 LG+C20 计算中; S3 midpoint ✅; S4 MAD ✅ |
| **4.2 AA vs 3Di tree** | **✅** | **nRF=0.74, 三大亚型单系一致 → QC3-YELLOW** |
| 4.3 Core ASR (S1) | ✅ | LogL=-2904285.44 |
| QC3 Root robustness | ⬜ | 需 S2 完成 |
| 5 结构验证 | ⬜ | 需 QC3 |
| 6 Core DCA | ⬜ | 6.1 可立即启动 |
| 7 论文蓝图 | ⬜ | — |

## Dataset

| Type | Raw Hits | PASS (QC1) | NR80 | Core Mapped |
|------|----------|------------|------|-------------|
| Iα   | 10,071 | 9,204 | 3,521 | — |
| Iβ   | 7,869 (post-KDOPS) | 6,401 | 3,073 | — |
| II   | 18,529 | 8,597 | 3,079 | — |
| **Total** | — | **24,202** | **9,673** | **9,393** |

Core MSA: 9,393 seqs × 521 cols (match-only), trimmed to 436 cols (tree) / 472 cols (ASR/DCA).

Structural panel: 35 structures (30 AFDB + 5 PDB). Skeleton trees: 46 seqs × 521 cols (AA: Q.PFAM+I+R4; 3Di: Q.3Di.AF+G4).

## Key Vulnerability Fixes

- **V2.0**: Type Iα/II internal insertions cannot be physically cut (Jiao 2020)
- **V2.0**: KDOPS sister clade removed by dual-HMM scoring (Yokoyama 2025)
- **V3.1**: Seed & Add replaces brute-force L-INS-i on 3k+ sequences
- **V3.1**: Mixed skeleton (PDB + stepping stones) prevents Mapping Cliff
- **V3.1**: ClipKIT replaces trimAl to preserve allosteric hinge evolution
- **V4.1**: FoldMason structural skeleton → LDDT-based core column definition
- **V4.1**: hmmalign + Stockholm strip-inserts eliminates MSA inflation
- **V4.1**: Module annotation decouples allosteric elements from core
- **V4.1**: Meff/L ≥ 3.0 hard gate prevents DCA artifacts [CHECK-03]
- **V5.0**: Assembly adjudication prevents wrong oligomeric state in ancestral structure prediction [CHECK-08]
- **V5.0**: Tree pruning ensures ASR tree–alignment tip-set consistency [Phase 4.3]
- **V5.0**: Apo-only prediction eliminates "modern ligand hallucination" risk
- **V5.0**: Core-only DCA focuses on information-rich layer; module DCA demoted to exploration
- **V5.0**: ICDC demoted to Discussion outlook — avoids overclaiming without full MD network
- **V5.1**: Multi root scenario framework prevents deep-root over-interpretation [CHECK-04]
- **V5.1**: QC3 mandatory gate separates working tree from narrative tree
- **V5.1**: Metrics manifest prevents document/numeric drift across PLAN/TASKS/README/paper [CHECK-06]
- **V5.1**: Phase 5 node selection tightened to four-way gating (support + root stability + annotation stability + ASR interpretability)
- **V5.1**: IQ-TREE 3.0.1 Q.3Di model case bug workaround: use `-m Q.3Di.AF+G4 --mdef` instead of `-mset`

## Environment

```bash
conda activate dah7ps_v4
# HMMER, MAFFT, IQ-TREE 2, ClipKIT, FoldMason, SeqKit, plmc, PastML, GROMACS
```

## References

- Jiao, W. et al. (2020). *Curr Opin Struct Biol* — DAH7PS allosteric evolution
- Yokoyama, R. et al. (2025). *Plant Direct* — KDOPS/DAH7PS phylogeny
- Gilchrist, C.L.M. et al. (2026). *Science* — FoldMason
- Steenwyk, J.L. et al. (2020). *PLoS Biology* — ClipKIT
- Seffernick, J.T. et al. (2025). *MBE* — 3Di substitution matrices
- Tria, F.D.K. et al. (2017). *Nat Ecol Evol* — MAD rooting
- Naser-Khdour, S. et al. (2022). *Syst Biol* — nonreversible rooting

## License

Academic research use.
