# :dna: DAH7PS Allosteric Evolution & Ancestral Sequence Reconstruction

> **Reconstructing the evolutionary dynamics of allosteric regulation in DAH7PS enzymes**

## Scientific Question

How did allosteric regulation evolve on top of the conserved TIM barrel catalytic core in DAH7PS? How did allosteric signals establish physical communication with the barrel?

## V5.0 Core Strategy

| Component | Approach |
|-----------|----------|
| **Sequence Mining** | HMM search of UniRef90 + KDOPS negative selection (dual-HMM) |
| **Structural Skeleton** | FoldMason easy-msa (AFDB panel + PDB anchors) → per-column LDDT core columns |
| **Core MSA** | hmmalign profile mapping → Stockholm → esl-alimask strip inserts → match-only AFA |
| **Module Annotation** | Coordinate-based (N_ext, α2β3, C_tail) + HMM-based (ACT, CM) → strict/relaxed matrices |
| **Smart Trimming** | Dual-version: ClipKIT kpic-smart-gap (tree) + minimal gap trim (ASR/DCA) |
| **Phylogeny & ASR** | KDOPS-rooted core tree (MFP + LG+C20 LBA check) + pruned ingroup ASR + nested full-length ASR + PastML module trait ASR |
| **Structural Validation** | Assembly adjudication (dimer vs tetramer) → AF3 Apo native-oligomer → validation MD ≥200 ns × 2 rep (V5.0: no Holo) |
| **Core DCA** | plmc on core MSA (Meff/L >> 5) → contact validation + functional site enrichment |
| **Integration** | ICDC demoted to Discussion outlook (V5.0): qualitative convergence of core DCA × ancestral structure × limited MD |

## V5.0 Key Changes (vs V4.1)

- **Phase 5 reduced**: 8 nodes × 2×2 factorial MD → 2–4 nodes × Apo-only × native-oligomer validation MD
- **Assembly adjudication (new)**: oligomeric state explicitly determined per ancestor node via literature survey + parallel AF3 dimer/tetramer + PISA scoring; no default tetramer [CHECK-08]
- **Tree–alignment tip consistency (new)**: KDOPS outgroup pruned from rooted tree before ASR; tip-set identity asserted [Phase 4.3]
- **Phase 6 focused**: core-layer DCA only as main evidence; module/joint DCA → optional exploration (ACT Meff/L ≈ 0.2–0.3)
- **ICDC demoted**: from Phase 6 main conclusion → Discussion "Cross-Evidence Convergence Outlook"
- **Phase 7 added**: paper narrative blueprint with figure list and ICDC positioning
- **Path discipline**: full-length stitched MSA output unified to `results/03_msa_full/`

## Dataset

| Type | Raw | After QC1 | NR80 | Core Mapped |
|------|-----|-----------|------|-------------|
| Iα   | 10,102 | 3,473 | 3,521 | — |
| Iβ   | 16,211 | 5,728 | 3,073 | — |
| II   | 9,862 | 3,064 | 3,079 | — |
| **Total** | — | — | **9,673** | **9,393** |

Core MSA: 9,393 seqs × 521 cols (match-only), trimmed to 436 cols (tree) / 472 cols (ASR/DCA).

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

## License

Academic research use.
