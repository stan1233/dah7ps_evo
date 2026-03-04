# :dna: DAH7PS Allosteric Evolution & Ancestral Sequence Reconstruction

> **Reconstructing the evolutionary dynamics of allosteric regulation in DAH7PS enzymes**

## Scientific Question

How did allosteric regulation evolve on top of the conserved TIM barrel catalytic core in DAH7PS? How did allosteric signals establish physical communication with the barrel?

## V4.1 Core Strategy

| Component | Approach |
|-----------|----------|
| **Sequence Mining** | HMM search of UniRef90 + KDOPS negative selection (dual-HMM) |
| **Structural Skeleton** | FoldMason easy-msa (AFDB panel + PDB anchors) → per-column LDDT core columns |
| **Core MSA** | hmmalign profile mapping → Stockholm → esl-alimask strip inserts → match-only AFA |
| **Module Annotation** | Coordinate-based (N_ext, α2β3, C_tail) + HMM-based (ACT, CM) → strict/relaxed matrices |
| **Smart Trimming** | Dual-version: ClipKIT kpic-smart-gap (tree) + minimal gap trim (ASR/DCA) |
| **ASR** | Nested ASR + AltAll ensemble sampling (PP₁=0.80, PP₂=0.20) |
| **Dynamics** | AlphaFold 3 → GROMACS MD ≥500ns×3 → RMSF/DCCM |
| **Integration** | ICDC: DCA co-evolution × DCCM dynamic correlation |

## Dataset

| Type | Raw | After QC1 | NR80 | Core Mapped |
|------|-----|-----------|------|-------------|
| Iα   | 10,102 | 3,473 | 3,521 | — |
| Iβ   | 16,211 | 5,728 | 3,073 | — |
| II   | 9,862 | 3,064 | 3,079 | — |
| **Total** | — | — | **9,673** | **9,393** |

Core MSA: 9,393 seqs × 521 cols (match-only), trimmed to 436 cols (tree) / 472 cols (ASR).

## Key Vulnerability Fixes

- **V2.0**: Type Iα/II internal insertions cannot be physically cut (Jiao 2020)
- **V2.0**: KDOPS sister clade removed by dual-HMM scoring (Yokoyama 2025)
- **V3.1**: Seed & Add replaces brute-force L-INS-i on 3k+ sequences
- **V3.1**: Mixed skeleton (PDB + stepping stones) prevents Mapping Cliff
- **V3.1**: ClipKIT replaces trimAl to preserve allosteric hinge evolution
- **V4.1**: FoldMason structural skeleton → LDDT-based core column definition
- **V4.1**: hmmalign + Stockholm strip-inserts eliminates MSA inflation
- **V4.1**: Module annotation decouples allosteric elements from core
- **V4.1**: Apo-first gating prevents AF3 "modern ligand hallucination" [CHECK-07]
- **V4.1**: Meff/L ≥ 3.0 hard gate prevents DCA artifacts [CHECK-03]

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