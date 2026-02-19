# 🧬 DAH7PS Allosteric Evolution & Ancestral Sequence Reconstruction

> **Reconstructing the evolutionary dynamics of allosteric regulation in DAH7PS enzymes**

## Scientific Question

How did allosteric regulation evolve on top of the conserved TIM barrel catalytic core in DAH7PS? How did allosteric signals establish physical communication with the barrel?

## V3.1 Core Strategy

| Component | Approach |
|-----------|----------|
| **Sequence Mining** | HMM search of UniRef90 + KDOPS negative selection (dual-HMM) |
| **Case A MSA** | Seed & Add: CD-HIT 60% seeds → L-INS-i backbone → incremental mapping |
| **Case B MSA** | PROMALS3D mixed skeleton (PDB anchors + evolutionary stepping stones) → MAFFT --add |
| **Smart Trimming** | ClipKIT kpi-smart-gap (preserves allosteric hinge gaps) |
| **ASR** | Nested ASR + AltAll ensemble sampling |
| **Dynamics** | AlphaFold 3 → GROMACS MD ≥500ns×3 → RMSF/DCCM |
| **Integration** | ICDC: DCA co-evolution × DCCM dynamic correlation |

## Dataset

| Type | Raw | After KDOPS Filter | After QC1 |
|------|-----|-------------------|-----------|
| Iα   | 10,102 | — | 3,473 |
| Iβ   | 16,211 | cleaned | 5,728 |
| II   | 9,862 | — | 3,064 |

## Key Vulnerability Fixes

- **V2.0**: Type Iα/II internal insertions cannot be physically cut (Jiao 2020)
- **V2.0**: KDOPS sister clade removed by dual-HMM scoring (Yokoyama 2025)
- **V3.1**: Seed & Add replaces brute-force L-INS-i on 3k+ sequences
- **V3.1**: Mixed skeleton (PDB + stepping stones) prevents Mapping Cliff
- **V3.1**: ClipKIT replaces trimAl to preserve allosteric hinge evolution

## Environment

```bash
conda activate dah7ps_evo
# HMMER, MAFFT, IQ-TREE 2, ClipKIT, CD-HIT, SeqKit, matplotlib
```

## References

- Jiao, W. et al. (2020). *Biochem J* — DAH7PS allosteric evolution
- Yokoyama, R. et al. (2025). *Plant Direct* — KDOPS/DAH7PS phylogeny
- Steenwyk, J.L. et al. (2020). *PLoS Biology* — ClipKIT

## License

Academic research use.