# DAH7PS V5.1 QC3 Root Robustness Gate Report

> Generated: 2026-04-08  
> Script: `scripts/qc_root_stability.py`  
> Status: **YELLOW**  

---

## 1. Scenarios Summary

| Scenario | Method | Tree File | Tips | Root Partition Ratio | Status |
|---|---|---|---|---|---|
| S1_MFP_KDOPS | outgroup | `CoreTree_rooted_MFP.treefile` | 9405 | 1.000 | ✅ |
| S2_LGC20_KDOPS | outgroup | `CoreTree_rooted_LGC20.treefile` | 9405 | 1.000 | ✅ |
| S3_MIDPOINT | midpoint | `CoreTree_rooted_midpoint_ingroup.treefile` | 9393 | 0.326 | ✅ |
| S4_MAD | MAD | `CoreTree_rooted_MAD_ingroup.treefile` | 9393 | 0.300 | ✅ |

## 2. Root Children Analysis

### S1_MFP_KDOPS — S1 MFP+KDOPS outgroup

- **Model**: `Q.PFAM+F+R10`  
- **LogL**: -2835379.61  
- **Total tips**: 9405  
- **KDOPS monophyletic**: ⚠️ NO (polyphyletic)  
- **KDOPS count**: 12  

**Root children:**

| Child | Total tips | KDOPS | Ingroup | Confidence | Example tips |
|---|---|---|---|---|---|
| #0 | 9403 | 10 | 9393 | 79 | UniRef90_A0A1Q3DYP8, UniRef90_G4TQP6, UniRef90_A0A0C2W198 |
| #1 | 1 | 1 | 0 | None | KDOPS_P0A715 |
| #2 | 1 | 1 | 0 | None | KDOPS_Q9ZFK4 |

**Deep node support values (depth, confidence):**

| Depth | Confidence |
|---|---|
| 0 | None |
| 1 | 79 |
| 2 | 100 |
| 3 | 17 |
| 4 | 67 |
| 4 | 100 |
| 5 | 69 |
| 6 | 68 |

### S2_LGC20_KDOPS — S2 LG+C20+KDOPS outgroup

- **Model**: `LG+C20+F+G`  
- **LogL**: -2820118.09  
- **Total tips**: 9405  
- **KDOPS monophyletic**: ⚠️ NO (polyphyletic)  
- **KDOPS count**: 12  

**Root children:**

| Child | Total tips | KDOPS | Ingroup | Confidence | Example tips |
|---|---|---|---|---|---|
| #0 | 9403 | 10 | 9393 | 100 | UniRef90_A0A1Q3DYP8, UniRef90_G4TQP6, UniRef90_A0A0C2W198 |
| #1 | 1 | 1 | 0 | None | KDOPS_P0A715 |
| #2 | 1 | 1 | 0 | None | KDOPS_Q9ZFK4 |

**Deep node support values (depth, confidence):**

| Depth | Confidence |
|---|---|
| 0 | None |
| 1 | 100 |
| 2 | 100 |
| 3 | 69 |
| 4 | 69 |
| 5 | 86 |
| 5 | 100 |
| 6 | 33 |

### S3_MIDPOINT — S3 Midpoint ingroup

- **Model**: `NA`  
- **Total tips**: 9393  

**Root children:**

| Child | Total tips | KDOPS | Ingroup | Confidence | Example tips |
|---|---|---|---|---|---|
| #0 | 3060 | 0 | 3060 | 100.0 | UniRef90_A0ACC1QTS7, UniRef90_A0AAJ0G0A8, UniRef90_A0ACC0UQB3 |
| #1 | 6333 | 0 | 6333 | 100.0 | UniRef90_A0ABZ2J1H6, UniRef90_A0A1G2G5S1, UniRef90_A0A3N5C5V9 |

**Deep node support values (depth, confidence):**

| Depth | Confidence |
|---|---|
| 0 | None |
| 1 | 100.0 |
| 1 | 100.0 |
| 2 | 92.0 |
| 2 | 99.0 |
| 2 | 100.0 |
| 3 | 92.0 |
| 3 | 100.0 |

### S4_MAD — S4 MAD ingroup

- **Model**: `NA`  
- **Total tips**: 9393  

**Root children:**

| Child | Total tips | KDOPS | Ingroup | Confidence | Example tips |
|---|---|---|---|---|---|
| #0 | 6573 | 0 | 6573 | None | UniRef90_A0A1Q3DYP8, UniRef90_G4TQP6, UniRef90_A0A0C2W198 |
| #1 | 2820 | 0 | 2820 | None | UniRef90_A0A6I4I514, UniRef90_A0A367GU25, UniRef90_A0A7K1YBY4 |

**Deep node support values (depth, confidence):**

| Depth | Confidence |
|---|---|
| 0 | None |
| 1 | None |
| 1 | None |
| 2 | None |
| 2 | None |
| 2 | None |
| 2 | None |
| 3 | None |

## 3. Pairwise RF Distances (Ingroup)

| Pair | RF | Max RF | nRF | Interpretation |
|---|---|---|---|---|
| S1_MFP_KDOPS vs S2_LGC20_KDOPS | 7748 | 18782 | 0.413 | moderately different |
| S1_MFP_KDOPS vs S3_MIDPOINT | 54 | 18782 | 0.003 | very similar |
| S1_MFP_KDOPS vs S4_MAD | 52 | 18782 | 0.003 | very similar |
| S2_LGC20_KDOPS vs S3_MIDPOINT | 7752 | 18782 | 0.413 | moderately different |
| S2_LGC20_KDOPS vs S4_MAD | 7750 | 18782 | 0.413 | moderately different |
| S3_MIDPOINT vs S4_MAD | 2 | 18782 | 0.000 | very similar |

## 4. KDOPS Outgroup Assessment

| Scenario | KDOPS in tree | Monophyletic | Note |
|---|---|---|---|
| S1_MFP_KDOPS | YES (12) | NO (polyphyletic) | outgroup-based rooting |
| S2_LGC20_KDOPS | YES (12) | NO (polyphyletic) | outgroup-based rooting |
| S3_MIDPOINT | NO | N/A | midpoint rooting |
| S4_MAD | NO | N/A | MAD rooting |

> **KDOPS polyphyly in ML trees (S1/S2)** is flagged in AGENTS.md §0.2. This does not invalidate the analysis — it means KDOPS cannot serve as a clean single outgroup and confirms that root uncertainty must be managed via multi-scenario approach (V5.1 strategy).

## 5. AA vs 3Di Tree (Phase 4.2)

*(From previous analysis — 2026-03-19)*

| Metric | Value |
|---|---|
| AA model | Q.PFAM+I+R4 |
| 3Di model | Q.3Di.AF+G4 |
| Normalized RF (nRF) | 0.7442 |
| Shared bipartitions | 11 of 43 |
| Ia monophyletic (both) | YES |
| Ib monophyletic (both) | YES |
| II monophyletic (both) | YES |
| Q2 (KDOPS placement) | N/A (no KDOPS in skeleton) |
| Q3 (fatal deep conflict) | NO |
| **QC3-4.2 verdict** | **YELLOW** (high nRF expected for AA vs 3Di) |

## 6. Module Event Root-Sensitivity

| Module | Strict % | Relaxed % | Event Type | Root Sensitivity |
|---|---|---|---|---|
| N_ext | 33.3% | 47.2% | partial_loss_candidate | ❌ root_sensitive |
| alpha2beta3_insert | 1.8% | 2.8% | rare_or_lineage_specific | ✅ root_robust |
| ACT_domain | 0.5% | 0.6% | rare_or_lineage_specific | ✅ root_robust |
| CM_domain | 4.3% | 4.4% | rare_or_lineage_specific | ❌ root_sensitive |
| C_tail | 3.8% | 10.8% | rare_or_lineage_specific | ❌ root_sensitive |

> **Note on root sensitivity**: Without running full ancestral state reconstruction on S2/S3/S4 trees, module gain/loss events at the deepest nodes cannot be precisely mapped. Labels here are heuristic estimates based on prevalence. Full ASR-based event mapping pending Phase 4.6 completion.

## 7. Root Stability Verdict

| Metric | Value |
|---|---|
| Scenarios loaded | 4 / 4 |
| Root partition ratios | S1=1.0, S2=1.0, S3=0.326, S4=0.3 |
| Tight agreement (±5%) | 2/4 |
| Loose agreement (±15%) | 2/4 |
| Root stability classification | **root_sensitive** |
| **Overall QC3 verdict** | **YELLOW** |

> 🟡 **YELLOW**: Root position shows partial consistency. Deep-root conclusions must be labelled `root_sensitive` in main text. Only root-robust events (confirmed in ≥2 outgroup scenarios) may enter primary claims.

## 8. Phase 5 Gating Recommendation

Phase 5 may proceed with the following constraints:

- **Working tree**: `CoreTree_rooted_ingroup.treefile` (S1 ingroup, MFP model)
- **Node eligibility**: nodes must be `root_robust` in ≥2 outgroup scenarios (S1 + S2) to enter `node_selection_registry.tsv` as `eligible`
- **Claim tier**: conclusions from root-sensitive events → `root_sensitive` tier (Supplement/Discussion only)

## 9. Required Follow-up Actions

- [ ] Run ASR on S2 (LGC20) ingroup tree — prune KDOPS first
- [ ] Run ASR on S4 (MAD) tree
- [ ] Run Phase 4.6 trait ASR on S1 + S2 trees (strict + relaxed)
- [ ] Cross-reference module events from S1 vs S2 ASR
- [ ] Update `node_selection_registry.tsv` with final eligibility
- [ ] Optional: S5 nonreversible reduced-set for rootstrap

---

> **Document hierarchy**: This file is subordinate to `PLAN.md`. All numbers must be reflected in `results/meta/metrics_manifest.tsv` and `results/meta/progress_snapshot.md`.
