# Module Boundary Robustness Report

**Date**: 2026-03-03  
**Total sequences**: 9393  
**Source**: `core_domain_coords.tsv` + `module_hits.domtbl`

---

## 1. Module Prevalence Summary

| Module | Strict (n) | Strict (%) | Relaxed (n) | Relaxed (%) | Δ (relaxed only) |
|--------|-----------|-----------|------------|------------|------------------|
| N_ext | 3130 | 33.3% | 4429 | 47.2% | 1299 |
| alpha2beta3_insert | 172 | 1.8% | 263 | 2.8% | 91 |
| ACT_domain | 47 | 0.5% | 60 | 0.6% | 13 |
| CM_domain | 408 | 4.3% | 412 | 4.4% | 4 |
| C_tail | 360 | 3.8% | 1018 | 10.8% | 658 |

## 2. Boundary Confidence Distribution

| Confidence | Count | % |
|-----------|-------|---|
| high | 7434 | 79.1% |
| medium | 1958 | 20.8% |
| low | 1 | 0.0% |

## 3. Coordinate Module Distributions

- **N_ext length**: n=9393 min=0 q1=1 median=8 q3=43 max=2189
- **α2β3 insert gap**: n=9393 min=0 q1=0 median=0 q3=0 max=368
- **C_tail length**: n=9393 min=0 q1=0 median=1 q3=5 max=1955

## 4. Thresholds Applied

| Module | Strict | Relaxed |
|--------|--------|--------|
| N_ext | ≥ 25 aa | ≥ 10 aa |
| alpha2beta3_insert | gap ≥ 5 aa | gap ≥ 1 aa |
| ACT_domain | i-Evalue ≤ 1e-05 | i-Evalue ≤ 0.001 |
| CM_domain | i-Evalue ≤ 1e-05 | i-Evalue ≤ 0.001 |
| C_tail | ≥ 25 aa & no HMM hit | ≥ 10 aa |

## 5. Strict ⊆ Relaxed Consistency

✅ **PASS**: All strict=1 entries are also relaxed=1 (strict ⊆ relaxed holds for all modules).

