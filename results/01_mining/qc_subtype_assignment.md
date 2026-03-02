# QC1b: Ia vs II Best-Hit Subtype Assignment Report

> Generated: 2026-02-25
> Input: results/01_mining/hits_Ia.domtbl vs results/01_mining/hits_II.domtbl

---

## 1. Overlap Summary

| Metric | Count |
|--------|-------|
| Ia total hits | 10071 |
| II total hits | 18529 |
| Ia ∩ II overlap | 8561 (85.0% of Ia) |
| Ia-only (no overlap) | 1510 |
| II-only (no overlap) | 9968 |

## 2. Assignment Results

| Assigned to | Count | % of overlap |
|-------------|-------|--------------|
| Ia | 8561 | 100.0% |
| II | 0 | 0.0% |

## 3. Confidence Distribution

| Confidence | Criterion | Count | % |
|------------|-----------|-------|---|
| HIGH | |Δ bits| ≥ 20 | 8561 | 100.0% |
| MED | 10 ≤ |Δ bits| < 20 | 0 | 0.0% |
| LOW | |Δ bits| < 10 | 0 | 0.0% |

## 4. Final Mutually Exclusive Sets

| Set | Count | Composition |
|-----|-------|-------------|
| Ia_final | 10071 | Ia-only(1510) + overlap→Ia(8561) |
| II_final | 9968 | II-only(9968) + overlap→II(0) |
| Ia∩II (should=0) | 0 | — |

## 5. Go/No-Go Evaluation

**✅ GO:** Ia_final (10071) and II_final (9968) are both substantial. LOW confidence = 0.0% (< 40%). Assignment is reliable.

## 6. Output Files

| File | Description |
|------|-------------|
| `besthit_Ia_vs_II.tsv` | Per-sequence assignment table (8561 rows) |
| `hits_Ia_final_ids.txt` | Mutually exclusive Ia IDs (10071) |
| `hits_II_final_ids.txt` | Mutually exclusive II IDs (9968) |
| `hits_IaII_lowconf_ids.txt` | LOW confidence IDs (0) |
