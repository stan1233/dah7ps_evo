# QC1: Phase 1 Mining Report (V4.1)

> Generated: 2026-02-24
> HMM models: V4.1 (expanded seeds: Ia=13, Ib=6, II=14)
> Database: UniRef90 (E-value threshold: 1e-10)
> KDOPS filter: dual-HMM competitive scoring (score_margin=0)

---

## 1. Hit Summary

| Subtype | V3.1 Hits | V4.1 Hits | Change | Seeds (V3.1→V4.1) |
|---------|-----------|-----------|--------|---------------------|
| Ia | 10,102 | 10,071 | -31 (−0.3%) | 3 → 13 |
| Ib (raw) | 16,211 | 15,608 | -603 (−3.7%) | 3 → 6 |
| Ib (KDOPS-clean) | 15,731 | 7,869 | -7,862 (−50.0%) | — |
| II | 9,862 | 18,529 | +8,667 (+87.9%) | 1 → 14 |

**Key observations:**

- Type II nearly doubled — V3.1 single-seed HMM severely undersampled Type II diversity (actinobacteria, plants, etc.)
- Ia stable — expanded seeds did not substantially change coverage (already well-represented by 3 E. coli seeds)
- Ib raw slightly decreased — new 6-seed HMM is more specific (shorter model L=334 vs V3.1 L=376)
- Ib KDOPS filtering removed 49.6% (7,739 sequences) — much higher than V3.1 (3.0%, 480 sequences). This is because V3.1 used a different E-value cutoff strategy (domain conditional E-value < 1e-10 on column 12), while V4.1 uses hmmsearch -E 1e-10 (full-sequence E-value). The V4.1 Ib raw set contains more borderline hits that are actually KDOPS.

---

## 2. E-value & Score Distributions

| Subtype | E-value min | E-value median | E-value max | Score min | Score median | Score max |
|---------|-------------|----------------|-------------|-----------|--------------|-----------|
| Ia | 0.0 | 2.7e-148 | 1.0e-10 | 54.1 | 506.6 | 1133.1 |
| Ib (raw) | 2.7e-177 | 8.5e-35 | 9.5e-11 | 54.5 | 133.6 | 602.0 |
| Ib (clean) | 2.7e-177 | 6.3e-110 | 9.4e-11 | 54.5 | 380.5 | 602.0 |
| II | 0.0 | 1.2e-54 | 1.0e-10 | 53.9 | 198.8 | 1603.1 |

**Note:** Ib raw median score (133.6) is much lower than Ib clean (380.5), confirming that KDOPS contaminants cluster at low DAH7PS scores. The bimodal score distribution cleanly separates DAH7PS from KDOPS (only 4 sequences with |delta| ≤ 20).

---

## 3. Length Distributions

### Type Ia (10,071 sequences)

| Range | Count | % |
|-------|-------|---|
| <200 | 444 | 4.4% |
| 200–300 | 335 | 3.3% |
| 300–350 | 1,480 | 14.7% |
| 350–400 | 7,221 | 71.7% |
| 400–450 | 368 | 3.7% |
| 450–500 | 58 | 0.6% |
| 500–600 | 38 | 0.4% |
| 600–800 | 50 | 0.5% |
| >800 | 70 | 0.7% |

Main peak: ~350–400 aa (consistent with V3.1)

### Type Ib clean (7,869 sequences)

| Range | Count | % |
|-------|-------|---|
| <200 | 240 | 3.0% |
| 200–300 | 1,039 | 13.2% |
| 300–350 | 3,338 | 42.4% |
| 350–400 | 2,990 | 38.0% |
| 400–450 | 18 | 0.2% |
| 450–500 | 5 | 0.1% |
| 500–600 | 11 | 0.1% |
| 600–800 | 221 | 2.8% |
| >800 | 5 | 0.1% |

Main peak: ~300–400 aa. The 600–800 range likely contains ACT/CM domain fusions.

### Type II (18,529 sequences)

| Range | Count | % |
|-------|-------|---|
| <200 | 691 | 3.7% |
| 200–300 | 579 | 3.1% |
| 300–350 | 911 | 4.9% |
| 350–400 | 8,142 | 43.9% |
| 400–450 | 2,343 | 12.6% |
| 450–500 | 4,705 | 25.4% |
| 500–600 | 865 | 4.7% |
| 600–800 | 148 | 0.8% |
| >800 | 136 | 0.7% |

Bimodal: main peak ~350–400 (bacterial), shoulder ~450–500 (plant chloroplast transit peptide). Consistent with V3.1 observation.

---

## 4. KDOPS Competitive Scoring (Type Ib)

| Metric | Count |
|--------|-------|
| Total Ib sequences | 15,608 |
| Only DAH7PS hit (no KDOPS) | 136 |
| Only KDOPS hit (no DAH7PS) | 0 |
| Hit both HMMs | 15,472 |
| DAH7PS wins | 7,733 |
| KDOPS wins (score ≥ DAH7PS) | 7,739 |
| **Total removed** | **7,739** |
| **Remaining (clean)** | **7,869** |

Score delta distribution (DAH7PS − KDOPS):

| Delta range | Count | Interpretation |
|-------------|-------|----------------|
| < −100 | 7,737 | Strong KDOPS → removed |
| −100 to −50 | 2 | KDOPS → removed |
| −50 to 0 | 0 | — |
| 0 to 20 | 4 | Borderline DAH7PS → kept |
| 20 to 100 | 296 | DAH7PS → kept |
| > 100 | 7,433 | Strong DAH7PS → kept |

The distribution is strongly bimodal with virtually no ambiguous cases. The filter is robust.

---

## 5. Taxonomy Coverage (Top 10 per subtype)

### Type Ia
| Taxon | Count |
|-------|-------|
| Pseudomonas | 96 |
| Vibrio | 95 |
| Escherichia coli | 60 |
| Shewanella | 60 |
| Acinetobacter | 49 |
| Streptococcus | 48 |
| Buchnera aphidicola | 45 |
| Bifidobacterium | 44 |
| marine metagenome | 42 |
| Halomonadaceae | 40 |

5,398 unique taxa. Dominated by Proteobacteria (γ-class), with Firmicutes and Actinobacteria representation.

### Type Ib (clean)
| Taxon | Count |
|-------|-------|
| Paenibacillus | 84 |
| Clostridium | 83 |
| Bacillaceae | 72 |
| marine sediment metagenome | 59 |
| Flavobacteriaceae | 57 |
| Flavobacterium | 44 |
| Hymenobacter | 36 |
| Eubacteriales | 34 |
| bioreactor metagenome | 33 |
| Shewanella | 32 |

4,510 unique taxa. Dominated by Firmicutes (Bacillales, Clostridiales) and Bacteroidetes.

### Type II
| Taxon | Count |
|-------|-------|
| Streptomyces | 454 |
| unclassified Streptomyces | 383 |
| Pseudomonas | 118 |
| Vibrio | 94 |
| marine metagenome | 87 |
| freshwater metagenome | 77 |
| Pseudomonas aeruginosa | 70 |
| Shewanella | 60 |
| Corynebacterium | 59 |
| hydrothermal vent metagenome | 54 |

9,215 unique taxa. Heavily enriched in Actinobacteria (Streptomyces alone = 4.5%). V4.1 expanded seeds successfully captured this diversity that V3.1 missed.

---

## 6. Output Files

| File | Description | Size |
|------|-------------|------|
| `hits_Ia.domtbl` | Ia domain hits | 2.9 MB |
| `hits_Ia_ids.txt` | Ia hit IDs | 10,071 lines |
| `hits_Ia_seqs.fasta` | Ia full-length sequences | 5.0 MB |
| `hits_Ib.domtbl` | Ib domain hits (raw) | 4.8 MB |
| `hits_Ib_ids.txt` | Ib hit IDs (raw) | 15,608 lines |
| `hits_Ib_seqs.fasta` | Ib full-length sequences (raw) | 6.9 MB |
| `hits_Ib_clean.fasta` | Ib after KDOPS removal | 3.7 MB (7,869 seqs) |
| `hits_Ib_vs_dah7ps_v41.domtbl` | Ib vs DAH7PS Ib HMM | 4.8 MB |
| `hits_Ib_vs_kdops_v41.domtbl` | Ib vs KDOPS HMM | 4.7 MB |
| `kdops_contaminants_v41.txt` | Removed KDOPS IDs | 7,739 lines |
| `kdops_filter_report_v41.tsv` | Per-sequence score comparison | 698 KB |
| `hits_II.domtbl` | II domain hits | 8.5 MB |
| `hits_II_ids.txt` | II hit IDs | 18,529 lines |
| `hits_II_seqs.fasta` | II full-length sequences | 10.0 MB |

---

## 7. Conclusion

Phase 1 mining with V4.1 expanded seeds is complete. Key outcomes:

1. Type II coverage dramatically improved (+88%), validating the seed expansion strategy
2. KDOPS filtering is clean and unambiguous (bimodal delta, no gray zone)
3. All three subtypes show expected taxonomic diversity
4. Ready for Phase 2 (length filtering + de-redundancy)