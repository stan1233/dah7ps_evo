# QC — Core alignment (Phase 3.6)

## Input

- Core HMM: `results/03_msa_core/core_global.hmm` (L=521 match states, `--symfrac 0.0`)
- Input sequences: `results/03_msa_core/nr80_all.fasta` (9,673 seqs = Ia 3,521 + Ib 3,073 + II 3,079)

## Core domain extraction (`extract_core_domains.py`)

- i-Evalue threshold: 1e-5
- HMM span minimum: 30 aa
- Merge gap (HMM space): 5
- Coverage threshold: 0.70 (from `params.json`)
- Hit stitching: **enabled** (CHECK-06 compliant)

| Metric | Count |
|---|---|
| Input sequences | 9,673 |
| Sequences with qualifying hits | 9,619 |
| No qualifying hit | 54 |
| Failed coverage (< 0.70) | 226 |
| **Passed (output)** | **9,393** |
| Stitched (>1 hit merged) | 1,069 |

Outputs:
- `all_core_only.fasta` (9,393 seqs)
- `core_domain_coords.tsv` (9,393 rows)

## Alignment pipeline (AGENT-mandated)

1. `hmmalign --amino --outformat Stockholm` → `core_global.sto`
2. `esl-alimask --rf-is-mask` → `core_global_matchonly.sto` (insert columns stripped)
3. `esl-reformat afa` → `core_global_matchonly.afa`

**No `hmmalign --outformat afa` shortcut used** — fully compliant with AGENT §0.6.

## Final alignment: `core_global_matchonly.afa`

| Metric | Value |
|---|---|
| Sequences | 9,393 |
| **Alignment columns (L)** | **521** |
| Non-gap residues min | 198 |
| Non-gap residues median | 345 |
| Non-gap residues max | 471 |
| Non-gap residues mean | 336.2 |

## Acceptance criteria

- **A** ✅ hmmbuild model length L = 521
- **B** ✅ 9,393 / 9,673 sequences passed coverage filter (97.1%)
- **C** ✅ core_global_matchonly.afa column count = 521 (no insert column inflation)

## Notes

- Column count 521 is within the target 400–600 range and exactly matches the core definition from Phase 3.4.
- No insert column inflation detected: pre-strip Stockholm contains insert columns; post-strip AFA is exactly 521 columns.
- 1,069 sequences (11.4%) required hit stitching, confirming CHECK-06 necessity (primarily Type II with α2β3 insertions).
- 54 sequences had no qualifying hit at the given thresholds — these likely represent heavily truncated or divergent entries in the nr80 set.
