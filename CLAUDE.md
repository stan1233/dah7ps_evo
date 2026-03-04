# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics research project reconstructing the evolutionary origin of allosteric regulation in DAH7PS enzymes. The workflow follows a strict phase-based SOP (V5.0) defined in `PLAN.md`. **`PLAN.md` is the single source of truth**; `AGENT.md` translates it into executable tasks; `TASKS.md` tracks completion status; `log.md` records every command executed.

## Environment

```bash
conda activate dah7ps_v4
# Tools: HMMER, MAFFT, IQ-TREE 2, ClipKIT, FoldMason, SeqKit, MMseqs2, CD-HIT, plmc, PastML
# External (HPC): GROMACS ≥2023, AlphaFold3, ESMFold
```

**Parallelism:** This machine has 28 cores. Always use **20 threads**: `--thread 20` / `--cpu 20` / `-T 20` / `-T AUTO` / `--threads 20` / `max_workers=20`. Never use `--thread -1`.

## Running Scripts

All scripts in `scripts/` must support `--help`, validate inputs, auto-create output dirs, and exit non-zero on failure.

```bash
# Example: core domain extraction (Phase 3.6)
python scripts/extract_core_domains.py --help

# Example: module annotation (Phase 3.8)
python scripts/annotate_modules.py --help
```

No top-level test suite exists. Correctness is validated via QC assertion scripts and the QC `.md` reports produced at each phase.

## Architecture

### Phase Execution Order (V5.0)

| Phase | Status | Key Output |
|-------|--------|------------|
| 0 — Environment | ✅ | `meta/params.json`, `results/meta/software_versions.tsv` |
| 1 — HMM Mining | ✅ | `results/01_mining/` — 9,673 NR80 sequences |
| 2 — QC & Dedup | ✅ | `results/02_qc/nr80_*.fasta` (Ia=3521, Ib=3073, II=3079) |
| 3.1–3.8 — Core MSA + Module Annotation | ✅ | `results/03_msa_core/core_asr.afa` (9393×472), `results/03_msa_modules/` |
| **3.9 — Full-length Stitching** | ⬜ **NEXT** | `results/03_msa_full/msa_full_Ib_v4.afa` |
| 4 — Phylogeny + ASR | ⬜ | `results/04_phylogeny_asr/` |
| 5 — Structural Validation | ⬜ | `results/05_struct_valid/` |
| 6 — Core DCA | ⬜ | `results/06_dca/` |
| 7 — Paper Blueprint | ⬜ | Narrative + Fig 1–6 |

### Core MSA Pipeline (Phases 3.1–3.6)

FoldMason `easy-msa` → per-column LDDT (knee detection) → core column mask (521 cols) → `hmmbuild` profile → `hmmalign` → **Stockholm → `esl-alimask --rf-is-mask` → AFA**. The Stockholm strip-insert path is mandatory; `hmmalign --outformat afa` is **prohibited** (causes insert column inflation).

Type II sequences require multi-domain hit stitching in `extract_core_domains.py` because the α2β3 insertion fragments HMM hits.

### Core–Module Decoupling

The core MSA (TIM-barrel) and module MSAs (ACT, CM, α2β3, N_ext, C_tail) are built independently:
- `results/03_msa_core/` — core-only alignments for phylogeny and DCA
- `results/03_msa_modules/` — per-module alignments + `module_presence_absence_strict/relaxed.tsv`
- `results/03_msa_full/` — subtype-scoped full-length stitched MSA (Phase 3.9, not yet executed); **never mix full-length MSAs into `03_msa_core/`**

### Parameters

All thresholds are externalized to `meta/params.json`. Never hardcode values in scripts; read from the params file instead.

### Hard Constraints (V5.0)

- **DCA gate**: Core DCA requires `Meff/L ≥ 3.0` (ideal ≥ 5.0). Module DCA (ACT Meff/L ≈ 0.2–0.3) is exploratory only, not main-line evidence.
- **Assembly adjudication**: Oligomeric state for ancestral nodes must be explicitly determined via Phase 5.0 (literature + parallel dimer/tetramer AF3 + PISA scoring → `assembly_adjudication.tsv`). Never hardcode "tetramer".
- **Tree–MSA tip consistency**: Before ASR, prune KDOPS outgroup tips from the rooted tree → `CoreTree_rooted_ingroup.treefile`; assert with `assert_tip_match.py`.
- **Ancestral sequences**: Full-length ancestors sent to AF3 must come from Phase 4.4 nested ASR output. Manual concatenation of core+module ASR fragments is prohibited.
- **Apo-only**: Phase 5 does not predict Holo structures (eliminates "modern ligand hallucination" risk).

### Result Preservation

Never overwrite results from a completed run. Either use a timestamped subdirectory (`results/run_YYYYMMDD/`) or keep the old files with a `.bak` suffix before producing new ones.

### QC Failures

Any QC assertion failure = **stop immediately and trace back to fix the root cause**. Never proceed with a known failure downstream.

### Logging

Every command, parameter, output file, result summary, and timestamp must be recorded in `log.md` **before moving on**, not retroactively.

### ACT Module — Low Prevalence Decision (2026-03-03)

ACT strict = 47 seqs (L=142), Meff/L ≈ 0.2–0.3. ACT DCA is excluded from the main evidence chain; exploratory appendix only.

## Common Failures (Quick Reference)

| Symptom | Cause | Fix |
|---------|-------|-----|
| Core MSA columns > 1000 | Insert columns not stripped, or `hmmalign --outformat afa` used | Return to Phase 3.6; use Stockholm → esl-alimask path |
| Type II core sequences cut in half | Hit stitching not applied | Check `extract_core_domains.py` stitching logic [CHECK-06] |
| Root position differs between MFP and LG+C20 | LBA risk | Declare root uncertainty; repeat sensitivity analysis under both root hypotheses (QC3) |
| IQ-TREE `-te` tip mismatch error | Rooted tree contains KDOPS outgroup but `core_asr.afa` does not | Run `prune_tree.py` first → `CoreTree_rooted_ingroup.treefile` |
| Module DCA looks good but Meff/L < 3 | Insufficient depth, spurious signal | Treat as noise; exclude from main-line conclusions [CHECK-03] |
| Ancestral structure interface collapse in MD | Wrong oligomeric state sent to AF3/MD | Check `assembly_adjudication.tsv`; consider alternative oligomer [CHECK-08] |
| Full-length stitched MSA not found | Wrong path | Must be in `results/03_msa_full/`, not `results/03_msa_core/` |
| Coordinate mapping mismatch | PDB insertion codes, ClipKIT trimming shifts, or GROMACS topology offset | Inspect each source of numbering change in `coordinate_mapper.py` |

## Scripts Pending Development (Phase 3.9+)

These scripts are referenced in `PLAN.md` but do not yet exist:
`extract_linkers.py`, `stitch_full_length_msa.py`, `select_sequences.py`, `prune_tree.py`, `assert_tip_match.py`, `adjudicate_assembly.py`, `prepare_dca_input.py`, `qc_root_stability.py`, `compare_trees.py`, `aggregate_gap_blocks.py`, `compute_meff.py`, `dca_significance.py`, `coordinate_mapper.py`
