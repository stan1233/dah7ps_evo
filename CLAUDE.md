# CLAUDE.md

This file provides guidance to Claude Code when working inside this repository.

---

## Project Overview

This repository contains a bioinformatics research workflow for reconstructing the evolutionary origin of allosteric regulation in the DAH7PS enzyme family.

The project now follows **V6 SOP rev7**.

### Source-of-truth hierarchy

Use the following order whenever documents disagree:

1. `PLAN.md` — strategy, scope, gates, interpretation boundaries
2. `results/meta/metrics_manifest.tsv` — canonical numbers and dimensions
3. `results/meta/progress_snapshot.md` — current execution snapshot
4. `AGENTS.md` — operational contract and acceptance rules
5. `TASKS.md` — live checklist / pending items
6. `log.md` — command-by-command execution record

### Current project status

Treat the repository as being in the following state unless a newer `progress_snapshot.md` says otherwise:

- Phases 0–3.8: completed
- **Phase 3.9: completed**
- Phase 4.1: rooted trees already generated
- Phase 4.3: KDOPS outgroup pruned and tip-match issue already resolved
- Main current bottleneck: **root robustness**, not core-MSA failure

### The most important V6 mindset shift

Do **not** assume that the current rooted ingroup tree is the final historical truth.
It is a **working tree** for computation.
Deep historical interpretation must wait until **QC3 root robustness** is complete.

---

## Environment

```bash
conda activate dah7ps_v4
```

### Core tools

- HMMER
- MAFFT
- IQ-TREE 2
- FoldMason
- MMseqs2
- CD-HIT
- SeqKit
- ClipKIT
- PastML
- plmc

### External / HPC tools

- AlphaFold3
- ESMFold
- GROMACS >= 2023
- PISA or equivalent interface assessment tools

### Parallelism

This machine has 28 CPU cores.
Use **20 threads** consistently unless there is an explicit reason not to.

Examples:

- `mafft --thread 20`
- `hmmsearch --cpu 20`
- `hmmalign --cpu 20`
- `iqtree -T 20`
- `mmseqs --threads 20`
- `cd-hit -T 20`
- Python: `max_workers=20`

Never use `--thread -1`.

---

## What this project is trying to prove

The main claim is **not** “a single perfect deep-rooted history for all DAH7PS modules.”

The main claim is:

- a structurally conserved TIM-barrel core can be aligned and compared family-wide,
- allosteric modules were recruited recurrently on top of that conserved core,
- some evolutionary conclusions are robust to root uncertainty,
- core DCA and limited ancestral structure validation can support those root-robust conclusions.

Whenever there is a conflict between “telling a complete story” and “staying inside the evidence boundary,” choose the evidence boundary.

---

## Current execution priorities

### Highest priority

1. Finish **QC3 root robustness gate**
2. Elevate **Phase 4.2 AA vs 3Di tree comparison** to a formal QC3 input
3. Initialize and maintain:
   - `results/04_phylogeny_asr/root_scenarios.tsv`
   - `results/04_phylogeny_asr/node_selection_registry.tsv`
   - `results/meta/metrics_manifest.tsv`
   - `results/meta/progress_snapshot.md`

### Safe to continue in parallel

These can proceed before QC3 is fully closed:

- subtype-local nested full-length ASR (especially `Iβ-ACT`)
- strict / relaxed module trait matrix preparation
- core DCA input preparation
- assembly adjudication literature/annotation preparation

### Must wait for QC3

Do **not** lock these before QC3:

- the deepest rooted historical narrative
- final Phase 5 ancestral node list
- any “first gain” or “earliest origin” claim that depends on a single root scenario

---

## Hard constraints

### 1) The core MSA must stay insert-free

The canonical path is:

```text
FoldMason skeleton -> hmmbuild -> hmmalign (Stockholm) -> RF mask -> esl-alimask -> AFA
```

Never use `hmmalign --outformat afa` to produce the core AFA directly.
That causes insert-column inflation.

### 2) Full-length ancestors must come from nested ASR

If a full-length ancestral sequence is sent to AF3, it must come from the formal Phase 4.4 subtype-local nested ASR output.
Manual concatenation of core and module fragments is prohibited.

### 3) Tree/MSA tip sets must match exactly

Before ASR:

- prune KDOPS outgroup from the rooted tree
- generate `CoreTree_rooted_ingroup.treefile`
- assert equality with `assert_tip_match.py`

The tip mismatch problem is already considered solved for the current Phase 4.3 state, but this rule remains mandatory for any regenerated tree.

### 4) The current rooted tree is only a working tree

You may use it for:

- ASR execution
- trait reconstruction preparation
- DCA input preparation

You may **not** use it as the sole basis for final deep-history claims until QC3 is complete.

### 5) Root robustness is a gate, not a footnote

QC3 must compare multiple root scenarios.
Minimum required scenarios:

- `MFP_KDOPS`
- `LGC20_KDOPS`
- `MIDPOINT_INGROUP`
- `MAD_INGROUP`

Optional strengthened scenarios:

- `NONREV_REDUCED`
- `ROOTSTRAP_REDUCED`

### 6) DCA has a hard evidence gate

Mainline DCA is **core-only**.

- minimum mainline gate: `Meff/L >= 3.0`
- preferred writing-grade target: `Meff/L >= 5.0`

Module DCA and joint cross-domain DCA are exploratory only.
Do not promote them into the main evidence chain.

### 7) Oligomeric state must be adjudicated explicitly

Never assume tetramer.
Every ancestral node entering structural validation must first pass Phase 5.0 assembly adjudication:

- literature / annotation scan
- dimer vs tetramer parallel AF3
- interface assessment / PISA-style evidence

All downstream structural work must read the copy number from `assembly_adjudication.tsv`.

### 8) Phase 5 remains Apo-only

Do not generate Holo ancestral structures for the main workflow.
No modern-ligand-driven “rescue” of ancestral conformations.

### 9) ICDC is not a main quantitative result

Treat DCA × structure × limited MD convergence as an outlook / discussion topic, not as a closed-loop quantitative proof.

### 10) Never overwrite formal results silently

Use timestamped run directories or preserve old files with `.bak` suffixes.

### 11) Every run must be logged immediately

Every significant command, parameter set, output path, result summary, and timestamp must be written to `log.md` before moving on.

---

## Key project architecture

### Core–module decoupling

The project is intentionally split into three alignment layers:

- `results/03_msa_core/` — core-only alignments for phylogeny, ASR, DCA
- `results/03_msa_modules/` — module-specific MSAs and strict/relaxed presence-absence matrices
- `results/03_msa_full/` — subtype-scoped stitched full-length alignments for nested ASR only

Never collapse these back into one family-wide full-length alignment.
That was the failure mode of earlier workflow generations.

### Important legacy warning

Old V3.1 `mafft --add` full-length artifacts may still exist in the repository.
They are archival only and must not be reused as Phase 4–6 inputs.

---

## Key files and canonical roles

### Core alignments

- `results/03_msa_core/core_global_matchonly.afa` — full core alignment with all match states
- `results/03_msa_core/core_tree.afa` — phylogeny-trimmed core alignment
- `results/03_msa_core/core_asr.afa` — ASR / DCA trimmed core alignment

### Module annotations

- `results/03_msa_modules/module_presence_absence_strict.tsv`
- `results/03_msa_modules/module_presence_absence_relaxed.tsv`
- module-specific MSAs under `results/03_msa_modules/`

### Full-length stitched alignments

- `results/03_msa_full/msa_full_Ib_v4.afa`
- `results/03_msa_full/msa_full_Ib_column_map.tsv`

### Rooting and ASR

- `results/04_phylogeny_asr/CoreTree_rooted_MFP.treefile`
- `results/04_phylogeny_asr/CoreTree_rooted_LGC20.treefile`
- `results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile`
- `results/04_phylogeny_asr/QC3_root_stability.md`
- `results/04_phylogeny_asr/root_scenarios.tsv`
- `results/04_phylogeny_asr/node_selection_registry.tsv`

### Structural validation

- `results/05_struct_valid/assembly_adjudication.tsv`
- `results/05_struct_valid/qc_struct_validation.md`

### DCA

- `results/06_dca/core_dca.afa`
- `results/06_dca/core_dca_stats.tsv`
- `results/06_dca/core_significant_couplings.tsv`
- `results/06_dca/qc_core_dca.md`

### Metadata / synchronization

- `results/meta/software_versions.tsv`
- `results/meta/model_files.tsv`
- `results/meta/metrics_manifest.tsv`
- `results/meta/progress_snapshot.md`

---

## Rooting policy and interpretation policy

### Root scenarios are first-class objects

Every serious historical interpretation must reference an explicit root scenario.
Do not speak as if there is only one root unless QC3 explicitly supports that.

### Use these claim tiers

Every result should be mentally tagged as one of:

- `root_robust`
- `root_sensitive`
- `exploratory`

Interpretation rules:

- `root_robust` — safe for main text
- `root_sensitive` — supplement or qualified discussion only
- `exploratory` — appendix, outlook, or method note only

### Node-selection policy for Phase 5

A candidate ancestral node is **not eligible** for mainline structural validation unless all of the following are true:

1. the associated event is stable in at least 2 root scenarios
2. `UFBoot >= 95`
3. if available, `SH-aLRT >= 80`
4. strict and relaxed trait reconstructions do not fundamentally contradict each other
5. the full-length ancestor can be obtained from formal nested ASR
6. assembly adjudication yields a defensible oligomeric state

Do not fall back to the old `bootstrap >= 70` style standard for mainline node selection.

---

## What to do when uncertainty appears

### If KDOPS is non-monophyletic in the rooted ML tree

Do **not** try to force KDOPS into a single “fixed” answer.
Treat this as a root-robustness problem.
Proceed via QC3.

### If MFP and LG+C20 disagree on root position

That is not a reason to discard the project.
It means the deepest historical polarity is uncertain.
Keep computing, but restrict strong claims to root-robust conclusions.

### If module DCA appears interesting despite poor depth

Do not promote it.
Low-depth module DCA stays exploratory.

### If an attractive ancestral node fails support criteria

Reject or hold it.
Do not keep it just because the story is biologically appealing.

### If documentation numbers disagree

Trust `results/meta/metrics_manifest.tsv`.
Then sync all derived documents.

---

## Script expectations

Every script in `scripts/` should:

- support `--help`
- validate inputs
- create output directories automatically
- fail with non-zero exit code on error
- write enough metadata to make the run reproducible

### High-priority scripts for the current stage

- `qc_root_stability.py`
- `compare_trees.py`
- `prepare_dca_input.py`
- `adjudicate_assembly.py`

### Additional important scripts

- `compute_meff.py`
- `dca_significance.py`
- `coordinate_mapper.py`
- `prune_tree.py`
- `assert_tip_match.py`
- `stitch_full_length_msa.py`

---

## Common failure modes

| Symptom                                | Likely cause                               | Correct response                                             |
| -------------------------------------- | ------------------------------------------ | ------------------------------------------------------------ |
| Core alignment suddenly becomes huge   | insert columns were retained               | go back to Stockholm -> RF-mask path                         |
| IQ-TREE `-te` tip mismatch             | KDOPS outgroup not pruned before ASR       | regenerate `CoreTree_rooted_ingroup.treefile` and rerun tip assertion |
| Root differs across models             | deep root instability                      | keep working, but downgrade interpretation to root-sensitive |
| A module DCA result looks exciting     | insufficient depth / phylogenetic artifact | keep as exploratory only                                     |
| AF3 ancestor requires assumed tetramer | oligomeric state was hardcoded             | return to assembly adjudication                              |
| Full-length ancestor is built manually | pipeline shortcut                          | reject and regenerate via nested ASR                         |
| README / TASKS / PLAN numbers diverge  | documentation drift                        | fix `metrics_manifest.tsv` first, then sync outward          |

---

## Practical working style for Claude Code

When editing or extending this repository:

1. check `progress_snapshot.md` first
2. confirm whether the task is root-sensitive or root-robust
3. do not change scientific scope silently
4. update manifest/snapshot when producing final numbers
5. log commands immediately
6. preserve old outputs
7. prefer smaller, well-audited advances over large undocumented leaps

---

## Bottom line

The current project is not “off track.”
It is at the point where the difference between a publishable story and an over-claimed story is whether root uncertainty is handled explicitly.

Your job in this repository is therefore:

- keep the workflow reproducible,
- keep the evidence layers separated,
- keep the documentation synchronized,
- and only make claims that survive the declared gates.
