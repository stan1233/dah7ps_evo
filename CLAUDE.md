# CLAUDE.md

Repository guidance for coding agents working on this project.

---

## Source Of Truth

1. `PLAN.md`
2. `results/meta/metrics_manifest.tsv`
3. `results/meta/progress_snapshot.md`
4. `AGENTS.md`
5. `TASKS.md`
6. `log.md`

---

## Current V6.1 State

- Phase 0–3.9: done
- New inserted stage: `provenance repair`
- `S2 prune + assert_tip_match + ASR` is now the absolute first priority
- `S4` is no longer singular:
  - `S4A_TOP500_PROXY`
  - `S4B_FULLSEARCH_PROXY`
- `QC3` is on `HOLD`, not passed
- `Phase 5` is frozen until:
  1. S2 ASR exists
  2. cross-scenario ASR sensitivity exists
  3. QC3 multidimensional gate is re-evaluated

---

## Hard Rules

### 1. Do not write deep-root conclusions yet

Before S2 ASR is available:

- no deep-root narrative lock
- no trait ASR finalization
- no Phase 5 node lock-in

### 2. Treat provenance as a gate

For S1–S4 artifacts, prefer:

- `artifact_manifest.tsv`
- `root_scenarios.tsv`
- explicit `treefile + summary + log`

### 3. S4 tie must stay split

If `S4a` and `S4b` share rho but differ in root identity, keep both.
Do not collapse them to a single official S4.

### 4. Trait ASR must use orthogonal features

Do not treat `C_tail` as a primary evolutionary character.
Use:

- terminal extension features
- internal insert feature
- ACT HMM support
- CM HMM support
- residual C-terminal bin only as secondary/conditional

### 5. QC3 is four-dimensional

QC3 now tracks:

- provenance
- topology/model sensitivity
- root tie/identity
- annotation sensitivity

No single partition-ratio shortcut.

### 6. S5 is conditional

Only start S5 nonreversible/rootstrap if repaired S4 + S2 ASR still leave the deepest conclusions unstable.

---

## Key Files Added In V6.1

- `results/04_phylogeny_asr/artifact_manifest.tsv`
- `results/04_phylogeny_asr/root_scenarios.tsv`
- `results/04_phylogeny_asr/QC3_root_stability.md`
- `results/03_msa_modules/module_feature_registry.tsv`
- `results/03_msa_modules/module_feature_matrix.tsv`
- `results/03_msa_modules/panel35_feature_calibration.tsv`
- `scripts/build_artifact_manifest.py`
- `scripts/cross_scenario_asr_sensitivity.py`
- `scripts/recode_module_features.py`

---

## Default Next Action

Unless the user explicitly says otherwise, the next operational step is:

1. prune S2
2. assert tip match
3. run S2 ASR
4. run cross-scenario ASR sensitivity
5. re-open QC3
