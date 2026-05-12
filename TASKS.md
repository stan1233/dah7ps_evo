# TASKS: Strategy A / Ib len255 representation

Status tags: `[ ]` pending, `[HOLD]` held, `[BLOCKED]` blocked, `[DONE]` done.

## Phase 0 - branch hygiene and docs-only reset

- [DONE] Create archive pointer for current audit branch.
- [DONE] Create new Strategy A branch from selected base.
- [DONE] Rewrite docs only: `TASKS.md`, `PLAN.md`, `AGENTS.md`, `README.md`.
- [DONE] Do not modify `log.md` in this phase.
- [DONE] Do not modify `meta/params.json` in this phase.
- [DONE] Do not modify `scripts/`, `results/`, `data/`, or `cache/` in this phase.

## Phase 1 - formal policy change

- [DONE] Change Type Ib `canonical_min` from 280 to 255 in `meta/params.json`.
- [DONE] Keep Type Ib `canonical_max` unchanged.
- [DONE] Keep `cov_min` unchanged.
- [DONE] Keep CD-HIT thresholds unchanged.
- [DONE] Record the policy change in the appropriate manifest or provenance file, not `log.md`.

## Phase 2 - formal Ib QC rerun

- [DONE] Rerun Type Ib QC.
- [DONE] Confirm `Q8U0A9` `PASS_CANONICAL`.
- [DONE] Confirm `Q9YEJ7` `PASS_CANONICAL`.
- [DONE] Audit `V8CS59` / `UniRef90_F9DH16` as the current Pni candidate.
- [DONE] Keep legacy `A0A0F2JEB6` unresolved unless accession/sequence adjudication resolves it.
- [DONE] Flag/adjudicate the one KDOPS-like rescued hit.
- [DONE] Flag/adjudicate hypothetical/uncharacterized rescued hits.
- [DONE] Flag/adjudicate other/ambiguous rescued hits.

## Phase 3 - nr80 representative audit

- [DONE] Rerun nr80.
- [DONE] Confirm `Q8U0A9` is represented by `UniRef90_UPI0002AF51CE` or updated equivalent.
- [DONE] Confirm `Q9YEJ7` is a direct nr80 representative.
- [DONE] Audit every expanded literature/anchor target for direct tip vs surrogate vs hold.
- [DONE] Produce `target_representation.tsv`.
- [DONE] Do not change CD-HIT thresholds.

## Round 3A - len255 QA handoff

- [DONE] Verify tracked Round 2 len255 outputs and validation scripts.
- [DONE] Record ignored/generated artifact provenance in `results/02_qc_len255/generated_artifacts_manifest.tsv`.
- [DONE] Review representative acceptability in `results/02_qc_len255/representative_acceptability_round3A.tsv`.
- [DONE] Flag rescued-risk records in `results/02_qc_len255/len255_flagged_rescue_adjudication.tsv`.
- [DONE] Prepare downstream handoff in `results/02_qc_len255/round3A_downstream_handoff.md`.
- [DONE] Keep final MSA/tree, IQ-TREE, QC3/rooting, and root-sensitive ASR out of Round 3A.

## Round 3B - all-subtype nr80 len255 handoff

- [DONE] Build versioned all-subtype nr80 FASTA at `results/02_qc_len255/nr80_all_len255.fasta`.
- [DONE] Audit counts, duplicate IDs, subtype composition, and target presence in `results/02_qc_len255/nr80_all_len255_target_presence.tsv`.
- [DONE] Carry forward flagged len255 rescue presence in `results/02_qc_len255/nr80_all_len255_flagged_presence.tsv`.
- [DONE] Record handoff provenance in `results/02_qc_len255/nr80_all_len255_handoff_manifest.tsv`.
- [DONE] Prepare Round 4 core-extraction handoff in `results/02_qc_len255/round3B_core_extraction_handoff.md`.
- [DONE] Keep MSA, IQ-TREE, final tree inference, QC3/rooting release, and ASR out of Round 3B.

## Phase 4 - final core-tree path

- [DONE] Rebuild versioned `nr80_all_len255`.
- [ ] Re-extract core domain sequences.
- [ ] Re-align core sequences.
- [ ] Build final core tree for representation.
- [ ] Display as unrooted/radial or visualization-rooted only.
- [ ] Do not release root-sensitive ASR or QC3.

## Phase 5 - reporting and labels

- [DONE] Update README target-status table.
- [ ] Update figure labels for direct tips and nr80 surrogates.
- [DONE] Label `Q8U0A9` as represented by nr80 surrogate, not direct tip, unless formal rerun changes that.
- [DONE] Label `Q9YEJ7` as direct nr80 tip if formal rerun confirms it.
- [DONE] Distinguish legacy `A0A0F2JEB6` from current Pni candidate `V8CS59` / `UniRef90_F9DH16`.
- [DONE] State root-sensitive ASR remains `HOLD`.

## HOLD / out of scope

- [HOLD] `noO66496 formal S1` completion.
- [HOLD] formal S2 noO66496.
- [HOLD] QC3 root stability release.
- [HOLD] root-sensitive ASR claims.
- [HOLD] exact-tip forced-retain policy.
- [HOLD] legacy `A0A0F2JEB6` unless sequence/accession adjudication resolves it.
