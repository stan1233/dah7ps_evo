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

- [ ] Change Type Ib `canonical_min` from 280 to 255 in `meta/params.json`.
- [ ] Keep Type Ib `canonical_max` unchanged.
- [ ] Keep `cov_min` unchanged.
- [ ] Keep CD-HIT thresholds unchanged.
- [ ] Record the policy change in the appropriate manifest or provenance file, not `log.md`.

## Phase 2 - formal Ib QC rerun

- [ ] Rerun Type Ib QC.
- [ ] Confirm `Q8U0A9` `PASS_CANONICAL`.
- [ ] Confirm `Q9YEJ7` `PASS_CANONICAL`.
- [ ] Audit `V8CS59` / `UniRef90_F9DH16` as the current Pni candidate.
- [ ] Keep legacy `A0A0F2JEB6` unresolved unless accession/sequence adjudication resolves it.
- [ ] Flag/adjudicate the one KDOPS-like rescued hit.
- [ ] Flag/adjudicate hypothetical/uncharacterized rescued hits.
- [ ] Flag/adjudicate other/ambiguous rescued hits.

## Phase 3 - nr80 representative audit

- [ ] Rerun nr80.
- [ ] Confirm `Q8U0A9` is represented by `UniRef90_UPI0002AF51CE` or updated equivalent.
- [ ] Confirm `Q9YEJ7` is a direct nr80 representative.
- [ ] Audit every expanded literature/anchor target for direct tip vs surrogate vs hold.
- [ ] Produce `target_representation.tsv`.
- [ ] Do not change CD-HIT thresholds.

## Phase 4 - final core-tree path

- [ ] Rebuild `nr80_all`.
- [ ] Re-extract core domain sequences.
- [ ] Re-align core sequences.
- [ ] Build final core tree for representation.
- [ ] Display as unrooted/radial or visualization-rooted only.
- [ ] Do not release root-sensitive ASR or QC3.

## Phase 5 - reporting and labels

- [ ] Update README target-status table.
- [ ] Update figure labels for direct tips and nr80 surrogates.
- [ ] Label `Q8U0A9` as represented by nr80 surrogate, not direct tip, unless formal rerun changes that.
- [ ] Label `Q9YEJ7` as direct nr80 tip if formal rerun confirms it.
- [ ] Distinguish legacy `A0A0F2JEB6` from current Pni candidate `V8CS59` / `UniRef90_F9DH16`.
- [ ] State root-sensitive ASR remains `HOLD`.

## HOLD / out of scope

- [HOLD] `noO66496 formal S1` completion.
- [HOLD] formal S2 noO66496.
- [HOLD] QC3 root stability release.
- [HOLD] root-sensitive ASR claims.
- [HOLD] exact-tip forced-retain policy.
- [HOLD] legacy `A0A0F2JEB6` unless sequence/accession adjudication resolves it.
