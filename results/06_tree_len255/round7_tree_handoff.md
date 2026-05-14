# Round 7A Tree Handoff

## Status

Round 7A is **blocked_or_pending_iqtree_completion**. The approved unrooted MFP command was started on the Round 5 masked nr80-based len255 alignment, but IQ-TREE was interrupted with SIGINT after prolonged ModelFinder runtime. No completed `.treefile` or `.iqtree` report was produced, so no placement audit or tree interpretation is released by this handoff.

This is a command/log handoff for an interrupted unrooted Strategy A len255 nr80 representation-tree attempt. No rooting, QC3, ASR, or root-sensitive evolutionary claim is released by this round. Round 6 manual-review calls remain active interpretation flags. No sequence exclusions were applied.

## Command

Alignment used:

```text
results/04_msa_len255/core_len255.masked.afa
```

The alignment is masked and contains 9,677 sequences and 521 columns, as reported by IQ-TREE.

Exact command run:

```bash
conda run -n dah7ps_v4 iqtree -s results/04_msa_len255/core_len255.masked.afa -m MFP -T 20 --seed 20260513 --prefix results/06_tree_len255/core_len255_masked_MFP_unrooted
```

The command used no `-o` option and added no KDOPS/O66496 outgroup tips.

## Completion

The run did not complete cleanly:

- `results/06_tree_len255/core_len255_masked_MFP_unrooted.treefile`: not produced
- `results/06_tree_len255/core_len255_masked_MFP_unrooted.iqtree`: not produced
- `results/06_tree_len255/core_len255_masked_MFP_unrooted.log`: interrupted-run log only
- `results/06_tree_len255/core_len255_masked_MFP_unrooted.model.gz`: interrupted ModelFinder checkpoint/diagnostic output only

ModelFinder had not completed, so no final model was selected. The latest observed log table had progressed through `DAYHOFF+R10` before interruption. These partial outputs must not be treated as formal completed tree results.

## Tip Set

The alignment has 9,677 unique FASTA headers. Because no treefile was produced, tree tip count, duplicate tree tips, missing tree tips, extra tree tips, and target placement cannot be evaluated.

No KDOPS/O66496 outgroup was added by the command. The alignment header check found zero `KDOPS_O66496` entries and zero headers containing `KDOPS`, but the tree-tip KDOPS status remains not evaluated because no treefile exists.

## Target Status

Q8U0A9 / PfuDAH7PS representation remains pending in tree context. The expected Round 7A tree tip would have been the nr80 surrogate `UniRef90_UPI0002AF51CE`, with label policy `PfuDAH7PS-like / represents Q8U0A9`; no treefile exists to confirm presence or placement.

Q9YEJ7 / ApeDAH7PS remains pending in tree context. The expected direct nr80 tip would have been `UniRef90_Q9YEJ7`; no treefile exists to confirm presence or placement.

V8CS59 / UniRef90_F9DH16 / PniDAH7PS candidate remains pending in tree context and still needs curation. The expected representative would have been `UniRef90_A0A379DXQ3`; no treefile exists to confirm placement, and Pni surrogate acceptance remains HOLD.

Legacy `A0A0F2JEB6` remains `unresolved_accession` and was not resolved by this round.

KDOPS-like `UniRef90_A0A5E4LPI9` remains a manual-review flagged hit. No treefile exists to evaluate its neighborhood, and it must not be called clean DAH7PS.

## Placement Audits

Nearest-neighbor audits for `UniRef90_UPI0002AF51CE`, `UniRef90_Q9YEJ7`, `UniRef90_A0A379DXQ3`, and `UniRef90_A0A5E4LPI9` were not generated because no completed unrooted tree exists. Round 6 manual-review calls remain active and were not converted into exclusions.

## HOLD

The following remain HOLD:

- noO66496 formal S1
- formal S2 noO66496
- QC3 root stability
- root-sensitive ASR
- Pni surrogate acceptance
- KDOPS-like flagged hit acceptance

## Must Not Be Claimed

This handoff does not support any of the following claims:

- rooting solved
- QC3 released
- ASR released
- all targets are direct tips
- Pni surrogate accepted
- a Type Ia placement claim for O66496

Supported O66496 wording remains only:

```text
In the original MFP/LGC20 O66496-containing trees, KDOPS_O66496 is nearest to Ib-labelled ingroup tips.
```

## Recommended Next Action

User decision is required before the next tree action. Options include resuming/continuing the interrupted full MFP run from checkpoint if feasible on this machine, rerunning on a longer-running compute context, or approving a revised model-search strategy. No completed Round 7A tree, tipset confirmation, target neighborhood audit, or final tree handoff exists yet.
