# DAH7PS evolution project

## Current branch purpose

This branch implements Strategy A: Type Ib short-diversity representation using a formal Type Ib `canonical_min=255` sampling policy.

Strategy A moves the main project direction away from treating `noO66496 formal S1` as the current release blocker. It preserves the O66496 and replacement-outgroup provenance, keeps rooting-sensitive claims on `HOLD`, and keeps the nr80-based final core phylogeny as the main final-tree route.

## What this branch does

- Applies a formal Type Ib lower canonical length change from 280 aa to 255 aa.
- Treats the len255 change as a broad Type Ib sampling/QC policy trial, not a one-off target rescue.
- Preserves nr80 as the final core-tree representation layer.
- Uses seeds60 only for diversity panels or structure-candidate preselection.
- Represents short and literature-relevant Type Ib diversity through direct nr80 tips or documented nr80 surrogates.
- Introduces a formal literature/anchor target representation policy.
- Keeps CD-HIT thresholds unchanged.

## What this branch does not claim

- It does not release QC3 or root stability.
- It does not treat `noO66496 formal S1` as complete.
- It does not make root-sensitive ASR claims.
- It does not claim replacement outgroup rooting is solved.
- It does not claim legacy `A0A0F2JEB6` has been rescued.
- It does not claim all literature targets will appear as direct tree tips.

O66496 is nearest to Ib-labelled ingroup tips in the original MFP/LGC20 O66496-containing trees and remains excluded from formal outgroup-rooting. This exclusion does not solve rooting. QC3, root stability, and root-sensitive ASR claims remain `HOLD` and are out of scope for Strategy A.

## Target representation summary

| Target | Strategy A status | Reporting policy |
|---|---|---|
| `Q8U0A9` / `UniRef90_Q8U0A9` | Formal len255 QC status is `PASS_CANONICAL`; formal nr80 status is represented by `UniRef90_UPI0002AF51CE` at 95.04%. | Label/report `UniRef90_UPI0002AF51CE` as representing PfuDAH7PS / `Q8U0A9`; do not force-retain exact tips in this round. |
| `Q9YEJ7` / `UniRef90_Q9YEJ7` | Formal len255 QC status is `PASS_CANONICAL`; formal nr80 status is direct representative. | Label as a direct ApeDAH7PS nr80 tip. |
| `V8CS59` / `UniRef90_F9DH16` | Formal len255 QC status is `PASS_CANONICAL`; formal nr80 status is represented by `UniRef90_A0A379DXQ3` at 81.92%; curation remains required. | Track as the current PniDAH7PS candidate via nr80 surrogate; do not mark the surrogate acceptable without curation and do not merge with legacy `A0A0F2JEB6`. |
| legacy `A0A0F2JEB6` | Unresolved accession; no formal exact sequence/accession adjudication in this round. | Do not claim represented. Keep as unresolved unless replaced by a validated current sequence/accession. |
| expanded literature/anchor targets | Formally audited in `results/02_qc_len255/target_representation.tsv`. | Use direct-tip, surrogate, absent, unresolved, curation, or hold status from that table. |

Important Pni distinction: legacy `A0A0F2JEB6` remains unresolved, while `V8CS59` / `UniRef90_F9DH16` is the current PniDAH7PS candidate to audit.

## Provenance status

Ib len255 is now formally applied for the Type Ib QC/CD-HIT rerun on this branch. Formal outputs are under `results/02_qc_len255/`, including `input_staging_manifest.tsv`, `target_representation.tsv`, and `len255_rescue_summary.tsv`.

The formal len255 rescue summary reports 479 sequences moving from previously `FRAG` to `PASS_CANONICAL`; 301 of those enter len255 nr80 representatives and 60 enter len255 seeds60 representatives. The KDOPS-like, hypothetical/uncharacterized, and other/ambiguous records are flagged for manual adjudication and are not silently filtered.

Round 3B adds the versioned all-subtype nr80 handoff at `results/02_qc_len255/nr80_all_len255.fasta`, with target presence recorded in `results/02_qc_len255/nr80_all_len255_target_presence.tsv`. Round 4 core extraction remains pending; QC3, rooting, and root-sensitive ASR remain `HOLD`.

## Next steps

Use [TASKS.md](TASKS.md) for the execution checklist and [PLAN.md](PLAN.md) for the formal Strategy A roadmap.
