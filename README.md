# DAH7PS evolution project

## Current branch purpose

This branch implements Strategy A: Type Ib short-diversity representation using a planned formal Type Ib `canonical_min=255` sampling policy.

Strategy A moves the main project direction away from treating `noO66496 formal S1` as the current release blocker. It preserves the O66496 and replacement-outgroup provenance, keeps rooting-sensitive claims on `HOLD`, and keeps the nr80-based final core phylogeny as the main final-tree route.

## What this branch does

- Plans a formal Type Ib lower canonical length change from 280 aa to 255 aa.
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
| `Q8U0A9` / `UniRef90_Q8U0A9` | Passes diagnostic len255 QC, but is absorbed at nr80 by `UniRef90_UPI0002AF51CE` at 95.04%. | Do not force-retain in round 1. Label/report `UniRef90_UPI0002AF51CE` as representing PfuDAH7PS / `Q8U0A9`, pending formal rerun confirmation. |
| `Q9YEJ7` / `UniRef90_Q9YEJ7` | Passes diagnostic len255 QC and is expected to be a direct nr80 representative because the final core tree uses nr80. | Label as ApeDAH7PS only if the formal rerun confirms direct nr80 tip status. |
| legacy `A0A0F2JEB6` | Deleted/unresolved accession; absent from the len255 CD-HIT diagnostic trial. | Do not claim rescued or represented. Keep as unresolved unless replaced by a validated current sequence/accession. |
| expanded literature/anchor targets | Not yet formally audited on this branch. | Generate `target_representation.tsv` before reporting direct-tip, surrogate, absent, unresolved, or curation status. |

Important Pni distinction: legacy `A0A0F2JEB6` remains unresolved, while `V8CS59` / `UniRef90_F9DH16` is the current PniDAH7PS candidate to audit.

## Provenance status

Diagnostic cache outputs support the Strategy A decision, but formal outputs must be regenerated on this branch. The diagnostic len255 evidence reports 479 sequences moving from previously `FRAG` to `PASS_CANONICAL`; 301 of those enter len255 nr80 representatives and 60 enter len255 seeds60 representatives.

The one KDOPS-like len255 hit and the small hypothetical, uncharacterized, other, or ambiguous set must be flagged for adjudication during the formal rerun.

## Next steps

Use [TASKS.md](TASKS.md) for the execution checklist and [PLAN.md](PLAN.md) for the formal Strategy A roadmap.
