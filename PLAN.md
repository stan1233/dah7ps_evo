# Strategy A plan: Ib len255 representation

## Rationale

The current Type Ib 280 aa lower bound excludes short but HMM-covered Type Ib candidates. Strategy A formally applies a Type Ib `canonical_min=255` policy for the branch-local Type Ib QC/CD-HIT rerun while keeping nr80 as the final core-tree representation layer.

Diagnostic len255 evidence motivating the formal trial:

- 479 sequences move from previously `FRAG` to `PASS_CANONICAL`.
- 301 of those 479 enter len255 nr80 representatives.
- 60 of those 479 enter len255 seeds60 representatives.
- The recovered length range is 255-279 aa, with median 268 aa.
- `cov_best` median is 0.7575.
- Annotation classes among the 479: 446 DAH7PS-like, 27 DAH7PS plus chorismate mutase, 3 hypothetical/uncharacterized, 1 KDOPS-like, and 2 other/ambiguous.

The formal core phylogeny continues to use the nr80 route, not seeds60. The goal is to represent Type Ib diversity by direct nr80 tips or explicit nr80 surrogates.

## Policy decisions

### P1. Rooting policy

- `KDOPS_O66496` remains excluded from formal outgroup-rooting inputs.
- Do not write that O66496 falls inside Type Ia.
- Supported wording: O66496 is nearest to Ib-labelled ingroup tips in the original MFP/LGC20 O66496-containing trees.
- `noO66496 formal S1` is not the Strategy A release blocker.
- This does not mean rooting is solved.
- Rooting remains unresolved because branch-local formal noO66496 S1 is not released as a completed gate; the diagnostic status to carry forward is that `S1_KDOPS11_noO66496` lacks a final branch-local `.iqtree`, UFBoot was not converged in the captured log with latest correlation 0.911, replacement outgroup `KDOPS_P0A715,KDOPS_Q9ZFK4` does not form a clean two-tip outgroup clade, and the IQ-TREE warning "branch separating outgroup not found" maps to the noO66496 MFP run using outgroup string `KDOPS_P0A715,KDOPS_Q9ZFK4`.
- QC3, root stability, and root-sensitive ASR claims remain `HOLD` and are out of scope for Strategy A.

### P2. Ib length policy

- Applied formal change: Type Ib `canonical_min` 280 -> 255.
- Type Ib `canonical_max` remains unchanged.
- `cov_min` remains unchanged.
- CD-HIT thresholds remain unchanged.
- The one KDOPS-like hit and the small hypothetical, uncharacterized, other, or ambiguous len255-rescued set must be adjudicated or explicitly flagged in the formal rerun.

### P3. Representative policy

- Final core phylogeny uses nr80, not seeds60.
- seeds60 may support diversity panels or structure candidate preselection, but seeds60 does not determine final core-tree tip presence.
- CD-HIT cluster representatives are expected behavior, not a pipeline error.
- Do not change CD-HIT thresholds just to make a literature target visible.
- Short literature targets may be direct nr80 tips, represented by acceptable nr80 surrogates, unresolved/absent, or requiring curation.
- If exact accession visibility is required later, define a separate curated forced-retain policy. It is not part of this first documentation reset.

Key Strategy A cases:

| Target | Desired label | Current Strategy A status | Policy |
|---|---|---|---|
| `Q8U0A9` / `UniRef90_Q8U0A9` | PfuDAH7PS | Formal len255 QC status is `PASS_CANONICAL`; formal nr80 status is represented by `UniRef90_UPI0002AF51CE` at 95.04%. The representative is an acceptable nr80 surrogate. | Do not force-retain in round 1. Label/report `UniRef90_UPI0002AF51CE` as representing PfuDAH7PS / `Q8U0A9`. |
| `Q9YEJ7` / `UniRef90_Q9YEJ7` | ApeDAH7PS | Formal len255 QC status is `PASS_CANONICAL`; formal nr80 status is direct representative. Because the final core tree uses nr80, it is expected to appear as a direct nr80 tip. | Label direct nr80 tip as ApeDAH7PS. |
| legacy `A0A0F2JEB6` | legacy PniDAH7PS accession | Deleted/unresolved accession; absent from the len255 CD-HIT diagnostic trial. | Do not claim rescued or represented. Keep as unresolved legacy accession unless replaced by a validated current sequence/accession. |

Important Pni rule: do not silently merge legacy `A0A0F2JEB6` with `V8CS59` / `UniRef90_F9DH16`. Treat `A0A0F2JEB6` as a legacy unresolved/deleted accession. Treat `V8CS59` / `UniRef90_F9DH16` as the current PniDAH7PS target to audit.

## Formal rerun stages

1. Done: update `meta/params.json` after the documentation reset.
2. Done: rerun Type Ib QC with `canonical_min=255`.
3. Done: rerun Type Ib CD-HIT using unchanged nr80 and seeds60 thresholds.
4. Done: rebuild versioned `nr80_all_len255` at `results/02_qc_len255/nr80_all_len255.fasta`.
5. Re-extract core domain sequences.
6. Re-align core sequences.
7. Build the core tree for representation as unrooted or visualization-rooted only.
8. Done: generate `results/02_qc_len255/target_representation.tsv`.
9. Update labels and legends according to the target representation policy.
10. Update QC notes while keeping root-sensitive claims `HOLD`.

## Target representation table

Create a formal `target_representation.tsv` covering the expanded literature and anchor target list. Do not assume direct nr80 tip status until the formal audit verifies it.

Required columns:

`abbrev`, `primary_accession`, `uniref_id`, `length`, `subtype_or_expected_group`, `current_qc_status`, `len255_qc_status`, `nr80_status`, `nr80_representative`, `nr80_identity_if_absorbed`, `nr80_local_target_coverage`, `seeds60_status`, `final_core_tree_expected_status`, `label_policy`, `evidence_source`, `curation_call`, `notes`

Allowed `final_core_tree_expected_status` values:

- `direct_nr80_tip`
- `represented_by_nr80_surrogate`
- `absent_from_input`
- `unresolved_accession`
- `needs_curation`
- `hold`

Expanded target audit list:

| Abbrev | Accession / target | Length |
|---|---|---|
| PfuDAH7PS | `Q8U0A9` / `UniRef90_Q8U0A9` | 262 aa |
| ApeDAH7PS | `Q9YEJ7` / `UniRef90_Q9YEJ7` | 270 aa |
| TmaDAH7PS | `Q9WYH8` / `UniRef90_Q9WYH8` | 338 aa |
| GspDAH7PS | `A0ACD6B8N4` | 360 aa |
| LmoDAH7PS | `Q8Y6T2` / `UniRef90_Q8Y6T2` | 361 aa |
| PniDAH7PS | `V8CS59` / `UniRef90_F9DH16` | 354 aa |
| FtuDAH7PS | `Q5NG89` / `UniRef90_Q5NG89` | 370 aa |
| NmeDAH7PS | `Q9K169` / `UniRef90_Q9K169` | 351 aa |
| SceDAH7PS ARO3 | `P14843` / `UniRef90_P14843` | 370 aa |
| SceDAH7PS ARO4 | `P32449` / `UniRef90_P32449` | 370 aa |
| EcoDAH7PS AroF | `P00888` / `UniRef90_P00888` | 356 aa |
| EcoDAH7PS AroG | `P0AB91` / `UniRef90_P0AB91` | 350 aa |
| MtuDAH7PS | `O53512` / `UniRef90_O53512` | 462 aa |
| HpyDAH7PS | `O24947` / `UniRef90_O24947` | 449 aa |
| CglDAH7PS | `P35170` / `UniRef90_P35170` | 366 aa |
| PaeDAH7PS PA2843 | `Q9I000` / `UniRef90_Q9I000` | 448 aa |
| PaeDAH7PS PA1901 | `Q7DC82` / `UniRef90_Q7DC82` | 405 aa |

## Release criteria

- Done: new Type Ib QC output exists at `results/02_qc_len255/qc_classification_Ib.tsv`.
- Done: new Type Ib nr80/seeds60 outputs exist under `results/02_qc_len255/`.
- Done: new `target_representation.tsv` exists on this branch.
- Done: versioned all-subtype nr80 handoff exists at `results/02_qc_len255/nr80_all_len255.fasta`.
- Done: post-merge target presence audit exists at `results/02_qc_len255/nr80_all_len255_target_presence.tsv`.
- Done: `Q8U0A9` and `Q9YEJ7` status is verified after the formal rerun.
- Done: the expanded target list has direct, surrogate, absent, unresolved, curation, or hold status.
- Done: KDOPS-like and ambiguous len255-rescued hits are flagged in `results/02_qc_len255/len255_rescue_summary.tsv`.
- Root-sensitive claims are absent or explicitly marked `HOLD`.
