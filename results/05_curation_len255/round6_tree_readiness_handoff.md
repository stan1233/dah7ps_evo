# Strategy A Round 6 Tree Readiness Handoff

## 1. Future Alignment

A future tree round would use:

```text
results/04_msa_len255/core_len255.masked.afa
```

This is the Round 5 RF-masked AFA alignment. Round 6 did not create a filtered
alignment and did not remove sequences.

## 2. Alignment Dimensions

- Sequences: 9677
- Columns: 521
- Masking: masked

## 3. High-Gap / Ambiguity Burden

- gap_or_ambiguity_fraction > 0.50: 1822
- gap_or_ambiguity_fraction > 0.70: 0
- candidate_exclude_before_tree: 0
- manual_review_before_tree: 373
- keep_for_tree: 9303

Round 5 reported 1,822 sequences above 0.50 gap/ambiguity and 0 sequences above
0.70 gap fraction; Round 6 confirms those counts from the formal tables.

## 4. Subtype Burden

- Ia high-gap count: 9
- Ib_len255 high-gap count: 1811
- II high-gap count: 2
- unknown high-gap count: 0

Ib_len255 has materially higher gap/ambiguity burden than Ia or II. The exact
rescued479 set contributes 274 high-gap rows among
1822 total high-gap rows, so the burden is
not mainly from newly rescued 255-279 aa sequences. This argues for user curation
review and cautious tree approval, not for invalidating Strategy A
`canonical_min=255` from Round 6 evidence alone.

Subtype details are in:

```text
results/05_curation_len255/round6_subtype_gap_review.tsv
```

## 5. Target Status

- Q8U0A9 / PfuDAH7PS: represented_by_nr80_surrogate; alignment_seq_id=UniRef90_UPI0002AF51CE; present_in_alignment=true; gap_fraction=0.4798; round6_call=keep_for_tree
- Q9YEJ7 / ApeDAH7PS: direct_nr80_tip; alignment_seq_id=UniRef90_Q9YEJ7; present_in_alignment=true; gap_fraction=0.5355; round6_call=manual_review_before_tree
- V8CS59 / UniRef90_F9DH16 / PniDAH7PS candidate: represented_by_nr80_surrogate; alignment_seq_id=UniRef90_A0A379DXQ3; present_in_alignment=true; gap_fraction=0.5470; round6_call=needs_curation
- legacy A0A0F2JEB6: unresolved_accession; not rescued, represented, solved, or merged with V8CS59 / UniRef90_F9DH16.
- KDOPS-like `UniRef90_A0A5E4LPI9`: KDOPS-like; present_in_alignment=true; alignment_seq_id=UniRef90_A0A5E4LPI9; gap_fraction=0.5355; round6_call=manual_review_before_tree

The Pni surrogate remains `needs_curation`. Presence of `UniRef90_A0A379DXQ3` in
the MSA is not sufficient to call it an acceptable Pni surrogate.

## 6. Candidate Exclusions

Candidate exclusions before tree inference:

- None under the Round 6 evidence-based advisory rules.

The advisory table is:

```text
results/05_curation_len255/round6_candidate_exclusion_list.tsv
```

No sequence removal was applied.

## 7. Manual Review

Manual review before tree inference is required for 373
sequence-level calls. Top severity examples:

- `UniRef90_A0A2S6GKF5` (Ib_len255; gap_or_ambiguity=0.6180; upper-tail high-gap Ib_len255 sequence near the observed masked-alignment maximum)
- `UniRef90_UPI00196046DF` (Ib_len255; gap_or_ambiguity=0.6123; exact rescued 255-279 aa set member has gap/ambiguity above 0.50)
- `UniRef90_UPI001FE533E7` (Ib_len255; gap_or_ambiguity=0.6084; upper-tail high-gap Ib_len255 sequence near the observed masked-alignment maximum)
- `UniRef90_UPI00080C3DB4` (Ib_len255; gap_or_ambiguity=0.6065; exact rescued 255-279 aa set member has gap/ambiguity above 0.50)
- `UniRef90_A0A090IXF0` (Ib_len255; gap_or_ambiguity=0.6046; upper-tail high-gap Ib_len255 sequence near the observed masked-alignment maximum)
- `UniRef90_UPI0036612F34` (Ib_len255; gap_or_ambiguity=0.6008; upper-tail high-gap Ib_len255 sequence near the observed masked-alignment maximum)
- `UniRef90_UPI001ED9971A` (Ib_len255; gap_or_ambiguity=0.6008; upper-tail high-gap Ib_len255 sequence near the observed masked-alignment maximum)
- `UniRef90_A0ABP7I939` (Ib_len255; gap_or_ambiguity=0.6008; upper-tail high-gap Ib_len255 sequence near the observed masked-alignment maximum)
- `UniRef90_A0A942TDW4` (Ib_len255; gap_or_ambiguity=0.6008; upper-tail high-gap Ib_len255 sequence near the observed masked-alignment maximum)
- `UniRef90_A0A7Y9UVM7` (Ib_len255; gap_or_ambiguity=0.6008; upper-tail high-gap Ib_len255 sequence near the observed masked-alignment maximum)

The full advisory keep/review table is:

```text
results/05_curation_len255/round6_candidate_keep_review_list.tsv
```

## 8. Recommendation

Do not run IQ-TREE yet. The next step is user approval of the Round 6 advisory
curation calls, especially the high-gap manual-review rows, the Pni surrogate,
and the KDOPS-like flagged hit. The Round 6 readiness call is:

```text
proceed_to_user_review_before_tree
```

## 9. Remains HOLD

- noO66496 formal S1
- formal S2 noO66496
- QC3 root stability
- root-sensitive ASR

In the original MFP/LGC20 O66496-containing trees, KDOPS_O66496 is nearest to
Ib-labelled ingroup tips.

## 10. Must Not Be Claimed

- solved rooting
- released QC3/root-stability status
- released ASR status
- direct-tip status for every target
- uncurated Pni surrogate acceptability
- Type Ia placement for O66496

## 11. Round 6 Outputs

- `results/05_curation_len255/round6_curation_summary.tsv`
- `results/05_curation_len255/round6_sequence_curation_calls.tsv`
- `results/05_curation_len255/round6_high_gap_review.tsv`
- `results/05_curation_len255/round6_subtype_gap_review.tsv`
- `results/05_curation_len255/round6_target_curation_review.tsv`
- `results/05_curation_len255/round6_flagged_rescue_review.tsv`
- `results/05_curation_len255/round6_pni_surrogate_review.tsv`
- `results/05_curation_len255/round6_candidate_exclusion_list.tsv`
- `results/05_curation_len255/round6_candidate_keep_review_list.tsv`
- `results/05_curation_len255/round6_artifact_manifest.tsv`
