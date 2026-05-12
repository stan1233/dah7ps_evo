# Strategy A Round 5 Tree Handoff

## 1. Alignment For Round 6

Round 6 should use:

```text
results/04_msa_len255/core_len255.masked.afa
```

## 2. Masking Status

This is the RF-masked alignment produced by `hmmalign --outformat Stockholm`, `esl-alimask --rf-is-mask`, and `esl-reformat afa`.

## 3. Alignment Dimensions

- Sequences: 9677
- Columns: 521

## 4. Target Presence

- Q8U0A9 representative `UniRef90_UPI0002AF51CE`: true
- Q9YEJ7 direct tip `UniRef90_Q9YEJ7`: true
- Pni representative `UniRef90_A0A379DXQ3`: true; remains `needs_curation`.

## 5. Flagged Rescue Review

Flagged rescue hits present in the alignment and requiring review before tree inference:

- `UniRef90_A0A5E4LPI9`
- `UniRef90_UPI00345644EC`
- `UniRef90_UPI0026F16738`
- `UniRef90_UPI003593E582`
- `UniRef90_A0A0E3W3M4`
- `UniRef90_A0A645BIV5`

## 6. Sequences Exceeding 50 Percent Gap Or Ambiguity

- Count: 1822
- Full list: `results/04_msa_len255/core_len255_gap_or_ambiguity_gt_0_50.tsv`

## 7. Round 6 Recommendation

Do not proceed directly to tree inference. Run a pre-tree curation gate first because flagged KDOPS-like, hypothetical/uncharacterized, and ambiguous len255 rescue hits are present in the alignment, and high-gap sequences require review.

## 8. HOLD Items

- `noO66496 formal S1`
- formal `S2 noO66496`
- QC3 root stability
- root-sensitive ASR

## 9. Forbidden Claims

- rooting solved
- QC3 released
- ASR released
- all targets direct tips
- Pni surrogate acceptable without curation
