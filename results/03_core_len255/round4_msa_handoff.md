# Strategy A Round 4 MSA Handoff

## 1. FASTA For Round 5

Round 5 should align:

```text
results/03_core_len255/all_core_only_len255.fasta
```

## 2. Candidate MSA Path

No repo script appears to wrap the first core MSA step. The prior workflow in `log.md` used HMMER/Easel directly:

```bash
conda run -n dah7ps_v4 hmmalign --amino --outformat Stockholm \
  results/03_msa_core/core_global.hmm \
  results/03_core_len255/all_core_only_len255.fasta \
  > results/03_core_len255/core_global_len255.sto

conda run -n dah7ps_v4 esl-alimask --rf-is-mask \
  results/03_core_len255/core_global_len255.sto \
  > results/03_core_len255/core_global_matchonly_len255.sto

conda run -n dah7ps_v4 esl-reformat afa \
  results/03_core_len255/core_global_matchonly_len255.sto \
  > results/03_core_len255/core_global_matchonly_len255.afa
```

Round 5 should get user confirmation before running this direct command path. Do not use `hmmalign --outformat afa` directly for the core AFA.

Candidate follow-on scripts after an MSA exists:

- `scripts/minimal_trim.py`
- `scripts/define_core_columns.py`
- `scripts/merge_alignments.py`

## 3. Versioned MSA Output Names

Use versioned outputs under `results/03_core_len255/`:

- `results/03_core_len255/core_global_len255.sto`
- `results/03_core_len255/core_global_matchonly_len255.sto`
- `results/03_core_len255/core_global_matchonly_len255.afa`
- later tree/asr-specific derivatives should also use `_len255` names.

Do not overwrite legacy `results/03_msa_core/` alignments.

## 4. Target Checks After MSA

Repeat target carry-through checks after MSA:

- `UniRef90_UPI0002AF51CE` remains present for PfuDAH7PS / `Q8U0A9`.
- `UniRef90_Q9YEJ7` remains present as the ApeDAH7PS direct tip.
- `UniRef90_A0A379DXQ3` remains present for the Pni candidate and remains `needs_curation`.
- legacy `A0A0F2JEB6` remains unresolved and separate from `V8CS59` / `UniRef90_F9DH16`.

## 5. Flagged Rescue Review Before Tree Inference

Carry these flagged records into manual review before any tree inference:

- `UniRef90_A0A5E4LPI9` KDOPS-like, `manual_review_before_final_tree`.
- `UniRef90_UPI00345644EC`, `UniRef90_UPI0026F16738`, and representative `UniRef90_UPI0026F16738` for `UniRef90_UPI003593E582`, all manual review.
- `UniRef90_A0A0E3W3M4` and `UniRef90_A0A645BIV5`, annotation check.

## 6. HOLD Items

- `noO66496 formal S1`
- formal `S2 noO66496`
- QC3 root stability
- root-sensitive ASR

## 7. Forbidden Claims

- rooting solved
- QC3 released
- ASR released
- all targets direct tips
- Pni surrogate acceptable without curation
