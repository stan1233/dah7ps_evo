# Strategy A Round 3A downstream handoff

## Formal Ib replacement files

- Use `results/02_qc_len255/nr80_Ib.fasta` instead of `results/02_qc/nr80_Ib.fasta` for any Round 3B `nr80_all` construction.
- Use `results/02_qc_len255/nr80_Ib.fasta.clstr` for target/surrogate checks.
- Keep `results/02_qc_len255/qc_classification_Ib.tsv` as the formal len255 Ib QC provenance table.

## Ia/II reference files

- Keep the current staged Ia/II references for Round 3B unless a later task explicitly reruns them:
  `results/02_qc/nr80_Ia.fasta`, `results/02_qc/nr80_Ia.fasta.clstr`,
  `results/02_qc/nr80_II.fasta`, and `results/02_qc/nr80_II.fasta.clstr`.
- The staged Ia/II files are reference artifacts copied from the authorized source worktree and recorded in `results/02_qc_len255/input_staging_manifest.tsv`.

## Ignored but required next-stage artifacts

- Required ignored artifacts: `results/02_qc_len255/nr80_Ib.fasta` and `results/02_qc_len255/nr80_Ib.fasta.clstr`.
- Required tracked artifact: `results/02_qc_len255/qc_classification_Ib.tsv`.
- MD5/file-size provenance is recorded in `results/02_qc_len255/generated_artifacts_manifest.tsv`.

## Round 3B nr80_all command

Do not overwrite legacy `results/02_qc/nr80_all` paths. Use a versioned len255 path:

```bash
cat \
  results/02_qc/nr80_Ia.fasta \
  results/02_qc_len255/nr80_Ib.fasta \
  results/02_qc/nr80_II.fasta \
  > results/02_qc_len255/nr80_all_len255.fasta
```

Then record sequence count and md5 for `results/02_qc_len255/nr80_all_len255.fasta` before any core extraction.

## Round 3B representation checks

- `Q8U0A9` should not be expected as a direct nr80 tip; check that `UniRef90_UPI0002AF51CE` is present in the len255 `nr80_all` output.
- `Q9YEJ7` should be present as `UniRef90_Q9YEJ7`.
- `V8CS59` / `UniRef90_F9DH16` should not be treated as resolved by exact tip; check that surrogate `UniRef90_A0A379DXQ3` is present and keep Pni labeling under curation.

Suggested checks after Round 3B builds `nr80_all_len255.fasta`:

```bash
rg -n '^>UniRef90_UPI0002AF51CE\b|^>UniRef90_Q9YEJ7\b|^>UniRef90_A0A379DXQ3\b' results/02_qc_len255/nr80_all_len255.fasta
rg -n '^>UniRef90_Q8U0A9\b|^>UniRef90_F9DH16\b|^>A0A0F2JEB6\b' results/02_qc_len255/nr80_all_len255.fasta
```

## HOLD status

- `noO66496 formal S1` remains HOLD.
- `S2 noO66496` remains HOLD.
- QC3/root stability remains HOLD.
- Root-sensitive ASR remains HOLD.

## Claims not allowed

- Do not claim rooting solved.
- Do not claim QC3 released.
- Do not claim ASR released.
- Do not claim all targets are direct tips.
- Do not introduce forced-retain policy or change CD-HIT thresholds.
