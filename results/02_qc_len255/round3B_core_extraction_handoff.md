# Strategy A Round 3B Core Extraction Handoff

## 1. All-subtype nr80 input

Round 4 should use:

```text
results/02_qc_len255/nr80_all_len255.fasta
```

This is the versioned Strategy A all-subtype nr80 FASTA built from staged Ia, formal Ib len255, and staged II nr80 representatives.

## 2. Files not to overwrite

Do not overwrite these legacy/non-versioned outputs:

- `results/02_qc/nr80_all.fasta`
- `results/03_msa_core/nr80_all.fasta`
- `results/03_msa_core/all_core_only.fasta`
- `results/03_msa_core/core_domain_coords.tsv`
- `results/03_msa_core/core_domain_coords_domtblout.txt`
- `results/03_msa_core/core_global_matchonly.afa`
- `results/03_msa_core/core_tree.afa`
- `results/03_msa_core/core_asr.afa`

## 3. Candidate core-extraction scripts

The maintained candidate is:

- `scripts/extract_core_domains.py`

Related later-stage scripts that are not the first extraction step:

- `scripts/define_core_columns.py`
- `scripts/minimal_trim.py`
- `scripts/merge_alignments.py`

A root-level `extract_core_domains.py` duplicate also exists, but formal workflow scripts should stay under `scripts/`.

## 4. Candidate Round 4 command

The prior Phase 3.6 command pattern in `log.md` used `scripts/extract_core_domains.py` with `results/03_msa_core/core_global.hmm`. Adapted for Strategy A versioned outputs, the candidate command is:

```bash
conda run -n dah7ps_v4 python scripts/extract_core_domains.py \
  --hmm results/03_msa_core/core_global.hmm \
  --fasta results/02_qc_len255/nr80_all_len255.fasta \
  --params meta/params.json \
  --out_fasta results/03_msa_core/all_core_only_len255.fasta \
  --out_tsv results/03_msa_core/core_domain_coords_len255.tsv \
  --ievalue 1e-5 --hmm_span_min 30 --merge_gap 5 --cpu 20
```

Round 4 requires user confirmation before running because `results/03_msa_core/core_global.hmm` is not currently available in this Strategy A worktree.

## 5. Versioned output names

Use versioned Strategy A names to avoid legacy overwrite:

- `results/03_msa_core/all_core_only_len255.fasta`
- `results/03_msa_core/core_domain_coords_len255.tsv`
- `results/03_msa_core/core_domain_coords_len255_domtblout.txt`

If Round 4 continues into alignment later, keep versioned names such as `core_global_matchonly_len255.afa`, `core_tree_len255.afa`, and `core_asr_len255.afa`.

## 6. Target checks after core extraction

Repeat these checks on `all_core_only_len255.fasta` and `core_domain_coords_len255.tsv`:

- `Q8U0A9` direct tip remains not expected; `UniRef90_UPI0002AF51CE` should be present.
- `Q9YEJ7` / `UniRef90_Q9YEJ7` should be present as a direct extracted core sequence.
- `V8CS59` / `UniRef90_F9DH16` direct tip remains not expected unless a formal table says otherwise; `UniRef90_A0A379DXQ3` should be present and still needs curation.
- legacy `A0A0F2JEB6` remains unresolved and must not be merged with `V8CS59` / `UniRef90_F9DH16`.
- flagged len255 rescue representatives, especially `UniRef90_A0A5E4LPI9`, remain manual-review items before final tree release.
- repeat duplicate-ID, empty-sequence, and subtype-map checks after extraction.

## 7. HOLD items

- `noO66496 formal S1`
- formal `S2 noO66496`
- QC3 root stability
- root-sensitive ASR

## 8. Forbidden claims

- rooting solved
- QC3 released
- ASR released
- all targets direct tips
- Pni surrogate acceptable without curation
