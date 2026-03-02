# QC — Core column definition

- MSA input: `results/03_msa_core/skeleton_aa.fa`
- LDDT report: `results/03_msa_core/skeleton_lddt_report.html`
- Reported msaLDDT (includes unscored cols in normalization): **0.2391**
- Alignment length (L): **1966**
- Sequences (N): **46**
- Valid scored columns (scores>=0): **1735**

## Parameters

- gap_fraction_max: **0.3**
- pad_residues: **20**
- lddt_min source: **from_params** (auto_inflection(knee))
- lddt_min used: **0.1814**

## Output

- Base core (before padding): **302** columns
- Core length (kept columns, after ±20 padding): **521**
- Core mask: `results/03_msa_core/core_columns.mask`
- Core AA MSA: `results/03_msa_core/skeleton_core_aa.fa`
- Core 3Di MSA: `results/03_msa_core/skeleton_core_3di.fa`

## Per-sequence non-gap length in core (sanity)

- min/median/max non-gap residues: **266 / 358 / 447**

## Notes

- Columns with score = -1 are treated as non-core and are not enabled by padding.
- If core_len is far outside ~400–600, adjust lddt_min or gap_fraction_max in params, then rerun.
