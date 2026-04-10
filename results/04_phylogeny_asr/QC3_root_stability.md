# QC3 Root Stability (V6.1 Multi-Dimensional Gate)

> This report supersedes the older single-ratio QC3 interpretation. QC3 is not closed until all four dimensions are addressed.

**Overall QC3 status:** `HOLD`

## 1. Provenance

Status: `PASS`

| Scenario | Formality | Artifact roles | Missing roles | Status |
|---|---|---|---|---|
| S1_MFP_KDOPS | formal | log,summary,treefile | - | PASS |
| S2_LGC20_KDOPS | formal | log,summary,treefile | - | PASS |
| S3_MIDPOINT_INGROUP | formal | log,summary,treefile | - | PASS |
| S4A_TOP500_PROXY | provisional | log,summary,treefile | - | PASS |
| S4B_FULLSEARCH_PROXY | provisional | log,summary,treefile | - | PASS |
| S5_NONREV_ROOTSTRAP_REDUCED | pending | none | - | PASS |

## 2. Topology / Model Sensitivity

Status: `HOLD`

| Comparison | nRF | Status | Note |
|---|---|---|---|
| S1_MFP_KDOPS vs S2_LGC20_KDOPS | 0.413 | SENSITIVE | meaningful topology/model shift between S1 and S2 |
| S3_MIDPOINT_INGROUP, S4A_TOP500_PROXY, S4B_FULLSEARCH_PROXY | reroot-only | CONDITIONAL | same ingroup base topology; sensitivity resides in root identity, not unrooted splits |
| AA vs 3Di skeleton | 0.7442 | ORTHOGONAL | structural panel check |

## 3. Root Tie / Identity

Status: `HOLD`

| Scenario | rho | Root split sizes | Root identity |
|---|---|---|---|
| S1_MFP_KDOPS | NA | 1,1,9403 | 1[KDOPS_P0A715] | 1[KDOPS_Q9ZFK4] |
| S2_LGC20_KDOPS | NA | 1,1,9403 | 1[KDOPS_P0A715] | 1[KDOPS_Q9ZFK4] |
| S3_MIDPOINT_INGROUP | NA | 3060,6333 | 3060[UniRef90_A0A022KZ14,UniRef90_A0A023J9D5,UniRef90_A0A059JHZ8] | 6333[UniRef90_A0A011R1U0,UniRef90_A0A017TFB1,UniRef90_A0A023CY89] |
| S4A_TOP500_PROXY | 0.163617 | 3060,6333 | 3060[UniRef90_A0A022KZ14,UniRef90_A0A023J9D5,UniRef90_A0A059JHZ8] | 6333[UniRef90_A0A011R1U0,UniRef90_A0A017TFB1,UniRef90_A0A023CY89] |
| S4B_FULLSEARCH_PROXY | 0.163617 | 2820,6573 | 2820[UniRef90_A0A017TFB1,UniRef90_A0A023CY89,UniRef90_A0A024QA01] | 6573[UniRef90_A0A011R1U0,UniRef90_A0A022KZ14,UniRef90_A0A023J9D5] |
| S5_NONREV_ROOTSTRAP_REDUCED | NA | NA | unavailable |
| rho tie | 0.163617 | S4A_TOP500_PROXY,S4B_FULLSEARCH_PROXY | same rho but different root identities |

## 4. Annotation Sensitivity

Status: `HOLD`

| Feature | Orthogonality | Trait ASR role | Panel calibration | Status |
|---|---|---|---|---|
| n_ext_len | primary | main | required | PASS |
| insert_len | primary | main | required | PASS |
| act_hmm | primary | main | required | PASS |
| cm_hmm | primary | main | required | PASS |
| c_ext_len | primary | main | required | PASS |
| c_residual | secondary | conditional_only | required | CONDITIONAL |
| panel35 calibration | 0/35 reviewed | - | - | HOLD |

## Gate Interpretation

- `PASS`: evidence bundle is complete and internally coherent.
- `HOLD`: missing provenance, unresolved root ties, or uncalibrated traits.
- `CONDITIONAL`: usable only for conditional conclusions, not main-text lock-in.
