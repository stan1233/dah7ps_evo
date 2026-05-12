# Agent instructions for Strategy A branch

## Branch purpose

This branch prioritizes Type Ib short-diversity representation using a planned Type Ib len255 policy and an nr80-based final core-tree path.

Strategy A stops treating `noO66496 formal S1` as the current release blocker. It preserves O66496 and replacement-outgroup provenance while keeping rooting, QC3, and root-sensitive ASR claims unresolved or `HOLD`.

## Evidence hierarchy

1. Formal rerun outputs generated on this branch.
2. Diagnostic cache summaries that motivated Strategy A.
3. Legacy or posthoc outputs.

Do not promote diagnostic cache outputs to formal results. The len255 diagnostic supports a formal policy trial, but formal outputs must be regenerated on this branch.

## Required behavior

- Keep rooting, QC3, root stability, and root-sensitive ASR claims `HOLD`.
- Treat `noO66496 formal S1` as incomplete unless final formal branch-local artifacts exist and the gate is explicitly released.
- Use nr80 as the final core-tree representation layer.
- Treat seeds60 as diversity-panel or preselection support, not final-tree tip presence.
- Maintain an explicit `target_representation.tsv` table.
- Separate direct nr80 tips from nr80 surrogate representation.
- Distinguish legacy `A0A0F2JEB6` from current Pni candidate `V8CS59` / `UniRef90_F9DH16`.
- Keep CD-HIT thresholds unchanged unless a later explicit task approves a separate policy change.
- Record formal rerun policy and result provenance in the appropriate manifests or representation tables, not as a strategy reset in `log.md`.

## Do not do

- Do not modify `log.md` for the strategy reset unless explicitly instructed.
- Do not claim O66496 falls inside Type Ia.
- Do not claim noO66496 solves rooting.
- Do not release QC3 or root stability.
- Do not make root-sensitive ASR claims.
- Do not change CD-HIT thresholds to force target visibility.
- Do not claim all literature targets are direct tree tips.
- Do not claim legacy `A0A0F2JEB6` / Pni is rescued.
- Do not silently merge legacy `A0A0F2JEB6` with `V8CS59` / `UniRef90_F9DH16`.
- Do not promote cache-only len255 diagnostics to formal branch results.

## Required wording

Use:

```text
O66496 is nearest to Ib-labelled ingroup tips in the original MFP/LGC20 O66496-containing trees and remains excluded from formal outgroup-rooting.
```

Use:

```text
Q8U0A9 is represented by an acceptable nr80 surrogate unless exact-tip forced retain is separately approved.
```

Use:

```text
Q9YEJ7 is expected as a direct nr80 representative under len255, pending formal rerun confirmation.
```

Use:

```text
legacy A0A0F2JEB6 remains unresolved; V8CS59 / UniRef90_F9DH16 is the current Pni candidate to audit.
```

Use:

```text
QC3/root stability/root-sensitive ASR claims remain HOLD and are out of scope for Strategy A.
```

## File hygiene

- `cache/` remains diagnostic.
- `results/` formal outputs must not be overwritten without an explicit task entry.
- `log.md` should not be used as the strategy reset mechanism.
- The first docs-only commit must modify only `README.md`, `PLAN.md`, `TASKS.md`, and `AGENTS.md`.
- Do not modify `meta/params.json`, `scripts/`, `results/`, `data/`, or `cache/` in the docs-only reset.

## Wording to avoid

Incorrect:

```text
O66496 falls inside Type Ia.
```

Use:

```text
O66496 is nearest to Ib-labelled ingroup tips in the original MFP/LGC20 O66496-containing trees.
```

Incorrect:

```text
noO66496 solves rooting.
```

Use:

```text
O66496 exclusion is justified, but replacement outgroup rooting remains unresolved; QC3 remains HOLD.
```

Incorrect:

```text
len255 rescues Pfu, Ape, and Pni.
```

Use:

```text
len255 rescues Q8U0A9 and Q9YEJ7 at the QC layer; Q8U0A9 is represented by an nr80 surrogate, Q9YEJ7 is a direct nr80 representative pending confirmation, legacy A0A0F2JEB6 remains unresolved, and V8CS59 / UniRef90_F9DH16 is the current Pni candidate to audit.
```

Incorrect:

```text
Change CD-HIT thresholds to keep the targets.
```

Use:

```text
Keep CD-HIT thresholds unchanged; use explicit target-representation policy.
```

Incorrect:

```text
All literature targets will appear as direct tree tips.
```

Use:

```text
Literature targets may be direct nr80 tips, represented by nr80 surrogates, absent, unresolved, or held for curation.
```

Incorrect:

```text
Posthoc noO66496 supports formal release.
```

Use:

```text
Posthoc noO66496 is diagnostic/provenance only and does not release QC3.
```
