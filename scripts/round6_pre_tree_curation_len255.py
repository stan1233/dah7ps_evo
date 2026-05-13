#!/usr/bin/env python3
"""Strategy A Round 6 pre-tree curation gate for the len255 core MSA.

This script is intentionally read-only outside results/05_curation_len255/.
It does not run tree inference, filter alignments, or remove sequences.
"""

from __future__ import annotations

import csv
import hashlib
import statistics
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results/05_curation_len255"

ALIGNMENT_AFA = "results/04_msa_len255/core_len255.masked.afa"
ALIGNMENT_STO = "results/04_msa_len255/core_len255.masked.sto"
ALIGNMENT_SUMMARY = "results/04_msa_len255/core_len255_alignment_summary.tsv"
ALIGNMENT_QC = "results/04_msa_len255/core_len255.alignment_qc.tsv"
SUBTYPE_QC = "results/04_msa_len255/core_len255_subtype_alignment_qc.tsv"
RESCUED479_QC = "results/04_msa_len255/core_len255_rescued479_alignment_qc.tsv"
MSA_TARGETS = "results/04_msa_len255/core_len255_target_carrythrough.tsv"
MSA_FLAGGED = "results/04_msa_len255/core_len255_flagged_carrythrough.tsv"
HIGH_GAP_INPUT = "results/04_msa_len255/core_len255_gap_or_ambiguity_gt_0_50.tsv"
ROUND5_MANIFEST = "results/04_msa_len255/core_len255_artifact_manifest.tsv"
ROUND5_HANDOFF = "results/04_msa_len255/round5_tree_handoff.md"
TARGET_REPRESENTATION = "results/02_qc_len255/target_representation.tsv"
SUBTYPE_MAP = "results/02_qc_len255/nr80_all_len255_subtype_map.tsv"
CORE_TARGETS = "results/03_core_len255/core_len255_target_carrythrough.tsv"

ROUND5_INPUTS = [
    ALIGNMENT_AFA,
    ALIGNMENT_STO,
    ALIGNMENT_SUMMARY,
    ALIGNMENT_QC,
    SUBTYPE_QC,
    RESCUED479_QC,
    MSA_TARGETS,
    MSA_FLAGGED,
    HIGH_GAP_INPUT,
    ROUND5_MANIFEST,
    ROUND5_HANDOFF,
    TARGET_REPRESENTATION,
    SUBTYPE_MAP,
    CORE_TARGETS,
]

GENERATED_OUTPUTS = [
    "results/05_curation_len255/round6_curation_summary.tsv",
    "results/05_curation_len255/round6_sequence_curation_calls.tsv",
    "results/05_curation_len255/round6_high_gap_review.tsv",
    "results/05_curation_len255/round6_subtype_gap_review.tsv",
    "results/05_curation_len255/round6_target_curation_review.tsv",
    "results/05_curation_len255/round6_flagged_rescue_review.tsv",
    "results/05_curation_len255/round6_pni_surrogate_review.tsv",
    "results/05_curation_len255/round6_candidate_exclusion_list.tsv",
    "results/05_curation_len255/round6_candidate_keep_review_list.tsv",
    "results/05_curation_len255/round6_tree_readiness_handoff.md",
    "results/05_curation_len255/round6_artifact_manifest.tsv",
]

SCRIPT_PATH = "scripts/round6_pre_tree_curation_len255.py"

CALLS = {
    "keep_for_tree",
    "manual_review_before_tree",
    "candidate_exclude_before_tree",
    "needs_curation",
    "needs_external_validation",
    "hold",
}


def read_tsv(relative_path: str) -> list[dict[str, str]]:
    path = ROOT / relative_path
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(relative_path: str, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path = ROOT / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean_value(row.get(field, "")) for field in fieldnames})


def clean_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def as_float(value: str) -> float | None:
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def as_bool(value: str) -> bool:
    return value.lower() == "true"


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.4f}"


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_status_bool(args: list[str], relative_path: str) -> bool:
    result = subprocess.run(
        ["git", *args, relative_path],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def is_tracked(relative_path: str) -> bool:
    return git_status_bool(["ls-files", "--error-unmatch"], relative_path)


def is_ignored(relative_path: str) -> bool:
    return git_status_bool(["check-ignore", "-q"], relative_path)


def parse_alignment_dimensions(relative_path: str) -> tuple[int, int]:
    count = 0
    current: list[str] = []
    lengths: list[int] = []
    with (ROOT / relative_path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if current:
                    lengths.append(len("".join(current)))
                    current = []
                count += 1
            else:
                current.append(line.strip())
        if current:
            lengths.append(len("".join(current)))
    unique_lengths = set(lengths)
    if len(unique_lengths) != 1:
        raise ValueError(f"alignment has non-uniform lengths: {sorted(unique_lengths)[:5]}")
    return count, lengths[0]


def by_key(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in rows if row.get(key)}


def append_map(target: dict[str, list[dict[str, str]]], key: str, row: dict[str, str]) -> None:
    if key:
        target.setdefault(key, []).append(row)


def unique_join(values: list[str]) -> str:
    seen: list[str] = []
    for value in values:
        if value and value not in seen:
            seen.append(value)
    return ";".join(seen)


def load_context() -> dict[str, object]:
    alignment_rows = read_tsv(ALIGNMENT_QC)
    alignment_by_id = by_key(alignment_rows, "seq_id")
    target_rows = read_tsv(MSA_TARGETS)
    target_rep_rows = read_tsv(TARGET_REPRESENTATION)
    flagged_rows = read_tsv(MSA_FLAGGED)
    rescued_rows = read_tsv(RESCUED479_QC)
    subtype_rows = read_tsv(SUBTYPE_QC)
    summary_rows = read_tsv(ALIGNMENT_SUMMARY)
    high_gap_input_rows = read_tsv(HIGH_GAP_INPUT)
    core_target_rows = read_tsv(CORE_TARGETS)

    target_by_alignment: dict[str, list[dict[str, str]]] = {}
    for row in target_rows:
        if as_bool(row.get("present_in_alignment", "")):
            append_map(target_by_alignment, row.get("alignment_seq_id", ""), row)

    flags_by_alignment: dict[str, list[dict[str, str]]] = {}
    for row in flagged_rows:
        if as_bool(row.get("present_in_alignment", "")):
            append_map(flags_by_alignment, row.get("alignment_seq_id", ""), row)

    rescued_present_by_id = {
        row["seq_id"]: row
        for row in rescued_rows
        if as_bool(row.get("present_in_alignment", "")) and row.get("seq_id") in alignment_by_id
    }

    return {
        "alignment_rows": alignment_rows,
        "alignment_by_id": alignment_by_id,
        "target_rows": target_rows,
        "target_rep_rows": target_rep_rows,
        "flagged_rows": flagged_rows,
        "rescued_rows": rescued_rows,
        "subtype_rows": subtype_rows,
        "summary_rows": summary_rows,
        "high_gap_input_rows": high_gap_input_rows,
        "core_target_rows": core_target_rows,
        "target_by_alignment": target_by_alignment,
        "flags_by_alignment": flags_by_alignment,
        "rescued_present_by_id": rescued_present_by_id,
    }


def sequence_call(
    row: dict[str, str],
    target_rows: list[dict[str, str]],
    flag_rows: list[dict[str, str]],
    is_rescued479: bool,
) -> tuple[str, str, str]:
    gap = as_float(row["gap_fraction"]) or 0.0
    gap_or_ambiguity = as_float(row["gap_or_ambiguity_fraction"]) or 0.0
    subtype = row["subtype_source"]
    is_high = gap_or_ambiguity > 0.50 or gap > 0.50
    has_pni = any(target.get("abbrev") == "PniDAH7PS" for target in target_rows)

    if has_pni:
        return (
            "needs_curation",
            "current Pni surrogate remains unresolved for annotation/taxonomy acceptability",
            "presence in the masked alignment does not make the Pni surrogate acceptable",
        )
    if (gap_or_ambiguity > 0.70 or gap > 0.70) and not target_rows and not flag_rows:
        return (
            "candidate_exclude_before_tree",
            "extreme gap/ambiguity burden above 0.70 without target or flagged-rescue protection",
            "advisory only; removal requires later user approval",
        )
    if flag_rows:
        return (
            "manual_review_before_tree",
            "flagged rescue hit carried forward for explicit manual review",
            "not silently filtered and not called clean DAH7PS",
        )
    if target_rows and is_high:
        return (
            "manual_review_before_tree",
            "target or nr80 representative has gap/ambiguity above 0.50",
            "target/representative is not automatically excluded by gap fraction alone",
        )
    if is_rescued479 and is_high:
        return (
            "manual_review_before_tree",
            "exact rescued 255-279 aa set member has gap/ambiguity above 0.50",
            "review before tree but do not remove without approval",
        )
    if is_high and subtype != "Ib_len255":
        return (
            "manual_review_before_tree",
            "non-Ib high-gap outlier relative to subtype burden",
            "review before tree; no sequence removed in Round 6",
        )
    if is_high and gap_or_ambiguity > 0.58:
        return (
            "manual_review_before_tree",
            "upper-tail high-gap Ib_len255 sequence near the observed masked-alignment maximum",
            "below 0.70 and advisory review only",
        )
    if is_high:
        return (
            "keep_for_tree",
            "moderate Ib_len255 high-gap burden without target, flagged, rescued, or extreme-tail signal",
            "kept for tree pending user review of the advisory curation tables",
        )
    return (
        "keep_for_tree",
        "below Round 6 high-gap review threshold and no special curation flag",
        "kept for future tree input if user approves Round 6 curation",
    )


def build_sequence_rows(context: dict[str, object]) -> list[dict[str, object]]:
    alignment_rows = context["alignment_rows"]
    target_by_alignment = context["target_by_alignment"]
    flags_by_alignment = context["flags_by_alignment"]
    rescued_present_by_id = context["rescued_present_by_id"]

    assert isinstance(alignment_rows, list)
    assert isinstance(target_by_alignment, dict)
    assert isinstance(flags_by_alignment, dict)
    assert isinstance(rescued_present_by_id, dict)

    output: list[dict[str, object]] = []
    for row in alignment_rows:
        seq_id = row["seq_id"]
        targets = target_by_alignment.get(seq_id, [])
        flags = flags_by_alignment.get(seq_id, [])
        is_rescued479 = seq_id in rescued_present_by_id
        call, reason, notes = sequence_call(row, targets, flags, is_rescued479)
        if call not in CALLS:
            raise ValueError(f"invalid curation call {call}")
        gap = as_float(row["gap_fraction"]) or 0.0
        gap_or_ambiguity = as_float(row["gap_or_ambiguity_fraction"]) or 0.0
        target_labels = unique_join([target.get("abbrev", "") for target in targets])
        flag_classes = unique_join([flag.get("flag_class", "") for flag in flags])
        output.append(
            {
                "seq_id": seq_id,
                "subtype_source": row["subtype_source"],
                "gap_fraction": row["gap_fraction"],
                "ambiguity_fraction": row["ambiguity_fraction"],
                "gap_or_ambiguity_fraction": row["gap_or_ambiguity_fraction"],
                "is_high_gap_gt_0_50": gap_or_ambiguity > 0.50 or gap > 0.50,
                "is_high_gap_gt_0_70": gap_or_ambiguity > 0.70 or gap > 0.70,
                "is_target_or_representative": bool(targets),
                "target_label_if_any": target_labels,
                "is_flagged_rescue": bool(flags),
                "flag_class_if_any": flag_classes,
                "is_rescued479": is_rescued479,
                "round6_curation_call": call,
                "reason": reason,
                "notes": notes,
            }
        )
    return output


def build_high_gap_review(
    context: dict[str, object], sequence_rows_by_id: dict[str, dict[str, object]]
) -> list[dict[str, object]]:
    alignment_rows = context["alignment_rows"]
    target_by_alignment = context["target_by_alignment"]
    flags_by_alignment = context["flags_by_alignment"]
    rescued_present_by_id = context["rescued_present_by_id"]

    assert isinstance(alignment_rows, list)
    assert isinstance(target_by_alignment, dict)
    assert isinstance(flags_by_alignment, dict)
    assert isinstance(rescued_present_by_id, dict)

    high_rows = [
        row
        for row in alignment_rows
        if (as_float(row["gap_or_ambiguity_fraction"]) or 0.0) > 0.50
        or (as_float(row["gap_fraction"]) or 0.0) > 0.50
    ]
    high_rows.sort(
        key=lambda row: (
            as_float(row["gap_or_ambiguity_fraction"]) or 0.0,
            as_float(row["gap_fraction"]) or 0.0,
            row["seq_id"],
        ),
        reverse=True,
    )
    output: list[dict[str, object]] = []
    for row in high_rows:
        seq_id = row["seq_id"]
        seq_call = sequence_rows_by_id[seq_id]
        output.append(
            {
                "seq_id": seq_id,
                "subtype_source": row["subtype_source"],
                "aligned_length": row["aligned_length"],
                "gap_fraction": row["gap_fraction"],
                "ambiguity_fraction": row["ambiguity_fraction"],
                "gap_or_ambiguity_fraction": row["gap_or_ambiguity_fraction"],
                "core_length_before_alignment": row["core_length_before_alignment"],
                "is_target_or_representative": bool(target_by_alignment.get(seq_id)),
                "is_flagged_rescue": bool(flags_by_alignment.get(seq_id)),
                "is_rescued479": seq_id in rescued_present_by_id,
                "recommended_curation_call": seq_call["round6_curation_call"],
                "reason": seq_call["reason"],
                "notes": seq_call["notes"],
            }
        )
    return output


def subtype_interpretation(
    subtype: str,
    row: dict[str, str],
    rescued_present_count: int,
    rescued_high_count: int,
    total_high_count: int,
) -> tuple[str, str, str]:
    if subtype == "Ib_len255":
        interpretation = (
            "Ib_len255 has materially higher gap/ambiguity burden than Ia or II; "
            "the exact rescued479 set contributes "
            f"{rescued_high_count}/{total_high_count} high-gap rows, so the burden is not mainly "
            "from newly rescued 255-279 aa sequences."
        )
        action = "proceed_to_user_review_before_tree; do not invalidate canonical_min=255 from this evidence alone"
        notes = (
            f"rescued479_present_in_alignment={rescued_present_count}; "
            f"rescued479_high_gap={rescued_high_count}; no gap_fraction_gt_0_70 rows"
        )
    elif subtype == "Ia":
        interpretation = "Ia burden is low; only a small outlier set exceeds 0.50."
        action = "review_high_gap_outliers; otherwise keep_for_tree"
        notes = "No evidence here that Ia broadly drives masked-alignment gap burden."
    elif subtype == "II":
        interpretation = "Type II burden is low overall with two high-gap outliers."
        action = "manual_review_before_tree for Type II outliers; otherwise keep_for_tree"
        notes = "The two Type II outliers should be checked before tree approval."
    else:
        interpretation = "No unknown-subtype sequences are present in the masked alignment."
        action = "no_action"
        notes = "Included for explicit Ia/Ib_len255/II/unknown comparison."
    if row:
        return interpretation, action, notes
    return interpretation, action, notes


def build_subtype_review(
    context: dict[str, object],
    rescued_present_count: int,
    rescued_high_count: int,
    total_high_count: int,
) -> list[dict[str, object]]:
    subtype_rows = context["subtype_rows"]
    assert isinstance(subtype_rows, list)
    subtype_by_name = {row["subtype_source"]: row for row in subtype_rows}
    output: list[dict[str, object]] = []
    for subtype in ["Ia", "Ib_len255", "II", "unknown"]:
        row = subtype_by_name.get(subtype, {})
        interpretation, action, notes = subtype_interpretation(
            subtype, row, rescued_present_count, rescued_high_count, total_high_count
        )
        output.append(
            {
                "subtype_source": subtype,
                "aligned_sequence_count": row.get("aligned_sequence_count", "0"),
                "mean_gap_fraction": row.get("mean_gap_fraction", ""),
                "median_gap_fraction": row.get("median_gap_fraction", ""),
                "max_gap_fraction": row.get("max_gap_fraction", ""),
                "gap_fraction_gt_0_50_count": row.get("gap_fraction_gt_0_50_count", "0"),
                "gap_fraction_gt_0_70_count": row.get("gap_fraction_gt_0_70_count", "0"),
                "mean_gap_or_ambiguity_fraction": row.get("mean_gap_or_ambiguity_fraction", ""),
                "median_gap_or_ambiguity_fraction": row.get(
                    "median_gap_or_ambiguity_fraction", ""
                ),
                "interpretation": interpretation,
                "recommended_action": action,
                "notes": notes,
            }
        )
    return output


def target_curation_call(row: dict[str, str]) -> tuple[str, str, str]:
    label = row["abbrev"]
    status = row["expected_status_from_target_representation"]
    gap_or_ambiguity = as_float(row.get("gap_or_ambiguity_fraction", ""))
    present = as_bool(row.get("present_in_alignment", ""))
    if label == "PniDAH7PS":
        return (
            "needs_curation",
            "user curation required before accepting Pni surrogate",
            "V8CS59 / UniRef90_F9DH16 remains represented by a surrogate but unresolved",
        )
    if label == "legacy A0A0F2JEB6" or status == "unresolved_accession":
        return (
            "needs_external_validation",
            "external accession/sequence validation required",
            "legacy A0A0F2JEB6 remains unresolved_accession and is not represented",
        )
    if status in {"absent_from_input", "hold"} or not present:
        return (
            "hold",
            "not present as an approved tree label in the masked alignment",
            "do not label as a direct final-tree tip",
        )
    if gap_or_ambiguity is not None and gap_or_ambiguity > 0.50:
        return (
            "manual_review_before_tree",
            "target or representative has gap/ambiguity above 0.50",
            "do not remove automatically; review before tree approval",
        )
    return (
        "keep_for_tree",
        "target or representative is present without Round 6 high-gap burden",
        "acceptable for future tree input if user approves curation handoff",
    )


def build_target_review(context: dict[str, object]) -> list[dict[str, object]]:
    target_rows = context["target_rows"]
    assert isinstance(target_rows, list)
    output: list[dict[str, object]] = []
    for row in target_rows:
        call, recommendation, reason = target_curation_call(row)
        output.append(
            {
                "target_label": row["abbrev"],
                "primary_accession": row["primary_accession"],
                "uniref_id": row["uniref_id"],
                "expected_status_from_target_representation": row[
                    "expected_status_from_target_representation"
                ],
                "representative_id": row["representative_id"],
                "alignment_seq_id": row["alignment_seq_id"],
                "present_in_alignment": row["present_in_alignment"],
                "gap_fraction": row["gap_fraction"],
                "gap_or_ambiguity_fraction": row["gap_or_ambiguity_fraction"],
                "label_policy": row["label_policy"],
                "curation_call": call,
                "round6_recommendation": recommendation,
                "reason": reason,
                "notes": row["notes"],
            }
        )
    return output


def build_flagged_review(
    context: dict[str, object], sequence_rows_by_id: dict[str, dict[str, object]]
) -> list[dict[str, object]]:
    flagged_rows = context["flagged_rows"]
    assert isinstance(flagged_rows, list)
    output: list[dict[str, object]] = []
    for row in flagged_rows:
        alignment_seq_id = row.get("alignment_seq_id", "")
        seq_call = sequence_rows_by_id.get(alignment_seq_id, {})
        round6_call = seq_call.get("round6_curation_call", "hold")
        if as_bool(row.get("present_in_alignment", "")):
            round6_call = "manual_review_before_tree"
        output.append(
            {
                "flag_class": row["flag_class"],
                "seq_id": row["seq_id"],
                "nr80_status": row["nr80_status"],
                "nr80_representative": row["nr80_representative"],
                "present_in_core": row["present_in_core"],
                "present_in_alignment": row["present_in_alignment"],
                "alignment_seq_id": row["alignment_seq_id"],
                "gap_fraction": row["gap_fraction"],
                "gap_or_ambiguity_fraction": row["gap_or_ambiguity_fraction"],
                "recommended_action": row["recommended_action"],
                "round6_curation_call": round6_call,
                "reason": "flagged rescue class requires explicit review before final tree",
                "notes": f"{row['notes']} Round 6 applies no silent filter.",
            }
        )
    return output


def find_target_row(rows: list[dict[str, str]], label: str) -> dict[str, str]:
    for row in rows:
        if row.get("abbrev") == label:
            return row
    return {}


def find_target_rep_row(rows: list[dict[str, str]], label: str) -> dict[str, str]:
    for row in rows:
        if row.get("abbrev") == label:
            return row
    return {}


def build_pni_review(context: dict[str, object]) -> list[dict[str, object]]:
    target_rows = context["target_rows"]
    target_rep_rows = context["target_rep_rows"]
    core_target_rows = context["core_target_rows"]
    alignment_by_id = context["alignment_by_id"]

    assert isinstance(target_rows, list)
    assert isinstance(target_rep_rows, list)
    assert isinstance(core_target_rows, list)
    assert isinstance(alignment_by_id, dict)

    pni_msa = find_target_row(target_rows, "PniDAH7PS")
    pni_rep = find_target_rep_row(target_rep_rows, "PniDAH7PS")
    pni_core = find_target_row(core_target_rows, "PniDAH7PS")
    legacy_msa = find_target_row(target_rows, "legacy A0A0F2JEB6")
    legacy_rep = find_target_rep_row(target_rep_rows, "legacy A0A0F2JEB6")
    legacy_core = find_target_row(core_target_rows, "legacy A0A0F2JEB6")

    representative_id = pni_msa.get("representative_id", "UniRef90_A0A379DXQ3")
    rep_alignment = alignment_by_id.get(representative_id, {})
    pni_present_in_core = pni_core.get(
        "present_in_core", pni_core.get("representative_present_in_core", "")
    )
    legacy_present_in_core = legacy_core.get(
        "present_in_core", legacy_core.get("representative_present_in_core", "false")
    )
    coverage_status = (
        "nr80_local_target_coverage="
        f"{pni_rep.get('nr80_local_target_coverage', '')}; "
        "nr80_identity_if_absorbed="
        f"{pni_rep.get('nr80_identity_if_absorbed', '')}"
    )
    evidence = unique_join(
        [
            "results/02_qc_len255/target_representation.tsv",
            "results/03_core_len255/core_len255_target_carrythrough.tsv",
            "results/04_msa_len255/core_len255_target_carrythrough.tsv",
            "results/04_msa_len255/core_len255.alignment_qc.tsv",
        ]
    )

    return [
        {
            "candidate_label": "V8CS59 / UniRef90_F9DH16",
            "primary_accession": "V8CS59",
            "representative_id": representative_id,
            "present_in_target_representation": bool(pni_rep),
            "present_in_core": pni_present_in_core,
            "present_in_alignment": pni_msa.get("present_in_alignment", ""),
            "gap_fraction": pni_msa.get("gap_fraction", ""),
            "gap_or_ambiguity_fraction": pni_msa.get("gap_or_ambiguity_fraction", ""),
            "annotation_status": "needs_curation",
            "taxonomy_status": "需要进一步验证",
            "coverage_status": coverage_status,
            "curation_call": "needs_curation",
            "evidence_source": evidence,
            "notes": "current Pni candidate remains represented_by_nr80_surrogate and is not accepted",
        },
        {
            "candidate_label": "UniRef90_A0A379DXQ3",
            "primary_accession": "A0A379DXQ3",
            "representative_id": representative_id,
            "present_in_target_representation": bool(pni_rep),
            "present_in_core": pni_present_in_core,
            "present_in_alignment": "true" if rep_alignment else "false",
            "gap_fraction": rep_alignment.get("gap_fraction", pni_msa.get("gap_fraction", "")),
            "gap_or_ambiguity_fraction": rep_alignment.get(
                "gap_or_ambiguity_fraction", pni_msa.get("gap_or_ambiguity_fraction", "")
            ),
            "annotation_status": "representative_for_Pni_candidate_needs_curation",
            "taxonomy_status": "需要进一步验证",
            "coverage_status": coverage_status,
            "curation_call": "needs_curation",
            "evidence_source": evidence,
            "notes": "presence as an nr80 representative does not resolve Pni surrogate acceptability",
        },
        {
            "candidate_label": "legacy A0A0F2JEB6",
            "primary_accession": "A0A0F2JEB6",
            "representative_id": legacy_msa.get("representative_id", ""),
            "present_in_target_representation": bool(legacy_rep),
            "present_in_core": legacy_present_in_core,
            "present_in_alignment": legacy_msa.get("present_in_alignment", "false"),
            "gap_fraction": legacy_msa.get("gap_fraction", ""),
            "gap_or_ambiguity_fraction": legacy_msa.get("gap_or_ambiguity_fraction", ""),
            "annotation_status": "unresolved_accession",
            "taxonomy_status": "需要进一步验证",
            "coverage_status": "需要进一步验证",
            "curation_call": "needs_external_validation",
            "evidence_source": evidence,
            "notes": "legacy accession is not rescued, represented, or merged with V8CS59 / UniRef90_F9DH16",
        },
    ]


def build_candidate_lists(
    sequence_rows: list[dict[str, object]]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    exclusions: list[dict[str, object]] = []
    keep_review: list[dict[str, object]] = []
    for row in sequence_rows:
        call = row["round6_curation_call"]
        if call == "candidate_exclude_before_tree":
            exclusions.append(
                {
                    "seq_id": row["seq_id"],
                    "subtype_source": row["subtype_source"],
                    "gap_fraction": row["gap_fraction"],
                    "gap_or_ambiguity_fraction": row["gap_or_ambiguity_fraction"],
                    "reason_for_candidate_exclusion": row["reason"],
                    "would_affect_target_representation": row["is_target_or_representative"],
                    "would_affect_flagged_rescue": row["is_flagged_rescue"],
                    "needs_user_approval_before_removal": "true",
                    "notes": "Advisory only; no removal or filtered alignment was created in Round 6.",
                }
            )
        elif (
            row["is_high_gap_gt_0_50"]
            or row["is_target_or_representative"]
            or row["is_flagged_rescue"]
            or row["is_rescued479"]
        ):
            status_bits = []
            if row["target_label_if_any"]:
                status_bits.append(f"target:{row['target_label_if_any']}")
            if row["flag_class_if_any"]:
                status_bits.append(f"flag:{row['flag_class_if_any']}")
            if row["is_rescued479"]:
                status_bits.append("rescued479:true")
            status_bits.append(f"call:{call}")
            keep_review.append(
                {
                    "seq_id": row["seq_id"],
                    "subtype_source": row["subtype_source"],
                    "gap_fraction": row["gap_fraction"],
                    "gap_or_ambiguity_fraction": row["gap_or_ambiguity_fraction"],
                    "reason_to_keep_or_review": row["reason"],
                    "target_or_flag_status": ";".join(status_bits),
                    "notes": row["notes"],
                }
            )
    exclusions.sort(
        key=lambda row: (as_float(str(row["gap_or_ambiguity_fraction"])) or 0.0, str(row["seq_id"])),
        reverse=True,
    )
    keep_review.sort(
        key=lambda row: (as_float(str(row["gap_or_ambiguity_fraction"])) or 0.0, str(row["seq_id"])),
        reverse=True,
    )
    return exclusions, keep_review


def metric_from_summary(summary_rows: list[dict[str, str]], metric: str) -> str:
    for row in summary_rows:
        if row["metric"] == metric:
            return row["value"]
    return ""


def summarize_rescued479(
    rescued_rows: list[dict[str, str]], sequence_rows_by_id: dict[str, dict[str, object]]
) -> dict[str, object]:
    present_rows = [row for row in rescued_rows if as_bool(row.get("present_in_alignment", ""))]
    gaps = [as_float(row["gap_fraction"]) for row in present_rows if row.get("gap_fraction")]
    gap_values = [gap for gap in gaps if gap is not None]
    calls = [
        sequence_rows_by_id[row["seq_id"]]["round6_curation_call"]
        for row in present_rows
        if row["seq_id"] in sequence_rows_by_id
    ]
    return {
        "total": len(rescued_rows),
        "present_in_nr80": sum(as_bool(row.get("present_in_nr80", "")) for row in rescued_rows),
        "present_in_core": sum(as_bool(row.get("present_in_core", "")) for row in rescued_rows),
        "present_in_alignment": len(present_rows),
        "mean_gap_fraction": statistics.mean(gap_values) if gap_values else None,
        "median_gap_fraction": statistics.median(gap_values) if gap_values else None,
        "max_gap_fraction": max(gap_values) if gap_values else None,
        "high_gap_count": sum(
            as_float(row.get("gap_or_ambiguity_fraction", "")) is not None
            and (as_float(row.get("gap_or_ambiguity_fraction", "")) or 0.0) > 0.50
            for row in present_rows
        ),
        "candidate_exclude_count": calls.count("candidate_exclude_before_tree"),
        "manual_review_count": calls.count("manual_review_before_tree"),
        "keep_count": calls.count("keep_for_tree"),
    }


def target_status_text(target_review_rows: list[dict[str, object]], label: str) -> str:
    for row in target_review_rows:
        if row["target_label"] == label:
            parts = [
                str(row["expected_status_from_target_representation"]),
                f"alignment_seq_id={row['alignment_seq_id'] or 'NA'}",
                f"present_in_alignment={row['present_in_alignment']}",
                f"gap_fraction={row['gap_fraction'] or 'NA'}",
                f"round6_call={row['curation_call']}",
            ]
            return "; ".join(parts)
    return "missing_from_round6_target_review"


def flagged_status_text(flagged_rows: list[dict[str, object]], seq_id: str) -> str:
    for row in flagged_rows:
        if row["seq_id"] == seq_id:
            parts = [
                row["flag_class"],
                f"present_in_alignment={row['present_in_alignment']}",
                f"alignment_seq_id={row['alignment_seq_id']}",
                f"gap_fraction={row['gap_fraction']}",
                f"round6_call={row['round6_curation_call']}",
            ]
            return "; ".join(str(part) for part in parts)
    return "missing_from_round6_flagged_review"


def build_summary(
    context: dict[str, object],
    sequence_rows: list[dict[str, object]],
    target_review_rows: list[dict[str, object]],
    flagged_review_rows: list[dict[str, object]],
    rescued_summary: dict[str, object],
    alignment_sequence_count: int,
    alignment_column_count: int,
) -> list[dict[str, object]]:
    alignment_rows = context["alignment_rows"]
    high_gap_input_rows = context["high_gap_input_rows"]
    summary_rows = context["summary_rows"]
    assert isinstance(alignment_rows, list)
    assert isinstance(high_gap_input_rows, list)
    assert isinstance(summary_rows, list)

    gap_or_amb_gt_050 = sum(
        (as_float(row["gap_or_ambiguity_fraction"]) or 0.0) > 0.50 for row in alignment_rows
    )
    gap_or_amb_gt_070 = sum(
        (as_float(row["gap_or_ambiguity_fraction"]) or 0.0) > 0.70 for row in alignment_rows
    )
    subtype_high = {
        subtype: sum(
            row["subtype_source"] == subtype
            and (as_float(row["gap_or_ambiguity_fraction"]) or 0.0) > 0.50
            for row in alignment_rows
        )
        for subtype in ["Ib_len255", "Ia", "II"]
    }
    call_counts = {
        call: sum(row["round6_curation_call"] == call for row in sequence_rows) for call in CALLS
    }
    notes_high = (
        f"input high-gap file rows={len(high_gap_input_rows)}; "
        f"Round 5 summary gap>0.50={metric_from_summary(summary_rows, 'sequences_gap_fraction_gt_0_50_count')}; "
        f"Round 5 summary gap>0.70={metric_from_summary(summary_rows, 'sequences_gap_fraction_gt_0_70_count')}"
    )
    return [
        {"metric": "alignment_used", "value": ALIGNMENT_AFA, "notes": "Round 5 masked AFA"},
        {"metric": "alignment_is_masked", "value": "true", "notes": "RF-masked alignment"},
        {
            "metric": "alignment_sequence_count",
            "value": alignment_sequence_count,
            "notes": "records parsed from masked AFA",
        },
        {
            "metric": "alignment_column_count",
            "value": alignment_column_count,
            "notes": "columns parsed from masked AFA",
        },
        {
            "metric": "total_sequences_reviewed",
            "value": len(sequence_rows),
            "notes": "all alignment rows received a Round 6 curation call",
        },
        {
            "metric": "gap_or_ambiguity_gt_0_50_count",
            "value": gap_or_amb_gt_050,
            "notes": notes_high,
        },
        {
            "metric": "gap_or_ambiguity_gt_0_70_count",
            "value": gap_or_amb_gt_070,
            "notes": "confirmed from per-sequence alignment QC",
        },
        {
            "metric": "candidate_exclude_before_tree_count",
            "value": call_counts["candidate_exclude_before_tree"],
            "notes": "advisory only; no sequence removals applied",
        },
        {
            "metric": "manual_review_before_tree_count",
            "value": call_counts["manual_review_before_tree"],
            "notes": "manual review rows in sequence-level calls",
        },
        {
            "metric": "keep_for_tree_count",
            "value": call_counts["keep_for_tree"],
            "notes": "kept for future tree input if user approves curation",
        },
        {
            "metric": "Ib_len255_high_gap_count",
            "value": subtype_high["Ib_len255"],
            "notes": "Ib_len255 rows with gap_or_ambiguity_fraction > 0.50",
        },
        {
            "metric": "Ia_high_gap_count",
            "value": subtype_high["Ia"],
            "notes": "Ia rows with gap_or_ambiguity_fraction > 0.50",
        },
        {
            "metric": "II_high_gap_count",
            "value": subtype_high["II"],
            "notes": "II rows with gap_or_ambiguity_fraction > 0.50",
        },
        {
            "metric": "rescued479_present_in_alignment_count",
            "value": rescued_summary["present_in_alignment"],
            "notes": "exact rescued set from Round 5 rescued479 table",
        },
        {
            "metric": "rescued479_high_gap_count",
            "value": rescued_summary["high_gap_count"],
            "notes": "present rescued479 rows with gap_or_ambiguity_fraction > 0.50",
        },
        {
            "metric": "Q8U0A9_representative_status",
            "value": target_status_text(target_review_rows, "PfuDAH7PS"),
            "notes": "surrogate label policy remains PfuDAH7PS-like / represents Q8U0A9",
        },
        {
            "metric": "Q9YEJ7_status",
            "value": target_status_text(target_review_rows, "ApeDAH7PS"),
            "notes": "direct nr80 tip, but high gap triggers manual review",
        },
        {
            "metric": "Pni_candidate_representative_status",
            "value": target_status_text(target_review_rows, "PniDAH7PS"),
            "notes": "V8CS59 / UniRef90_F9DH16 remains needs_curation",
        },
        {
            "metric": "legacy_A0A0F2JEB6_status",
            "value": "unresolved_accession",
            "notes": "not rescued, represented, solved, or merged with current Pni candidate",
        },
        {
            "metric": "KDOPS_like_flagged_hit_status",
            "value": flagged_status_text(flagged_review_rows, "UniRef90_A0A5E4LPI9"),
            "notes": "manual review before final tree; not called clean DAH7PS",
        },
        {
            "metric": "round6_tree_readiness_call",
            "value": "proceed_to_user_review_before_tree",
            "notes": "IQ-TREE should wait for user approval of advisory curation decisions",
        },
        {
            "metric": "rescued479_mean_gap_fraction",
            "value": format_float(rescued_summary["mean_gap_fraction"]),
            "notes": "present rescued479 alignment rows only",
        },
        {
            "metric": "rescued479_median_gap_fraction",
            "value": format_float(rescued_summary["median_gap_fraction"]),
            "notes": "present rescued479 alignment rows only",
        },
        {
            "metric": "rescued479_max_gap_fraction",
            "value": format_float(rescued_summary["max_gap_fraction"]),
            "notes": "present rescued479 alignment rows only",
        },
    ]


def render_handoff(
    summary_rows: list[dict[str, object]],
    subtype_rows: list[dict[str, object]],
    target_review_rows: list[dict[str, object]],
    flagged_review_rows: list[dict[str, object]],
    candidate_exclusions: list[dict[str, object]],
    sequence_rows: list[dict[str, object]],
    rescued_summary: dict[str, object],
) -> str:
    summary = {str(row["metric"]): str(row["value"]) for row in summary_rows}
    subtype_by_name = {str(row["subtype_source"]): row for row in subtype_rows}
    manual_rows = [row for row in sequence_rows if row["round6_curation_call"] == "manual_review_before_tree"]
    manual_rows.sort(
        key=lambda row: (as_float(str(row["gap_or_ambiguity_fraction"])) or 0.0, str(row["seq_id"])),
        reverse=True,
    )
    top_manual = "\n".join(
        f"- `{row['seq_id']}` ({row['subtype_source']}; gap_or_ambiguity={row['gap_or_ambiguity_fraction']}; "
        f"{row['reason']})"
        for row in manual_rows[:10]
    )
    if not top_manual:
        top_manual = "- None."
    exclusion_text = "\n".join(
        f"- `{row['seq_id']}` ({row['subtype_source']}; gap_or_ambiguity={row['gap_or_ambiguity_fraction']})"
        for row in candidate_exclusions
    )
    if not exclusion_text:
        exclusion_text = "- None under the Round 6 evidence-based advisory rules."

    pfu = target_status_text(target_review_rows, "PfuDAH7PS")
    ape = target_status_text(target_review_rows, "ApeDAH7PS")
    pni = target_status_text(target_review_rows, "PniDAH7PS")
    kdops = flagged_status_text(flagged_review_rows, "UniRef90_A0A5E4LPI9")

    return f"""# Strategy A Round 6 Tree Readiness Handoff

## 1. Future Alignment

A future tree round would use:

```text
{ALIGNMENT_AFA}
```

This is the Round 5 RF-masked AFA alignment. Round 6 did not create a filtered
alignment and did not remove sequences.

## 2. Alignment Dimensions

- Sequences: {summary['alignment_sequence_count']}
- Columns: {summary['alignment_column_count']}
- Masking: masked

## 3. High-Gap / Ambiguity Burden

- gap_or_ambiguity_fraction > 0.50: {summary['gap_or_ambiguity_gt_0_50_count']}
- gap_or_ambiguity_fraction > 0.70: {summary['gap_or_ambiguity_gt_0_70_count']}
- candidate_exclude_before_tree: {summary['candidate_exclude_before_tree_count']}
- manual_review_before_tree: {summary['manual_review_before_tree_count']}
- keep_for_tree: {summary['keep_for_tree_count']}

Round 5 reported 1,822 sequences above 0.50 gap/ambiguity and 0 sequences above
0.70 gap fraction; Round 6 confirms those counts from the formal tables.

## 4. Subtype Burden

- Ia high-gap count: {summary['Ia_high_gap_count']}
- Ib_len255 high-gap count: {summary['Ib_len255_high_gap_count']}
- II high-gap count: {summary['II_high_gap_count']}
- unknown high-gap count: 0

Ib_len255 has materially higher gap/ambiguity burden than Ia or II. The exact
rescued479 set contributes {rescued_summary['high_gap_count']} high-gap rows among
{summary['gap_or_ambiguity_gt_0_50_count']} total high-gap rows, so the burden is
not mainly from newly rescued 255-279 aa sequences. This argues for user curation
review and cautious tree approval, not for invalidating Strategy A
`canonical_min=255` from Round 6 evidence alone.

Subtype details are in:

```text
results/05_curation_len255/round6_subtype_gap_review.tsv
```

## 5. Target Status

- Q8U0A9 / PfuDAH7PS: {pfu}
- Q9YEJ7 / ApeDAH7PS: {ape}
- V8CS59 / UniRef90_F9DH16 / PniDAH7PS candidate: {pni}
- legacy A0A0F2JEB6: unresolved_accession; not rescued, represented, solved, or merged with V8CS59 / UniRef90_F9DH16.
- KDOPS-like `UniRef90_A0A5E4LPI9`: {kdops}

The Pni surrogate remains `needs_curation`. Presence of `UniRef90_A0A379DXQ3` in
the MSA is not sufficient to call it an acceptable Pni surrogate.

## 6. Candidate Exclusions

Candidate exclusions before tree inference:

{exclusion_text}

The advisory table is:

```text
results/05_curation_len255/round6_candidate_exclusion_list.tsv
```

No sequence removal was applied.

## 7. Manual Review

Manual review before tree inference is required for {summary['manual_review_before_tree_count']}
sequence-level calls. Top severity examples:

{top_manual}

The full advisory keep/review table is:

```text
results/05_curation_len255/round6_candidate_keep_review_list.tsv
```

## 8. Recommendation

Do not run IQ-TREE yet. The next step is user approval of the Round 6 advisory
curation calls, especially the high-gap manual-review rows, the Pni surrogate,
and the KDOPS-like flagged hit. The Round 6 readiness call is:

```text
{summary['round6_tree_readiness_call']}
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
"""


def build_artifact_manifest() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    output_source = ";".join(ROUND5_INPUTS)
    paths = [(path, "round6_output") for path in GENERATED_OUTPUTS]
    paths.append((SCRIPT_PATH, "round6_helper_script"))
    paths.extend((path, "round5_or_source_input") for path in ROUND5_INPUTS)

    for relative_path, role in paths:
        path = ROOT / relative_path
        if path.exists() and relative_path != "results/05_curation_len255/round6_artifact_manifest.tsv":
            size = path.stat().st_size
            md5 = md5_file(path)
        elif relative_path == "results/05_curation_len255/round6_artifact_manifest.tsv":
            size = path.stat().st_size if path.exists() else ""
            md5 = "SELF_REFERENTIAL"
        else:
            size = ""
            md5 = "MISSING"
        if role == "round6_output":
            required = "true"
            source_inputs = output_source
            recreate = "python scripts/round6_pre_tree_curation_len255.py"
            notes = "Round 6 curation output; advisory only"
        elif role == "round6_helper_script":
            required = "true"
            source_inputs = output_source
            recreate = "tracked helper script"
            notes = "Deterministic; does not run tree inference"
        else:
            required = "true"
            source_inputs = ""
            recreate = "source artifact from Round 5 or earlier Strategy A stage"
            notes = "Input used by Round 6 curation gate"
        if relative_path == "results/05_curation_len255/round6_artifact_manifest.tsv":
            notes = "Self-referential manifest row; md5 cannot include its own final md5"
        rows.append(
            {
                "relative_path": relative_path,
                "artifact_role": role,
                "file_size_bytes": size,
                "md5": md5,
                "tracked_by_git": is_tracked(relative_path),
                "ignored_by_git": is_ignored(relative_path),
                "required_for_next_stage": required,
                "source_inputs": source_inputs,
                "recreate_command_or_source": recreate,
                "notes": notes,
            }
        )
    return rows


def main() -> int:
    for relative_path in ROUND5_INPUTS:
        if not (ROOT / relative_path).is_file():
            raise FileNotFoundError(relative_path)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    context = load_context()
    alignment_sequence_count, alignment_column_count = parse_alignment_dimensions(ALIGNMENT_AFA)

    sequence_rows = build_sequence_rows(context)
    sequence_rows_by_id = {str(row["seq_id"]): row for row in sequence_rows}
    high_gap_rows = build_high_gap_review(context, sequence_rows_by_id)
    rescued_rows = context["rescued_rows"]
    assert isinstance(rescued_rows, list)
    rescued_summary = summarize_rescued479(rescued_rows, sequence_rows_by_id)
    subtype_rows = build_subtype_review(
        context,
        int(rescued_summary["present_in_alignment"]),
        int(rescued_summary["high_gap_count"]),
        len(high_gap_rows),
    )
    target_review_rows = build_target_review(context)
    flagged_review_rows = build_flagged_review(context, sequence_rows_by_id)
    pni_review_rows = build_pni_review(context)
    candidate_exclusions, candidate_keep_review = build_candidate_lists(sequence_rows)
    summary_rows = build_summary(
        context,
        sequence_rows,
        target_review_rows,
        flagged_review_rows,
        rescued_summary,
        alignment_sequence_count,
        alignment_column_count,
    )

    write_tsv(
        "results/05_curation_len255/round6_sequence_curation_calls.tsv",
        [
            "seq_id",
            "subtype_source",
            "gap_fraction",
            "ambiguity_fraction",
            "gap_or_ambiguity_fraction",
            "is_high_gap_gt_0_50",
            "is_high_gap_gt_0_70",
            "is_target_or_representative",
            "target_label_if_any",
            "is_flagged_rescue",
            "flag_class_if_any",
            "is_rescued479",
            "round6_curation_call",
            "reason",
            "notes",
        ],
        sequence_rows,
    )
    write_tsv(
        "results/05_curation_len255/round6_high_gap_review.tsv",
        [
            "seq_id",
            "subtype_source",
            "aligned_length",
            "gap_fraction",
            "ambiguity_fraction",
            "gap_or_ambiguity_fraction",
            "core_length_before_alignment",
            "is_target_or_representative",
            "is_flagged_rescue",
            "is_rescued479",
            "recommended_curation_call",
            "reason",
            "notes",
        ],
        high_gap_rows,
    )
    write_tsv(
        "results/05_curation_len255/round6_subtype_gap_review.tsv",
        [
            "subtype_source",
            "aligned_sequence_count",
            "mean_gap_fraction",
            "median_gap_fraction",
            "max_gap_fraction",
            "gap_fraction_gt_0_50_count",
            "gap_fraction_gt_0_70_count",
            "mean_gap_or_ambiguity_fraction",
            "median_gap_or_ambiguity_fraction",
            "interpretation",
            "recommended_action",
            "notes",
        ],
        subtype_rows,
    )
    write_tsv(
        "results/05_curation_len255/round6_target_curation_review.tsv",
        [
            "target_label",
            "primary_accession",
            "uniref_id",
            "expected_status_from_target_representation",
            "representative_id",
            "alignment_seq_id",
            "present_in_alignment",
            "gap_fraction",
            "gap_or_ambiguity_fraction",
            "label_policy",
            "curation_call",
            "round6_recommendation",
            "reason",
            "notes",
        ],
        target_review_rows,
    )
    write_tsv(
        "results/05_curation_len255/round6_flagged_rescue_review.tsv",
        [
            "flag_class",
            "seq_id",
            "nr80_status",
            "nr80_representative",
            "present_in_core",
            "present_in_alignment",
            "alignment_seq_id",
            "gap_fraction",
            "gap_or_ambiguity_fraction",
            "recommended_action",
            "round6_curation_call",
            "reason",
            "notes",
        ],
        flagged_review_rows,
    )
    write_tsv(
        "results/05_curation_len255/round6_pni_surrogate_review.tsv",
        [
            "candidate_label",
            "primary_accession",
            "representative_id",
            "present_in_target_representation",
            "present_in_core",
            "present_in_alignment",
            "gap_fraction",
            "gap_or_ambiguity_fraction",
            "annotation_status",
            "taxonomy_status",
            "coverage_status",
            "curation_call",
            "evidence_source",
            "notes",
        ],
        pni_review_rows,
    )
    write_tsv(
        "results/05_curation_len255/round6_candidate_exclusion_list.tsv",
        [
            "seq_id",
            "subtype_source",
            "gap_fraction",
            "gap_or_ambiguity_fraction",
            "reason_for_candidate_exclusion",
            "would_affect_target_representation",
            "would_affect_flagged_rescue",
            "needs_user_approval_before_removal",
            "notes",
        ],
        candidate_exclusions,
    )
    write_tsv(
        "results/05_curation_len255/round6_candidate_keep_review_list.tsv",
        [
            "seq_id",
            "subtype_source",
            "gap_fraction",
            "gap_or_ambiguity_fraction",
            "reason_to_keep_or_review",
            "target_or_flag_status",
            "notes",
        ],
        candidate_keep_review,
    )
    write_tsv(
        "results/05_curation_len255/round6_curation_summary.tsv",
        ["metric", "value", "notes"],
        summary_rows,
    )
    handoff = render_handoff(
        summary_rows,
        subtype_rows,
        target_review_rows,
        flagged_review_rows,
        candidate_exclusions,
        sequence_rows,
        rescued_summary,
    )
    (OUT_DIR / "round6_tree_readiness_handoff.md").write_text(handoff, encoding="utf-8")

    manifest_rows = build_artifact_manifest()
    write_tsv(
        "results/05_curation_len255/round6_artifact_manifest.tsv",
        [
            "relative_path",
            "artifact_role",
            "file_size_bytes",
            "md5",
            "tracked_by_git",
            "ignored_by_git",
            "required_for_next_stage",
            "source_inputs",
            "recreate_command_or_source",
            "notes",
        ],
        manifest_rows,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
