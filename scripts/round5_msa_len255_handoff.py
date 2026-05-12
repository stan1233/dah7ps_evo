#!/usr/bin/env python3
"""Generate Strategy A Round 5 MSA QC and tree-handoff artifacts."""

from __future__ import annotations

import csv
import hashlib
import statistics
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/04_msa_len255"

CORE_FASTA = ROOT / "results/03_core_len255/all_core_only_len255.fasta"
CORE_HMM = ROOT / "results/03_msa_core/core_global.hmm"
UNMASKED_STO = OUTDIR / "core_len255.sto"
MASKED_STO = OUTDIR / "core_len255.masked.sto"
MASKED_AFA = OUTDIR / "core_len255.masked.afa"

SUBTYPE_MAP = ROOT / "results/02_qc_len255/nr80_all_len255_subtype_map.tsv"
TARGET_REPRESENTATION = ROOT / "results/02_qc_len255/target_representation.tsv"
CORE_TARGET_CARRY = ROOT / "results/03_core_len255/core_len255_target_carrythrough.tsv"
CORE_FLAGGED_CARRY = ROOT / "results/03_core_len255/core_len255_flagged_carrythrough.tsv"
CORE_SUBTYPE_CARRY = ROOT / "results/03_core_len255/core_len255_subtype_carrythrough.tsv"
OLD_IB_QC = ROOT / "results/02_qc/qc_classification_Ib.tsv"
LEN255_IB_QC = ROOT / "results/02_qc_len255/qc_classification_Ib.tsv"
NR80_IB_CLSTR = ROOT / "results/02_qc_len255/nr80_Ib.fasta.clstr"
FLAGGED_ADJUDICATION = ROOT / "results/02_qc_len255/len255_flagged_rescue_adjudication.tsv"

SUMMARY_TSV = OUTDIR / "core_len255_alignment_summary.tsv"
ALIGNMENT_QC_TSV = OUTDIR / "core_len255.alignment_qc.tsv"
SUBTYPE_QC_TSV = OUTDIR / "core_len255_subtype_alignment_qc.tsv"
RESCUED479_QC_TSV = OUTDIR / "core_len255_rescued479_alignment_qc.tsv"
TARGET_CARRY_TSV = OUTDIR / "core_len255_target_carrythrough.tsv"
FLAGGED_CARRY_TSV = OUTDIR / "core_len255_flagged_carrythrough.tsv"
HIGH_GAP_TSV = OUTDIR / "core_len255_gap_or_ambiguity_gt_0_50.tsv"
ARTIFACT_MANIFEST = OUTDIR / "core_len255_artifact_manifest.tsv"
TREE_HANDOFF_MD = OUTDIR / "round5_tree_handoff.md"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
GAP_CHARS = set("-.")
ALLOWED_SUBTYPES = {"Ia", "Ib_len255", "II"}

HMMALIGN_CMD = (
    "conda run -n dah7ps_v4 hmmalign --amino --outformat Stockholm "
    "-o results/04_msa_len255/core_len255.sto "
    "results/03_msa_core/core_global.hmm "
    "results/03_core_len255/all_core_only_len255.fasta"
)
ALIMASK_CMD = (
    "conda run -n dah7ps_v4 esl-alimask --rf-is-mask "
    "-o results/04_msa_len255/core_len255.masked.sto "
    "results/04_msa_len255/core_len255.sto"
)
REFORMAT_CMD = (
    "conda run -n dah7ps_v4 esl-reformat "
    "-o results/04_msa_len255/core_len255.masked.afa "
    "afa results/04_msa_len255/core_len255.masked.sto"
)


@dataclass(frozen=True)
class ClusterHit:
    representative: str
    is_representative: bool


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def require_file(path: Path) -> None:
    if not path.is_file():
        raise SystemExit(f"required file missing: {path}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    require_file(path)
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    if path.parent != OUTDIR:
        raise SystemExit(f"refusing to write outside {rel(OUTDIR)}: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=columns,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            output_row = {column: row.get(column, "") for column in columns}
            if output_row.get(columns[-1], "") == "":
                output_row[columns[-1]] = "NA"
            writer.writerow(output_row)


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    require_file(path)
    records: list[tuple[str, str]] = []
    current_id = ""
    seq_parts: list[str] = []

    def flush() -> None:
        if current_id:
            records.append((current_id, "".join(seq_parts)))

    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                flush()
                current_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    flush()
    return records


def read_fasta(path: Path) -> dict[str, str]:
    return dict(read_fasta_records(path))


def duplicate_ids_from_records(records: list[tuple[str, str]]) -> list[str]:
    counts = Counter(seq_id for seq_id, _ in records)
    return sorted(seq_id for seq_id, count in counts.items() if count > 1)


def read_stockholm(path: Path) -> dict[str, str]:
    require_file(path)
    fragments: dict[str, list[str]] = defaultdict(list)
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or line == "//":
                continue
            parts = line.split()
            if len(parts) >= 2:
                fragments[parts[0]].append(parts[1])
    return {seq_id: "".join(parts) for seq_id, parts in fragments.items()}


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_bool(*args: str) -> bool:
    return (
        subprocess.run(
            ["git", *args],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        ).returncode
        == 0
    )


def git_tracked(path: Path) -> str:
    return str(git_bool("ls-files", "--error-unmatch", rel(path))).lower()


def git_ignored(path: Path) -> str:
    return str(git_bool("check-ignore", "-q", rel(path))).lower()


def fmt(value: object) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def mean(values: list[float]) -> float:
    return statistics.mean(values) if values else 0.0


def median(values: list[float]) -> float:
    return statistics.median(values) if values else 0.0


def max_or_zero(values: list[float]) -> float:
    return max(values) if values else 0.0


def load_subtypes() -> dict[str, str]:
    subtypes: dict[str, str] = {}
    for row in read_tsv(SUBTYPE_MAP):
        subtype = row.get("subtype_source", "")
        subtypes[row["seq_id"]] = subtype if subtype in ALLOWED_SUBTYPES else "unknown"
    return subtypes


def compute_alignment_qc(
    alignment_records: list[tuple[str, str]],
    core_records: dict[str, str],
    subtypes: dict[str, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for seq_id, sequence in alignment_records:
        aligned_length = len(sequence)
        gap_count = 0
        ambiguity_count = 0
        standard_count = 0
        for char in sequence:
            upper = char.upper()
            if char in GAP_CHARS:
                gap_count += 1
            elif upper in STANDARD_AA:
                standard_count += 1
            else:
                ambiguity_count += 1
        subtype = subtypes.get(seq_id, "unknown")
        notes = "" if subtype != "unknown" else "subtype_source not found in nr80 map"
        rows.append(
            {
                "seq_id": seq_id,
                "subtype_source": subtype,
                "aligned_length": aligned_length,
                "non_gap_non_ambiguity_count": standard_count,
                "gap_count": gap_count,
                "ambiguity_count": ambiguity_count,
                "gap_fraction": gap_count / aligned_length if aligned_length else 0.0,
                "ambiguity_fraction": ambiguity_count / aligned_length
                if aligned_length
                else 0.0,
                "gap_or_ambiguity_fraction": (gap_count + ambiguity_count)
                / aligned_length
                if aligned_length
                else 0.0,
                "core_length_before_alignment": len(core_records.get(seq_id, "")),
                "notes": notes,
            }
        )
    return rows


def write_alignment_qc(rows: list[dict[str, object]]) -> None:
    columns = [
        "seq_id",
        "subtype_source",
        "aligned_length",
        "non_gap_non_ambiguity_count",
        "gap_count",
        "ambiguity_count",
        "gap_fraction",
        "ambiguity_fraction",
        "gap_or_ambiguity_fraction",
        "core_length_before_alignment",
        "notes",
    ]
    write_tsv(
        ALIGNMENT_QC_TSV,
        columns,
        [{key: fmt(value) for key, value in row.items()} for row in rows],
    )


def write_summary(
    core_records: dict[str, str],
    alignment_records: list[tuple[str, str]],
    unmasked_records: dict[str, str],
    qc_rows: list[dict[str, object]],
) -> None:
    duplicates = duplicate_ids_from_records(alignment_records)
    aligned_count = len(alignment_records)
    core_input_count = len(core_records)
    masked_columns = len(alignment_records[0][1]) if alignment_records else 0
    unmasked_columns = len(next(iter(unmasked_records.values()))) if unmasked_records else 0
    gap_fracs = [float(row["gap_fraction"]) for row in qc_rows]
    empty_aligned = sum(
        1
        for row in qc_rows
        if int(row["non_gap_non_ambiguity_count"]) + int(row["ambiguity_count"]) == 0
    )
    rows = [
        {
            "metric": "core_input_count",
            "value": core_input_count,
            "notes": "records in results/03_core_len255/all_core_only_len255.fasta",
        },
        {
            "metric": "aligned_sequence_count",
            "value": aligned_count,
            "notes": "records in results/04_msa_len255/core_len255.masked.afa",
        },
        {
            "metric": "alignment_columns_unmasked",
            "value": unmasked_columns,
            "notes": "columns in results/04_msa_len255/core_len255.sto",
        },
        {
            "metric": "alignment_columns_masked",
            "value": masked_columns,
            "notes": "RF-mask columns in results/04_msa_len255/core_len255.masked.afa",
        },
        {
            "metric": "sequence_count_matches_core_input",
            "value": aligned_count == core_input_count,
            "notes": "true if aligned sequence count equals core FASTA count",
        },
        {
            "metric": "duplicate_ids_count",
            "value": len(duplicates),
            "notes": "unique sequence IDs repeated in masked AFA",
        },
        {
            "metric": "duplicate_ids",
            "value": ",".join(duplicates),
            "notes": "comma-separated; empty means none",
        },
        {
            "metric": "empty_aligned_sequences_count",
            "value": empty_aligned,
            "notes": "aligned records with no non-gap residues",
        },
        {
            "metric": "mean_gap_fraction",
            "value": mean(gap_fracs),
            "notes": "mean per-sequence gap fraction on masked AFA",
        },
        {
            "metric": "median_gap_fraction",
            "value": median(gap_fracs),
            "notes": "median per-sequence gap fraction on masked AFA",
        },
        {
            "metric": "max_gap_fraction",
            "value": max_or_zero(gap_fracs),
            "notes": "maximum per-sequence gap fraction on masked AFA",
        },
        {
            "metric": "sequences_gap_fraction_gt_0_50_count",
            "value": sum(1 for value in gap_fracs if value > 0.50),
            "notes": "masked AFA records with gap_fraction > 0.50",
        },
        {
            "metric": "sequences_gap_fraction_gt_0_70_count",
            "value": sum(1 for value in gap_fracs if value > 0.70),
            "notes": "masked AFA records with gap_fraction > 0.70",
        },
    ]
    write_tsv(
        SUMMARY_TSV,
        ["metric", "value", "notes"],
        [
            {"metric": row["metric"], "value": fmt(row["value"]), "notes": row["notes"]}
            for row in rows
        ],
    )


def write_subtype_qc(qc_rows: list[dict[str, object]]) -> None:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in qc_rows:
        grouped[str(row["subtype_source"])].append(row)
    subtype_means = {
        subtype: mean([float(row["gap_or_ambiguity_fraction"]) for row in rows])
        for subtype, rows in grouped.items()
    }
    ia_ii_mean = mean(
        [subtype_means.get("Ia", 0.0), subtype_means.get("II", 0.0)]
    )
    ib_delta = subtype_means.get("Ib_len255", 0.0) - ia_ii_mean
    ib_note = (
        "Ib_len255 has materially higher gap/ambiguity burden than Ia/II mean."
        if ib_delta > 0.05
        else "Ib_len255 gap/ambiguity burden is not materially higher than Ia/II mean."
    )
    rows_out: list[dict[str, object]] = []
    for subtype in ["Ia", "Ib_len255", "II", "unknown"]:
        rows = grouped.get(subtype, [])
        if not rows and subtype == "unknown":
            continue
        gap_fracs = [float(row["gap_fraction"]) for row in rows]
        gap_amb_fracs = [float(row["gap_or_ambiguity_fraction"]) for row in rows]
        rows_out.append(
            {
                "subtype_source": subtype,
                "aligned_sequence_count": len(rows),
                "mean_gap_fraction": mean(gap_fracs),
                "median_gap_fraction": median(gap_fracs),
                "max_gap_fraction": max_or_zero(gap_fracs),
                "gap_fraction_gt_0_50_count": sum(
                    1 for value in gap_fracs if value > 0.50
                ),
                "gap_fraction_gt_0_70_count": sum(
                    1 for value in gap_fracs if value > 0.70
                ),
                "mean_gap_or_ambiguity_fraction": mean(gap_amb_fracs),
                "median_gap_or_ambiguity_fraction": median(gap_amb_fracs),
                "notes": ib_note if subtype == "Ib_len255" else "",
            }
        )
    columns = [
        "subtype_source",
        "aligned_sequence_count",
        "mean_gap_fraction",
        "median_gap_fraction",
        "max_gap_fraction",
        "gap_fraction_gt_0_50_count",
        "gap_fraction_gt_0_70_count",
        "mean_gap_or_ambiguity_fraction",
        "median_gap_or_ambiguity_fraction",
        "notes",
    ]
    write_tsv(
        SUBTYPE_QC_TSV,
        columns,
        [{key: fmt(value) for key, value in row.items()} for row in rows_out],
    )


def parse_cdhit_clusters(path: Path) -> dict[str, ClusterHit]:
    require_file(path)
    hits: dict[str, ClusterHit] = {}
    current_members: list[dict[str, object]] = []

    def flush() -> None:
        if not current_members:
            return
        representatives = [
            str(member["seq_id"]) for member in current_members if member["is_rep"]
        ]
        representative = representatives[0] if representatives else ""
        for member in current_members:
            hits[str(member["seq_id"])] = ClusterHit(
                representative=representative,
                is_representative=bool(member["is_rep"]),
            )

    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip()
            if line.startswith(">Cluster "):
                flush()
                current_members = []
                continue
            if "\t" not in line or ">" not in line or "..." not in line:
                continue
            seq_id = line.split(">", 1)[1].split("...", 1)[0]
            current_members.append({"seq_id": seq_id, "is_rep": line.endswith("*")})
    flush()
    return hits


def reconstruct_rescued_rows() -> list[dict[str, str]]:
    old_rows = {row["seq_id"]: row for row in read_tsv(OLD_IB_QC)}
    rescued: list[dict[str, str]] = []
    for row in read_tsv(LEN255_IB_QC):
        seq_id = row["seq_id"]
        old = old_rows.get(seq_id, {})
        length = int(row.get("L_seq") or 0)
        if (
            old.get("bin") == "FRAG"
            and row.get("bin") == "PASS_CANONICAL"
            and 255 <= length <= 279
        ):
            rescued.append(row)
    return rescued


def write_rescued479_qc(
    alignment_qc: dict[str, dict[str, str]],
    core_records: dict[str, str],
    nr80_ids: set[str],
) -> None:
    cluster_hits = parse_cdhit_clusters(NR80_IB_CLSTR)
    flagged = {row["seq_id"]: row for row in read_tsv(FLAGGED_ADJUDICATION)}
    rows_out: list[dict[str, object]] = []
    rescued_rows = reconstruct_rescued_rows()
    exact_note = (
        "exact rescued set reconstructed from old/new Ib QC tables"
        if len(rescued_rows) == 479
        else f"rescued set reconstruction count={len(rescued_rows)}; expected 479"
    )
    for row in rescued_rows:
        seq_id = row["seq_id"]
        hit = cluster_hits.get(seq_id)
        if hit is None:
            nr80_status = "absent_from_nr80"
        elif hit.is_representative:
            nr80_status = "direct_nr80_representative"
        else:
            nr80_status = "absorbed_by_nr80_surrogate"
        flagged_row = flagged.get(seq_id, {})
        qc = alignment_qc.get(seq_id, {})
        rows_out.append(
            {
                "seq_id": seq_id,
                "length": row.get("L_seq", ""),
                "cov_best": row.get("cov_best", ""),
                "present_in_nr80": seq_id in nr80_ids,
                "nr80_status": nr80_status,
                "present_in_core": seq_id in core_records,
                "present_in_alignment": seq_id in alignment_qc,
                "gap_fraction": qc.get("gap_fraction", ""),
                "ambiguity_fraction": qc.get("ambiguity_fraction", ""),
                "gap_or_ambiguity_fraction": qc.get(
                    "gap_or_ambiguity_fraction", ""
                ),
                "flag_class": flagged_row.get("flag_class", ""),
                "recommended_action": flagged_row.get("recommended_action", ""),
                "notes": exact_note,
            }
        )
    columns = [
        "seq_id",
        "length",
        "cov_best",
        "present_in_nr80",
        "nr80_status",
        "present_in_core",
        "present_in_alignment",
        "gap_fraction",
        "ambiguity_fraction",
        "gap_or_ambiguity_fraction",
        "flag_class",
        "recommended_action",
        "notes",
    ]
    write_tsv(
        RESCUED479_QC_TSV,
        columns,
        [{key: fmt(value) for key, value in row.items()} for row in rows_out],
    )


def target_label_policy(core_row: dict[str, str]) -> str:
    abbrev = core_row.get("abbrev", "")
    if abbrev == "PfuDAH7PS":
        return "PfuDAH7PS-like / represents Q8U0A9"
    if abbrev == "ApeDAH7PS":
        return "label_direct_nr80_tip"
    return core_row.get("label_policy", "")


def write_target_carrythrough(alignment_qc: dict[str, dict[str, str]]) -> None:
    target_rows = {row["abbrev"]: row for row in read_tsv(TARGET_REPRESENTATION)}
    rows_out: list[dict[str, object]] = []
    missing_required: list[str] = []
    for row in read_tsv(CORE_TARGET_CARRY):
        abbrev = row["abbrev"]
        label_id = row.get("core_sequence_id_used_for_label", "")
        present = bool(label_id and label_id in alignment_qc)
        qc = alignment_qc.get(label_id, {})
        curation_call = row.get("curation_call", "") or target_rows.get(
            abbrev, {}
        ).get("curation_call", "")
        if abbrev == "PniDAH7PS":
            curation_call = "needs_curation"
        if abbrev in {"PfuDAH7PS", "ApeDAH7PS", "PniDAH7PS"} and not present:
            missing_required.append(f"{abbrev}:{label_id}")
        notes = row.get("notes", "")
        if abbrev == "PfuDAH7PS":
            notes = (
                f"{notes} Representative remains PfuDAH7PS-like / represents Q8U0A9."
            )
        if abbrev == "PniDAH7PS":
            notes = f"{notes} Pni surrogate remains needs_curation."
        rows_out.append(
            {
                "abbrev": abbrev,
                "primary_accession": row.get("primary_accession", ""),
                "uniref_id": row.get("uniref_id", ""),
                "expected_status_from_target_representation": row.get(
                    "expected_status_from_target_representation", ""
                ),
                "representative_id": row.get("representative_id", ""),
                "core_sequence_id_used_for_label": label_id,
                "present_in_core": row.get("representative_present_in_core", ""),
                "present_in_alignment": present,
                "alignment_seq_id": label_id if present else "",
                "gap_fraction": qc.get("gap_fraction", ""),
                "gap_or_ambiguity_fraction": qc.get(
                    "gap_or_ambiguity_fraction", ""
                ),
                "label_policy": target_label_policy(row),
                "curation_call": curation_call,
                "notes": notes,
            }
        )
    columns = [
        "abbrev",
        "primary_accession",
        "uniref_id",
        "expected_status_from_target_representation",
        "representative_id",
        "core_sequence_id_used_for_label",
        "present_in_core",
        "present_in_alignment",
        "alignment_seq_id",
        "gap_fraction",
        "gap_or_ambiguity_fraction",
        "label_policy",
        "curation_call",
        "notes",
    ]
    write_tsv(
        TARGET_CARRY_TSV,
        columns,
        [{key: fmt(value) for key, value in row.items()} for row in rows_out],
    )
    if missing_required:
        raise SystemExit(
            "required target IDs missing from alignment: " + ",".join(missing_required)
        )


def write_flagged_carrythrough(alignment_qc: dict[str, dict[str, str]]) -> None:
    rows_out: list[dict[str, object]] = []
    for row in read_tsv(CORE_FLAGGED_CARRY):
        seq_id = row.get("core_sequence_id", "") or row.get("nr80_representative", "")
        present = bool(seq_id and seq_id in alignment_qc)
        qc = alignment_qc.get(seq_id, {})
        gap_amb = float(qc.get("gap_or_ambiguity_fraction") or 0.0)
        notes = row.get("notes", "")
        if present:
            notes = f"{notes} Present in masked alignment."
        if gap_amb > 0.50:
            notes = f"{notes} gap_or_ambiguity_fraction exceeds 0.50."
        action = row.get("recommended_action", "")
        if row.get("seq_id") == "UniRef90_A0A5E4LPI9":
            action = "manual_review_before_final_tree"
        rows_out.append(
            {
                "flag_class": row.get("flag_class", ""),
                "seq_id": row.get("seq_id", ""),
                "nr80_status": row.get("nr80_status", ""),
                "nr80_representative": row.get("nr80_representative", ""),
                "present_in_core": row.get("present_in_core", ""),
                "present_in_alignment": present,
                "alignment_seq_id": seq_id if present else "",
                "gap_fraction": qc.get("gap_fraction", ""),
                "gap_or_ambiguity_fraction": qc.get(
                    "gap_or_ambiguity_fraction", ""
                ),
                "recommended_action": action,
                "notes": notes.strip(),
            }
        )
    columns = [
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
        "notes",
    ]
    write_tsv(
        FLAGGED_CARRY_TSV,
        columns,
        [{key: fmt(value) for key, value in row.items()} for row in rows_out],
    )


def write_high_gap_table(qc_rows: list[dict[str, object]]) -> None:
    rows = [
        row
        for row in qc_rows
        if float(row["gap_or_ambiguity_fraction"]) > 0.50
    ]
    columns = [
        "seq_id",
        "subtype_source",
        "gap_fraction",
        "ambiguity_fraction",
        "gap_or_ambiguity_fraction",
        "core_length_before_alignment",
        "notes",
    ]
    write_tsv(
        HIGH_GAP_TSV,
        columns,
        [{key: fmt(value) for key, value in row.items()} for row in rows],
    )


def artifact_row(
    path: Path,
    role: str,
    required: bool,
    source_inputs: str,
    recreate: str,
    notes: str,
) -> dict[str, object]:
    return {
        "relative_path": rel(path),
        "artifact_role": role,
        "file_size_bytes": path.stat().st_size if path.exists() else "",
        "md5": md5(path) if path.exists() else "",
        "tracked_by_git": git_tracked(path) if path.exists() else "false",
        "ignored_by_git": git_ignored(path) if path.exists() else "false",
        "required_for_next_stage": str(required).lower(),
        "source_inputs": source_inputs,
        "recreate_command_or_source": recreate,
        "notes": notes if path.exists() else f"MISSING: {notes}",
    }


def write_artifact_manifest() -> None:
    source_inputs = (
        "results/03_core_len255/all_core_only_len255.fasta;"
        "results/03_msa_core/core_global.hmm"
    )
    msa_cmds = f"{HMMALIGN_CMD}; {ALIMASK_CMD}; {REFORMAT_CMD}"
    rows = [
        artifact_row(
            CORE_FASTA,
            "strategyA_len255_core_extracted_fasta",
            True,
            "",
            "conda run -n dah7ps_v4 python scripts/extract_core_domains.py ...",
            "Round 4 ignored FASTA input for Round 5 MSA.",
        ),
        artifact_row(
            CORE_HMM,
            "staged_core_global_hmm",
            True,
            "/home/luogu/dah7ps_evo/results/03_msa_core/core_global.hmm",
            "cp /home/luogu/dah7ps_evo/results/03_msa_core/core_global.hmm results/03_msa_core/core_global.hmm",
            "Staged HMM input; no legacy MSA outputs copied.",
        ),
        artifact_row(
            UNMASKED_STO,
            "strategyA_len255_core_unmasked_stockholm",
            True,
            source_inputs,
            HMMALIGN_CMD,
            "Unmasked hmmalign Stockholm output.",
        ),
        artifact_row(
            MASKED_STO,
            "strategyA_len255_core_rf_masked_stockholm",
            True,
            "results/04_msa_len255/core_len255.sto",
            ALIMASK_CMD,
            "RF-masked Stockholm output from esl-alimask.",
        ),
        artifact_row(
            MASKED_AFA,
            "strategyA_len255_core_rf_masked_afa",
            True,
            "results/04_msa_len255/core_len255.masked.sto",
            REFORMAT_CMD,
            "Round 6 tree-inference candidate alignment.",
        ),
        artifact_row(
            ALIGNMENT_QC_TSV,
            "strategyA_len255_alignment_per_sequence_qc",
            True,
            f"{source_inputs};results/04_msa_len255/core_len255.masked.afa",
            "python scripts/round5_msa_len255_handoff.py",
            "Per-sequence gap and ambiguity QC.",
        ),
        artifact_row(
            SUMMARY_TSV,
            "strategyA_len255_alignment_summary",
            True,
            f"{source_inputs};results/04_msa_len255/core_len255.masked.afa",
            "python scripts/round5_msa_len255_handoff.py",
            "Normalized alignment summary.",
        ),
        artifact_row(
            SUBTYPE_QC_TSV,
            "strategyA_len255_subtype_alignment_qc",
            True,
            "results/02_qc_len255/nr80_all_len255_subtype_map.tsv;"
            "results/04_msa_len255/core_len255.masked.afa",
            "python scripts/round5_msa_len255_handoff.py",
            "Subtype gap and ambiguity burden summary.",
        ),
        artifact_row(
            RESCUED479_QC_TSV,
            "strategyA_len255_rescued479_alignment_qc",
            True,
            "results/02_qc/qc_classification_Ib.tsv;"
            "results/02_qc_len255/qc_classification_Ib.tsv;"
            "results/04_msa_len255/core_len255.masked.afa",
            "python scripts/round5_msa_len255_handoff.py",
            "QC for reconstructed 479 len255 Ib rescued sequences.",
        ),
        artifact_row(
            TARGET_CARRY_TSV,
            "strategyA_len255_target_alignment_carrythrough",
            True,
            "results/03_core_len255/core_len255_target_carrythrough.tsv;"
            "results/04_msa_len255/core_len255.alignment_qc.tsv",
            "python scripts/round5_msa_len255_handoff.py",
            "Target carry-through after MSA.",
        ),
        artifact_row(
            FLAGGED_CARRY_TSV,
            "strategyA_len255_flagged_alignment_carrythrough",
            True,
            "results/03_core_len255/core_len255_flagged_carrythrough.tsv;"
            "results/04_msa_len255/core_len255.alignment_qc.tsv",
            "python scripts/round5_msa_len255_handoff.py",
            "Flagged len255 rescue carry-through after MSA.",
        ),
        artifact_row(
            HIGH_GAP_TSV,
            "strategyA_len255_high_gap_or_ambiguity_table",
            True,
            "results/04_msa_len255/core_len255.alignment_qc.tsv",
            "python scripts/round5_msa_len255_handoff.py",
            "Sequences with gap_or_ambiguity_fraction > 0.50.",
        ),
        artifact_row(
            TREE_HANDOFF_MD,
            "strategyA_round5_tree_handoff_plan",
            True,
            "Round 5 alignment and QC outputs",
            "python scripts/round5_msa_len255_handoff.py",
            "Round 6 tree handoff; no IQ-TREE run.",
        ),
    ]
    columns = [
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
    ]
    write_tsv(ARTIFACT_MANIFEST, columns, rows)


def write_tree_handoff(
    qc_rows: list[dict[str, object]],
    target_rows: list[dict[str, str]],
    flagged_rows: list[dict[str, str]],
    masked_columns: int,
) -> None:
    high_gap_ids = [
        str(row["seq_id"])
        for row in qc_rows
        if float(row["gap_or_ambiguity_fraction"]) > 0.50
    ]
    target_lookup = {row["abbrev"]: row for row in target_rows}
    flagged_present = [
        row["seq_id"] for row in flagged_rows if row.get("present_in_alignment") == "true"
    ]
    text = f"""# Strategy A Round 5 Tree Handoff

## 1. Alignment For Round 6

Round 6 should use:

```text
results/04_msa_len255/core_len255.masked.afa
```

## 2. Masking Status

This is the RF-masked alignment produced by `hmmalign --outformat Stockholm`, `esl-alimask --rf-is-mask`, and `esl-reformat afa`.

## 3. Alignment Dimensions

- Sequences: {len(qc_rows)}
- Columns: {masked_columns}

## 4. Target Presence

- Q8U0A9 representative `UniRef90_UPI0002AF51CE`: {target_lookup.get('PfuDAH7PS', {}).get('present_in_alignment', '')}
- Q9YEJ7 direct tip `UniRef90_Q9YEJ7`: {target_lookup.get('ApeDAH7PS', {}).get('present_in_alignment', '')}
- Pni representative `UniRef90_A0A379DXQ3`: {target_lookup.get('PniDAH7PS', {}).get('present_in_alignment', '')}; remains `needs_curation`.

## 5. Flagged Rescue Review

Flagged rescue hits present in the alignment and requiring review before tree inference:

{chr(10).join(f'- `{seq_id}`' for seq_id in flagged_present)}

## 6. Sequences Exceeding 50 Percent Gap Or Ambiguity

- Count: {len(high_gap_ids)}
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
"""
    TREE_HANDOFF_MD.write_text(text)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for path in [
        CORE_FASTA,
        CORE_HMM,
        UNMASKED_STO,
        MASKED_STO,
        MASKED_AFA,
        SUBTYPE_MAP,
        TARGET_REPRESENTATION,
        CORE_TARGET_CARRY,
        CORE_FLAGGED_CARRY,
        CORE_SUBTYPE_CARRY,
        OLD_IB_QC,
        LEN255_IB_QC,
        NR80_IB_CLSTR,
        FLAGGED_ADJUDICATION,
    ]:
        require_file(path)

    core_records = read_fasta(CORE_FASTA)
    alignment_records = read_fasta_records(MASKED_AFA)
    alignment = dict(alignment_records)
    unmasked = read_stockholm(UNMASKED_STO)
    subtypes = load_subtypes()
    nr80_ids = set(subtypes)

    qc_rows = compute_alignment_qc(alignment_records, core_records, subtypes)
    write_alignment_qc(qc_rows)
    alignment_qc = {
        str(row["seq_id"]): {key: fmt(value) for key, value in row.items()}
        for row in qc_rows
    }
    write_summary(core_records, alignment_records, unmasked, qc_rows)
    write_subtype_qc(qc_rows)
    write_rescued479_qc(alignment_qc, core_records, nr80_ids)
    write_target_carrythrough(alignment_qc)
    write_flagged_carrythrough(alignment_qc)
    write_high_gap_table(qc_rows)
    target_rows = read_tsv(TARGET_CARRY_TSV)
    flagged_rows = read_tsv(FLAGGED_CARRY_TSV)
    masked_columns = len(next(iter(alignment.values()))) if alignment else 0
    write_tree_handoff(qc_rows, target_rows, flagged_rows, masked_columns)
    write_artifact_manifest()

    print(f"wrote {rel(SUMMARY_TSV)}")
    print(f"wrote {rel(ALIGNMENT_QC_TSV)}")
    print(f"wrote {rel(SUBTYPE_QC_TSV)}")
    print(f"wrote {rel(RESCUED479_QC_TSV)}")
    print(f"wrote {rel(TARGET_CARRY_TSV)}")
    print(f"wrote {rel(FLAGGED_CARRY_TSV)}")
    print(f"wrote {rel(HIGH_GAP_TSV)}")
    print(f"wrote {rel(ARTIFACT_MANIFEST)}")
    print(f"wrote {rel(TREE_HANDOFF_MD)}")


if __name__ == "__main__":
    main()
