#!/usr/bin/env python3
"""Build Strategy A all-subtype nr80 len255 handoff artifacts."""

from __future__ import annotations

import csv
import hashlib
import statistics
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/02_qc_len255"

FASTA_INPUTS = [
    ("Ia", ROOT / "results/02_qc/nr80_Ia.fasta"),
    ("Ib_len255", OUTDIR / "nr80_Ib.fasta"),
    ("II", ROOT / "results/02_qc/nr80_II.fasta"),
]

TARGET_REPRESENTATION = OUTDIR / "target_representation.tsv"
ACCEPTABILITY = OUTDIR / "representative_acceptability_round3A.tsv"
FLAGGED_ADJUDICATION = OUTDIR / "len255_flagged_rescue_adjudication.tsv"
GENERATED_MANIFEST = OUTDIR / "generated_artifacts_manifest.tsv"
INPUT_STAGING_MANIFEST = OUTDIR / "input_staging_manifest.tsv"

ALL_FASTA = OUTDIR / "nr80_all_len255.fasta"
SUMMARY_TSV = OUTDIR / "nr80_all_len255_summary.tsv"
SUBTYPE_MAP_TSV = OUTDIR / "nr80_all_len255_subtype_map.tsv"
TARGET_PRESENCE_TSV = OUTDIR / "nr80_all_len255_target_presence.tsv"
FLAGGED_PRESENCE_TSV = OUTDIR / "nr80_all_len255_flagged_presence.tsv"
HANDOFF_MANIFEST_TSV = OUTDIR / "nr80_all_len255_handoff_manifest.tsv"
CORE_HANDOFF_MD = OUTDIR / "round3B_core_extraction_handoff.md"


@dataclass(frozen=True)
class FastaRecord:
    seq_id: str
    header: str
    sequence: str
    subtype_source: str
    source_fasta: Path


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def require_file(path: Path) -> None:
    if not path.is_file():
        raise SystemExit(f"required input missing: {rel(path)}")


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
            writer.writerow({column: row.get(column, "") for column in columns})


def parse_fasta(path: Path, subtype_source: str) -> list[FastaRecord]:
    require_file(path)
    records: list[FastaRecord] = []
    current_header = ""
    current_id = ""
    seq_parts: list[str] = []

    def flush() -> None:
        if not current_id:
            return
        records.append(
            FastaRecord(
                seq_id=current_id,
                header=current_header,
                sequence="".join(seq_parts),
                subtype_source=subtype_source,
                source_fasta=path,
            )
        )

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                current_header = line[1:].strip()
                current_id = current_header.split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    flush()
    return records


def copy_fasta_inputs() -> None:
    if ALL_FASTA.parent != OUTDIR:
        raise SystemExit(f"refusing to write outside {rel(OUTDIR)}: {ALL_FASTA}")
    with ALL_FASTA.open("wb") as out_handle:
        for _, path in FASTA_INPUTS:
            require_file(path)
            data = path.read_bytes()
            out_handle.write(data)
            if data and not data.endswith(b"\n"):
                out_handle.write(b"\n")


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_git(*args: str) -> bool:
    return subprocess.run(
        ["git", *args],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    ).returncode == 0


def git_tracked(path: Path) -> str:
    return str(run_git("ls-files", "--error-unmatch", rel(path))).lower()


def git_ignored(path: Path) -> str:
    return str(run_git("check-ignore", "-q", rel(path))).lower()


def format_value(value: object) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def build_record_indexes(
    records_by_subtype: dict[str, list[FastaRecord]],
) -> tuple[dict[str, FastaRecord], dict[str, list[FastaRecord]]]:
    by_id: dict[str, FastaRecord] = {}
    occurrences: dict[str, list[FastaRecord]] = defaultdict(list)
    for records in records_by_subtype.values():
        for record in records:
            occurrences[record.seq_id].append(record)
            by_id.setdefault(record.seq_id, record)
    return by_id, occurrences


def write_summary(
    records_by_subtype: dict[str, list[FastaRecord]],
    occurrences: dict[str, list[FastaRecord]],
) -> list[str]:
    counts = {subtype: len(records) for subtype, records in records_by_subtype.items()}
    all_records = [record for records in records_by_subtype.values() for record in records]
    expected_sum = counts["Ia"] + counts["Ib_len255"] + counts["II"]
    duplicate_ids = sorted(seq_id for seq_id, hits in occurrences.items() if len(hits) > 1)
    lengths = [len(record.sequence) for record in all_records]
    empty_count = sum(1 for length in lengths if length == 0)
    median = statistics.median(lengths) if lengths else ""

    rows = [
        {
            "metric": "nr80_Ia_count",
            "value": counts["Ia"],
            "notes": "records in results/02_qc/nr80_Ia.fasta",
        },
        {
            "metric": "nr80_Ib_len255_count",
            "value": counts["Ib_len255"],
            "notes": "records in results/02_qc_len255/nr80_Ib.fasta",
        },
        {
            "metric": "nr80_II_count",
            "value": counts["II"],
            "notes": "records in results/02_qc/nr80_II.fasta",
        },
        {
            "metric": "nr80_all_len255_count",
            "value": len(all_records),
            "notes": "records in results/02_qc_len255/nr80_all_len255.fasta",
        },
        {
            "metric": "expected_sum",
            "value": expected_sum,
            "notes": "Ia + Ib_len255 + II",
        },
        {
            "metric": "count_matches_expected_sum",
            "value": len(all_records) == expected_sum,
            "notes": "true if concatenated count equals component sum",
        },
        {
            "metric": "duplicate_sequence_ids_count",
            "value": len(duplicate_ids),
            "notes": "unique sequence IDs observed more than once",
        },
        {
            "metric": "duplicate_sequence_ids",
            "value": ",".join(duplicate_ids),
            "notes": "comma-separated; empty means none",
        },
        {
            "metric": "empty_sequences_count",
            "value": empty_count,
            "notes": "records with zero residues after parsing",
        },
        {
            "metric": "min_length",
            "value": min(lengths) if lengths else "",
            "notes": "minimum parsed sequence length",
        },
        {
            "metric": "median_length",
            "value": median,
            "notes": "median parsed sequence length",
        },
        {
            "metric": "max_length",
            "value": max(lengths) if lengths else "",
            "notes": "maximum parsed sequence length",
        },
    ]
    columns = ["metric", "value", "notes"]
    write_tsv(
        SUMMARY_TSV,
        columns,
        [{key: format_value(value) for key, value in row.items()} for row in rows],
    )
    return duplicate_ids


def write_subtype_map(records_by_subtype: dict[str, list[FastaRecord]]) -> None:
    columns = ["seq_id", "subtype_source", "source_fasta", "length"]
    rows: list[dict[str, object]] = []
    for records in records_by_subtype.values():
        for record in records:
            rows.append(
                {
                    "seq_id": record.seq_id,
                    "subtype_source": record.subtype_source,
                    "source_fasta": rel(record.source_fasta),
                    "length": len(record.sequence),
                }
            )
    write_tsv(SUBTYPE_MAP_TSV, columns, rows)


def acceptability_by_target() -> dict[str, dict[str, str]]:
    rows = read_tsv(ACCEPTABILITY)
    by_target: dict[str, dict[str, str]] = {}
    for row in rows:
        by_target[row["target"]] = row
        if row.get("target_id"):
            by_target[row["target_id"]] = row
    return by_target


def representative_source(
    row: dict[str, str],
    by_id: dict[str, FastaRecord],
    representative: str,
) -> str:
    direct_id = row.get("uniref_id") or row.get("primary_accession", "")
    if representative in by_id:
        return by_id[representative].subtype_source
    if direct_id in by_id:
        return by_id[direct_id].subtype_source
    return ""


def resolve_representative_id(
    representative_id: str,
    by_id: dict[str, FastaRecord],
) -> tuple[str, bool, str]:
    if not representative_id:
        return "", False, ""
    if representative_id in by_id:
        return representative_id, True, ""
    prefix_matches = sorted(
        seq_id for seq_id in by_id if seq_id.startswith(representative_id)
    )
    if len(prefix_matches) == 1:
        resolved_id = prefix_matches[0]
        return (
            resolved_id,
            True,
            "target_representation representative_id "
            f"{representative_id} resolved to full FASTA ID {resolved_id} "
            "by unique prefix match.",
        )
    if len(prefix_matches) > 1:
        return (
            representative_id,
            False,
            "target_representation representative_id "
            f"{representative_id} has multiple FASTA prefix matches: "
            + ",".join(prefix_matches),
        )
    return representative_id, False, ""


def target_notes(
    row: dict[str, str],
    direct_present: bool,
    representative_present: bool,
    acceptability: dict[str, str] | None,
    resolution_note: str,
) -> str:
    notes = row.get("notes", "")
    status = row.get("final_core_tree_expected_status", "")
    primary = row.get("primary_accession", "")
    uniref = row.get("uniref_id", "")

    additions: list[str] = []
    if primary == "Q8U0A9":
        additions.append(
            "Direct tip is not expected; surrogate UniRef90_UPI0002AF51CE is expected."
        )
    if primary == "Q9YEJ7":
        additions.append("Direct nr80 tip UniRef90_Q9YEJ7 is expected.")
    if primary == "V8CS59" or uniref == "UniRef90_F9DH16":
        additions.append(
            "Current Pni candidate remains represented_by_nr80_surrogate and needs curation."
        )
    if primary == "A0A0F2JEB6":
        additions.append("Legacy accession remains unresolved; not merged with V8CS59.")
    if status in {"direct_nr80_tip", "represented_by_nr80_surrogate"}:
        additions.append(
            f"presence_check direct={str(direct_present).lower()}; "
            f"representative={str(representative_present).lower()}"
        )
    if acceptability and acceptability.get("acceptability_call"):
        additions.append(
            "round3A_acceptability="
            f"{acceptability.get('acceptability_call')}; "
            f"round3A_curation={acceptability.get('curation_call')}"
        )
    if resolution_note:
        additions.append(resolution_note)
    return " ".join(part for part in [notes, *additions] if part)


def write_target_presence(by_id: dict[str, FastaRecord]) -> None:
    acceptability = acceptability_by_target()
    columns = [
        "abbrev",
        "primary_accession",
        "uniref_id",
        "expected_status_from_target_representation",
        "direct_tip_present_in_nr80_all_len255",
        "representative_present_in_nr80_all_len255",
        "representative_id",
        "subtype_source",
        "label_policy",
        "curation_call",
        "notes",
    ]
    rows: list[dict[str, object]] = []
    absent_expected: list[str] = []
    for row in read_tsv(TARGET_REPRESENTATION):
        direct_id = row.get("uniref_id") or row.get("primary_accession", "")
        representative_id = row.get("nr80_representative", "")
        resolved_representative_id, representative_present, resolution_note = (
            resolve_representative_id(representative_id, by_id)
        )
        status = row.get("final_core_tree_expected_status", "")
        direct_present = bool(direct_id and direct_id in by_id)
        acc_row = (
            acceptability.get(row.get("abbrev", ""))
            or acceptability.get(row.get("uniref_id", ""))
            or acceptability.get(row.get("primary_accession", ""))
        )
        curation_call = row.get("curation_call", "")
        if acc_row and acc_row.get("curation_call"):
            curation_call = acc_row["curation_call"]
        if (
            status in {"direct_nr80_tip", "represented_by_nr80_surrogate"}
            and representative_id
            and not representative_present
        ):
            absent_expected.append(f"{row.get('abbrev')}:{representative_id}")
        rows.append(
            {
                "abbrev": row.get("abbrev", ""),
                "primary_accession": row.get("primary_accession", ""),
                "uniref_id": row.get("uniref_id", ""),
                "expected_status_from_target_representation": status,
                "direct_tip_present_in_nr80_all_len255": direct_present,
                "representative_present_in_nr80_all_len255": representative_present,
                "representative_id": resolved_representative_id,
                "subtype_source": representative_source(
                    row, by_id, resolved_representative_id
                ),
                "label_policy": row.get("label_policy", ""),
                "curation_call": curation_call,
                "notes": target_notes(
                    row,
                    direct_present,
                    representative_present,
                    acc_row,
                    resolution_note,
                ),
            }
        )
    write_tsv(
        TARGET_PRESENCE_TSV,
        columns,
        [{key: format_value(value) for key, value in row.items()} for row in rows],
    )
    if absent_expected:
        raise SystemExit(
            "expected nr80 representatives absent from nr80_all_len255: "
            + ", ".join(absent_expected)
        )


def write_flagged_presence(by_id: dict[str, FastaRecord]) -> None:
    columns = [
        "flag_class",
        "seq_id",
        "nr80_status",
        "nr80_representative",
        "present_in_nr80_all_len255",
        "subtype_source",
        "recommended_action",
        "notes",
    ]
    rows: list[dict[str, object]] = []
    for row in read_tsv(FLAGGED_ADJUDICATION):
        seq_id = row.get("seq_id", "")
        representative = row.get("nr80_representative", "")
        lookup_id = representative or seq_id
        present = bool(lookup_id and lookup_id in by_id)
        source = by_id[lookup_id].subtype_source if lookup_id in by_id else ""
        notes = row.get("notes", "")
        if row.get("nr80_status") == "direct_nr80_representative" and present:
            notes = (
                f"{notes} Direct nr80 representative enters nr80_all_len255; "
                "carry forward for manual review before final tree."
            ).strip()
        rows.append(
            {
                "flag_class": row.get("flag_class", ""),
                "seq_id": seq_id,
                "nr80_status": row.get("nr80_status", ""),
                "nr80_representative": representative,
                "present_in_nr80_all_len255": present,
                "subtype_source": source,
                "recommended_action": row.get("recommended_action", ""),
                "notes": notes,
            }
        )
    write_tsv(
        FLAGGED_PRESENCE_TSV,
        columns,
        [{key: format_value(value) for key, value in row.items()} for row in rows],
    )


def source_lookup() -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for path in [INPUT_STAGING_MANIFEST, GENERATED_MANIFEST]:
        if not path.is_file():
            continue
        for row in read_tsv(path):
            lookup[row.get("relative_path", "")] = row
    return lookup


def source_command_or_note(
    path: Path,
    lookup: dict[str, dict[str, str]],
) -> tuple[str, str]:
    row = lookup.get(rel(path), {})
    if "recreate_command_or_source" in row and row.get("recreate_command_or_source"):
        return row["recreate_command_or_source"], row.get("notes", "")
    if "source_worktree" in row:
        source = (
            f"copied from {row.get('source_worktree')} "
            f"{row.get('source_branch')}@{row.get('source_head')}"
        )
        note = (
            "Staged reference artifact copied from old worktree; "
            "not rerun on strategyA-ib-len255-qc."
        )
        return source, note
    return "", ""


def manifest_row(
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
        "tracked_by_git": git_tracked(path),
        "ignored_by_git": git_ignored(path),
        "required_for_next_stage": str(required).lower(),
        "source_inputs": source_inputs,
        "recreate_command_or_source": recreate,
        "notes": notes if path.exists() else f"MISSING: {notes}",
    }


def write_handoff_manifest() -> None:
    lookup = source_lookup()
    ia_cmd, ia_note = source_command_or_note(ROOT / "results/02_qc/nr80_Ia.fasta", lookup)
    ib_cmd, ib_note = source_command_or_note(OUTDIR / "nr80_Ib.fasta", lookup)
    ii_cmd, ii_note = source_command_or_note(ROOT / "results/02_qc/nr80_II.fasta", lookup)
    inputs = (
        "results/02_qc/nr80_Ia.fasta;"
        "results/02_qc_len255/nr80_Ib.fasta;"
        "results/02_qc/nr80_II.fasta"
    )
    script_cmd = "python scripts/build_nr80_all_len255.py"
    cat_cmd = (
        "cat results/02_qc/nr80_Ia.fasta "
        "results/02_qc_len255/nr80_Ib.fasta "
        "results/02_qc/nr80_II.fasta "
        "> results/02_qc_len255/nr80_all_len255.fasta"
    )
    rows = [
        manifest_row(
            ROOT / "results/02_qc/nr80_Ia.fasta",
            "staged_ia_nr80_representatives",
            True,
            "",
            ia_cmd,
            ia_note,
        ),
        manifest_row(
            OUTDIR / "nr80_Ib.fasta",
            "formal_len255_ib_nr80_representatives",
            True,
            "results/02_qc_len255/Ib_len255_pass_plus_long.fasta",
            ib_cmd,
            ib_note,
        ),
        manifest_row(
            ROOT / "results/02_qc/nr80_II.fasta",
            "staged_ii_nr80_representatives",
            True,
            "",
            ii_cmd,
            ii_note,
        ),
        manifest_row(
            ALL_FASTA,
            "strategyA_all_subtype_nr80_len255_fasta",
            True,
            inputs,
            cat_cmd,
            "Versioned all-subtype nr80 FASTA; legacy "
            "results/02_qc/nr80_all.fasta not overwritten.",
        ),
        manifest_row(
            SUMMARY_TSV,
            "strategyA_all_subtype_nr80_len255_summary",
            True,
            inputs,
            script_cmd,
            "Counts, duplicate-ID audit, and length summary for nr80_all_len255.",
        ),
        manifest_row(
            SUBTYPE_MAP_TSV,
            "strategyA_all_subtype_nr80_len255_subtype_map",
            True,
            inputs,
            script_cmd,
            "Per-tip subtype/source map for the versioned all-subtype nr80 FASTA.",
        ),
        manifest_row(
            TARGET_PRESENCE_TSV,
            "strategyA_all_subtype_nr80_len255_target_presence",
            True,
            f"{inputs};results/02_qc_len255/target_representation.tsv",
            script_cmd,
            "Target direct-tip and representative presence audit after nr80_all_len255 merge.",
        ),
        manifest_row(
            FLAGGED_PRESENCE_TSV,
            "strategyA_all_subtype_nr80_len255_flagged_presence",
            True,
            f"{inputs};results/02_qc_len255/len255_flagged_rescue_adjudication.tsv",
            script_cmd,
            "Flagged len255 rescue carry-forward presence audit.",
        ),
        manifest_row(
            CORE_HANDOFF_MD,
            "strategyA_round3b_core_extraction_handoff_plan",
            True,
            "repo scripts/docs inspection",
            "manual markdown prepared in Round 3B",
            "Exact next-stage core extraction handoff plan; no extraction run.",
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
    write_tsv(HANDOFF_MANIFEST_TSV, columns, rows)


def main() -> None:
    for _, path in FASTA_INPUTS:
        require_file(path)
    for path in [
        TARGET_REPRESENTATION,
        GENERATED_MANIFEST,
        ACCEPTABILITY,
        FLAGGED_ADJUDICATION,
    ]:
        require_file(path)

    records_by_subtype = {
        subtype: parse_fasta(path, subtype) for subtype, path in FASTA_INPUTS
    }
    by_id, occurrences = build_record_indexes(records_by_subtype)

    copy_fasta_inputs()
    duplicate_ids = write_summary(records_by_subtype, occurrences)
    write_subtype_map(records_by_subtype)
    if duplicate_ids:
        raise SystemExit(
            "duplicate IDs detected across nr80 inputs; stopping before handoff: "
            + ", ".join(duplicate_ids)
        )

    write_target_presence(by_id)
    write_flagged_presence(by_id)
    write_handoff_manifest()

    print(f"wrote {rel(ALL_FASTA)}")
    print(f"wrote {rel(SUMMARY_TSV)}")
    print(f"wrote {rel(SUBTYPE_MAP_TSV)}")
    print(f"wrote {rel(TARGET_PRESENCE_TSV)}")
    print(f"wrote {rel(FLAGGED_PRESENCE_TSV)}")
    print(f"wrote {rel(HANDOFF_MANIFEST_TSV)}")


if __name__ == "__main__":
    main()
