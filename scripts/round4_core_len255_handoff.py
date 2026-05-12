#!/usr/bin/env python3
"""Generate Strategy A Round 4 core-extraction handoff artifacts."""

from __future__ import annotations

import csv
import hashlib
import statistics
import subprocess
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/03_core_len255"
SOURCE_WORKTREE = Path("/home/luogu/dah7ps_evo")

NR80_FASTA = ROOT / "results/02_qc_len255/nr80_all_len255.fasta"
SUBTYPE_MAP = ROOT / "results/02_qc_len255/nr80_all_len255_subtype_map.tsv"
TARGET_REPRESENTATION = ROOT / "results/02_qc_len255/target_representation.tsv"
TARGET_PRESENCE = ROOT / "results/02_qc_len255/nr80_all_len255_target_presence.tsv"
FLAGGED_PRESENCE = ROOT / "results/02_qc_len255/nr80_all_len255_flagged_presence.tsv"
FLAGGED_ADJUDICATION = ROOT / "results/02_qc_len255/len255_flagged_rescue_adjudication.tsv"
CORE_HMM = ROOT / "results/03_msa_core/core_global.hmm"
CORE_FASTA = OUTDIR / "all_core_only_len255.fasta"
CORE_TSV = OUTDIR / "core_extraction_len255.tsv"
CORE_DOMTBLOUT = OUTDIR / "core_extraction_len255_domtblout.txt"
CORE_LOG = OUTDIR / "core_extraction_len255.log"

INPUT_STAGING_MANIFEST = OUTDIR / "input_staging_manifest.tsv"
SUMMARY_TSV = OUTDIR / "core_extraction_len255_summary.tsv"
SUBTYPE_CARRY_TSV = OUTDIR / "core_len255_subtype_carrythrough.tsv"
TARGET_CARRY_TSV = OUTDIR / "core_len255_target_carrythrough.tsv"
FLAGGED_CARRY_TSV = OUTDIR / "core_len255_flagged_carrythrough.tsv"
ARTIFACT_MANIFEST = OUTDIR / "core_len255_artifact_manifest.tsv"
MSA_HANDOFF_MD = OUTDIR / "round4_msa_handoff.md"

CORE_COMMAND = (
    "conda run -n dah7ps_v4 python scripts/extract_core_domains.py "
    "--hmm results/03_msa_core/core_global.hmm "
    "--fasta results/02_qc_len255/nr80_all_len255.fasta "
    "--params meta/params.json "
    "--out_fasta results/03_core_len255/all_core_only_len255.fasta "
    "--out_tsv results/03_core_len255/core_extraction_len255.tsv "
    "--ievalue 1e-5 --hmm_span_min 30 --merge_gap 5 --cpu 20"
)


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
            writer.writerow({column: row.get(column, "") for column in columns})


def read_fasta(path: Path) -> dict[str, str]:
    require_file(path)
    records: dict[str, str] = {}
    current_id = ""
    seq_parts: list[str] = []

    def flush() -> None:
        if current_id:
            records[current_id] = "".join(seq_parts)

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


def duplicate_ids(path: Path) -> list[str]:
    counts: Counter[str] = Counter()
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                counts[line[1:].split()[0]] += 1
    return sorted(seq_id for seq_id, count in counts.items() if count > 1)


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_text(worktree: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=worktree,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else ""


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


def write_input_staging_manifest() -> None:
    source_branch = git_text(SOURCE_WORKTREE, "branch", "--show-current")
    source_head = git_text(SOURCE_WORKTREE, "rev-parse", "HEAD")
    current_branch = git_text(ROOT, "branch", "--show-current")
    current_head = git_text(ROOT, "rev-parse", "HEAD")
    rows = [
        {
            "relative_path": "results/03_msa_core/core_global.hmm",
            "source_worktree": SOURCE_WORKTREE.as_posix(),
            "source_branch": source_branch,
            "source_head": source_head,
            "destination_path": CORE_HMM.as_posix(),
            "file_size_bytes": CORE_HMM.stat().st_size,
            "md5": md5(CORE_HMM),
            "purpose": "core_extraction_hmm",
            "status": "copied",
        },
        {
            "relative_path": "results/02_qc_len255/nr80_all_len255.fasta",
            "source_worktree": ROOT.as_posix(),
            "source_branch": current_branch,
            "source_head": current_head,
            "destination_path": NR80_FASTA.as_posix(),
            "file_size_bytes": NR80_FASTA.stat().st_size,
            "md5": md5(NR80_FASTA),
            "purpose": "core_extraction_input_fasta",
            "status": "existing_round3b_artifact",
        },
        {
            "relative_path": "results/02_qc_len255/target_representation.tsv",
            "source_worktree": ROOT.as_posix(),
            "source_branch": current_branch,
            "source_head": current_head,
            "destination_path": TARGET_REPRESENTATION.as_posix(),
            "file_size_bytes": TARGET_REPRESENTATION.stat().st_size,
            "md5": md5(TARGET_REPRESENTATION),
            "purpose": "target_audit_reference",
            "status": "existing_round3a_artifact",
        },
        {
            "relative_path": "results/02_qc_len255/nr80_all_len255_target_presence.tsv",
            "source_worktree": ROOT.as_posix(),
            "source_branch": current_branch,
            "source_head": current_head,
            "destination_path": TARGET_PRESENCE.as_posix(),
            "file_size_bytes": TARGET_PRESENCE.stat().st_size,
            "md5": md5(TARGET_PRESENCE),
            "purpose": "target_audit_reference",
            "status": "existing_round3b_artifact",
        },
        {
            "relative_path": "results/02_qc_len255/nr80_all_len255_subtype_map.tsv",
            "source_worktree": ROOT.as_posix(),
            "source_branch": current_branch,
            "source_head": current_head,
            "destination_path": SUBTYPE_MAP.as_posix(),
            "file_size_bytes": SUBTYPE_MAP.stat().st_size,
            "md5": md5(SUBTYPE_MAP),
            "purpose": "target_audit_reference",
            "status": "existing_round3b_artifact",
        },
        {
            "relative_path": "results/02_qc_len255/nr80_all_len255_flagged_presence.tsv",
            "source_worktree": ROOT.as_posix(),
            "source_branch": current_branch,
            "source_head": current_head,
            "destination_path": FLAGGED_PRESENCE.as_posix(),
            "file_size_bytes": FLAGGED_PRESENCE.stat().st_size,
            "md5": md5(FLAGGED_PRESENCE),
            "purpose": "target_audit_reference",
            "status": "existing_round3b_artifact",
        },
    ]
    columns = [
        "relative_path",
        "source_worktree",
        "source_branch",
        "source_head",
        "destination_path",
        "file_size_bytes",
        "md5",
        "purpose",
        "status",
    ]
    write_tsv(INPUT_STAGING_MANIFEST, columns, rows)


def write_summary(nr80_records: dict[str, str], core_records: dict[str, str]) -> None:
    core_lengths = [len(seq) for seq in core_records.values()]
    dupes = duplicate_ids(CORE_FASTA)
    input_count = len(nr80_records)
    core_count = len(core_records)
    rows = [
        {
            "metric": "nr80_all_len255_input_count",
            "value": input_count,
            "notes": "records in results/02_qc_len255/nr80_all_len255.fasta",
        },
        {
            "metric": "core_extracted_count",
            "value": core_count,
            "notes": "records in results/03_core_len255/all_core_only_len255.fasta",
        },
        {
            "metric": "core_failed_or_missing_count",
            "value": input_count - core_count,
            "notes": "input records absent from core FASTA after HMM hit/coverage filters",
        },
        {
            "metric": "core_extracted_fraction",
            "value": core_count / input_count if input_count else 0,
            "notes": "core_extracted_count / nr80_all_len255_input_count",
        },
        {
            "metric": "min_core_length",
            "value": min(core_lengths) if core_lengths else "",
            "notes": "minimum extracted core sequence length",
        },
        {
            "metric": "median_core_length",
            "value": statistics.median(core_lengths) if core_lengths else "",
            "notes": "median extracted core sequence length",
        },
        {
            "metric": "max_core_length",
            "value": max(core_lengths) if core_lengths else "",
            "notes": "maximum extracted core sequence length",
        },
        {
            "metric": "empty_core_sequences_count",
            "value": sum(1 for length in core_lengths if length == 0),
            "notes": "core FASTA records with zero residues",
        },
        {
            "metric": "duplicate_core_ids_count",
            "value": len(dupes),
            "notes": "unique sequence IDs observed more than once in core FASTA",
        },
        {
            "metric": "duplicate_core_ids",
            "value": ",".join(dupes),
            "notes": "comma-separated; empty means none",
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


def write_subtype_carrythrough(core_records: dict[str, str]) -> None:
    rows_by_subtype: dict[str, list[str]] = {"Ia": [], "Ib_len255": [], "II": []}
    for row in read_tsv(SUBTYPE_MAP):
        subtype = row["subtype_source"]
        if subtype in rows_by_subtype:
            rows_by_subtype[subtype].append(row["seq_id"])

    rows: list[dict[str, object]] = []
    for subtype in ["Ia", "Ib_len255", "II"]:
        ids = rows_by_subtype[subtype]
        extracted = sum(1 for seq_id in ids if seq_id in core_records)
        total = len(ids)
        fraction = extracted / total if total else 0.0
        rows.append(
            {
                "subtype_source": subtype,
                "nr80_input_count": total,
                "core_extracted_count": extracted,
                "core_failed_or_missing_count": total - extracted,
                "core_extracted_fraction": fraction,
            }
        )
    columns = [
        "subtype_source",
        "nr80_input_count",
        "core_extracted_count",
        "core_failed_or_missing_count",
        "core_extracted_fraction",
    ]
    write_tsv(
        SUBTYPE_CARRY_TSV,
        columns,
        [{key: fmt(value) for key, value in row.items()} for row in rows],
    )


def by_key(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in rows if row.get(key)}


def target_note(
    row: dict[str, str],
    direct_core: bool,
    representative_core: bool,
) -> str:
    primary = row.get("primary_accession", "")
    uniref = row.get("uniref_id", "")
    status = row.get("expected_status_from_target_representation", "")
    notes = [row.get("notes", "")]
    if primary == "Q8U0A9":
        notes.append(
            "Direct tip absent as expected; label representative as "
            "PfuDAH7PS-like / represents Q8U0A9 if present in core."
        )
    if primary == "Q9YEJ7":
        notes.append("Direct core sequence should be labeled ApeDAH7PS if present.")
    if primary == "V8CS59" or uniref == "UniRef90_F9DH16":
        notes.append(
            "Pni candidate surrogate remains needs_curation; do not call accepted."
        )
    if primary == "A0A0F2JEB6":
        notes.append("Legacy accession remains unresolved and separate from V8CS59.")
    if status in {"direct_nr80_tip", "represented_by_nr80_surrogate"}:
        notes.append(
            f"core_presence direct={str(direct_core).lower()}; "
            f"representative={str(representative_core).lower()}"
        )
    return " ".join(note for note in notes if note)


def write_target_carrythrough(core_records: dict[str, str]) -> None:
    target_rows = by_key(read_tsv(TARGET_REPRESENTATION), "abbrev")
    rows: list[dict[str, object]] = []
    missing_expected: list[str] = []
    columns = [
        "abbrev",
        "primary_accession",
        "uniref_id",
        "expected_status_from_target_representation",
        "direct_tip_present_in_nr80_all_len255",
        "representative_id",
        "representative_present_in_nr80_all_len255",
        "direct_tip_present_in_core",
        "representative_present_in_core",
        "core_sequence_id_used_for_label",
        "core_length",
        "label_policy",
        "curation_call",
        "notes",
    ]
    for presence in read_tsv(TARGET_PRESENCE):
        abbrev = presence["abbrev"]
        target = target_rows.get(abbrev, {})
        direct_id = presence.get("uniref_id") or presence.get("primary_accession", "")
        representative_id = presence.get("representative_id", "")
        status = presence.get("expected_status_from_target_representation", "")
        direct_core = bool(direct_id and direct_id in core_records)
        representative_core = bool(representative_id and representative_id in core_records)
        if status == "direct_nr80_tip":
            label_id = direct_id if direct_core else representative_id
        elif status == "represented_by_nr80_surrogate":
            label_id = representative_id
        else:
            label_id = ""
        core_length = len(core_records[label_id]) if label_id in core_records else ""
        if (
            status in {"direct_nr80_tip", "represented_by_nr80_surrogate"}
            and presence.get("representative_present_in_nr80_all_len255") == "true"
            and not representative_core
        ):
            missing_expected.append(f"{abbrev}:{representative_id}")
        curation_call = presence.get("curation_call", "") or target.get("curation_call", "")
        if abbrev == "PniDAH7PS":
            curation_call = "needs_curation"
        rows.append(
            {
                "abbrev": abbrev,
                "primary_accession": presence.get("primary_accession", ""),
                "uniref_id": presence.get("uniref_id", ""),
                "expected_status_from_target_representation": status,
                "direct_tip_present_in_nr80_all_len255": presence.get(
                    "direct_tip_present_in_nr80_all_len255", ""
                ),
                "representative_id": representative_id,
                "representative_present_in_nr80_all_len255": presence.get(
                    "representative_present_in_nr80_all_len255", ""
                ),
                "direct_tip_present_in_core": direct_core,
                "representative_present_in_core": representative_core,
                "core_sequence_id_used_for_label": label_id,
                "core_length": core_length,
                "label_policy": presence.get("label_policy", ""),
                "curation_call": curation_call,
                "notes": target_note(presence, direct_core, representative_core),
            }
        )
    write_tsv(
        TARGET_CARRY_TSV,
        columns,
        [{key: fmt(value) for key, value in row.items()} for row in rows],
    )
    if missing_expected:
        raise SystemExit(
            "expected target representatives missing from core output: "
            + ", ".join(missing_expected)
        )


def write_flagged_carrythrough(core_records: dict[str, str]) -> None:
    adjudication = by_key(read_tsv(FLAGGED_ADJUDICATION), "seq_id")
    rows: list[dict[str, object]] = []
    columns = [
        "flag_class",
        "seq_id",
        "nr80_status",
        "nr80_representative",
        "present_in_nr80_all_len255",
        "present_in_core",
        "core_sequence_id",
        "core_length",
        "recommended_action",
        "notes",
    ]
    for row in read_tsv(FLAGGED_PRESENCE):
        seq_id = row.get("seq_id", "")
        representative = row.get("nr80_representative", "")
        lookup_id = representative or seq_id
        present = lookup_id in core_records
        note = row.get("notes", "")
        if present:
            note = (
                f"{note} Present in core output; carry forward before MSA/tree review."
            ).strip()
        if seq_id == "UniRef90_A0A5E4LPI9":
            action = "manual_review_before_final_tree"
        else:
            action = (
                row.get("recommended_action", "")
                or adjudication.get(seq_id, {}).get("recommended_action", "")
            )
        rows.append(
            {
                "flag_class": row.get("flag_class", ""),
                "seq_id": seq_id,
                "nr80_status": row.get("nr80_status", ""),
                "nr80_representative": representative,
                "present_in_nr80_all_len255": row.get("present_in_nr80_all_len255", ""),
                "present_in_core": present,
                "core_sequence_id": lookup_id if present else "",
                "core_length": len(core_records[lookup_id]) if present else "",
                "recommended_action": action,
                "notes": note,
            }
        )
    write_tsv(
        FLAGGED_CARRY_TSV,
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
        "results/02_qc_len255/nr80_all_len255.fasta;"
        "results/03_msa_core/core_global.hmm"
    )
    rows = [
        artifact_row(
            NR80_FASTA,
            "strategyA_all_subtype_nr80_len255_fasta",
            True,
            "",
            "python scripts/build_nr80_all_len255.py",
            "Round 3B ignored FASTA input for core extraction.",
        ),
        artifact_row(
            CORE_HMM,
            "staged_core_global_hmm",
            True,
            "/home/luogu/dah7ps_evo/results/03_msa_core/core_global.hmm",
            "cp /home/luogu/dah7ps_evo/results/03_msa_core/core_global.hmm results/03_msa_core/core_global.hmm",
            "Staged HMM only; no old MSA/tree outputs copied.",
        ),
        artifact_row(
            CORE_FASTA,
            "strategyA_len255_core_extracted_fasta",
            True,
            source_inputs,
            CORE_COMMAND,
            "Ignored FASTA for Round 5 MSA input.",
        ),
        artifact_row(
            CORE_TSV,
            "strategyA_len255_core_extraction_table",
            True,
            source_inputs,
            CORE_COMMAND,
            "Versioned core extraction coordinate table.",
        ),
        artifact_row(
            CORE_DOMTBLOUT,
            "strategyA_len255_core_extraction_domtblout",
            False,
            source_inputs,
            CORE_COMMAND,
            "HMMER domtblout intermediate from core extraction.",
        ),
        artifact_row(
            CORE_LOG,
            "strategyA_len255_core_extraction_log",
            False,
            source_inputs,
            CORE_COMMAND,
            "Ignored command log from core extraction.",
        ),
        artifact_row(
            SUMMARY_TSV,
            "strategyA_len255_core_extraction_summary",
            True,
            f"{source_inputs};results/03_core_len255/all_core_only_len255.fasta",
            "python scripts/round4_core_len255_handoff.py",
            "Normalized Round 4 extraction summary.",
        ),
        artifact_row(
            SUBTYPE_CARRY_TSV,
            "strategyA_len255_core_subtype_carrythrough",
            True,
            "results/02_qc_len255/nr80_all_len255_subtype_map.tsv;"
            "results/03_core_len255/all_core_only_len255.fasta",
            "python scripts/round4_core_len255_handoff.py",
            "Subtype carry-through after core extraction.",
        ),
        artifact_row(
            TARGET_CARRY_TSV,
            "strategyA_len255_core_target_carrythrough",
            True,
            "results/02_qc_len255/target_representation.tsv;"
            "results/02_qc_len255/nr80_all_len255_target_presence.tsv;"
            "results/03_core_len255/all_core_only_len255.fasta",
            "python scripts/round4_core_len255_handoff.py",
            "Target carry-through after core extraction.",
        ),
        artifact_row(
            FLAGGED_CARRY_TSV,
            "strategyA_len255_core_flagged_carrythrough",
            True,
            "results/02_qc_len255/nr80_all_len255_flagged_presence.tsv;"
            "results/02_qc_len255/len255_flagged_rescue_adjudication.tsv;"
            "results/03_core_len255/all_core_only_len255.fasta",
            "python scripts/round4_core_len255_handoff.py",
            "Flagged len255 rescue carry-through after core extraction.",
        ),
        artifact_row(
            MSA_HANDOFF_MD,
            "strategyA_round4_msa_handoff_plan",
            True,
            "repo scripts/docs inspection;Round 4 core carry-through tables",
            "python scripts/round4_core_len255_handoff.py",
            "Round 5 MSA handoff; no MSA run.",
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


def write_msa_handoff() -> None:
    text = """# Strategy A Round 4 MSA Handoff

## 1. FASTA For Round 5

Round 5 should align:

```text
results/03_core_len255/all_core_only_len255.fasta
```

## 2. Candidate MSA Path

No repo script appears to wrap the first core MSA step. The prior workflow in `log.md` used HMMER/Easel directly:

```bash
conda run -n dah7ps_v4 hmmalign --amino --outformat Stockholm \\
  results/03_msa_core/core_global.hmm \\
  results/03_core_len255/all_core_only_len255.fasta \\
  > results/03_core_len255/core_global_len255.sto

conda run -n dah7ps_v4 esl-alimask --rf-is-mask \\
  results/03_core_len255/core_global_len255.sto \\
  > results/03_core_len255/core_global_matchonly_len255.sto

conda run -n dah7ps_v4 esl-reformat afa \\
  results/03_core_len255/core_global_matchonly_len255.sto \\
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
"""
    MSA_HANDOFF_MD.write_text(text)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for path in [
        NR80_FASTA,
        SUBTYPE_MAP,
        TARGET_REPRESENTATION,
        TARGET_PRESENCE,
        FLAGGED_PRESENCE,
        FLAGGED_ADJUDICATION,
        CORE_HMM,
        CORE_FASTA,
        CORE_TSV,
        CORE_LOG,
    ]:
        require_file(path)

    nr80_records = read_fasta(NR80_FASTA)
    core_records = read_fasta(CORE_FASTA)

    write_input_staging_manifest()
    write_summary(nr80_records, core_records)
    write_subtype_carrythrough(core_records)
    write_target_carrythrough(core_records)
    write_flagged_carrythrough(core_records)
    write_msa_handoff()
    write_artifact_manifest()

    print(f"wrote {rel(INPUT_STAGING_MANIFEST)}")
    print(f"wrote {rel(SUMMARY_TSV)}")
    print(f"wrote {rel(SUBTYPE_CARRY_TSV)}")
    print(f"wrote {rel(TARGET_CARRY_TSV)}")
    print(f"wrote {rel(FLAGGED_CARRY_TSV)}")
    print(f"wrote {rel(ARTIFACT_MANIFEST)}")
    print(f"wrote {rel(MSA_HANDOFF_MD)}")


if __name__ == "__main__":
    main()
