#!/usr/bin/env python3
"""Build a provenance-focused artifact manifest from a spec TSV.

The input spec must contain at least:
  artifact_id, scenario_id, artifact_role, file_path, source_script, command,
  input_paths, generated_at, formal_status

The script computes:
  - output_md5
  - input_md5
  - git_commit

Usage:
  python scripts/build_artifact_manifest.py \
    --spec results/04_phylogeny_asr/artifact_manifest.spec.tsv \
    --output results/04_phylogeny_asr/artifact_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import subprocess
import sys
from pathlib import Path


REQUIRED_COLUMNS = [
    "artifact_id",
    "scenario_id",
    "artifact_role",
    "file_path",
    "source_script",
    "command",
    "input_paths",
    "generated_at",
    "formal_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build artifact_manifest.tsv with MD5, commit, and input MD5."
    )
    parser.add_argument("--spec", required=True, help="Input spec TSV")
    parser.add_argument("--output", required=True, help="Output manifest TSV")
    return parser.parse_args()


def fail(message: str) -> None:
    print(f"[build_artifact_manifest] ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def md5_of_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except Exception:
        return "UNKNOWN"


def normalize_paths(raw: str) -> list[str]:
    if not raw.strip():
        return []
    return [part.strip() for part in raw.split(",") if part.strip()]


def main() -> None:
    args = parse_args()
    spec_path = Path(args.spec)
    output_path = Path(args.output)

    if not spec_path.is_file():
        fail(f"spec not found: {spec_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    commit = git_commit()

    rows: list[dict[str, str]] = []
    with spec_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail("spec file has no header")
        for column in REQUIRED_COLUMNS:
            if column not in reader.fieldnames:
                fail(f"missing required column: {column}")

        for row in reader:
            file_path = Path(row["file_path"])
            if not file_path.is_file():
                fail(f"artifact file not found: {file_path}")

            input_paths = normalize_paths(row["input_paths"])
            missing_inputs = [path for path in input_paths if not Path(path).is_file()]
            if missing_inputs:
                fail(
                    "missing input files for "
                    f"{row['artifact_id']}: {', '.join(missing_inputs)}"
                )

            row["output_md5"] = md5_of_file(file_path)
            row["input_md5"] = ",".join(md5_of_file(Path(path)) for path in input_paths)
            row["git_commit"] = commit
            rows.append(row)

    fieldnames = [
        "artifact_id",
        "scenario_id",
        "artifact_role",
        "file_path",
        "output_md5",
        "source_script",
        "git_commit",
        "command",
        "input_paths",
        "input_md5",
        "generated_at",
        "formal_status",
        "notes",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})

    print(f"[build_artifact_manifest] Wrote {len(rows)} rows to {output_path}")


if __name__ == "__main__":
    main()
