#!/usr/bin/env python3
"""Drop selected FASTA/AFA records by sequence ID.

This is used for provenance-preserving sensitivity inputs where a known
outlier record must be excluded without changing the existing alignment
columns.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Drop selected FASTA records by ID.")
    parser.add_argument("--input", required=True, help="Input FASTA/AFA file")
    parser.add_argument("--output", required=True, help="Output FASTA/AFA file")
    parser.add_argument(
        "--remove-id",
        action="append",
        required=True,
        help="Sequence ID to remove. Can be supplied more than once.",
    )
    parser.add_argument("--expect-input-count", type=int, help="Expected input record count")
    parser.add_argument("--expect-output-count", type=int, help="Expected output record count")
    return parser.parse_args()


def fail(message: str) -> None:
    print(f"[drop_fasta_records] ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def record_id(header: str) -> str:
    return header[1:].strip().split()[0]


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    remove_ids = set(args.remove_id)

    if not input_path.is_file():
        fail(f"input not found: {input_path}")
    if output_path.exists():
        fail(f"output already exists; refusing to overwrite: {output_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    input_count = 0
    output_count = 0
    removed: set[str] = set()
    current_header: str | None = None
    current_seq: list[str] = []

    def flush_record(handle) -> None:
        nonlocal input_count, output_count, current_header, current_seq
        if current_header is None:
            return
        input_count += 1
        seq_id = record_id(current_header)
        if seq_id in remove_ids:
            removed.add(seq_id)
        else:
            output_count += 1
            handle.write(current_header)
            if not current_header.endswith("\n"):
                handle.write("\n")
            for line in current_seq:
                handle.write(line)
                if not line.endswith("\n"):
                    handle.write("\n")

    with input_path.open() as source, output_path.open("w") as target:
        for line in source:
            if line.startswith(">"):
                flush_record(target)
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        flush_record(target)

    missing = remove_ids - removed
    if missing:
        output_path.unlink(missing_ok=True)
        fail(f"requested IDs not found: {','.join(sorted(missing))}")

    if args.expect_input_count is not None and input_count != args.expect_input_count:
        output_path.unlink(missing_ok=True)
        fail(f"input count {input_count} != expected {args.expect_input_count}")

    if args.expect_output_count is not None and output_count != args.expect_output_count:
        output_path.unlink(missing_ok=True)
        fail(f"output count {output_count} != expected {args.expect_output_count}")

    print(f"[drop_fasta_records] Input records: {input_count}")
    print(f"[drop_fasta_records] Removed records: {len(removed)} ({','.join(sorted(removed))})")
    print(f"[drop_fasta_records] Output records: {output_count}")
    print(f"[drop_fasta_records] Output: {output_path}")


if __name__ == "__main__":
    main()
