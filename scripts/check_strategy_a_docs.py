#!/usr/bin/env python3
"""Lightweight Strategy A documentation checks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


REQUIRED_SNIPPETS = [
    "QC3",
    "HOLD",
    "CD-HIT thresholds unchanged",
    "target_representation.tsv",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="Markdown documents to check")
    args = parser.parse_args()

    errors: list[str] = []
    for raw_path in args.paths:
        path = Path(raw_path)
        if not path.is_file():
            errors.append(f"missing document: {path}")
            continue
        text = path.read_text()
        if path.name == "log.md":
            errors.append("log.md must not be part of the Strategy A docs check")
        if not text.strip():
            errors.append(f"empty document: {path}")

    combined = "\n".join(Path(path).read_text() for path in args.paths if Path(path).is_file())
    for snippet in REQUIRED_SNIPPETS:
        if snippet not in combined:
            errors.append(f"required Strategy A snippet missing: {snippet}")

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"OK: checked {len(args.paths)} Strategy A docs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
