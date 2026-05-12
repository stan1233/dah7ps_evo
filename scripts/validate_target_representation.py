#!/usr/bin/env python3
"""Validate the Strategy A target representation TSV."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = [
    "abbrev",
    "primary_accession",
    "uniref_id",
    "length",
    "subtype_or_expected_group",
    "current_qc_status",
    "len255_qc_status",
    "nr80_status",
    "nr80_representative",
    "nr80_identity_if_absorbed",
    "nr80_local_target_coverage",
    "seeds60_status",
    "final_core_tree_expected_status",
    "label_policy",
    "evidence_source",
    "curation_call",
    "notes",
]

REQUIRED_TARGETS = {
    "PfuDAH7PS",
    "ApeDAH7PS",
    "TmaDAH7PS",
    "GspDAH7PS",
    "LmoDAH7PS",
    "PniDAH7PS",
    "FtuDAH7PS",
    "NmeDAH7PS",
    "SceDAH7PS ARO3",
    "SceDAH7PS ARO4",
    "EcoDAH7PS AroF",
    "EcoDAH7PS AroG",
    "MtuDAH7PS",
    "HpyDAH7PS",
    "CglDAH7PS",
    "PaeDAH7PS PA2843",
    "PaeDAH7PS PA1901",
    "legacy A0A0F2JEB6",
}

ALLOWED_FINAL_STATUS = {
    "direct_nr80_tip",
    "represented_by_nr80_surrogate",
    "absent_from_input",
    "unresolved_accession",
    "needs_curation",
    "hold",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", help="target_representation.tsv path")
    args = parser.parse_args()
    path = Path(args.path)
    errors: list[str] = []
    if not path.is_file():
        print(f"ERROR: not a file: {path}", file=sys.stderr)
        return 1

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != REQUIRED_COLUMNS:
            errors.append("header does not match required columns exactly")
        rows = list(reader)

    seen = {row.get("abbrev", "") for row in rows}
    missing = REQUIRED_TARGETS - seen
    extra = seen - REQUIRED_TARGETS
    if missing:
        errors.append(f"missing targets: {', '.join(sorted(missing))}")
    if extra:
        errors.append(f"unexpected targets: {', '.join(sorted(extra))}")

    for line_no, row in enumerate(rows, start=2):
        status = row.get("final_core_tree_expected_status", "")
        if status not in ALLOWED_FINAL_STATUS:
            errors.append(f"line {line_no}: invalid final status {status!r}")
        if status == "direct_nr80_tip" and row.get("nr80_status") != "direct_nr80_representative":
            errors.append(f"line {line_no}: direct final status without direct nr80 status")
        if status == "represented_by_nr80_surrogate" and not row.get("nr80_representative"):
            errors.append(f"line {line_no}: surrogate status without nr80 representative")
        if not row.get("evidence_source"):
            errors.append(f"line {line_no}: missing evidence_source")

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"OK: validated {len(rows)} target rows in {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
