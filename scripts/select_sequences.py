#!/usr/bin/env python3
"""
Filter sequences by subtype and required module presence.

Phase 3.9 of the DAH7PS V5.0 SOP.

Given a module presence/absence matrix and an optional FASTA (or ID list)
defining the subtype subset, output the IDs that satisfy all required module
conditions.

Usage:
    python scripts/select_sequences.py \
        --presence_table results/03_msa_modules/module_presence_absence_strict.tsv \
        --subtype_fasta  results/02_qc/nr80_Ib.fasta \
        --require_module ACT_domain \
        --output         results/03_msa_full/Ib_ACT.ids

    # Multiple required modules (all must be 1):
    python scripts/select_sequences.py \
        --presence_table results/03_msa_modules/module_presence_absence_strict.tsv \
        --subtype_fasta  results/02_qc/nr80_Ib.fasta \
        --require_core 1 \
        --require_module ACT_domain CM_domain \
        --output results/03_msa_full/Ib_ACT_CM.ids
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Filter sequences by subtype FASTA and required module flags."
    )
    p.add_argument("--presence_table", required=True,
                   help="module_presence_absence_strict/relaxed.tsv")
    p.add_argument("--subtype_fasta",
                   help="FASTA file whose sequence IDs define the allowed subtype subset "
                        "(first whitespace-delimited token of each header). "
                        "If omitted, no subtype filter is applied.")
    p.add_argument("--subtype_ids",
                   help="Plain text file of allowed seq IDs (one per line), "
                        "alternative to --subtype_fasta.")
    p.add_argument("--require_core", type=int, choices=[0, 1],
                   help="If set, require the 'core' column to equal this value "
                        "(column must be named 'core' in the table).")
    p.add_argument("--require_module", nargs="+", default=[],
                   metavar="MODULE",
                   help="Module column name(s) that must equal 1 (e.g. ACT_domain).")
    p.add_argument("--exclude_module", nargs="+", default=[],
                   metavar="MODULE",
                   help="Module column name(s) that must equal 0.")
    p.add_argument("--output", required=True,
                   help="Output file: one seq_id per line.")
    return p.parse_args()


def read_ids_from_fasta(path):
    ids = set()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split()[0].strip())
    return ids


def read_ids_from_file(path):
    ids = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                ids.add(line)
    return ids


def main():
    args = parse_args()

    # Validate inputs
    if not os.path.isfile(args.presence_table):
        print(f"[ERROR] presence_table not found: {args.presence_table}", file=sys.stderr)
        sys.exit(1)

    # Load subtype filter
    subtype_ids = None
    if args.subtype_fasta:
        if not os.path.isfile(args.subtype_fasta):
            print(f"[ERROR] subtype_fasta not found: {args.subtype_fasta}", file=sys.stderr)
            sys.exit(1)
        subtype_ids = read_ids_from_fasta(args.subtype_fasta)
        print(f"[select_sequences] Loaded {len(subtype_ids)} subtype IDs from {args.subtype_fasta}",
              file=sys.stderr)
    elif args.subtype_ids:
        if not os.path.isfile(args.subtype_ids):
            print(f"[ERROR] subtype_ids not found: {args.subtype_ids}", file=sys.stderr)
            sys.exit(1)
        subtype_ids = read_ids_from_file(args.subtype_ids)
        print(f"[select_sequences] Loaded {len(subtype_ids)} subtype IDs from {args.subtype_ids}",
              file=sys.stderr)

    # Parse presence/absence table and apply filters
    selected = []
    header = None
    n_total = 0
    with open(args.presence_table) as f:
        for line in f:
            line = line.rstrip("\n")
            if header is None:
                header = line.split("\t")
                # Validate requested columns
                for mod in args.require_module + args.exclude_module:
                    if mod not in header:
                        print(f"[ERROR] Module '{mod}' not found in table columns: {header}",
                              file=sys.stderr)
                        sys.exit(1)
                if args.require_core is not None and "core" not in header:
                    print(f"[WARNING] --require_core set but 'core' column not in table; "
                          f"ignoring core filter.", file=sys.stderr)
                    args.require_core = None
                continue

            fields = line.split("\t")
            row = dict(zip(header, fields))
            seq_id = row["seq_id"]
            n_total += 1

            # Subtype filter
            if subtype_ids is not None and seq_id not in subtype_ids:
                continue

            # Core filter
            if args.require_core is not None:
                if int(row.get("core", -1)) != args.require_core:
                    continue

            # Required module filter (all must be 1)
            skip = False
            for mod in args.require_module:
                if int(row.get(mod, 0)) != 1:
                    skip = True
                    break
            if skip:
                continue

            # Excluded module filter (all must be 0)
            for mod in args.exclude_module:
                if int(row.get(mod, 1)) != 0:
                    skip = True
                    break
            if skip:
                continue

            selected.append(seq_id)

    print(f"[select_sequences] {n_total} sequences in table → {len(selected)} pass all filters",
          file=sys.stderr)

    if len(selected) == 0:
        print("[ERROR] No sequences passed the filters. Check module names and subtype FASTA.",
              file=sys.stderr)
        sys.exit(1)

    # Write output
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as f:
        for sid in selected:
            f.write(sid + "\n")
    print(f"[select_sequences] Wrote {len(selected)} IDs to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
