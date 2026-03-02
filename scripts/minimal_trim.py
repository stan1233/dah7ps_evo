#!/usr/bin/env python3
"""
Minimal MSA trimming: remove only columns that exceed a gap-fraction threshold.

Phase 3.7 of the DAH7PS V4.1 SOP — produces the ASR/DCA version of the
core MSA by removing near-empty columns while preserving maximum site
information for model estimation.

Unlike ClipKIT (which also considers parsimony-informativeness), this script
applies a single, transparent criterion: gap fraction per column.

Outputs:
  - Trimmed alignment (FASTA)
  - Column retention report (TSV): original_col, gap_fraction, kept
  - Summary statistics to stderr
"""

import argparse
import os
import sys


def parse_fasta(path):
    """Parse aligned FASTA, return list of (seqid, sequence) tuples."""
    records = []
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        records.append((current_id, "".join(current_seq)))
    return records


def compute_gap_fractions(records):
    """Compute per-column gap fraction."""
    if not records:
        return []
    n_seqs = len(records)
    n_cols = len(records[0][1])
    gap_counts = [0] * n_cols
    for _, seq in records:
        for i, ch in enumerate(seq):
            if ch == "-" or ch == ".":
                gap_counts[i] += 1
    return [c / n_seqs for c in gap_counts]


def main():
    parser = argparse.ArgumentParser(
        description="Minimal MSA trimming: remove columns exceeding gap-fraction threshold (Phase 3.7)."
    )
    parser.add_argument(
        "--input", required=True, help="Input aligned FASTA"
    )
    parser.add_argument(
        "--output", required=True, help="Output trimmed aligned FASTA"
    )
    parser.add_argument(
        "--gap_col_threshold",
        type=float,
        default=0.95,
        help="Remove columns with gap fraction > this value (default: 0.95)",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Output column report TSV (optional; default: <output>.cols.tsv)",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"[ERROR] File not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    report_path = args.report or args.output.rsplit(".", 1)[0] + ".cols.tsv"

    # Parse
    records = parse_fasta(args.input)
    if not records:
        print("[ERROR] No sequences found in input.", file=sys.stderr)
        sys.exit(1)

    n_seqs = len(records)
    n_cols = len(records[0][1])
    print(f"[minimal_trim] Input: {n_seqs} seqs × {n_cols} cols", file=sys.stderr)
    print(f"[minimal_trim] Gap threshold: > {args.gap_col_threshold}", file=sys.stderr)

    # Validate alignment
    for sid, seq in records:
        if len(seq) != n_cols:
            print(
                f"[ERROR] Sequence {sid} has {len(seq)} cols, expected {n_cols}.",
                file=sys.stderr,
            )
            sys.exit(1)

    # Compute gap fractions
    gap_fracs = compute_gap_fractions(records)

    # Determine kept columns
    kept_cols = [
        i for i, gf in enumerate(gap_fracs) if gf <= args.gap_col_threshold
    ]
    n_kept = len(kept_cols)
    n_removed = n_cols - n_kept
    print(
        f"[minimal_trim] Columns kept: {n_kept} / {n_cols} "
        f"(removed {n_removed}, {n_removed/n_cols*100:.1f}%)",
        file=sys.stderr,
    )

    # Extract kept columns
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    kept_set = set(kept_cols)
    with open(args.output, "w") as f:
        for sid, seq in records:
            trimmed = "".join(seq[i] for i in kept_cols)
            f.write(f">{sid}\n{trimmed}\n")
    print(f"[minimal_trim] Wrote {args.output}", file=sys.stderr)

    # Write column report
    os.makedirs(os.path.dirname(report_path) or ".", exist_ok=True)
    with open(report_path, "w") as f:
        f.write("original_col_1based\tgap_fraction\tkept\n")
        for i, gf in enumerate(gap_fracs):
            f.write(f"{i+1}\t{gf:.6f}\t{'Y' if i in kept_set else 'N'}\n")
    print(f"[minimal_trim] Wrote {report_path}", file=sys.stderr)

    # Summary
    if kept_cols:
        kept_gap_fracs = [gap_fracs[i] for i in kept_cols]
        mean_gap = sum(kept_gap_fracs) / len(kept_gap_fracs)
        max_gap = max(kept_gap_fracs)
        print(
            f"[minimal_trim] Kept columns — mean gap: {mean_gap:.4f}, max gap: {max_gap:.4f}",
            file=sys.stderr,
        )
    else:
        print("[WARNING] All columns removed!", file=sys.stderr)
        sys.exit(1)

    print(f"[minimal_trim] Done. {n_seqs} seqs × {n_kept} cols.", file=sys.stderr)


if __name__ == "__main__":
    main()
