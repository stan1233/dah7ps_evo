#!/usr/bin/env python3
"""
KDOPS Negative Selection Filter (V4.1)

Compare hmmsearch scores: KDOPS HMM vs DAH7PS Ib HMM on Type Iβ sequences.
Remove sequences where KDOPS score >= DAH7PS Ib score (or only hit KDOPS).

Usage:
  python filter_kdops.py \
    --dah7ps_domtbl hits_Ib_vs_dah7ps.domtbl \
    --kdops_domtbl hits_Ib_vs_kdops.domtbl \
    --input hits_Ib_seqs.fasta \
    --output hits_Ib_clean.fasta \
    --contaminants kdops_contaminants.txt \
    --report kdops_filter_report.tsv
"""

import argparse
import os
import sys


def parse_domtbl(path):
    """Parse domtblout, return dict: seqid -> best full-sequence score (column 7, 0-indexed)."""
    best = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 22:
                continue
            seqid = cols[0]
            score = float(cols[7])  # full sequence score
            if seqid not in best or score > best[seqid]:
                best[seqid] = score
    return best


def write_filtered_fasta(input_path, output_path, remove_ids):
    """Write FASTA excluding remove_ids. Returns (kept, removed) counts."""
    kept = 0
    removed = 0
    skip = False
    with open(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                seqid = line[1:].split()[0]
                if seqid in remove_ids:
                    skip = True
                    removed += 1
                else:
                    skip = False
                    kept += 1
                    fout.write(line)
            else:
                if not skip:
                    fout.write(line)
    return kept, removed


def main():
    parser = argparse.ArgumentParser(
        description="KDOPS negative selection filter for DAH7PS Type Iβ candidates"
    )
    parser.add_argument(
        "--dah7ps_domtbl",
        required=True,
        help="domtblout from hmmsearch of Ib seqs vs DAH7PS Ib HMM",
    )
    parser.add_argument(
        "--kdops_domtbl",
        required=True,
        help="domtblout from hmmsearch of Ib seqs vs KDOPS HMM",
    )
    parser.add_argument(
        "--input", required=True, help="Input FASTA of Ib candidate sequences"
    )
    parser.add_argument(
        "--output", required=True, help="Output FASTA after KDOPS removal"
    )
    parser.add_argument(
        "--contaminants",
        default=None,
        help="Output file listing removed sequence IDs (optional)",
    )
    parser.add_argument(
        "--report", default=None, help="Per-sequence score comparison TSV (optional)"
    )
    parser.add_argument(
        "--score_margin",
        type=float,
        default=0.0,
        help="Remove if kdops_score >= dah7ps_score - margin (default: 0)",
    )
    args = parser.parse_args()

    # Validate inputs
    for path in [args.dah7ps_domtbl, args.kdops_domtbl, args.input]:
        if not os.path.isfile(path):
            print(f"ERROR: Input file not found: {path}", file=sys.stderr)
            sys.exit(1)

    # Ensure output directories exist
    for path in [args.output, args.contaminants, args.report]:
        if path:
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    # Parse scores
    dah7ps_scores = parse_domtbl(args.dah7ps_domtbl)
    kdops_scores = parse_domtbl(args.kdops_domtbl)

    # Classify sequences
    all_seqs = set(dah7ps_scores.keys()) | set(kdops_scores.keys())
    only_kdops = set(kdops_scores.keys()) - set(dah7ps_scores.keys())
    only_dah7ps = set(dah7ps_scores.keys()) - set(kdops_scores.keys())
    both = set(kdops_scores.keys()) & set(dah7ps_scores.keys())

    remove_ids = set()
    # Remove sequences only hit by KDOPS
    remove_ids |= only_kdops
    # Remove sequences where KDOPS wins (with margin)
    for seqid in both:
        if kdops_scores[seqid] >= dah7ps_scores[seqid] - args.score_margin:
            remove_ids.add(seqid)

    # Write per-sequence report if requested
    if args.report:
        with open(args.report, "w") as f:
            f.write("seq_id\tdah7ps_score\tkdops_score\tdelta\tverdict\n")
            for seqid in sorted(all_seqs):
                ds = dah7ps_scores.get(seqid, float("nan"))
                ks = kdops_scores.get(seqid, float("nan"))
                delta = ds - ks if seqid in both else float("nan")
                verdict = "REMOVE" if seqid in remove_ids else "KEEP"
                f.write(f"{seqid}\t{ds:.1f}\t{ks:.1f}\t{delta:.1f}\t{verdict}\n")

    # Write contaminant IDs
    if args.contaminants:
        with open(args.contaminants, "w") as f:
            for seqid in sorted(remove_ids):
                f.write(seqid + "\n")

    # Filter FASTA
    kept, removed = write_filtered_fasta(args.input, args.output, remove_ids)

    # Summary
    kdops_wins = {
        s for s in both if kdops_scores[s] >= dah7ps_scores[s] - args.score_margin
    }
    print(f"=== KDOPS Negative Selection Report (V4.1) ===")
    print(f"DAH7PS domtbl:  {args.dah7ps_domtbl}")
    print(f"KDOPS domtbl:   {args.kdops_domtbl}")
    print(f"Score margin:   {args.score_margin}")
    print(f"")
    print(f"Total unique sequences with hits:  {len(all_seqs)}")
    print(f"  Only DAH7PS (no KDOPS hit):      {len(only_dah7ps)}")
    print(f"  Only KDOPS (no DAH7PS hit):      {len(only_kdops)}")
    print(f"  Hit both HMMs:                   {len(both)}")
    print(f"    DAH7PS wins:                   {len(both) - len(kdops_wins)}")
    print(f"    KDOPS wins (>= margin):        {len(kdops_wins)}")
    print(f"")
    print(f"Total removed:                     {len(remove_ids)}")
    print(f"  (only_kdops={len(only_kdops)} + kdops_wins={len(kdops_wins)})")
    print(f"FASTA kept:                        {kept}")
    print(f"FASTA removed:                     {removed}")

    if args.contaminants:
        print(f"\nContaminant IDs → {args.contaminants}")
    if args.report:
        print(f"Score comparison → {args.report}")
    print(f"Filtered FASTA  → {args.output}")


if __name__ == "__main__":
    main()
