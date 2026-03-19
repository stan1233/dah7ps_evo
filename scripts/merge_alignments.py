#!/usr/bin/env python3
"""
Merge outgroup sequences into an existing trimmed core MSA.

The outgroup sequences are aligned to the full (untrimmed) core HMM profile
and therefore have the same number of columns as core_global_matchonly.afa.
This script identifies which columns of the full alignment were retained in
the trimmed alignment (core_tree.afa or core_asr.afa) by comparing a
reference sequence present in both, then subsets the outgroup to the same
columns.

Phase 4.1 of the DAH7PS V5.0 SOP.

Usage:
    python scripts/merge_alignments.py \
        --core       results/03_msa_core/core_tree.afa \
        --core_full  results/03_msa_core/core_global_matchonly.afa \
        --outgroup   results/04_phylogeny_asr/kdops_core_aligned.afa \
        --output     results/04_phylogeny_asr/core_with_outgroup.afa
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Merge outgroup into a trimmed core MSA by column mapping."
    )
    p.add_argument("--core", required=True,
                   help="Trimmed core MSA (e.g. core_tree.afa, 436 cols)")
    p.add_argument("--core_full", required=True,
                   help="Full (untrimmed) core MSA (core_global_matchonly.afa, 521 cols)")
    p.add_argument("--outgroup", required=True,
                   help="Outgroup MSA aligned to same HMM (521 cols)")
    p.add_argument("--output", required=True,
                   help="Output merged MSA (core + outgroup, trimmed cols)")
    return p.parse_args()


def parse_msa(path):
    seqs = {}
    order = []
    cur_id = None
    cur_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_seq)
                    order.append(cur_id)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_id is not None:
        seqs[cur_id] = "".join(cur_seq)
        order.append(cur_id)
    return seqs, order


def find_kept_columns(full_seq, trimmed_seq):
    """Find which columns of full_seq map to trimmed_seq.

    Returns list of 0-based column indices in the full alignment
    that correspond to the trimmed alignment columns.
    """
    kept = []
    t_pos = 0
    t_len = len(trimmed_seq)
    for f_pos in range(len(full_seq)):
        if t_pos >= t_len:
            break
        if full_seq[f_pos].upper() == trimmed_seq[t_pos].upper():
            kept.append(f_pos)
            t_pos += 1
    if t_pos != t_len:
        return None  # matching failed
    return kept


def write_fasta(records, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        for sid, seq in records:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def main():
    args = parse_args()

    for path in [args.core, args.core_full, args.outgroup]:
        if not os.path.isfile(path):
            print(f"[ERROR] File not found: {path}", file=sys.stderr)
            sys.exit(1)

    core_seqs, core_order = parse_msa(args.core)
    full_seqs, _ = parse_msa(args.core_full)
    og_seqs, og_order = parse_msa(args.outgroup)

    core_width = len(next(iter(core_seqs.values())))
    full_width = len(next(iter(full_seqs.values())))
    og_width = len(next(iter(og_seqs.values())))

    print(f"[merge_alignments] core: {len(core_seqs)} seqs × {core_width} cols", file=sys.stderr)
    print(f"[merge_alignments] full: {len(full_seqs)} seqs × {full_width} cols", file=sys.stderr)
    print(f"[merge_alignments] outgroup: {len(og_seqs)} seqs × {og_width} cols", file=sys.stderr)

    if og_width != full_width:
        print(f"[ERROR] Outgroup ({og_width} cols) and full ({full_width} cols) "
              f"must have the same number of columns.", file=sys.stderr)
        sys.exit(1)

    # Find column mapping: use first sequence in core as reference
    kept_cols = None
    ref_id = None
    for sid in core_order:
        if sid in full_seqs:
            ref_id = sid
            kept_cols = find_kept_columns(full_seqs[sid], core_seqs[sid])
            if kept_cols is not None and len(kept_cols) == core_width:
                break
            kept_cols = None

    if kept_cols is None:
        print(f"[ERROR] Could not determine column mapping between full and trimmed MSA. "
              f"No reference sequence matched.", file=sys.stderr)
        sys.exit(1)

    print(f"[merge_alignments] Column mapping: {full_width} → {len(kept_cols)} cols "
          f"(reference: {ref_id})", file=sys.stderr)

    # Validate mapping with 2 more reference sequences
    n_validated = 0
    for sid in core_order:
        if sid == ref_id or sid not in full_seqs:
            continue
        trimmed = "".join(full_seqs[sid][c] for c in kept_cols)
        if trimmed.upper() == core_seqs[sid].upper():
            n_validated += 1
        else:
            print(f"[WARNING] Column mapping validation failed for {sid}", file=sys.stderr)
        if n_validated >= 2:
            break

    if n_validated < 2:
        print(f"[WARNING] Could only validate column mapping with {n_validated} "
              f"additional sequences.", file=sys.stderr)

    print(f"[merge_alignments] Column mapping validated with {n_validated + 1} sequences",
          file=sys.stderr)

    # Apply column mapping to outgroup sequences
    og_trimmed = {}
    for sid in og_order:
        og_trimmed[sid] = "".join(og_seqs[sid][c] for c in kept_cols)

    # Check for ID conflicts
    overlap = set(core_order) & set(og_order)
    if overlap:
        # If outgroup IDs are from --mapali skeleton, remove them
        print(f"[merge_alignments] {len(overlap)} IDs overlap between core and outgroup "
              f"(skeleton sequences from --mapali). Removing duplicates from outgroup.",
              file=sys.stderr)
        og_order = [sid for sid in og_order if sid not in overlap]

    # Write merged output: core first, then outgroup
    records = [(sid, core_seqs[sid]) for sid in core_order]
    records += [(sid, og_trimmed[sid]) for sid in og_order]

    write_fasta(records, args.output)
    total = len(records)
    print(f"[merge_alignments] Output: {total} seqs × {core_width} cols → {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
