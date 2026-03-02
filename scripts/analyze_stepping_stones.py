#!/usr/bin/env python3
"""Analyze stepping stones cluster composition and generate panel subset.

Outputs:
  1. Per-subtype distribution of 258 stepping stone reps
  2. Cross-subtype mixing analysis
  3. Secondary clustering at --min-seq-id 0.25 for panel candidates (target 20-50)
"""

import argparse
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


def load_ids_from_fasta(fasta_path):
    """Extract sequence IDs from FASTA file."""
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                ids.add(line[1:].strip().split()[0])
    return ids


def get_subtype(seqid, subtype_map):
    return subtype_map.get(seqid, 'Unknown')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', default='/home/tynan/0218')
    parser.add_argument('--seeds60_ia', default='results/02_qc/seeds60_Ia.fasta')
    parser.add_argument('--seeds60_ib', default='results/02_qc/seeds60_Ib.fasta')
    parser.add_argument('--seeds60_ii', default='results/02_qc/seeds60_II.fasta')
    parser.add_argument('--cluster_tsv', default='results/02_qc/stepping_stones_cluster.tsv')
    parser.add_argument('--rep_fasta', default='results/02_qc/stepping_stones_rep_seq.fasta')
    parser.add_argument('--outdir', default='results/02_qc')
    args = parser.parse_args()

    wd = Path(args.workdir)

    # Build subtype map from seeds60 files
    subtype_map = {}
    for st, f in [('Ia', args.seeds60_ia), ('Ib', args.seeds60_ib), ('II', args.seeds60_ii)]:
        for sid in load_ids_from_fasta(wd / f):
            subtype_map[sid] = st

    print(f"Loaded subtype map: Ia={sum(1 for v in subtype_map.values() if v=='Ia')}, "
          f"Ib={sum(1 for v in subtype_map.values() if v=='Ib')}, "
          f"II={sum(1 for v in subtype_map.values() if v=='II')}")

    # Parse cluster.tsv
    clusters = defaultdict(list)
    with open(wd / args.cluster_tsv) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                rep, member = parts
                clusters[rep].append(member)

    print(f"\nTotal clusters (stepping stone reps): {len(clusters)}")

    # Rep subtype distribution
    rep_subtypes = Counter(get_subtype(rep, subtype_map) for rep in clusters)
    print("\n=== Stepping stones reps by subtype ===")
    for st in ['Ia', 'Ib', 'II', 'Unknown']:
        if rep_subtypes[st]:
            print(f"  {st}: {rep_subtypes[st]}")
    print(f"  Total: {sum(rep_subtypes.values())}")

    # Cross-subtype mixing
    mixed_count = 0
    mixed_details = []
    for rep, members in clusters.items():
        subtypes = set(get_subtype(m, subtype_map) for m in members)
        if len(subtypes) > 1:
            mixed_count += 1
            rep_st = get_subtype(rep, subtype_map)
            member_sts = Counter(get_subtype(m, subtype_map) for m in members)
            mixed_details.append((rep, rep_st, dict(member_sts), len(members)))

    print(f"\n=== Cross-subtype cluster analysis ===")
    print(f"  Pure (single-subtype) clusters: {len(clusters) - mixed_count}")
    print(f"  Mixed (cross-subtype) clusters: {mixed_count}")

    if mixed_details:
        print("\n  Mixed clusters detail (up to 15):")
        for rep, rep_st, sts, size in mixed_details[:15]:
            st_str = ", ".join(f"{k}:{v}" for k, v in sorted(sts.items()))
            print(f"    rep={rep} (subtype={rep_st}), size={size}, composition: {st_str}")
        if len(mixed_details) > 15:
            print(f"    ... ({len(mixed_details) - 15} more)")

    # Cluster size distribution
    sizes = [len(m) for m in clusters.values()]
    print(f"\n=== Cluster size distribution ===")
    print(f"  Min: {min(sizes)}, Max: {max(sizes)}, "
          f"Mean: {sum(sizes)/len(sizes):.1f}, Median: {sorted(sizes)[len(sizes)//2]}")
    size_hist = Counter()
    for s in sizes:
        if s == 1:
            size_hist['singleton'] += 1
        elif s <= 5:
            size_hist['2-5'] += 1
        elif s <= 10:
            size_hist['6-10'] += 1
        elif s <= 20:
            size_hist['11-20'] += 1
        else:
            size_hist['>20'] += 1
    for k in ['singleton', '2-5', '6-10', '11-20', '>20']:
        print(f"  {k}: {size_hist.get(k, 0)}")


if __name__ == '__main__':
    main()
