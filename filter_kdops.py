#!/usr/bin/env python3
"""
Compare hmmsearch scores: KDOPS HMM vs DAH7PS Ib HMM on Type Iβ sequences.
Remove sequences where KDOPS score >= DAH7PS Ib score (or only hit KDOPS).
"""
import sys

def parse_domtbl(path):
    """Parse domtblout, return dict: seqid -> best full-sequence score (column 7, 0-indexed)."""
    best = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) < 22:
                continue
            seqid = cols[0]
            score = float(cols[7])  # full sequence score
            if seqid not in best or score > best[seqid]:
                best[seqid] = score
    return best

# Parse both result tables
kdops_scores = parse_domtbl('domhits_Ib_vs_kdops.tbl')
dah7ps_scores = parse_domtbl('domhits_Ib_vs_dah7ps.tbl')

# Find sequences to remove: KDOPS score >= DAH7PS score, or only in KDOPS
remove_ids = set()
for seqid, kscore in kdops_scores.items():
    dscore = dah7ps_scores.get(seqid, -999999)
    if kscore >= dscore:
        remove_ids.add(seqid)

# Report
all_seqs = set(dah7ps_scores.keys()) | set(kdops_scores.keys())
only_kdops = set(kdops_scores.keys()) - set(dah7ps_scores.keys())
both = set(kdops_scores.keys()) & set(dah7ps_scores.keys())
kdops_wins = {s for s in both if kdops_scores[s] >= dah7ps_scores[s]}

print(f"=== KDOPS Negative Selection Report ===")
print(f"Total unique sequences with hits:  {len(all_seqs)}")
print(f"Sequences with KDOPS hits:         {len(kdops_scores)}")
print(f"Sequences with DAH7PS Ib hits:     {len(dah7ps_scores)}")
print(f"Only hit KDOPS (no DAH7PS hit):    {len(only_kdops)}")
print(f"Hit both, KDOPS score >= DAH7PS:   {len(kdops_wins)}")
print(f"Total to remove:                   {len(remove_ids)}")
print(f"Remaining after filtering:         {len(all_seqs) - len(remove_ids)}")

# Write IDs to remove
with open('kdops_contaminants.txt', 'w') as f:
    for seqid in sorted(remove_ids):
        f.write(seqid + '\n')

print(f"\nContaminant IDs written to: kdops_contaminants.txt")

# Also write the list of IDs to keep (from raw_full_Ib.fasta)
# We need all IDs from the original FASTA, then subtract contaminants
print(f"\nNow use seqkit to filter:")
print(f"  seqkit grep -v -f kdops_contaminants.txt raw_full_Ib.fasta > raw_full_Ib_clean.fasta")
