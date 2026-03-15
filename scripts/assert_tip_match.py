#!/usr/bin/env python3
"""Assert that tree tip names match MSA sequence IDs exactly.

V5.0 hard constraint: IQ-TREE -te requires tree and alignment tip sets
to be identical. This script verifies that constraint.

Usage:
    python scripts/assert_tip_match.py \
        --tree results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
        --msa results/03_msa_core/core_asr.afa \
        --assert_identical
"""

import argparse
import sys
import re
from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser(
        description="Assert tree tip set matches MSA sequence ID set."
    )
    p.add_argument("--tree", required=True, help="Newick tree file")
    p.add_argument("--msa", required=True, help="FASTA/AFA alignment file")
    p.add_argument(
        "--assert_identical",
        action="store_true",
        help="Exit with error if tip sets differ",
    )
    return p.parse_args()


def extract_tips_from_newick(path: str) -> set[str]:
    """Extract tip names from a Newick file using regex (fast, no tree loading)."""
    with open(path) as f:
        nwk = f.read().strip()
    tips = re.findall(r"[\(,]\s*([A-Za-z0-9_.\-|]+)\s*:", nwk)
    return set(tips)


def extract_ids_from_msa(path: str) -> set[str]:
    """Extract sequence IDs from a FASTA/AFA file."""
    ids = set()
    for rec in SeqIO.parse(path, "fasta"):
        ids.add(rec.id)
    return ids


def main():
    args = parse_args()

    print(f"[assert_tip_match.py] Tree: {args.tree}")
    print(f"[assert_tip_match.py] MSA:  {args.msa}")

    tree_tips = extract_tips_from_newick(args.tree)
    msa_ids = extract_ids_from_msa(args.msa)

    print(f"  Tree tips: {len(tree_tips)}")
    print(f"  MSA seqs:  {len(msa_ids)}")

    only_tree = tree_tips - msa_ids
    only_msa = msa_ids - tree_tips

    if only_tree:
        print(f"\n  In TREE but not in MSA ({len(only_tree)}):")
        for name in sorted(only_tree)[:20]:
            print(f"    - {name}")
        if len(only_tree) > 20:
            print(f"    ... and {len(only_tree) - 20} more")

    if only_msa:
        print(f"\n  In MSA but not in TREE ({len(only_msa)}):")
        for name in sorted(only_msa)[:20]:
            print(f"    - {name}")
        if len(only_msa) > 20:
            print(f"    ... and {len(only_msa) - 20} more")

    if not only_tree and not only_msa:
        print("\n  ✅ PASS: Tree tips and MSA sequence IDs are IDENTICAL.")
        print(f"     Matched: {len(tree_tips)} / {len(msa_ids)}")
    else:
        print(f"\n  ❌ FAIL: Tip sets differ.")
        print(f"     Only in tree: {len(only_tree)}")
        print(f"     Only in MSA:  {len(only_msa)}")

        if args.assert_identical:
            print(
                "\n  ERROR: --assert_identical flag set. Aborting.",
                file=sys.stderr,
            )
            sys.exit(1)


if __name__ == "__main__":
    main()
