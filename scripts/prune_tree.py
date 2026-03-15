#!/usr/bin/env python3
"""Prune tips from a rooted Newick tree by name prefix.

Designed for Phase 4.3: removing KDOPS outgroup tips from the rooted
ML tree to create an ingroup-only tree for downstream ASR.

Usage:
    python scripts/prune_tree.py \
        --input results/04_phylogeny_asr/CoreTree_rooted_MFP.treefile \
        --remove_prefix KDOPS_ \
        --output results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
        --assert_rooted
"""

import argparse
import sys
import re


def parse_args():
    p = argparse.ArgumentParser(
        description="Prune tips matching a prefix from a rooted Newick tree."
    )
    p.add_argument("--input", required=True, help="Input rooted Newick tree file")
    p.add_argument("--output", required=True, help="Output pruned Newick tree file")
    p.add_argument(
        "--remove_prefix",
        required=True,
        help="Prefix of tip names to remove (e.g. KDOPS_)",
    )
    p.add_argument(
        "--assert_rooted",
        action="store_true",
        help="Assert the output tree is still rooted (bifurcating root)",
    )
    return p.parse_args()


def read_newick(path: str) -> str:
    """Read a Newick string from a file."""
    with open(path) as f:
        return f.read().strip()


def extract_tip_names(nwk: str) -> list[str]:
    """Extract all tip names from a Newick string using regex."""
    # Tip names are tokens that are not preceded by ')' and are followed by
    # ':' (branch length) or ',' or ')'. We match any word-like token before ':'
    # But Newick tip names can contain various chars. We use a robust pattern.
    # In IQ-TREE output, tip names are alphanumeric with underscores.
    tips = re.findall(r"[\(,]\s*([A-Za-z0-9_.\-|]+)\s*:", nwk)
    return tips


def prune_newick_string(nwk: str, remove_names: set[str]) -> str:
    """Prune tips from a Newick string by direct string manipulation.

    This is much faster than loading a full tree object for large trees.
    Strategy: iteratively remove leaf tokens and clean up the resulting
    tree structure (empty clades, dangling commas, single-child internal nodes).
    """
    from ete3 import Tree

    print("  Loading tree with ete3 (this may take a minute for large trees)...")
    sys.stdout.flush()
    # format=1 includes internal node names and branch lengths
    t = Tree(nwk, format=1)

    total_tips = len(t.get_leaves())
    print(f"  Tree loaded: {total_tips} tips")

    # Find tips to remove
    tips_to_remove = [leaf for leaf in t.get_leaves() if leaf.name in remove_names]
    found_names = {leaf.name for leaf in tips_to_remove}
    missing = remove_names - found_names
    if missing:
        print(f"  WARNING: {len(missing)} names not found in tree: {missing}")

    print(f"  Tips to prune: {len(tips_to_remove)}")

    # --- Monophyly diagnostic ---
    if len(found_names) >= 2:
        try:
            mono, clade_type, violating = t.check_monophyly(
                values=found_names, target_attr="name"
            )
            if mono:
                print(f"  KDOPS monophyly: YES (monophyletic)")
            else:
                print(f"  KDOPS monophyly: NO ({clade_type})")
                if violating:
                    viol_names = [n.name for n in list(violating)[:5]]
                    print(f"    Violating taxa (first 5): {viol_names}")
        except Exception as e:
            print(f"  Monophyly check failed: {e}")

    # --- Check root position relative to KDOPS ---
    root = t.get_tree_root()
    root_children = root.get_children()
    print(f"  Root children: {len(root_children)}")
    for i, child in enumerate(root_children):
        child_leaves = child.get_leaf_names()
        kdops_in_child = [n for n in child_leaves if n in found_names]
        print(
            f"    Child {i}: {len(child_leaves)} tips "
            f"({len(kdops_in_child)} KDOPS)"
        )

    # --- Prune by keeping only non-matching tips ---
    keep_names = [
        leaf.name for leaf in t.get_leaves() if leaf.name not in remove_names
    ]
    print(f"  Keeping {len(keep_names)} tips, pruning {len(tips_to_remove)}...")
    sys.stdout.flush()
    t.prune(keep_names, preserve_branch_length=True)

    remaining = len(t.get_leaves())
    print(f"  Pruned tree: {remaining} tips")

    return t, remaining


def main():
    args = parse_args()

    print(f"[prune_tree.py] Input: {args.input}")
    print(f"[prune_tree.py] Remove prefix: {args.remove_prefix}")

    nwk = read_newick(args.input)

    # Quick scan for tip names to remove
    all_tips = extract_tip_names(nwk)
    remove_names = {name for name in all_tips if name.startswith(args.remove_prefix)}

    print(f"  Total tips (regex scan): {len(all_tips)}")
    print(f"  Tips matching '{args.remove_prefix}*': {len(remove_names)}")
    for name in sorted(remove_names):
        print(f"    - {name}")

    if not remove_names:
        print("  No tips to remove. Copying input to output.")
        with open(args.output, "w") as f:
            f.write(nwk + "\n")
        return

    # Prune
    pruned_tree, n_remaining = prune_newick_string(nwk, remove_names)

    # Assert rooted
    root_children = pruned_tree.get_tree_root().get_children()
    is_rooted = len(root_children) == 2
    print(f"  Output tree rooted (bifurcating root): {is_rooted}")

    if args.assert_rooted and not is_rooted:
        print(
            f"  ERROR: --assert_rooted failed. "
            f"Root has {len(root_children)} children (expected 2).",
            file=sys.stderr,
        )
        sys.exit(1)

    # Write output
    # format=1 preserves internal node names + branch lengths + support values
    pruned_tree.write(outfile=args.output, format=1)

    print(f"[prune_tree.py] Output: {args.output}")
    print(f"[prune_tree.py] Done. {len(remove_names)} tips pruned, {n_remaining} remain.")


if __name__ == "__main__":
    main()
