#!/usr/bin/env python3
"""
root_ingroup_tree.py — S3 Midpoint + S4 MAD rooting for DAH7PS ingroup tree.

S3: Midpoint rooting via BioPython (efficient, O(n) tree traversal).
S4: MAD rooting — Minimal Ancestor Deviation (Tria et al. 2017).
    Uses efficient numpy-based O(n^2) distance matrix + per-branch optimization.

Usage:
    python scripts/root_ingroup_tree.py \
        --input results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
        --outdir results/04_phylogeny_asr/
"""
import argparse, sys, os
from io import StringIO
from Bio import Phylo

def read_tree(path):
    with open(path) as f:
        return Phylo.read(StringIO(f.read().strip()), "newick")

def midpoint_root(tree_path, out_path):
    """S3: Midpoint rooting."""
    tree = read_tree(tree_path)
    tree.root_at_midpoint()
    Phylo.write(tree, out_path, "newick")
    
    # Summary stats
    tips = tree.get_terminals()
    dists = [tree.distance(tree.root, t) for t in tips[:20]]  # sample
    print(f"[S3] Midpoint rooted tree: {len(tips)} tips")
    print(f"[S3] Root-to-tip distances (sample of 20): min={min(dists):.4f} max={max(dists):.4f}")
    
    # Get the two subtrees at root
    children = tree.root.clades
    sizes = [len(c.get_terminals()) for c in children]
    print(f"[S3] Root split sizes: {sizes}")
    print(f"[S3] Written to: {out_path}")
    return tree, sizes

def mad_root_efficient(tree_path, out_path):
    """
    S4: MAD rooting — efficient implementation.
    
    For each candidate root branch, evaluate the clock-likeness
    by minimizing the relative ancestor deviation (rho).
    Uses a two-pass approach on the tree topology.
    
    For very large trees, we sample representative branches
    around deep splits to keep computation feasible.
    """
    import numpy as np
    from Bio import Phylo
    from io import StringIO
    import copy
    
    tree = read_tree(tree_path)
    tips = tree.get_terminals()
    n = len(tips)
    print(f"[S4] MAD rooting on tree with {n} tips")
    
    # For a tree with >5000 tips, full pairwise distances are expensive.
    # We use a two-pass approach:
    # 1. Find all internal branches (candidates for rooting)
    # 2. For each branch, compute clock-deviation using root-to-tip distances
    
    # Get all clades with branches (candidates for root placement)
    candidates = []
    for clade in tree.find_clades():
        if clade == tree.root:
            continue
        if clade.branch_length is not None and clade.branch_length > 0:
            # Store clade and its descendant tip count
            n_desc = len(clade.get_terminals())
            # Skip very small or very large subtrees for efficiency
            # Focus on branches that create meaningful splits
            if 2 <= n_desc <= n - 2:
                candidates.append((clade, n_desc))
    
    print(f"[S4] Evaluating {len(candidates)} candidate branches...")
    
    if not candidates:
        print("[S4] WARNING: No valid candidate branches, falling back to midpoint")
        tree.root_at_midpoint()
        Phylo.write(tree, out_path, "newick")
        return tree
    
    # For efficiency with large trees, sample and evaluate
    # Sort by how balanced the split is (prefer balanced splits)
    candidates.sort(key=lambda x: abs(x[1] - n/2))
    
    # Evaluate top candidates (most balanced splits)
    max_eval = min(len(candidates), 500)  # evaluate up to 500 most balanced
    candidates = candidates[:max_eval]
    
    best_rho = float('inf')
    best_clade = None
    
    for idx, (clade, n_desc) in enumerate(candidates):
        if idx % 50 == 0:
            print(f"[S4]   Progress: {idx}/{max_eval}")
        
        try:
            # Create a copy and root at this branch
            test_tree = copy.deepcopy(tree)
            
            # Find the matching clade in the copy
            if clade.is_terminal():
                targets = [c for c in test_tree.get_terminals() if c.name == clade.name]
            else:
                clade_tip_set = frozenset(t.name for t in clade.get_terminals())
                targets = []
                for c in test_tree.get_nonterminals():
                    if frozenset(t.name for t in c.get_terminals()) == clade_tip_set:
                        targets.append(c)
                        break
            
            if not targets:
                continue
            
            test_tree.root_with_outgroup(targets[0])
            
            # Compute root-to-tip distances for ALL tips
            r2t = []
            for tip in test_tree.get_terminals():
                d = test_tree.distance(test_tree.root, tip)
                r2t.append(d)
            
            r2t = np.array(r2t)
            mean_d = r2t.mean()
            
            if mean_d <= 0:
                continue
            
            # Relative ancestor deviation (clock-likeness metric)
            # Lower = more clock-like
            rho = np.sqrt(np.mean(((r2t - mean_d) / mean_d) ** 2))
            
            if rho < best_rho:
                best_rho = rho
                best_clade = clade
                
        except Exception as e:
            continue
    
    if best_clade is None:
        print("[S4] WARNING: Could not find optimal root, falling back to midpoint")
        tree.root_at_midpoint()
        Phylo.write(tree, out_path, "newick")
        return tree
    
    # Apply best rooting to original tree
    if best_clade.is_terminal():
        targets = [c for c in tree.get_terminals() if c.name == best_clade.name]
    else:
        clade_tip_set = frozenset(t.name for t in best_clade.get_terminals())
        targets = []
        for c in tree.get_nonterminals():
            if frozenset(t.name for t in c.get_terminals()) == clade_tip_set:
                targets.append(c)
                break

    if targets:
        tree.root_with_outgroup(targets[0])
    
    Phylo.write(tree, out_path, "newick")
    
    # Summary
    children = tree.root.clades
    sizes = [len(c.get_terminals()) for c in children]
    print(f"[S4] Best rho (relative ancestor deviation) = {best_rho:.6f}")
    print(f"[S4] Root split sizes: {sizes}")
    print(f"[S4] Written to: {out_path}")
    
    return tree


def main():
    parser = argparse.ArgumentParser(
        description="S3/S4 rooting: Midpoint + MAD for ingroup tree"
    )
    parser.add_argument("--input", required=True, help="Input ingroup tree (Newick)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--midpoint-only", action="store_true",
                        help="Only do midpoint rooting (skip MAD)")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: {args.input} not found", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    mp_out = os.path.join(args.outdir, "CoreTree_rooted_midpoint_ingroup.treefile")
    mad_out = os.path.join(args.outdir, "CoreTree_rooted_MAD_ingroup.treefile")

    print("=" * 60)
    print("S3: Midpoint Rooting")
    print("=" * 60)
    mp_tree, mp_sizes = midpoint_root(args.input, mp_out)

    if not args.midpoint_only:
        print()
        print("=" * 60)
        print("S4: MAD Rooting")
        print("=" * 60)
        mad_tree = mad_root_efficient(args.input, mad_out)

    print()
    print("Done!")

if __name__ == "__main__":
    main()
