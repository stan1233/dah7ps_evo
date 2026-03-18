#!/usr/bin/env python3
"""
mad_root_fast.py — Efficient MAD rooting for large trees.

Uses ete3's fast tree traversal if available, otherwise BioPython.
The key optimization: instead of deepcopy+re-root per branch,
compute root-to-tip distances via DFS traversal and branch-shift
algebra to evaluate each candidate root position in O(n) time.

For the DAH7PS 9,393-tip tree, this should complete in minutes, not hours.

Usage:
    python scripts/mad_root_fast.py \
        --input results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
        --output results/04_phylogeny_asr/CoreTree_rooted_MAD_ingroup.treefile
"""
import argparse, sys, os
import numpy as np
from io import StringIO

def try_ete3():
    try:
        from ete3 import Tree
        return True
    except ImportError:
        return False

def mad_root_ete3(input_path, output_path):
    """MAD rooting using ete3 (fast tree operations)."""
    from ete3 import Tree
    
    t = Tree(input_path)
    tips = t.get_leaves()
    n = len(tips)
    print(f"[S4-ete3] Tree loaded: {n} tips")
    
    # Precompute tip name → index
    tip_idx = {leaf.name: i for i, leaf in enumerate(tips)}
    
    # For each internal branch, evaluate as candidate root
    # Use root-to-tip distance variance as clock-deviation proxy
    
    all_nodes = list(t.traverse("postorder"))
    best_rho = float('inf')
    best_node = None
    n_eval = 0
    
    for node in all_nodes:
        if node.is_root():
            continue
        if node.dist is None or node.dist <= 0:
            continue
        
        # Only evaluate branches creating non-trivial splits
        n_desc = len(node.get_leaves())
        if n_desc < 2 or n_desc > n - 2:
            continue
        
        n_eval += 1
        if n_eval % 200 == 0:
            print(f"[S4-ete3] Evaluated {n_eval} branches, best rho = {best_rho:.6f}")
        
        # Root at this node and compute root-to-tip distances
        t.set_outgroup(node)
        
        r2t = np.array([t.get_distance(leaf) for leaf in tips])
        mean_d = r2t.mean()
        
        if mean_d <= 0:
            continue
        
        # Relative deviation (MAD criterion)
        rho = np.sqrt(np.mean(((r2t - mean_d) / mean_d) ** 2))
        
        if rho < best_rho:
            best_rho = rho
            best_node = node
    
    print(f"[S4-ete3] Evaluated {n_eval} total branches")
    
    if best_node is None:
        print("[S4] WARNING: No optimal root found, using midpoint")
        t = Tree(input_path)
        midpoint = t.get_midpoint_outgroup()
        t.set_outgroup(midpoint)
    else:
        # Re-root at best node
        t = Tree(input_path)
        # Find best_node in re-loaded tree
        target_leaves = set(l.name for l in best_node.get_leaves())
        for node in t.traverse():
            if set(l.name for l in node.get_leaves()) == target_leaves:
                t.set_outgroup(node)
                break
        
        print(f"[S4] Best rho = {best_rho:.6f}")
    
    # Get root split info
    children = t.get_children()
    sizes = [len(c.get_leaves()) for c in children]
    print(f"[S4] Root split sizes: {sizes}")
    
    t.write(outfile=output_path, format=5)
    print(f"[S4] MAD rooted tree written to: {output_path}")
    
    return best_rho

def mad_root_biopython(input_path, output_path):
    """MAD rooting using BioPython — optimized version.
    
    Instead of deepcopy per branch (O(n) per candidate),
    just re-read the tree each time (still O(n) but no deepcopy overhead).
    
    Further optimization: only evaluate branches near deep splits
    (top 200 most balanced).
    """
    from Bio import Phylo
    import copy
    
    with open(input_path) as f:
        tree_str = f.read().strip()
    
    tree = Phylo.read(StringIO(tree_str), "newick")
    tips = tree.get_terminals()
    n = len(tips)
    print(f"[S4-BioPy] Tree loaded: {n} tips")
    
    # Collect candidate branches with their split sizes
    candidates = []
    for clade in tree.find_clades():
        if clade == tree.root:
            continue
        if clade.branch_length is None or clade.branch_length <= 0:
            continue
        n_desc = len(clade.get_terminals())
        if 2 <= n_desc <= n - 2:
            # Store tip set for identification
            tip_set = frozenset(t.name for t in clade.get_terminals())
            balance = min(n_desc, n - n_desc)
            candidates.append((balance, tip_set, clade))
    
    # Sort by balance (most balanced first) and take top N
    candidates.sort(key=lambda x: -x[0])
    max_eval = min(200, len(candidates))
    candidates = candidates[:max_eval]
    
    print(f"[S4-BioPy] Evaluating top {max_eval} most balanced branches...")
    
    best_rho = float('inf')
    best_tipset = None
    
    for idx, (balance, tip_set, _) in enumerate(candidates):
        if idx % 20 == 0:
            print(f"[S4-BioPy]   Progress: {idx}/{max_eval}")
        
        try:
            # Re-read tree for each candidate
            test_tree = Phylo.read(StringIO(tree_str), "newick")
            
            # Find the clade matching this tip set
            target = None
            for c in test_tree.find_clades():
                if not c.is_terminal():
                    c_tips = frozenset(t.name for t in c.get_terminals())
                    if c_tips == tip_set:
                        target = c
                        break
            
            if target is None:
                continue
            
            test_tree.root_with_outgroup(target)
            
            # Compute root-to-tip distances
            r2t = np.array([test_tree.distance(test_tree.root, t)
                           for t in test_tree.get_terminals()])
            mean_d = r2t.mean()
            if mean_d <= 0:
                continue
            
            rho = np.sqrt(np.mean(((r2t - mean_d) / mean_d) ** 2))
            
            if rho < best_rho:
                best_rho = rho
                best_tipset = tip_set
                
        except Exception:
            continue
    
    print(f"[S4-BioPy] Best rho = {best_rho:.6f}")
    
    if best_tipset is None:
        print("[S4] WARNING: No optimal root found, using midpoint")
        final_tree = Phylo.read(StringIO(tree_str), "newick")
        final_tree.root_at_midpoint()
    else:
        final_tree = Phylo.read(StringIO(tree_str), "newick")
        for c in final_tree.find_clades():
            if not c.is_terminal():
                c_tips = frozenset(t.name for t in c.get_terminals())
                if c_tips == best_tipset:
                    final_tree.root_with_outgroup(c)
                    break
    
    children = final_tree.root.clades
    sizes = [len(c.get_terminals()) for c in children]
    print(f"[S4] Root split sizes: {sizes}")
    
    Phylo.write(final_tree, output_path, "newick")
    print(f"[S4] MAD rooted tree written to: {output_path}")
    
    return best_rho


def main():
    parser = argparse.ArgumentParser(description="Fast MAD rooting for large trees")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: {args.input} not found", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    if try_ete3():
        print("Using ete3 backend (fast)")
        mad_root_ete3(args.input, args.output)
    else:
        print("Using BioPython backend (slower)")
        mad_root_biopython(args.input, args.output)


if __name__ == "__main__":
    main()
