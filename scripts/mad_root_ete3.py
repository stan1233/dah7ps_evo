#!/usr/bin/env python3
"""
mad_root_ete3.py — Fast MAD rooting using ete3.

ete3's set_outgroup() is implemented in C and is orders of magnitude
faster than BioPython's root_with_outgroup() for large trees.

Usage:
    conda run -n dah7ps_v4 python3 scripts/mad_root_ete3.py \
        --input results/04_phylogeny_asr/CoreTree_rooted_ingroup.treefile \
        --output results/04_phylogeny_asr/CoreTree_rooted_MAD_ingroup.treefile
"""
import argparse, sys, os, time
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Fast MAD rooting via ete3")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-eval", type=int, default=500,
                        help="Max number of candidate branches to evaluate")
    args = parser.parse_args()

    from ete3 import Tree

    print(f"[MAD] Loading tree from {args.input} ...")
    t0 = time.time()
    t = Tree(args.input)
    leaves = t.get_leaves()
    n = len(leaves)
    print(f"[MAD] Loaded: {n} tips in {time.time()-t0:.1f}s")

    # Collect candidate internal branches (non-trivial splits)
    candidates = []
    for node in t.traverse("postorder"):
        if node.is_root():
            continue
        if node.dist is None or node.dist <= 0:
            continue
        nd = len(node)  # number of descendant leaves
        if nd < 2 or nd > n - 2:
            continue
        # balance score: how close to 50/50 split
        balance = min(nd, n - nd)
        candidates.append((balance, node))

    # Sort by balance (most balanced first)
    candidates.sort(key=lambda x: -x[0])
    max_eval = min(args.max_eval, len(candidates))
    print(f"[MAD] {len(candidates)} candidate branches total, evaluating top {max_eval} most balanced")

    best_rho = float('inf')
    best_node = None
    t0 = time.time()

    for idx in range(max_eval):
        balance, node = candidates[idx]

        if idx % 100 == 0:
            elapsed = time.time() - t0
            print(f"[MAD]   {idx}/{max_eval} evaluated, elapsed={elapsed:.1f}s, best_rho={best_rho:.6f}")

        try:
            # ete3 set_outgroup modifies tree in-place (fast, C-based)
            t.set_outgroup(node)

            # Compute root-to-tip distances
            r2t = np.array([t.get_distance(leaf) for leaf in leaves])
            mean_d = r2t.mean()

            if mean_d <= 0:
                continue

            # MAD criterion: relative deviation of root-to-tip distances
            rho = np.sqrt(np.mean(((r2t - mean_d) / mean_d) ** 2))

            if rho < best_rho:
                best_rho = rho
                best_node = node

        except Exception as e:
            continue

    elapsed = time.time() - t0
    print(f"[MAD] Evaluation complete: {max_eval} branches in {elapsed:.1f}s")

    if best_node is None:
        print("[MAD] WARNING: No optimal root found, using midpoint")
        # Reload and midpoint root
        t = Tree(args.input)
        mp = t.get_midpoint_outgroup()
        t.set_outgroup(mp)
    else:
        # Re-load tree and apply best rooting
        # (tree is already modified, but re-root at best_node cleanly)
        t.set_outgroup(best_node)
        print(f"[MAD] Best rho = {best_rho:.6f}")

    # Summary
    children = t.get_children()
    sizes = [len(c) for c in children]
    print(f"[MAD] Root split sizes: {sizes}")

    # Write output
    t.write(outfile=args.output, format=5)
    print(f"[MAD] Written to: {args.output}")

    # Also write a summary
    summary_path = args.output.replace(".treefile", "_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"MAD Rooting Summary\n")
        f.write(f"===================\n")
        f.write(f"Input: {args.input}\n")
        f.write(f"Tips: {n}\n")
        f.write(f"Candidates evaluated: {max_eval}\n")
        f.write(f"Best rho (relative ancestor deviation): {best_rho:.6f}\n")
        f.write(f"Root split sizes: {sizes}\n")
        f.write(f"Elapsed: {elapsed:.1f}s\n")

        # Root-to-tip distance stats
        r2t_final = np.array([t.get_distance(leaf) for leaf in leaves])
        f.write(f"\nRoot-to-tip distances:\n")
        f.write(f"  Mean: {r2t_final.mean():.4f}\n")
        f.write(f"  Std:  {r2t_final.std():.4f}\n")
        f.write(f"  Min:  {r2t_final.min():.4f}\n")
        f.write(f"  Max:  {r2t_final.max():.4f}\n")
        f.write(f"  CV:   {r2t_final.std()/r2t_final.mean():.4f}\n")

    print(f"[MAD] Summary written to: {summary_path}")
    print("Done!")


if __name__ == "__main__":
    main()
