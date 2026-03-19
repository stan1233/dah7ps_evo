#!/usr/bin/env python3
"""compare_trees_lite.py — Lightweight tree comparison (no external deps).

Parses newick trees with a minimal recursive descent parser,
computes RF distance and monophyly checks.
"""
import sys, os, re, datetime

# ── Minimal Newick parser ────────────────────────────────────────────────────
class Node:
    __slots__ = ['name', 'children', 'branch_length']
    def __init__(self, name="", children=None, branch_length=None):
        self.name = name
        self.children = children or []
        self.branch_length = branch_length

    def is_leaf(self):
        return len(self.children) == 0

    def get_leaves(self):
        if self.is_leaf():
            return [self.name]
        leaves = []
        for c in self.children:
            leaves.extend(c.get_leaves())
        return leaves


def parse_newick(s):
    """Parse a newick string into a tree of Node objects."""
    s = s.strip().rstrip(';').strip()
    pos = [0]  # mutable index

    def _parse():
        children = []
        if pos[0] < len(s) and s[pos[0]] == '(':
            pos[0] += 1  # skip '('
            children.append(_parse())
            while pos[0] < len(s) and s[pos[0]] == ',':
                pos[0] += 1  # skip ','
                children.append(_parse())
            if pos[0] < len(s) and s[pos[0]] == ')':
                pos[0] += 1  # skip ')'

        # Parse name and branch length
        name = ""
        while pos[0] < len(s) and s[pos[0]] not in (',', ')', ':', ';'):
            name += s[pos[0]]
            pos[0] += 1
        bl = None
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            bl_str = ""
            while pos[0] < len(s) and s[pos[0]] not in (',', ')', ';'):
                bl_str += s[pos[0]]
                pos[0] += 1
            try:
                bl = float(bl_str)
            except ValueError:
                pass

        return Node(name=name.strip(), children=children, branch_length=bl)

    return _parse()


def get_bipartitions(root):
    """Get all bipartitions (as frozenset pairs) from an unrooted tree."""
    all_leaves = frozenset(root.get_leaves())
    bips = set()

    def _collect(node):
        if node.is_leaf():
            return frozenset([node.name])
        clade_set = frozenset()
        for c in node.children:
            clade_set = clade_set | _collect(c)
        # Add bipartition (skip root-level trivial ones)
        complement = all_leaves - clade_set
        if len(clade_set) > 0 and len(complement) > 0 and \
           len(clade_set) > 1 and len(complement) > 1:
            bip = frozenset([clade_set, complement])
            bips.add(bip)
        return clade_set

    _collect(root)
    return bips


def check_monophyly(root, target_names):
    """Check if target_names form a monophyletic clade."""
    all_leaves = set(root.get_leaves())
    present = set(n for n in target_names if n in all_leaves)
    if len(present) < 2:
        return {"n": len(present), "mono": "NA", "type": "NA"}

    # Find the smallest clade containing all target leaves
    def _find_mrca(node):
        if node.is_leaf():
            return (frozenset([node.name]), node)
        child_results = [_find_mrca(c) for c in node.children]
        combined = frozenset()
        for cr, _ in child_results:
            combined = combined | cr
        # Check if any single child contains all targets
        for cr, cn in child_results:
            if present.issubset(cr):
                return (cr, cn)
        # This node is the MRCA
        return (combined, node)

    mrca_leaves, mrca_node = _find_mrca(root)
    mrca_leaf_set = set(mrca_node.get_leaves())

    if mrca_leaf_set == present:
        return {"n": len(present), "mono": True, "type": "monophyletic"}
    extra = mrca_leaf_set - present
    ratio = len(present) / len(mrca_leaf_set)
    if ratio > 0.9:
        return {"n": len(present), "mono": False,
                "type": f"paraphyletic ({len(extra)} intruders, {ratio:.0%})"}
    return {"n": len(present), "mono": False,
            "type": f"polyphyletic ({len(extra)} intruders, {ratio:.0%})"}


# ── PDB subtype mapping ─────────────────────────────────────────────────────
PDB_SUBTYPE = {
    "PDB-1KFL_A": "Ia", "PDB-1KFL_B": "Ia", "PDB-1KFL_C": "Ia",
    "PDB-1KFL_D": "Ia", "PDB-1KFL_E": "Ia", "PDB-1KFL_F": "Ia",
    "PDB-1KFL_G": "Ia", "PDB-1KFL_H": "Ia",
    "PDB-1RZM_A": "Ib", "PDB-1RZM_B": "Ib",
    "PDB-3NV8_A": "II", "PDB-3NV8_B": "II",
    "PDB-5CKV_A": "II", "PDB-5CKV_B": "II",
    "PDB-2B7O_A": "II", "PDB-2B7O_B": "II",
}


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Compare AA vs 3Di skeleton trees")
    parser.add_argument("--aa", required=True)
    parser.add_argument("--threedi", required=True)
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument("--out_md", required=True)
    parser.add_argument("--aa_model", default="Q.PFAM+I+R4")
    parser.add_argument("--threedi_model", default="Q.3Di.AF+G4")
    args = parser.parse_args()

    for p in [args.aa, args.threedi]:
        if not os.path.isfile(p):
            print(f"ERROR: {p} not found", file=sys.stderr); sys.exit(1)

    print(f"Loading AA tree: {args.aa}")
    with open(args.aa) as f:
        aa_tree = parse_newick(f.read())
    aa_tips = aa_tree.get_leaves()
    print(f"  Tips: {len(aa_tips)}")

    print(f"Loading 3Di tree: {args.threedi}")
    with open(args.threedi) as f:
        di_tree = parse_newick(f.read())
    di_tips = di_tree.get_leaves()
    print(f"  Tips: {len(di_tips)}")

    # Tip match
    aa_set, di_set = set(aa_tips), set(di_tips)
    if aa_set != di_set:
        print(f"WARNING: mismatch! AA-only={aa_set-di_set}, 3Di-only={di_set-aa_set}",
              file=sys.stderr)

    # RF distance
    print("Computing RF distance ...")
    b_aa = get_bipartitions(aa_tree)
    b_di = get_bipartitions(di_tree)
    shared = b_aa & b_di
    only_aa = b_aa - b_di
    only_di = b_di - b_aa
    rf = len(only_aa) + len(only_di)
    max_rf = len(b_aa) + len(b_di)
    nrf = rf / max_rf if max_rf > 0 else 0
    print(f"  RF = {rf}/{max_rf} (nRF = {nrf:.4f})")
    print(f"  Shared = {len(shared)}, AA-only = {len(only_aa)}, 3Di-only = {len(only_di)}")

    # Monophyly
    ia_tips = [t for t in aa_tips if PDB_SUBTYPE.get(t) == "Ia"]
    ib_tips = [t for t in aa_tips if PDB_SUBTYPE.get(t) == "Ib"]
    ii_tips = [t for t in aa_tips if PDB_SUBTYPE.get(t) == "II"]
    mono_results = []
    for label, tips in [("Type_Ia (PDB-1KFL)", ia_tips),
                        ("Type_Ib (PDB-1RZM)", ib_tips),
                        ("Type_II (PDB-3NV8/5CKV/2B7O)", ii_tips)]:
        aa_m = check_monophyly(aa_tree, tips)
        di_m = check_monophyly(di_tree, tips)
        mono_results.append((label, aa_m, di_m))
        print(f"  {label}: AA={aa_m['mono']}, 3Di={di_m['mono']}")

    # ── Write TSV ──
    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)
    with open(args.out_tsv, "w") as f:
        f.write("metric\tvalue\tnote\n")
        f.write(f"aa_tree_bipartitions\t{len(b_aa)}\tInternal bipartitions in AA tree\n")
        f.write(f"threedi_tree_bipartitions\t{len(b_di)}\tInternal bipartitions in 3Di tree\n")
        f.write(f"shared_bipartitions\t{len(shared)}\tIdentical in both\n")
        f.write(f"aa_only_bipartitions\t{len(only_aa)}\tUnique to AA\n")
        f.write(f"threedi_only_bipartitions\t{len(only_di)}\tUnique to 3Di\n")
        f.write(f"RF_distance\t{rf}\tRobinson-Foulds\n")
        f.write(f"RF_max\t{max_rf}\tMax possible\n")
        f.write(f"RF_normalized\t{nrf:.4f}\t0=identical 1=max different\n")
    print(f"TSV → {args.out_tsv}")

    # ── Write MD ──
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    q1_consistent = all(aa_m['mono'] == di_m['mono'] for _, aa_m, di_m in mono_results)
    with open(args.out_md, "w") as f:
        f.write(f"# Phase 4.2: AA vs 3Di Skeleton Tree Comparison\n\n")
        f.write(f"> Generated: {now}\n\n---\n\n")
        f.write("## Overview\n\n| Property | AA Tree | 3Di Tree |\n|---|---|---|\n")
        f.write(f"| Input | `skeleton_core_aa.fa` | `skeleton_core_3di.fa` |\n")
        f.write(f"| Model | {args.aa_model} | {args.threedi_model} |\n")
        f.write(f"| Tips | {len(aa_tips)} | {len(di_tips)} |\n\n")

        f.write("## Topological Comparison\n\n| Metric | Value |\n|---|---|\n")
        f.write(f"| Shared bipartitions | {len(shared)} |\n")
        f.write(f"| AA-only bipartitions | {len(only_aa)} |\n")
        f.write(f"| 3Di-only bipartitions | {len(only_di)} |\n")
        f.write(f"| Robinson-Foulds distance | {rf} / {max_rf} |\n")
        f.write(f"| Normalized RF | {nrf:.4f} |\n\n")

        f.write("## Q1: Are subtype deep divergence directions consistent?\n\n")
        f.write("| Clade | AA monophyletic? | 3Di monophyletic? | Consistent? |\n")
        f.write("|---|---|---|---|\n")
        for label, aa_m, di_m in mono_results:
            c = "✅" if aa_m['mono'] == di_m['mono'] else "⚠️"
            f.write(f"| {label} | {aa_m['mono']} ({aa_m['type']}) "
                    f"| {di_m['mono']} ({di_m['type']}) | {c} |\n")
        f.write("\n")
        if q1_consistent:
            f.write("> **Conclusion:** PDB-based subtype groupings show consistent "
                    "monophyly patterns between AA and 3Di trees.\n\n")
        else:
            f.write("> **Conclusion:** Some subtype groupings differ between AA and 3Di trees.\n\n")

        f.write("## Q2: KDOPS/outgroup separation compatibility\n\n")
        f.write("Skeleton trees contain the structural panel only (no KDOPS). "
                "Direct KDOPS compatibility test not applicable.\n\n")

        f.write("## Q3: Focal deep split conflicts\n\n")
        f.write(f"With nRF = **{nrf:.4f}**, ")
        if nrf < 0.3:
            f.write("the two trees are broadly consistent.\n\n")
        elif nrf < 0.6:
            f.write("the trees show moderate differences.\n\n")
        else:
            f.write("the trees show substantial disagreement at deep nodes.\n\n")

        verdict = "GREEN" if nrf < 0.5 and q1_consistent else \
                  "YELLOW" if nrf < 0.7 else "RED"
        f.write(f"---\n\n## Summary for QC3\n\n**4.2 verdict: QC3-{verdict}**\n\n")
        if verdict == "GREEN":
            f.write("AA and 3Di trees do not produce fatal contradictions at the deep "
                    "topology level.\n")
        elif verdict == "YELLOW":
            f.write("Some discrepancies exist. Deep-level conclusions need caution.\n")
        else:
            f.write("Serious conflicts detected. Deep inferences undermined.\n")

    print(f"MD  → {args.out_md}")
    print("Done.")


if __name__ == "__main__":
    main()
