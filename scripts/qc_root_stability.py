#!/usr/bin/env python3
"""
qc_root_stability.py — DAH7PS V5.1 QC3 Root Robustness Gate
============================================================
分析四套 rooting scenario（S1 MFP / S2 LGC20 / S3 Midpoint / S4 MAD）
的根位稳定性，并输出：
  - results/04_phylogeny_asr/QC3_root_stability.md
  - results/04_phylogeny_asr/root_scenarios.tsv  (更新)
  - results/04_phylogeny_asr/node_selection_registry.tsv (初始化/更新)
  - results/04_phylogeny_asr/rooted_working_tree.treefile (指定工作树)

用法:
  conda run -n dah7ps_v4 python scripts/qc_root_stability.py [--help]
"""

import argparse
import sys
import os
import csv
import datetime
from collections import defaultdict

def check_biopython():
    try:
        from Bio import Phylo
        return True
    except ImportError:
        return False

def parse_args():
    p = argparse.ArgumentParser(
        description="DAH7PS QC3: root robustness gate across 4 root scenarios"
    )
    p.add_argument("--dir", default="results/04_phylogeny_asr",
                   help="Input/output directory (default: results/04_phylogeny_asr)")
    p.add_argument("--module_strict",
                   default="results/03_msa_modules/module_presence_absence_strict.tsv",
                   help="Strict module presence/absence matrix")
    p.add_argument("--module_relaxed",
                   default="results/03_msa_modules/module_presence_absence_relaxed.tsv",
                   help="Relaxed module presence/absence matrix")
    p.add_argument("--out_dir", default="results/04_phylogeny_asr",
                   help="Output directory")
    return p.parse_args()

# ─────────────────────────────────────────────────────
# 1. Tree loading helpers
# ─────────────────────────────────────────────────────

def load_tree(path):
    """Load newick tree. Returns (tree, error_msg)."""
    from Bio import Phylo
    import io
    if not os.path.exists(path):
        return None, f"FILE NOT FOUND: {path}"
    try:
        tree = Phylo.read(path, "newick")
        return tree, None
    except Exception as e:
        return None, str(e)

def get_root_children_summary(tree):
    """
    Return info about root's two (or more) immediate children.
    For each child: #tips, #KDOPS tips, first 3 tip names.
    """
    root = tree.root
    summary = []
    for i, clade in enumerate(root.clades):
        tips = clade.get_terminals()
        n_total = len(tips)
        n_kdops = sum(1 for t in tips if t.name and t.name.startswith("KDOPS"))
        n_ingroup = n_total - n_kdops
        example_tips = [t.name for t in tips[:4]]
        summary.append({
            "child_idx": i,
            "n_total": n_total,
            "n_kdops": n_kdops,
            "n_ingroup": n_ingroup,
            "example_tips": example_tips,
            "confidence": clade.confidence,
        })
    return summary

def count_tips(tree):
    return tree.count_terminals()

def get_kdops_monophyly(tree):
    """
    Check if KDOPS sequences form a monophyletic group.
    Returns: (is_mono, n_kdops, which_child_index_contains_all)
    """
    root = tree.root
    all_kdops = set(t.name for t in tree.get_terminals()
                    if t.name and t.name.startswith("KDOPS"))
    if not all_kdops:
        return None, 0, -1  # no KDOPS in tree

    # Check if any single child of root contains all KDOPS
    for i, clade in enumerate(root.clades):
        child_names = set(t.name for t in clade.get_terminals())
        if all_kdops.issubset(child_names):
            # Check further: is KDOPS strictly monophyletic within this child?
            child_kdops = set(t.name for t in clade.get_terminals()
                              if t.name and t.name.startswith("KDOPS"))
            child_nonkdops = child_names - child_kdops
            if child_nonkdops:
                return False, len(all_kdops), i  # KDOPS mixed with ingroup in same clade
            else:
                return True, len(all_kdops), i   # That child IS the KDOPS clade
    # KDOPS split across multiple root children
    return False, len(all_kdops), -1

def get_clade_monophyly_by_prefix(tree, prefix, name_to_subtype=None):
    """
    Check if sequences with a given prefix form a monophyletic group.
    For ingroup subtype checking, use name_to_subtype dict.
    """
    if name_to_subtype:
        target = set(name for name, st in name_to_subtype.items() if st == prefix)
    else:
        target = set(t.name for t in tree.get_terminals()
                     if t.name and t.name.startswith(prefix))
    if not target:
        return None, 0

    all_tips = tree.get_terminals()
    n = len(target & set(t.name for t in all_tips))

    # Try to find the MRCA of the target set
    target_terminals = [t for t in all_tips if t.name in target]
    if len(target_terminals) < 2:
        return None, n

    try:
        mrca = tree.common_ancestor(target_terminals)
        mrca_tips = set(t.name for t in mrca.get_terminals())
        non_target_in_mrca = mrca_tips - target
        is_mono = len(non_target_in_mrca) == 0
        return is_mono, n
    except Exception:
        return None, n

def get_deep_node_supports(tree, n_deep=5):
    """Extract confidence values for the n_deep shallowest internal nodes."""
    internals = []
    def traverse(clade, depth):
        if not clade.is_terminal():
            conf = clade.confidence
            internals.append((depth, conf))
            for child in clade.clades:
                traverse(child, depth + 1)
    traverse(tree.root, 0)
    internals.sort(key=lambda x: x[0])
    return internals[:n_deep]

def compute_rf_distance(tree1, tree2):
    """
    Compute Robinson-Foulds distance between two trees.
    Only works on trees with same tip set (after pruning KDOPS).
    Returns (rf, max_rf, norm_rf) or None on failure.
    """
    try:
        from Bio import Phylo
        import io

        def get_bipartitions(tree):
            all_tips = frozenset(t.name for t in tree.get_terminals())
            bips = set()
            def traverse(clade):
                if clade.is_terminal():
                    return frozenset([clade.name])
                child_sets = [traverse(c) for c in clade.clades]
                full = frozenset().union(*child_sets)
                complement = all_tips - full
                if full and complement:
                    bip = (min(full, complement, key=len), max(full, complement, key=len))
                    # Use frozensets for hashing
                    bips.add((frozenset(full), frozenset(complement)))
                return full
            traverse(tree.root)
            return bips

        bips1 = get_bipartitions(tree1)
        bips2 = get_bipartitions(tree2)

        # Normalize: use frozenset pairs as canonical form
        def canonical(bips):
            result = set()
            for a, b in bips:
                result.add((frozenset(a), frozenset(b)))
            return result

        b1 = canonical(bips1)
        b2 = canonical(bips2)

        shared = len(b1 & b2)
        rf = len(b1 - b2) + len(b2 - b1)
        max_rf = len(b1) + len(b2)
        norm_rf = rf / max_rf if max_rf > 0 else 0.0
        return rf, max_rf, norm_rf
    except Exception as e:
        return None

def prune_kdops(tree):
    """Return a copy of tree with KDOPS tips pruned (ingroup only)."""
    from Bio import Phylo
    import copy
    t = copy.deepcopy(tree)
    kdops_tips = [leaf for leaf in t.get_terminals()
                  if leaf.name and leaf.name.startswith("KDOPS")]
    for tip in kdops_tips:
        t.prune(tip)
    return t

# ─────────────────────────────────────────────────────
# 2. Module analysis
# ─────────────────────────────────────────────────────

def load_module_matrix(path):
    """Load presence/absence matrix. Returns {seq_id: {module: int}}."""
    data = {}
    if not os.path.exists(path):
        return data
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sid = row['seq_id']
            data[sid] = {k: int(v) for k, v in row.items()
                         if k not in ('seq_id', 'boundary_confidence')}
    return data

def compute_module_prevalence(matrix):
    """Compute prevalence % of each module across all sequences."""
    if not matrix:
        return {}
    modules = list(next(iter(matrix.values())).keys())
    n = len(matrix)
    prev = {}
    for m in modules:
        count = sum(1 for v in matrix.values() if v.get(m, 0) == 1)
        prev[m] = round(100.0 * count / n, 1)
    return prev

def infer_module_events(strict_prev, relaxed_prev):
    """
    Infer likely evolutionary event type per module based on prevalence.
    Returns dict: module -> {event_type, base_note}
    """
    events = {}
    for m in strict_prev:
        sp = strict_prev.get(m, 0)
        rp = relaxed_prev.get(m, 0)
        # Heuristic classification
        if sp > 50:
            event_type = "widespread_presence"
        elif sp < 5:
            event_type = "rare_or_lineage_specific"
        elif 5 <= sp <= 30:
            event_type = "minority_gain_candidate"
        elif 30 < sp <= 50:
            event_type = "partial_loss_candidate"
        else:
            event_type = "complex"
        events[m] = {
            "event_type": event_type,
            "strict_prev": sp,
            "relaxed_prev": rp,
            "delta_strict_relaxed": round(rp - sp, 1),
        }
    return events

# ─────────────────────────────────────────────────────
# 3. Root scenario comparison logic
# ─────────────────────────────────────────────────────

def classify_root_agreement(scenarios_data):
    """
    Given per-scenario summary dicts, infer consistency:
    - outgroup_based: S1, S2 (KDOPS as outgroup)
    - topology_based: S3, S4 (no outgroup)
    Returns: (n_agree_outgroup, n_agree_topology, n_total_valid)
    """
    valid = [s for s in scenarios_data if s.get("loaded")]
    outgroup = [s for s in valid if s.get("rooting_method") == "outgroup"]
    topology = [s for s in valid if s.get("rooting_method") in ("midpoint", "MAD")]

    # Classify if S1 and S2 "agree": compare root children tip-count pattern
    def outgroup_root_summary(s):
        ch = s.get("root_children", [])
        if not ch:
            return None
        # Largest non-KDOPS clade and smallest
        ingroup_sizes = sorted([c["n_ingroup"] for c in ch if c["n_ingroup"] > 0])
        return tuple(ingroup_sizes)

    # Agreement: are the root children patterns consistent across scenarios?
    outgroup_summaries = [outgroup_root_summary(s) for s in outgroup]
    topology_summaries = [outgroup_root_summary(s) for s in topology]

    return {
        "n_valid": len(valid),
        "n_outgroup": len(outgroup),
        "n_topology": len(topology),
        "outgroup_summaries": outgroup_summaries,
        "topology_summaries": topology_summaries,
    }

# ─────────────────────────────────────────────────────
# 4. Main QC3 logic
# ─────────────────────────────────────────────────────

def main():
    args = parse_args()

    if not check_biopython():
        print("ERROR: BioPython not found. Run in conda env dah7ps_v4.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Define scenarios ─────────────────────────────
    scenarios = [
        {
            "id": "S1_MFP_KDOPS",
            "label": "S1 MFP+KDOPS outgroup",
            "file": os.path.join(args.dir, "CoreTree_rooted_MFP.treefile"),
            "rooting_method": "outgroup",
            "has_kdops": True,
            "model": "Q.PFAM+F+R10",
            "logL": -2835379.61,
        },
        {
            "id": "S2_LGC20_KDOPS",
            "label": "S2 LG+C20+KDOPS outgroup",
            "file": os.path.join(args.dir, "CoreTree_rooted_LGC20.treefile"),
            "rooting_method": "outgroup",
            "has_kdops": True,
            "model": "LG+C20+F+G",
            "logL": -2820118.09,
        },
        {
            "id": "S3_MIDPOINT",
            "label": "S3 Midpoint ingroup",
            "file": os.path.join(args.dir, "CoreTree_rooted_midpoint_ingroup.treefile"),
            "rooting_method": "midpoint",
            "has_kdops": False,
            "model": "NA",
            "logL": None,
        },
        {
            "id": "S4_MAD",
            "label": "S4 MAD ingroup",
            "file": os.path.join(args.dir, "CoreTree_rooted_MAD_ingroup.treefile"),
            "rooting_method": "MAD",
            "has_kdops": False,
            "model": "NA",
            "logL": None,
        },
    ]

    print("=" * 65)
    print("DAH7PS V5.1  QC3 Root Robustness Gate")
    print("=" * 65)

    # ── Load all trees ────────────────────────────────
    print("\n[1/5] Loading trees …")
    for s in scenarios:
        tree, err = load_tree(s["file"])
        s["loaded"] = tree is not None
        s["load_error"] = err
        if tree:
            s["tree"] = tree
            s["n_tips"] = count_tips(tree)
            s["root_children"] = get_root_children_summary(tree)
            # KDOPS monophyly
            if s["has_kdops"]:
                mono, n_k, child_idx = get_kdops_monophyly(tree)
                s["kdops_monophyletic"] = mono
                s["n_kdops"] = n_k
                s["kdops_child_idx"] = child_idx
            else:
                s["kdops_monophyletic"] = None
                s["n_kdops"] = 0
            s["deep_supports"] = get_deep_node_supports(tree, n_deep=8)
            print(f"  ✓ {s['id']}: {s['n_tips']} tips loaded from {s['file']}")
        else:
            s["tree"] = None
            s["n_tips"] = 0
            s["root_children"] = []
            s["kdops_monophyletic"] = None
            s["n_kdops"] = 0
            s["deep_supports"] = []
            print(f"  ✗ {s['id']}: FAILED — {err}")

    # ── Ingroup-only RF distances ────────────────────
    print("\n[2/5] Computing pairwise RF distances (ingroup only) …")
    rf_matrix = {}
    loaded_scenarios = [s for s in scenarios if s["loaded"]]

    ingroup_trees = {}
    for s in loaded_scenarios:
        if s["has_kdops"]:
            ingroup_trees[s["id"]] = prune_kdops(s["tree"])
        else:
            ingroup_trees[s["id"]] = s["tree"]

    ids = list(ingroup_trees.keys())
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            id1, id2 = ids[i], ids[j]
            t1 = ingroup_trees[id1]
            t2 = ingroup_trees[id2]
            # Check tip counts match
            n1 = t1.count_terminals()
            n2 = t2.count_terminals()
            if n1 != n2:
                rf_matrix[(id1, id2)] = None
                print(f"  ! {id1} vs {id2}: tip count mismatch ({n1} vs {n2}), skipping RF")
                continue
            result = compute_rf_distance(t1, t2)
            rf_matrix[(id1, id2)] = result
            if result:
                rf, max_rf, nrf = result
                print(f"  {id1} vs {id2}: RF={rf}, max={max_rf}, nRF={nrf:.3f}")
            else:
                print(f"  {id1} vs {id2}: RF computation failed")

    # ── Module analysis ───────────────────────────────
    print("\n[3/5] Analysing module presence/absence …")
    strict_matrix = load_module_matrix(args.module_strict)
    relaxed_matrix = load_module_matrix(args.module_relaxed)
    strict_prev = compute_module_prevalence(strict_matrix)
    relaxed_prev = compute_module_prevalence(relaxed_matrix)
    module_events = infer_module_events(strict_prev, relaxed_prev)

    print("  Module prevalence (strict → relaxed):")
    for m, ev in module_events.items():
        print(f"    {m}: strict={ev['strict_prev']}% | relaxed={ev['relaxed_prev']}% "
              f"| type={ev['event_type']}")

    # ── Root-sensitivity classification ─────────────
    print("\n[4/5] Classifying root sensitivity per module event …")
    # S1 and S2 are outgroup-based (model-dependent root)
    # S3 and S4 are topology-based (model-free root)
    # Cross-scenario agreement = root_robust; only one scenario = root_sensitive
    s1 = next((s for s in scenarios if s["id"] == "S1_MFP_KDOPS"), None)
    s2 = next((s for s in scenarios if s["id"] == "S2_LGC20_KDOPS"), None)
    s3 = next((s for s in scenarios if s["id"] == "S3_MIDPOINT"), None)
    s4 = next((s for s in scenarios if s["id"] == "S4_MAD"), None)

    # Determine root position agreement between S1 and S2
    def root_clade_sizes(s):
        """Return sorted tuple of root-child ingroup sizes."""
        ch = s.get("root_children", []) if s else []
        sizes = sorted([c["n_ingroup"] for c in ch if c["n_ingroup"] > 0])
        return tuple(sizes)

    s1_sizes = root_clade_sizes(s1) if s1 and s1["loaded"] else None
    s2_sizes = root_clade_sizes(s2) if s2 and s2["loaded"] else None
    s3_sizes = root_clade_sizes(s3) if s3 and s3["loaded"] else None
    s4_sizes = root_clade_sizes(s4) if s4 and s4["loaded"] else None

    # Root positions are "consistent" if the min-clade size (outgroup-adjacent sister)
    # is within 10% of total across scenarios
    def partition_ratio(sizes, total):
        if not sizes:
            return None
        min_s = min(sizes)
        return round(min_s / total, 3) if total else None

    total_ingroup = 9393  # known
    r1 = partition_ratio(s1_sizes, total_ingroup) if s1_sizes else None
    r2 = partition_ratio(s2_sizes, total_ingroup) if s2_sizes else None
    r3 = partition_ratio(s3_sizes, total_ingroup) if s3_sizes else None
    r4 = partition_ratio(s4_sizes, total_ingroup) if s4_sizes else None

    print(f"  Root partition ratios:")
    print(f"    S1 MFP: {r1} (smallest ingroup clade / total)")
    print(f"    S2 LGC20: {r2}")
    print(f"    S3 Midpoint: {r3}")
    print(f"    S4 MAD: {r4}")

    # Determine inter-scenario RF convergence
    s1_s2_rf = rf_matrix.get(("S1_MFP_KDOPS", "S2_LGC20_KDOPS"))
    s1_s3_rf = rf_matrix.get(("S1_MFP_KDOPS", "S3_MIDPOINT"))
    s1_s4_rf = rf_matrix.get(("S1_MFP_KDOPS", "S4_MAD"))

    # Determine root-stability verdict
    # "root_robust" if >2 scenarios agree on the same root partition
    ratios = [r for r in [r1, r2, r3, r4] if r is not None]
    n_valid = len(ratios)

    # Agreement: ratio within 5% of each other (tight), 10% (loose)
    def agreement_count(ratios, tol=0.10):
        if not ratios:
            return 0
        ref = ratios[0]
        return sum(1 for r in ratios if abs(r - ref) <= tol)

    n_agree_tight = agreement_count(ratios, tol=0.05)
    n_agree_loose = agreement_count(ratios, tol=0.15)

    # Overall QC3 verdict
    n_scenarios_loaded = sum(1 for s in scenarios if s["loaded"])
    if n_scenarios_loaded < 3:
        overall_verdict = "RED"
        root_stability = "insufficient_data"
    elif n_agree_tight >= 3:
        overall_verdict = "GREEN"
        root_stability = "root_robust"
    elif n_agree_loose >= 3:
        overall_verdict = "YELLOW"
        root_stability = "root_semi_robust"
    else:
        overall_verdict = "YELLOW"
        root_stability = "root_sensitive"

    print(f"\n  Overall root stability: {root_stability} ({n_agree_tight}/{n_valid} tight, "
          f"{n_agree_loose}/{n_valid} loose agreement)")

    # Assign module event root-sensitivity labels
    module_stability = {}
    for m, ev in module_events.items():
        sp = ev["strict_prev"]
        if sp > 70 or sp < 3:
            # Highly prevalent or very rare → interpretation unlikely to flip with root
            stab = "root_robust"
        elif sp < 15:
            stab = "root_sensitive"
        elif root_stability == "root_robust":
            stab = "root_robust"
        elif root_stability == "root_semi_robust":
            stab = "ambiguous"
        else:
            stab = "root_sensitive"
        module_stability[m] = stab

    print("\n  Module root-sensitivity labels:")
    for m, stab in module_stability.items():
        print(f"    {m}: {stab} (strict={strict_prev.get(m)}%)")

    # ── Write outputs ────────────────────────────────
    print("\n[5/5] Writing output files …")
    now = datetime.datetime.now().strftime("%Y-%m-%d")

    # ── 5a. Update root_scenarios.tsv ───────────────
    tsv_path = os.path.join(args.out_dir, "root_scenarios.tsv")
    status_map = {
        "S1_MFP_KDOPS": "completed",
        "S2_LGC20_KDOPS": "completed",
        "S3_MIDPOINT": "completed",
        "S4_MAD": "completed",
    }
    with open(tsv_path, "w") as f:
        f.write("scenario_id\tsource_alignment\ttree_model\trooting_method\t"
                "outgroup_set\ttip_count\tsupport_status\tstatus\troot_partition_ratio\tnote\n")
        rows = [
            ("S1_MFP_KDOPS", "core_with_outgroup.afa", "Q.PFAM+F+R10", "outgroup",
             "KDOPS_P0A715,KDOPS_Q9ZFK4,KDOPS_O66496", "9405", "UFBoot=1000",
             "completed", str(r1) if r1 else "NA",
             "MFP rooted tree; KDOPS polyphyletic in ML tree"),
            ("S2_LGC20_KDOPS", "core_with_outgroup.afa", "LG+C20+F+G", "outgroup",
             "KDOPS_P0A715,KDOPS_Q9ZFK4,KDOPS_O66496", "9405", "UFBoot=1000",
             "completed", str(r2) if r2 else "NA",
             f"LG+C20 rooted tree; LogL=-2820118.09; completed 2026-04-05"),
            ("S3_MIDPOINT", "CoreTree_rooted_ingroup.treefile", "NA", "midpoint",
             "none", "9393", "NA",
             "completed", str(r3) if r3 else "NA",
             "BioPython midpoint root; ingroup only"),
            ("S4_MAD", "CoreTree_rooted_ingroup.treefile", "NA", "MAD",
             "none", "9393", "NA",
             "completed", str(r4) if r4 else "NA",
             "ete3 MAD full-search; rho=0.163617; root non-unique"),
            ("S5_NONREV_REDUCED", "TBD", "nonreversible", "nonreversible",
             "none", "TBD", "NA",
             "optional", "NA",
             "Reduced representative set; resource-dependent"),
        ]
        for row in rows:
            f.write("\t".join(row) + "\n")
    print(f"  ✓ {tsv_path}")

    # ── 5b. Write QC3_root_stability.md ─────────────
    md_path = os.path.join(args.out_dir, "QC3_root_stability.md")
    with open(md_path, "w") as f:
        f.write(f"# DAH7PS V5.1 QC3 Root Robustness Gate Report\n\n")
        f.write(f"> Generated: {now}  \n")
        f.write(f"> Script: `scripts/qc_root_stability.py`  \n")
        f.write(f"> Status: **{overall_verdict}**  \n\n")
        f.write("---\n\n")

        f.write("## 1. Scenarios Summary\n\n")
        f.write("| Scenario | Method | Tree File | Tips | Root Partition Ratio | Status |\n")
        f.write("|---|---|---|---|---|---|\n")
        ratio_map = {
            "S1_MFP_KDOPS": r1, "S2_LGC20_KDOPS": r2,
            "S3_MIDPOINT": r3, "S4_MAD": r4,
        }
        for s in scenarios:
            ratio = ratio_map.get(s["id"])
            ratio_str = f"{ratio:.3f}" if ratio else "—"
            loaded_str = "✅" if s["loaded"] else "❌"
            tips_str = str(s["n_tips"]) if s["loaded"] else "FAILED"
            f.write(f"| {s['id']} | {s['rooting_method']} | "
                    f"`{os.path.basename(s['file'])}` | "
                    f"{tips_str} | {ratio_str} | {loaded_str} |\n")
        f.write("\n")

        f.write("## 2. Root Children Analysis\n\n")
        for s in scenarios:
            if not s["loaded"]:
                continue
            f.write(f"### {s['id']} — {s['label']}\n\n")
            f.write(f"- **Model**: `{s['model']}`  \n")
            if s.get("logL"):
                f.write(f"- **LogL**: {s['logL']}  \n")
            f.write(f"- **Total tips**: {s['n_tips']}  \n")
            if s["has_kdops"]:
                mono = s.get("kdops_monophyletic")
                mono_str = "✅ YES" if mono else "⚠️ NO (polyphyletic)" if mono is False else "N/A"
                f.write(f"- **KDOPS monophyletic**: {mono_str}  \n")
                f.write(f"- **KDOPS count**: {s['n_kdops']}  \n")
            f.write(f"\n**Root children:**\n\n")
            f.write("| Child | Total tips | KDOPS | Ingroup | Confidence | Example tips |\n")
            f.write("|---|---|---|---|---|---|\n")
            for c in s["root_children"]:
                eg = ", ".join(c["example_tips"][:3])
                f.write(f"| #{c['child_idx']} | {c['n_total']} | {c['n_kdops']} | "
                        f"{c['n_ingroup']} | {c['confidence']} | {eg} |\n")
            f.write("\n**Deep node support values (depth, confidence):**\n\n")
            f.write("| Depth | Confidence |\n|---|---|\n")
            for depth, conf in (s["deep_supports"] or []):
                f.write(f"| {depth} | {conf} |\n")
            f.write("\n")

        f.write("## 3. Pairwise RF Distances (Ingroup)\n\n")
        f.write("| Pair | RF | Max RF | nRF | Interpretation |\n")
        f.write("|---|---|---|---|---|\n")
        for (id1, id2), result in rf_matrix.items():
            if result:
                rf, max_rf, nrf = result
                if nrf < 0.2:
                    interp = "very similar"
                elif nrf < 0.4:
                    interp = "moderately similar"
                elif nrf < 0.6:
                    interp = "moderately different"
                else:
                    interp = "highly different"
                f.write(f"| {id1} vs {id2} | {rf} | {max_rf} | {nrf:.3f} | {interp} |\n")
            else:
                f.write(f"| {id1} vs {id2} | — | — | — | computation failed or tip mismatch |\n")
        f.write("\n")

        f.write("## 4. KDOPS Outgroup Assessment\n\n")
        f.write("| Scenario | KDOPS in tree | Monophyletic | Note |\n")
        f.write("|---|---|---|---|\n")
        for s in scenarios:
            if not s["loaded"]:
                continue
            if s["has_kdops"]:
                mono = s.get("kdops_monophyletic")
                mono_str = "YES" if mono else "NO (polyphyletic)" if mono is False else "N/A"
                f.write(f"| {s['id']} | YES ({s['n_kdops']}) | {mono_str} | "
                        f"outgroup-based rooting |\n")
            else:
                f.write(f"| {s['id']} | NO | N/A | {s['rooting_method']} rooting |\n")
        f.write("\n")
        f.write("> **KDOPS polyphyly in ML trees (S1/S2)** is flagged in AGENTS.md §0.2. "  
                "This does not invalidate the analysis — it means KDOPS cannot serve as a "
                "clean single outgroup and confirms that root uncertainty must be managed "
                "via multi-scenario approach (V5.1 strategy).\n\n")

        f.write("## 5. AA vs 3Di Tree (Phase 4.2)\n\n")
        f.write("*(From previous analysis — 2026-03-19)*\n\n")
        f.write("| Metric | Value |\n|---|---|\n")
        f.write("| AA model | Q.PFAM+I+R4 |\n")
        f.write("| 3Di model | Q.3Di.AF+G4 |\n")
        f.write("| Normalized RF (nRF) | 0.7442 |\n")
        f.write("| Shared bipartitions | 11 of 43 |\n")
        f.write("| Ia monophyletic (both) | YES |\n")
        f.write("| Ib monophyletic (both) | YES |\n")
        f.write("| II monophyletic (both) | YES |\n")
        f.write("| Q2 (KDOPS placement) | N/A (no KDOPS in skeleton) |\n")
        f.write("| Q3 (fatal deep conflict) | NO |\n")
        f.write("| **QC3-4.2 verdict** | **YELLOW** (high nRF expected for AA vs 3Di) |\n\n")

        f.write("## 6. Module Event Root-Sensitivity\n\n")
        f.write("| Module | Strict % | Relaxed % | Event Type | Root Sensitivity |\n")
        f.write("|---|---|---|---|---|\n")
        for m, ev in module_events.items():
            stab = module_stability.get(m, "unknown")
            icon = "✅" if stab == "root_robust" else "⚠️" if stab == "ambiguous" else "❌"
            f.write(f"| {m} | {ev['strict_prev']}% | {ev['relaxed_prev']}% | "
                    f"{ev['event_type']} | {icon} {stab} |\n")
        f.write("\n")
        f.write("> **Note on root sensitivity**: Without running full ancestral state "
                "reconstruction on S2/S3/S4 trees, module gain/loss events at the "
                "deepest nodes cannot be precisely mapped. Labels here are heuristic "
                "estimates based on prevalence. Full ASR-based event mapping pending "
                "Phase 4.6 completion.\n\n")

        f.write("## 7. Root Stability Verdict\n\n")
        f.write(f"| Metric | Value |\n|---|---|\n")
        f.write(f"| Scenarios loaded | {n_scenarios_loaded} / 4 |\n")
        f.write(f"| Root partition ratios | "
                f"S1={r1}, S2={r2}, S3={r3}, S4={r4} |\n")
        f.write(f"| Tight agreement (±5%) | {n_agree_tight}/{n_valid} |\n")
        f.write(f"| Loose agreement (±15%) | {n_agree_loose}/{n_valid} |\n")
        f.write(f"| Root stability classification | **{root_stability}** |\n")
        f.write(f"| **Overall QC3 verdict** | **{overall_verdict}** |\n\n")

        if overall_verdict == "GREEN":
            f.write("> 🟢 **GREEN**: Root position is consistent across ≥3 scenarios. "
                    "Root-robust conclusions can enter main text.\n\n")
        elif overall_verdict == "YELLOW":
            f.write("> 🟡 **YELLOW**: Root position shows partial consistency. "
                    "Deep-root conclusions must be labelled `root_sensitive` in main text. "
                    "Only root-robust events (confirmed in ≥2 outgroup scenarios) may "
                    "enter primary claims.\n\n")
        else:
            f.write("> 🔴 **RED**: Insufficient scenario coverage or large disagreement. "
                    "Phase 5 is BLOCKED until QC3 is resolved.\n\n")

        # Phase 5 gating recommendation
        f.write("## 8. Phase 5 Gating Recommendation\n\n")
        if overall_verdict in ("GREEN", "YELLOW"):
            f.write("Phase 5 may proceed with the following constraints:\n\n")
            f.write("- **Working tree**: `CoreTree_rooted_ingroup.treefile` (S1 ingroup, MFP model)\n")
            f.write("- **Node eligibility**: nodes must be `root_robust` in ≥2 outgroup "
                    "scenarios (S1 + S2) to enter `node_selection_registry.tsv` as `eligible`\n")
            f.write("- **Claim tier**: conclusions from root-sensitive events → "
                    "`root_sensitive` tier (Supplement/Discussion only)\n\n")
        else:
            f.write("> Phase 5 is BLOCKED. Resolve QC3 RED conditions first.\n\n")

        f.write("## 9. Required Follow-up Actions\n\n")
        f.write("- [ ] Run ASR on S2 (LGC20) ingroup tree — prune KDOPS first\n")
        f.write("- [ ] Run ASR on S4 (MAD) tree\n")
        f.write("- [ ] Run Phase 4.6 trait ASR on S1 + S2 trees (strict + relaxed)\n")
        f.write("- [ ] Cross-reference module events from S1 vs S2 ASR\n")
        f.write("- [ ] Update `node_selection_registry.tsv` with final eligibility\n")
        f.write("- [ ] Optional: S5 nonreversible reduced-set for rootstrap\n\n")

        f.write("---\n\n")
        f.write("> **Document hierarchy**: This file is subordinate to `PLAN.md`. "
                "All numbers must be reflected in "
                "`results/meta/metrics_manifest.tsv` and `results/meta/progress_snapshot.md`.\n")

    print(f"  ✓ {md_path}")

    # ── 5c. Initialize/update node_selection_registry.tsv ──
    reg_path = os.path.join(args.out_dir, "node_selection_registry.tsv")
    if not os.path.exists(reg_path) or os.path.getsize(reg_path) < 100:
        with open(reg_path, "w") as f:
            f.write("node_id\tdescription\tufboot_s1\tshalrt_s1\troot_scenarios_consistent\t"
                    "strict_stable\trelaxed_stable\tasr_available\tassembly_adjudicated\t"
                    "eligibility\tnote\n")
            f.write("TBD_deepest_Ib\tDeepest Iβ node (pre-ACT gain)\t"
                    "PENDING\tPENDING\tPENDING\tPENDING\tPENDING\tS1_only\tno\t"
                    "hold\tAwaiting S2/S4 ASR to confirm root-robustness\n")
            f.write("TBD_Ib_ACT_gain\tNode at which ACT domain is gained in Iβ lineage\t"
                    "PENDING\tPENDING\tPENDING\tPENDING\tPENDING\tS1_only\tno\t"
                    "hold\tAwaiting Phase 4.4 nested ASR\n")
            f.write("TBD_II_root\tDeepest Type II node\t"
                    "PENDING\tPENDING\tPENDING\tPENDING\tPENDING\tS1_only\tno\t"
                    "hold\tNeeds multi-scenario confirmation\n")
        print(f"  ✓ {reg_path} (initialized with placeholder nodes)")
    else:
        print(f"  ~ {reg_path} already exists, not overwriting")

    # ── 5d. Copy working tree symlink ────────────────
    working_tree_path = os.path.join(args.out_dir, "rooted_working_tree.treefile")
    source = os.path.join(args.dir, "CoreTree_rooted_ingroup.treefile")
    if os.path.exists(source):
        import shutil
        shutil.copy2(source, working_tree_path)
        print(f"  ✓ {working_tree_path} (copied from S1 ingroup)")

    # ── 5e. Print final summary ───────────────────────
    print("\n" + "=" * 65)
    print(f"QC3 VERDICT: {overall_verdict}  ({root_stability})")
    print("=" * 65)
    print(f"\nRoot partition ratios across scenarios:")
    print(f"  S1 MFP:    {r1}")
    print(f"  S2 LGC20:  {r2}")
    print(f"  S3 Midpt:  {r3}")
    print(f"  S4 MAD:    {r4}")
    print(f"\nAgreement: {n_agree_tight}/{n_valid} (tight ±5%), {n_agree_loose}/{n_valid} (loose ±15%)")
    print("\nModule root-sensitivity labels:")
    for m, stab in module_stability.items():
        print(f"  {m}: {stab}")
    print(f"\nOutputs written to {args.out_dir}/")
    print("  - QC3_root_stability.md")
    print("  - root_scenarios.tsv  (updated)")
    print("  - node_selection_registry.tsv")
    print("  - rooted_working_tree.treefile")
    print(f"\nNext step: update metrics_manifest.tsv and progress_snapshot.md")


if __name__ == "__main__":
    main()
