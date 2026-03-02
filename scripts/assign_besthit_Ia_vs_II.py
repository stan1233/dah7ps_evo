#!/usr/bin/env python3
"""
Best-hit subtype assignment for Ia/II overlapping sequences.

Compares full-sequence bit scores from hmmsearch domtblout to assign each
overlapping sequence to its best-matching subtype, with 3-tier confidence.

Usage:
  python scripts/assign_besthit_Ia_vs_II.py \
    --ia_domtbl results/01_mining/hits_Ia.domtbl \
    --ii_domtbl results/01_mining/hits_II.domtbl \
    --overlap_ids results/01_mining/overlap_Ia_II.txt \
    --ia_all_ids results/01_mining/hits_Ia_ids.txt \
    --ii_all_ids results/01_mining/hits_II_ids.txt \
    --outdir results/01_mining

Outputs:
  besthit_Ia_vs_II.tsv         Per-sequence assignment table
  hits_Ia_final_ids.txt        Mutually exclusive Ia IDs
  hits_II_final_ids.txt        Mutually exclusive II IDs
  hits_IaII_lowconf_ids.txt    LOW confidence IDs (exclude from seeds)
  qc_subtype_assignment.md     QC1b report
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_domtbl_scores(path):
    """Parse domtblout, return dict: seqid -> best full-sequence score (col 7, 0-indexed)."""
    best = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 22:
                continue
            seqid = cols[0]
            score = float(cols[7])  # full-sequence bit score
            if seqid not in best or score > best[seqid]:
                best[seqid] = score
    return best


def load_ids(path):
    """Load IDs from text file."""
    ids = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                ids.add(line)
    return ids


def classify_confidence(delta_bits):
    """Classify confidence based on delta bits."""
    abs_delta = abs(delta_bits)
    if abs_delta >= 20:
        return "HIGH"
    elif abs_delta >= 10:
        return "MED"
    else:
        return "LOW"


def main():
    parser = argparse.ArgumentParser(
        description="Best-hit subtype assignment for Ia/II overlap"
    )
    parser.add_argument("--ia_domtbl", required=True, help="Ia domtblout")
    parser.add_argument("--ii_domtbl", required=True, help="II domtblout")
    parser.add_argument("--overlap_ids", required=True, help="Overlap IDs file")
    parser.add_argument("--ia_all_ids", required=True, help="All Ia hit IDs")
    parser.add_argument("--ii_all_ids", required=True, help="All II hit IDs")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    # Validate inputs
    for path in [args.ia_domtbl, args.ii_domtbl, args.overlap_ids, args.ia_all_ids, args.ii_all_ids]:
        if not os.path.isfile(path):
            print(f"ERROR: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # Parse scores
    print("Parsing domtblout files...")
    ia_scores = parse_domtbl_scores(args.ia_domtbl)
    ii_scores = parse_domtbl_scores(args.ii_domtbl)

    # Load IDs
    overlap_ids = load_ids(args.overlap_ids)
    ia_all = load_ids(args.ia_all_ids)
    ii_all = load_ids(args.ii_all_ids)

    print(f"  Ia scores: {len(ia_scores)} sequences")
    print(f"  II scores: {len(ii_scores)} sequences")
    print(f"  Overlap:   {len(overlap_ids)} sequences")

    # --- Assignment ---
    assignments = []
    conf_counts = defaultdict(int)
    subtype_counts = defaultdict(int)
    low_conf_ids = set()

    for seqid in sorted(overlap_ids):
        score_ia = ia_scores.get(seqid, float("-inf"))
        score_ii = ii_scores.get(seqid, float("-inf"))
        delta = score_ia - score_ii  # positive = Ia wins

        if delta >= 0:
            assigned = "Ia"
        else:
            assigned = "II"

        conf = classify_confidence(delta)
        conf_counts[conf] += 1
        subtype_counts[assigned] += 1

        if conf == "LOW":
            low_conf_ids.add(seqid)

        assignments.append({
            "seq_id": seqid,
            "score_Ia": score_ia,
            "score_II": score_ii,
            "delta_bits": delta,
            "assigned_subtype": assigned,
            "confidence": conf,
        })

    # --- Build final mutually exclusive ID sets ---
    # Non-overlap Ia-only and II-only sequences stay in their original subtype
    ia_only = ia_all - overlap_ids
    ii_only = ii_all - overlap_ids

    # From overlap, assign to winner
    overlap_to_ia = {a["seq_id"] for a in assignments if a["assigned_subtype"] == "Ia"}
    overlap_to_ii = {a["seq_id"] for a in assignments if a["assigned_subtype"] == "II"}

    ia_final = ia_only | overlap_to_ia
    ii_final = ii_only | overlap_to_ii

    # Sanity check: mutual exclusivity
    assert len(ia_final & ii_final) == 0, "FATAL: final ID sets are not mutually exclusive!"

    # --- Write outputs ---

    # 1. Assignment table
    tsv_path = os.path.join(args.outdir, "besthit_Ia_vs_II.tsv")
    with open(tsv_path, "w") as f:
        f.write("seq_id\tscore_Ia\tscore_II\tdelta_bits\tassigned_subtype\tconfidence\n")
        for a in assignments:
            f.write(f"{a['seq_id']}\t{a['score_Ia']:.1f}\t{a['score_II']:.1f}\t"
                    f"{a['delta_bits']:.1f}\t{a['assigned_subtype']}\t{a['confidence']}\n")

    # 2. Final ID lists
    ia_final_path = os.path.join(args.outdir, "hits_Ia_final_ids.txt")
    with open(ia_final_path, "w") as f:
        for sid in sorted(ia_final):
            f.write(sid + "\n")

    ii_final_path = os.path.join(args.outdir, "hits_II_final_ids.txt")
    with open(ii_final_path, "w") as f:
        for sid in sorted(ii_final):
            f.write(sid + "\n")

    # 3. Low-confidence IDs
    low_path = os.path.join(args.outdir, "hits_IaII_lowconf_ids.txt")
    with open(low_path, "w") as f:
        for sid in sorted(low_conf_ids):
            f.write(sid + "\n")

    # --- Summary ---
    print()
    print("=" * 60)
    print("Best-hit Ia vs II Assignment Results")
    print("=" * 60)
    print()
    print(f"Overlap sequences:      {len(overlap_ids)}")
    print(f"  → Assigned to Ia:     {subtype_counts['Ia']}")
    print(f"  → Assigned to II:     {subtype_counts['II']}")
    print()
    print(f"Confidence breakdown:")
    print(f"  HIGH (|Δ| ≥ 20):     {conf_counts['HIGH']}  ({conf_counts['HIGH']/len(overlap_ids)*100:.1f}%)")
    print(f"  MED  (10 ≤ |Δ| < 20): {conf_counts['MED']}  ({conf_counts['MED']/len(overlap_ids)*100:.1f}%)")
    print(f"  LOW  (|Δ| < 10):     {conf_counts['LOW']}  ({conf_counts['LOW']/len(overlap_ids)*100:.1f}%)")
    print()
    print(f"Final mutually exclusive sets:")
    print(f"  Ia_final: {len(ia_final):>6}  (Ia-only={len(ia_only)} + overlap→Ia={len(overlap_to_ia)})")
    print(f"  II_final: {len(ii_final):>6}  (II-only={len(ii_only)} + overlap→II={len(overlap_to_ii)})")
    print(f"  Ia∩II:    {len(ia_final & ii_final):>6}  (should be 0)")
    print()
    print(f"LOW confidence IDs (exclude from seeds): {len(low_conf_ids)}")
    print()

    # --- Evaluation of go/no-go ---
    low_pct = conf_counts["LOW"] / len(overlap_ids) * 100 if overlap_ids else 0
    print("=" * 60)
    print("Go/No-Go Evaluation")
    print("=" * 60)

    warnings = []
    if len(ia_final) < 1500:
        warnings.append(f"  ⚠ Ia_final={len(ia_final)} is very small (< 1,500)")
    if low_pct > 40:
        warnings.append(f"  ⚠ LOW confidence = {low_pct:.1f}% (> 40%) — HMMs may lack discriminating power")

    if warnings:
        print("  🚨 POTENTIAL ISSUES DETECTED:")
        for w in warnings:
            print(w)
        print()
        print("  Consider fallback strategies:")
        print("    1. Merge Ia+II into joint Phase 2; split by tree clade in Phase 4")
        print("    2. Build module-specific HMMs (α2β3 / β-hairpin) for subtype discrimination")
    else:
        print("  ✅ GO: Ia_final and II_final are both substantial")
        print(f"     LOW confidence < 40% ({low_pct:.1f}%), assignment is reliable")
    print()

    # --- Write QC1b report ---
    qc_path = os.path.join(args.outdir, "qc_subtype_assignment.md")
    with open(qc_path, "w") as f:
        f.write("# QC1b: Ia vs II Best-Hit Subtype Assignment Report\n\n")
        f.write(f"> Generated: 2026-02-25\n")
        f.write(f"> Input: {args.ia_domtbl} vs {args.ii_domtbl}\n\n")
        f.write("---\n\n")
        f.write("## 1. Overlap Summary\n\n")
        f.write(f"| Metric | Count |\n")
        f.write(f"|--------|-------|\n")
        f.write(f"| Ia total hits | {len(ia_all)} |\n")
        f.write(f"| II total hits | {len(ii_all)} |\n")
        f.write(f"| Ia ∩ II overlap | {len(overlap_ids)} ({len(overlap_ids)/len(ia_all)*100:.1f}% of Ia) |\n")
        f.write(f"| Ia-only (no overlap) | {len(ia_only)} |\n")
        f.write(f"| II-only (no overlap) | {len(ii_only)} |\n\n")

        f.write("## 2. Assignment Results\n\n")
        f.write(f"| Assigned to | Count | % of overlap |\n")
        f.write(f"|-------------|-------|--------------|\n")
        f.write(f"| Ia | {subtype_counts['Ia']} | {subtype_counts['Ia']/len(overlap_ids)*100:.1f}% |\n")
        f.write(f"| II | {subtype_counts['II']} | {subtype_counts['II']/len(overlap_ids)*100:.1f}% |\n\n")

        f.write("## 3. Confidence Distribution\n\n")
        f.write(f"| Confidence | Criterion | Count | % |\n")
        f.write(f"|------------|-----------|-------|---|\n")
        f.write(f"| HIGH | |Δ bits| ≥ 20 | {conf_counts['HIGH']} | {conf_counts['HIGH']/len(overlap_ids)*100:.1f}% |\n")
        f.write(f"| MED | 10 ≤ |Δ bits| < 20 | {conf_counts['MED']} | {conf_counts['MED']/len(overlap_ids)*100:.1f}% |\n")
        f.write(f"| LOW | |Δ bits| < 10 | {conf_counts['LOW']} | {conf_counts['LOW']/len(overlap_ids)*100:.1f}% |\n\n")

        f.write("## 4. Final Mutually Exclusive Sets\n\n")
        f.write(f"| Set | Count | Composition |\n")
        f.write(f"|-----|-------|-------------|\n")
        f.write(f"| Ia_final | {len(ia_final)} | Ia-only({len(ia_only)}) + overlap→Ia({len(overlap_to_ia)}) |\n")
        f.write(f"| II_final | {len(ii_final)} | II-only({len(ii_only)}) + overlap→II({len(overlap_to_ii)}) |\n")
        f.write(f"| Ia∩II (should=0) | {len(ia_final & ii_final)} | — |\n\n")

        f.write("## 5. Go/No-Go Evaluation\n\n")
        if warnings:
            f.write("**🚨 Issues detected:**\n\n")
            for w in warnings:
                f.write(f"- {w.strip()}\n")
            f.write("\n**Consider fallback strategies** (see log.md for details)\n")
        else:
            f.write(f"**✅ GO:** Ia_final ({len(ia_final)}) and II_final ({len(ii_final)}) are both substantial. ")
            f.write(f"LOW confidence = {low_pct:.1f}% (< 40%). Assignment is reliable.\n")

        f.write("\n## 6. Output Files\n\n")
        f.write("| File | Description |\n")
        f.write("|------|-------------|\n")
        f.write(f"| `besthit_Ia_vs_II.tsv` | Per-sequence assignment table ({len(overlap_ids)} rows) |\n")
        f.write(f"| `hits_Ia_final_ids.txt` | Mutually exclusive Ia IDs ({len(ia_final)}) |\n")
        f.write(f"| `hits_II_final_ids.txt` | Mutually exclusive II IDs ({len(ii_final)}) |\n")
        f.write(f"| `hits_IaII_lowconf_ids.txt` | LOW confidence IDs ({len(low_conf_ids)}) |\n")

    print(f"Outputs written to {args.outdir}/")
    print(f"  besthit_Ia_vs_II.tsv")
    print(f"  hits_Ia_final_ids.txt ({len(ia_final)} IDs)")
    print(f"  hits_II_final_ids.txt ({len(ii_final)} IDs)")
    print(f"  hits_IaII_lowconf_ids.txt ({len(low_conf_ids)} IDs)")
    print(f"  qc_subtype_assignment.md")


if __name__ == "__main__":
    main()
