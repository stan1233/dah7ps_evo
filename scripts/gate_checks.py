#!/usr/bin/env python3
"""
Pre-Phase 2 Gate Checks

Gate A: Cross-subtype overlap in hits IDs
Gate B: Isolate borderline KDOPS sequences (delta 0-20)

Usage:
  python scripts/gate_checks.py --workdir /home/tynan/0218
"""

import argparse
import os
import sys


def load_ids_from_file(path):
    """Load IDs from a text file (one per line)."""
    ids = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                ids.add(line)
    return ids


def load_ids_from_fasta(path):
    """Load sequence IDs from a FASTA file."""
    ids = set()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                seqid = line[1:].split()[0]
                ids.add(seqid)
    return ids


def gate_a(workdir):
    """Check cross-subtype overlap."""
    print("=" * 60)
    print("Gate A: Cross-subtype hits overlap")
    print("=" * 60)

    ia_ids = load_ids_from_file(os.path.join(workdir, "results/01_mining/hits_Ia_ids.txt"))
    # For Ib, use clean FASTA IDs (post-KDOPS)
    ib_ids = load_ids_from_fasta(os.path.join(workdir, "results/01_mining/hits_Ib_clean.fasta"))
    ii_ids = load_ids_from_file(os.path.join(workdir, "results/01_mining/hits_II_ids.txt"))

    print(f"  Ia:       {len(ia_ids):>6}")
    print(f"  Ib clean: {len(ib_ids):>6}")
    print(f"  II:       {len(ii_ids):>6}")
    print()

    ia_ib = ia_ids & ib_ids
    ia_ii = ia_ids & ii_ids
    ib_ii = ib_ids & ii_ids
    triple = ia_ids & ib_ids & ii_ids

    print(f"  Ia ∩ Ib:       {len(ia_ib):>6}  ({len(ia_ib) / min(len(ia_ids), len(ib_ids)) * 100:.2f}% of smaller set)")
    print(f"  Ia ∩ II:       {len(ia_ii):>6}  ({len(ia_ii) / min(len(ia_ids), len(ii_ids)) * 100:.2f}% of smaller set)")
    print(f"  Ib ∩ II:       {len(ib_ii):>6}  ({len(ib_ii) / min(len(ib_ids), len(ii_ids)) * 100:.2f}% of smaller set)")
    print(f"  Ia ∩ Ib ∩ II:  {len(triple):>6}")
    print()

    # Decision
    max_pct = max(
        len(ia_ib) / min(len(ia_ids), len(ib_ids)),
        len(ia_ii) / min(len(ia_ids), len(ii_ids)),
        len(ib_ii) / min(len(ib_ids), len(ii_ids)),
    )

    if max_pct > 0.05:
        print("  ⚠ OVERLAP > 5% detected — best-hit subtype assignment recommended")
        print("  Writing overlap details...")
        # Save overlap IDs for potential best-hit assignment
        overlap_dir = os.path.join(workdir, "results/01_mining")
        with open(os.path.join(overlap_dir, "overlap_Ia_Ib.txt"), "w") as f:
            for sid in sorted(ia_ib):
                f.write(sid + "\n")
        with open(os.path.join(overlap_dir, "overlap_Ia_II.txt"), "w") as f:
            for sid in sorted(ia_ii):
                f.write(sid + "\n")
        with open(os.path.join(overlap_dir, "overlap_Ib_II.txt"), "w") as f:
            for sid in sorted(ib_ii):
                f.write(sid + "\n")
        return True  # needs best-hit assignment
    else:
        print("  ✅ All overlaps < 5% — no deduplication needed")
        return False


def gate_b(workdir):
    """Isolate borderline KDOPS sequences."""
    print()
    print("=" * 60)
    print("Gate B: Ib borderline KDOPS isolation")
    print("=" * 60)

    report_path = os.path.join(workdir, "results/01_mining/kdops_filter_report_v41.tsv")
    if not os.path.isfile(report_path):
        print(f"  ERROR: Report file not found: {report_path}")
        return []

    borderline = []
    with open(report_path) as f:
        header = f.readline()  # skip header
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 5:
                continue
            seq_id = cols[0]
            verdict = cols[4]
            try:
                delta = float(cols[3])
            except ValueError:
                continue
            # Borderline: kept (KEEP) but delta is small (0 <= delta <= 20)
            if verdict == "KEEP" and 0 <= delta <= 20:
                borderline.append((seq_id, delta))

    output_path = os.path.join(workdir, "results/01_mining/kdops_borderline_ids.txt")
    with open(output_path, "w") as f:
        for seq_id, delta in sorted(borderline):
            f.write(f"{seq_id}\n")

    print(f"  Found {len(borderline)} borderline sequences (delta 0–20, KEEP):")
    for seq_id, delta in sorted(borderline, key=lambda x: x[1]):
        print(f"    {seq_id}  delta={delta:.1f}")
    print(f"  → Written to {output_path}")

    return [b[0] for b in borderline]


def main():
    parser = argparse.ArgumentParser(description="Pre-Phase 2 gate checks")
    parser.add_argument("--workdir", default="/home/tynan/0218", help="Project root")
    args = parser.parse_args()

    needs_bha = gate_a(args.workdir)
    borderline_ids = gate_b(args.workdir)

    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Gate A: {'⚠ needs best-hit assignment' if needs_bha else '✅ pass (overlap < 5%)'}")
    print(f"  Gate B: {len(borderline_ids)} borderline IDs isolated")
    print()

    if needs_bha:
        sys.exit(1)  # signal that Phase 2 needs extra step
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
