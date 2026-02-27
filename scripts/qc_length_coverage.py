#!/usr/bin/env python3
"""
Phase 2.1: Length + HMM Coverage Triple-Bin Filtering

Classifies sequences into PASS_CANONICAL / PASS_LONG / FRAG based on:
  1. HMM coverage (best single domain or merged multi-domain)
  2. Sequence length within subtype-specific windows

Supports two coverage modes:
  --cov_mode best    : use best single domain hit coverage (default)
  --cov_mode merged  : stitch all qualifying domain hits and compute union
                       coverage (required for Type II due to α2β3 fragmentation)

Usage:
  # Ia/Ib (single domain is sufficient)
  python scripts/qc_length_coverage.py \
    --fasta results/01_mining/hits_Ia_seqs.fasta \
    --domtbl results/01_mining/hits_Ia.domtbl \
    --hmm_length 355 --subtype Ia \
    --canonical_min 300 --canonical_max 450 \
    --cov_min 0.70 --outdir results/02_qc

  # Type II (multi-domain stitching for α2β3 insertion)
  python scripts/qc_length_coverage.py \
    --fasta results/01_mining/hits_II_final_seqs.fasta \
    --domtbl results/01_mining/hits_II.domtbl \
    --hmm_length 471 --subtype II \
    --canonical_min 320 --canonical_max 600 \
    --cov_min 0.70 --cov_mode merged --merge_gap 5 \
    --outdir results/02_qc

Outputs per subtype:
  qc_pass_{subtype}.fasta         PASS_CANONICAL sequences
  qc_long_{subtype}.fasta         PASS_LONG sequences
  fragments_{subtype}.fasta       Fragment sequences
  qc_classification_{subtype}.tsv Per-sequence classification table
"""

import argparse
import os
import sys

# Import from hmmer_utils (same directory)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hmmer_utils import (
    parse_domtbl_all_domains,
    parse_domtbl_best_domain,
    union_coverage_hmm,
)


def load_fasta_records(path):
    """Load FASTA, return ordered list of (seqid, header_line, sequence)."""
    records = []
    current_id = None
    current_header = None
    current_seq = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, current_header, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_header = line.rstrip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id is not None:
            records.append((current_id, current_header, "".join(current_seq)))

    return records


def classify_sequence(l_seq, cov, canonical_min, canonical_max, cov_min):
    """Classify a sequence into one of three bins."""
    if cov < cov_min:
        return "FRAG"
    if l_seq < canonical_min:
        return "FRAG"
    if l_seq > canonical_max:
        return "PASS_LONG"
    return "PASS_CANONICAL"


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2.1: Length + HMM coverage triple-bin filtering"
    )
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--domtbl", required=True, help="domtblout for this subtype")
    parser.add_argument("--hmm_length", type=int, required=True, help="HMM model length (match states)")
    parser.add_argument("--subtype", required=True, help="Subtype label (Ia/Ib/II)")
    parser.add_argument("--canonical_min", type=int, required=True, help="Min length for PASS_CANONICAL")
    parser.add_argument("--canonical_max", type=int, required=True, help="Max length for PASS_CANONICAL")
    parser.add_argument("--cov_min", type=float, default=0.70, help="Min HMM coverage for PASS (default: 0.70)")
    parser.add_argument("--cov_mode", choices=["best", "merged"], default="best",
                        help="Coverage mode: 'best' (single domain) or 'merged' (multi-domain stitching)")
    parser.add_argument("--merge_gap", type=int, default=0,
                        help="Gap tolerance for merging adjacent HMM intervals (default: 0)")
    parser.add_argument("--ievalue_max", type=float, default=1e-5,
                        help="Max i-Evalue for domain hits in merged mode (default: 1e-5)")
    parser.add_argument("--min_hmm_span", type=int, default=30,
                        help="Min HMM span for domain hits in merged mode (default: 30)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    # Validate inputs
    for path in [args.fasta, args.domtbl]:
        if not os.path.isfile(path):
            print(f"ERROR: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # Print config
    print(f"Processing {args.subtype}...")
    print(f"  FASTA: {args.fasta}")
    print(f"  domtbl: {args.domtbl}")
    print(f"  HMM length: {args.hmm_length}")
    print(f"  Canonical window: [{args.canonical_min}, {args.canonical_max}]")
    print(f"  Coverage threshold: {args.cov_min}")
    print(f"  Coverage mode: {args.cov_mode}")
    if args.cov_mode == "merged":
        print(f"  Merge gap: {args.merge_gap}")
        print(f"  i-Evalue max: {args.ievalue_max}")
        print(f"  Min HMM span: {args.min_hmm_span}")
    print()

    # Parse domain hits
    if args.cov_mode == "merged":
        all_domains = parse_domtbl_all_domains(
            args.domtbl,
            ievalue_max=args.ievalue_max,
            min_hmm_span=args.min_hmm_span,
        )
    # Always parse best for comparison
    best_domains = parse_domtbl_best_domain(args.domtbl)

    # Load FASTA records
    records = load_fasta_records(args.fasta)

    # Classify each sequence
    classifications = {}
    bins = {"PASS_CANONICAL": [], "PASS_LONG": [], "FRAG": []}
    rescued_count = 0

    for seqid, header, seq in records:
        l_seq = len(seq)

        # Compute cov_best (single best domain)
        if seqid in best_domains:
            _, hmm_from, hmm_to, _, _ = best_domains[seqid]
            cov_best = (hmm_to - hmm_from + 1) / args.hmm_length
        else:
            cov_best = 0.0

        # Compute cov_merged (multi-domain stitching)
        if args.cov_mode == "merged" and seqid in all_domains:
            cov_merged, merged_intervals, n_doms = union_coverage_hmm(
                all_domains[seqid], args.hmm_length, args.merge_gap
            )
        elif args.cov_mode == "merged":
            cov_merged = 0.0
            n_doms = 0
        else:
            cov_merged = cov_best
            n_doms = 1 if seqid in best_domains else 0

        # Use the appropriate coverage for classification
        if args.cov_mode == "merged":
            cov_used = cov_merged
        else:
            cov_used = cov_best

        bin_label = classify_sequence(l_seq, cov_used, args.canonical_min, args.canonical_max, args.cov_min)

        # Track rescues (would be FRAG by best, but PASS by merged)
        rescued = (cov_best < args.cov_min and cov_merged >= args.cov_min and l_seq >= args.canonical_min)
        if rescued:
            rescued_count += 1

        classifications[seqid] = {
            "bin": bin_label,
            "l_seq": l_seq,
            "cov_best": cov_best,
            "cov_merged": cov_merged,
            "n_domains": n_doms,
            "rescued": rescued,
        }
        bins[bin_label].append((seqid, header, seq))

    # Write output FASTAs
    for bin_label, prefix in [("PASS_CANONICAL", "qc_pass"), ("PASS_LONG", "qc_long"), ("FRAG", "fragments")]:
        out_path = os.path.join(args.outdir, f"{prefix}_{args.subtype}.fasta")
        with open(out_path, "w") as f:
            for seqid, header, seq in bins[bin_label]:
                f.write(f">{header[1:]}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i + 80] + "\n")

    # Write classification TSV
    tsv_path = os.path.join(args.outdir, f"qc_classification_{args.subtype}.tsv")
    with open(tsv_path, "w") as f:
        f.write("seq_id\tL_seq\tcov_best\tcov_merged\tn_domains\trescued_by_stitching\tbin\n")
        for seqid, header, seq in records:
            c = classifications[seqid]
            f.write(f"{seqid}\t{c['l_seq']}\t{c['cov_best']:.4f}\t{c['cov_merged']:.4f}\t"
                    f"{c['n_domains']}\t{c['rescued']}\t{c['bin']}\n")

    # Summary
    total = len(records)
    n_pass = len(bins["PASS_CANONICAL"])
    n_long = len(bins["PASS_LONG"])
    n_frag = len(bins["FRAG"])

    print(f"{'=' * 60}")
    print(f"Phase 2.1 Results — {args.subtype} (cov_mode={args.cov_mode})")
    print(f"{'=' * 60}")
    print(f"  Total input:       {total:>7}")
    print(f"  PASS_CANONICAL:    {n_pass:>7}  ({n_pass / total * 100:.1f}%)")
    print(f"  PASS_LONG:         {n_long:>7}  ({n_long / total * 100:.1f}%)")
    print(f"  FRAG:              {n_frag:>7}  ({n_frag / total * 100:.1f}%)")
    if args.cov_mode == "merged":
        print(f"  Rescued by stitching: {rescued_count:>4}")
    print()

    # Per-bin stats
    for bin_label in ["PASS_CANONICAL", "PASS_LONG", "FRAG"]:
        seqs = bins[bin_label]
        if seqs:
            lens = [len(s[2]) for s in seqs]
            covs_best = [classifications[s[0]]["cov_best"] for s in seqs]
            covs_merged = [classifications[s[0]]["cov_merged"] for s in seqs]
            print(f"  {bin_label}:")
            print(f"    Length: min={min(lens)}, median={sorted(lens)[len(lens) // 2]}, max={max(lens)}")
            print(f"    cov_best:   min={min(covs_best):.3f}, median={sorted(covs_best)[len(covs_best) // 2]:.3f}, max={max(covs_best):.3f}")
            if args.cov_mode == "merged":
                print(f"    cov_merged: min={min(covs_merged):.3f}, median={sorted(covs_merged)[len(covs_merged) // 2]:.3f}, max={max(covs_merged):.3f}")
        else:
            print(f"  {bin_label}: (empty)")
    print()

    # Sanity check
    assert n_pass + n_long + n_frag == total, \
        f"FATAL: bins don't sum to total! {n_pass}+{n_long}+{n_frag} != {total}"
    print(f"  ✅ Sanity check passed: {n_pass}+{n_long}+{n_frag} = {total}")
    print()
    print(f"Outputs:")
    print(f"  {args.outdir}/qc_pass_{args.subtype}.fasta ({n_pass})")
    print(f"  {args.outdir}/qc_long_{args.subtype}.fasta ({n_long})")
    print(f"  {args.outdir}/fragments_{args.subtype}.fasta ({n_frag})")
    print(f"  {args.outdir}/qc_classification_{args.subtype}.tsv")


if __name__ == "__main__":
    main()
