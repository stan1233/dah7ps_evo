#!/usr/bin/env python3
"""
Phase 2.1: Length + HMM Coverage Triple-Bin Filtering

Classifies sequences into PASS_CANONICAL / PASS_LONG / FRAG based on:
  1. HMM coverage: cov_hmm = (hmm_to - hmm_from + 1) / hmm_length
  2. Sequence length: L_seq within subtype-specific windows

Usage:
  python scripts/qc_length_coverage.py \
    --fasta results/01_mining/hits_Ia_seqs.fasta \
    --domtbl results/01_mining/hits_Ia.domtbl \
    --hmm_length 355 \
    --subtype Ia \
    --canonical_min 300 --canonical_max 450 \
    --cov_min 0.70 \
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


def parse_domtbl_best_domain(path):
    """Parse domtblout, return dict: seqid -> (best_score, hmm_from, hmm_to, env_from, env_to).

    For each sequence, selects the domain hit with the highest individual domain score
    (column 13, 0-indexed). Also records HMM coordinates for coverage calculation.
    """
    best = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 22:
                continue
            seqid = cols[0]
            dom_score = float(cols[13])  # domain score (this domain)
            hmm_from = int(cols[15])     # hmm coord from
            hmm_to = int(cols[16])       # hmm coord to
            env_from = int(cols[19])     # envelope from (on sequence)
            env_to = int(cols[20])       # envelope to (on sequence)

            if seqid not in best or dom_score > best[seqid][0]:
                best[seqid] = (dom_score, hmm_from, hmm_to, env_from, env_to)
    return best


def load_fasta_lengths(path):
    """Load FASTA, return dict: seqid -> length, and ordered list of (seqid, header, seq)."""
    lengths = {}
    records = []
    current_id = None
    current_header = None
    current_seq = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    seq = "".join(current_seq)
                    lengths[current_id] = len(seq)
                    records.append((current_id, current_header, seq))
                current_id = line[1:].split()[0]
                current_header = line.rstrip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id is not None:
            seq = "".join(current_seq)
            lengths[current_id] = len(seq)
            records.append((current_id, current_header, seq))

    return lengths, records


def classify_sequence(l_seq, cov_hmm, canonical_min, canonical_max, cov_min):
    """Classify a sequence into one of three bins."""
    if cov_hmm < cov_min:
        return "FRAG"
    if l_seq < canonical_min:
        # Short but high coverage — could be a truncated but valid core
        # Still classify as FRAG if below canonical minimum
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
    parser.add_argument("--canonical_max", type=int, required=True, help="Max length for PASS_CANONICAL (above = PASS_LONG)")
    parser.add_argument("--cov_min", type=float, default=0.70, help="Min HMM coverage for PASS (default: 0.70)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    # Validate inputs
    for path in [args.fasta, args.domtbl]:
        if not os.path.isfile(path):
            print(f"ERROR: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # Parse data
    print(f"Processing {args.subtype}...")
    print(f"  FASTA: {args.fasta}")
    print(f"  domtbl: {args.domtbl}")
    print(f"  HMM length: {args.hmm_length}")
    print(f"  Canonical window: [{args.canonical_min}, {args.canonical_max}]")
    print(f"  Coverage threshold: {args.cov_min}")
    print()

    dom_info = parse_domtbl_best_domain(args.domtbl)
    lengths, records = load_fasta_lengths(args.fasta)

    # Classify
    classifications = {}  # seqid -> (bin, l_seq, cov_hmm)
    bins = {"PASS_CANONICAL": [], "PASS_LONG": [], "FRAG": []}

    for seqid, header, seq in records:
        l_seq = len(seq)

        if seqid in dom_info:
            _, hmm_from, hmm_to, _, _ = dom_info[seqid]
            cov_hmm = (hmm_to - hmm_from + 1) / args.hmm_length
        else:
            # No domain hit in domtbl — treat as fragment
            cov_hmm = 0.0

        bin_label = classify_sequence(l_seq, cov_hmm, args.canonical_min, args.canonical_max, args.cov_min)
        classifications[seqid] = (bin_label, l_seq, cov_hmm)
        bins[bin_label].append((seqid, header, seq))

    # Write output FASTAs
    for bin_label, prefix in [("PASS_CANONICAL", "qc_pass"), ("PASS_LONG", "qc_long"), ("FRAG", "fragments")]:
        out_path = os.path.join(args.outdir, f"{prefix}_{args.subtype}.fasta")
        with open(out_path, "w") as f:
            for seqid, header, seq in bins[bin_label]:
                f.write(f">{header[1:]}\n")
                # Write sequence in 80-char lines
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")

    # Write classification TSV
    tsv_path = os.path.join(args.outdir, f"qc_classification_{args.subtype}.tsv")
    with open(tsv_path, "w") as f:
        f.write("seq_id\tL_seq\tcov_hmm\tbin\n")
        for seqid, header, seq in records:
            bin_label, l_seq, cov_hmm = classifications[seqid]
            f.write(f"{seqid}\t{l_seq}\t{cov_hmm:.4f}\t{bin_label}\n")

    # Summary
    total = len(records)
    n_pass = len(bins["PASS_CANONICAL"])
    n_long = len(bins["PASS_LONG"])
    n_frag = len(bins["FRAG"])

    print(f"{'='*60}")
    print(f"Phase 2.1 Results — {args.subtype}")
    print(f"{'='*60}")
    print(f"  Total input:       {total:>7}")
    print(f"  PASS_CANONICAL:    {n_pass:>7}  ({n_pass/total*100:.1f}%)")
    print(f"  PASS_LONG:         {n_long:>7}  ({n_long/total*100:.1f}%)")
    print(f"  FRAG:              {n_frag:>7}  ({n_frag/total*100:.1f}%)")
    print()

    # Quick stats per bin
    for bin_label in ["PASS_CANONICAL", "PASS_LONG", "FRAG"]:
        seqs = bins[bin_label]
        if seqs:
            lens = [len(s[2]) for s in seqs]
            covs = [classifications[s[0]][2] for s in seqs]
            print(f"  {bin_label}:")
            print(f"    Length: min={min(lens)}, median={sorted(lens)[len(lens)//2]}, max={max(lens)}")
            print(f"    HMM cov: min={min(covs):.3f}, median={sorted(covs)[len(covs)//2]:.3f}, max={max(covs):.3f}")
        else:
            print(f"  {bin_label}: (empty)")
    print()

    # Sanity check
    assert n_pass + n_long + n_frag == total, f"FATAL: bins don't sum to total! {n_pass}+{n_long}+{n_frag} != {total}"
    print(f"  ✅ Sanity check passed: {n_pass}+{n_long}+{n_frag} = {total}")
    print()
    print(f"Outputs:")
    print(f"  {args.outdir}/qc_pass_{args.subtype}.fasta ({n_pass})")
    print(f"  {args.outdir}/qc_long_{args.subtype}.fasta ({n_long})")
    print(f"  {args.outdir}/fragments_{args.subtype}.fasta ({n_frag})")
    print(f"  {args.outdir}/qc_classification_{args.subtype}.tsv")


if __name__ == "__main__":
    main()
