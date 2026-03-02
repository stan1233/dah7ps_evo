#!/usr/bin/env python3
"""
Extract core-domain subsequences from full-length DAH7PS sequences using
an HMM profile, with mandatory hit stitching (SOP CHECK-06) and boundary
padding (SOP CHECK-01).

Phase 3.6 of the DAH7PS V4.1 SOP.

Pipeline:
  1. Run hmmsearch --domtblout against the core HMM
  2. Parse domain hits, filter by i-Evalue and HMM span
  3. Stitch multi-domain hits in HMM coordinate space (merge_gap tolerance)
  4. Compute coverage = stitched_hmm_span / model_length
  5. Filter by coverage >= threshold
  6. Extract padded envelope region from each sequence (CHECK-01: --pad N)
  7. Output FASTA + TSV (with per-segment env coordinates for traceability)

V4.1-fix changes vs original:
  - [Fix B / CHECK-01] Added --pad N: extends extraction by N residues beyond
    the outermost envelope boundaries, giving hmmalign context to recover
    terminal match states. Extra residues become insert states in hmmalign
    and are stripped by esl-alimask.
  - [Fix A] TSV now includes 'env_segments' column with per-segment
    (env_from-env_to) coordinates, enabling downstream linker extraction
    (Phase 3.9) and coordinate_mapper.py (Phase 6.3) to trace multi-segment
    hits without re-parsing domtblout.
  - Extraction changed from "concatenate segments" to "single continuous
    region (padded_start..padded_end)". For multi-segment sequences, the
    internal insertion (e.g. α2β3) is now included; hmmalign handles it
    as insert states → stripped by esl-alimask. This is safer than the
    previous concat approach which could confuse hmmalign at junction points.

Reuses hmmer_utils.py for domtblout parsing and interval merging.
"""

import argparse
import json
import os
import subprocess
import sys

# Add scripts directory to path for hmmer_utils
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hmmer_utils import parse_domtbl_all_domains, merge_intervals, union_coverage_hmm


def parse_fasta(path):
    """Parse FASTA file, return dict of seqid -> sequence."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(current_seq)
    return seqs


def get_hmm_length(hmm_path):
    """Extract model length from HMM file header (LENG line)."""
    with open(hmm_path) as f:
        for line in f:
            if line.startswith("LENG"):
                return int(line.split()[1])
    raise ValueError(f"Could not find LENG in {hmm_path}")


def run_hmmsearch(hmm_path, fasta_path, domtblout_path, cpu=4):
    """Run hmmsearch and produce domtblout."""
    cmd = [
        "hmmsearch",
        "--domtblout", domtblout_path,
        "--noali",
        "--cpu", str(cpu),
        hmm_path,
        fasta_path,
    ]
    print(f"[extract_core_domains] Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] hmmsearch failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    print(f"[extract_core_domains] hmmsearch completed.", file=sys.stderr)


def format_env_segments(env_merged):
    """Format merged envelope segments as semicolon-separated ranges.

    Example: [(5, 200), (250, 399)] -> "5-200;250-399"
    """
    return ";".join(f"{s}-{e}" for s, e in env_merged)


def extract_core_domains(
    hmm_path,
    fasta_path,
    out_fasta,
    out_tsv,
    ievalue=1e-5,
    hmm_span_min=30,
    merge_gap=5,
    coverage_min=0.70,
    pad=0,
    cpu=4,
):
    """Main extraction pipeline."""
    # 0. Get model length
    hmm_len = get_hmm_length(hmm_path)
    print(f"[extract_core_domains] HMM model length: {hmm_len}", file=sys.stderr)
    print(f"[extract_core_domains] Boundary padding: {pad} aa", file=sys.stderr)

    # 1. Run hmmsearch
    domtbl_path = out_tsv.replace(".tsv", "_domtblout.txt")
    run_hmmsearch(hmm_path, fasta_path, domtbl_path, cpu=cpu)

    # 2. Parse domain hits
    all_domains = parse_domtbl_all_domains(
        domtbl_path, ievalue_max=ievalue, min_hmm_span=hmm_span_min
    )
    print(
        f"[extract_core_domains] Sequences with qualifying hits: {len(all_domains)}",
        file=sys.stderr,
    )

    # 3. Parse input FASTA
    seqs = parse_fasta(fasta_path)
    print(f"[extract_core_domains] Input sequences: {len(seqs)}", file=sys.stderr)

    # 4. Stitch, filter, extract
    n_pass = 0
    n_fail_cov = 0
    n_no_hit = 0
    n_stitched = 0

    tsv_rows = []
    fasta_records = []

    for seqid, seq in seqs.items():
        seq_len = len(seq)

        if seqid not in all_domains:
            n_no_hit += 1
            continue

        domains = all_domains[seqid]
        n_hits = len(domains)

        # Compute stitched HMM coverage (NOT affected by padding)
        coverage, hmm_merged, _ = union_coverage_hmm(
            domains, hmm_len, merge_gap=merge_gap
        )

        stitched = n_hits > 1 and len(hmm_merged) < n_hits
        if stitched:
            n_stitched += 1

        # Filter by coverage
        if coverage < coverage_min:
            n_fail_cov += 1
            continue

        # Compute envelope segments (before padding, for traceability)
        env_intervals = [(d["env_from"], d["env_to"]) for d in domains]
        env_merged = merge_intervals(env_intervals, merge_gap=0)

        # Raw envelope boundaries (before padding)
        raw_env_start = env_merged[0][0]
        raw_env_end = env_merged[-1][1]

        # [CHECK-01] Apply boundary padding: extend extraction beyond envelope
        # Padding gives hmmalign context to recover terminal match states.
        # Extra residues become insert states → stripped by esl-alimask.
        # Padding is applied to outer boundaries only, not between segments.
        padded_start = max(1, raw_env_start - pad)
        padded_end = min(seq_len, raw_env_end + pad)

        # Extract as single continuous region (padded_start..padded_end).
        # For multi-segment sequences, the internal gap (e.g. α2β3 insert)
        # is now included; hmmalign will assign it to insert states.
        # This is safer than concatenating segments which could confuse
        # hmmalign at splice junctions.
        core_seq = seq[padded_start - 1 : padded_end]

        # Actual padding applied (may be less than requested at seq boundaries)
        pad_left = raw_env_start - padded_start
        pad_right = padded_end - raw_env_end

        # HMM coverage range
        hmm_from = hmm_merged[0][0]
        hmm_to = hmm_merged[-1][1]

        n_pass += 1
        fasta_records.append((seqid, core_seq))
        tsv_rows.append(
            {
                "seq_id": seqid,
                "seq_len": seq_len,
                "core_start": padded_start,
                "core_end": padded_end,
                "raw_env_start": raw_env_start,
                "raw_env_end": raw_env_end,
                "pad_left": pad_left,
                "pad_right": pad_right,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "coverage": f"{coverage:.4f}",
                "n_hits": n_hits,
                "n_env_segments": len(env_merged),
                "env_segments": format_env_segments(env_merged),
                "stitched_flag": "Y" if stitched else "N",
                "core_residues": len(core_seq),
            }
        )

    # 5. Print summary
    print(f"\n[extract_core_domains] === Summary ===", file=sys.stderr)
    print(f"  Input sequences:    {len(seqs)}", file=sys.stderr)
    print(f"  No qualifying hit:  {n_no_hit}", file=sys.stderr)
    print(f"  Failed coverage:    {n_fail_cov} (< {coverage_min})", file=sys.stderr)
    print(f"  Passed (output):    {n_pass}", file=sys.stderr)
    print(f"  Stitched (>1 hit merged): {n_stitched}", file=sys.stderr)
    print(f"  Boundary padding:   {pad} aa (CHECK-01)", file=sys.stderr)

    # 6. Write outputs
    os.makedirs(os.path.dirname(out_fasta) or ".", exist_ok=True)
    with open(out_fasta, "w") as f:
        for seqid, seq in fasta_records:
            f.write(f">{seqid}\n{seq}\n")
    print(f"[extract_core_domains] Wrote {out_fasta} ({n_pass} seqs)", file=sys.stderr)

    tsv_header = [
        "seq_id",
        "seq_len",
        "core_start",
        "core_end",
        "raw_env_start",
        "raw_env_end",
        "pad_left",
        "pad_right",
        "hmm_from",
        "hmm_to",
        "coverage",
        "n_hits",
        "n_env_segments",
        "env_segments",
        "stitched_flag",
        "core_residues",
    ]
    with open(out_tsv, "w") as f:
        f.write("\t".join(tsv_header) + "\n")
        for row in tsv_rows:
            f.write("\t".join(str(row[k]) for k in tsv_header) + "\n")
    print(
        f"[extract_core_domains] Wrote {out_tsv} ({n_pass} rows)", file=sys.stderr
    )

    return n_pass


def main():
    parser = argparse.ArgumentParser(
        description="Extract core-domain subsequences using HMM profile with hit stitching (Phase 3.6)."
    )
    parser.add_argument(
        "--hmm", required=True, help="Path to core HMM profile (core_global.hmm)"
    )
    parser.add_argument(
        "--fasta", required=True, help="Path to input FASTA (nr80_all.fasta)"
    )
    parser.add_argument(
        "--params",
        default=None,
        help="Path to params.json (reads qc.hmm_coverage_min and core_definition.pad_residues)",
    )
    parser.add_argument(
        "--out_fasta",
        required=True,
        help="Output FASTA of core-domain subsequences",
    )
    parser.add_argument(
        "--out_tsv", required=True, help="Output TSV with domain coordinates"
    )
    parser.add_argument(
        "--ievalue",
        type=float,
        default=1e-5,
        help="Max individual domain E-value (default: 1e-5)",
    )
    parser.add_argument(
        "--hmm_span_min",
        type=int,
        default=30,
        help="Min HMM span per hit (default: 30)",
    )
    parser.add_argument(
        "--merge_gap",
        type=int,
        default=5,
        help="Max HMM gap for stitching (default: 5)",
    )
    parser.add_argument(
        "--coverage_min",
        type=float,
        default=None,
        help="Min HMM coverage to retain (default: from params or 0.70)",
    )
    parser.add_argument(
        "--pad",
        type=int,
        default=None,
        help="Boundary padding in residues (CHECK-01; default: from params or 20)",
    )
    parser.add_argument(
        "--cpu", type=int, default=4, help="CPUs for hmmsearch (default: 4)"
    )

    args = parser.parse_args()

    # Validate inputs
    for p in [args.hmm, args.fasta]:
        if not os.path.isfile(p):
            print(f"[ERROR] File not found: {p}", file=sys.stderr)
            sys.exit(1)

    # Read defaults from params.json
    coverage_min = 0.70
    pad = 20  # CHECK-01 default
    if args.params and os.path.isfile(args.params):
        with open(args.params) as f:
            params = json.load(f)
        coverage_min = params.get("qc", {}).get("hmm_coverage_min", 0.70)
        pad = params.get("core_definition", {}).get("pad_residues", 20)
        print(
            f"[extract_core_domains] From params.json: coverage_min={coverage_min}, pad={pad}",
            file=sys.stderr,
        )

    # CLI overrides
    if args.coverage_min is not None:
        coverage_min = args.coverage_min
        print(
            f"[extract_core_domains] coverage_min override: {coverage_min}",
            file=sys.stderr,
        )
    if args.pad is not None:
        pad = args.pad
        print(
            f"[extract_core_domains] pad override: {pad}",
            file=sys.stderr,
        )

    n_pass = extract_core_domains(
        hmm_path=args.hmm,
        fasta_path=args.fasta,
        out_fasta=args.out_fasta,
        out_tsv=args.out_tsv,
        ievalue=args.ievalue,
        hmm_span_min=args.hmm_span_min,
        merge_gap=args.merge_gap,
        coverage_min=coverage_min,
        pad=pad,
        cpu=args.cpu,
    )

    if n_pass == 0:
        print("[ERROR] No sequences passed filters!", file=sys.stderr)
        sys.exit(1)

    print(f"\n[extract_core_domains] Done. {n_pass} sequences extracted.", file=sys.stderr)


if __name__ == "__main__":
    main()
