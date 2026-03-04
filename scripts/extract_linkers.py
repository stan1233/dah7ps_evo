#!/usr/bin/env python3
"""
Extract linker (C-flank) regions from full-length sequences for Phase 3.9.

For each sequence in seq_ids, the linker is defined as the region OUTSIDE the
core envelope (raw_env_end+1 .. seq_len).  The ACT domain is an insert state
embedded within the core HMM envelope and is already handled by the separate
ACT_domain_msa.afa; it is NOT part of the linker.

Optionally also writes the N-flank (1 .. raw_env_start-1) to a separate FASTA
for inspection, but the primary output used for E-INS-i alignment is the C-flank.

Phase 3.9 of the DAH7PS V5.0 SOP.

Usage:
    python scripts/extract_linkers.py \
        --full_length_fasta results/02_qc/nr80_Ib.fasta \
        --seq_ids           results/03_msa_full/Ib_ACT.ids \
        --core_coords       results/03_msa_core/core_domain_coords.tsv \
        --output            results/03_msa_full/Ib_ACT_linkers.fasta \
        [--output_nflank    results/03_msa_full/Ib_ACT_nflanks.fasta]
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract C-flank linker regions for Phase 3.9 profile-anchored stitching."
    )
    p.add_argument("--full_length_fasta", required=True,
                   help="Full-length FASTA for the relevant subtype (e.g. nr80_Ib.fasta)")
    p.add_argument("--seq_ids", required=True,
                   help="Text file of sequence IDs to process (one per line)")
    p.add_argument("--core_coords", required=True,
                   help="core_domain_coords.tsv from Phase 3.6")
    p.add_argument("--output", required=True,
                   help="Output FASTA: C-flank (raw_env_end+1 .. seq_len) per sequence")
    p.add_argument("--output_nflank",
                   help="Optional output FASTA: N-flank (1 .. raw_env_start-1) per sequence")
    p.add_argument("--min_cflank", type=int, default=1,
                   help="Minimum C-flank length to include (default: 1). "
                        "Sequences with shorter C-flanks get an empty/gap placeholder.")
    return p.parse_args()


def read_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                ids.append(line)
    return ids


def parse_fasta(path, keep_ids=None):
    """Parse FASTA → dict id -> sequence (uppercase). keep_ids: set of IDs to load."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    if keep_ids is None or current_id in keep_ids:
                        seqs[current_id] = "".join(current_seq).upper()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        if keep_ids is None or current_id in keep_ids:
            seqs[current_id] = "".join(current_seq).upper()
    return seqs


def load_core_coords(path, keep_ids):
    """Load core_domain_coords.tsv → dict id -> {raw_env_start, raw_env_end, seq_len}."""
    coords = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {col: i for i, col in enumerate(header)}
        required = {"seq_id", "raw_env_start", "raw_env_end", "seq_len"}
        missing = required - set(idx)
        if missing:
            print(f"[ERROR] core_coords missing columns: {missing}", file=sys.stderr)
            sys.exit(1)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            sid = fields[idx["seq_id"]]
            if sid not in keep_ids:
                continue
            coords[sid] = {
                "raw_env_start": int(fields[idx["raw_env_start"]]),
                "raw_env_end": int(fields[idx["raw_env_end"]]),
                "seq_len": int(fields[idx["seq_len"]]),
            }
    return coords


def write_fasta(records, path):
    """Write list of (id, seq) to FASTA."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        for sid, seq in records:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def main():
    args = parse_args()

    for path in [args.full_length_fasta, args.seq_ids, args.core_coords]:
        if not os.path.isfile(path):
            print(f"[ERROR] Input not found: {path}", file=sys.stderr)
            sys.exit(1)

    seq_ids = read_ids(args.seq_ids)
    keep_set = set(seq_ids)
    print(f"[extract_linkers] Processing {len(seq_ids)} sequences", file=sys.stderr)

    # Load full-length sequences
    fl_seqs = parse_fasta(args.full_length_fasta, keep_ids=keep_set)
    missing_seqs = keep_set - set(fl_seqs)
    if missing_seqs:
        print(f"[ERROR] {len(missing_seqs)} seq IDs not found in {args.full_length_fasta}: "
              f"{sorted(missing_seqs)[:5]}...", file=sys.stderr)
        sys.exit(1)

    # Load core coordinates
    coords = load_core_coords(args.core_coords, keep_set)
    missing_coords = keep_set - set(coords)
    if missing_coords:
        print(f"[ERROR] {len(missing_coords)} seq IDs not found in core_coords: "
              f"{sorted(missing_coords)[:5]}...", file=sys.stderr)
        sys.exit(1)

    cflank_records = []
    nflank_records = []
    stats = {"cflank_len": [], "nflank_len": []}

    for sid in seq_ids:
        seq = fl_seqs[sid]
        c = coords[sid]
        raw_start = c["raw_env_start"]  # 1-based
        raw_end = c["raw_env_end"]      # 1-based
        seq_len = len(seq)

        # Sanity check
        if seq_len != c["seq_len"]:
            print(f"[WARNING] {sid}: FASTA length {seq_len} != coords seq_len {c['seq_len']}; "
                  f"using FASTA length", file=sys.stderr)

        # C-flank: raw_env_end+1 to seq_len (0-based: raw_end .. seq_len)
        cflank = seq[raw_end:]   # raw_end is 1-based end → 0-based index = raw_end
        stats["cflank_len"].append(len(cflank))

        if len(cflank) >= args.min_cflank:
            cflank_records.append((sid, cflank))
        else:
            print(f"[WARNING] {sid}: C-flank too short ({len(cflank)} aa < min {args.min_cflank}); "
                  f"skipping", file=sys.stderr)

        # N-flank: 1 to raw_env_start-1 (0-based: 0 .. raw_start-1)
        nflank = seq[:raw_start - 1]
        stats["nflank_len"].append(len(nflank))
        if args.output_nflank:
            nflank_records.append((sid, nflank))

    # Summary stats
    cf = sorted(stats["cflank_len"])
    nf = sorted(stats["nflank_len"])
    n = len(cf)
    print(f"[extract_linkers] C-flank stats: n={n} "
          f"min={cf[0]} median={cf[n//2]} max={cf[-1]}", file=sys.stderr)
    print(f"[extract_linkers] N-flank stats: n={n} "
          f"min={nf[0]} median={nf[n//2]} max={nf[-1]}", file=sys.stderr)

    if len(cflank_records) == 0:
        print("[ERROR] No C-flank sequences extracted. Check core_coords raw_env_end values.",
              file=sys.stderr)
        sys.exit(1)

    write_fasta(cflank_records, args.output)
    print(f"[extract_linkers] Wrote {len(cflank_records)} C-flank sequences → {args.output}",
          file=sys.stderr)

    if args.output_nflank and nflank_records:
        write_fasta(nflank_records, args.output_nflank)
        print(f"[extract_linkers] Wrote {len(nflank_records)} N-flank sequences → {args.output_nflank}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
