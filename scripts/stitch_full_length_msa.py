#!/usr/bin/env python3
"""
Stitch core + module + linker alignments into a full-length MSA.

Column order in output:
    [core columns (fixed, from core_asr.afa)]
    [module columns (fixed, from module_msa.afa)]
    [linker/C-flank columns (free, from linker_einsi.afa)]

Hard assertion:
    The core segment of the output exactly matches the corresponding rows in
    core_msa (i.e. core columns are never re-aligned).

Phase 3.9 of the DAH7PS V5.0 SOP.

Usage:
    python scripts/stitch_full_length_msa.py \
        --seq_ids        results/03_msa_full/Ib_ACT.ids \
        --core_msa       results/03_msa_core/core_asr.afa \
        --module_msa     results/03_msa_modules/ACT_domain_msa.afa \
        --module_name    ACT \
        --linker_msa     results/03_msa_full/Ib_ACT_linkers_einsi.afa \
        --output         results/03_msa_full/msa_full_Ib_v4.afa \
        --emit_column_map results/03_msa_full/msa_full_Ib_column_map.tsv \
        --assert_core_columns_unchanged
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Stitch core + module + linker MSAs into a full-length alignment."
    )
    p.add_argument("--seq_ids", required=True,
                   help="Text file: one seq_id per line (defines the subset and order)")
    p.add_argument("--core_msa", required=True,
                   help="Core MSA (e.g. core_asr.afa). Rows for seq_ids are extracted.")
    p.add_argument("--module_msa", required=True,
                   help="Module MSA (e.g. ACT_domain_msa.afa). Must contain all seq_ids.")
    p.add_argument("--module_name", required=True,
                   help="Module label used in column map (e.g. ACT)")
    p.add_argument("--linker_msa", required=True,
                   help="Aligned linker/C-flank MSA (from E-INS-i). Must contain all seq_ids.")
    p.add_argument("--output", required=True,
                   help="Output full-length MSA (AFA format)")
    p.add_argument("--emit_column_map", required=True,
                   help="Output TSV: col_start, col_end, segment_type, source_file")
    p.add_argument("--assert_core_columns_unchanged", action="store_true",
                   help="Verify that core columns in output exactly match core_msa (hard assertion).")
    return p.parse_args()


def read_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                ids.append(line)
    return ids


def parse_msa(path, keep_ids=None):
    """Parse aligned FASTA → dict id -> aligned_sequence (gaps preserved).

    keep_ids: if provided, only load those IDs (set).
    All sequences must have the same length; validated after loading.
    """
    seqs = {}
    order = []
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    if keep_ids is None or current_id in keep_ids:
                        seqs[current_id] = "".join(current_seq).upper()
                        order.append(current_id)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        if keep_ids is None or current_id in keep_ids:
            seqs[current_id] = "".join(current_seq).upper()
            order.append(current_id)

    # Validate alignment width
    lengths = set(len(s) for s in seqs.values())
    if len(lengths) > 1:
        print(f"[ERROR] {path}: sequences have unequal lengths: {lengths}", file=sys.stderr)
        sys.exit(1)

    width = lengths.pop() if lengths else 0
    return seqs, order, width


def write_fasta(records, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        for sid, seq in records:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def write_column_map(segments, path):
    """segments: list of (label, source, col_start, col_end) (1-based, inclusive)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        f.write("segment\tsource\tcol_start\tcol_end\tn_cols\n")
        for label, source, cs, ce in segments:
            f.write(f"{label}\t{source}\t{cs}\t{ce}\t{ce - cs + 1}\n")
    print(f"[stitch_full_length_msa] Column map → {path}", file=sys.stderr)


def main():
    args = parse_args()

    for path in [args.seq_ids, args.core_msa, args.module_msa, args.linker_msa]:
        if not os.path.isfile(path):
            print(f"[ERROR] Input not found: {path}", file=sys.stderr)
            sys.exit(1)

    seq_ids = read_ids(args.seq_ids)
    keep_set = set(seq_ids)
    print(f"[stitch_full_length_msa] Stitching {len(seq_ids)} sequences", file=sys.stderr)

    # Load MSAs
    core_seqs, _, core_width = parse_msa(args.core_msa, keep_ids=keep_set)
    mod_seqs, _, mod_width   = parse_msa(args.module_msa, keep_ids=keep_set)
    lnk_seqs, _, lnk_width   = parse_msa(args.linker_msa, keep_ids=keep_set)

    print(f"[stitch_full_length_msa] core={core_width} cols, "
          f"module({args.module_name})={mod_width} cols, "
          f"linker={lnk_width} cols", file=sys.stderr)

    # Check all seq_ids are present in every MSA
    for label, msa_dict in [("core", core_seqs), ("module", mod_seqs), ("linker", lnk_seqs)]:
        missing = keep_set - set(msa_dict)
        if missing:
            print(f"[ERROR] {len(missing)} seq IDs missing from {label} MSA: "
                  f"{sorted(missing)[:5]}...", file=sys.stderr)
            sys.exit(1)

    # Stitch: core | module | linker
    stitched = []
    for sid in seq_ids:
        stitched.append((sid, core_seqs[sid] + mod_seqs[sid] + lnk_seqs[sid]))

    total_cols = core_width + mod_width + lnk_width
    print(f"[stitch_full_length_msa] Output: {len(stitched)} seqs × {total_cols} cols",
          file=sys.stderr)

    # Hard assertion: core columns unchanged
    if args.assert_core_columns_unchanged:
        n_fail = 0
        for sid in seq_ids:
            stitched_core = stitched[seq_ids.index(sid)][1][:core_width]
            if stitched_core != core_seqs[sid]:
                print(f"[ASSERT FAIL] {sid}: core columns differ from core_msa", file=sys.stderr)
                n_fail += 1
        if n_fail > 0:
            print(f"[ERROR] Core column assertion failed for {n_fail} sequences. "
                  f"Aborting.", file=sys.stderr)
            sys.exit(1)
        print(f"[stitch_full_length_msa] ✅ ASSERT PASS: core columns unchanged "
              f"for all {len(seq_ids)} sequences ({core_width} cols)", file=sys.stderr)

    # Write output MSA
    write_fasta(stitched, args.output)
    print(f"[stitch_full_length_msa] Wrote stitched MSA → {args.output}", file=sys.stderr)

    # Write column map (1-based, inclusive)
    segments = [
        ("core",            os.path.basename(args.core_msa),   1,             core_width),
        (f"module:{args.module_name}", os.path.basename(args.module_msa),
         core_width + 1,   core_width + mod_width),
        ("linker:C-flank",  os.path.basename(args.linker_msa),
         core_width + mod_width + 1, total_cols),
    ]
    write_column_map(segments, args.emit_column_map)


if __name__ == "__main__":
    main()
