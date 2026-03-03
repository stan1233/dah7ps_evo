#!/usr/bin/env python3
"""
Extract module subsequences from full-length DAH7PS sequences.

Two modes:
  1. --extract-tails   Extract C-terminal tails (residues after raw_env_end)
                       for downstream HMM scanning.
  2. --extract-modules Extract per-module sequences using the presence/absence
                       matrix and coordinate rules.

Phase 3.8 of the DAH7PS V4.1 SOP.

Usage:
    # Mode 1: C-tail extraction for HMM scan
    python scripts/extract_module_seqs.py --extract-tails \\
        --coords results/03_msa_core/core_domain_coords.tsv \\
        --sequences results/02_qc/nr80_Ia.fasta results/02_qc/nr80_Ib.fasta results/02_qc/nr80_II.fasta \\
        --output results/03_msa_modules/c_tails.fasta \\
        --min-tail 10

    # Mode 2: Full module extraction
    python scripts/extract_module_seqs.py --extract-modules \\
        --coords results/03_msa_core/core_domain_coords.tsv \\
        --matrix results/03_msa_modules/module_presence_absence_strict.tsv \\
        --sequences results/02_qc/nr80_Ia.fasta results/02_qc/nr80_Ib.fasta results/02_qc/nr80_II.fasta \\
        --domtbl results/03_msa_modules/module_hits.domtbl \\
        --outdir results/03_msa_modules
"""

import argparse
import os
import sys
import csv


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract module subsequences from full-length DAH7PS sequences."
    )
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--extract-tails", action="store_true",
                       help="Mode 1: extract C-terminal tails for HMM scanning")
    mode.add_argument("--extract-modules", action="store_true",
                       help="Mode 2: extract per-module sequences using presence matrix")

    p.add_argument("--coords", required=True,
                   help="core_domain_coords.tsv from Phase 3.6")
    p.add_argument("--sequences", nargs="+", required=True,
                   help="Full-length FASTA file(s) (nr80_*.fasta)")
    p.add_argument("--output", help="Output FASTA for --extract-tails mode")
    p.add_argument("--min-tail", type=int, default=10,
                   help="Minimum tail length to include (default: 10)")

    # Mode 2 specific
    p.add_argument("--matrix",
                   help="module_presence_absence_*.tsv (required for --extract-modules)")
    p.add_argument("--domtbl",
                   help="module_hits.domtbl from hmmsearch (for ACT/CM coords)")
    p.add_argument("--outdir",
                   help="Output directory for per-module files (required for --extract-modules)")
    return p.parse_args()


def parse_fasta_multi(paths):
    """Parse multiple FASTA files, return dict of seqid -> sequence."""
    seqs = {}
    for path in paths:
        if not os.path.isfile(path):
            print(f"[ERROR] FASTA not found: {path}", file=sys.stderr)
            sys.exit(1)
        sid = None
        chunks = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if sid is not None:
                        seqs[sid] = "".join(chunks)
                    sid = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
        if sid is not None:
            seqs[sid] = "".join(chunks)
    return seqs


def load_coords(path):
    """Load core_domain_coords.tsv, return dict of seqid -> row dict."""
    if not os.path.isfile(path):
        print(f"[ERROR] Coords file not found: {path}", file=sys.stderr)
        sys.exit(1)
    coords = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            coords[row["seq_id"]] = row
    return coords


def load_domtbl(path):
    """Parse hmmsearch domtblout, return dict of seqid -> list of hit dicts."""
    hits = {}
    if path is None or not os.path.isfile(path):
        return hits
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue
            seq_id = fields[0]
            hit = {
                "target": seq_id,
                "query": fields[3],   # HMM name (ACT, CM_1, CM_2)
                "acc": fields[4],     # HMM accession
                "ievalue": float(fields[12]),
                "score": float(fields[13]),
                "ali_from": int(fields[17]),  # alignment start on target
                "ali_to": int(fields[18]),    # alignment end on target
                "env_from": int(fields[19]),
                "env_to": int(fields[20]),
            }
            if seq_id not in hits:
                hits[seq_id] = []
            hits[seq_id].append(hit)
    return hits


def extract_tails(coords, seqs, min_tail, output):
    """Extract C-terminal tails (residues after raw_env_end)."""
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    n_written = 0
    n_short = 0
    with open(output, "w") as fout:
        for seq_id, row in coords.items():
            raw_env_end = int(row["raw_env_end"])
            seq_len = int(row["seq_len"])
            tail_len = seq_len - raw_env_end
            if tail_len < min_tail:
                n_short += 1
                continue
            if seq_id not in seqs:
                continue
            full_seq = seqs[seq_id]
            tail_seq = full_seq[raw_env_end:]  # 0-indexed: raw_env_end onwards
            if len(tail_seq) < min_tail:
                n_short += 1
                continue
            fout.write(f">{seq_id} tail_start={raw_env_end + 1} tail_end={seq_len}\n")
            # Write in 80-char lines
            for i in range(0, len(tail_seq), 80):
                fout.write(tail_seq[i:i+80] + "\n")
            n_written += 1

    print(f"[extract_module_seqs] C-tail extraction: {n_written} sequences written, "
          f"{n_short} skipped (< {min_tail} aa)", file=sys.stderr)
    print(f"[extract_module_seqs] Output: {output}", file=sys.stderr)
    return n_written


def parse_env_segments(seg_str):
    """Parse env_segments like '5-271;274-399' into list of (start, end) tuples."""
    segments = []
    for part in seg_str.split(";"):
        part = part.strip()
        if "-" in part:
            s, e = part.split("-", 1)
            segments.append((int(s), int(e)))
    return segments


def extract_modules(coords, seqs, matrix_path, domtbl_path, outdir):
    """Extract per-module sequences based on presence/absence matrix."""
    os.makedirs(outdir, exist_ok=True)

    # Load presence/absence matrix
    if not os.path.isfile(matrix_path):
        print(f"[ERROR] Matrix not found: {matrix_path}", file=sys.stderr)
        sys.exit(1)

    with open(matrix_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        matrix_cols = reader.fieldnames
        matrix = {row["seq_id"]: row for row in reader}

    # Load HMM hits for ACT/CM coordinates
    hmm_hits = load_domtbl(domtbl_path)

    # Determine which modules are present in the matrix
    module_cols = [c for c in matrix_cols if c not in ("seq_id", "boundary_confidence")]

    module_stats = {}
    for mod in module_cols:
        fasta_path = os.path.join(outdir, f"{mod}_seqs.fasta")
        coords_path = os.path.join(outdir, f"{mod}_domain_coords.tsv")
        n_written = 0

        with open(fasta_path, "w") as fasta_out, \
             open(coords_path, "w") as coords_out:
            coords_out.write("seq_id\tseq_len\tmodule\tfrom\tto\tmodule_len\n")

            for seq_id, row in matrix.items():
                if row.get(mod, "0") != "1":
                    continue
                if seq_id not in seqs:
                    continue
                if seq_id not in coords:
                    continue

                full_seq = seqs[seq_id]
                coord_row = coords[seq_id]
                seq_len = int(coord_row["seq_len"])
                raw_env_start = int(coord_row["raw_env_start"])
                raw_env_end = int(coord_row["raw_env_end"])

                start, end = None, None

                if mod == "N_ext":
                    start = 1
                    end = raw_env_start - 1
                elif mod == "alpha2beta3_insert":
                    # Extract the gap between env_segments
                    segs = parse_env_segments(coord_row["env_segments"])
                    if len(segs) >= 2:
                        # Take the largest internal gap
                        gaps = []
                        for i in range(len(segs) - 1):
                            gap_start = segs[i][1] + 1
                            gap_end = segs[i + 1][0] - 1
                            if gap_end >= gap_start:
                                gaps.append((gap_start, gap_end, gap_end - gap_start + 1))
                        if gaps:
                            # Use the largest gap
                            gaps.sort(key=lambda x: x[2], reverse=True)
                            start, end = gaps[0][0], gaps[0][1]
                elif mod == "C_tail":
                    start = raw_env_end + 1
                    end = seq_len
                elif mod in ("ACT_domain", "CM_domain"):
                    # Use HMM envelope coordinates
                    if seq_id in hmm_hits:
                        relevant = []
                        for h in hmm_hits[seq_id]:
                            if mod == "ACT_domain" and "ACT" in h["query"]:
                                relevant.append(h)
                            elif mod == "CM_domain" and "CM" in h["query"]:
                                relevant.append(h)
                        if relevant:
                            # Take best hit (lowest i-Evalue)
                            best = min(relevant, key=lambda x: x["ievalue"])
                            start = best["env_from"]
                            end = best["env_to"]

                if start is not None and end is not None and end >= start:
                    # Convert to 0-indexed for slicing
                    subseq = full_seq[start - 1:end]
                    if len(subseq) > 0:
                        fasta_out.write(f">{seq_id} {mod}={start}-{end}\n")
                        for i in range(0, len(subseq), 80):
                            fasta_out.write(subseq[i:i+80] + "\n")
                        coords_out.write(f"{seq_id}\t{seq_len}\t{mod}\t{start}\t{end}\t{end - start + 1}\n")
                        n_written += 1

        module_stats[mod] = n_written
        print(f"[extract_module_seqs] {mod}: {n_written} sequences → {fasta_path}",
              file=sys.stderr)

    return module_stats


def main():
    args = parse_args()

    coords = load_coords(args.coords)
    print(f"[extract_module_seqs] Loaded {len(coords)} coordinate records",
          file=sys.stderr)

    seqs = parse_fasta_multi(args.sequences)
    print(f"[extract_module_seqs] Loaded {len(seqs)} full-length sequences",
          file=sys.stderr)

    if args.extract_tails:
        if not args.output:
            print("[ERROR] --output required for --extract-tails mode",
                  file=sys.stderr)
            sys.exit(1)
        extract_tails(coords, seqs, args.min_tail, args.output)

    elif args.extract_modules:
        if not args.matrix:
            print("[ERROR] --matrix required for --extract-modules mode",
                  file=sys.stderr)
            sys.exit(1)
        if not args.outdir:
            print("[ERROR] --outdir required for --extract-modules mode",
                  file=sys.stderr)
            sys.exit(1)
        extract_modules(coords, seqs, args.matrix, args.domtbl, args.outdir)


if __name__ == "__main__":
    main()
