#!/usr/bin/env python3
"""
Extract a structural subset from an MSA, with ID remapping for FoldMason.

Reads an aligned FASTA and a panel manifest, extracts sequences whose IDs
appear in the manifest, remaps headers from UniRef90_<ACC> to
AF-<ACC>-F1-model_v4 (matching panelDb.lookup keys), and strips any
non-FASTA lines (e.g. seqkit [INFO] pollution).

Phase 3.8 engineering debt — formalises the one-off fix from Phase 3.7.

Usage:
    python scripts/extract_struct_subset.py \
        --input  results/03_msa_core/core_tree.afa \
        --manifest results/03_msa_core/panel_manifest.tsv \
        --output results/03_msa_core/core_tree_struct_subset.foldmason.fa \
        [--id_col rep_id] [--no-remap]
"""

import argparse
import os
import re
import sys


ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
UNIREF_RE = re.compile(r"UniRef\d+_([A-Za-z0-9]+)")


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract structural subset from MSA with FoldMason-compatible ID remapping."
    )
    p.add_argument("--input", required=True, help="Input aligned FASTA")
    p.add_argument("--manifest", required=True,
                   help="Panel manifest TSV (must have header row)")
    p.add_argument("--output", required=True, help="Output FASTA")
    p.add_argument("--id_col", default="rep_id",
                   help="Column name in manifest for sequence IDs (default: rep_id)")
    p.add_argument("--no-remap", action="store_true",
                   help="Skip UniRef→AF header remapping (keep original IDs)")
    return p.parse_args()


def load_panel_ids(manifest_path, id_col):
    """Read panel IDs from manifest TSV."""
    ids = set()
    with open(manifest_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if id_col not in header:
            print(f"[ERROR] Column '{id_col}' not found in manifest. "
                  f"Available: {header}", file=sys.stderr)
            sys.exit(1)
        col_idx = header.index(id_col)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if col_idx < len(fields) and fields[col_idx].strip():
                ids.add(fields[col_idx].strip())
    return ids


def is_noise_line(line):
    """Detect seqkit/conda log lines that pollute FASTA output."""
    cleaned = ANSI_RE.sub("", line).strip()
    if not cleaned:
        return True
    if cleaned.startswith("[INFO]"):
        return True
    if "patterns loaded from file" in cleaned:
        return True
    if cleaned.startswith("# >>>") or cleaned.startswith("ERROR conda"):
        return True
    return False


def remap_header(sid):
    """UniRef90_<ACC> → AF-<ACC>-F1-model_v4."""
    m = UNIREF_RE.match(sid)
    if m:
        return f"AF-{m.group(1)}-F1-model_v4"
    return sid


def main():
    args = parse_args()

    if not os.path.isfile(args.input):
        print(f"[ERROR] Input not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.manifest):
        print(f"[ERROR] Manifest not found: {args.manifest}", file=sys.stderr)
        sys.exit(1)

    panel_ids = load_panel_ids(args.manifest, args.id_col)
    print(f"[extract_struct_subset] Loaded {len(panel_ids)} panel IDs from "
          f"{args.manifest} (col={args.id_col})", file=sys.stderr)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    n_found = 0
    writing = False
    with open(args.input) as fin, open(args.output, "w") as fout:
        for raw_line in fin:
            line = ANSI_RE.sub("", raw_line.rstrip("\n"))

            if is_noise_line(line):
                continue

            if line.startswith(">"):
                sid = line[1:].split()[0]
                writing = sid in panel_ids
                if writing:
                    n_found += 1
                    out_id = sid if args.no_remap else remap_header(sid)
                    fout.write(f">{out_id}\n")
            elif writing:
                # Strip spaces (shouldn't be in aligned FASTA)
                fout.write(line.replace(" ", "") + "\n")

    print(f"[extract_struct_subset] Extracted {n_found}/{len(panel_ids)} "
          f"sequences → {args.output}", file=sys.stderr)

    missing = panel_ids - set()  # We don't track found IDs, just count
    if n_found < len(panel_ids):
        print(f"[WARNING] {len(panel_ids) - n_found} panel IDs not found "
              f"in input MSA", file=sys.stderr)

    if n_found == 0:
        print("[ERROR] No sequences extracted!", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
