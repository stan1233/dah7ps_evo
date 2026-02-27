#!/usr/bin/env python3
# scripts/define_core_columns.py
"""
Define core columns from FoldMason per-column LDDT scores + gap fraction.
- Input: skeleton_aa.fa (MSA), skeleton_lddt_report.html (msa2lddtreport output)
- Output:
  - results/03_msa_core/core_columns.mask
  - results/03_msa_core/skeleton_core_aa.fa
  - optionally skeleton_core_3di.fa
  - core_columns.tsv + qc_core_definition.md

Implements meta/params.json:
  core_definition.lddt_min = "auto_inflection" or numeric
  core_definition.gap_fraction_max
  core_definition.pad_residues
"""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional


GAP_CHARS = set(["-", "."])


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    seq_chunks: List[str] = []
    with path.open("r", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_chunks)))
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if name is not None:
        records.append((name, "".join(seq_chunks)))
    if not records:
        raise SystemExit(f"[ERROR] No FASTA records read from: {path}")
    return records


def extract_json_array_after_key(text: str, key: str) -> str:
    """
    Extract the JSON array string following a key like '"scores"'.
    Robust bracket matching.
    """
    idx = text.find(key)
    if idx == -1:
        raise ValueError(f"Key not found in report: {key}")
    # find first '[' after key
    start = text.find("[", idx)
    if start == -1:
        raise ValueError(f"Could not find '[' after key {key}")
    depth = 0
    end = None
    for i in range(start, len(text)):
        c = text[i]
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                end = i
                break
    if end is None:
        raise ValueError(f"Could not find matching ']' for array after {key}")
    return text[start : end + 1]


def parse_lddt_from_html(html_path: Path) -> Tuple[List[float], Optional[float]]:
    text = html_path.read_text(errors="replace")
    msa_lddt = None
    m = re.search(r'"msaLDDT"\s*:\s*([0-9]+(?:\.[0-9]+)?)', text)
    if m:
        try:
            msa_lddt = float(m.group(1))
        except ValueError:
            msa_lddt = None

    arr_str = extract_json_array_after_key(text, '"scores"')
    scores_raw = json.loads(arr_str)
    if not isinstance(scores_raw, list):
        raise ValueError("Parsed scores is not a list")
    scores = [float(x) for x in scores_raw]
    return scores, msa_lddt


def knee_threshold(scores: List[float]) -> float:
    """
    Simple 'max distance to diagonal' knee on sorted valid scores.
    scores: per-column, with -1 meaning 'unscored'
    Returns threshold on original score scale.
    """
    valid = [s for s in scores if s >= 0.0]
    if not valid:
        raise SystemExit("[ERROR] No valid (>=0) LDDT scores found.")
    valid_sorted = sorted(valid, reverse=True)

    if len(valid_sorted) < 3:
        return min(valid_sorted)

    s_max = max(valid_sorted)
    s_min = min(valid_sorted)
    if s_max == s_min:
        return s_max

    # normalize to [0,1]
    y = [(s - s_min) / (s_max - s_min) for s in valid_sorted]
    n = len(y)
    best_i = 0
    best_dist = -1.0
    for i, yi in enumerate(y):
        xi = i / (n - 1)
        line_y = 1.0 - xi  # diagonal from (0,1) to (1,0)
        dist = line_y - yi
        if dist > best_dist:
            best_dist = dist
            best_i = i
    return valid_sorted[best_i]


def contiguous_blocks(mask: List[bool]) -> List[Tuple[int, int]]:
    blocks: List[Tuple[int, int]] = []
    in_block = False
    start = 0
    for i, keep in enumerate(mask):
        if keep and not in_block:
            in_block = True
            start = i
        elif (not keep) and in_block:
            in_block = False
            blocks.append((start, i - 1))
    if in_block:
        blocks.append((start, len(mask) - 1))
    return blocks


def apply_block_padding(
    keep: List[bool],
    pad: int,
    scores: List[float],
) -> List[bool]:
    """
    Expand each kept contiguous block by ±pad columns.
    Safety: do NOT turn on columns whose score is -1 (unscored), even if padded.
    """
    if pad <= 0:
        return keep[:]
    L = len(keep)
    out = keep[:]
    for s, e in contiguous_blocks(keep):
        s2 = max(0, s - pad)
        e2 = min(L - 1, e + pad)
        for j in range(s2, e2 + 1):
            if scores[j] >= 0.0:  # exclude unscored columns
                out[j] = True
    return out


def write_fasta(path: Path, records: List[Tuple[str, str]]) -> None:
    with path.open("w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            # wrap for readability
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Define core columns from FoldMason per-column LDDT + gap fraction."
    )
    ap.add_argument("--msa", required=True,
                     help="AA MSA FASTA (e.g. results/03_msa_core/skeleton_aa.fa)")
    ap.add_argument("--lddt-report", required=True,
                     help="FoldMason msa2lddtreport HTML (e.g. skeleton_lddt_report.html)")
    ap.add_argument("--msa-3di", default=None,
                     help="Optional 3Di MSA FASTA to mask with the same columns")
    ap.add_argument("--params", default="meta/params.json",
                     help="meta/params.json")
    ap.add_argument("--outdir", default="results/03_msa_core",
                     help="Output directory")
    ap.add_argument("--prefix", default="skeleton",
                     help="Output prefix for core FASTA files")
    ap.add_argument("--lddt-min", default=None,
                     help="Override core_definition.lddt_min (float), e.g. 0.70")
    args = ap.parse_args()

    msa_path = Path(args.msa)
    html_path = Path(args.lddt_report)
    outdir = Path(args.outdir)

    # --- input validation ---
    if not msa_path.is_file():
        raise SystemExit(f"[ERROR] MSA file not found: {msa_path}")
    if not html_path.is_file():
        raise SystemExit(f"[ERROR] LDDT report not found: {html_path}")
    params_path = Path(args.params)
    if not params_path.is_file():
        raise SystemExit(f"[ERROR] Params file not found: {params_path}")

    outdir.mkdir(parents=True, exist_ok=True)

    params = json.loads(params_path.read_text())
    core_params = params.get("core_definition", {})
    gap_max = float(core_params.get("gap_fraction_max", 0.30))
    pad = int(core_params.get("pad_residues", 20))

    if args.lddt_min is not None:
        lddt_min = float(args.lddt_min)
        lddt_method = "manual_override"
    else:
        lddt_min = core_params.get("lddt_min", "auto_inflection")
        lddt_method = "from_params"

    # --- read MSA ---
    records = read_fasta(msa_path)
    names = [r[0] for r in records]
    seqs = [r[1] for r in records]

    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise SystemExit("[ERROR] MSA is not rectangular: sequences have different alignment lengths.")

    print(f"[INFO] MSA: {len(seqs)} sequences, alignment length {L}")

    # --- parse LDDT scores ---
    scores, msa_lddt = parse_lddt_from_html(html_path)
    if len(scores) != L:
        raise SystemExit(f"[ERROR] LDDT scores length ({len(scores)}) != MSA length ({L}). "
                         f"Make sure the report was generated from this exact MSA.")

    print(f"[INFO] LDDT scores parsed: {len(scores)} columns, msaLDDT={msa_lddt}")

    # --- determine threshold ---
    if isinstance(lddt_min, str) and lddt_min == "auto_inflection":
        threshold = knee_threshold(scores)
        lddt_min_used = threshold
        lddt_min_note = "auto_inflection(knee)"
    elif isinstance(lddt_min, (int, float)):
        lddt_min_used = float(lddt_min)
        lddt_min_note = "numeric_from_params"
    else:
        # allow strings like "0.70"
        try:
            lddt_min_used = float(lddt_min)
            lddt_min_note = "coerced_string"
        except Exception as e:
            raise SystemExit(f"[ERROR] Unrecognized lddt_min in params: {lddt_min!r} ({e})")

    print(f"[INFO] LDDT threshold: {lddt_min_used:.4f} ({lddt_min_note})")

    # --- compute gap fractions & base mask ---
    n = len(seqs)
    gap_fracs: List[float] = []
    keep0: List[bool] = []

    for i in range(L):
        gaps = 0
        for s in seqs:
            if s[i] in GAP_CHARS:
                gaps += 1
        gf = gaps / n
        gap_fracs.append(gf)
        score = scores[i]
        # base core rule
        keep = (score >= 0.0) and (score >= lddt_min_used) and (gf <= gap_max)
        keep0.append(keep)

    base_core_len = sum(1 for x in keep0 if x)
    print(f"[INFO] Base core (before padding): {base_core_len} columns")

    # --- apply padding ---
    keep = apply_block_padding(keep0, pad=pad, scores=scores)

    core_len = sum(1 for x in keep if x)
    valid_scores = [s for s in scores if s >= 0.0]
    n_valid = len(valid_scores)

    print(f"[INFO] Core after padding (±{pad}): {core_len} columns")

    # --- write mask ---
    mask_path = outdir / "core_columns.mask"
    mask_str = "".join("1" if k else "0" for k in keep)
    mask_path.write_text(mask_str + "\n")

    # --- write masked AA fasta ---
    core_records = []
    per_seq_non_gap = []
    for name, seq in zip(names, seqs):
        core_seq = "".join(ch for ch, k in zip(seq, keep) if k)
        core_records.append((name, core_seq))
        nongap = sum(1 for ch in core_seq if ch not in GAP_CHARS)
        per_seq_non_gap.append(nongap)

    aa_out = outdir / f"{args.prefix}_core_aa.fa"
    write_fasta(aa_out, core_records)

    # --- optionally mask 3di alignment with same columns ---
    if args.msa_3di:
        msa3_path = Path(args.msa_3di)
        if not msa3_path.is_file():
            raise SystemExit(f"[ERROR] 3Di MSA file not found: {msa3_path}")
        rec3 = read_fasta(msa3_path)
        seq3 = [r[1] for r in rec3]
        if len(seq3) != n:
            raise SystemExit("[ERROR] 3Di MSA record count != AA MSA record count; cannot apply same mask safely.")
        if any(len(s) != L for s in seq3):
            raise SystemExit("[ERROR] 3Di MSA length != AA MSA length; cannot apply same mask safely.")
        core3 = []
        for (name3, s3) in rec3:
            core_seq3 = "".join(ch for ch, k in zip(s3, keep) if k)
            core3.append((name3, core_seq3))
        out3 = outdir / f"{args.prefix}_core_3di.fa"
        write_fasta(out3, core3)

    # --- write per-column table ---
    tsv_path = outdir / "core_columns.tsv"
    with tsv_path.open("w") as fh:
        fh.write("col_index\tlddt\tgap_fraction\tkeep_base\tkeep_padded\n")
        for i in range(L):
            fh.write(
                f"{i+1}\t{scores[i]:.6f}\t{gap_fracs[i]:.6f}\t{int(keep0[i])}\t{int(keep[i])}\n"
            )

    # --- QC report ---
    qc_path = outdir / "qc_core_definition.md"
    non_gap_min = min(per_seq_non_gap) if per_seq_non_gap else 0
    non_gap_med = sorted(per_seq_non_gap)[len(per_seq_non_gap)//2] if per_seq_non_gap else 0
    non_gap_max = max(per_seq_non_gap) if per_seq_non_gap else 0

    qc_lines = []
    qc_lines.append("# QC — Core column definition\n\n")
    qc_lines.append(f"- MSA input: `{msa_path}`\n")
    qc_lines.append(f"- LDDT report: `{html_path}`\n")
    if msa_lddt is not None:
        qc_lines.append(f"- Reported msaLDDT (includes unscored cols in normalization): **{msa_lddt:.4f}**\n")
    qc_lines.append(f"- Alignment length (L): **{L}**\n")
    qc_lines.append(f"- Sequences (N): **{n}**\n")
    qc_lines.append(f"- Valid scored columns (scores>=0): **{n_valid}**\n")
    qc_lines.append("\n## Parameters\n\n")
    qc_lines.append(f"- gap_fraction_max: **{gap_max}**\n")
    qc_lines.append(f"- pad_residues: **{pad}**\n")
    qc_lines.append(f"- lddt_min source: **{lddt_method}** ({lddt_min_note})\n")
    qc_lines.append(f"- lddt_min used: **{lddt_min_used:.4f}**\n")
    qc_lines.append("\n## Output\n\n")
    qc_lines.append(f"- Base core (before padding): **{base_core_len}** columns\n")
    qc_lines.append(f"- Core length (kept columns, after ±{pad} padding): **{core_len}**\n")
    qc_lines.append(f"- Core mask: `{mask_path}`\n")
    qc_lines.append(f"- Core AA MSA: `{aa_out}`\n")
    if args.msa_3di:
        qc_lines.append(f"- Core 3Di MSA: `{outdir / (args.prefix + '_core_3di.fa')}`\n")
    qc_lines.append("\n## Per-sequence non-gap length in core (sanity)\n\n")
    qc_lines.append(f"- min/median/max non-gap residues: **{non_gap_min} / {non_gap_med} / {non_gap_max}**\n")
    qc_lines.append("\n## Notes\n\n")
    qc_lines.append("- Columns with score = -1 are treated as non-core and are not enabled by padding.\n")
    qc_lines.append("- If core_len is far outside ~400–600, adjust lddt_min or gap_fraction_max in params, then rerun.\n")

    qc_path.write_text("".join(qc_lines))

    print(f"\n[OK] Wrote:")
    print(f"  - {mask_path}")
    print(f"  - {aa_out}")
    if args.msa_3di:
        print(f"  - {outdir / (args.prefix + '_core_3di.fa')}")
    print(f"  - {tsv_path}")
    print(f"  - {qc_path}")
    print(f"[SUMMARY] core_len={core_len}, base_core={base_core_len}, lddt_min_used={lddt_min_used:.4f}, gap_max={gap_max}, pad={pad}")


if __name__ == "__main__":
    main()
