#!/usr/bin/env python3
"""
Annotate DAH7PS module presence/absence per sequence.

Combines coordinate-based evidence (N_ext, alpha2beta3_insert, C_tail) with
HMM-based evidence (ACT_domain, CM_domain) to produce strict and relaxed
presence/absence matrices.

Phase 3.8 of the DAH7PS V4.1 SOP.

Usage:
    python scripts/annotate_modules.py \\
        --coords results/03_msa_core/core_domain_coords.tsv \\
        --domtbl results/03_msa_modules/module_hits.domtbl \\
        --outdir results/03_msa_modules \\
        [--params meta/params.json]
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict


# ── Default thresholds (overridable via params.json) ─────────────────────────
DEFAULT_THRESHOLDS = {
    "N_ext": {"strict": 25, "relaxed": 10},
    "alpha2beta3_insert": {"strict": 5, "relaxed": 1},
    "ACT_domain": {"strict": 1e-5, "relaxed": 1e-3},
    "CM_domain": {"strict": 1e-5, "relaxed": 1e-3},
    "C_tail": {"strict": 25, "relaxed": 10},
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate DAH7PS module presence/absence."
    )
    p.add_argument("--coords", required=True,
                   help="core_domain_coords.tsv from Phase 3.6")
    p.add_argument("--domtbl",
                   help="module_hits.domtbl from hmmsearch (ACT/CM HMM hits)")
    p.add_argument("--outdir", required=True,
                   help="Output directory for matrices and report")
    p.add_argument("--params",
                   help="params.json (optional, for threshold overrides)")
    return p.parse_args()


def load_coords(path):
    """Load core_domain_coords.tsv → list of row dicts (preserving order)."""
    if not os.path.isfile(path):
        print(f"[ERROR] Coords not found: {path}", file=sys.stderr)
        sys.exit(1)
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def load_domtbl(path):
    """Parse hmmsearch domtblout → dict of seq_id → list of hit dicts."""
    hits = defaultdict(list)
    if path is None or not os.path.isfile(path):
        return hits
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue
            hit = {
                "target": fields[0],
                "query": fields[3],
                "acc": fields[4],
                "ievalue": float(fields[12]),
                "score": float(fields[13]),
                "env_from": int(fields[19]),
                "env_to": int(fields[20]),
            }
            hits[fields[0]].append(hit)
    return hits


def load_thresholds(params_path):
    """Load thresholds from params.json, falling back to defaults."""
    thresholds = dict(DEFAULT_THRESHOLDS)
    if params_path and os.path.isfile(params_path):
        with open(params_path) as f:
            params = json.load(f)
        if "module_annotation" in params:
            ma = params["module_annotation"]
            for mod, vals in ma.items():
                if mod in thresholds:
                    thresholds[mod].update(vals)
    return thresholds


def parse_env_segments(seg_str):
    """Parse '5-271;274-399' → [(5,271), (274,399)]."""
    segments = []
    for part in seg_str.split(";"):
        part = part.strip()
        if "-" in part:
            s, e = part.split("-", 1)
            segments.append((int(s), int(e)))
    return segments


def get_largest_insert_gap(env_segments_str):
    """Get the largest gap between env_segments (α2β3 insert)."""
    segs = parse_env_segments(env_segments_str)
    if len(segs) < 2:
        return 0
    max_gap = 0
    for i in range(len(segs) - 1):
        gap = segs[i + 1][0] - segs[i][1] - 1
        if gap > max_gap:
            max_gap = gap
    return max_gap


def best_hmm_hit(hits_list, module_name):
    """Get best HMM hit for a module (lowest i-Evalue)."""
    relevant = []
    for h in hits_list:
        if module_name == "ACT_domain" and "ACT" in h["query"]:
            relevant.append(h)
        elif module_name == "CM_domain" and "CM" in h["query"]:
            relevant.append(h)
    if relevant:
        return min(relevant, key=lambda x: x["ievalue"])
    return None


def annotate_all(coord_rows, hmm_hits, thresholds):
    """Annotate all sequences, return list of annotation dicts."""
    results = []
    modules = ["N_ext", "alpha2beta3_insert", "ACT_domain", "CM_domain", "C_tail"]

    for row in coord_rows:
        seq_id = row["seq_id"]
        seq_len = int(row["seq_len"])
        raw_env_start = int(row["raw_env_start"])
        raw_env_end = int(row["raw_env_end"])
        n_ext_len = raw_env_start - 1
        c_tail_len = seq_len - raw_env_end
        insert_gap = get_largest_insert_gap(row["env_segments"])

        # Get HMM hits for this sequence
        seq_hits = hmm_hits.get(seq_id, [])
        act_hit = best_hmm_hit(seq_hits, "ACT_domain")
        cm_hit = best_hmm_hit(seq_hits, "CM_domain")

        strict = {}
        relaxed = {}

        # N_ext
        strict["N_ext"] = 1 if n_ext_len >= thresholds["N_ext"]["strict"] else 0
        relaxed["N_ext"] = 1 if n_ext_len >= thresholds["N_ext"]["relaxed"] else 0

        # alpha2beta3_insert
        strict["alpha2beta3_insert"] = 1 if insert_gap >= thresholds["alpha2beta3_insert"]["strict"] else 0
        relaxed["alpha2beta3_insert"] = 1 if insert_gap >= thresholds["alpha2beta3_insert"]["relaxed"] else 0

        # ACT_domain
        if act_hit and act_hit["ievalue"] <= thresholds["ACT_domain"]["strict"]:
            strict["ACT_domain"] = 1
        else:
            strict["ACT_domain"] = 0

        if act_hit and act_hit["ievalue"] <= thresholds["ACT_domain"]["relaxed"]:
            relaxed["ACT_domain"] = 1
        else:
            relaxed["ACT_domain"] = 0

        # CM_domain
        if cm_hit and cm_hit["ievalue"] <= thresholds["CM_domain"]["strict"]:
            strict["CM_domain"] = 1
        else:
            strict["CM_domain"] = 0

        if cm_hit and cm_hit["ievalue"] <= thresholds["CM_domain"]["relaxed"]:
            relaxed["CM_domain"] = 1
        else:
            relaxed["CM_domain"] = 0

        # C_tail: ≥ threshold AND no HMM hit (i.e. not ACT/CM)
        # Use relaxed-level HMM exclusion for BOTH matrices to preserve strict ⊆ relaxed
        has_hmm_hit_relaxed = (relaxed["ACT_domain"] == 1 or relaxed["CM_domain"] == 1)
        strict["C_tail"] = 1 if (c_tail_len >= thresholds["C_tail"]["strict"] and not has_hmm_hit_relaxed) else 0
        relaxed["C_tail"] = 1 if (c_tail_len >= thresholds["C_tail"]["relaxed"] and not has_hmm_hit_relaxed) else 0

        # Boundary confidence: count how many evidence sources agree
        n_evidence = 0
        if n_ext_len > 0:
            n_evidence += 1
        if insert_gap > 0:
            n_evidence += 1
        if act_hit is not None:
            n_evidence += 1
        if cm_hit is not None:
            n_evidence += 1
        if c_tail_len > 0:
            n_evidence += 1
        # Confidence: high if strict == relaxed for all modules
        mismatches = sum(1 for m in modules if strict[m] != relaxed[m])
        if mismatches == 0:
            confidence = "high"
        elif mismatches <= 2:
            confidence = "medium"
        else:
            confidence = "low"

        results.append({
            "seq_id": seq_id,
            "strict": strict,
            "relaxed": relaxed,
            "confidence": confidence,
            "raw_values": {
                "n_ext_len": n_ext_len,
                "insert_gap": insert_gap,
                "c_tail_len": c_tail_len,
                "act_ievalue": act_hit["ievalue"] if act_hit else None,
                "cm_ievalue": cm_hit["ievalue"] if cm_hit else None,
            }
        })

    return results, modules


def write_matrix(results, modules, path, version):
    """Write presence/absence matrix TSV."""
    with open(path, "w") as f:
        f.write("seq_id\t" + "\t".join(modules) + "\tboundary_confidence\n")
        for r in results:
            vals = r[version]
            row = [r["seq_id"]] + [str(vals[m]) for m in modules] + [r["confidence"]]
            f.write("\t".join(row) + "\n")
    print(f"[annotate_modules] Wrote {version} matrix: {path} ({len(results)} rows)",
          file=sys.stderr)


def write_robustness_report(results, modules, thresholds, outdir):
    """Write boundary_robustness.md report."""
    path = os.path.join(outdir, "boundary_robustness.md")

    # Compute summary statistics
    total = len(results)
    strict_counts = {m: 0 for m in modules}
    relaxed_counts = {m: 0 for m in modules}
    delta_counts = {m: 0 for m in modules}

    confidence_dist = {"high": 0, "medium": 0, "low": 0}

    for r in results:
        for m in modules:
            if r["strict"][m] == 1:
                strict_counts[m] += 1
            if r["relaxed"][m] == 1:
                relaxed_counts[m] += 1
            if r["relaxed"][m] == 1 and r["strict"][m] == 0:
                delta_counts[m] += 1
        confidence_dist[r["confidence"]] += 1

    # Raw value distributions for coordinate modules
    n_ext_vals = [r["raw_values"]["n_ext_len"] for r in results]
    insert_vals = [r["raw_values"]["insert_gap"] for r in results]
    c_tail_vals = [r["raw_values"]["c_tail_len"] for r in results]

    def percentiles(vals):
        s = sorted(vals)
        n = len(s)
        if n == 0:
            return "n=0"
        return (f"n={n} min={s[0]} q1={s[n//4]} median={s[n//2]} "
                f"q3={s[3*n//4]} max={s[-1]}")

    with open(path, "w") as f:
        f.write("# Module Boundary Robustness Report\n\n")
        f.write(f"**Date**: 2026-03-03  \n")
        f.write(f"**Total sequences**: {total}  \n")
        f.write(f"**Source**: `core_domain_coords.tsv` + `module_hits.domtbl`\n\n")
        f.write("---\n\n")

        f.write("## 1. Module Prevalence Summary\n\n")
        f.write("| Module | Strict (n) | Strict (%) | Relaxed (n) | Relaxed (%) | Δ (relaxed only) |\n")
        f.write("|--------|-----------|-----------|------------|------------|------------------|\n")
        for m in modules:
            sp = strict_counts[m] / total * 100
            rp = relaxed_counts[m] / total * 100
            f.write(f"| {m} | {strict_counts[m]} | {sp:.1f}% | {relaxed_counts[m]} "
                    f"| {rp:.1f}% | {delta_counts[m]} |\n")
        f.write("\n")

        f.write("## 2. Boundary Confidence Distribution\n\n")
        f.write("| Confidence | Count | % |\n")
        f.write("|-----------|-------|---|\n")
        for level in ["high", "medium", "low"]:
            pct = confidence_dist[level] / total * 100
            f.write(f"| {level} | {confidence_dist[level]} | {pct:.1f}% |\n")
        f.write("\n")

        f.write("## 3. Coordinate Module Distributions\n\n")
        f.write(f"- **N_ext length**: {percentiles(n_ext_vals)}\n")
        f.write(f"- **α2β3 insert gap**: {percentiles(insert_vals)}\n")
        f.write(f"- **C_tail length**: {percentiles(c_tail_vals)}\n\n")

        f.write("## 4. Thresholds Applied\n\n")
        f.write("| Module | Strict | Relaxed |\n")
        f.write("|--------|--------|--------|\n")
        for m in modules:
            s_val = thresholds[m]["strict"]
            r_val = thresholds[m]["relaxed"]
            if m in ("ACT_domain", "CM_domain"):
                f.write(f"| {m} | i-Evalue ≤ {s_val} | i-Evalue ≤ {r_val} |\n")
            elif m == "C_tail":
                f.write(f"| {m} | ≥ {s_val} aa & no HMM hit | ≥ {r_val} aa |\n")
            else:
                if m == "alpha2beta3_insert":
                    f.write(f"| {m} | gap ≥ {s_val} aa | gap ≥ {r_val} aa |\n")
                else:
                    f.write(f"| {m} | ≥ {s_val} aa | ≥ {r_val} aa |\n")
        f.write("\n")

        f.write("## 5. Strict ⊆ Relaxed Consistency\n\n")
        violations = 0
        for r in results:
            for m in modules:
                if r["strict"][m] == 1 and r["relaxed"][m] == 0:
                    violations += 1
        if violations == 0:
            f.write("✅ **PASS**: All strict=1 entries are also relaxed=1 "
                    "(strict ⊆ relaxed holds for all modules).\n")
        else:
            f.write(f"⚠️ **FAIL**: {violations} violations of strict ⊆ relaxed.\n")
        f.write("\n")

    print(f"[annotate_modules] Wrote robustness report: {path}", file=sys.stderr)


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    coord_rows = load_coords(args.coords)
    print(f"[annotate_modules] Loaded {len(coord_rows)} coordinate records",
          file=sys.stderr)

    hmm_hits = load_domtbl(args.domtbl)
    if hmm_hits:
        n_seqs_with_hits = len(hmm_hits)
        n_total_hits = sum(len(v) for v in hmm_hits.values())
        print(f"[annotate_modules] Loaded {n_total_hits} HMM hits for "
              f"{n_seqs_with_hits} sequences", file=sys.stderr)
    else:
        print("[annotate_modules] No HMM hits loaded (ACT/CM will be 0 for all)",
              file=sys.stderr)

    thresholds = load_thresholds(args.params)

    results, modules = annotate_all(coord_rows, hmm_hits, thresholds)

    # Write strict matrix
    strict_path = os.path.join(args.outdir, "module_presence_absence_strict.tsv")
    write_matrix(results, modules, strict_path, "strict")

    # Write relaxed matrix
    relaxed_path = os.path.join(args.outdir, "module_presence_absence_relaxed.tsv")
    write_matrix(results, modules, relaxed_path, "relaxed")

    # Write robustness report
    write_robustness_report(results, modules, thresholds, args.outdir)


if __name__ == "__main__":
    main()
