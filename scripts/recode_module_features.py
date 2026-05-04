#!/usr/bin/env python3
"""Recode module annotation into orthogonal features and panel calibration tables.

This keeps the legacy 5-module view available, but downstream trait ASR should
operate on orthogonal features rather than treating C_tail as a primary
evolutionary character.

Outputs:
  - module_feature_registry.tsv
  - module_feature_matrix.tsv
  - panel35_feature_calibration.tsv
  - module_feature_encoding.md
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import defaultdict


DEFAULT_THRESHOLDS = {
    "n_ext_relaxed": 10,
    "n_ext_strict": 25,
    "insert_relaxed": 1,
    "insert_strict": 5,
    "c_ext_relaxed": 10,
    "c_ext_strict": 25,
    "act_relaxed": 1e-3,
    "act_strict": 1e-5,
    "cm_relaxed": 1e-3,
    "cm_strict": 1e-5,
}

PDB_ANCHORS = [
    ("PDB-1KFL", "Ia"),
    ("PDB-1RZM", "Ib"),
    ("PDB-3NV8", "II"),
    ("PDB-5CKV", "II"),
    ("PDB-2B7O", "II"),
]

PDB_ANCHOR_NOTES = {
    "PDB-3NV8": (
        "O53512_AROG_MYCTU state variant; same M. tuberculosis AroG/DAHP "
        "synthase core enzyme as PDB-2B7O and PDB-5CKV."
    ),
    "PDB-5CKV": (
        "O53512_AROG_MYCTU state variant; His-tagged/full feedback-inhibited "
        "crystal state of the same core enzyme as PDB-2B7O and PDB-3NV8."
    ),
    "PDB-2B7O": (
        "O53512_AROG_MYCTU state variant; same M. tuberculosis AroG/DAHP "
        "synthase core enzyme as PDB-3NV8 and PDB-5CKV."
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recode module annotations into orthogonal features."
    )
    parser.add_argument("--coords", required=True, help="core_domain_coords.tsv")
    parser.add_argument("--domtbl", required=True, help="module_hits.domtbl")
    parser.add_argument("--panel_manifest", required=True, help="panel_manifest.tsv")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--params", help="Optional params.json")
    return parser.parse_args()


def fail(message: str) -> None:
    print(f"[recode_module_features] ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def load_thresholds(params_path: str | None) -> dict[str, float]:
    thresholds = dict(DEFAULT_THRESHOLDS)
    if params_path and os.path.isfile(params_path):
        with open(params_path) as handle:
            params = json.load(handle)
        feature_params = params.get("module_feature_encoding", {})
        thresholds.update(feature_params)
    return thresholds


def load_coords(path: str) -> list[dict[str, str]]:
    if not os.path.isfile(path):
        fail(f"coords not found: {path}")
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_domtbl(path: str) -> dict[str, list[dict[str, object]]]:
    if not os.path.isfile(path):
        fail(f"domtbl not found: {path}")
    hits: dict[str, list[dict[str, object]]] = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue
            hits[fields[0]].append(
                {
                    "query": fields[3],
                    "ievalue": float(fields[12]),
                    "env_from": int(fields[19]),
                    "env_to": int(fields[20]),
                }
            )
    return hits


def load_panel_ids(path: str) -> dict[str, dict[str, str]]:
    if not os.path.isfile(path):
        fail(f"panel manifest not found: {path}")
    panel = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            panel[row["rep_id"]] = row
    return panel


def parse_env_segments(raw: str) -> list[tuple[int, int]]:
    result = []
    for part in raw.split(";"):
        part = part.strip()
        if not part or "-" not in part:
            continue
        left, right = part.split("-", 1)
        result.append((int(left), int(right)))
    return result


def largest_gap(raw: str) -> int:
    segments = parse_env_segments(raw)
    if len(segments) < 2:
        return 0
    gaps = [segments[i + 1][0] - segments[i][1] - 1 for i in range(len(segments) - 1)]
    return max(gaps) if gaps else 0


def best_hit(seq_hits: list[dict[str, object]], token: str) -> dict[str, object] | None:
    candidates = [hit for hit in seq_hits if token in str(hit["query"])]
    if not candidates:
        return None
    return min(candidates, key=lambda hit: float(hit["ievalue"]))


def classify_binary(value: bool) -> str:
    return "1" if value else "0"


def main() -> None:
    args = parse_args()
    thresholds = load_thresholds(args.params)
    coords = load_coords(args.coords)
    domtbl = load_domtbl(args.domtbl)
    panel = load_panel_ids(args.panel_manifest)

    os.makedirs(args.outdir, exist_ok=True)

    registry_path = os.path.join(args.outdir, "module_feature_registry.tsv")
    matrix_path = os.path.join(args.outdir, "module_feature_matrix.tsv")
    panel_path = os.path.join(args.outdir, "panel35_feature_calibration.tsv")
    report_path = os.path.join(args.outdir, "module_feature_encoding.md")

    matrix_rows = []
    panel_rows = []

    for row in coords:
        seq_id = row["seq_id"]
        seq_len = int(row["seq_len"])
        raw_env_start = int(row["raw_env_start"])
        raw_env_end = int(row["raw_env_end"])
        n_ext_len = raw_env_start - 1
        c_ext_len = seq_len - raw_env_end
        insert_len = largest_gap(row["env_segments"])

        seq_hits = domtbl.get(seq_id, [])
        act_hit = best_hit(seq_hits, "ACT")
        cm_hit = best_hit(seq_hits, "CM")
        act_relaxed = act_hit is not None and float(act_hit["ievalue"]) <= thresholds["act_relaxed"]
        act_strict = act_hit is not None and float(act_hit["ievalue"]) <= thresholds["act_strict"]
        cm_relaxed = cm_hit is not None and float(cm_hit["ievalue"]) <= thresholds["cm_relaxed"]
        cm_strict = cm_hit is not None and float(cm_hit["ievalue"]) <= thresholds["cm_strict"]

        c_residual_relaxed = c_ext_len >= thresholds["c_ext_relaxed"] and not (act_relaxed or cm_relaxed)
        c_residual_strict = c_ext_len >= thresholds["c_ext_strict"] and not (act_relaxed or cm_relaxed)

        matrix_row = {
            "seq_id": seq_id,
            "n_ext_len": str(n_ext_len),
            "n_ext_relaxed": classify_binary(n_ext_len >= thresholds["n_ext_relaxed"]),
            "n_ext_strict": classify_binary(n_ext_len >= thresholds["n_ext_strict"]),
            "insert_len": str(insert_len),
            "insert_relaxed": classify_binary(insert_len >= thresholds["insert_relaxed"]),
            "insert_strict": classify_binary(insert_len >= thresholds["insert_strict"]),
            "act_hmm_relaxed": classify_binary(act_relaxed),
            "act_hmm_strict": classify_binary(act_strict),
            "cm_hmm_relaxed": classify_binary(cm_relaxed),
            "cm_hmm_strict": classify_binary(cm_strict),
            "c_ext_len": str(c_ext_len),
            "c_ext_relaxed": classify_binary(c_ext_len >= thresholds["c_ext_relaxed"]),
            "c_ext_strict": classify_binary(c_ext_len >= thresholds["c_ext_strict"]),
            "c_residual_relaxed": classify_binary(c_residual_relaxed),
            "c_residual_strict": classify_binary(c_residual_strict),
            "trait_asr_priority": (
                "secondary_residual"
                if c_residual_relaxed or c_residual_strict
                else "primary_orthogonal"
            ),
        }
        matrix_rows.append(matrix_row)

        if seq_id in panel:
            panel_rows.append(
                {
                    "rep_id": seq_id,
                    "subtype": panel[seq_id]["subtype"],
                    "structure_source": panel[seq_id]["structure_source"],
                    "seq_len": panel[seq_id]["seq_len"],
                    "n_ext_len": matrix_row["n_ext_len"],
                    "insert_len": matrix_row["insert_len"],
                    "act_hmm_relaxed": matrix_row["act_hmm_relaxed"],
                    "cm_hmm_relaxed": matrix_row["cm_hmm_relaxed"],
                    "c_ext_len": matrix_row["c_ext_len"],
                    "c_residual_relaxed": matrix_row["c_residual_relaxed"],
                    "calibration_status": "pending_manual_review",
                    "panel_note": "",
                }
            )

    existing_panel_ids = {row["rep_id"] for row in panel_rows}
    for rep_id, subtype in PDB_ANCHORS:
        if rep_id in existing_panel_ids:
            continue
        panel_rows.append(
            {
                "rep_id": rep_id,
                "subtype": subtype,
                "structure_source": "PDB",
                "seq_len": "NA",
                "n_ext_len": "NA",
                "insert_len": "NA",
                "act_hmm_relaxed": "NA",
                "cm_hmm_relaxed": "NA",
                "c_ext_len": "NA",
                "c_residual_relaxed": "NA",
                "calibration_status": "pending_pdb_anchor_mapping",
                "panel_note": PDB_ANCHOR_NOTES.get(
                    rep_id,
                    "Anchor structure retained for manual feature calibration.",
                ),
            }
        )

    with open(registry_path, "w", newline="") as handle:
        fieldnames = [
            "feature_id",
            "feature_class",
            "orthogonality",
            "trait_asr_role",
            "panel35_calibration",
            "note",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            [
                {
                    "feature_id": "n_ext_len",
                    "feature_class": "terminal_extension",
                    "orthogonality": "primary",
                    "trait_asr_role": "main",
                    "panel35_calibration": "required",
                    "note": "Length-based N-terminal extension feature.",
                },
                {
                    "feature_id": "insert_len",
                    "feature_class": "internal_insert",
                    "orthogonality": "primary",
                    "trait_asr_role": "main",
                    "panel35_calibration": "required",
                    "note": "alpha2beta3 insert should be treated independently.",
                },
                {
                    "feature_id": "act_hmm",
                    "feature_class": "domain_hit",
                    "orthogonality": "primary",
                    "trait_asr_role": "main",
                    "panel35_calibration": "required",
                    "note": "ACT-domain HMM support.",
                },
                {
                    "feature_id": "cm_hmm",
                    "feature_class": "domain_hit",
                    "orthogonality": "primary",
                    "trait_asr_role": "main",
                    "panel35_calibration": "required",
                    "note": "CM-domain HMM support.",
                },
                {
                    "feature_id": "c_ext_len",
                    "feature_class": "terminal_extension",
                    "orthogonality": "primary",
                    "trait_asr_role": "main",
                    "panel35_calibration": "required",
                    "note": "C-terminal extension length irrespective of HMM identity.",
                },
                {
                    "feature_id": "c_residual",
                    "feature_class": "residual_bin",
                    "orthogonality": "secondary",
                    "trait_asr_role": "conditional_only",
                    "panel35_calibration": "required",
                    "note": "Residual no-HMM C-tail bucket; not a primary evolutionary character.",
                },
            ]
        )

    with open(matrix_path, "w", newline="") as handle:
        fieldnames = list(matrix_rows[0].keys()) if matrix_rows else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(matrix_rows)

    with open(panel_path, "w", newline="") as handle:
        fieldnames = [
            "rep_id",
            "subtype",
            "structure_source",
            "seq_len",
            "n_ext_len",
            "insert_len",
            "act_hmm_relaxed",
            "cm_hmm_relaxed",
            "c_ext_len",
            "c_residual_relaxed",
            "calibration_status",
            "panel_note",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(panel_rows)

    with open(report_path, "w") as handle:
        handle.write("# Orthogonal Module Feature Encoding\n\n")
        handle.write(
            "Trait ASR should operate on orthogonal features rather than the legacy "
            "`C_tail` residual bin.\n\n"
        )
        handle.write(f"- Total sequences recoded: {len(matrix_rows)}\n")
        handle.write(f"- Structure panel entries for calibration: {len(panel_rows)}\n")
        handle.write("- `c_residual` is explicitly marked secondary/conditional.\n")

    print(
        f"[recode_module_features] wrote {len(matrix_rows)} feature rows and "
        f"{len(panel_rows)} panel calibration rows"
    )


if __name__ == "__main__":
    main()
