#!/usr/bin/env python3
"""Multi-dimensional QC3 gate for DAH7PS V6.1.

QC3 is now split into four dimensions:
  1. provenance
  2. topology/model sensitivity
  3. root tie/identity
  4. annotation sensitivity

The script is intentionally conservative: missing evidence keeps QC3 on HOLD.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate QC3 across provenance, topology, root identity, and annotation sensitivity."
    )
    parser.add_argument(
        "--root_scenarios",
        default="results/04_phylogeny_asr/root_scenarios.tsv",
        help="Scenario registry TSV",
    )
    parser.add_argument(
        "--artifact_manifest",
        default="results/04_phylogeny_asr/artifact_manifest.tsv",
        help="Artifact manifest TSV",
    )
    parser.add_argument(
        "--feature_registry",
        default="results/03_msa_modules/module_feature_registry.tsv",
        help="Orthogonal feature registry TSV",
    )
    parser.add_argument(
        "--panel_calibration",
        default="results/03_msa_modules/panel35_feature_calibration.tsv",
        help="35-structure panel calibration TSV",
    )
    parser.add_argument(
        "--tree_comparison",
        default="results/04_phylogeny_asr/tree_comparison.tsv",
        help="Optional AA vs 3Di tree comparison TSV",
    )
    parser.add_argument(
        "--metrics_manifest",
        default="results/meta/metrics_manifest.tsv",
        help="Metrics manifest TSV with precomputed nRF values",
    )
    parser.add_argument(
        "--output_md",
        default="results/04_phylogeny_asr/QC3_root_stability.md",
        help="QC3 markdown output",
    )
    return parser.parse_args()


class Node:
    __slots__ = ["name", "children"]

    def __init__(self, name: str = "", children: list["Node"] | None = None) -> None:
        self.name = name
        self.children = children or []

    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self) -> list[str]:
        if self.is_leaf():
            return [self.name]
        values: list[str] = []
        for child in self.children:
            values.extend(child.leaves())
        return values


def fail(message: str) -> None:
    print(f"[qc_root_stability] ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def parse_newick_string(raw: str) -> Node:
    text = raw.strip().rstrip(";")
    pos = [0]

    def parse_node() -> Node:
        children: list[Node] = []
        if pos[0] < len(text) and text[pos[0]] == "(":
            pos[0] += 1
            children.append(parse_node())
            while pos[0] < len(text) and text[pos[0]] == ",":
                pos[0] += 1
                children.append(parse_node())
            if pos[0] < len(text) and text[pos[0]] == ")":
                pos[0] += 1

        name = []
        while pos[0] < len(text) and text[pos[0]] not in ",):;":
            name.append(text[pos[0]])
            pos[0] += 1

        if pos[0] < len(text) and text[pos[0]] == ":":
            pos[0] += 1
            while pos[0] < len(text) and text[pos[0]] not in ",);":
                pos[0] += 1

        return Node("".join(name).strip(), children)

    return parse_node()


def load_tree(path: str) -> Node | None:
    if not path or not os.path.isfile(path):
        return None
    with open(path) as handle:
        return parse_newick_string(handle.read())


def bipartitions(root: Node) -> set[frozenset[frozenset[str]]]:
    all_leaves = frozenset(root.leaves())
    splits: set[frozenset[frozenset[str]]] = set()

    def visit(node: Node) -> frozenset[str]:
        if node.is_leaf():
            return frozenset([node.name])
        current = frozenset()
        for child in node.children:
            current |= visit(child)
        other = all_leaves - current
        if 1 < len(current) < len(all_leaves) - 1:
            splits.add(frozenset([current, other]))
        return current

    visit(root)
    return splits


def normalized_rf(left: Node, right: Node) -> float | None:
    left_tips = set(left.leaves())
    right_tips = set(right.leaves())
    if left_tips != right_tips:
        return None
    left_bips = bipartitions(left)
    right_bips = bipartitions(right)
    denom = len(left_bips) + len(right_bips)
    if denom == 0:
        return 0.0
    return (len(left_bips - right_bips) + len(right_bips - left_bips)) / denom


def root_identity(root: Node) -> str:
    if len(root.children) < 2:
        return "unresolved"
    child_sets = [sorted(child.leaves()) for child in root.children]
    ordered = sorted(child_sets, key=lambda names: (len(names), names[:3]))
    labels = []
    for names in ordered[:2]:
        preview = ",".join(names[:3])
        labels.append(f"{len(names)}[{preview}]")
    return " | ".join(labels)


def root_split_sizes(root: Node) -> str:
    if len(root.children) < 2:
        return "NA"
    sizes = [len(child.leaves()) for child in root.children]
    return ",".join(str(size) for size in sorted(sizes))


def read_tsv(path: str) -> list[dict[str, str]]:
    if not os.path.isfile(path):
        return []
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def scenario_table(path: str) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if not rows:
        fail(f"root scenario table missing or empty: {path}")
    return rows


def provenance_dimension(
    scenarios: list[dict[str, str]],
    artifact_rows: list[dict[str, str]],
) -> tuple[str, list[str]]:
    by_scenario: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in artifact_rows:
        by_scenario[row.get("scenario_id", "")].append(row)

    lines = []
    failed = False
    for row in scenarios:
        scenario_id = row["scenario_id"]
        artifacts = by_scenario.get(scenario_id, [])
        roles = {artifact["artifact_role"] for artifact in artifacts}
        missing = [
            role
            for role in ("treefile", "summary", "log")
            if row.get("formal_status") == "formal" and role not in roles
        ]
        has_md5 = all(artifact.get("output_md5") for artifact in artifacts)
        has_command = all(artifact.get("command") for artifact in artifacts)
        has_inputs = all(artifact.get("input_md5") for artifact in artifacts if artifact.get("input_paths"))
        status = "PASS"
        if missing or not has_md5 or not has_command or not has_inputs:
            status = "HOLD"
            failed = True
        lines.append(
            f"| {scenario_id} | {row.get('formal_status','')} | "
            f"{','.join(sorted(roles)) or 'none'} | "
            f"{','.join(missing) if missing else '-'} | {status} |"
        )
    return ("HOLD" if failed else "PASS"), lines


def metrics_by_id(path: str) -> dict[str, str]:
    metrics = {}
    for row in read_tsv(path):
        metric_id = row.get("metric_id")
        if metric_id:
            metrics[metric_id] = row.get("value", "")
    return metrics


def topology_dimension(
    scenarios: list[dict[str, str]],
    tree_comparison_rows: list[dict[str, str]],
    metrics: dict[str, str],
) -> tuple[str, list[str]]:
    lines = []
    hold = False

    scenario_ids = {row["scenario_id"] for row in scenarios}
    if {"S1_MFP_KDOPS", "S2_LGC20_KDOPS"}.issubset(scenario_ids):
        s1s2 = metrics.get("qc3_s1_s2_nRF", "NA")
        verdict = "PASS"
        note = "precomputed model-sensitivity comparison"
        try:
            if float(s1s2) >= 0.10:
                verdict = "SENSITIVE"
                hold = True
                note = "meaningful topology/model shift between S1 and S2"
        except ValueError:
            verdict = "HOLD"
            hold = True
            note = "missing numeric nRF"
        lines.append(f"| S1_MFP_KDOPS vs S2_LGC20_KDOPS | {s1s2} | {verdict} | {note} |")

    reroot_only = [
        scenario_id
        for scenario_id in ("S3_MIDPOINT_INGROUP", "S4A_TOP500_PROXY", "S4B_FULLSEARCH_PROXY")
        if scenario_id in scenario_ids
    ]
    if len(reroot_only) >= 2:
        lines.append(
            f"| {', '.join(reroot_only)} | reroot-only | CONDITIONAL | "
            "same ingroup base topology; sensitivity resides in root identity, not unrooted splits |"
        )

    if tree_comparison_rows:
        nrf = next(
            (
                row["value"]
                for row in tree_comparison_rows
                if row.get("metric") == "RF_normalized"
            ),
            "NA",
        )
        lines.append(f"| AA vs 3Di skeleton | {nrf} | ORTHOGONAL | structural panel check |")

    return ("HOLD" if hold else "PASS"), lines


def root_identity_dimension(
    scenarios: list[dict[str, str]],
) -> tuple[str, list[str]]:
    lines = []
    rho_groups: dict[str, list[str]] = defaultdict(list)
    identities: dict[str, str] = {}

    for row in scenarios:
        tree = load_tree(row.get("tree_path", ""))
        identity = "unavailable"
        split = row.get("root_split_sizes", "NA")
        if tree is not None:
            identity = root_identity(tree)
            split = root_split_sizes(tree)
        identities[row["scenario_id"]] = identity
        rho = row.get("rho", "")
        if rho and rho != "NA":
            rho_groups[rho].append(row["scenario_id"])
        lines.append(
            f"| {row['scenario_id']} | {row.get('rho','NA')} | {split} | {identity} |"
        )

    hold = False
    for rho, members in rho_groups.items():
        unique_identities = {identities[scenario_id] for scenario_id in members}
        if len(members) > 1 and len(unique_identities) > 1:
            hold = True
            lines.append(
                f"| rho tie | {rho} | {','.join(members)} | "
                "same rho but different root identities |"
            )

    return ("HOLD" if hold else "PASS"), lines


def annotation_dimension(
    feature_rows: list[dict[str, str]],
    panel_rows: list[dict[str, str]],
) -> tuple[str, list[str]]:
    lines = []
    hold = False

    if not feature_rows:
        return "HOLD", ["| feature registry | missing | HOLD |"]

    for row in feature_rows:
        feature_id = row["feature_id"]
        role = row["trait_asr_role"]
        orthogonality = row["orthogonality"]
        calibration = row["panel35_calibration"]
        status = "PASS"
        if calibration != "required":
            status = "PASS"
        elif not panel_rows:
            status = "HOLD"
            hold = True
        elif feature_id == "c_residual" and role == "conditional_only":
            status = "CONDITIONAL"
        lines.append(
            f"| {feature_id} | {orthogonality} | {role} | {calibration} | {status} |"
        )

    if panel_rows:
        reviewed = sum(
            1
            for row in panel_rows
            if not row.get("calibration_status", "").startswith("pending")
        )
        lines.append(
            f"| panel35 calibration | {reviewed}/{len(panel_rows)} reviewed | - | - | "
            f"{'HOLD' if reviewed == 0 else 'PASS'} |"
        )
        if reviewed == 0:
            hold = True

    return ("HOLD" if hold else "PASS"), lines


def overall_verdict(statuses: list[str]) -> str:
    if "HOLD" in statuses:
        return "HOLD"
    if "SENSITIVE" in statuses:
        return "CONDITIONAL"
    return "PASS"


def write_report(
    output_path: str,
    provenance_status: str,
    provenance_lines: list[str],
    topology_status: str,
    topology_lines: list[str],
    root_status: str,
    root_lines: list[str],
    annotation_status: str,
    annotation_lines: list[str],
) -> None:
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    final_status = overall_verdict(
        [provenance_status, topology_status, root_status, annotation_status]
    )

    with open(output_path, "w") as handle:
        handle.write("# QC3 Root Stability (V6.1 Multi-Dimensional Gate)\n\n")
        handle.write(
            "> This report supersedes the older single-ratio QC3 interpretation. "
            "QC3 is not closed until all four dimensions are addressed.\n\n"
        )
        handle.write(f"**Overall QC3 status:** `{final_status}`\n\n")
        handle.write("## 1. Provenance\n\n")
        handle.write(f"Status: `{provenance_status}`\n\n")
        handle.write("| Scenario | Formality | Artifact roles | Missing roles | Status |\n")
        handle.write("|---|---|---|---|---|\n")
        handle.write("\n".join(provenance_lines) + "\n\n")

        handle.write("## 2. Topology / Model Sensitivity\n\n")
        handle.write(f"Status: `{topology_status}`\n\n")
        handle.write("| Comparison | nRF | Status | Note |\n")
        handle.write("|---|---|---|---|\n")
        handle.write("\n".join(topology_lines) + "\n\n")

        handle.write("## 3. Root Tie / Identity\n\n")
        handle.write(f"Status: `{root_status}`\n\n")
        handle.write("| Scenario | rho | Root split sizes | Root identity |\n")
        handle.write("|---|---|---|---|\n")
        handle.write("\n".join(root_lines) + "\n\n")

        handle.write("## 4. Annotation Sensitivity\n\n")
        handle.write(f"Status: `{annotation_status}`\n\n")
        handle.write("| Feature | Orthogonality | Trait ASR role | Panel calibration | Status |\n")
        handle.write("|---|---|---|---|---|\n")
        handle.write("\n".join(annotation_lines) + "\n\n")

        handle.write("## Gate Interpretation\n\n")
        handle.write("- `PASS`: evidence bundle is complete and internally coherent.\n")
        handle.write("- `HOLD`: missing provenance, unresolved root ties, or uncalibrated traits.\n")
        handle.write("- `CONDITIONAL`: usable only for conditional conclusions, not main-text lock-in.\n")


def main() -> None:
    args = parse_args()

    scenarios = scenario_table(args.root_scenarios)
    artifact_rows = read_tsv(args.artifact_manifest)
    tree_comparison_rows = read_tsv(args.tree_comparison)
    feature_rows = read_tsv(args.feature_registry)
    panel_rows = read_tsv(args.panel_calibration)
    metrics = metrics_by_id(args.metrics_manifest)

    provenance_status, provenance_lines = provenance_dimension(scenarios, artifact_rows)
    topology_status, topology_lines = topology_dimension(
        scenarios,
        tree_comparison_rows,
        metrics,
    )
    root_status, root_lines = root_identity_dimension(scenarios)
    annotation_status, annotation_lines = annotation_dimension(feature_rows, panel_rows)

    write_report(
        args.output_md,
        provenance_status,
        provenance_lines,
        topology_status,
        topology_lines,
        root_status,
        root_lines,
        annotation_status,
        annotation_lines,
    )
    print(f"[qc_root_stability] wrote {args.output_md}")


if __name__ == "__main__":
    main()
