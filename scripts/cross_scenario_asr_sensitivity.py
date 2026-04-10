#!/usr/bin/env python3
"""Quantify cross-scenario ASR sensitivity from IQ-TREE .state files.

Outputs:
  - site-level table: MAP agreement, posterior deltas, information range
  - node-level table: agreement rates and conflict counts
  - scenario coverage summary

Usage:
  python scripts/cross_scenario_asr_sensitivity.py \
    --state S1=results/04_phylogeny_asr/ASR_core.state \
    --state S2=results/04_phylogeny_asr/ASR_core_S2.state \
    --out_prefix results/04_phylogeny_asr/asr_cross_scenario
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict


AA_COLUMNS = [f"p_{aa}" for aa in "ARNDCQEGHILKMFPSTWYV"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare IQ-TREE ASR state files across root/model scenarios."
    )
    parser.add_argument(
        "--state",
        action="append",
        required=True,
        help="Scenario assignment in the form SCENARIO_ID=path/to/file.state",
    )
    parser.add_argument(
        "--node_map",
        help="Optional TSV with columns reference_node,scenario_id,node_id",
    )
    parser.add_argument("--out_prefix", required=True, help="Output prefix")
    parser.add_argument(
        "--min_pp_high_conf",
        type=float,
        default=0.90,
        help="Posterior threshold for counting high-confidence conflicts",
    )
    return parser.parse_args()


def fail(message: str) -> None:
    print(f"[cross_scenario_asr_sensitivity] ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def parse_state_arg(item: str) -> tuple[str, str]:
    if "=" not in item:
        fail(f"--state must be SCENARIO=PATH, got: {item}")
    scenario, path = item.split("=", 1)
    scenario = scenario.strip()
    path = path.strip()
    if not scenario or not path:
        fail(f"invalid --state assignment: {item}")
    if not os.path.isfile(path):
        fail(f"state file not found: {path}")
    return scenario, path


def read_node_map(path: str | None, scenarios: list[str]) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    if path is None:
        return mapping
    if not os.path.isfile(path):
        fail(f"node map not found: {path}")
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = {"reference_node", "scenario_id", "node_id"}
        if reader.fieldnames is None or not expected.issubset(reader.fieldnames):
            fail("node map must contain reference_node, scenario_id, node_id")
        for row in reader:
            ref = row["reference_node"]
            scenario = row["scenario_id"]
            if scenario not in scenarios:
                continue
            mapping.setdefault(ref, {})[scenario] = row["node_id"]
    return mapping


def entropy(probs: list[float]) -> float:
    value = 0.0
    for p in probs:
        if p > 0:
            value -= p * math.log2(p)
    return value


def information_bits(probs: list[float]) -> float:
    return math.log2(len(probs)) - entropy(probs)


def read_state_table(path: str) -> dict[str, dict[int, dict[str, object]]]:
    node_site: dict[str, dict[int, dict[str, object]]] = defaultdict(dict)
    with open(path) as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("#")),
            delimiter="\t",
        )
        expected = {"Node", "Site", "State", *AA_COLUMNS}
        if reader.fieldnames is None or not expected.issubset(reader.fieldnames):
            fail(f"unexpected .state header in {path}")
        for row in reader:
            probs = [float(row[column]) for column in AA_COLUMNS]
            node_site[row["Node"]][int(row["Site"])] = {
                "map_state": row["State"],
                "top_pp": max(probs),
                "entropy_bits": entropy(probs),
                "information_bits": information_bits(probs),
                "probs": probs,
            }
    return node_site


def resolve_reference_nodes(
    state_tables: dict[str, dict[str, dict[int, dict[str, object]]]],
    node_map: dict[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    if node_map:
        return node_map

    common_nodes = None
    for table in state_tables.values():
        nodes = set(table.keys())
        common_nodes = nodes if common_nodes is None else common_nodes & nodes
    if not common_nodes:
        fail("no common node IDs found across scenarios; provide --node_map")

    resolved: dict[str, dict[str, str]] = {}
    for node_id in sorted(common_nodes):
        resolved[node_id] = {scenario: node_id for scenario in state_tables}
    return resolved


def site_union(
    state_tables: dict[str, dict[str, dict[int, dict[str, object]]]],
    per_scenario_nodes: dict[str, str],
) -> list[int]:
    sites: set[int] = set()
    for scenario, node_id in per_scenario_nodes.items():
        sites.update(state_tables[scenario][node_id].keys())
    return sorted(sites)


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def main() -> None:
    args = parse_args()
    scenario_paths = [parse_state_arg(item) for item in args.state]
    scenarios = [scenario for scenario, _ in scenario_paths]
    state_tables = {
        scenario: read_state_table(path) for scenario, path in scenario_paths
    }
    node_map = read_node_map(args.node_map, scenarios)
    reference_nodes = resolve_reference_nodes(state_tables, node_map)

    out_prefix = args.out_prefix
    out_dir = os.path.dirname(out_prefix) or "."
    os.makedirs(out_dir, exist_ok=True)

    site_path = f"{out_prefix}_site.tsv"
    node_path = f"{out_prefix}_node.tsv"
    scenario_path = f"{out_prefix}_scenario.tsv"

    site_rows: list[dict[str, object]] = []
    node_rows: list[dict[str, object]] = []

    for ref_node, per_scenario_nodes in sorted(reference_nodes.items()):
        missing = [scenario for scenario in scenarios if scenario not in per_scenario_nodes]
        if missing:
            continue

        agreement_flags: list[int] = []
        pp_deltas: list[float] = []
        info_ranges: list[float] = []
        high_conflict_sites = 0

        for site in site_union(state_tables, per_scenario_nodes):
            per_scenario = {}
            for scenario, node_id in per_scenario_nodes.items():
                record = state_tables[scenario][node_id].get(site)
                if record is not None:
                    per_scenario[scenario] = record

            if len(per_scenario) < 2:
                continue

            map_states = {scenario: rec["map_state"] for scenario, rec in per_scenario.items()}
            top_pps = {scenario: float(rec["top_pp"]) for scenario, rec in per_scenario.items()}
            info_bits = {
                scenario: float(rec["information_bits"]) for scenario, rec in per_scenario.items()
            }
            prob_vectors = [rec["probs"] for rec in per_scenario.values()]

            map_consistent = int(len(set(map_states.values())) == 1)
            max_pp_delta = max(top_pps.values()) - min(top_pps.values())
            l1_deltas = []
            scenarios_present = sorted(per_scenario)
            for i, left in enumerate(scenarios_present):
                for right in scenarios_present[i + 1 :]:
                    left_probs = per_scenario[left]["probs"]
                    right_probs = per_scenario[right]["probs"]
                    l1_deltas.append(
                        sum(abs(a - b) for a, b in zip(left_probs, right_probs)) / 2.0
                    )
            mean_l1_delta = mean(l1_deltas)
            info_range = max(info_bits.values()) - min(info_bits.values())

            agreement_flags.append(map_consistent)
            pp_deltas.append(max_pp_delta)
            info_ranges.append(info_range)
            if (
                not map_consistent
                and max(top_pps.values()) >= args.min_pp_high_conf
                and min(top_pps.values()) >= args.min_pp_high_conf
            ):
                high_conflict_sites += 1

            site_rows.append(
                {
                    "reference_node": ref_node,
                    "site": site,
                    "scenario_count": len(per_scenario),
                    "scenario_ids": ",".join(scenarios_present),
                    "map_states": ";".join(
                        f"{scenario}:{map_states[scenario]}" for scenario in scenarios_present
                    ),
                    "map_consistent": map_consistent,
                    "max_top_pp_delta": f"{max_pp_delta:.6f}",
                    "mean_l1_posterior_delta": f"{mean_l1_delta:.6f}",
                    "information_bits": ";".join(
                        f"{scenario}:{info_bits[scenario]:.6f}"
                        for scenario in scenarios_present
                    ),
                    "information_range_bits": f"{info_range:.6f}",
                }
            )

        node_rows.append(
            {
                "reference_node": ref_node,
                "scenario_count": len(per_scenario_nodes),
                "site_count": len(agreement_flags),
                "map_consistency_rate": (
                    f"{sum(agreement_flags) / len(agreement_flags):.6f}"
                    if agreement_flags
                    else "NA"
                ),
                "mean_max_top_pp_delta": f"{mean(pp_deltas):.6f}" if pp_deltas else "NA",
                "mean_information_range_bits": (
                    f"{mean(info_ranges):.6f}" if info_ranges else "NA"
                ),
                "high_conflict_sites": high_conflict_sites,
            }
        )

    with open(site_path, "w", newline="") as handle:
        fieldnames = [
            "reference_node",
            "site",
            "scenario_count",
            "scenario_ids",
            "map_states",
            "map_consistent",
            "max_top_pp_delta",
            "mean_l1_posterior_delta",
            "information_bits",
            "information_range_bits",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(site_rows)

    with open(node_path, "w", newline="") as handle:
        fieldnames = [
            "reference_node",
            "scenario_count",
            "site_count",
            "map_consistency_rate",
            "mean_max_top_pp_delta",
            "mean_information_range_bits",
            "high_conflict_sites",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(node_rows)

    with open(scenario_path, "w", newline="") as handle:
        fieldnames = ["scenario_id", "state_path", "node_count"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for scenario, path in scenario_paths:
            writer.writerow(
                {
                    "scenario_id": scenario,
                    "state_path": path,
                    "node_count": len(state_tables[scenario]),
                }
            )

    print(
        "[cross_scenario_asr_sensitivity] "
        f"wrote {len(site_rows)} site rows and {len(node_rows)} node rows"
    )


if __name__ == "__main__":
    main()
