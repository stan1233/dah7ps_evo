#!/usr/bin/env python3
"""Compare two rooted Newick trees with identical tip sets."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Node:
    name: str = ""
    children: list["Node"] = field(default_factory=list)

    def is_leaf(self) -> bool:
        return not self.children


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare two rooted Newick trees.")
    parser.add_argument("--left", required=True, help="Left/reference Newick tree")
    parser.add_argument("--right", required=True, help="Right/query Newick tree")
    parser.add_argument("--left-label", default="left", help="Label for left tree")
    parser.add_argument("--right-label", default="right", help="Label for right tree")
    parser.add_argument("--output-tsv", required=True, help="Comparison TSV output")
    parser.add_argument("--output-md", required=True, help="Markdown summary output")
    parser.add_argument(
        "--membership-fasta",
        action="append",
        default=[],
        metavar="GROUP=PATH",
        help="Optional group membership FASTA, e.g. Ia=results/02_qc/nr80_Ia.fasta",
    )
    return parser.parse_args()


def fail(message: str) -> None:
    print(f"[compare_rooted_tree_pairs] ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def read_text(path: Path) -> str:
    if not path.is_file():
        fail(f"tree not found: {path}")
    return path.read_text().strip()


def parse_newick(raw: str) -> Node:
    text = raw.strip().rstrip(";")
    pos = 0

    def parse_node() -> Node:
        nonlocal pos
        children: list[Node] = []
        if pos < len(text) and text[pos] == "(":
            pos += 1
            children.append(parse_node())
            while pos < len(text) and text[pos] == ",":
                pos += 1
                children.append(parse_node())
            if pos >= len(text) or text[pos] != ")":
                fail("malformed Newick: missing closing parenthesis")
            pos += 1

        name_chars: list[str] = []
        while pos < len(text) and text[pos] not in ",):;":
            name_chars.append(text[pos])
            pos += 1

        if pos < len(text) and text[pos] == ":":
            pos += 1
            while pos < len(text) and text[pos] not in ",);":
                pos += 1

        return Node("".join(name_chars).strip(), children)

    return parse_node()


def leaves(node: Node) -> list[str]:
    if node.is_leaf():
        return [node.name]
    result: list[str] = []
    for child in node.children:
        result.extend(leaves(child))
    return result


def bipartitions(root: Node) -> set[frozenset[frozenset[str]]]:
    all_leaves = frozenset(leaves(root))
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


def normalized_rf(left: Node, right: Node) -> tuple[int, int, float]:
    left_bips = bipartitions(left)
    right_bips = bipartitions(right)
    rf = len(left_bips - right_bips) + len(right_bips - left_bips)
    denom = len(left_bips) + len(right_bips)
    nrf = 0.0 if denom == 0 else rf / denom
    return rf, denom, nrf


def root_child_sets(root: Node) -> list[list[str]]:
    return [sorted(leaves(child)) for child in root.children]


def root_identity(root: Node) -> str:
    child_sets = sorted(root_child_sets(root), key=lambda names: (len(names), names[:3]))
    labels = []
    for names in child_sets:
        labels.append(f"{len(names)}[{','.join(names[:3])}]")
    return " | ".join(labels) if labels else "unresolved"


def root_split_sizes(root: Node) -> str:
    sizes = sorted(len(child) for child in root_child_sets(root))
    return ",".join(str(size) for size in sizes) if sizes else "NA"


def read_fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


def load_membership(raw_items: list[str]) -> dict[str, str]:
    membership: dict[str, str] = {}
    for item in raw_items:
        if "=" not in item:
            fail(f"membership item must be GROUP=PATH: {item}")
        group, raw_path = item.split("=", 1)
        for seq_id in read_fasta_ids(Path(raw_path)):
            membership[seq_id] = group
    return membership


def group_counts(names: list[str], membership: dict[str, str]) -> str:
    if not membership:
        return "NA"
    counts: dict[str, int] = {}
    for name in names:
        group = membership.get(name, "Other")
        counts[group] = counts.get(group, 0) + 1
    return ";".join(f"{group}:{counts[group]}" for group in sorted(counts))


def write_outputs(
    args: argparse.Namespace,
    left: Node,
    right: Node,
    left_tips: set[str],
    right_tips: set[str],
    membership: dict[str, str],
) -> None:
    if left_tips != right_tips:
        only_left = sorted(left_tips - right_tips)
        only_right = sorted(right_tips - left_tips)
        fail(
            "tip sets differ: "
            f"only_left={len(only_left)} only_right={len(only_right)}"
        )

    rf, denom, nrf = normalized_rf(left, right)
    left_root_children = root_child_sets(left)
    right_root_children = root_child_sets(right)
    root_identity_same = root_identity(left) == root_identity(right)

    output_tsv = Path(args.output_tsv)
    output_md = Path(args.output_md)
    if output_tsv.exists():
        fail(f"output already exists; refusing to overwrite: {output_tsv}")
    if output_md.exists():
        fail(f"output already exists; refusing to overwrite: {output_md}")
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    rows = [
        {"metric": "left_label", "value": args.left_label, "note": ""},
        {"metric": "right_label", "value": args.right_label, "note": ""},
        {"metric": "tips", "value": str(len(left_tips)), "note": "identical tip set"},
        {"metric": "rf_distance", "value": str(rf), "note": f"denominator={denom}"},
        {"metric": "normalized_rf", "value": f"{nrf:.6f}", "note": ""},
        {"metric": "left_root_split_sizes", "value": root_split_sizes(left), "note": ""},
        {"metric": "right_root_split_sizes", "value": root_split_sizes(right), "note": ""},
        {"metric": "left_root_identity", "value": root_identity(left), "note": ""},
        {"metric": "right_root_identity", "value": root_identity(right), "note": ""},
        {"metric": "root_identity_same", "value": str(root_identity_same), "note": ""},
    ]
    for label, child_sets in ((args.left_label, left_root_children), (args.right_label, right_root_children)):
        for index, names in enumerate(sorted(child_sets, key=lambda item: (len(item), item[:3])), start=1):
            rows.append(
                {
                    "metric": f"{label}_root_child_{index}_group_counts",
                    "value": group_counts(names, membership),
                    "note": f"n={len(names)}",
                }
            )

    with output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value", "note"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    output_md.write_text(
        "\n".join(
            [
                "# Rooted Tree Pair Comparison",
                "",
                f"- Left: `{args.left_label}`",
                f"- Right: `{args.right_label}`",
                f"- Tips: `{len(left_tips)}` identical tips",
                f"- RF distance: `{rf} / {denom}`",
                f"- Normalized RF: `{nrf:.6f}`",
                f"- Left root split: `{root_split_sizes(left)}`",
                f"- Right root split: `{root_split_sizes(right)}`",
                f"- Root identity same: `{root_identity_same}`",
                "",
                "Root identities:",
                "",
                f"- Left: `{root_identity(left)}`",
                f"- Right: `{root_identity(right)}`",
                "",
            ]
        )
    )

    print(f"[compare_rooted_tree_pairs] Wrote {output_tsv}")
    print(f"[compare_rooted_tree_pairs] Wrote {output_md}")


def main() -> None:
    args = parse_args()
    left = parse_newick(read_text(Path(args.left)))
    right = parse_newick(read_text(Path(args.right)))
    left_tips = set(leaves(left))
    right_tips = set(leaves(right))
    membership = load_membership(args.membership_fasta)
    write_outputs(args, left, right, left_tips, right_tips, membership)


if __name__ == "__main__":
    main()
