#!/usr/bin/env python3
"""Render one circular tree per completed Phase 4 rooting scenario."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import math
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "figures" / ".mplconfig"))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Wedge


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from figures.project_figures_lib import (  # noqa: E402
    FIGURES_DIR,
    SUBTYPE_COLORS,
    TreeNode,
    circular_layout,
    ensure_dir,
    load_tsv,
    parse_newick,
    publication_style,
    read_text,
    save_matplotlib_figure,
    thousands,
    tree_label_to_logical_id,
)


KDOPS_COLOR = "#B6423A"
OTHER_COLOR = "#B8B8B8"
MIXED_COLOR = "#262626"
GROUP_COLORS = {
    "Ia": SUBTYPE_COLORS["Ia"],
    "Ib": SUBTYPE_COLORS["Ib"],
    "II": SUBTYPE_COLORS["II"],
    "KDOPS": KDOPS_COLOR,
    "Other": OTHER_COLOR,
    "Mixed": MIXED_COLOR,
}

FIGURE_ORDER = [
    ("F16", "S1_MFP_KDOPS"),
    ("F17", "S2_LGC20_KDOPS"),
    ("F18", "S3_MIDPOINT_INGROUP"),
    ("F19", "S4A_TOP500_PROXY"),
    ("F20", "S4B_FULLSEARCH_PROXY"),
]

SCENARIO_SHORT = {
    "S1_MFP_KDOPS": "S1",
    "S2_LGC20_KDOPS": "S2",
    "S3_MIDPOINT_INGROUP": "S3",
    "S4A_TOP500_PROXY": "S4a",
    "S4B_FULLSEARCH_PROXY": "S4b",
}

SCENARIO_TITLE = {
    "S1_MFP_KDOPS": "KDOPS + MFP rooted tree",
    "S2_LGC20_KDOPS": "KDOPS + LG+C20 rooted tree",
    "S3_MIDPOINT_INGROUP": "Ingroup midpoint rooted tree",
    "S4A_TOP500_PROXY": "Top-500 reroot proxy tree",
    "S4B_FULLSEARCH_PROXY": "Full-search reroot proxy tree",
}

OUTPUT_MARKDOWN = "PHASE4_SCENARIO_TREES.md"
LEGEND_GROUPS = ["Ia", "Ib", "II", "KDOPS"]
PURE_SHADE_GROUPS = {"Ia", "Ib", "II", "KDOPS"}
PURE_WEDGE_MIN_LEAVES = 40


@dataclass(frozen=True)
class NodeSummary:
    group: str
    n_leaves: int


def read_fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


def backup_if_exists(path: Path) -> Path | None:
    if not path.exists():
        return None
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = path.with_name(f"{path.name}.bak_{timestamp}")
    path.replace(backup_path)
    return backup_path


def load_group_membership() -> dict[str, str]:
    membership: dict[str, str] = {}
    for subtype in ["Ia", "Ib", "II"]:
        fasta_path = REPO_ROOT / "results" / "02_qc" / f"nr80_{subtype}.fasta"
        for seq_id in read_fasta_ids(fasta_path):
            membership[seq_id] = subtype
    return membership


def classify_leaf(label: str, membership: dict[str, str]) -> str:
    logical_id = tree_label_to_logical_id(label)
    if logical_id in membership:
        return membership[logical_id]
    if logical_id.startswith("KDOPS_"):
        return "KDOPS"
    return "Other"


def parse_tree(path: Path) -> TreeNode:
    return parse_newick(read_text(path))


def polar_xy(angle: float, radius: float) -> tuple[float, float]:
    return radius * math.cos(angle), radius * math.sin(angle)


def annotate_groups(node: TreeNode, membership: dict[str, str], cache: dict[int, NodeSummary]) -> NodeSummary:
    if node.is_leaf():
        summary = NodeSummary(classify_leaf(node.name, membership), 1)
        cache[id(node)] = summary
        return summary

    child_summaries = [annotate_groups(child, membership, cache) for child in node.children]
    n_leaves = sum(item.n_leaves for item in child_summaries)
    child_groups = {item.group for item in child_summaries}
    if len(child_groups) == 1:
        group = next(iter(child_groups))
    else:
        group = "Mixed"
    summary = NodeSummary(group, n_leaves)
    cache[id(node)] = summary
    return summary


def draw_outer_tip_ring(ax: plt.Axes, leaf_angles: list[float], leaf_groups: list[str], radius_outer: float, width: float) -> None:
    if not leaf_angles:
        return
    step = (2.0 * math.pi) / len(leaf_angles)
    for angle, group in zip(leaf_angles, leaf_groups):
        theta1 = math.degrees(angle - 0.5 * step)
        theta2 = math.degrees(angle + 0.5 * step)
        ax.add_patch(
            Wedge(
                center=(0, 0),
                r=radius_outer,
                theta1=theta1,
                theta2=theta2,
                width=width,
                facecolor=GROUP_COLORS[group],
                edgecolor="none",
                alpha=0.95,
                zorder=1,
            )
        )


def collect_pure_wedges(
    node: TreeNode,
    parent_group: str,
    span_lookup: dict[int, tuple[float, float]],
    depth_lookup: dict[int, float],
    max_depth: int,
    group_cache: dict[int, NodeSummary],
    wedges: list[tuple[str, float, float, float]],
) -> None:
    summary = group_cache[id(node)]
    if summary.group in PURE_SHADE_GROUPS and summary.group != parent_group and summary.n_leaves >= PURE_WEDGE_MIN_LEAVES:
        radius_inner = 0.10 + 0.78 * (depth_lookup[id(node)] / max(max_depth, 1))
        angle_min, angle_max = span_lookup[id(node)]
        wedges.append((summary.group, radius_inner, angle_min, angle_max))
        next_parent = summary.group
    else:
        next_parent = parent_group
    for child in node.children:
        if not child.is_leaf():
            collect_pure_wedges(child, next_parent, span_lookup, depth_lookup, max_depth, group_cache, wedges)


def add_pure_clade_wedges(
    ax: plt.Axes,
    root: TreeNode,
    span_lookup: dict[int, tuple[float, float]],
    depth_lookup: dict[int, float],
    max_depth: int,
    group_cache: dict[int, NodeSummary],
) -> None:
    wedges: list[tuple[str, float, float, float]] = []
    collect_pure_wedges(root, "Mixed", span_lookup, depth_lookup, max_depth, group_cache, wedges)
    for group, radius_inner, angle_min, angle_max in sorted(wedges, key=lambda item: item[1]):
        ax.add_patch(
            Wedge(
                center=(0, 0),
                r=1.00,
                theta1=math.degrees(angle_min),
                theta2=math.degrees(angle_max),
                width=max(1.00 - radius_inner, 0.01),
                facecolor=GROUP_COLORS[group],
                edgecolor="none",
                alpha=0.08,
                zorder=0,
            )
        )


def draw_tree_lines(
    ax: plt.Axes,
    root: TreeNode,
    angle_lookup: dict[int, float],
    depth_lookup: dict[int, float],
    max_depth: int,
    group_cache: dict[int, NodeSummary],
) -> None:
    radial_segments: list[list[tuple[float, float]]] = []
    radial_colors: list[str] = []
    arc_segments: list[list[tuple[float, float]]] = []
    arc_colors: list[str] = []

    def radius(node: TreeNode) -> float:
        return 0.10 + 0.78 * (depth_lookup[id(node)] / max(max_depth, 1))

    def walk(node: TreeNode) -> None:
        if node.is_leaf():
            return
        node_radius = radius(node)
        child_angles = [angle_lookup[id(child)] for child in node.children]
        start_angle = min(child_angles)
        end_angle = max(child_angles)
        arc_group = group_cache[id(node)].group
        arc_color = GROUP_COLORS.get(arc_group, MIXED_COLOR) if arc_group in PURE_SHADE_GROUPS else MIXED_COLOR
        arc_steps = max(8, min(32, int((end_angle - start_angle) / 0.08)))
        arc_angles = np.linspace(start_angle, end_angle, arc_steps)
        for left, right in zip(arc_angles[:-1], arc_angles[1:]):
            arc_segments.append([polar_xy(left, node_radius), polar_xy(right, node_radius)])
            arc_colors.append(arc_color)

        for child in node.children:
            child_radius = radius(child)
            child_angle = angle_lookup[id(child)]
            child_group = group_cache[id(child)].group
            line_color = GROUP_COLORS.get(child_group, MIXED_COLOR) if child_group in {"Ia", "Ib", "II", "KDOPS"} else MIXED_COLOR
            radial_segments.append([polar_xy(child_angle, node_radius), polar_xy(child_angle, child_radius)])
            radial_colors.append(line_color)
            walk(child)

    walk(root)

    ax.add_collection(LineCollection(arc_segments, colors=arc_colors, linewidths=0.20, zorder=2))
    ax.add_collection(LineCollection(radial_segments, colors=radial_colors, linewidths=0.20, zorder=2))


def add_kdops_labels(
    ax: plt.Axes,
    leaves: list[TreeNode],
    angle_lookup: dict[int, float],
    membership: dict[str, str],
) -> int:
    kdops_count = 0
    for leaf in leaves:
        if classify_leaf(leaf.name, membership) != "KDOPS":
            continue
        kdops_count += 1
        angle = angle_lookup[id(leaf)]
        x_tip, y_tip = polar_xy(angle, 1.01)
        x_mark, y_mark = polar_xy(angle, 1.05)
        x_text, y_text = polar_xy(angle, 1.16)
        angle_deg = math.degrees(angle)
        if -90 <= angle_deg <= 90:
            rotation = angle_deg
            ha = "left"
        else:
            rotation = angle_deg + 180
            ha = "right"
        ax.plot([x_tip, x_mark], [y_tip, y_mark], color=KDOPS_COLOR, lw=0.55, zorder=4)
        ax.scatter([x_mark], [y_mark], s=12, marker="*", color=KDOPS_COLOR, edgecolors="white", linewidths=0.25, zorder=5)
        ax.text(
            x_text,
            y_text,
            leaf.name.replace("KDOPS_", ""),
            fontsize=5.0,
            color=KDOPS_COLOR,
            rotation=rotation,
            rotation_mode="anchor",
            ha=ha,
            va="center",
            zorder=5,
        )
    return kdops_count


def counts_by_group(leaves: list[TreeNode], membership: dict[str, str]) -> dict[str, int]:
    counts = {"Ia": 0, "Ib": 0, "II": 0, "KDOPS": 0, "Other": 0}
    for leaf in leaves:
        counts[classify_leaf(leaf.name, membership)] += 1
    return counts


def scenario_paths() -> list[tuple[str, object, Path]]:
    table = load_tsv(REPO_ROOT / "results" / "04_phylogeny_asr" / "root_scenarios.tsv")
    rows_by_id = {row.scenario_id: row for row in table.itertuples(index=False)}
    paths: list[tuple[str, object, Path]] = []
    for figure_id, scenario_id in FIGURE_ORDER:
        row = rows_by_id[scenario_id]
        tree_path = str(row.tree_path).strip()
        if tree_path.upper() in {"NA", "NAN", ""}:
            raise FileNotFoundError(f"No tree_path registered for {scenario_id}")
        paths.append((figure_id, row, REPO_ROOT / tree_path))
    return paths


def render_scenario_figure(
    figure_id: str,
    row: object,
    tree_path: Path,
    membership: dict[str, str],
) -> dict[str, str]:
    root = parse_tree(tree_path)
    group_cache: dict[int, NodeSummary] = {}
    annotate_groups(root, membership, group_cache)
    angle_lookup, depth_lookup, span_lookup, leaves, max_depth = circular_layout(root)

    fig, ax = plt.subplots(figsize=(10.0, 10.0))
    fig.subplots_adjust(top=0.93, bottom=0.06, left=0.05, right=0.95)

    leaf_angles = [angle_lookup[id(leaf)] for leaf in leaves]
    leaf_groups = [classify_leaf(leaf.name, membership) for leaf in leaves]
    add_pure_clade_wedges(ax, root, span_lookup, depth_lookup, max_depth, group_cache)
    draw_outer_tip_ring(ax, leaf_angles, leaf_groups, radius_outer=1.04, width=0.05)
    draw_tree_lines(ax, root, angle_lookup, depth_lookup, max_depth, group_cache)

    kdops_count = 0
    if str(row.scenario_id) in {"S1_MFP_KDOPS", "S2_LGC20_KDOPS"}:
        kdops_count = add_kdops_labels(ax, leaves, angle_lookup, membership)

    counts = counts_by_group(leaves, membership)
    total = sum(counts.values())
    rho_text = "rho = NA" if str(row.rho) == "NA" else f"rho = {row.rho}"
    counts_text = (
        f"n = {thousands(total)}\n"
        f"Ia {thousands(counts['Ia'])}  Ib {thousands(counts['Ib'])}\n"
        f"II {thousands(counts['II'])}  KDOPS {thousands(counts['KDOPS'])}"
    )
    center_text = [
        figure_id,
        SCENARIO_SHORT[str(row.scenario_id)],
        str(row.label),
        rho_text,
        counts_text,
    ]
    if kdops_count:
        center_text.append(f"KDOPS labeled: {kdops_count}")
    ax.text(
        0.0,
        0.0,
        "\n".join(center_text),
        ha="center",
        va="center",
        fontsize=8.1,
        color="#222222",
        linespacing=1.25,
        bbox={"boxstyle": "round,pad=0.4", "fc": "white", "ec": "#DDDDDD", "alpha": 0.9},
        zorder=6,
    )

    legend_handles = [
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=GROUP_COLORS[group], markersize=7, label=group)
        for group in LEGEND_GROUPS
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.03), ncol=4)
    ax.set_title(f"{figure_id}. {SCENARIO_TITLE[str(row.scenario_id)]}", fontsize=11, pad=14)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.27, 1.27)
    ax.set_ylim(-1.27, 1.27)

    stem = f"{figure_id}_{SCENARIO_SHORT[str(row.scenario_id)]}_circular_rooted_tree"
    backup_if_exists(FIGURES_DIR / f"{stem}.png")
    backup_if_exists(FIGURES_DIR / f"{stem}.pdf")
    return save_matplotlib_figure(fig, stem)


def build_markdown(entries: list[tuple[str, object, dict[str, str]]]) -> str:
    lines = [
        "# Phase 4 Scenario Trees",
        "",
        "This file is autogenerated by `scripts/render_phase4_scenario_trees.py`.",
        "Each completed scenario is rendered as a standalone circular tree using the stored computed `.treefile`.",
        "",
    ]
    for figure_id, row, path_map in entries:
        png_name = Path(path_map["png"]).name
        pdf_name = Path(path_map["pdf"]).name
        kdops_note = ""
        if str(row.scenario_id) in {"S1_MFP_KDOPS", "S2_LGC20_KDOPS"}:
            kdops_note = " KDOPS tips are marked individually in red."
        lines.extend(
            [
                f"## {figure_id}. {SCENARIO_TITLE[str(row.scenario_id)]}",
                f"Files: [PNG](./{png_name}) | [PDF](./{pdf_name})",
                "",
                f"![{figure_id}](./{png_name})",
                "",
                (
                    f"{figure_id} shows the actual computed tree for {SCENARIO_SHORT[str(row.scenario_id)]} "
                    f"({row.label}) in circular layout. Translucent branch-sector shading marks clades dominated "
                    f"by a single subtype, the outer ring preserves leaf-order subtype composition, and the center "
                    f"panel records scenario metadata directly from `root_scenarios.tsv`.{kdops_note}"
                ),
                "",
            ]
        )
    lines.extend(
        [
            "## Pending",
            "",
            "S5 is not rendered because `root_scenarios.tsv` still registers no computed `tree_path` for that scenario.",
            "",
        ]
    )
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    publication_style()
    ensure_dir(FIGURES_DIR)

    # Retire the old multi-panel output from the previous version of this script.
    backup_if_exists(FIGURES_DIR / "F16_phase4_actual_root_scenario_trees.png")
    backup_if_exists(FIGURES_DIR / "F16_phase4_actual_root_scenario_trees.pdf")

    membership = load_group_membership()
    entries: list[tuple[str, object, dict[str, str]]] = []
    for figure_id, row, tree_path in scenario_paths():
        path_map = render_scenario_figure(figure_id, row, tree_path, membership)
        entries.append((figure_id, row, path_map))
        print(figure_id, path_map["png"], path_map["pdf"])

    markdown_path = FIGURES_DIR / OUTPUT_MARKDOWN
    backup_if_exists(markdown_path)
    markdown_path.write_text(build_markdown(entries))
    print("MARKDOWN", markdown_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
