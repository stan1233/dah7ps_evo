#!/usr/bin/env python3
"""Render Phase 1 HMM profile and KDOPS decoy diagnostics.

Outputs three publication-style figure sets into figures/:
  F21_phase1_kdops_profile_logo.{pdf,png}
  F22_phase1_hmm_emission_heatmap.{pdf,png}
  F23_phase1_kdops_competitive_scoring.{pdf,png}
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
import math
import os
from pathlib import Path
import sys

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parents[1] / "figures" / ".mplconfig"),
)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D


REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = REPO_ROOT / "figures"

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
AA_GROUP_COLORS = {
    **{aa: "#3D6FA3" for aa in "AVLIMFWY"},
    **{aa: "#3F8F5B" for aa in "STNQ"},
    **{aa: "#7B5BA7" for aa in "KRH"},
    **{aa: "#C44E52" for aa in "DE"},
    "G": "#777777",
    "P": "#A66A2B",
    "C": "#D9A441",
}
MODEL_COLORS = {
    "Ia": "#355C8A",
    "Ib": "#C46E2B",
    "II": "#3B8B52",
    "KDOPS": "#B6423A",
}


@dataclass(frozen=True)
class HMMProfile:
    label: str
    path: Path
    name: str
    length: int
    nseq: int | None
    amino_acids: list[str]
    emissions: np.ndarray
    positions: np.ndarray

    @property
    def information_bits(self) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            entropy = -np.nansum(
                np.where(self.emissions > 0, self.emissions * np.log2(self.emissions), 0.0),
                axis=1,
            )
        return np.log2(len(self.amino_acids)) - entropy

    @property
    def weighted_emissions(self) -> np.ndarray:
        return self.emissions * self.information_bits[:, None]


def rel(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def parse_neglog_probability(token: str) -> float:
    if token == "*":
        return 0.0
    return math.exp(-float(token))


def parse_hmmer_profile(path: Path, label: str) -> HMMProfile:
    if not path.is_file():
        raise FileNotFoundError(path)

    name = path.stem
    length = None
    nseq = None
    amino_acids: list[str] | None = None
    positions: list[int] = []
    emissions: list[list[float]] = []

    with path.open() as handle:
        in_profile = False
        for raw in handle:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped == "//":
                break
            if stripped.startswith("NAME"):
                parts = stripped.split()
                if len(parts) > 1:
                    name = parts[1]
                continue
            if stripped.startswith("LENG"):
                length = int(stripped.split()[1])
                continue
            if stripped.startswith("NSEQ"):
                nseq = int(stripped.split()[1])
                continue
            if stripped.startswith("HMM "):
                parts = stripped.split()
                amino_acids = parts[1:]
                if amino_acids != AA_ORDER:
                    raise ValueError(
                        f"Unexpected amino-acid order in {path}: {' '.join(amino_acids)}"
                    )
                in_profile = True
                continue
            if not in_profile:
                continue

            parts = stripped.split()
            if not parts or not parts[0].isdigit():
                continue
            if amino_acids is None:
                raise ValueError(f"Missing HMM alphabet line in {path}")
            if len(parts) < 1 + len(amino_acids):
                raise ValueError(f"Malformed match-state line in {path}: {stripped}")
            positions.append(int(parts[0]))
            emissions.append(
                [parse_neglog_probability(token) for token in parts[1 : 1 + len(amino_acids)]]
            )

    if amino_acids is None or length is None:
        raise ValueError(f"Could not parse HMMER profile header from {path}")
    if len(emissions) != length:
        raise ValueError(
            f"Parsed {len(emissions)} match states from {path}, expected LENG={length}"
        )

    matrix = np.asarray(emissions, dtype=float)
    row_sums = matrix.sum(axis=1, keepdims=True)
    matrix = np.divide(matrix, row_sums, out=np.zeros_like(matrix), where=row_sums > 0)

    return HMMProfile(
        label=label,
        path=path,
        name=name,
        length=length,
        nseq=nseq,
        amino_acids=amino_acids,
        emissions=matrix,
        positions=np.asarray(positions, dtype=int),
    )


def backup_if_exists(path: Path) -> None:
    if not path.exists():
        return
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = path.with_name(f"{path.name}.bak_{timestamp}")
    path.replace(backup)


def save_figure(fig: mpl.figure.Figure, stem: str) -> dict[str, Path]:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    outputs = {
        "pdf": FIGURES_DIR / f"{stem}.pdf",
        "png": FIGURES_DIR / f"{stem}.png",
    }
    for output in outputs.values():
        backup_if_exists(output)
    fig.savefig(outputs["pdf"], bbox_inches="tight", pad_inches=0.04)
    fig.savefig(outputs["png"], dpi=600, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return outputs


def set_publication_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 600,
        }
    )


def draw_letter(
    ax: mpl.axes.Axes,
    letter: str,
    x_left: float,
    y_bottom: float,
    width: float,
    height: float,
    color: str,
    font: FontProperties,
) -> None:
    if height <= 0:
        return
    path = TextPath((0, 0), letter, size=1.0, prop=font)
    bbox = path.get_extents()
    if bbox.width <= 0 or bbox.height <= 0:
        return
    transform = (
        Affine2D()
        .translate(-bbox.x0, -bbox.y0)
        .scale(width / bbox.width, height / bbox.height)
        .translate(x_left, y_bottom)
        + ax.transData
    )
    patch = PathPatch(path, transform=transform, facecolor=color, edgecolor="none")
    ax.add_patch(patch)


def draw_logo_panel(
    ax: mpl.axes.Axes,
    profile: HMMProfile,
    start: int,
    stop: int,
    max_bits: float,
    min_letter_bits: float,
    show_xticklabels: bool,
) -> None:
    font = FontProperties(family="DejaVu Sans", weight="bold")
    weighted = profile.weighted_emissions[start:stop]
    positions = profile.positions[start:stop]

    for column_index, heights in enumerate(weighted):
        x_left = column_index + 0.075
        y_bottom = 0.0
        for aa_index in np.argsort(heights):
            height = float(heights[aa_index])
            if height < min_letter_bits:
                continue
            aa = profile.amino_acids[aa_index]
            draw_letter(
                ax,
                aa,
                x_left=x_left,
                y_bottom=y_bottom,
                width=0.85,
                height=height,
                color=AA_GROUP_COLORS.get(aa, "#333333"),
                font=font,
            )
            y_bottom += height

    ax.set_xlim(0, max(1, stop - start))
    ax.set_ylim(0, max_bits)
    ax.set_ylabel("bits")
    tick_positions = []
    tick_labels = []
    for offset, position in enumerate(positions):
        if position == positions[0] or position == positions[-1] or position % 10 == 0:
            tick_positions.append(offset + 0.5)
            tick_labels.append(str(position))
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels if show_xticklabels else [])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.4, alpha=0.7)
    ax.text(
        0.01,
        0.94,
        f"HMM states {positions[0]}-{positions[-1]}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="#333333",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
    )


def render_kdops_profile_logo(profile: HMMProfile, stem: str) -> dict[str, Path]:
    window = 70
    n_rows = int(math.ceil(profile.length / window))
    max_bits = math.log2(len(profile.amino_acids))
    fig_height = 1.95 * n_rows + 1.05
    fig, axes = plt.subplots(n_rows, 1, figsize=(12.5, fig_height), sharey=True)
    axes_array = np.asarray(axes).reshape(-1)

    for row, ax in enumerate(axes_array):
        start = row * window
        stop = min((row + 1) * window, profile.length)
        if start >= profile.length:
            ax.set_visible(False)
            continue
        draw_logo_panel(
            ax,
            profile=profile,
            start=start,
            stop=stop,
            max_bits=max_bits,
            min_letter_bits=0.025,
            show_xticklabels=row == len(axes_array) - 1,
        )

    axes_array[-1].set_xlabel("KDOPS HMM match-state position")
    fig.subplots_adjust(left=0.055, right=0.995, top=0.91, bottom=0.085, hspace=0.38)
    fig.suptitle(
        f"KDOPS decoy HMM profile logo ({profile.name}; L={profile.length}, NSEQ={profile.nseq})",
        y=0.975,
        fontsize=11,
        fontweight="bold",
    )
    return save_figure(fig, stem)


def render_emission_heatmap(profiles: list[HMMProfile], stem: str) -> dict[str, Path]:
    max_value = max(float(profile.weighted_emissions.max()) for profile in profiles)
    cmap = LinearSegmentedColormap.from_list(
        "hmm_emission",
        ["#FFFFFF", "#E7EEF5", "#8FB5CF", "#2F5F8A", "#1E2F4A"],
    )

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2), sharey=False)
    axes_array = axes.reshape(-1)
    image = None
    for index, profile in enumerate(profiles):
        ax = axes_array[index]
        matrix = profile.weighted_emissions.T
        image = ax.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            vmin=0,
            vmax=max_value,
            extent=(1, profile.length, len(profile.amino_acids) - 0.5, -0.5),
        )
        ax.set_yticks(np.arange(len(profile.amino_acids)))
        ax.set_yticklabels(profile.amino_acids)
        ax.tick_params(axis="y", labelsize=6, pad=1)
        ax.set_ylabel("")
        ax.set_title(
            f"{profile.label} ({profile.name}; L={profile.length}, NSEQ={profile.nseq})",
            loc="left",
            fontsize=9,
            color=MODEL_COLORS.get(profile.label, "#333333"),
            fontweight="bold",
            pad=4,
        )
        ax.set_xlim(1, profile.length)
        ax.set_xlabel("HMM match-state position")
        ax.tick_params(axis="both", length=2.5)
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

    for ax in axes_array[len(profiles) :]:
        ax.set_visible(False)

    fig.subplots_adjust(left=0.065, right=0.985, top=0.9, bottom=0.14, wspace=0.12, hspace=0.48)
    cax = fig.add_axes([0.28, 0.06, 0.44, 0.025])
    if image is not None:
        colorbar = fig.colorbar(image, cax=cax, orientation="horizontal")
        colorbar.set_label("information-weighted emission (bits)", labelpad=3)
    fig.suptitle(
        "Phase 1 HMM emission heatmaps",
        y=0.975,
        fontsize=11,
        fontweight="bold",
    )
    return save_figure(fig, stem)


def render_kdops_competitive_scoring(report_path: Path, stem: str) -> dict[str, Path]:
    if not report_path.is_file():
        raise FileNotFoundError(report_path)
    data = pd.read_csv(report_path, sep="\t")
    required = {"seq_id", "dah7ps_score", "kdops_score", "delta", "verdict"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing columns in {report_path}: {', '.join(sorted(missing))}")

    for column in ["dah7ps_score", "kdops_score", "delta"]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna(subset=["dah7ps_score", "kdops_score", "delta", "verdict"]).copy()
    data["verdict"] = data["verdict"].astype(str)

    keep = data[data["verdict"] == "KEEP"]
    remove = data[data["verdict"] == "REMOVE"]
    borderline = data[(data["delta"] >= 0) & (data["delta"] <= 20)]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.1), gridspec_kw={"width_ratios": [1, 1.05]})
    ax_hist, ax_scatter = axes

    bins = np.linspace(float(data["delta"].min()), float(data["delta"].max()), 90)
    ax_hist.hist(
        [remove["delta"], keep["delta"]],
        bins=bins,
        stacked=True,
        color=[MODEL_COLORS["KDOPS"], MODEL_COLORS["Ib"]],
        label=[f"REMOVE ({len(remove):,})", f"KEEP ({len(keep):,})"],
        linewidth=0,
        alpha=0.92,
    )
    ax_hist.axvline(0, color="#222222", linewidth=1.0, linestyle="--")
    ax_hist.axvspan(0, 20, color="#E3C34D", alpha=0.22, linewidth=0)
    ax_hist.set_yscale("log")
    ax_hist.set_xlabel("Delta bits: DAH7PS Ib score - KDOPS score")
    ax_hist.set_ylabel("Sequences (log scale)")
    ax_hist.set_title("Competitive score delta", loc="left", fontweight="bold")
    ax_hist.legend(frameon=False, loc="upper left")
    ax_hist.text(
        0.98,
        0.92,
        f"0-20 borderline: {len(borderline):,}",
        transform=ax_hist.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "#BDBDBD", "linewidth": 0.6, "pad": 2.5},
    )

    scatter_order = [("REMOVE", remove, MODEL_COLORS["KDOPS"]), ("KEEP", keep, MODEL_COLORS["Ib"])]
    for verdict, subset, color in scatter_order:
        ax_scatter.scatter(
            subset["dah7ps_score"],
            subset["kdops_score"],
            s=5,
            alpha=0.32,
            linewidths=0,
            color=color,
            label=verdict,
            rasterized=True,
        )
    score_max = float(max(data["dah7ps_score"].max(), data["kdops_score"].max()))
    ax_scatter.plot([0, score_max], [0, score_max], color="#222222", linestyle="--", linewidth=1.0)
    ax_scatter.fill_between(
        [0, score_max],
        [0, score_max],
        [20, score_max + 20],
        color="#E3C34D",
        alpha=0.18,
        linewidth=0,
    )
    ax_scatter.set_xlim(0, score_max * 1.03)
    ax_scatter.set_ylim(0, score_max * 1.03)
    ax_scatter.set_aspect("equal", adjustable="box")
    ax_scatter.set_xlabel("DAH7PS Ib full-sequence bit score")
    ax_scatter.set_ylabel("KDOPS full-sequence bit score")
    ax_scatter.set_title("Ib HMM vs KDOPS decoy HMM", loc="left", fontweight="bold")
    ax_scatter.legend(frameon=False, loc="upper left")

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(color="#D9D9D9", linewidth=0.45, alpha=0.7)

    fig.suptitle(
        "KDOPS negative-selection filter for Type Ib candidates",
        y=0.99,
        fontsize=11,
        fontweight="bold",
    )
    return save_figure(fig, stem)


def default_profiles() -> list[tuple[str, Path]]:
    return [
        ("Ia", REPO_ROOT / "results" / "01_mining" / "model_Ia.hmm"),
        ("Ib", REPO_ROOT / "results" / "01_mining" / "model_Ib.hmm"),
        ("II", REPO_ROOT / "results" / "01_mining" / "model_II.hmm"),
        ("KDOPS", REPO_ROOT / "results" / "01_mining" / "kdo8ps.hmm"),
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Render Phase 1 HMM profile logo, emission heatmap, and KDOPS scoring plots."
    )
    parser.add_argument(
        "--kdops_hmm",
        default="results/01_mining/kdo8ps.hmm",
        help="KDOPS decoy HMM used for the profile logo.",
    )
    parser.add_argument(
        "--ia_hmm",
        default="results/01_mining/model_Ia.hmm",
        help="Type Ia subtype HMM.",
    )
    parser.add_argument(
        "--ib_hmm",
        default="results/01_mining/model_Ib.hmm",
        help="Type Ib subtype HMM.",
    )
    parser.add_argument(
        "--ii_hmm",
        default="results/01_mining/model_II.hmm",
        help="Type II subtype HMM.",
    )
    parser.add_argument(
        "--kdops_report",
        default="results/01_mining/kdops_filter_report_v41.tsv",
        help="KDOPS competitive scoring report TSV.",
    )
    parser.add_argument(
        "--logo_stem",
        default="F21_phase1_kdops_profile_logo",
        help="Output stem for the KDOPS profile logo.",
    )
    parser.add_argument(
        "--heatmap_stem",
        default="F22_phase1_hmm_emission_heatmap",
        help="Output stem for the four-HMM emission heatmap.",
    )
    parser.add_argument(
        "--scoring_stem",
        default="F23_phase1_kdops_competitive_scoring",
        help="Output stem for the KDOPS competitive scoring figure.",
    )
    parser.add_argument(
        "--skip_scoring",
        action="store_true",
        help="Only render the HMM profile logo and heatmap; leave the scoring figure untouched.",
    )
    return parser


def main() -> int:
    set_publication_style()
    args = build_parser().parse_args()

    paths = {
        "Ia": rel(Path(args.ia_hmm)),
        "Ib": rel(Path(args.ib_hmm)),
        "II": rel(Path(args.ii_hmm)),
        "KDOPS": rel(Path(args.kdops_hmm)),
    }
    profiles = [parse_hmmer_profile(path, label) for label, path in paths.items()]
    kdops_profile = next(profile for profile in profiles if profile.label == "KDOPS")

    outputs = []
    outputs.append(render_kdops_profile_logo(kdops_profile, args.logo_stem))
    outputs.append(render_emission_heatmap(profiles, args.heatmap_stem))
    if not args.skip_scoring:
        outputs.append(
            render_kdops_competitive_scoring(rel(Path(args.kdops_report)), args.scoring_stem)
        )

    print("Rendered Phase 1 HMM figures:")
    for output_pair in outputs:
        print(f"  {output_pair['pdf'].relative_to(REPO_ROOT)}")
        print(f"  {output_pair['png'].relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
