#!/usr/bin/env python3
"""Plot sequence length distribution histograms for four raw_full_*.fasta files."""

import matplotlib.pyplot as plt
from pathlib import Path


def parse_fasta_lengths(fasta_path: str) -> list[int]:
    """Parse a FASTA file and return a list of sequence lengths."""
    lengths = []
    current_len = 0
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line)
    if current_len > 0:
        lengths.append(current_len)
    return lengths


def main():
    base_dir = Path(__file__).parent
    files = [
        "raw_full_Ia.fasta",
        "raw_full_Ib.fasta",
        "raw_full_Ib_clean.fasta",
        "raw_full_II.fasta",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("Sequence Length Distribution", fontsize=16, fontweight="bold")

    colors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

    for ax, fname, color in zip(axes.flat, files, colors):
        path = base_dir / fname
        lengths = parse_fasta_lengths(str(path))

        label = fname.replace("raw_full_", "").replace(".fasta", "")
        ax.hist(lengths, bins=50, color=color, edgecolor="white", alpha=0.85)
        ax.set_title(f"{label}  (n={len(lengths):,})", fontsize=13)
        ax.set_xlabel("Sequence Length (aa)")
        ax.set_ylabel("Count")
        ax.set_xlim(0, 1000)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = base_dir / "length_distribution.png"
    plt.savefig(out_path, dpi=200)
    print(f"Saved to {out_path}")
    plt.show()


if __name__ == "__main__":
    main()
