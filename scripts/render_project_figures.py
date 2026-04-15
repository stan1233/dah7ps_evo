#!/usr/bin/env python3
"""Render all DAH7PS project figures into figures/."""

from __future__ import annotations

from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from figures.project_figures_lib import FIGURE_SPECS, render_all_figures, write_caption_markdown


def main() -> int:
    outputs = render_all_figures(display=False)
    markdown_path = write_caption_markdown(outputs=outputs)
    for spec, path_map in zip(FIGURE_SPECS, outputs):
        print(spec.figure_id, path_map.get("pdf", ""), path_map.get("png", ""))
    print("CAPTIONS", markdown_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
