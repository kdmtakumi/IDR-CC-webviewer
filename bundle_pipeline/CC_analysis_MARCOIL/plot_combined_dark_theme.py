#!/usr/bin/env python3
"""
Create a dark-theme combined plot (Disorder + Coiled-coil) from a plot_data CSV.

Inputs:
  - CSV with columns:
      Position, Disorder_Score_Original, Disorder_Score_3res_MA,
      CC_Probability_Original, CC_Probability_3res_MA
Outputs:
  - <output_prefix>_dark.png                (normal dark theme)
  - <output_prefix>_dark_threshold.png      (threshold-applied, dark blend)

Features:
  - Black background / white text
  - Legends and layout aligned to the standard combined plot
  - X ticks configurable (default 200)
  - Threshold processing darkens (<50%) regions toward black instead of white
"""

import argparse
import csv
import io
import os
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from PIL import Image


# Palette / style constants
COLOR_DISORDER = "#ff595e"
COLOR_CC = "#1982c4"
COLOR_THRESH = "#ffca3a"
COLOR_GRID = "#666666"
COLOR_TEXT = "white"
COLOR_BG = "black"


def ensure_mplconfig():
    """Use a local matplotlib config dir to avoid permission issues."""
    local_cfg = Path(__file__).resolve().parent / ".mplconfig"
    local_cfg.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(local_cfg))


def load_plot_data(csv_path: Path):
    positions = []
    disorder_orig = []
    disorder_ma = []
    cc_orig = []
    cc_ma = []

    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            positions.append(int(row["Position"]))
            disorder_orig.append(float(row["Disorder_Score_Original"]))
            disorder_ma.append(float(row["Disorder_Score_3res_MA"]))
            cc_orig.append(float(row["CC_Probability_Original"]))
            cc_ma.append(float(row["CC_Probability_3res_MA"]))

    return positions, disorder_orig, disorder_ma, cc_orig, cc_ma


def style_axes(ax, twin, add_grid=True):
    ax.set_facecolor(COLOR_BG)
    twin.set_facecolor("none")
    ax.set_ylim(-5, 105)
    twin.set_ylim(-5, 105)

    for axis in (ax, twin):
        axis.tick_params(axis="y", colors=COLOR_TEXT, labelsize=13, width=1, length=4)
        for label in axis.get_yticklabels():
            label.set_fontweight("bold")
        for spine in axis.spines.values():
            spine.set_color(COLOR_TEXT)
            spine.set_linewidth(1)

    ax.axhline(0, color=COLOR_GRID, linewidth=0.8, alpha=0.7)
    ax.grid(add_grid, alpha=0.3, linestyle=":")


def build_figure(
    positions: List[int],
    disorder_orig: List[float],
    disorder_ma: List[float],
    cc_orig: List[float],
    cc_ma: List[float],
    seq_name: str,
    xtick_step: int,
) -> Tuple[plt.Figure, List[Tuple[plt.Axes, plt.Axes]]]:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))
    fig.patch.set_facecolor(COLOR_BG)

    # Top panel
    ax1_twin = ax1.twinx()
    ax1.plot(
        positions,
        disorder_orig,
        linewidth=2.0,
        color=COLOR_DISORDER,
        alpha=0.8,
        label="Disorder Score",
    )
    ax1_twin.plot(
        positions,
        cc_orig,
        linewidth=2.0,
        color=COLOR_CC,
        alpha=0.8,
        label="Coiled-coil Probability",
    )
    ax1_twin.axhline(y=50, color=COLOR_THRESH, linestyle="--", linewidth=1.2, alpha=0.9)
    ax1_twin.axhline(y=90, color=COLOR_GRID, linestyle="--", linewidth=1.0, alpha=0.8)
    ax1.set_ylabel("Disorder Score (%)", fontsize=12, fontweight="bold", color=COLOR_TEXT)
    ax1_twin.set_ylabel(
        "Coiled-coil Probability (%)", fontsize=12, fontweight="bold", color=COLOR_TEXT
    )
    ax1.set_title(
        f"{seq_name} - Original Data (Disorder vs Coiled-coil)",
        fontsize=14,
        fontweight="bold",
        color=COLOR_TEXT,
    )
    ax1.set_xlim(min(positions), max(positions))
    style_axes(ax1, ax1_twin, add_grid=True)
    ax1.set_xticks(list(range(0, max(positions)+1, xtick_step)))

    legend_elements = [
        Line2D([0], [0], color=COLOR_DISORDER, linewidth=2, alpha=0.8, label="Disorder Score"),
        Line2D([0], [0], color=COLOR_CC, linewidth=2, alpha=0.8, label="Coiled-coil Probability"),
    ]
    leg1 = ax1.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=2,
        frameon=True,
        facecolor="white",
        edgecolor="black",
        fontsize=11,
    )
    for text in leg1.get_texts():
        text.set_color(COLOR_TEXT)
        text.set_fontweight("bold")
    ax1.set_xlabel("Position in Sequence", fontsize=12, fontweight="bold", color=COLOR_TEXT, labelpad=30)

    # Bottom panel
    ax2_twin = ax2.twinx()
    ax2.plot(
        positions,
        disorder_ma,
        linewidth=2.0,
        color=COLOR_DISORDER,
        alpha=0.8,
        label="Disorder Score (3-res MA)",
    )
    ax2_twin.plot(
        positions,
        cc_ma,
        linewidth=2.0,
        color=COLOR_CC,
        alpha=0.8,
        label="Coiled-coil Probability (3-res MA)",
    )
    ax2_twin.axhline(y=50, color=COLOR_THRESH, linestyle="--", linewidth=1.2, alpha=0.9)
    ax2_twin.axhline(y=90, color=COLOR_GRID, linestyle="--", linewidth=1.0, alpha=0.8)
    ax2.set_ylabel("Disorder Score (%)", fontsize=12, fontweight="bold", color=COLOR_TEXT)
    ax2_twin.set_ylabel(
        "Coiled-coil Probability (%)", fontsize=12, fontweight="bold", color=COLOR_TEXT
    )
    ax2.set_title(
        f"{seq_name} - 3-residue Moving Average (Disorder vs Coiled-coil)",
        fontsize=14,
        fontweight="bold",
        color=COLOR_TEXT,
    )
    ax2.set_xlim(min(positions), max(positions))
    style_axes(ax2, ax2_twin, add_grid=True)
    ax2.set_xticks(list(range(0, max(positions)+1, xtick_step)))

    legend_elements2 = [
        Line2D([0], [0], color=COLOR_DISORDER, linewidth=2, alpha=0.8, label="Disorder Score (3-res MA)"),
        Line2D([0], [0], color=COLOR_CC, linewidth=2, alpha=0.8, label="Coiled-coil Probability (3-res MA)"),
    ]
    leg2 = ax2.legend(
        handles=legend_elements2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=2,
        frameon=True,
        facecolor="white",
        edgecolor="black",
        fontsize=11,
    )
    for text in leg2.get_texts():
        text.set_color(COLOR_TEXT)
        text.set_fontweight("bold")
    ax2.set_xlabel("Position in Sequence", fontsize=12, fontweight="bold", color=COLOR_TEXT, labelpad=35)

    xticks = list(range(0, max(positions) + 1, xtick_step))
    for ax in (ax1, ax2):
        ax.set_xticks(xticks)
        ax.tick_params(axis="x", labelbottom=True)
    for lbl in ax1.get_xticklabels() + ax2.get_xticklabels():
        lbl.set_color(COLOR_TEXT)
        lbl.set_fontweight("bold")
        lbl.set_fontsize(13)

    plt.tight_layout(rect=(0, 0.05, 1, 0.98))
    return fig, [(ax1, ax1_twin), (ax2, ax2_twin)]


def darken_threshold(fig: plt.Figure, axes_targets, threshold=50.0, alpha_blend=0.2, dpi=300):
    """
    Darken pixels below the threshold line toward black for the given axes/line colors.
    axes_targets: list of (axis, target_rgb_0to1)
    """
    buf = io.BytesIO()
    fig.savefig(buf, dpi=dpi, facecolor=fig.get_facecolor())
    buf.seek(0)
    img = Image.open(buf)
    img_array = np.array(img).astype(float) / 255.0
    buf.close()

    height, width = img_array.shape[:2]

    for axis, target in axes_targets:
        bbox = axis.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        left_px = int(bbox.x0 * dpi)
        right_px = int(bbox.x1 * dpi)
        top_px = int((fig.get_figheight() - bbox.y1) * dpi)
        bottom_px = int((fig.get_figheight() - bbox.y0) * dpi)

        ylim = axis.get_ylim()
        thr_norm = (threshold - ylim[0]) / (ylim[1] - ylim[0])
        thr_px = top_px + (bottom_px - top_px) * (1 - thr_norm)

        start_x = max(0, left_px - 5)
        target_vec = np.array(target)

        for px_x in range(start_x, min(right_px + 1, width)):
            for px_y in range(max(0, top_px), min(bottom_px + 1, height)):
                if px_y >= thr_px:
                    pixel = img_array[px_y, px_x, :3]
                    if np.linalg.norm(pixel - target_vec) < 0.6:  # tolerance in 0-1 scale
                        img_array[px_y, px_x, :3] = pixel * alpha_blend  # blend toward black

    return Image.fromarray((img_array * 255).astype(np.uint8))


def parse_args():
    parser = argparse.ArgumentParser(description="Generate dark-theme combined plot with threshold darkening.")
    parser.add_argument("--csv", required=True, help="Plot data CSV (exported via export_plot_data.py)")
    parser.add_argument("--seq-name", default=None, help="Sequence name for titles (default: CSV stem)")
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Output prefix (default: <csv_stem>) -> generates <prefix>_dark.png and _dark_threshold.png",
    )
    parser.add_argument("--xtick-step", type=int, default=200, help="X tick interval (default: 200)")
    parser.add_argument("--alpha-blend", type=float, default=0.2, help="Blend ratio to darken (<50%%) (default: 0.2)")
    return parser.parse_args()


def main():
    ensure_mplconfig()
    args = parse_args()

    csv_path = Path(args.csv)
    seq_name = args.seq_name or csv_path.stem
    output_prefix = args.output_prefix or csv_path.stem

    positions, disorder_orig, disorder_ma, cc_orig, cc_ma = load_plot_data(csv_path)
    fig, axes_pairs = build_figure(
        positions,
        disorder_orig,
        disorder_ma,
        cc_orig,
        cc_ma,
        seq_name=seq_name,
        xtick_step=args.xtick_step,
    )

    normal_png = Path(f"{output_prefix}_dark.png")
    fig.savefig(normal_png, dpi=300, facecolor=fig.get_facecolor())
    print(f"Saved dark theme plot: {normal_png}")

    axes_targets = [
        (axes_pairs[0][0], np.array([1.0, 0.0, 0.0])),  # top disorder (red)
        (axes_pairs[0][1], np.array([0.0, 0.0, 1.0])),  # top CC (blue)
        (axes_pairs[1][0], np.array([1.0, 0.0, 0.0])),  # bottom disorder (red)
        (axes_pairs[1][1], np.array([0.0, 0.0, 1.0])),  # bottom CC (blue)
    ]
    threshold_img = darken_threshold(
        fig,
        axes_targets,
        threshold=50.0,
        alpha_blend=args.alpha_blend,
        dpi=300,
    )
    threshold_png = Path(f"{output_prefix}_dark_threshold.png")
    threshold_img.save(threshold_png, dpi=(300, 300))
    print(f"Saved dark theme threshold plot: {threshold_png}")

    plt.close(fig)


if __name__ == "__main__":
    main()
