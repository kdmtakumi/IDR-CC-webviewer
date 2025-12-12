#!/usr/bin/env python3
"""
Overlay DeepTMHMM-predicted TM helices on Disorder/CC plots with threshold styling.

Inputs:
  - Plot data CSV (exported via export_plot_data.py)
  - DeepTMHMM GFF3 (prediction.gff3) containing TMhelix features
Outputs:
  - <out_prefix>_white.png (light theme, <50% lightened toward white)
  - <out_prefix>_dark.png  (dark theme, <50% darkened toward black)

Design matches combined_dark_threshold (legends, xtick spacing).
"""

from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path
from typing import List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from PIL import Image

from CC_analysis_MARCOIL.plot_combined_dark_theme import (
    COLOR_BG,
    COLOR_CC,
    COLOR_DISORDER,
    COLOR_GRID,
    COLOR_TEXT,
    COLOR_THRESH,
)

COLOR_TOL = 230  # tolerance to catch antialiased line pixels


def load_plot_data(csv_path: Path):
    pos: List[int] = []
    residues: List[str] = []
    heptads: List[str] = []
    d0: List[float] = []
    d3: List[float] = []
    c0: List[float] = []
    c3: List[float] = []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            pos.append(int(row["Position"]))
            residues.append(row.get("Residue", ""))
            heptads.append(row.get("Heptad_Phase", ""))
            d0.append(float(row["Disorder_Score_Original"]))
            d3.append(float(row["Disorder_Score_3res_MA"]))
            c0.append(float(row["CC_Probability_Original"]))
            c3.append(float(row["CC_Probability_3res_MA"]))
    return np.array(pos), residues, np.array(d0), np.array(d3), np.array(c0), np.array(c3), heptads


def parse_gff_tm_spans(gff_path: Path) -> List[Tuple[int, int]]:
    spans: List[Tuple[int, int]] = []
    tm_labels = {"TMhelix", "TMH", "TM", "TRANSMEM"}
    with gff_path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            feature = None
            start = None
            end = None
            if len(parts) >= 5:
                feature = parts[2]
                start = parts[3]
                end = parts[4]
            elif len(parts) >= 4:
                # Some DeepTMHMM outputs use seqid, feature, start, end only
                feature = parts[1]
                start = parts[2]
                end = parts[3]
            if feature is None or start is None or end is None:
                continue
            if feature not in tm_labels:
                continue
            spans.append((int(start), int(end)))
    return spans


def flags_to_intervals(flags: List[bool]) -> List[Tuple[int, int]]:
    intervals: List[Tuple[int, int]] = []
    if not flags:
        return intervals
    start = None
    for i, val in enumerate(flags, start=1):
        if val and start is None:
            start = i
        if start is not None and (not val or i == len(flags)):
            end = i if val and i == len(flags) else i - 1
            intervals.append((start, end))
            start = None
    return intervals


def add_tm_spans(ax, spans: List[Tuple[int, int]], dark: bool):
    color = "#8CFF66" if dark else "limegreen"  # brighter for dark background
    for s, e in spans:
        ax.axvspan(s, e, color=color, alpha=0.25)
    ax.set_ylim(-5, 105)


def build_base_fig(
    name: str,
    pos: np.ndarray,
    residues: List[str],
    d0: np.ndarray,
    d3: np.ndarray,
    c0: np.ndarray,
    c3: np.ndarray,
    heptads: List[str],
    tm_spans: List[Tuple[int, int]],
    dark: bool,
    xtick_step: int,
    threshold: float,
    out_prefix: Path,
):
    if dark:
        plt.style.use("dark_background")
        facecolor = COLOR_BG
    else:
        plt.style.use("default")
        facecolor = "white"

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

    xticks = list(range(0, int(pos.max()) + xtick_step, xtick_step))
    c_dis = COLOR_DISORDER if dark else "red"
    c_cc = COLOR_CC if dark else "blue"

    def style_axes(ax, twin):
        ax.set_ylim(-5, 105)
        twin.set_ylim(-5, 105)
        ax.set_xlim(pos.min(), pos.max())
        ax.grid(True, alpha=0.3, linestyle=":")
        ax.set_xticks(xticks)
        ax.tick_params(axis="x", labelbottom=True)
        if dark:
            ax.set_facecolor(COLOR_BG)
            twin.set_facecolor("none")
            for axis in (ax, twin):
                axis.tick_params(axis="both", colors=COLOR_TEXT, labelsize=13, width=1, length=4)
                for lbl in axis.get_yticklabels():
                    lbl.set_fontweight("bold")
                for spine in axis.spines.values():
                    spine.set_color(COLOR_TEXT)
                    spine.set_linewidth(1)
            for lbl in ax.get_xticklabels():
                lbl.set_color(COLOR_TEXT)
                lbl.set_fontweight("bold")
                lbl.set_fontsize(13)

    # Top panel
    ax1_twin = ax1.twinx()
    ax1.plot(pos, d0, linewidth=2.0, color=c_dis, alpha=0.8 if dark else 0.7)
    ax1_twin.plot(pos, c0, linewidth=2.0, color=c_cc, alpha=0.8 if dark else 0.7)
    ax1_twin.axhline(y=50, color=COLOR_THRESH if dark else "orange", ls="--", lw=1.2 if dark else 1, alpha=0.9 if dark else 0.7)
    ax1_twin.axhline(y=90, color=COLOR_GRID if dark else "gray", ls="--", lw=1.0, alpha=0.8 if dark else 0.7)
    add_tm_spans(ax1, tm_spans, dark)
    ax1.set_ylabel("Disorder Score (%)", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    ax1_twin.set_ylabel("Coiled-coil Probability (%)", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    ax1.set_title(f"{name} - Original (threshold + DeepTMHMM)", fontsize=14, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    style_axes(ax1, ax1_twin)
    legend_tm = Patch(facecolor="limegreen", edgecolor="none", alpha=0.25, label="TM helix (DeepTMHMM)")
    legend_top = [
        Line2D([0], [0], color=c_dis, lw=2, alpha=0.8 if dark else 0.7, label="Disorder"),
        Line2D([0], [0], color=c_cc, lw=2, alpha=0.8 if dark else 0.7, label="Coiled-coil"),
        legend_tm,
    ]
    leg1 = ax1.legend(
        handles=legend_top,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=3,
        framealpha=0.95,
        fontsize=11,
        fancybox=True,
        facecolor=("none" if dark else "white"),
        edgecolor=("white" if dark else "black"),
    )
    if dark:
        for txt in leg1.get_texts():
            txt.set_color(COLOR_TEXT)
            txt.set_fontweight("bold")

    # Bottom panel
    ax2_twin = ax2.twinx()
    ax2.plot(pos, d3, linewidth=2.0, color=c_dis, alpha=0.8 if dark else 0.7)
    ax2_twin.plot(pos, c3, linewidth=2.0, color=c_cc, alpha=0.8 if dark else 0.7)
    ax2_twin.axhline(y=50, color=COLOR_THRESH if dark else "orange", ls="--", lw=1.2 if dark else 1, alpha=0.9 if dark else 0.7)
    ax2_twin.axhline(y=90, color=COLOR_GRID if dark else "gray", ls="--", lw=1.0, alpha=0.8 if dark else 0.7)
    add_tm_spans(ax2, tm_spans, dark)
    ax2.set_ylabel("Disorder Score (%)", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    ax2_twin.set_ylabel("Coiled-coil Probability (%)", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    ax2.set_xlabel("Position in Sequence", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"), labelpad=30)
    ax2.set_title(f"{name} - 3-residue MA (threshold + DeepTMHMM)", fontsize=14, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    style_axes(ax2, ax2_twin)
    legend_bottom = [
        Line2D([0], [0], color=c_dis, lw=2, alpha=0.8 if dark else 0.7, label="Disorder (3-res MA)"),
        Line2D([0], [0], color=c_cc, lw=2, alpha=0.8 if dark else 0.7, label="Coiled-coil (3-res MA)"),
        legend_tm,
    ]
    leg2 = ax2.legend(
        handles=legend_bottom,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=3,
        framealpha=0.95,
        fontsize=11,
        fancybox=True,
        facecolor=("none" if dark else "white"),
        edgecolor=("white" if dark else "black"),
    )
    if dark:
        for txt in leg2.get_texts():
            txt.set_color(COLOR_TEXT)
            txt.set_fontweight("bold")

    fig.tight_layout(rect=(0, 0.05, 1, 0.98))
    fig.patch.set_facecolor(facecolor)
    return fig, (ax1, ax1_twin, ax2, ax2_twin), facecolor


def apply_threshold_to_image(img_array: np.ndarray, fig, ax, color_rgb, dpi=300, threshold=50.0, alpha_blend=0.2):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    left_px = int(bbox.x0 * dpi)
    right_px = int(bbox.x1 * dpi)
    top_px = int((fig.get_figheight() - bbox.y1) * dpi)
    bottom_px = int((fig.get_figheight() - bbox.y0) * dpi)
    ylim = ax.get_ylim()
    threshold_y_normalized = (threshold - ylim[0]) / (ylim[1] - ylim[0])
    threshold_y_pixel = top_px + (bottom_px - top_px) * (1 - threshold_y_normalized)

    target = np.array(color_rgb) * 255
    start_x = max(0, left_px - 5)

    sub = img_array[max(0, int(top_px)) : int(bottom_px) + 1, start_x : int(right_px) + 1, :3]
    if sub.size == 0:
        return
    y_indices = np.arange(max(0, int(top_px)), int(bottom_px) + 1)
    mask_y = (y_indices[:, None] >= threshold_y_pixel)
    diff = np.linalg.norm(sub - target, axis=2)
    mask_color = diff < COLOR_TOL
    mask = mask_y & mask_color
    sub[mask] = sub[mask] * alpha_blend + 255 * (1 - alpha_blend)
    img_array[max(0, int(top_px)) : int(bottom_px) + 1, start_x : int(right_px) + 1, :3] = sub


def apply_light_threshold(fig, axes, facecolor: str, outfile: Path, alpha_blend: float = 0.2, threshold: float = 50.0):
    dpi = 300
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor=facecolor)
    buf.seek(0)
    img_array = np.array(Image.open(buf)).astype(float)
    buf.close()

    ax1, ax1_twin, ax2, ax2_twin = axes
    apply_threshold_to_image(img_array, fig, ax1, (1.0, 0.0, 0.0), dpi=dpi, threshold=threshold, alpha_blend=alpha_blend)
    apply_threshold_to_image(img_array, fig, ax1_twin, (0.0, 0.0, 1.0), dpi=dpi, threshold=threshold, alpha_blend=alpha_blend)
    apply_threshold_to_image(img_array, fig, ax2, (1.0, 0.0, 0.0), dpi=dpi, threshold=threshold, alpha_blend=alpha_blend)
    apply_threshold_to_image(img_array, fig, ax2_twin, (0.0, 0.0, 1.0), dpi=dpi, threshold=threshold, alpha_blend=alpha_blend)

    img_array = np.clip(img_array, 0, 255)
    Image.fromarray(img_array.astype(np.uint8)).save(outfile, dpi=(dpi, dpi))
    print(f"Saved {outfile}")


def apply_dark_threshold(fig, axes, outfile: Path, alpha_blend: float = 0.2, threshold: float = 50.0):
    dpi = 300
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor=COLOR_BG)
    buf.seek(0)
    img_array = np.array(Image.open(buf)).astype(float)
    buf.close()

    def darken(ax, color_rgb):
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        left_px = int(bbox.x0 * dpi)
        right_px = int(bbox.x1 * dpi)
        top_px = int((fig.get_figheight() - bbox.y1) * dpi)
        bottom_px = int((fig.get_figheight() - bbox.y0) * dpi)
        ylim = ax.get_ylim()
        threshold_y_normalized = (threshold - ylim[0]) / (ylim[1] - ylim[0])
        threshold_y_pixel = top_px + (bottom_px - top_px) * (1 - threshold_y_normalized)

        target = np.array(color_rgb) * 255
        start_x = max(0, left_px - 5)
        sub = img_array[max(0, int(top_px)) : int(bottom_px) + 1, start_x : int(right_px) + 1, :3]
        if sub.size == 0:
            return
        y_indices = np.arange(max(0, int(top_px)), int(bottom_px) + 1)
        mask_y = (y_indices[:, None] >= threshold_y_pixel)
        diff = np.linalg.norm(sub - target, axis=2)
        mask_color = diff < COLOR_TOL
        mask = mask_y & mask_color
        sub[mask] = sub[mask] * alpha_blend
        img_array[max(0, int(top_px)) : int(bottom_px) + 1, start_x : int(right_px) + 1, :3] = sub

    ax1, ax1_twin, ax2, ax2_twin = axes
    darken(ax1, (1.0, 0.0, 0.0))
    darken(ax1_twin, (0.0, 0.0, 1.0))
    darken(ax2, (1.0, 0.0, 0.0))
    darken(ax2_twin, (0.0, 0.0, 1.0))

    img_array = np.clip(img_array, 0, 255)
    Image.fromarray(img_array.astype(np.uint8)).save(outfile, dpi=(dpi, dpi))
    print(f"Saved {outfile}")


def export_regions(
    out_prefix: Path,
    pos: List[int],
    residues: List[str],
    heptads: List[str],
    d0: np.ndarray,
    d3: np.ndarray,
    c0: np.ndarray,
    c3: np.ndarray,
    tm_spans: List[Tuple[int, int]],
    threshold: float,
):
    idr_flags = [v >= threshold for v in d0]
    cc_flags = [v >= threshold for v in c0]
    tm_flags = [any(s <= p <= e for s, e in tm_spans) for p in pos] if tm_spans else [False] * len(pos)

    pos_csv = out_prefix.with_name(out_prefix.name + "_regions_positions.csv")
    with pos_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Position",
            "Residue",
            "Disorder_Score_Original",
            "Disorder_Score_3res_MA",
            "In_IDR_>=50",
            "CC_Probability_Original",
            "CC_Probability_3res_MA",
            "In_CC_>=50",
            "Heptad_Phase",
            "In_TM",
        ])
        for p, r, hp, d_orig, d_ma, c_orig, c_ma, fi, fc, ft in zip(pos, residues, heptads, d0, d3, c0, c3, idr_flags, cc_flags, tm_flags):
            writer.writerow([
                p,
                r,
                f"{d_orig:.2f}",
                f"{d_ma:.2f}",
                "Yes" if fi else "No",
                f"{c_orig:.2f}",
                f"{c_ma:.2f}",
                "Yes" if fc else "No",
                hp,
                "Yes" if ft else "No",
            ])

    interval_csv = out_prefix.with_name(out_prefix.name + "_regions_intervals.csv")
    with interval_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Type", "Start", "End", "Length"])
        for t, flags in (("IDR", idr_flags), ("CC", cc_flags)):
            for s, e in flags_to_intervals(flags):
                writer.writerow([t, s, e, e - s + 1])
        for s, e in tm_spans:
            writer.writerow(["TM", s, e, e - s + 1])

    print(f"Saved region tables: {pos_csv}, {interval_csv}")


def main():
    parser = argparse.ArgumentParser(description="Plot DeepTMHMM overlay with threshold styling (pixel-based).")
    parser.add_argument("--name", required=True, help="Sequence name (titles/legends)")
    parser.add_argument("--plot-csv", required=True, help="Plot data CSV from export_plot_data.py")
    parser.add_argument("--gff", required=False, help="DeepTMHMM prediction.gff3 (optional)")
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix (writes <prefix>_white.png and <prefix>_dark.png)",
    )
    parser.add_argument("--xtick-step", type=int, default=200, help="X tick interval (default 200)")
    parser.add_argument("--threshold", type=float, default=50.0, help="Threshold for IDR/CC flags (default 50%)")
    args = parser.parse_args()

    pos, residues, d0, d3, c0, c3, heptads = load_plot_data(Path(args.plot_csv))
    spans = parse_gff_tm_spans(Path(args.gff)) if args.gff else []

    fig_l, axes_l, face_l = build_base_fig(args.name, pos, residues, d0, d3, c0, c3, heptads, spans, dark=False, xtick_step=args.xtick_step, threshold=args.threshold, out_prefix=Path(args.out_prefix))
    apply_light_threshold(fig_l, axes_l, face_l, Path(f"{args.out_prefix}_white.png"), threshold=args.threshold)
    plt.close(fig_l)

    fig_d, axes_d, face_d = build_base_fig(args.name, pos, residues, d0, d3, c0, c3, heptads, spans, dark=True, xtick_step=args.xtick_step, threshold=args.threshold, out_prefix=Path(args.out_prefix))
    apply_dark_threshold(fig_d, axes_d, Path(f"{args.out_prefix}_dark.png"), threshold=args.threshold)
    plt.close(fig_d)

    export_regions(Path(args.out_prefix), list(pos), residues, heptads, d0, d3, c0, c3, spans, args.threshold)


if __name__ == "__main__":
    main()
