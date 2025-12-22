#!/usr/bin/env python3
"""
Threshold-styled Disorder/CC plot without TM overlay.
Outputs light/dark PNGs and region tables (IDR/CC) at 50% by default.
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

COLOR_TOL = 230


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


def flags_to_intervals(flags: List[bool]) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []
    start = None
    for i, val in enumerate(flags, start=1):
        if val and start is None:
            start = i
        if start is not None and (not val or i == len(flags)):
            end = i if val and i == len(flags) else i - 1
            out.append((start, end))
            start = None
    return out


def build_base_fig(name, pos, d0, d3, c0, c3, dark, xtick_step):
    if dark:
        plt.style.use("dark_background")
        face = COLOR_BG
    else:
        plt.style.use("default")
        face = "white"
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
    xticks = list(range(0, int(pos.max()) + xtick_step, xtick_step))
    c_dis = COLOR_DISORDER if dark else "red"
    c_cc = COLOR_CC if dark else "blue"

    def style(ax, twin):
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

    # top
    ax1_twin = ax1.twinx()
    ax1.plot(pos, d0, lw=2, color=c_dis, alpha=0.8 if dark else 0.7)
    ax1_twin.plot(pos, c0, lw=2, color=c_cc, alpha=0.8 if dark else 0.7)
    ax1_twin.axhline(y=50, color=COLOR_THRESH if dark else "orange", ls="--", lw=1.2 if dark else 1, alpha=0.9 if dark else 0.7)
    ax1_twin.axhline(y=90, color=COLOR_GRID if dark else "gray", ls="--", lw=1.0, alpha=0.8 if dark else 0.7)
    ax1.set_ylabel("Disorder Score (%)", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    ax1_twin.set_ylabel("Coiled-coil Probability (%)", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    ax1.set_title(f"{name} - Original (threshold)", fontsize=14, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    style(ax1, ax1_twin)
    legend_top = [
        Line2D([0], [0], color=c_dis, lw=2, alpha=0.8 if dark else 0.7, label="Disorder"),
        Line2D([0], [0], color=c_cc, lw=2, alpha=0.8 if dark else 0.7, label="Coiled-coil"),
    ]
    ax1.legend(handles=legend_top, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2,
               framealpha=0.95, fontsize=11, fancybox=True,
               facecolor=("none" if dark else "white"), edgecolor=("white" if dark else "black"))
    ax1.set_xlabel("Position in Sequence", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"), labelpad=30)

    # bottom
    ax2_twin = ax2.twinx()
    ax2.plot(pos, d3, lw=2, color=c_dis, alpha=0.8 if dark else 0.7)
    ax2_twin.plot(pos, c3, lw=2, color=c_cc, alpha=0.8 if dark else 0.7)
    ax2_twin.axhline(y=50, color=COLOR_THRESH if dark else "orange", ls="--", lw=1.2 if dark else 1, alpha=0.9 if dark else 0.7)
    ax2_twin.axhline(y=90, color=COLOR_GRID if dark else "gray", ls="--", lw=1.0, alpha=0.8 if dark else 0.7)
    ax2.set_ylabel("Disorder Score (%)", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    ax2_twin.set_ylabel("Coiled-coil Probability (%)", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    ax2.set_xlabel("Position in Sequence", fontsize=12, fontweight="bold", color=(COLOR_TEXT if dark else "black"), labelpad=30)
    ax2.set_title(f"{name} - 3-residue MA (threshold)", fontsize=14, fontweight="bold", color=(COLOR_TEXT if dark else "black"))
    style(ax2, ax2_twin)
    legend_bottom = [
        Line2D([0], [0], color=c_dis, lw=2, alpha=0.8 if dark else 0.7, label="Disorder (3-res MA)"),
        Line2D([0], [0], color=c_cc, lw=2, alpha=0.8 if dark else 0.7, label="Coiled-coil (3-res MA)"),
    ]
    ax2.legend(handles=legend_bottom, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2,
               framealpha=0.95, fontsize=11, fancybox=True,
               facecolor=("none" if dark else "white"), edgecolor=("white" if dark else "black"))

    fig.tight_layout(rect=(0, 0.05, 1, 0.98))
    fig.patch.set_facecolor(face)
    return fig, (ax1, ax1_twin, ax2, ax2_twin), face


def apply_threshold_to_image(img_array, fig, ax, color_rgb, dpi, threshold, alpha_blend):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    left_px = int(bbox.x0 * dpi)
    right_px = int(bbox.x1 * dpi)
    top_px = int((fig.get_figheight() - bbox.y1) * dpi)
    bottom_px = int((fig.get_figheight() - bbox.y0) * dpi)
    ylim = ax.get_ylim()
    thr_norm = (threshold - ylim[0]) / (ylim[1] - ylim[0])
    thr_px = top_px + (bottom_px - top_px) * (1 - thr_norm)
    target = np.array(color_rgb) * 255
    start_x = max(0, left_px - 5)
    sub = img_array[max(0, int(top_px)) : int(bottom_px) + 1, start_x : int(right_px) + 1, :3]
    if sub.size == 0:
        return
    y_indices = np.arange(max(0, int(top_px)), int(bottom_px) + 1)
    mask_y = (y_indices[:, None] >= thr_px)
    diff = np.linalg.norm(sub - target, axis=2)
    mask_color = diff < COLOR_TOL
    mask = mask_y & mask_color
    sub[mask] = sub[mask] * alpha_blend + 255 * (1 - alpha_blend)
    img_array[max(0, int(top_px)) : int(bottom_px) + 1, start_x : int(right_px) + 1, :3] = sub


def apply_light_threshold(fig, axes, face, out_png, threshold, alpha_blend=0.2):
    dpi = 300
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor=face)
    buf.seek(0)
    img = np.array(Image.open(buf)).astype(float)
    buf.close()
    ax1, ax1_twin, ax2, ax2_twin = axes
    apply_threshold_to_image(img, fig, ax1, (1.0, 0.0, 0.0), dpi, threshold, alpha_blend)
    apply_threshold_to_image(img, fig, ax1_twin, (0.0, 0.0, 1.0), dpi, threshold, alpha_blend)
    apply_threshold_to_image(img, fig, ax2, (1.0, 0.0, 0.0), dpi, threshold, alpha_blend)
    apply_threshold_to_image(img, fig, ax2_twin, (0.0, 0.0, 1.0), dpi, threshold, alpha_blend)
    img = np.clip(img, 0, 255)
    Image.fromarray(img.astype(np.uint8)).save(out_png, dpi=(dpi, dpi))
    print(f"Saved {out_png}")


def apply_dark_threshold(fig, axes, out_png, threshold, alpha_blend=0.2):
    dpi = 300
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor=COLOR_BG)
    buf.seek(0)
    img = np.array(Image.open(buf)).astype(float)
    buf.close()
    ax1, ax1_twin, ax2, ax2_twin = axes

    def darken(ax, color_rgb):
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        left_px = int(bbox.x0 * dpi)
        right_px = int(bbox.x1 * dpi)
        top_px = int((fig.get_figheight() - bbox.y1) * dpi)
        bottom_px = int((fig.get_figheight() - bbox.y0) * dpi)
        ylim = ax.get_ylim()
        thr_norm = (threshold - ylim[0]) / (ylim[1] - ylim[0])
        thr_px = top_px + (bottom_px - top_px) * (1 - thr_norm)
        target = np.array(color_rgb) * 255
        start_x = max(0, left_px - 5)
        sub = img[max(0, int(top_px)) : int(bottom_px) + 1, start_x : int(right_px) + 1, :3]
        if sub.size == 0:
            return
        y_idx = np.arange(max(0, int(top_px)), int(bottom_px) + 1)
        mask_y = (y_idx[:, None] >= thr_px)
        diff = np.linalg.norm(sub - target, axis=2)
        mask_color = diff < COLOR_TOL
        mask = mask_y & mask_color
        sub[mask] = sub[mask] * alpha_blend
        img[max(0, int(top_px)) : int(bottom_px) + 1, start_x : int(right_px) + 1, :3] = sub

    darken(ax1, (1.0, 0.0, 0.0))
    darken(ax1_twin, (0.0, 0.0, 1.0))
    darken(ax2, (1.0, 0.0, 0.0))
    darken(ax2_twin, (0.0, 0.0, 1.0))
    img = np.clip(img, 0, 255)
    Image.fromarray(img.astype(np.uint8)).save(out_png, dpi=(dpi, dpi))
    print(f"Saved {out_png}")


def export_regions(out_prefix: Path, pos, residues, heptads, d0, d3, c0, c3, threshold):
    idr_flags = [v >= threshold for v in d0]
    cc_flags = [v >= threshold for v in c0]
    tm_flags = [False] * len(pos)

    pos_csv = out_prefix.with_name(out_prefix.name + "_regions_positions.csv")
    with pos_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "Position","Residue",
            "Disorder_Score_Original","Disorder_Score_3res_MA","In_IDR_>=50",
            "CC_Probability_Original","CC_Probability_3res_MA","In_CC_>=50",
            "Heptad_Phase","In_TM"
        ])
        for p, r, hp, d_orig, d_ma, c_orig, c_ma, fi, fc, ft in zip(pos, residues, heptads, d0, d3, c0, c3, idr_flags, cc_flags, tm_flags):
            w.writerow([p, r, f"{d_orig:.2f}", f"{d_ma:.2f}", "Yes" if fi else "No",
                        f"{c_orig:.2f}", f"{c_ma:.2f}", "Yes" if fc else "No",
                        hp, "Yes" if ft else "No"])

    interval_csv = out_prefix.with_name(out_prefix.name + "_regions_intervals.csv")
    with interval_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Type","Start","End","Length"])
        for t, flags in (("IDR", idr_flags), ("CC", cc_flags)):
            for s, e in flags_to_intervals(flags):
                w.writerow([t, s, e, e - s + 1])
    print(f"Saved region tables: {pos_csv}, {interval_csv}")


def main():
    ap = argparse.ArgumentParser(description="Plot thresholded Disorder/CC (no TM) and export region tables.")
    ap.add_argument("--name", required=True)
    ap.add_argument("--plot-csv", required=True, help="Plot data CSV from export_plot_data.py")
    ap.add_argument("--out-prefix", required=True, help="Output prefix for PNG and CSV")
    ap.add_argument("--xtick-step", type=int, default=200)
    ap.add_argument("--threshold", type=float, default=50.0)
    args = ap.parse_args()

    pos, residues, d0, d3, c0, c3, heptads = load_plot_data(Path(args.plot_csv))

    fig_l, axes_l, face_l = build_base_fig(args.name, pos, d0, d3, c0, c3, dark=False, xtick_step=args.xtick_step)
    apply_light_threshold(fig_l, axes_l, face_l, Path(f"{args.out_prefix}_white.png"), threshold=args.threshold)
    plt.close(fig_l)

    fig_d, axes_d, face_d = build_base_fig(args.name, pos, d0, d3, c0, c3, dark=True, xtick_step=args.xtick_step)
    apply_dark_threshold(fig_d, axes_d, Path(f"{args.out_prefix}_dark.png"), threshold=args.threshold)
    plt.close(fig_d)

    export_regions(Path(args.out_prefix), pos, residues, heptads, d0, d3, c0, c3, args.threshold)


if __name__ == "__main__":
    main()
