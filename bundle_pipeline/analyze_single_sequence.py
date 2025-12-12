#!/usr/bin/env python3
"""
End-to-end IDR + coiled-coil analysis for a single sequence.
Usage example:
  python analyze_single_sequence.py --name TANGO1_KEGG --sequence "MAAAP..."
or
  python analyze_single_sequence.py --name TANGO1_KEGG --fasta path/to/input.fasta

Outputs are collected under <output_root>/<name>/ including:
  - FASTA, disorder plot/CSV/summary
  - MARCOIL domain CSV
  - Plot data CSV
  - Combined plots (light + threshold) and combined CSV
  - Dark theme plots (normal + threshold)
  - Replot variants (original, MA, 2-panel)
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

# Safer matplotlib cache configuration before importing pyplot
DEFAULT_CACHE_DIR = Path(".cache")
DEFAULT_MPL_DIR = Path(".matplotlib_cache")
os.environ.setdefault("XDG_CACHE_HOME", str(DEFAULT_CACHE_DIR))
os.environ.setdefault("MPLCONFIGDIR", str(DEFAULT_MPL_DIR))

import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import metapredict as meta  # type: ignore

from CC_analysis_MARCOIL import run_marcoil as rm
from CC_analysis_MARCOIL.export_plot_data import export_plot_data
from CC_analysis_MARCOIL.plot_combined_analysis import plot_combined_analysis
from CC_analysis_MARCOIL.plot_combined_dark_theme import (
    build_figure as build_dark_fig,
    darken_threshold,
    ensure_mplconfig,
    load_plot_data as load_dark_data,
)
from CC_analysis_MARCOIL.replot_from_csv import (
    load_plot_data as load_replot_data,
    plot_combined_2panel,
    plot_integrated_ma,
    plot_integrated_original,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run IDR + MARCOIL pipeline for one sequence.")
    parser.add_argument("--name", required=True, help="Sequence name (used for folder and titles)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sequence", help="Raw amino acid sequence (whitespace allowed)")
    group.add_argument("--fasta", help="Input FASTA file path")
    parser.add_argument(
        "--output-root",
        default=".",
        help="Root directory to place the per-sequence folder (default: current directory)",
    )
    parser.add_argument(
        "--marcoil-dir",
        default="CC_analysis_MARCOIL/MARCOIL",
        help="Path to MARCOIL directory containing the executable (default: CC_analysis_MARCOIL/MARCOIL)",
    )
    parser.add_argument(
        "--mode",
        default="H",
        choices=["H", "L"],
        help="MARCOIL mode: H (high sensitivity) or L (low) (default: H)",
    )
    return parser.parse_args()


def write_fasta(name: str, sequence: str, out_path: Path) -> None:
    sequence = "".join(sequence.split())
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(f">{name}\n{sequence}\n")


def run_disorder(name: str, sequence: str, out_dir: Path) -> Tuple[Path, Path, Path]:
    scores = meta.predict_disorder(sequence, version=3)
    domains = meta.predict_disorder_domains(sequence, version=3, disorder_threshold=0.5).disordered_domain_boundaries
    positions = np.arange(1, len(sequence) + 1)

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(positions, scores, linewidth=1.5, color="#2E86AB", label="Disorder score")
    ax.axhline(y=0.5, color="red", linestyle="--", linewidth=1.5, alpha=0.7, label="Threshold (0.5)")
    for i, (start, end) in enumerate(domains):
        ax.axvspan(start, end, alpha=0.2, color="orange", label="IDR" if i == 0 else "")
    ax.set_xlabel("Residue Position", fontsize=12, fontweight="bold")
    ax.set_ylabel("Disorder Score", fontsize=12, fontweight="bold")
    ax.set_title(f"{name} Protein Disorder Prediction (MetaPredict v3)", fontsize=14, fontweight="bold")
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(0, len(sequence))
    ax.grid(True, alpha=0.3, linestyle=":", linewidth=0.5)
    ax.legend(loc="upper right", fontsize=10)

    mean_disorder = float(np.mean(scores))
    idr_residues = sum((end - start + 1) for start, end in domains)
    idr_percentage = (idr_residues / len(sequence) * 100) if sequence else 0.0
    ax.text(
        0.02,
        0.98,
        f"Mean disorder: {mean_disorder:.3f}\nIDR coverage: {idr_percentage:.1f}%",
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    fig.tight_layout()
    plot_path = out_dir / f"{name}_disorder_plot.png"
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    csv_path = out_dir / f"{name}_disorder_scores.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Position", "Residue", "Disorder_Score", "In_IDR"])
        for pos, residue, score in zip(positions, sequence, scores):
            in_idr = any(start <= pos <= end for start, end in domains)
            writer.writerow([pos, residue, f"{score:.4f}", "Yes" if in_idr else "No"])

    summary_path = out_dir / f"{name}_disorder_summary.txt"
    with summary_path.open("w") as f:
        f.write(f"{name} Disorder Analysis Summary\n")
        f.write("=" * 50 + "\n")
        f.write(f"Sequence length: {len(sequence)} residues\n")
        f.write(f"Mean disorder score: {mean_disorder:.4f}\n")
        f.write(f"Max disorder score: {float(np.max(scores)):.4f}\n")
        f.write(f"Min disorder score: {float(np.min(scores)):.4f}\n")
        f.write(f"Number of IDRs: {len(domains)}\n")
        f.write(f"IDR coverage: {idr_percentage:.2f}%\n")
        f.write(f"Total IDR residues: {idr_residues}\n")
        if domains:
            f.write("IDR boundaries:\n")
            for i, (start, end) in enumerate(domains, 1):
                f.write(f"  IDR {i}: {start}-{end} (length {end-start+1})\n")

    return plot_path, csv_path, summary_path


def run_marcoil_pipeline(
    name: str, fasta_path: Path, out_dir: Path, marcoil_dir: Path, mode: str
) -> Tuple[Path, Path, Path]:
    outputs_dir = rm.run_marcoil(marcoil_dir, fasta_path, mode=mode)
    if outputs_dir is None:
        sys.exit("MARCOIL failed; aborting.")

    domains_file = outputs_dir / "Domains"
    problist_file = outputs_dir / "ProbList"
    domains_data = rm.parse_domains_file(domains_file)
    domain_csv = out_dir / f"{name}_coiled_coil_results.csv"
    if domains_data:
        rm.write_csv_output(domains_data, domain_csv)
    else:
        domain_csv.write_text("sequence_name,threshold,domain_number,start,end,length,max_probability\n")
    return problist_file, domains_file, domain_csv


def run_combined_plots(
    name: str,
    problist_file: Path,
    domains_file: Path,
    disorder_csv: Path,
    plot_data_csv: Path,
    out_dir: Path,
) -> Tuple[Path, Path, Path, Path]:
    # Plot data export (aligned disorder + CC with MA)
    export_plot_data(problist_file, disorder_csv, plot_data_csv, name)

    # Light theme combined + threshold + combined CSV
    combined_png = out_dir / f"{name}_combined_analysis.png"
    plot_combined_analysis(
        problist_file,
        disorder_csv,
        domains_file,
        str(combined_png),
        name,
        create_threshold_version=True,
        alpha_blend=0.2,
        export_csv=True,
        csv_output_dir=out_dir,
    )
    combined_threshold_png = out_dir / f"{name}_combined_analysis_threshold.png"
    combined_csv = out_dir / f"{name}_combined_scores.csv"

    # Dark theme combined + threshold
    ensure_mplconfig()
    positions, disorder_orig, disorder_ma, cc_orig, cc_ma = load_dark_data(plot_data_csv)
    fig, axes_pairs = build_dark_fig(
        positions,
        disorder_orig,
        disorder_ma,
        cc_orig,
        cc_ma,
        seq_name=name,
        xtick_step=200,
    )
    dark_png = out_dir / f"{name}_dark.png"
    fig.savefig(dark_png, dpi=300, facecolor=fig.get_facecolor())

    axes_targets = [
        (axes_pairs[0][0], np.array([1.0, 0.0, 0.0])),  # top disorder (red)
        (axes_pairs[0][1], np.array([0.0, 0.0, 1.0])),  # top CC (blue)
        (axes_pairs[1][0], np.array([1.0, 0.0, 0.0])),  # bottom disorder (red)
        (axes_pairs[1][1], np.array([0.0, 0.0, 1.0])),  # bottom CC (blue)
    ]
    dark_thresh_img = darken_threshold(
        fig,
        axes_targets,
        threshold=50.0,
        alpha_blend=0.2,
        dpi=300,
    )
    dark_threshold_png = out_dir / f"{name}_dark_threshold.png"
    dark_thresh_img.save(dark_threshold_png, dpi=(300, 300))
    plt.close(fig)

    return combined_png, combined_threshold_png, dark_png, dark_threshold_png


def run_replots(plot_data_csv: Path, out_prefix: Path, title: str) -> List[Path]:
    data = load_replot_data(plot_data_csv)
    original_png = Path(f"{out_prefix}_original.png")
    ma_png = Path(f"{out_prefix}_ma.png")
    panel_png = Path(f"{out_prefix}_2panel.png")
    plot_integrated_original(data, original_png, title)
    plot_integrated_ma(data, ma_png, title)
    plot_combined_2panel(data, panel_png, title)
    return [original_png, ma_png, panel_png]


def main() -> None:
    args = parse_args()
    name = args.name
    output_root = Path(args.output_root).resolve()
    out_dir = output_root / name
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prepare FASTA
    fasta_path = out_dir / f"{name}.fasta"
    if args.sequence:
        write_fasta(name, args.sequence, fasta_path)
    else:
        src = Path(args.fasta)
        if not src.exists():
            sys.exit(f"FASTA not found: {src}")
        seq_lines = [line.strip() for line in src.read_text().splitlines() if not line.startswith(">")]
        write_fasta(name, "".join(seq_lines), fasta_path)

    sequence = "".join(line.strip() for line in fasta_path.read_text().splitlines() if not line.startswith(">"))
    print(f"[INFO] Sequence length: {len(sequence)}")
    print(f"[INFO] Output directory: {out_dir}")

    # Disorder
    print("[STEP] Running disorder prediction...")
    disorder_plot, disorder_csv, disorder_summary = run_disorder(name, sequence, out_dir)
    print(f"  - Disorder plot: {disorder_plot}")
    print(f"  - Disorder scores: {disorder_csv}")
    print(f"  - Disorder summary: {disorder_summary}")

    # MARCOIL + domains CSV
    print("[STEP] Running MARCOIL...")
    problist_file, domains_file, domain_csv = run_marcoil_pipeline(
        name, fasta_path, out_dir, Path(args.marcoil_dir).resolve(), args.mode
    )
    print(f"  - Domains CSV: {domain_csv}")

    # Combined outputs (light + dark)
    print("[STEP] Building combined plots...")
    plot_data_csv = out_dir / f"{name}_plot_data.csv"
    combined_png, combined_threshold_png, dark_png, dark_threshold_png = run_combined_plots(
        name, problist_file, domains_file, disorder_csv, plot_data_csv, out_dir
    )
    combined_csv = out_dir / f"{name}_combined_scores.csv"
    print(f"  - Plot data CSV: {plot_data_csv}")
    print(f"  - Combined (light): {combined_png}")
    print(f"  - Combined threshold (light): {combined_threshold_png}")
    print(f"  - Combined CSV: {combined_csv}")
    print(f"  - Combined (dark): {dark_png}")
    print(f"  - Combined threshold (dark): {dark_threshold_png}")

    # Replots
    print("[STEP] Generating replot variants...")
    replot_prefix = out_dir / f"{name}_plot_data_replot"
    replot_paths = run_replots(plot_data_csv, replot_prefix, f"{name} analysis")
    for p in replot_paths:
        print(f"  - {p}")

    print("\n✓ Done.")


if __name__ == "__main__":
    main()
