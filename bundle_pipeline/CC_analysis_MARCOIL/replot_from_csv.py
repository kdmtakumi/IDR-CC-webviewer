#!/usr/bin/env python3
"""
Re-plot graphs from exported CSV data
This script demonstrates how to create plots from the CSV data exported by export_plot_data.py
"""

import csv
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path


def load_plot_data(csv_file):
    """Load plot data from CSV"""
    positions = []
    residues = []
    disorder_orig = []
    disorder_ma = []
    cc_orig = []
    cc_ma = []
    heptad_phase = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            positions.append(int(row['Position']))
            residues.append(row['Residue'])
            disorder_orig.append(float(row['Disorder_Score_Original']))
            disorder_ma.append(float(row['Disorder_Score_3res_MA']))
            cc_orig.append(float(row['CC_Probability_Original']))
            cc_ma.append(float(row['CC_Probability_3res_MA']))
            heptad_phase.append(row['Heptad_Phase'])

    return {
        'positions': positions,
        'residues': residues,
        'disorder_orig': disorder_orig,
        'disorder_ma': disorder_ma,
        'cc_orig': cc_orig,
        'cc_ma': cc_ma,
        'heptad_phase': heptad_phase
    }


def plot_integrated_original(data, output_file, title='Disorder vs Coiled-coil'):
    """
    Create integrated plot with ORIGINAL data (no moving average)

    Args:
        data: Dictionary containing plot data
        output_file: Output PNG file path
        title: Plot title
    """
    fig, ax = plt.subplots(1, 1, figsize=(16, 6))
    ax_twin = ax.twinx()

    # Plot disorder score on left axis (red)
    ax.plot(data['positions'], data['disorder_orig'], linewidth=2.0, color='red',
            alpha=0.7, label='Disorder Score')
    ax.set_ylabel('Disorder Score (%)', fontsize=12, fontweight='bold', color='black')
    ax.tick_params(axis='y', labelcolor='black')
    ax.set_ylim(-5, 105)

    # Plot coiled-coil on right axis (blue)
    ax_twin.plot(data['positions'], data['cc_orig'], linewidth=2.0, color='blue',
                 alpha=0.7, label='Coiled-coil Probability')
    ax_twin.axhline(y=50, color='orange', linestyle='--', linewidth=1, alpha=0.7)
    ax_twin.axhline(y=90, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax_twin.set_ylabel('Coiled-coil Probability (%)', fontsize=12, fontweight='bold', color='black')
    ax_twin.tick_params(axis='y', labelcolor='black')
    ax_twin.set_ylim(-5, 105)

    ax.set_title(f'{title} - Original Data', fontsize=14, fontweight='bold')
    ax.set_xlim(min(data['positions']), max(data['positions']))
    ax.grid(True, alpha=0.3, linestyle=':')

    # Combined legend
    legend_elements = [
        Line2D([0], [0], color='red', linewidth=2, alpha=0.7, label='Disorder Score'),
        Line2D([0], [0], color='blue', linewidth=2, alpha=0.7, label='Coiled-coil Probability')
    ]
    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05),
              ncol=2, framealpha=0.95, fontsize=11, fancybox=True)
    ax.set_xlabel('Position in Sequence', fontsize=12, fontweight='bold', labelpad=35)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_integrated_ma(data, output_file, title='Disorder vs Coiled-coil'):
    """
    Create integrated plot with 3-RESIDUE MOVING AVERAGE data

    Args:
        data: Dictionary containing plot data
        output_file: Output PNG file path
        title: Plot title
    """
    fig, ax = plt.subplots(1, 1, figsize=(16, 6))
    ax_twin = ax.twinx()

    # Plot disorder score MA on left axis (red)
    ax.plot(data['positions'], data['disorder_ma'], linewidth=2.0, color='red',
            alpha=0.7, label='Disorder Score (3-res MA)')
    ax.set_ylabel('Disorder Score (%)', fontsize=12, fontweight='bold', color='black')
    ax.tick_params(axis='y', labelcolor='black')
    ax.set_ylim(-5, 105)

    # Plot coiled-coil MA on right axis (blue)
    ax_twin.plot(data['positions'], data['cc_ma'], linewidth=2.0, color='blue',
                 alpha=0.7, label='Coiled-coil Probability (3-res MA)')
    ax_twin.axhline(y=50, color='orange', linestyle='--', linewidth=1, alpha=0.7)
    ax_twin.axhline(y=90, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax_twin.set_ylabel('Coiled-coil Probability (%)', fontsize=12, fontweight='bold', color='black')
    ax_twin.tick_params(axis='y', labelcolor='black')
    ax_twin.set_ylim(-5, 105)

    ax.set_title(f'{title} - 3-residue Moving Average', fontsize=14, fontweight='bold')
    ax.set_xlim(min(data['positions']), max(data['positions']))
    ax.grid(True, alpha=0.3, linestyle=':')

    # Combined legend
    legend_elements = [
        Line2D([0], [0], color='red', linewidth=2, alpha=0.7, label='Disorder Score (3-res MA)'),
        Line2D([0], [0], color='blue', linewidth=2, alpha=0.7, label='Coiled-coil Probability (3-res MA)')
    ]
    ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05),
              ncol=2, framealpha=0.95, fontsize=11, fancybox=True)
    ax.set_xlabel('Position in Sequence', fontsize=12, fontweight='bold', labelpad=35)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_combined_2panel(data, output_file, title='Disorder vs Coiled-coil'):
    """
    Create 2-panel plot: Original (top) and 3-res MA (bottom)

    Args:
        data: Dictionary containing plot data
        output_file: Output PNG file path
        title: Plot title
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))

    # === TOP PANEL: Original data ===
    ax1_twin = ax1.twinx()

    # Disorder (red, left axis)
    ax1.plot(data['positions'], data['disorder_orig'], linewidth=2.0, color='red',
             alpha=0.7, label='Disorder Score')
    ax1.set_ylabel('Disorder Score (%)', fontsize=12, fontweight='bold', color='black')
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.set_ylim(-5, 105)

    # Coiled-coil (blue, right axis)
    ax1_twin.plot(data['positions'], data['cc_orig'], linewidth=2.0, color='blue',
                  alpha=0.7, label='Coiled-coil Probability')
    ax1_twin.axhline(y=50, color='orange', linestyle='--', linewidth=1, alpha=0.7)
    ax1_twin.axhline(y=90, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax1_twin.set_ylabel('Coiled-coil Probability (%)', fontsize=12, fontweight='bold', color='black')
    ax1_twin.tick_params(axis='y', labelcolor='black')
    ax1_twin.set_ylim(-5, 105)

    ax1.set_title(f'{title} - Original Data', fontsize=14, fontweight='bold')
    ax1.set_xlim(min(data['positions']), max(data['positions']))
    ax1.grid(True, alpha=0.3, linestyle=':')

    # Legend for top panel
    legend_elements = [
        Line2D([0], [0], color='red', linewidth=2, alpha=0.7, label='Disorder Score'),
        Line2D([0], [0], color='blue', linewidth=2, alpha=0.7, label='Coiled-coil Probability')
    ]
    ax1.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05),
               ncol=2, framealpha=0.95, fontsize=11, fancybox=True)
    ax1.set_xlabel('Position in Sequence', fontsize=12, fontweight='bold', labelpad=35)

    # === BOTTOM PANEL: 3-residue moving average ===
    ax2_twin = ax2.twinx()

    # Disorder MA (red, left axis)
    ax2.plot(data['positions'], data['disorder_ma'], linewidth=2.0, color='red',
             alpha=0.7, label='Disorder Score (3-res MA)')
    ax2.set_ylabel('Disorder Score (%)', fontsize=12, fontweight='bold', color='black')
    ax2.tick_params(axis='y', labelcolor='black')
    ax2.set_ylim(-5, 105)

    # Coiled-coil MA (blue, right axis)
    ax2_twin.plot(data['positions'], data['cc_ma'], linewidth=2.0, color='blue',
                  alpha=0.7, label='Coiled-coil Probability (3-res MA)')
    ax2_twin.axhline(y=50, color='orange', linestyle='--', linewidth=1, alpha=0.7)
    ax2_twin.axhline(y=90, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax2_twin.set_ylabel('Coiled-coil Probability (%)', fontsize=12, fontweight='bold', color='black')
    ax2_twin.tick_params(axis='y', labelcolor='black')
    ax2_twin.set_ylim(-5, 105)

    ax2.set_title(f'{title} - 3-residue Moving Average', fontsize=14, fontweight='bold')
    ax2.set_xlim(min(data['positions']), max(data['positions']))
    ax2.grid(True, alpha=0.3, linestyle=':')

    # Legend for bottom panel
    legend_elements = [
        Line2D([0], [0], color='red', linewidth=2, alpha=0.7, label='Disorder Score (3-res MA)'),
        Line2D([0], [0], color='blue', linewidth=2, alpha=0.7, label='Coiled-coil Probability (3-res MA)')
    ]
    ax2.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05),
               ncol=2, framealpha=0.95, fontsize=11, fancybox=True)
    ax2.set_xlabel('Position in Sequence', fontsize=12, fontweight='bold', labelpad=35)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 replot_from_csv.py <csv_file> [output_prefix] [title]")
        print("\nExample:")
        print("  python3 replot_from_csv.py ELKS2_plot_data.csv ELKS2_replot 'ELKS2 Analysis'")
        print("\nThis will create:")
        print("  - {prefix}_original.png       (original data)")
        print("  - {prefix}_ma.png            (3-residue moving average)")
        print("  - {prefix}_2panel.png        (both in 2-panel layout)")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else Path(csv_file).stem + "_replot"
    title = sys.argv[3] if len(sys.argv) > 3 else Path(csv_file).stem

    if not Path(csv_file).exists():
        print(f"ERROR: CSV file not found: {csv_file}")
        sys.exit(1)

    # Load data
    print(f"Loading data from {csv_file}...")
    data = load_plot_data(csv_file)
    print(f"  Loaded {len(data['positions'])} data points")

    # Create plots
    print("\nCreating plots...")
    plot_integrated_original(data, f"{output_prefix}_original.png", title)
    plot_integrated_ma(data, f"{output_prefix}_ma.png", title)
    plot_combined_2panel(data, f"{output_prefix}_2panel.png", title)

    print("\n✓ Done! Created 3 plot files:")
    print(f"  - {output_prefix}_original.png")
    print(f"  - {output_prefix}_ma.png")
    print(f"  - {output_prefix}_2panel.png")
    print("\nYou can modify this script to customize colors, sizes, labels, etc.")


if __name__ == "__main__":
    main()
