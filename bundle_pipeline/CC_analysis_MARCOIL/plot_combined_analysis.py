#!/usr/bin/env python3
"""
Plot coiled-coil probability with disorder scores on dual y-axes
Generates both normal and threshold-based visualizations
"""

import re
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import csv
from PIL import Image
import io


def parse_problist(problist_file):
    """Parse MARCOIL ProbList output file"""
    results = {}
    current_seq = None

    with open(problist_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                seq_info = line[1:].strip()
                if '##' in seq_info:
                    current_seq = seq_info.split('##')[0].strip()
                else:
                    current_seq = seq_info
                results[current_seq] = []
                continue

            match = re.match(r'\s*(\d+)\s+([A-Z])\s+([\d.]+)\s+([a-g])', line)
            if match and current_seq:
                pos, aa, prob, phase = match.groups()
                results[current_seq].append({
                    'position': int(pos),
                    'amino_acid': aa,
                    'probability': float(prob),
                    'heptad_phase': phase
                })

    return results


def moving_average(data, window=3):
    """
    Calculate a centered moving average with edge padding so plot outputs stay
    aligned with the CSV export routine (3-residue MA with edge extension).
    """
    if not data:
        return []

    if window <= 0:
        raise ValueError("window size must be positive")

    if window % 2 == 0:
        raise ValueError("centered moving average requires an odd window size")

    if len(data) < window:
        avg = sum(data) / len(data)
        return [avg] * len(data)

    half_window = window // 2
    valid_averages = [
        sum(data[i:i + window]) / window
        for i in range(len(data) - window + 1)
    ]

    result = []
    result.extend([valid_averages[0]] * half_window)
    result.extend(valid_averages)
    result.extend([valid_averages[-1]] * half_window)
    return result


def parse_disorder_scores(disorder_csv):
    """Parse disorder scores CSV file"""
    positions = []
    scores = []

    with open(disorder_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            positions.append(int(row['Position']))
            # Convert to percentage if needed (0-1 range)
            score = float(row['Disorder_Score'])
            if score <= 1.0:
                score = score * 100
            scores.append(score)

    return positions, scores


def parse_domains_for_plot(domains_file):
    """Parse domains file to get domain regions"""
    domains = []
    current_seq = None

    with open(domains_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                seq_info = line[1:].strip()
                if '##' in seq_info:
                    current_seq = seq_info.split('##')[0].strip()
                else:
                    current_seq = seq_info
                continue

            if 'THRESHOLD 50.0' in line:
                threshold = 50.0
                continue

            domain_match = re.match(r'\s*(\d+)\.\s+from\s+(\d+)\s+to\s+(\d+)', line)
            if domain_match and current_seq:
                domain_num, start, end = domain_match.groups()
                domains.append({
                    'sequence': current_seq,
                    'start': int(start),
                    'end': int(end),
                    'threshold': 50.0
                })

    return domains


COLOR_TOL = 230  # tolerance for matching line colors during threshold lightening


def apply_threshold_lightening(fig, ax, color_rgb, threshold=50.0, alpha_blend=0.2):
    """
    Apply Y-coordinate based threshold lightening
    Pixels below the threshold line are lightened

    Args:
        fig: matplotlib figure
        ax: matplotlib axis
        color_rgb: target color as (R, G, B) in 0-1 scale
        threshold: threshold value (default 50.0 for 50%)
        alpha_blend: blend ratio for lightening (default 0.2 = 20% original, 80% white)

    Returns:
        Modified image as PIL Image
    """
    # Render the figure to a numpy array
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    img = Image.open(buf)
    img_array = np.array(img).astype(float)

    # Get the axes position in figure coordinates
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

    # Convert to pixel coordinates (at 300 dpi)
    dpi = 300
    left_px = int(bbox.x0 * dpi)
    right_px = int(bbox.x1 * dpi)
    top_px = int((fig.get_figheight() - bbox.y1) * dpi)
    bottom_px = int((fig.get_figheight() - bbox.y0) * dpi)

    # Get y-axis limits
    ylim = ax.get_ylim()

    # Calculate the pixel Y-coordinate of the threshold line
    threshold_y_normalized = (threshold - ylim[0]) / (ylim[1] - ylim[0])
    threshold_y_pixel = top_px + (bottom_px - top_px) * (1 - threshold_y_normalized)

    # Process pixels in the plot area
    height, width = img_array.shape[:2]
    pixels_modified = 0

    # Extend left_px slightly to catch edge pixels
    start_x = max(0, left_px - 5)

    for px_x in range(start_x, min(right_px + 1, width)):
        for px_y in range(max(0, top_px), min(bottom_px + 1, height)):
            # Check if this pixel is below the threshold line (>= to include boundary)
            if px_y >= threshold_y_pixel:
                pixel = img_array[px_y, px_x, :3]

                # Check if this pixel is close to the line color
                target_color = np.array(color_rgb) * 255
                color_diff = np.linalg.norm(pixel - target_color)

                # If pixel is close to the line color, make it lighter
                if color_diff < COLOR_TOL:  # Tolerance for color matching
                    # Blend with white to create lighter effect
                    img_array[px_y, px_x, :3] = pixel * alpha_blend + 255 * (1 - alpha_blend)
                    pixels_modified += 1

    buf.close()
    return Image.fromarray(img_array.astype(np.uint8)), pixels_modified


def export_combined_csv(sequence_name, disorder_csv, disorder_pos, disorder_scores,
                        disorder_ma, cc_entries, cc_ma_values, output_path):
    """
    Export combined Disorder/Coiled-coil scores (original + MA) to CSV using
    central moving averages and zero-filling for missing CC probabilities.
    """
    # Collect residues from disorder CSV if available
    residues = []
    try:
        with open(disorder_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                residues.append(row.get('Residue', ''))
    except FileNotFoundError:
        residues = []

    if len(residues) != len(disorder_pos):
        residues = [''] * len(disorder_pos)

    # Build lookup map for MARCOIL data (only available positions)
    cc_map = {}
    for entry, ma_val in zip(cc_entries, cc_ma_values):
        cc_map[entry['position']] = {
            'prob_orig': entry['probability'],
            'prob_ma': ma_val,
            'heptad': entry['heptad_phase']
        }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Position',
            'Residue',
            'Disorder_Score_Original',
            'Disorder_Score_3res_MA',
            'CC_Probability_Original',
            'CC_Probability_3res_MA',
            'Heptad_Phase',
            'Data_Source'
        ])

        for idx, pos in enumerate(disorder_pos):
            residue = residues[idx] if idx < len(residues) else ''
            disorder_orig = disorder_scores[idx]
            disorder_ma_val = disorder_ma[idx]

            if pos in cc_map:
                cc_orig = cc_map[pos]['prob_orig']
                cc_ma_val = cc_map[pos]['prob_ma']
                heptad = cc_map[pos]['heptad']
                source = 'MARCOIL'
            else:
                cc_orig = 0.0
                cc_ma_val = 0.0
                heptad = 'N/A'
                source = 'Filled (0.0%)'

            writer.writerow([
                pos,
                residue,
                f"{disorder_orig:.4f}",
                f"{disorder_ma_val:.4f}",
                f"{cc_orig:.4f}",
                f"{cc_ma_val:.4f}",
                heptad,
                source
            ])

    print(f"  ➜ Exported combined CSV to {output_path}")


def plot_combined_analysis(problist_file, disorder_csv, domains_file=None,
                          output_file=None, sequence_name=None,
                          create_threshold_version=True, alpha_blend=0.2,
                          export_csv=False, csv_output_dir=None):
    """
    Create a two-panel plot with disorder scores and coiled-coil probability

    Args:
        problist_file: Path to MARCOIL ProbList file
        disorder_csv: Path to disorder scores CSV
        domains_file: Path to MARCOIL Domains file (optional)
        output_file: Output filename (optional)
        sequence_name: Name of sequence to plot (optional)
        create_threshold_version: If True, also creates threshold-based visualization (default: True)
        alpha_blend: Blend ratio for threshold version (default: 0.2)
        export_csv: If True, export combined Disorder/CC CSV (default: False)
        csv_output_dir: Directory to place CSV files (default: current directory)
    """

    # Parse probability data
    prob_data = parse_problist(problist_file)

    if not prob_data:
        print("No probability data found in file")
        return

    # Parse disorder scores
    disorder_pos, disorder_scores = parse_disorder_scores(disorder_csv)

    # Parse domain data if available
    domains = []
    if domains_file and Path(domains_file).exists():
        domains = parse_domains_for_plot(domains_file)

    # Determine which sequence to plot
    if sequence_name:
        sequences_to_plot = [sequence_name] if sequence_name in prob_data else []
    else:
        sequences_to_plot = list(prob_data.keys())

    if not sequences_to_plot:
        print(f"No sequences found to plot")
        return

    # Create plot for each sequence
    for seq_name in sequences_to_plot:
        data = prob_data[seq_name]

        if not data:
            continue

        # Extract positions and probabilities
        positions = [d['position'] for d in data]
        probabilities = [d['probability'] for d in data]

        # Calculate 3-residue moving average
        prob_ma = moving_average(probabilities, window=3)
        disorder_ma = moving_average(disorder_scores, window=3)

        # Create figure with 2 subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))

        # Get domain regions for this sequence
        seq_domains = [d for d in domains if d['sequence'] == seq_name]

        # ===== TOP PANEL: Original data =====
        # Create twin axis for disorder
        ax1_twin = ax1.twinx()

        # Plot disorder score on left axis
        line1 = ax1.plot(disorder_pos, disorder_scores, linewidth=2.0, color='red',
                        alpha=0.7, label='Disorder Score')
        ax1.set_ylabel('Disorder Score (%)', fontsize=12, fontweight='bold', color='black')
        ax1.tick_params(axis='y', labelcolor='black')
        ax1.set_ylim(-5, 105)

        # Plot coiled-coil on right axis
        line2 = ax1_twin.plot(positions, probabilities, linewidth=2.0, color='blue',
                             alpha=0.7, label='Coiled-coil Probability')
        ax1_twin.axhline(y=50, color='orange', linestyle='--', linewidth=1, alpha=0.7)
        ax1_twin.axhline(y=90, color='gray', linestyle='--', linewidth=1, alpha=0.7)
        ax1_twin.set_ylabel('Coiled-coil Probability (%)', fontsize=12, fontweight='bold', color='black')
        ax1_twin.tick_params(axis='y', labelcolor='black')
        ax1_twin.set_ylim(-5, 105)

        ax1.set_title(f'{seq_name} - Original Data (Disorder vs Coiled-coil)', fontsize=14, fontweight='bold')
        ax1.set_xlim(min(positions), max(positions))
        ax1.grid(True, alpha=0.3, linestyle=':')
        ax1.set_xticks(list(range(0, max(positions)+1, 200)))
        ax1.tick_params(axis='x', labelbottom=True)

        # Combined legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='red', linewidth=2, alpha=0.7, label='Disorder Score'),
            Line2D([0], [0], color='blue', linewidth=2, alpha=0.7, label='Coiled-coil Probability')
        ]
        ax1.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.08),
                  ncol=2, framealpha=0.95, fontsize=11, fancybox=True)
        ax1.set_xlabel('Position in Sequence', fontsize=12, fontweight='bold', labelpad=25)

        # ===== BOTTOM PANEL: 3-residue moving average =====
        # Create twin axis for disorder
        ax2_twin = ax2.twinx()

        # Plot disorder score MA on left axis
        line3 = ax2.plot(disorder_pos, disorder_ma, linewidth=2.0, color='red',
                        alpha=0.7, label='Disorder Score (3-res MA)')
        ax2.set_ylabel('Disorder Score (%)', fontsize=12, fontweight='bold', color='black')
        ax2.tick_params(axis='y', labelcolor='black')
        ax2.set_ylim(-5, 105)

        # Plot coiled-coil MA on right axis
        line4 = ax2_twin.plot(positions, prob_ma, linewidth=2.0, color='blue',
                             alpha=0.7, label='Coiled-coil Probability (3-res MA)')
        ax2_twin.axhline(y=50, color='orange', linestyle='--', linewidth=1, alpha=0.7)
        ax2_twin.axhline(y=90, color='gray', linestyle='--', linewidth=1, alpha=0.7)
        ax2_twin.set_ylabel('Coiled-coil Probability (%)', fontsize=12, fontweight='bold', color='black')
        ax2_twin.tick_params(axis='y', labelcolor='black')
        ax2_twin.set_ylim(-5, 105)

        ax2.set_title(f'{seq_name} - 3-residue Moving Average (Disorder vs Coiled-coil)',
                     fontsize=14, fontweight='bold')
        ax2.set_xlim(min(positions), max(positions))
        ax2.grid(True, alpha=0.3, linestyle=':')
        ax2.set_xticks(list(range(0, max(positions)+1, 200)))
        ax2.tick_params(axis='x', labelbottom=True)

        # Combined legend
        legend_elements = [
            Line2D([0], [0], color='red', linewidth=2, alpha=0.7, label='Disorder Score (3-res MA)'),
            Line2D([0], [0], color='blue', linewidth=2, alpha=0.7, label='Coiled-coil Probability (3-res MA)')
        ]
        ax2.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.08),
                  ncol=2, framealpha=0.95, fontsize=11, fancybox=True)
        ax2.set_xlabel('Position in Sequence', fontsize=12, fontweight='bold', labelpad=30)

        plt.tight_layout(rect=(0, 0.02, 1, 0.98))

        # Save normal version
        if output_file:
            normal_output = output_file
        else:
            normal_output = f"{seq_name.replace(' ', '_')}_combined_analysis.png"

        plt.savefig(normal_output, dpi=300, bbox_inches='tight')
        print(f"Normal plot saved to {normal_output}")

        # Create threshold version if requested
        if create_threshold_version:
            print("\nCreating threshold-based visualization...")

            # Determine threshold output filename
            if output_file:
                base = output_file.rsplit('.', 1)[0]
                ext = output_file.rsplit('.', 1)[1] if '.' in output_file else 'png'
                threshold_output = f"{base}_threshold.{ext}"
            else:
                threshold_output = f"{seq_name.replace(' ', '_')}_combined_analysis_threshold.png"

            # Process each panel and axis
            print("  Processing top panel - Disorder (red)...")
            img_modified, count = apply_threshold_lightening(fig, ax1, (1.0, 0.0, 0.0),
                                                             threshold=50.0, alpha_blend=alpha_blend)
            print(f"    Modified {count} pixels")

            # Convert to array for further processing
            img_array = np.array(img_modified).astype(float)

            # Process top panel - Coiled-coil (blue)
            print("  Processing top panel - Coiled-coil (blue)...")
            bbox = ax1_twin.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            dpi = 300
            left_px = int(bbox.x0 * dpi)
            right_px = int(bbox.x1 * dpi)
            top_px = int((fig.get_figheight() - bbox.y1) * dpi)
            bottom_px = int((fig.get_figheight() - bbox.y0) * dpi)
            ylim = ax1_twin.get_ylim()

            threshold_y_normalized = (50.0 - ylim[0]) / (ylim[1] - ylim[0])
            threshold_y_pixel = top_px + (bottom_px - top_px) * (1 - threshold_y_normalized)

            height, width = img_array.shape[:2]
            pixels_modified = 0
            start_x = max(0, left_px - 5)

            for px_x in range(start_x, min(right_px + 1, width)):
                for px_y in range(max(0, top_px), min(bottom_px + 1, height)):
                    if px_y >= threshold_y_pixel:
                        pixel = img_array[px_y, px_x, :3]
                        target_color = np.array([0.0, 0.0, 1.0]) * 255
                        color_diff = np.linalg.norm(pixel - target_color)

                        if color_diff < 150:
                            img_array[px_y, px_x, :3] = pixel * alpha_blend + 255 * (1 - alpha_blend)
                            pixels_modified += 1

            print(f"    Modified {pixels_modified} pixels")
            img_modified = Image.fromarray(img_array.astype(np.uint8))

            # Process bottom panel - Disorder (red)
            print("  Processing bottom panel - Disorder (red)...")
            img_array = np.array(img_modified).astype(float)

            bbox = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            left_px = int(bbox.x0 * dpi)
            right_px = int(bbox.x1 * dpi)
            top_px = int((fig.get_figheight() - bbox.y1) * dpi)
            bottom_px = int((fig.get_figheight() - bbox.y0) * dpi)
            ylim = ax2.get_ylim()

            threshold_y_normalized = (50.0 - ylim[0]) / (ylim[1] - ylim[0])
            threshold_y_pixel = top_px + (bottom_px - top_px) * (1 - threshold_y_normalized)

            pixels_modified = 0
            start_x = max(0, left_px - 5)

            for px_x in range(start_x, min(right_px + 1, width)):
                for px_y in range(max(0, top_px), min(bottom_px + 1, height)):
                    if px_y >= threshold_y_pixel:
                        pixel = img_array[px_y, px_x, :3]
                        target_color = np.array([1.0, 0.0, 0.0]) * 255
                        color_diff = np.linalg.norm(pixel - target_color)

                        if color_diff < 150:
                            img_array[px_y, px_x, :3] = pixel * alpha_blend + 255 * (1 - alpha_blend)
                            pixels_modified += 1

            print(f"    Modified {pixels_modified} pixels")
            img_modified = Image.fromarray(img_array.astype(np.uint8))

            # Process bottom panel - Coiled-coil (blue)
            print("  Processing bottom panel - Coiled-coil (blue)...")
            img_array = np.array(img_modified).astype(float)

            bbox = ax2_twin.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            left_px = int(bbox.x0 * dpi)
            right_px = int(bbox.x1 * dpi)
            top_px = int((fig.get_figheight() - bbox.y1) * dpi)
            bottom_px = int((fig.get_figheight() - bbox.y0) * dpi)
            ylim = ax2_twin.get_ylim()

            threshold_y_normalized = (50.0 - ylim[0]) / (ylim[1] - ylim[0])
            threshold_y_pixel = top_px + (bottom_px - top_px) * (1 - threshold_y_normalized)

            pixels_modified = 0
            start_x = max(0, left_px - 5)

            for px_x in range(start_x, min(right_px + 1, width)):
                for px_y in range(max(0, top_px), min(bottom_px + 1, height)):
                    if px_y >= threshold_y_pixel:
                        pixel = img_array[px_y, px_x, :3]
                        target_color = np.array([0.0, 0.0, 1.0]) * 255
                        color_diff = np.linalg.norm(pixel - target_color)

                        if color_diff < 150:
                            img_array[px_y, px_x, :3] = pixel * alpha_blend + 255 * (1 - alpha_blend)
                            pixels_modified += 1

            print(f"    Modified {pixels_modified} pixels")
            img_modified = Image.fromarray(img_array.astype(np.uint8))

            # Save threshold version
            img_modified.save(threshold_output, dpi=(300, 300))
            print(f"\nThreshold plot saved to {threshold_output}")

        if export_csv:
            csv_dir = Path(csv_output_dir) if csv_output_dir else Path('.')
            csv_filename = f"{seq_name.replace(' ', '_')}_combined_scores.csv"
            csv_path = csv_dir / csv_filename
            export_combined_csv(
                seq_name,
                disorder_csv,
                disorder_pos,
                disorder_scores,
                disorder_ma,
                data,
                prob_ma,
                csv_path
            )

        plt.close()


if __name__ == "__main__":
    from pathlib import Path

    problist_file = Path('MARCOIL/Outputs/ProbList')
    domains_file = Path('MARCOIL/Outputs/Domains')
    disorder_csv = Path('ELKS2_disorder_scores.csv')

    plot_combined_analysis(
        problist_file,
        disorder_csv,
        domains_file,
        'ELKS2_SS306_combined_analysis_20251016.png',
        'ELKS2_SS306',
        create_threshold_version=True,  # Set to False to skip threshold version
        alpha_blend=0.2,  # Adjust lightening (0.1=very light, 0.3=darker)
        export_csv=True,
        csv_output_dir=Path('.')
    )
