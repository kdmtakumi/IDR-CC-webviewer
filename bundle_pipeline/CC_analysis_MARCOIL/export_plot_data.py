#!/usr/bin/env python3
"""
Export plot data to CSV for re-plotting
Exports both original and 3-residue moving average data for any sequence
"""

import re
import csv
import sys
from pathlib import Path


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

            match = re.match(r'\s*(\d+)\s+([A-Z*])\s+([\d.]+)\s+([a-g])', line)
            if match and current_seq:
                pos, aa, prob, phase = match.groups()
                results[current_seq].append({
                    'position': int(pos),
                    'amino_acid': aa,
                    'probability': float(prob),
                    'heptad_phase': phase
                })

    return results


def parse_disorder_scores(disorder_csv):
    """Parse disorder scores CSV file"""
    positions = []
    residues = []
    scores = []

    with open(disorder_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            positions.append(int(row['Position']))
            residues.append(row['Residue'])
            # Convert to percentage if needed (0-1 range)
            score = float(row['Disorder_Score'])
            if score <= 1.0:
                score = score * 100
            scores.append(score)

    return positions, residues, scores


def moving_average(data, window=3):
    """
    Calculate a centered moving average with edge padding so exported data stay
    aligned with the comprehensive CSV (3-residue MA with edge extension).
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


def export_plot_data(problist_file, disorder_csv, output_csv, sequence_name):
    """
    Export plot data to CSV for re-plotting

    Args:
        problist_file: Path to MARCOIL ProbList file
        disorder_csv: Path to disorder scores CSV
        output_csv: Path to output CSV file
        sequence_name: Name of sequence to export
    """

    print(f"Exporting plot data for {sequence_name}...")

    # Parse probability data
    prob_data = parse_problist(problist_file)

    if sequence_name not in prob_data:
        print(f"ERROR: Sequence '{sequence_name}' not found in ProbList")
        print(f"Available sequences: {', '.join(prob_data.keys())}")
        return False

    data = prob_data[sequence_name]

    if not data:
        print("ERROR: No coiled-coil data found")
        return False

    # Extract positions and probabilities
    cc_positions = [d['position'] for d in data]
    cc_amino_acids = [d['amino_acid'] for d in data]
    cc_probabilities = [d['probability'] for d in data]

    # Parse disorder scores
    disorder_pos, disorder_residues, disorder_scores = parse_disorder_scores(disorder_csv)

    # Create mapping of CC data by position for quick lookup
    cc_data_map = {}
    for i, pos in enumerate(cc_positions):
        cc_data_map[pos] = {
            'amino_acid': cc_amino_acids[i],
            'prob_orig': cc_probabilities[i],
            'heptad': data[i]['heptad_phase']
        }

    # Build full-length coiled-coil probability series aligned to disorder positions
    cc_probs_full = [cc_data_map.get(pos, {}).get('prob_orig', 0.0) for pos in disorder_pos]

    # Calculate 3-residue moving average on full-length arrays (central + edge padding)
    disorder_ma = moving_average(disorder_scores, window=3)
    cc_prob_ma_full = moving_average(cc_probs_full, window=3)

    # Create comprehensive CSV aligned to full sequence length
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Position',
            'Residue',
            'Disorder_Score_Original',
            'Disorder_Score_3res_MA',
            'CC_Probability_Original',
            'CC_Probability_3res_MA',
            'Heptad_Phase'
        ])

        for i, pos in enumerate(disorder_pos):
            residue = disorder_residues[i]
            disorder_orig = disorder_scores[i]
            disorder_ma_val = disorder_ma[i]
            cc_orig = cc_probs_full[i]
            cc_ma_val = cc_prob_ma_full[i]

            if pos in cc_data_map:
                heptad = cc_data_map[pos]['heptad']
            else:
                heptad = '-'

            writer.writerow([
                pos,
                residue,
                f"{disorder_orig:.4f}",
                f"{disorder_ma_val:.4f}",
                f"{cc_orig:.4f}",
                f"{cc_ma_val:.4f}",
                heptad
            ])

    print(f"✓ Exported plot data to {output_csv}")
    print(f"  Total positions: {len(disorder_pos)}")
    print(f"  Coiled-coil positions: {len(cc_positions)}")
    print(f"\nColumns in CSV:")
    print(f"  - Position: Residue position (1-based)")
    print(f"  - Residue: Amino acid 1-letter code")
    print(f"  - Disorder_Score_Original: Original disorder score (%)")
    print(f"  - Disorder_Score_3res_MA: 3-residue moving average disorder score (%)")
    print(f"  - CC_Probability_Original: Original coiled-coil probability (%)")
    print(f"  - CC_Probability_3res_MA: 3-residue moving average CC probability (%)")
    print(f"  - Heptad_Phase: Heptad register (a-g, or '-' if not coiled-coil)")

    return True


def main():
    if len(sys.argv) < 4:
        print("Usage: python3 export_plot_data.py <sequence_name> <disorder_csv> <output_csv>")
        print("\nExample:")
        print("  python3 export_plot_data.py ELKS2_SS306 ELKS2_disorder_scores.csv ELKS2_plot_data.csv")
        print("\nNote: This script uses the most recent MARCOIL output in MARCOIL/Outputs/")
        print("      Make sure you have run MARCOIL analysis first.")
        sys.exit(1)

    sequence_name = sys.argv[1]
    disorder_csv = sys.argv[2]
    output_csv = sys.argv[3]

    # Use default MARCOIL output location
    script_dir = Path(__file__).parent
    problist_file = script_dir / "MARCOIL" / "Outputs" / "ProbList"

    if not problist_file.exists():
        print(f"ERROR: ProbList file not found at {problist_file}")
        print("Please run MARCOIL analysis first using run_marcoil.py")
        sys.exit(1)

    if not Path(disorder_csv).exists():
        print(f"ERROR: Disorder CSV file not found: {disorder_csv}")
        sys.exit(1)

    success = export_plot_data(problist_file, disorder_csv, output_csv, sequence_name)

    if success:
        print(f"\n✓ Success! You can now use this CSV to re-plot the data.")
        print(f"  See replot_from_csv.py for an example re-plotting script.")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
