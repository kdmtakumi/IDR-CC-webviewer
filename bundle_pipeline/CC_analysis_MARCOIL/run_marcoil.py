#!/usr/bin/env python3
"""
MARCOIL Coiled-Coil Prediction Wrapper
Runs coiled-coil prediction on FASTA sequences and outputs results to CSV
"""

import os
import sys
import subprocess
import csv
from pathlib import Path
import re


def check_marcoil_executable(marcoil_dir):
    """Check if marcoil executable exists"""
    marcoil_path = Path(marcoil_dir) / "marcoil"
    return marcoil_path.exists() and os.access(marcoil_path, os.X_OK)


def compile_marcoil(marcoil_dir):
    """Attempt to compile MARCOIL (requires fixing old C++ headers)"""
    print("Attempting to compile MARCOIL...")
    try:
        result = subprocess.run(
            ["make"],
            cwd=marcoil_dir,
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.returncode == 0:
            print("Compilation successful!")
            return True
        else:
            print(f"Compilation failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"Compilation error: {e}")
        return False


def run_marcoil(marcoil_dir, fasta_file, mode="H", options="+dlsS"):
    """
    Run MARCOIL on a FASTA file

    Args:
        marcoil_dir: Path to MARCOIL directory
        fasta_file: Path to input FASTA file
        mode: "H" (high sensitivity) or "L" (low sensitivity)
        options: MARCOIL output options

    Returns:
        Path to output directory
    """
    marcoil_path = Path(marcoil_dir) / "marcoil"

    if not marcoil_path.exists():
        raise FileNotFoundError(f"MARCOIL executable not found at {marcoil_path}")

    # Set mode flag
    mode_flag = f"-{mode}"

    # Convert fasta_file to absolute path
    fasta_file = Path(fasta_file).resolve()

    # Run MARCOIL
    cmd = [str(marcoil_path), mode_flag, options, str(fasta_file)]

    print(f"Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            cwd=marcoil_dir,
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            print(f"MARCOIL error: {result.stderr}")
            return None

        print(result.stdout)
        return Path(marcoil_dir) / "Outputs"

    except subprocess.TimeoutExpired:
        print("MARCOIL execution timed out")
        return None
    except Exception as e:
        print(f"Error running MARCOIL: {e}")
        return None


def parse_domains_file(domains_file):
    """
    Parse MARCOIL Domains output file

    Returns list of dicts with domain information
    """
    results = []
    current_seq = None
    current_threshold = None

    with open(domains_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Sequence name line starts with '>'
            if line.startswith('>'):
                # Extract sequence name (before ##)
                seq_info = line[1:].strip()
                if '##' in seq_info:
                    current_seq = seq_info.split('##')[0].strip()
                else:
                    current_seq = seq_info
                continue

            # Threshold line: "NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 50.0 : 1"
            threshold_match = re.match(r'NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD ([\d.]+)\s*:\s*(\d+)', line)
            if threshold_match:
                current_threshold = float(threshold_match.group(1))
                continue

            # Domain line format: "  1. from 2 to 281 (length = 281) with max = 100.0"
            domain_match = re.match(r'\s*(\d+)\.\s+from\s+(\d+)\s+to\s+(\d+)\s+\(length\s*=\s*(\d+)\)\s+with\s+max\s*=\s*([\d.]+)', line)
            if domain_match and current_seq and current_threshold:
                domain_num, start, end, length, max_prob = domain_match.groups()
                results.append({
                    'sequence_name': current_seq,
                    'threshold': current_threshold,
                    'domain_number': int(domain_num),
                    'start': int(start),
                    'end': int(end),
                    'length': int(length),
                    'max_probability': float(max_prob)
                })

    return results


def parse_compact_profile(compact_file):
    """
    Parse MARCOIL CompactProfile output

    Returns dict with sequence-level statistics
    """
    results = {}
    current_seq = None

    with open(compact_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                current_seq = line[1:].strip()
                results[current_seq] = {'positions': []}
                continue

            # Parse position data
            if current_seq and line[0].isdigit():
                parts = line.split()
                if len(parts) >= 2:
                    pos = int(parts[0])
                    prob = float(parts[1])
                    results[current_seq]['positions'].append({
                        'position': pos,
                        'probability': prob
                    })

    return results


def write_csv_output(domains_data, output_csv):
    """Write domain predictions to CSV file"""

    if not domains_data:
        print("No domains found to write")
        return

    fieldnames = ['sequence_name', 'threshold', 'domain_number', 'start', 'end', 'length', 'max_probability']

    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(domains_data)

    print(f"Results written to {output_csv}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python run_marcoil.py <fasta_file> [output.csv] [mode:H/L]")
        print("Example: python run_marcoil.py sequence.fasta results.csv H")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else "marcoil_results.csv"
    mode = sys.argv[3] if len(sys.argv) > 3 else "H"

    # Get MARCOIL directory (in same directory as this script)
    script_dir = Path(__file__).parent
    marcoil_dir = script_dir / "MARCOIL"

    if not marcoil_dir.exists():
        print(f"MARCOIL directory not found at {marcoil_dir}")
        sys.exit(1)

    # Check if MARCOIL executable exists
    if not check_marcoil_executable(marcoil_dir):
        print("MARCOIL executable not found. Attempting to compile...")
        if not compile_marcoil(marcoil_dir):
            print("\nERROR: MARCOIL executable is not available.")
            print("Please compile MARCOIL manually first.")
            print("This requires fixing old C++ headers (iostream.h -> iostream, etc.)")
            sys.exit(1)

    # Run MARCOIL
    output_dir = run_marcoil(marcoil_dir, fasta_file, mode=mode)

    if not output_dir:
        print("MARCOIL execution failed")
        sys.exit(1)

    # Parse results
    domains_file = output_dir / "Domains"

    if not domains_file.exists():
        print(f"Output file not found: {domains_file}")
        sys.exit(1)

    print(f"Parsing results from {domains_file}")
    domains_data = parse_domains_file(domains_file)

    # Write CSV
    write_csv_output(domains_data, output_csv)

    print(f"\nSummary: Found {len(domains_data)} coiled-coil domain(s)")


if __name__ == "__main__":
    main()
