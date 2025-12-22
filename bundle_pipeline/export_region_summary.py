#!/usr/bin/env python3
"""
Create integrated per-residue table and region intervals for IDR/CC/TM.

Inputs:
  --combined-csv: combined_scores.csv (from plot_combined_analysis.py)
  --gff: TM prediction GFF (e.g., DeepTMHMM TMRs.gff3). Optional.
  --out-prefix: output prefix (e.g., TANGO1_LAB002/TANGO1_LAB002_regions)

Outputs:
  <out-prefix>_positions.csv  : per-residue scores/flags (IDR>=50, CC>=50, TM)
  <out-prefix>_intervals.csv  : contiguous intervals for IDR/CC/TM
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, List, Tuple


def load_combined(combined_csv: Path):
    pos: List[int] = []
    res: List[str] = []
    d0: List[float] = []
    d3: List[float] = []
    c0: List[float] = []
    c3: List[float] = []

    with combined_csv.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            pos.append(int(row["Position"]))
            res.append(row.get("Residue", ""))
            d0.append(float(row["Disorder_Score_Original"]) * (100 if float(row["Disorder_Score_Original"]) <= 1 else 1))
            d3.append(float(row["Disorder_Score_3res_MA"]) * (100 if float(row["Disorder_Score_3res_MA"]) <= 1 else 1))
            c0.append(float(row["CC_Probability_Original"]) * (100 if float(row["CC_Probability_Original"]) <= 1 else 1))
            c3.append(float(row["CC_Probability_3res_MA"]) * (100 if float(row["CC_Probability_3res_MA"]) <= 1 else 1))
    return pos, res, d0, d3, c0, c3


def load_tm_gff(gff_path: Path) -> List[Tuple[int, int]]:
    spans: List[Tuple[int, int]] = []
    tm_labels = {"TMhelix", "TMH", "TM", "TRANSMEM"}
    if not gff_path:
        return spans
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


def main():
    parser = argparse.ArgumentParser(description="Export integrated IDR/CC/TM positions and intervals.")
    parser.add_argument("--combined-csv", required=True, help="combined_scores.csv")
    parser.add_argument("--gff", help="TM prediction GFF (DeepTMHMM/TMHMM)")
    parser.add_argument("--out-prefix", required=True, help="Output prefix for CSVs")
    parser.add_argument("--threshold", type=float, default=50.0, help="Threshold for IDR/CC (default 50%)")
    args = parser.parse_args()

    combined_csv = Path(args.combined_csv)
    gff_path = Path(args.gff) if args.gff else None
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    pos, res, d0, d3, c0, c3 = load_combined(combined_csv)
    thr = args.threshold

    idr_flags = [x >= thr for x in d0]
    cc_flags = [x >= thr for x in c0]
    tm_spans = load_tm_gff(gff_path) if gff_path else []
    tm_flags = [any(s <= p <= e for s, e in tm_spans) for p in pos] if tm_spans else [False] * len(pos)

    # Per-residue output
    pos_csv = out_prefix.with_name(out_prefix.name + "_positions.csv")
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
            "In_TM",
        ])
        for p, r, d_orig, d_ma, c_orig, c_ma, fi, fc, ft in zip(pos, res, d0, d3, c0, c3, idr_flags, cc_flags, tm_flags):
            writer.writerow([
                p,
                r,
                f"{d_orig:.2f}",
                f"{d_ma:.2f}",
                "Yes" if fi else "No",
                f"{c_orig:.2f}",
                f"{c_ma:.2f}",
                "Yes" if fc else "No",
                "Yes" if ft else "No",
            ])

    # Intervals summary
    interval_csv = out_prefix.with_name(out_prefix.name + "_intervals.csv")
    with interval_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Type", "Start", "End", "Length"])
        for t, flags in (("IDR", idr_flags), ("CC", cc_flags)):
            for s, e in flags_to_intervals(flags):
                writer.writerow([t, s, e, e - s + 1])
        for s, e in tm_spans:
            writer.writerow(["TM", s, e, e - s + 1])

    print(f"Saved per-residue: {pos_csv}")
    print(f"Saved intervals  : {interval_csv}")


if __name__ == "__main__":
    main()
