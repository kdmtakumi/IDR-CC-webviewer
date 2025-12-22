Pipeline bundle for IDR / coiled-coil / membrane plotting
=========================================================

This folder packages the scripts and MARCOIL assets needed to run the single-sequence pipeline in any location. Only external dependencies are Python libraries (and DeepTMHMM if you want membrane spans). DeepTMHMM is optional; when omitted, you still get IDR/CC plots and region CSVs.

Contents
- `analyze_single_sequence.py` — end-to-end IDR (metapredict) + MARCOIL + combined plots (light/dark, threshold).
- `plot_deeptmhmm_overlay_threshold.py` — overlays TM spans from DeepTMHMM/TMHMM GFF on the combined plots and exports region tables.
- `plot_overlay_threshold_basic.py` — same plotting without TM input (use when you skip DeepTMHMM).
- `export_region_summary.py` — standalone region table exporter (IDR/CC + optional TM).
- `CC_analysis_MARCOIL/` — MARCOIL binary and plotting helpers (`run_marcoil.py`, `export_plot_data.py`, `plot_combined_analysis.py`, `plot_combined_dark_theme.py`, `replot_from_csv.py`, docs).

Requirements
- Python 3.9+ with `pip install metapredict matplotlib numpy pillow biolib` (biolib only if you run DeepTMHMM).
- MARCOIL binary is included under `CC_analysis_MARCOIL/MARCOIL/marcoil`; ensure it is executable: `chmod +x CC_analysis_MARCOIL/MARCOIL/marcoil`.
- Run commands from inside this folder and keep `PYTHONPATH=.` so the local `CC_analysis_MARCOIL` package resolves.
- DeepTMHMM via `biolib.load('DTU/DeepTMHMM')` downloads the model the first time; internet access is required for that initial fetch.

Quick start (from this folder)
1) Analyze one sequence (IDR + MARCOIL)
```
PYTHONPATH=. MPLCONFIGDIR=.matplotlib_cache XDG_CACHE_HOME=.cache \
  python analyze_single_sequence.py \
    --name SAMPLE \
    --fasta /path/to/input.fasta \
    --marcoil-dir CC_analysis_MARCOIL/MARCOIL \
    --mode H
```
Outputs go to `./SAMPLE/` (FASTA, disorder/CC plot data CSV, combined plots, thresholds, region CSVs).

2) (Optional) Run DeepTMHMM to get TM spans (needs internet on first run)
```
python - <<'PY'
import biolib
deeptmhmm = biolib.load('DTU/DeepTMHMM')
job = deeptmhmm.cli(args='--fasta SAMPLE/SAMPLE.fasta')
job.save_files('deeptmhmm_SAMPLE')
print('DeepTMHMM outputs in deeptmhmm_SAMPLE/')
PY
```
Use `deeptmhmm_SAMPLE/prediction.gff3` or `deeptmhmm_SAMPLE/TMRs.gff3` for plotting.

3a) Overlay DeepTMHMM on the combined plots (recommended over legacy TMHMM)
```
PYTHONPATH=. MPLCONFIGDIR=.matplotlib_cache XDG_CACHE_HOME=.cache \
  python plot_deeptmhmm_overlay_threshold.py \
    --name SAMPLE \
    --plot-csv SAMPLE/SAMPLE_plot_data.csv \
    --gff deeptmhmm_SAMPLE/prediction.gff3 \
    --out-prefix SAMPLE/SAMPLE_deeptmhmm_overlay_threshold \
    --xtick-step 200 \
    --threshold 50
```
Generates light/dark PNGs and region summary CSVs (`*_regions_positions.csv`, `*_regions_intervals.csv`).

3b) If you skip DeepTMHMM (no TM data), make the threshold plots anyway
```
PYTHONPATH=. MPLCONFIGDIR=.matplotlib_cache XDG_CACHE_HOME=.cache \
  python plot_overlay_threshold_basic.py \
    --name SAMPLE \
    --plot-csv SAMPLE/SAMPLE_plot_data.csv \
    --out-prefix SAMPLE/SAMPLE_overlay_threshold_noTM \
    --xtick-step 200 \
    --threshold 50
```

4) Export region tables separately (IDR/CC threshold=50%, TM from GFF if provided)
```
PYTHONPATH=. \
  python export_region_summary.py \
    --combined-csv SAMPLE/SAMPLE_combined_scores.csv \
    --out-prefix SAMPLE/SAMPLE_region_summary \
    --threshold 50 \
    --gff deeptmhmm_SAMPLE/prediction.gff3   # omit this flag if no TM data
```
Outputs:
- `<out_prefix>_regions_positions.csv` (per-residue scores + calls)
- `<out_prefix>_regions_intervals.csv` (IDR/CC/TM contiguous ranges)

Notes
- Keep `--xtick-step` consistent if you want axis design identical across plots (use 200 to match defaults).
- You can switch MARCOIL sensitivity with `--mode L` (default is `H`).
- DeepTMHMM is preferred over legacy TMHMM for membrane prediction; TM input is optional, and plots still work without it.
