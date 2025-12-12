# CC_analysis_MARCOIL - Coiled-Coil Analysis Toolkit

このフォルダには、タンパク質配列のコイルドコイル解析とDisorder Score解析を統合的に行うために必要なすべてのファイルが含まれています。

## 📚 ドキュメント一覧

- **README.md** (このファイル) - 全体概要と各スクリプトの詳細
- **HOW_TO_USE_YOUR_DATA.md** - 🌟自分のデータで解析する手順（ステップバイステップ）
- **QUICK_START.md** - 最短実行手順
- **FILE_FORMATS.md** - ファイル形式の詳細仕様
- **TEST_RUN.sh** - サンプルを使った動作確認スクリプト

## 🚀 初めての方へ

**自分のデータで解析したい場合**:
→ **[HOW_TO_USE_YOUR_DATA.md](HOW_TO_USE_YOUR_DATA.md)** を読んでください

**まず動作確認したい場合**:
```bash
./TEST_RUN.sh
```

## フォルダ構成

```
CC_analysis_MARCOIL/
├── README.md                          # このファイル
├── MARCOIL/                           # MARCOILプログラムディレクトリ
│   ├── marcoil                        # MARCOIL実行ファイル
│   ├── Inputs/                        # MARCOILパラメータファイル
│   └── Outputs/                       # MARCOIL出力ディレクトリ
├── run_marcoil.py                     # MARCOIL実行スクリプト
├── plot_coiled_coil.py                # 基本プロットスクリプト
├── plot_with_moving_average.py        # 移動平均プロットスクリプト
├── plot_combined_analysis.py          # Disorder+CC統合プロットスクリプト
├── export_plot_data.py                # プロットデータCSVエクスポートスクリプト
├── replot_from_csv.py                 # CSVからの再プロットスクリプト
├── example_input.fasta                # サンプル入力FASTAファイル
└── example_disorder_scores.csv        # サンプルDisorderスコアCSV
```

## 必要な環境

- Python 3.x
- 必要なPythonパッケージ:
  ```bash
  pip install matplotlib numpy
  ```

## 使い方

### 1. MARCOIL解析の実行

配列からコイルドコイル領域を予測します。

```bash
python3 run_marcoil.py <入力FASTA> [出力CSV] [モード:H/L]
```

**パラメータ:**
- `<入力FASTA>`: 解析するFASTAファイル（必須）
- `[出力CSV]`: 結果CSVファイル名（省略可、デフォルト: `marcoil_results.csv`）
- `[モード]`: `H`（高感度、デフォルト）または `L`（低感度）

**実行例:**
```bash
python3 run_marcoil.py example_input.fasta results.csv H
```

**出力:**
- CSVファイル: ドメイン情報（開始位置、終了位置、長さ、最大確率）
- `MARCOIL/Outputs/`ディレクトリ内にも詳細結果が保存されます

---

### 2. 基本的なコイルドコイルプロット

MARCOILの結果を可視化します。

```bash
python3 plot_coiled_coil.py <配列名> [出力PNG]
```

**実行例:**
```bash
python3 -c "
from pathlib import Path
import sys
from plot_coiled_coil import plot_coiled_coil_probability

plot_coiled_coil_probability(
    Path('MARCOIL/Outputs/ProbList'),
    Path('MARCOIL/Outputs/Domains'),
    'output_plot.png',
    '配列名'
)
"
```

**出力:**
- コイルドコイル確率のグラフ（PNG）
- 50%と90%の閾値ライン表示

---

### 3. 移動平均プロット

3残基移動平均を適用したプロットを作成します。すべての移動平均処理は
「中央寄せ＋エッジ拡張（端は隣接ウィンドウ平均を複製）」で統一しており、
`ELKS2_SS306+M_MA_data_corrected.csv` と同じ値になります。

**使い方:**

スクリプトを編集して配列名を指定:

```python
# plot_with_moving_average.py の最後の部分を編集
plot_with_moving_average(
    problist_file,
    domains_file,
    '出力ファイル名.png',
    '配列名'
)
```

実行:
```bash
python3 plot_with_moving_average.py
```

**出力:**
- 2パネルのグラフ
  - 上: オリジナルデータ
  - 下: 3残基移動平均
- `export_csv=True` にすると `ELKS2_SS306_combined_scores.csv` のようなCSVも生成され、
  Disorder/Coiled-coilのオリジナル値と3残基移動平均を同時に確認できます。CC確率は
  0.5%未満でレポートされない残基も自動的に0.0で補完されます。

---

### 4. Disorder Score + Coiled-Coil統合解析

Disorder ScoreとCoiled-Coil確率を重ね合わせて表示します。

**前提条件:**
- MARCOILを実行済み（`MARCOIL/Outputs/`に結果がある）
- Disorder ScoreのCSVファイルがある

**Disorder Score CSVの形式:**
```csv
Position,Residue,Disorder_Score,In_IDR
1,M,0.95,Yes
2,A,0.92,Yes
...
```

**使い方:**

スクリプトを編集してファイル名と配列名を指定:

```python
# plot_combined_analysis.py の最後の部分を編集
problist_file = Path('MARCOIL/Outputs/ProbList')
domains_file = Path('MARCOIL/Outputs/Domains')
disorder_csv = Path('あなたのDisorderスコア.csv')

plot_combined_analysis(
    problist_file,
    disorder_csv,
    domains_file,
    '出力ファイル名.png',
    '配列名'
)
```

実行:
```bash
python3 plot_combined_analysis.py
```

**出力:**
- 2パネルの統合グラフ
  - 上: オリジナルデータ（Disorder vs CC）
  - 下: 3残基移動平均（Disorder vs CC）
  - 左軸（赤）: Disorder Score (%)
  - 右軸（青）: Coiled-Coil Probability (%)

---

### 5. プロットデータのエクスポート（再プロット用）

解析結果をCSV形式でエクスポートして、他の人が再プロットできるようにします。この時、MARCOILの結果は0.5%などの閾値未満の残基がスキップされる可能性があるので、実際の位置とズレる可能性があります。

**使い方:**

```bash
python3 export_plot_data.py <配列名> <disorder_csv> <出力CSV>
```

**実行例:**
```bash
python3 export_plot_data.py ELKS2_SS306 ELKS2_disorder_scores.csv ELKS2_plot_data.csv
```

**出力:**
- CSVファイル（以下のカラムを含む）:
  - `Position`: 残基位置
  - `Residue`: アミノ酸1文字コード
  - `Disorder_Score_Original`: オリジナルDisorderスコア(%)
  - `Disorder_Score_3res_MA`: 3残基移動平均Disorderスコア(%)
  - `CC_Probability_Original`: オリジナルCoiled-Coil確率(%)
  - `CC_Probability_3res_MA`: 3残基移動平均CC確率(%)
  - `Heptad_Phase`: ヘプタッドレジスター(a-g、またはn/a)

---

### 6. CSVデータから再プロット

エクスポートされたCSVファイルからグラフを再作成します。

**使い方:**

```bash
python3 replot_from_csv.py <CSVファイル> [出力プレフィックス] [タイトル]
```

**実行例:**
```bash
python3 replot_from_csv.py ELKS2_plot_data.csv ELKS2_replot "ELKS2 Analysis"
```

**出力:**
- `{prefix}_original.png` - オリジナルデータプロット
- `{prefix}_ma.png` - 3残基移動平均プロット
- `{prefix}_2panel.png` - 2パネル統合プロット（上:オリジナル、下:移動平均）

**カスタマイズ:**
`replot_from_csv.py`を編集して、色、線の太さ、フォントサイズなどを変更できます。

---

### 7. 黒背景＋閾値フェード版の生成（汎用スクリプト）

`export_plot_data.py` で得た `*_plot_data.csv` から、黒背景・白文字の統合グラフと、50%未満を黒側へブレンドした閾値版を自動生成します。凡例・タイトル・軸位置は通常版と同じ配置です。

```bash
# 例: ELKS2_SS306+M_plot_data.csv から黒背景版を生成
python3 plot_combined_dark_theme.py \
  --csv ../ELKS2_SS306+M/ELKS2_SS306+M_plot_data.csv \
  --seq-name "ELKS2_SS306+M" \
  --output-prefix ../ELKS2_SS306+M/ELKS2_SS306+M_combined_analysis \
  --xtick-step 200
```

**出力:**
- `<prefix>_dark.png` … 黒背景の通常版（2パネル）
- `<prefix>_dark_threshold.png` … 50%未満を黒側へブレンドした閾値版

**補足:**
- 目盛り間隔は `--xtick-step`（デフォルト200）。
- matplotlibのキャッシュ書き込みに失敗する環境向けに、スクリプト内でローカル`.mplconfig`を自動作成しています。

## 完全ワークフロー例

ELKS2(SS306)タンパク質の解析例:

```bash
# 1. MARCOIL解析を実行
python3 run_marcoil.py example_input.fasta ELKS2_results.csv H

# 2. 統合プロットを作成（Disorder + CC）
# plot_combined_analysis.py を編集:
# - disorder_csv = Path('example_disorder_scores.csv')
# - 配列名 = 'ELKS2_SS306'
# - 出力ファイル名 = 'ELKS2_combined.png'

python3 plot_combined_analysis.py

# 3. プロットデータをエクスポート（再プロット用）
python3 export_plot_data.py ELKS2_SS306 example_disorder_scores.csv ELKS2_plot_data.csv

# 4. CSVから再プロット（他の人がカスタマイズして使用可能）
python3 replot_from_csv.py ELKS2_plot_data.csv ELKS2_replot "ELKS2 Analysis"
```

---

## 各スクリプトの詳細

### run_marcoil.py
- MARCOILプログラムのPythonラッパー
- FASTAファイルを入力として受け取る
- コイルドコイルドメインをCSV形式で出力

### plot_coiled_coil.py
- MARCOILの`ProbList`と`Domains`ファイルを読み込み
- コイルドコイル確率を位置ごとにプロット
- 閾値ラインとドメイン領域を表示

### plot_with_moving_average.py
- `plot_coiled_coil.py`の拡張版
- 3残基移動平均を計算して表示
- オリジナルとスムージングデータを2パネルで比較
- 移動平均は中央寄せ＋エッジ拡張で `ELKS2_SS306+M_MA_data_corrected.csv` と一致
- MARCOILが閾値未満で省略した残基も0.0で補完したうえで移動平均を計算する点が重要
- ターミナルから一括解析する手順は「カスタム解析ワークフロー」を参照

### plot_combined_analysis.py
- Disorder ScoreとCoiled-Coil確率を統合
- 双軸グラフで2つのデータを重ね合わせ
- 移動平均も適用可能
- 左軸: Disorder Score (赤)
- 右軸: Coiled-Coil Probability (青)
- 3残基移動平均を使う場合も同じ中央寄せ処理を適用
- `export_csv=True` を指定すると、各残基についてDisorder/CCのオリジナル値と
  3残基移動平均をまとめたCSVを出力（CC確率が0.5%未満でスキップされた位置も0.0で補完）

### export_plot_data.py
- 解析結果（オリジナルと3残基移動平均）をCSV形式でエクスポート
- 他の人が再プロットやカスタマイズできるようにデータを共有
- DisorderスコアとCoiled-Coil確率の両方を含む
- ヘプタッドレジスター情報も含む
- エクスポートする移動平均列は中央寄せ＋エッジ拡張で計算
- MARCOIL欠番は0.0で埋め戻してから移動平均を取っているため、手動処理時もこの順序を守る
- ターミナル操作コマンドは下部の手順例を参照

### replot_from_csv.py
- エクスポートされたCSVファイルからグラフを再作成
- 3種類のプロット（オリジナル、移動平均、2パネル）を自動生成
- スクリプトを編集してカスタマイズ可能
- Pythonの知識があれば、色、サイズ、ラベルなど自由に変更可能

---

## カスタム解析ワークフロー（ターミナル実行例）

任意の配列をターミナルから解析する際の例です。`MyProtein` や `YOUR_AMINO_ACID_SEQUENCE` を目的の名称・配列に置き換えてください。

```bash
cd "/path/to/IDR+CC coanalysis"
test -d venv || python3 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install metapredict matplotlib pillow
```

```bash
python3 - <<'PY'
from pathlib import Path
seq = "YOUR_AMINO_ACID_SEQUENCE"
fasta = ">MyProtein\n" + "\n".join(seq[i:i+60] for i in range(0, len(seq), 60)) + "\n"
Path("MyProtein").mkdir(parents=True, exist_ok=True)
Path("MyProtein/MyProtein.fasta").write_text(fasta)
PY
```

```bash
python3 - <<'PY'
from pathlib import Path
import csv
import metapredict as meta
sequence = "".join(line.strip() for line in Path("MyProtein/MyProtein.fasta").read_text().splitlines() if not line.startswith(">"))
scores = meta.predict_disorder(sequence, version=3)
domains = meta.predict_disorder_domains(sequence, version=3, disorder_threshold=0.5).disordered_domain_boundaries
with open("MyProtein/MyProtein_disorder_scores.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Position", "Residue", "Disorder_Score", "In_IDR"])
    for idx, (res, score) in enumerate(zip(sequence, scores), 1):
        in_idr = any(start <= idx <= end for start, end in domains)
        writer.writerow([idx, res, f"{score:.4f}", "Yes" if in_idr else "No"])
PY
```

```bash
python3 CC_analysis_MARCOIL/run_marcoil.py MyProtein/MyProtein.fasta MyProtein/MyProtein_coiled_coil_results.csv H
python3 CC_analysis_MARCOIL/export_plot_data.py MyProtein MyProtein/MyProtein_disorder_scores.csv MyProtein/MyProtein_plot_data.csv
```

```bash
mkdir -p .matplotlib_cache
MPLCONFIGDIR=$PWD/.matplotlib_cache MPLBACKEND=Agg PYTHONPATH=CC_analysis_MARCOIL python3 - <<'PY'
from pathlib import Path
from plot_combined_analysis import plot_combined_analysis
plot_combined_analysis(
    Path("CC_analysis_MARCOIL/MARCOIL/Outputs/ProbList"),
    Path("MyProtein/MyProtein_disorder_scores.csv"),
    Path("CC_analysis_MARCOIL/MARCOIL/Outputs/Domains"),
    "MyProtein/MyProtein_combined_analysis.png",
    sequence_name="MyProtein",
    create_threshold_version=True,
    alpha_blend=0.2,
    export_csv=True,
    csv_output_dir=Path("MyProtein")
)
PY
```

```bash
MPLCONFIGDIR=$PWD/.matplotlib_cache MPLBACKEND=Agg python3 CC_analysis_MARCOIL/replot_from_csv.py MyProtein/MyProtein_plot_data.csv
```

---

## トラブルシューティング

### MARCOILが実行できない

```bash
cd MARCOIL
chmod +x marcoil
```

### パスエラーが発生する

スクリプトは`CC_analysis_MARCOIL`ディレクトリから実行してください:
```bash
cd CC_analysis_MARCOIL
python3 run_marcoil.py ...
```

### グラフが表示されない

- matplotlibがインストールされているか確認
- バックエンドエラーの場合、スクリプト冒頭に以下を追加:
  ```python
  import matplotlib
  matplotlib.use('Agg')
  ```

### Disorder Scoreのスケールが合わない

`plot_combined_analysis.py`は自動的に0-1スケールを0-100%に変換します。
すでに%表示の場合はそのまま使用されます。

---

## 入力ファイルの準備

### FASTAファイル

```fasta
>配列名
MKLAVLVLSLVGALAVGQG...
```

### Disorder Score CSV

```csv
Position,Residue,Disorder_Score,In_IDR
1,M,0.9066,Yes
2,K,0.9583,Yes
3,L,0.9599,Yes
...
```

必須カラム:
- `Position`: 残基位置（1から開始）
- `Residue`: アミノ酸1文字コード
- `Disorder_Score`: 0-1または0-100のスコア

---

## 参考情報

- **MARCOIL論文**: Delorenzi & Speed (2002) BMC Bioinformatics
- **推奨モード**:
  - `-H`: 高感度（偽陽性が多い可能性あり、スクリーニングに適する）
  - `-L`: 低感度（偽陰性が多い可能性あり、高信頼度ドメインのみ）

---

## よくある質問

**Q: 複数の配列を一度に解析できますか？**
A: MARCOILはmultiFASTAに対応していますが、プロットは配列ごとに個別に実行する必要があります。

**Q: グラフの色を変更できますか？**
A: 各スクリプト内の`color`パラメータを編集してください。

**Q: 閾値を変更したい**
A: `axhline(y=50, ...)`の値を変更してください。

**Q: Disorder ScoreがないがCC解析だけしたい**
A: `plot_coiled_coil.py`または`plot_with_moving_average.py`を使用してください。

---

## サポート

問題が発生した場合は、以下を確認してください:
1. Pythonバージョン（3.6以上推奨）
2. 必要なパッケージのインストール
3. MARCOILの実行権限
4. ファイルパスの正確性
