# クイックスタートガイド

このフォルダで同じ解析を実行する最短手順です。

## 準備

1. このフォルダ(`CC_analysis_MARCOIL`)に移動
```bash
cd CC_analysis_MARCOIL
```

2. 必要なパッケージをインストール
```bash
pip install matplotlib numpy
```

## 最速実行手順

### パターン1: Disorder Score + Coiled-Coil統合解析

```bash
# 1. FASTAファイルを準備（your_sequence.fasta）

# 2. Disorder Score CSVを準備（your_disorder_scores.csv）
#    形式: Position,Residue,Disorder_Score,In_IDR

# 3. MARCOIL解析
python3 run_marcoil.py your_sequence.fasta results.csv H

# 4. plot_combined_analysis.pyを編集
# 最後の部分を以下のように変更:
#   disorder_csv = Path('your_disorder_scores.csv')
#   '配列名' を FASTAファイル内の>の後の名前に変更
#   出力ファイル名を指定

# 5. 統合グラフ作成
python3 plot_combined_analysis.py
```

### パターン2: Coiled-Coilのみの解析

```bash
# 1. FASTAファイルを準備（your_sequence.fasta）

# 2. MARCOIL解析
python3 run_marcoil.py your_sequence.fasta results.csv H

# 3. plot_with_moving_average.pyを編集
# 最後の部分で配列名と出力ファイル名を指定

# 4. グラフ作成
python3 plot_with_moving_average.py
```

## サンプルファイルを使った動作確認

```bash
# サンプルを使ってすぐに試せます
python3 run_marcoil.py example_input.fasta test_results.csv H

# 統合グラフを作成（サンプルのDisorderスコア使用）
python3 plot_combined_analysis.py
# → ELKS2_SS306_combined_analysis.png が生成されます
```

## ファイル編集が必要な箇所

### plot_combined_analysis.py の編集箇所

ファイルの最後（`if __name__ == "__main__":`以降）:

```python
if __name__ == "__main__":
    from pathlib import Path
    
    problist_file = Path('MARCOIL/Outputs/ProbList')
    domains_file = Path('MARCOIL/Outputs/Domains')
    disorder_csv = Path('あなたのDisorderファイル.csv')  # ← 変更
    
    plot_combined_analysis(
        problist_file,
        disorder_csv,
        domains_file,
        '出力ファイル名.png',  # ← 変更
        'FASTAの配列名'        # ← 変更（>の後の名前）
    )
```

### plot_with_moving_average.py の編集箇所

ファイルの最後:

```python
if __name__ == "__main__":
    from pathlib import Path
    
    problist_file = Path('MARCOIL/Outputs/ProbList')
    domains_file = Path('MARCOIL/Outputs/Domains')
    
    plot_with_moving_average(
        problist_file,
        domains_file,
        '出力ファイル名.png',  # ← 変更
        'FASTAの配列名'        # ← 変更
    )
```

## 完了！

これで同じ解析が実行できます。詳細は`README.md`を参照してください。
