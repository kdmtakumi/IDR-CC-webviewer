# ファイルフォーマット仕様

このドキュメントでは、解析に使用する各ファイルの形式を詳しく説明します。

## 入力ファイル

### 1. FASTAファイル（配列ファイル）

**ファイル名**: 任意（例: `sequence.fasta`, `protein.fa`）

**形式**:
```fasta
>配列名
MKLAVLVLSLVGALAVGQGYGSARTISNPEGSPSRSPRLPRSPRLGHRRTSSGGGGGTG
KTLSMENIQSLNAAYATSGPMYLSDHEGVASTTYPKGTMTLGRATNRAVYGGRVTAMGS
SPNIASAGLSHTDVLSYTDQHGGLGGSSHHHHHQVPSMLRQ
```

**ルール**:
- 1行目: `>`で始まり、配列名を記述（スペースを含まない推奨）
- 2行目以降: アミノ酸配列（1文字コード）
- 大文字・小文字どちらでも可（内部で大文字に変換）
- 改行で分割されていてもOK

**サンプル（複数配列）**:
```fasta
>Protein1
MKLAVLVLSLVGALAVGQG
>Protein2
YGSARTISNPEGSPSRSPRL
```

---

### 2. Disorder Score CSV

**ファイル名**: 任意（例: `disorder_scores.csv`）

**必須カラム**:
- `Position`: 残基位置（1から開始の整数）
- `Residue`: アミノ酸1文字コード
- `Disorder_Score`: Disorderスコア（0-1または0-100）

**オプションカラム**:
- `In_IDR`: IDR領域かどうか（Yes/No）
- その他任意のカラム

**形式例1（0-1スケール）**:
```csv
Position,Residue,Disorder_Score,In_IDR
1,M,0.9066,Yes
2,K,0.9583,Yes
3,L,0.9599,Yes
4,A,0.8234,Yes
5,V,0.2145,No
```

**形式例2（0-100スケール）**:
```csv
Position,Residue,Disorder_Score,In_IDR
1,M,90.66,Yes
2,K,95.83,Yes
3,L,95.99,Yes
```

**注意**:
- スコアが0-1の場合、自動的に100倍されてパーセント表示になります
- すでに0-100の場合はそのまま使用されます
- カラム名は大文字小文字を区別します

---

## 出力ファイル

### 1. MARCOIL結果CSV

**ファイル名**: `run_marcoil.py`の第2引数で指定（デフォルト: `marcoil_results.csv`）

**形式**:
```csv
sequence_name,threshold,domain_number,start,end,length,max_probability
Protein1,50.0,1,148,168,21,99.8
Protein1,50.0,2,226,250,25,92.5
Protein1,80.0,1,382,613,232,100.0
```

**カラムの説明**:
- `sequence_name`: 配列名
- `threshold`: 予測に使用した確率閾値（%）
- `domain_number`: ドメイン番号
- `start`: 開始残基位置
- `end`: 終了残基位置
- `length`: ドメインの長さ（残基数）
- `max_probability`: ドメイン内の最大コイルドコイル確率（%）

---

### 2. MARCOILの生出力ファイル

保存場所: `MARCOIL/Outputs/`

#### ProbList
各残基のコイルドコイル確率を記載

```
>Protein1   ## 1
YGSARTISN PEGSPSRSP...

  cc-probability in percent and best heptad phase
   1 Y   0.0 e
   2 G   0.1 f
 136 S   0.8 e
 137 M   1.5 f
 148 L  81.8 c
```

#### Domains
予測されたドメイン領域

```
>Protein1   ## 1

NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 50.0 : 9
  1. from 148 to 168 (length = 21) with max = 99.8
  2. from 226 to 250 (length = 25) with max = 92.5
```

#### ProbPerState
各隠れマルコフモデル状態の確率（詳細解析用）

---

### 3. グラフ出力（PNG）

各スクリプトで指定したファイル名で保存されます。

**推奨設定**:
- DPI: 300（高解像度）
- フォーマット: PNG
- サイズ: 14×6インチ（基本）、16×10インチ（2パネル）

---

## ファイル配置例

```
CC_analysis_MARCOIL/
├── my_protein.fasta              # 入力配列
├── my_disorder_scores.csv        # Disorderスコア
├── results.csv                   # MARCOIL結果
├── my_protein_combined.png       # 統合グラフ
└── MARCOIL/
    └── Outputs/
        ├── ProbList              # 確率データ
        ├── Domains               # ドメイン情報
        └── ProbPerState          # 状態別確率
```

---

## よくある質問

**Q: FASTAファイルの配列名にスペースを含めてもいいですか？**
A: 可能ですが、`##`より前の部分が配列名として使われます。スペースなしを推奨します。

**Q: Disorder ScoreのCSVに追加のカラムがあっても大丈夫ですか？**
A: はい、必須カラムがあれば他のカラムは無視されます。

**Q: 複数配列のFASTAファイルを使った場合、プロットはどうなりますか？**
A: MARCOILは全配列を解析しますが、プロット時に配列名を指定する必要があります。

**Q: グラフのフォーマットをPDF等に変更できますか？**
A: スクリプト内の`plt.savefig()`の拡張子を`.pdf`等に変更してください。
