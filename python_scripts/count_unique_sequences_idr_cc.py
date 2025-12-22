#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IDR+CCファイルのユニーク配列数をカウント
"""
import csv
from collections import Counter

def count_unique_sequences():
    """
    proteins_with_both_IDR_and_CC.csv内のユニーク配列数を計算
    """
    input_file = "../IDR+CC_DB/proteins_with_both_IDR_and_CC.csv"

    print(f"ファイル読み込み中: {input_file}")

    sequences = []
    protein_data = []

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row['Sequence']
            sequences.append(seq)
            protein_data.append({
                'UniProt_ID': row['UniProt_ID'],
                'Gene_Name': row['Gene_Name'],
                'Sequence': seq
            })

    total_proteins = len(sequences)
    unique_sequences = len(set(sequences))

    print(f"\n{'='*70}")
    print(f"総タンパク質数: {total_proteins:,}")
    print(f"ユニーク配列数: {unique_sequences:,}")
    print(f"重複タンパク質数: {total_proteins - unique_sequences:,}")
    print(f"ユニーク率: {unique_sequences/total_proteins*100:.2f}%")
    print(f"{'='*70}")

    # 重複配列の分析
    seq_counts = Counter(sequences)
    duplicates = {seq: count for seq, count in seq_counts.items() if count > 1}

    if duplicates:
        print(f"\n重複配列パターン数: {len(duplicates):,}")
        print(f"重複配列を持つタンパク質数: {sum(duplicates.values()):,}")

        # 重複数が多い順にソート
        sorted_duplicates = sorted(duplicates.items(), key=lambda x: x[1], reverse=True)

        print(f"\n重複数トップ10:")
        print("-" * 70)
        for i, (seq, count) in enumerate(sorted_duplicates[:10], 1):
            # この配列を持つタンパク質を取得
            proteins_with_seq = [p for p in protein_data if p['Sequence'] == seq]
            gene_names = ', '.join([p['Gene_Name'][:15] for p in proteins_with_seq[:3]])
            if len(proteins_with_seq) > 3:
                gene_names += f", ... ({len(proteins_with_seq)-3}件)"

            print(f"{i:2}. 重複数: {count:3}件 | 配列長: {len(seq):5} | 遺伝子: {gene_names}")

if __name__ == "__main__":
    count_unique_sequences()
