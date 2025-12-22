#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
アミノ酸配列の重複チェック
"""
import csv
from collections import Counter

def check_sequence_duplicates(input_file="human_protein_details_all.csv"):
    """
    同じアミノ酸配列を持つタンパク質を検出
    """
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
                'Protein_Name': row['Protein_Name'],
                'Sequence': seq,
                'Sequence_Length': row['Sequence_Length']
            })

    print(f"総タンパク質数: {len(sequences)}")

    # 配列の出現回数をカウント
    seq_counts = Counter(sequences)

    # 重複している配列を抽出
    duplicates = {seq: count for seq, count in seq_counts.items() if count > 1}

    print(f"ユニークな配列数: {len(seq_counts)}")
    print(f"重複配列数: {len(duplicates)}")

    if duplicates:
        print(f"\n重複している配列を持つタンパク質:")
        print("=" * 80)

        # 重複回数が多い順にソート
        sorted_duplicates = sorted(duplicates.items(), key=lambda x: x[1], reverse=True)

        for seq, count in sorted_duplicates[:20]:  # 上位20件表示
            print(f"\n配列長: {len(seq)} | 重複数: {count}件")
            print("-" * 80)

            # この配列を持つタンパク質をリストアップ
            proteins_with_seq = [p for p in protein_data if p['Sequence'] == seq]

            for p in proteins_with_seq:
                print(f"  {p['UniProt_ID']:<15} {p['Gene_Name']:<15} {p['Protein_Name'][:50]}")

        if len(duplicates) > 20:
            print(f"\n... 他 {len(duplicates) - 20} 件の重複配列")

    else:
        print("\n配列の重複はありません")

    # 統計
    print("\n" + "=" * 80)
    print("統計:")
    print(f"  総タンパク質数: {len(sequences)}")
    print(f"  ユニーク配列数: {len(seq_counts)}")
    print(f"  重複配列を持つタンパク質数: {sum(duplicates.values()) if duplicates else 0}")
    print(f"  重複配列パターン数: {len(duplicates)}")
    print("=" * 80)

if __name__ == "__main__":
    check_sequence_duplicates()
