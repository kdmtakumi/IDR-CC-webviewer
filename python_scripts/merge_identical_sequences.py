#!/usr/bin/env python3
"""
同じアミノ酸配列を持つタンパク質を統合
IDR/CC情報は同じなので、UniProt_ID, Gene_Name, Protein_Nameのみを連結
"""
import csv
from pathlib import Path
from collections import defaultdict

def merge_identical_sequences():
    """同じ配列のタンパク質を統合"""

    input_file = Path('../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation_ver3.csv')
    output_file = Path('../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation_ver4.csv')

    print("=" * 70)
    print("同一配列タンパク質統合プログラム")
    print("=" * 70)

    # ver3ファイルを読み込み、配列ごとにグループ化
    print(f"\nver3ファイルを読み込み中: {input_file}")
    sequence_groups = defaultdict(list)

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        for row in reader:
            sequence = row['Sequence']
            sequence_groups[sequence].append(row)

    total_proteins = sum(len(group) for group in sequence_groups.values())
    unique_sequences = len(sequence_groups)
    duplicates = sum(1 for group in sequence_groups.values() if len(group) > 1)

    print(f"  総タンパク質数: {total_proteins}")
    print(f"  ユニーク配列数: {unique_sequences}")
    print(f"  重複配列数: {duplicates}")

    # 統合処理
    print("\n統合処理中...")
    merged_rows = []

    for sequence, group in sequence_groups.items():
        if len(group) == 1:
            # 重複なし - そのまま使用
            merged_rows.append(group[0])
        else:
            # 重複あり - 統合
            base_row = group[0].copy()

            # UniProt_ID, Gene_Name, Protein_Nameを連結
            uniprot_ids = [row['UniProt_ID'] for row in group]
            gene_names = [row['Gene_Name'] for row in group]
            protein_names = [row['Protein_Name'] for row in group]

            base_row['UniProt_ID'] = '; '.join(uniprot_ids)
            base_row['Gene_Name'] = '; '.join(gene_names)
            base_row['Protein_Name'] = '; '.join(protein_names)

            merged_rows.append(base_row)

    print(f"  統合後のタンパク質数: {len(merged_rows)}")

    # ver4ファイルに書き込み
    print(f"\nver4ファイルに書き込み中: {output_file}")
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(merged_rows)

    print(f"  書き込み完了: {len(merged_rows)} 行")

    # 統計情報
    print("\n" + "=" * 70)
    print("統計情報")
    print("=" * 70)
    print(f"統合前: {total_proteins} タンパク質")
    print(f"統合後: {len(merged_rows)} ユニーク配列")
    print(f"削減数: {total_proteins - len(merged_rows)} ({(total_proteins - len(merged_rows)) / total_proteins * 100:.2f}%)")

    # 最も重複が多い配列のトップ5
    print("\n最も重複が多い配列 (トップ5):")
    sorted_groups = sorted(sequence_groups.items(), key=lambda x: len(x[1]), reverse=True)
    for i, (seq, group) in enumerate(sorted_groups[:5], 1):
        merged_row = next((r for r in merged_rows if r['Sequence'] == seq), None)
        if merged_row:
            print(f"\n{i}. 重複数: {len(group)}")
            print(f"   配列長: {len(seq)} aa")
            print(f"   UniProt_ID: {merged_row['UniProt_ID'][:100]}...")
            print(f"   Gene_Name: {merged_row['Gene_Name'][:100]}...")

    print("\n" + "=" * 70)
    print("完了!")
    print("=" * 70)

if __name__ == '__main__':
    merge_identical_sequences()
