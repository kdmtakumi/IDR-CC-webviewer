#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
全ファイルを統合して元データと同じ順序にソート
"""
import csv

def merge_and_sort():
    """
    4つのファイルを統合し、元データの順序に並び替え
    """
    # 元データのIDリストを読み込み（順序を保持）
    print("元データの順序を読み込み中...")
    original_order = []
    with open("human_protein_ids_separated.csv", 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            original_order.append(row['UniProt_ID'])

    print(f"元データ総数: {len(original_order)} 件")

    # 全ファイルからデータを読み込み（辞書形式で格納）
    print("\n全ファイルを読み込み中...")
    all_data = {}

    files = [
        "protein_details_50k_complete.csv",
        "protein_details_50k_to_100k_complete.csv",
        "protein_details_100k_to_150k_complete.csv",
        "protein_details_150k_to_end_complete.csv",
    ]

    for filename in files:
        print(f"  {filename} を読み込み中...")
        with open(filename, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                all_data[row['UniProt_ID']] = row

    print(f"\n読み込んだデータ総数: {len(all_data)} 件")

    # 元データの順序に従ってソート
    print("\n元データの順序でソート中...")
    sorted_data = []
    missing_ids = []

    for protein_id in original_order:
        if protein_id in all_data:
            sorted_data.append(all_data[protein_id])
        else:
            missing_ids.append(protein_id)

    print(f"ソート完了: {len(sorted_data)} 件")

    if missing_ids:
        print(f"⚠️  漏れ: {len(missing_ids)} 件")
        print(f"  漏れID例: {missing_ids[:5]}")
    else:
        print("✓ 漏れなし")

    # CSV出力
    output_file = "human_protein_details_all.csv"
    fieldnames = ['UniProt_ID', 'Gene_Name', 'Protein_Name', 'Sequence_Length', 'Sequence', 'Status']

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted_data)

    print(f"\n{'='*70}")
    print(f"完了: {output_file} に保存")
    print(f"総件数: {len(sorted_data)} 件")
    print(f"成功: {sum(1 for row in sorted_data if row['Status'] == 'Success')} 件")
    print(f"{'='*70}")

if __name__ == "__main__":
    merge_and_sort()
