#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proteins_with_both_IDR_and_CC.csv内のセレノプロテインと不明アミノ酸含有タンパク質を確認
"""
import csv

def check_special_aa():
    """
    UやXを含むタンパク質の確認
    """
    input_file = "../proteins_with_both_IDR_and_CC.csv"

    print(f"ファイル読み込み中: {input_file}")

    total = 0
    selenoproteins = []
    unknown_aa_proteins = []

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)

        for row in reader:
            total += 1

            # セレノプロテイン
            if row['Is_Selenoprotein'] == 'Yes':
                selenoproteins.append({
                    'UniProt_ID': row['UniProt_ID'],
                    'Gene_Name': row['Gene_Name'],
                    'Protein_Name': row['Protein_Name'],
                    'Sequence_Length': row['Sequence_Length']
                })

            # 不明アミノ酸含有
            try:
                unknown_count = int(row['Unknown_AA_Count'])
                if unknown_count > 0:
                    unknown_aa_proteins.append({
                        'UniProt_ID': row['UniProt_ID'],
                        'Gene_Name': row['Gene_Name'],
                        'Protein_Name': row['Protein_Name'],
                        'Sequence_Length': row['Sequence_Length'],
                        'Unknown_AA_Count': row['Unknown_AA_Count'],
                        'Unknown_AA_Percentage': row['Unknown_AA_Percentage']
                    })
            except (ValueError, KeyError):
                pass

    print(f"\n総タンパク質数: {total:,}")
    print(f"\nセレノプロテイン（U含む）: {len(selenoproteins)} 件")
    if selenoproteins:
        print("=" * 80)
        for p in selenoproteins:
            print(f"  {p['UniProt_ID']:<15} {p['Gene_Name']:<15} {p['Protein_Name'][:40]:<40} 長さ:{p['Sequence_Length']}")

    print(f"\n不明アミノ酸含有（X含む）: {len(unknown_aa_proteins)} 件")
    if unknown_aa_proteins:
        print("=" * 80)
        # 上位20件表示
        for p in unknown_aa_proteins[:20]:
            print(f"  {p['UniProt_ID']:<15} {p['Gene_Name']:<15} X数:{p['Unknown_AA_Count']:<5} ({p['Unknown_AA_Percentage']}%) {p['Protein_Name'][:30]}")

        if len(unknown_aa_proteins) > 20:
            print(f"\n  ... 他 {len(unknown_aa_proteins) - 20} 件")

if __name__ == "__main__":
    check_special_aa()
