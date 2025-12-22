#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IDRとCoiled coilを両方持つタンパク質を抽出
"""
import csv

def extract_idr_cc_proteins(
    input_file="../human_protein_database_with_IDR-CCinformation.csv",
    output_file="../proteins_with_both_IDR_and_CC.csv"
):
    """
    IDRとCCを両方持つタンパク質を抽出
    条件：
    - Has_IDR == "Yes"
    - Num_CC_Domains >= 1
    """
    print(f"ファイル読み込み中: {input_file}")

    idr_cc_proteins = []
    total_count = 0
    idr_only = 0
    cc_only = 0
    both = 0
    neither = 0

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        for row in reader:
            total_count += 1

            has_idr = row['Has_IDR'] == 'Yes'
            has_cc = False

            # Num_CC_Domainsをチェック
            try:
                num_cc = int(row['Num_CC_Domains'])
                has_cc = num_cc >= 1
            except (ValueError, KeyError):
                has_cc = False

            # 統計
            if has_idr and has_cc:
                both += 1
                idr_cc_proteins.append(row)
            elif has_idr:
                idr_only += 1
            elif has_cc:
                cc_only += 1
            else:
                neither += 1

    print(f"\n統計:")
    print(f"  総タンパク質数: {total_count:,}")
    print(f"  IDRのみ: {idr_only:,} ({idr_only/total_count*100:.2f}%)")
    print(f"  CCのみ: {cc_only:,} ({cc_only/total_count*100:.2f}%)")
    print(f"  両方持つ: {both:,} ({both/total_count*100:.2f}%)")
    print(f"  どちらも持たない: {neither:,} ({neither/total_count*100:.2f}%)")

    # CSV出力
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(idr_cc_proteins)

    print(f"\n{'='*70}")
    print(f"完了: {output_file} に保存")
    print(f"抽出件数: {len(idr_cc_proteins):,}")
    print(f"{'='*70}")

if __name__ == "__main__":
    extract_idr_cc_proteins()
