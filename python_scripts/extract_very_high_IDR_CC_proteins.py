#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IDRとCCが両方とも50%以上あるタンパク質を抽出
"""
import csv

def extract_very_high_idr_cc():
    """
    IDR_Percentage >= 50% かつ CC_Percentage >= 50% のタンパク質を抽出
    """
    input_file = "../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation.csv"
    output_file = "../proteins_with_very_high_IDR_and_CC.csv"

    print(f"ファイル読み込み中: {input_file}")

    very_high_idr_cc_proteins = []
    total_count = 0

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        for row in reader:
            total_count += 1

            try:
                idr_pct = float(row['IDR_Percentage'])
                cc_pct = float(row['CC_Percentage'])

                if idr_pct >= 50.0 and cc_pct >= 50.0:
                    very_high_idr_cc_proteins.append(row)

            except (ValueError, KeyError):
                # 数値に変換できない場合はスキップ
                continue

    # CSV出力
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(very_high_idr_cc_proteins)

    # 統計
    print(f"\n{'='*70}")
    print(f"総タンパク質数: {total_count:,}")
    print(f"抽出件数: {len(very_high_idr_cc_proteins):,} ({len(very_high_idr_cc_proteins)/total_count*100:.2f}%)")
    print(f"\n条件:")
    print(f"  IDR_Percentage >= 50%")
    print(f"  AND")
    print(f"  CC_Percentage >= 50%")
    print(f"\n完了: {output_file} に保存")
    print(f"{'='*70}")

    # 上位10件を表示
    if very_high_idr_cc_proteins:
        # IDR+CC合計が高い順にソート
        sorted_proteins = sorted(
            very_high_idr_cc_proteins,
            key=lambda x: float(x['IDR_Percentage']) + float(x['CC_Percentage']),
            reverse=True
        )

        print(f"\nIDR+CC合計が高い順 トップ10:")
        print("-" * 90)
        print(f"{'UniProt_ID':<15} {'Gene':<12} {'IDR%':>6} {'CC%':>6} {'合計%':>7} {'タンパク質名'}")
        print("-" * 90)

        for p in sorted_proteins[:10]:
            idr_pct = float(p['IDR_Percentage'])
            cc_pct = float(p['CC_Percentage'])
            total_pct = idr_pct + cc_pct
            protein_name = p['Protein_Name'][:40]

            print(f"{p['UniProt_ID']:<15} {p['Gene_Name']:<12} {idr_pct:6.1f} {cc_pct:6.1f} {total_pct:7.1f} {protein_name}")

if __name__ == "__main__":
    extract_very_high_idr_cc()
