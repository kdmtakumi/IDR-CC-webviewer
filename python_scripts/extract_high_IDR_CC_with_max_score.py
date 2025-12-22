#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IDRとCCが両方とも30%以上、かつmax scoreが両方とも0.5以上のタンパク質を抽出
"""
import csv

def extract_high_idr_cc_with_max_score():
    """
    条件:
    - IDR_Percentage >= 30%
    - CC_Percentage >= 30%
    - Max_Disorder_Score >= 0.5
    - CC_Max_Score >= 0.5
    """
    input_file = "../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation.csv"
    output_file = "../proteins_with_high_IDR_CC_and_max_scores.csv"

    print(f"ファイル読み込み中: {input_file}")

    filtered_proteins = []
    total_count = 0

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        for row in reader:
            total_count += 1

            try:
                idr_pct = float(row['IDR_Percentage'])
                cc_pct = float(row['CC_Percentage'])
                max_disorder_score = float(row['Max_Disorder_Score'])
                cc_max_score = float(row['CC_Max_Score'])

                if (idr_pct >= 30.0 and cc_pct >= 30.0 and
                    max_disorder_score >= 0.5 and cc_max_score >= 0.5):
                    filtered_proteins.append(row)

            except (ValueError, KeyError):
                # 数値に変換できない場合はスキップ
                continue

    # CSV出力
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(filtered_proteins)

    # 統計
    print(f"\n{'='*70}")
    print(f"総タンパク質数: {total_count:,}")
    print(f"抽出件数: {len(filtered_proteins):,} ({len(filtered_proteins)/total_count*100:.2f}%)")
    print(f"\n条件:")
    print(f"  IDR_Percentage >= 30%")
    print(f"  AND CC_Percentage >= 30%")
    print(f"  AND Max_Disorder_Score >= 0.5")
    print(f"  AND CC_Max_Score >= 0.5")
    print(f"\n完了: {output_file} に保存")
    print(f"{'='*70}")

    # 上位10件を表示
    if filtered_proteins:
        # IDR+CC合計が高い順にソート
        sorted_proteins = sorted(
            filtered_proteins,
            key=lambda x: float(x['IDR_Percentage']) + float(x['CC_Percentage']),
            reverse=True
        )

        print(f"\nIDR+CC合計が高い順 トップ10:")
        print("-" * 110)
        print(f"{'UniProt_ID':<15} {'Gene':<12} {'IDR%':>6} {'CC%':>6} {'IDR_Max':>8} {'CC_Max':>7} {'合計%':>7} {'タンパク質名'}")
        print("-" * 110)

        for p in sorted_proteins[:10]:
            idr_pct = float(p['IDR_Percentage'])
            cc_pct = float(p['CC_Percentage'])
            max_disorder = float(p['Max_Disorder_Score'])
            cc_max = float(p['CC_Max_Score'])
            total_pct = idr_pct + cc_pct
            protein_name = p['Protein_Name'][:35]

            print(f"{p['UniProt_ID']:<15} {p['Gene_Name']:<12} {idr_pct:6.1f} {cc_pct:6.1f} {max_disorder:8.3f} {cc_max:7.3f} {total_pct:7.1f} {protein_name}")

if __name__ == "__main__":
    extract_high_idr_cc_with_max_score()
