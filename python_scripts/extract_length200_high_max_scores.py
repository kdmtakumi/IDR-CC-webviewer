#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
アミノ酸長200以上、IDRとCCのmax scoreが両方とも0.5以上のタンパク質を抽出
"""
import csv

def extract_length200_high_max_scores():
    """
    条件:
    - Sequence_Length >= 200
    - Max_Disorder_Score >= 0.5
    - CC_Max_Score >= 0.5
    """
    input_file = "../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation.csv"
    output_file = "../proteins_length200_high_max_scores.csv"

    print(f"ファイル読み込み中: {input_file}")

    filtered_proteins = []
    total_count = 0

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        for row in reader:
            total_count += 1

            try:
                seq_length = int(row['Sequence_Length'])
                max_disorder_score = float(row['Max_Disorder_Score'])
                cc_max_score = float(row['CC_Max_Score'])

                if (seq_length >= 200 and
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
    print(f"  Sequence_Length >= 200")
    print(f"  AND Max_Disorder_Score >= 0.5")
    print(f"  AND CC_Max_Score >= 0.5")
    print(f"\n完了: {output_file} に保存")
    print(f"{'='*70}")

    # 統計情報
    if filtered_proteins:
        # IDRとCC両方を持つタンパク質をカウント
        has_both_idr_cc = 0
        for p in filtered_proteins:
            try:
                if p['Has_IDR'] == 'Yes' and int(p.get('Num_CC_Domains', 0) or 0) >= 1:
                    has_both_idr_cc += 1
            except:
                pass

        print(f"\n追加統計:")
        print(f"  IDRとCCドメイン両方保有: {has_both_idr_cc:,} ({has_both_idr_cc/len(filtered_proteins)*100:.1f}%)")

        # 配列長でソート
        sorted_by_length = sorted(
            filtered_proteins,
            key=lambda x: int(x['Sequence_Length']),
            reverse=True
        )

        print(f"\n配列長が長い順 トップ10:")
        print("-" * 110)
        print(f"{'UniProt_ID':<15} {'Gene':<12} {'Length':>7} {'IDR_Max':>8} {'CC_Max':>7} {'IDR%':>6} {'CC%':>6} {'タンパク質名'}")
        print("-" * 110)

        for p in sorted_by_length[:10]:
            seq_len = int(p['Sequence_Length'])
            max_disorder = float(p['Max_Disorder_Score'])
            cc_max = float(p['CC_Max_Score'])
            idr_pct = float(p.get('IDR_Percentage', 0))
            cc_pct = float(p.get('CC_Percentage', 0))
            protein_name = p['Protein_Name'][:35]

            print(f"{p['UniProt_ID']:<15} {p['Gene_Name']:<12} {seq_len:7} {max_disorder:8.3f} {cc_max:7.3f} {idr_pct:6.1f} {cc_pct:6.1f} {protein_name}")

if __name__ == "__main__":
    extract_length200_high_max_scores()
