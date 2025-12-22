#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
human_protein_details_all.csvにIDRとCoiled coilの全情報を統合
（セレノプロテイン、不明アミノ酸情報含む）
"""
import csv

def merge_all_data_complete():
    """
    3つのCSVファイルを統合（全情報を含む）
    - human_protein_details_all.csv (ベース)
    - human_proteins_disorder_analysis_complete.csv (IDR情報)
    - final_results_all_proteins_summary.csv (Coiled coil情報)
    """
    print("ファイル読み込み中...")

    # IDR情報を読み込み（全フィールド）
    print("  IDR情報読み込み中...")
    idr_data = {}
    with open("human_proteins_disorder_analysis_complete.csv", 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            idr_data[row['UniProt_ID']] = {
                'Is_Selenoprotein': row['Is_Selenoprotein'],
                'Processing_Status': row['Processing_Status'],
                'Has_IDR': row['Has_IDR'],
                'Num_IDRs': row['Num_IDRs'],
                'IDR_Boundaries': row['IDR_Boundaries'],
                'IDR_Residues': row['IDR_Residues'],
                'IDR_Percentage': row['IDR_Percentage'],
                'Mean_Disorder_Score': row['Mean_Disorder_Score'],
                'Max_Disorder_Score': row['Max_Disorder_Score'],
                'Unknown_AA_Count': row['Unknown_AA_Count'],
                'Unknown_AA_Percentage': row['Unknown_AA_Percentage']
            }
    print(f"    IDR情報: {len(idr_data)} 件")

    # Coiled coil情報を読み込み
    print("  Coiled coil情報読み込み中...")
    cc_data = {}
    with open("final_results_all_proteins_summary.csv", 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cc_data[row['UniProt_ID']] = {
                'Num_CC_Domains': row['Num_CC_Domains'],
                'Total_CC_Length': row['Total_CC_Length'],
                'CC_Percentage': row['CC_Percentage'],
                'Longest_CC_Domain_Length': row['Longest_Domain_Length'],
                'Mean_CC_Domain_Length': row['Mean_Domain_Length'],
                'CC_Mean_Score': row['Overall_Mean_Score'],
                'CC_Max_Score': row['Overall_Max_Score']
            }
    print(f"    Coiled coil情報: {len(cc_data)} 件")

    # ベースファイルを読み込んで統合
    print("  ベースファイル読み込み & 統合中...")
    merged_data = []
    match_idr = 0
    match_cc = 0
    selenoprotein_count = 0
    unknown_aa_count = 0

    with open("human_protein_details_all.csv", 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)

        for row in reader:
            uniprot_id = row['UniProt_ID']

            # IDR情報を追加
            if uniprot_id in idr_data:
                row.update(idr_data[uniprot_id])
                match_idr += 1

                # 統計カウント
                if idr_data[uniprot_id]['Is_Selenoprotein'] == 'Yes':
                    selenoprotein_count += 1
                if idr_data[uniprot_id]['Unknown_AA_Count'] != '0':
                    unknown_aa_count += 1
            else:
                row.update({
                    'Is_Selenoprotein': 'N/A',
                    'Processing_Status': 'N/A',
                    'Has_IDR': 'N/A',
                    'Num_IDRs': 'N/A',
                    'IDR_Boundaries': 'N/A',
                    'IDR_Residues': 'N/A',
                    'IDR_Percentage': 'N/A',
                    'Mean_Disorder_Score': 'N/A',
                    'Max_Disorder_Score': 'N/A',
                    'Unknown_AA_Count': 'N/A',
                    'Unknown_AA_Percentage': 'N/A'
                })

            # Coiled coil情報を追加
            if uniprot_id in cc_data:
                row.update(cc_data[uniprot_id])
                match_cc += 1
            else:
                row.update({
                    'Num_CC_Domains': 'N/A',
                    'Total_CC_Length': 'N/A',
                    'CC_Percentage': 'N/A',
                    'Longest_CC_Domain_Length': 'N/A',
                    'Mean_CC_Domain_Length': 'N/A',
                    'CC_Mean_Score': 'N/A',
                    'CC_Max_Score': 'N/A'
                })

            merged_data.append(row)

    print(f"    統合データ: {len(merged_data)} 件")
    print(f"    IDRマッチ: {match_idr} 件")
    print(f"    CCマッチ: {match_cc} 件")
    print(f"    セレノプロテイン: {selenoprotein_count} 件")
    print(f"    不明アミノ酸含有: {unknown_aa_count} 件")

    # CSV出力
    output_file = "human_protein_details_with_IDR_CC_complete.csv"

    fieldnames = [
        'UniProt_ID',
        'Gene_Name',
        'Protein_Name',
        'Sequence_Length',
        'Sequence',
        'Status',
        # IDR解析ステータス
        'Is_Selenoprotein',
        'Processing_Status',
        'Unknown_AA_Count',
        'Unknown_AA_Percentage',
        # IDR情報
        'Has_IDR',
        'Num_IDRs',
        'IDR_Boundaries',
        'IDR_Residues',
        'IDR_Percentage',
        'Mean_Disorder_Score',
        'Max_Disorder_Score',
        # Coiled coil情報
        'Num_CC_Domains',
        'Total_CC_Length',
        'CC_Percentage',
        'Longest_CC_Domain_Length',
        'Mean_CC_Domain_Length',
        'CC_Mean_Score',
        'CC_Max_Score'
    ]

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(merged_data)

    print(f"\n{'='*70}")
    print(f"完了: {output_file} に保存")
    print(f"総件数: {len(merged_data)}")
    print(f"{'='*70}")

if __name__ == "__main__":
    merge_all_data_complete()
