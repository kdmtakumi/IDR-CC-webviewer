#!/usr/bin/env python3
"""
全チャンクファイルの細胞内局在情報をver2ファイルに統合してver3を作成
"""
import csv
from pathlib import Path

def merge_subcellular_location():
    """細胞内局在情報を統合"""

    # ファイルパス
    input_file = Path('../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation_ver2.csv')
    output_file = Path('../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation_ver3.csv')

    # チャンクファイルのリスト
    chunk_files = [
        '../all_human_protein_database_with_IDR-CCinformation_10k_with_location.csv',
        '../all_human_protein_database_10001-20000_with_location.csv',
        '../all_human_protein_database_20001-30000_with_location.csv',
        '../all_human_protein_database_30001-40000_with_location.csv',
        '../all_human_protein_database_40001-50000_with_location.csv',
        '../all_human_protein_database_50001-60000_with_location.csv',
        '../all_human_protein_database_60001-70000_with_location.csv',
        '../all_human_protein_database_70001-80000_with_location.csv',
        '../all_human_protein_database_80001-90000_with_location.csv',
        '../all_human_protein_database_90001-100000_with_location.csv',
        '../all_human_protein_database_100001-110000_with_location.csv',
        '../all_human_protein_database_110001-120000_with_location.csv',
        '../all_human_protein_database_120001-130000_with_location.csv',
        '../all_human_protein_database_130001-140000_with_location.csv',
        '../all_human_protein_database_140001-150000_with_location.csv',
        '../all_human_protein_database_150001-160000_with_location.csv',
        '../all_human_protein_database_160001-180000_with_location.csv',
        '../all_human_protein_database_180001-200000_with_location.csv',
        '../all_human_protein_database_200001-205294_with_location.csv',
    ]

    print("=" * 70)
    print("細胞内局在情報統合プログラム")
    print("=" * 70)

    # 細胞内局在情報をUniProt_IDをキーとした辞書に格納
    print("\nチャンクファイルから細胞内局在情報を読み込み中...")
    location_dict = {}

    for chunk_file in chunk_files:
        chunk_path = Path(chunk_file)
        if not chunk_path.exists():
            print(f"警告: {chunk_file} が見つかりません")
            continue

        with open(chunk_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                uniprot_id = row['UniProt_ID']
                subcellular_location = row.get('Subcellular_Location', 'N/A')
                location_dict[uniprot_id] = subcellular_location

        print(f"  読み込み完了: {chunk_path.name} ({len(location_dict)} 件)")

    print(f"\n細胞内局在情報: {len(location_dict)} 件")

    # ver2ファイルを読み込んで、細胞内局在情報を追加
    print(f"\nver2ファイルを読み込み中: {input_file}")
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames + ['Subcellular_Location']

        rows = []
        matched = 0
        not_matched = 0

        for row in reader:
            uniprot_id = row['UniProt_ID']

            # 細胞内局在情報を取得
            if uniprot_id in location_dict:
                row['Subcellular_Location'] = location_dict[uniprot_id]
                matched += 1
            else:
                row['Subcellular_Location'] = 'N/A'
                not_matched += 1

            rows.append(row)

    print(f"  読み込み完了: {len(rows)} 行")
    print(f"  マッチ: {matched} 件")
    print(f"  未マッチ: {not_matched} 件")

    # ver3ファイルに書き込み
    print(f"\nver3ファイルに書き込み中: {output_file}")
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"  書き込み完了: {len(rows)} 行")

    # 統計情報
    print("\n" + "=" * 70)
    print("統計情報")
    print("=" * 70)

    # Subcellular_Locationの内容別カウント
    location_counts = {}
    for row in rows:
        loc = row['Subcellular_Location']
        location_counts[loc] = location_counts.get(loc, 0) + 1

    print(f"\n細胞内局在情報別カウント:")
    print(f"  N/A (データなし): {location_counts.get('N/A', 0)} 件")
    print(f"  データあり: {sum(v for k, v in location_counts.items() if k != 'N/A')} 件")

    if 'Timeout' in location_counts:
        print(f"  Timeout (警告!): {location_counts['Timeout']} 件")

    print("\n" + "=" * 70)
    print("完了!")
    print("=" * 70)

if __name__ == '__main__':
    merge_subcellular_location()
