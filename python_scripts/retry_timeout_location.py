#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Subcellular Locationのタイムアウトエントリを再取得
"""
import requests
import csv
import time

def get_subcellular_location(protein_id, session):
    """
    UniProt APIから細胞内局在情報を取得
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"

    try:
        response = session.get(url, timeout=60)  # 60秒に延長
        response.raise_for_status()

        data = response.json()

        # Subcellular Location情報を抽出
        locations = []

        comments = data.get('comments', [])
        for comment in comments:
            if comment.get('commentType') == 'SUBCELLULAR LOCATION':
                subcellular_locations = comment.get('subcellularLocations', [])
                for loc in subcellular_locations:
                    location = loc.get('location', {})
                    location_value = location.get('value', '')
                    if location_value and location_value not in locations:
                        locations.append(location_value)

        # 複数の局在がある場合はカンマ区切りで結合
        return ', '.join(locations) if locations else 'N/A'

    except requests.exceptions.Timeout:
        return 'Timeout'
    except Exception as e:
        return f'Error: {str(e)[:50]}'

def retry_timeout_location():
    """
    Timeoutエントリを再取得してCSVを更新
    """
    input_file = "../all_human_protein_database_with_IDR-CCinformation_10k_with_location.csv"
    output_file = "../all_human_protein_database_with_IDR-CCinformation_10k_with_location_complete.csv"

    print(f"ファイル読み込み中: {input_file}")

    # データを読み込み
    all_data = []
    timeout_entries = []

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        for row in reader:
            all_data.append(row)
            if row['Subcellular_Location'] == 'Timeout':
                timeout_entries.append(row)

    print(f"総データ数: {len(all_data):,}")
    print(f"タイムアウトエントリ: {len(timeout_entries)}")

    if not timeout_entries:
        print("タイムアウトエントリがありません")
        return

    print(f"\nタイムアウトエントリを再取得中...")
    print("=" * 70)

    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
    })

    retry_results = {}

    for i, entry in enumerate(timeout_entries, 1):
        protein_id = entry['UniProt_ID']
        location = get_subcellular_location(protein_id, session)
        retry_results[protein_id] = location

        status = "成功" if location not in ['N/A', 'Timeout'] else location
        print(f"{i}/{len(timeout_entries)}: {protein_id} -> {status}")

        time.sleep(1)  # レート制限対策

    # データを更新
    print("\nデータを統合中...")
    for row in all_data:
        if row['UniProt_ID'] in retry_results:
            row['Subcellular_Location'] = retry_results[row['UniProt_ID']]

    # CSV出力
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_data)

    # 統計
    success_count = sum(1 for row in all_data
                       if row['Subcellular_Location'] not in ['N/A', 'Timeout']
                       and not row['Subcellular_Location'].startswith('Error'))
    na_count = sum(1 for row in all_data if row['Subcellular_Location'] == 'N/A')
    timeout_count = sum(1 for row in all_data if row['Subcellular_Location'] == 'Timeout')

    print(f"\n{'='*70}")
    print(f"完了: {output_file} に保存")
    print(f"総件数: {len(all_data):,}")
    print(f"局在情報あり: {success_count:,} ({success_count/len(all_data)*100:.1f}%)")
    print(f"N/A（局在情報なし）: {na_count:,} ({na_count/len(all_data)*100:.1f}%)")
    print(f"タイムアウト残: {timeout_count}")
    print(f"{'='*70}")

if __name__ == "__main__":
    retry_timeout_location()
