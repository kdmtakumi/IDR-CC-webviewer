#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import requests
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

# スレッドローカルストレージ
thread_local = threading.local()

def get_session():
    if not hasattr(thread_local, "session"):
        thread_local.session = requests.Session()
    return thread_local.session

def get_subcellular_location(protein_id):
    """UniProtから細胞内局在情報を取得（90秒タイムアウト）"""
    session = get_session()
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"

    try:
        response = session.get(url, timeout=90)
        response.raise_for_status()
        data = response.json()

        # 細胞内局在情報を抽出
        locations = []
        comments = data.get('comments', [])
        for comment in comments:
            if comment.get('commentType') == 'SUBCELLULAR LOCATION':
                subcellular_locations = comment.get('subcellularLocations', [])
                for loc in subcellular_locations:
                    location_value = loc.get('location', {}).get('value', '')
                    if location_value and location_value not in locations:
                        locations.append(location_value)

        return ', '.join(locations) if locations else 'N/A'

    except requests.exceptions.Timeout:
        return 'Timeout'
    except Exception as e:
        return f'Error: {str(e)}'

def retry_timeouts_in_file(input_file):
    """ファイル内のタイムアウトエントリーを再試行（90秒タイムアウト、5並列）"""
    print(f"\n処理中: {input_file}")

    # データ読み込み
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    # タイムアウトエントリーを抽出
    timeout_indices = []
    for idx, row in enumerate(rows):
        if row.get('Subcellular_Location') == 'Timeout':
            timeout_indices.append(idx)

    total_timeouts = len(timeout_indices)
    print(f"タイムアウト件数: {total_timeouts}")

    if total_timeouts == 0:
        print("タイムアウトなし。スキップします。")
        return

    # 再試行処理
    success_count = 0
    still_timeout_count = 0
    na_count = 0

    def process_entry(idx):
        row = rows[idx]
        protein_id = row['UniProt_ID']
        location = get_subcellular_location(protein_id)
        time.sleep(0.2)
        return idx, location

    print("再試行中（タイムアウト: 90秒、並列数: 5）...")
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = {executor.submit(process_entry, idx): idx for idx in timeout_indices}

        for i, future in enumerate(as_completed(futures), 1):
            idx, location = future.result()
            rows[idx]['Subcellular_Location'] = location

            if location == 'Timeout':
                still_timeout_count += 1
            elif location == 'N/A':
                na_count += 1
            else:
                success_count += 1

            if i % 10 == 0 or i == total_timeouts:
                print(f"進捗: {i}/{total_timeouts} ({i*100//total_timeouts}%) | "
                      f"成功: {success_count} | N/A: {na_count} | タイムアウト: {still_timeout_count}")

    # 結果を書き込み
    with open(input_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n完了: {input_file}")
    print(f"成功: {success_count} | N/A: {na_count} | 残タイムアウト: {still_timeout_count}")

def main():
    files = [
        '../all_human_protein_database_60001-70000_with_location.csv',
        '../all_human_protein_database_70001-80000_with_location.csv',
        '../all_human_protein_database_80001-90000_with_location.csv',
        '../all_human_protein_database_90001-100000_with_location.csv',
        '../all_human_protein_database_100001-110000_with_location.csv'
    ]

    print("残タイムアウト分を再実行します（90秒タイムアウト、5並列）")
    print("=" * 70)

    for file in files:
        retry_timeouts_in_file(file)

    print("\n" + "=" * 70)
    print("全ファイルの処理が完了しました")

if __name__ == '__main__':
    main()
