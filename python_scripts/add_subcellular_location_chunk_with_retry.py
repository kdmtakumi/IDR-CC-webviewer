#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import requests
import time
import threading
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

# スレッドローカルストレージ
thread_local = threading.local()

def get_session():
    if not hasattr(thread_local, "session"):
        thread_local.session = requests.Session()
    return thread_local.session

def get_subcellular_location(protein_id, timeout=30):
    """UniProtから細胞内局在情報を取得"""
    session = get_session()
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"

    try:
        response = session.get(url, timeout=timeout)
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

def add_subcellular_location_chunk(start_idx, end_idx):
    """指定範囲にSubcellular_Location列を追加（初回実行）"""
    input_file = '../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation.csv'
    output_file = f'../all_human_protein_database_{start_idx+1}-{end_idx}_with_location.csv'

    print(f"\nファイル読み込み中: {input_file}")
    print(f"処理範囲: {start_idx+1}〜{end_idx}件目")

    # データ読み込み
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)
        fieldnames = reader.fieldnames

    # 指定範囲を抽出
    target_rows = all_rows[start_idx:end_idx]
    print(f"処理対象: {len(target_rows):,} 件")

    # Subcellular_Location列を追加
    new_fieldnames = list(fieldnames) + ['Subcellular_Location']

    print("\nUniProtから細胞内局在情報を取得中...")
    print("並列数: 10")
    print("=" * 70)

    success_count = 0
    timeout_count = 0
    na_count = 0

    def process_protein(row):
        protein_id = row['UniProt_ID']
        location = get_subcellular_location(protein_id, timeout=30)
        time.sleep(0.1)
        return location

    # 並列処理
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(process_protein, row): idx for idx, row in enumerate(target_rows)}

        for i, future in enumerate(as_completed(futures), 1):
            idx = futures[future]
            location = future.result()
            target_rows[idx]['Subcellular_Location'] = location

            if location == 'Timeout':
                timeout_count += 1
            elif location == 'N/A':
                na_count += 1
            else:
                success_count += 1

            if i % 100 == 0 or i == len(target_rows):
                print(f"進捗: {i:,}/{len(target_rows):,} ({i*100//len(target_rows)}%) | "
                      f"成功: {success_count:,} | タイムアウト: {timeout_count:,}")

    # 結果を書き込み
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=new_fieldnames)
        writer.writeheader()
        writer.writerows(target_rows)

    print("\n" + "=" * 70)
    print(f"初回処理完了: {output_file}")
    print(f"成功: {success_count:,} | N/A: {na_count:,} | タイムアウト: {timeout_count:,}")

    return output_file, timeout_count

def retry_timeouts(output_file):
    """タイムアウトエントリーを再試行"""
    print(f"\n{'=' * 70}")
    print(f"タイムアウト分を再実行中: {output_file}")
    print("=" * 70)

    # データ読み込み
    with open(output_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    # タイムアウトエントリーを抽出
    timeout_indices = []
    for idx, row in enumerate(rows):
        if row.get('Subcellular_Location') == 'Timeout':
            timeout_indices.append(idx)

    total_timeouts = len(timeout_indices)
    print(f"タイムアウト件数: {total_timeouts:,}")

    if total_timeouts == 0:
        print("タイムアウトなし。完了しました。")
        return

    # 再試行処理（60秒タイムアウト、10並列）
    success_count = 0
    still_timeout_count = 0
    na_count = 0

    def process_entry(idx):
        row = rows[idx]
        protein_id = row['UniProt_ID']
        location = get_subcellular_location(protein_id, timeout=60)
        time.sleep(0.1)
        return idx, location

    print("再試行中（タイムアウト: 60秒、並列数: 10）...")
    with ThreadPoolExecutor(max_workers=10) as executor:
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

            if i % 100 == 0 or i == total_timeouts:
                print(f"進捗: {i:,}/{total_timeouts:,} ({i*100//total_timeouts}%) | "
                      f"成功: {success_count:,} | N/A: {na_count:,} | タイムアウト: {still_timeout_count:,}")

    # 結果を書き込み
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("\n" + "=" * 70)
    print(f"再実行完了: {output_file}")
    print(f"成功: {success_count:,} | N/A: {na_count:,} | 残タイムアウト: {still_timeout_count:,}")

def main():
    if len(sys.argv) != 3:
        print("使用方法: python add_subcellular_location_chunk_with_retry.py <start_idx> <end_idx>")
        sys.exit(1)

    start_idx = int(sys.argv[1])
    end_idx = int(sys.argv[2])

    print("=" * 70)
    print(f"細胞内局在情報取得プログラム（自動再試行付き）")
    print(f"範囲: {start_idx+1}〜{end_idx}件目")
    print("=" * 70)

    # 初回実行
    output_file, timeout_count = add_subcellular_location_chunk(start_idx, end_idx)

    # タイムアウトがあれば自動再試行
    if timeout_count > 0:
        retry_timeouts(output_file)
    else:
        print(f"\nタイムアウトなし。処理完了しました。")

    print("\n" + "=" * 70)
    print("全処理完了")
    print("=" * 70)

if __name__ == '__main__':
    main()
