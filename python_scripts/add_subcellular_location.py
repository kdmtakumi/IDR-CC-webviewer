#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
all_human_protein_database_with_IDR-CCinformation.csvにSubcellular Location情報を追加
"""
import requests
import csv
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

def get_subcellular_location(protein_id, session):
    """
    UniProt APIから細胞内局在情報を取得
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"

    try:
        response = session.get(url, timeout=30)
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

def add_subcellular_location():
    """
    CSVファイルにSubcellular Location列を追加
    """
    input_file = "../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation.csv"
    output_file = "../IDR+CC_DB/all_human_protein_database_with_IDR-CCinformation_with_location.csv"

    print(f"ファイル読み込み中: {input_file}")

    # 既存データを読み込み
    proteins = []
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames)
        for row in reader:
            proteins.append(row)

    total_proteins = len(proteins)
    print(f"総タンパク質数: {total_proteins:,}")

    # Subcellular Location列を追加
    fieldnames.append('Subcellular_Location')

    print(f"\nUniProtから細胞内局在情報を取得中...")
    print(f"並列数: 10")
    print("=" * 70)

    # スレッドローカルなセッション
    thread_local = threading.local()

    def get_session():
        if not hasattr(thread_local, "session"):
            thread_local.session = requests.Session()
            thread_local.session.headers.update({
                'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
            })
        return thread_local.session

    # 進捗管理
    lock = threading.Lock()
    completed = 0
    success_count = 0
    timeout_count = 0

    def fetch_location(protein):
        nonlocal completed, success_count, timeout_count

        session = get_session()
        location = get_subcellular_location(protein['UniProt_ID'], session)
        protein['Subcellular_Location'] = location

        with lock:
            completed += 1
            if location not in ['N/A', 'Timeout'] and not location.startswith('Error'):
                success_count += 1
            elif location == 'Timeout':
                timeout_count += 1

            if completed % 100 == 0 or completed == total_proteins:
                print(f"進捗: {completed:,}/{total_proteins:,} ({completed/total_proteins*100:.1f}%) | "
                      f"成功: {success_count:,} | タイムアウト: {timeout_count:,}", flush=True)

        time.sleep(0.1)  # レート制限対策
        return protein

    # 並列で取得
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(fetch_location, p): p for p in proteins}

        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"  エラー: {e}")

    # CSV出力
    print(f"\nCSV出力中...")
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(proteins)

    print(f"\n{'='*70}")
    print(f"完了: {output_file} に保存")
    print(f"総件数: {total_proteins:,}")
    print(f"局在情報取得成功: {success_count:,} ({success_count/total_proteins*100:.1f}%)")
    print(f"タイムアウト: {timeout_count:,}")
    print(f"{'='*70}")

    # サンプル表示
    print(f"\n局在情報の例（最初の10件）:")
    print("-" * 90)
    for i, p in enumerate(proteins[:10]):
        loc = p['Subcellular_Location'][:60] if len(p['Subcellular_Location']) > 60 else p['Subcellular_Location']
        print(f"{p['UniProt_ID']:<15} {p['Gene_Name']:<12} {loc}")

if __name__ == "__main__":
    add_subcellular_location()
