#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
タイムアウトしたエントリを再取得してCSVに統合
"""
import requests
import csv
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

def get_protein_details_individual(protein_id, session):
    """
    個別のタンパク質IDから詳細情報を取得（タイムアウト対策で60秒）
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"

    try:
        response = session.get(url, timeout=60)
        response.raise_for_status()

        data = response.json()

        # UniProt ID
        uniprot_id = data.get('primaryAccession', '')

        # Gene name
        genes = data.get('genes', [])
        gene_name = ''
        if genes:
            gene_name = genes[0].get('geneName', {}).get('value', '')

        # Protein name
        protein_desc = data.get('proteinDescription', {})
        recommended_name = protein_desc.get('recommendedName', {}).get('fullName', {}).get('value', '')

        # Sequence
        sequence_data = data.get('sequence', {})
        sequence = sequence_data.get('value', '')
        length = sequence_data.get('length', 0)

        return {
            'UniProt_ID': uniprot_id,
            'Gene_Name': gene_name,
            'Protein_Name': recommended_name,
            'Sequence_Length': length,
            'Sequence': sequence,
            'Status': 'Success'
        }

    except requests.exceptions.Timeout:
        return {
            'UniProt_ID': protein_id,
            'Gene_Name': '',
            'Protein_Name': '',
            'Sequence_Length': 0,
            'Sequence': '',
            'Status': 'Timeout'
        }
    except Exception as e:
        return {
            'UniProt_ID': protein_id,
            'Gene_Name': '',
            'Protein_Name': '',
            'Sequence_Length': 0,
            'Sequence': '',
            'Status': f'Error: {str(e)}'
        }

def retry_timeout_entries(
    input_file="protein_details_50k_optimized.csv",
    output_file="protein_details_50k_complete.csv",
    max_workers=10
):
    """
    タイムアウトエントリを再取得してCSVを更新
    """
    print(f"入力ファイル読み込み中: {input_file}")

    # 既存データを読み込み
    all_data = []
    timeout_ids = []

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            all_data.append(row)
            if row['Status'] == 'Timeout':
                timeout_ids.append(row['UniProt_ID'])

    print(f"全データ数: {len(all_data)} 件")
    print(f"タイムアウトエントリ: {len(timeout_ids)} 件")

    if not timeout_ids:
        print("タイムアウトエントリがありません")
        return

    print(f"\n{max_workers}並列でタイムアウトエントリを再取得中...")
    print("=" * 60)

    # スレッドローカルなセッション
    thread_local = threading.local()

    def get_session():
        if not hasattr(thread_local, "session"):
            thread_local.session = requests.Session()
            thread_local.session.headers.update({
                'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
            })
        return thread_local.session

    # 再取得結果を格納
    retry_results = {}
    lock = threading.Lock()
    completed_count = 0

    def fetch_retry(protein_id):
        nonlocal completed_count
        session = get_session()
        result = get_protein_details_individual(protein_id, session)

        with lock:
            retry_results[protein_id] = result
            completed_count += 1

            if completed_count % 10 == 0 or completed_count == len(timeout_ids):
                print(f"再取得進捗: {completed_count}/{len(timeout_ids)} 件", flush=True)

        time.sleep(0.5)  # レート制限対策
        return result

    # 並列で再取得
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fetch_retry, pid): pid for pid in timeout_ids}

        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                protein_id = futures[future]
                print(f"  エラー: {protein_id} - {e}")

    # データを更新
    print("\nデータを統合中...")
    for i, row in enumerate(all_data):
        if row['UniProt_ID'] in retry_results:
            all_data[i] = retry_results[row['UniProt_ID']]

    # 成功数をカウント
    success_count = sum(1 for row in all_data if row['Status'] == 'Success')
    timeout_count = sum(1 for row in all_data if row['Status'] == 'Timeout')

    # CSV出力
    fieldnames = ['UniProt_ID', 'Gene_Name', 'Protein_Name', 'Sequence_Length', 'Sequence', 'Status']

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_data)

    print("=" * 60)
    print(f"完了: {output_file} に保存")
    print(f"成功: {success_count} 件 ({success_count/len(all_data)*100:.1f}%)")
    print(f"タイムアウト残: {timeout_count} 件")
    print("=" * 60)

if __name__ == "__main__":
    retry_timeout_entries()
