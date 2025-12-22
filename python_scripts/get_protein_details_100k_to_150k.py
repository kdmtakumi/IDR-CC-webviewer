#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ヒトタンパク質詳細情報取得（100,001〜150,000件目）- バッチAPI最適化版
"""
import requests
import csv
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

def get_proteins_batch(protein_ids_batch, session):
    """
    バッチAPIで複数のタンパク質情報を一度に取得
    """
    ids_str = ','.join(protein_ids_batch)
    url = "https://rest.uniprot.org/uniprotkb/accessions"

    params = {
        'accessions': ids_str,
        'format': 'json',
        'size': len(protein_ids_batch)
    }

    try:
        response = session.get(url, params=params, timeout=30)
        response.raise_for_status()

        results_data = response.json()
        results = []

        if 'results' in results_data:
            for data in results_data['results']:
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

                results.append({
                    'UniProt_ID': uniprot_id,
                    'Gene_Name': gene_name,
                    'Protein_Name': recommended_name,
                    'Sequence_Length': length,
                    'Sequence': sequence,
                    'Status': 'Success'
                })

        return results

    except requests.exceptions.Timeout:
        return [{'UniProt_ID': pid, 'Gene_Name': '', 'Protein_Name': '',
                'Sequence_Length': 0, 'Sequence': '', 'Status': 'Timeout'}
                for pid in protein_ids_batch]
    except Exception as e:
        return [{'UniProt_ID': pid, 'Gene_Name': '', 'Protein_Name': '',
                'Sequence_Length': 0, 'Sequence': '', 'Status': f'Error: {str(e)}'}
                for pid in protein_ids_batch]

def get_protein_details_batch(
    input_file="human_protein_ids_separated.csv",
    output_file="protein_details_100k_to_150k.csv",
    start_index=100000,
    end_index=150000,
    batch_size=100,
    max_workers=20
):
    """
    指定範囲のタンパク質詳細情報を取得
    """
    print(f"タンパク質ID読み込み中: {input_file}")
    print(f"取得範囲: {start_index+1}〜{end_index}件目")

    # IDを読み込み
    protein_ids = []
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if start_index <= i < end_index:
                protein_ids.append(row['UniProt_ID'])

    total_proteins = len(protein_ids)
    print(f"取得対象: {total_proteins} 件")

    # バッチに分割
    batches = [protein_ids[i:i+batch_size] for i in range(0, len(protein_ids), batch_size)]
    total_batches = len(batches)
    print(f"バッチ数: {total_batches} (各バッチ {batch_size} ID)")
    print(f"並列数: {max_workers}")
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

    # 結果を格納
    all_results = []
    lock = threading.Lock()
    completed_batches = 0

    def fetch_batch(batch):
        nonlocal completed_batches
        session = get_session()
        results = get_proteins_batch(batch, session)

        with lock:
            all_results.extend(results)
            completed_batches += 1

            processed = completed_batches * batch_size
            if processed > total_proteins:
                processed = total_proteins

            if completed_batches % 10 == 0 or completed_batches == total_batches:
                success_count = sum(1 for r in all_results if r['Status'] == 'Success')
                print(f"進捗: {processed}/{total_proteins} 件処理完了 ({processed/total_proteins*100:.1f}%) | "
                      f"成功: {success_count} 件", flush=True)

        time.sleep(0.2)
        return results

    # 並列でバッチ処理
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fetch_batch, batch): batch for batch in batches}

        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"  バッチ処理エラー: {e}")

    # 結果をCSVに保存
    fieldnames = ['UniProt_ID', 'Gene_Name', 'Protein_Name', 'Sequence_Length', 'Sequence', 'Status']

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    # 統計
    success_count = sum(1 for r in all_results if r['Status'] == 'Success')
    timeout_count = sum(1 for r in all_results if r['Status'] == 'Timeout')

    print("=" * 60)
    print(f"完了: {output_file} に保存")
    print(f"処理: {len(all_results)} 件")
    print(f"成功: {success_count} 件 ({success_count/len(all_results)*100:.1f}%)")
    print(f"タイムアウト: {timeout_count} 件")
    print("=" * 60)

if __name__ == "__main__":
    get_protein_details_batch()
