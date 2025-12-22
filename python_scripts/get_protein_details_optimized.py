#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UniProtからタンパク質の詳細情報を高速取得（最適化版）
- バッチAPIを使用して複数IDを一度に取得
- 並列処理の最適化
"""
import requests
import csv
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# スレッドローカルストレージでセッションを管理
thread_local = threading.local()

def get_session():
    """スレッドごとにセッションを取得"""
    if not hasattr(thread_local, "session"):
        thread_local.session = requests.Session()
        thread_local.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
        })
    return thread_local.session

def get_proteins_batch(protein_ids_batch, session):
    """
    複数のタンパク質情報を一度に取得（バッチAPI使用）
    UniProtのバッチAPIは accessions パラメータで複数IDを一度に取得可能
    """
    # IDをカンマ区切りで結合
    ids_str = ','.join(protein_ids_batch)
    url = "https://rest.uniprot.org/uniprotkb/accessions"

    params = {
        'accessions': ids_str,
        'format': 'json',
        'size': len(protein_ids_batch)
    }

    results = []

    try:
        response = session.get(url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        # 結果を辞書に変換（高速検索用）
        results_dict = {}
        if 'results' in data:
            for entry in data['results']:
                protein_id = entry.get('primaryAccession', '')

                # アミノ酸配列を取得
                sequence = entry.get('sequence', {}).get('value', '')

                # 遺伝子名を取得
                gene_name = ''
                if 'genes' in entry and len(entry['genes']) > 0:
                    gene = entry['genes'][0]
                    if 'geneName' in gene:
                        gene_name = gene['geneName'].get('value', '')

                # タンパク質名を取得
                protein_name = ''
                if 'proteinDescription' in entry:
                    rec_name = entry['proteinDescription'].get('recommendedName', {})
                    if 'fullName' in rec_name:
                        protein_name = rec_name['fullName'].get('value', '')
                    else:
                        sub_names = entry['proteinDescription'].get('submittedName', [])
                        if sub_names and len(sub_names) > 0:
                            protein_name = sub_names[0].get('fullName', {}).get('value', '')

                results_dict[protein_id] = {
                    'UniProt_ID': protein_id,
                    'Gene_Name': gene_name,
                    'Protein_Name': protein_name,
                    'Sequence': sequence,
                    'Sequence_Length': len(sequence),
                    'Status': 'Success'
                }

        # リクエストしたすべてのIDについて結果を作成
        for protein_id in protein_ids_batch:
            if protein_id in results_dict:
                results.append(results_dict[protein_id])
            else:
                # 取得できなかった場合
                results.append({
                    'UniProt_ID': protein_id,
                    'Gene_Name': '',
                    'Protein_Name': '',
                    'Sequence': '',
                    'Sequence_Length': 0,
                    'Status': 'Not Found'
                })

    except requests.exceptions.Timeout:
        # タイムアウトの場合、全IDをエラーとして返す
        for protein_id in protein_ids_batch:
            results.append({
                'UniProt_ID': protein_id,
                'Gene_Name': '',
                'Protein_Name': '',
                'Sequence': '',
                'Sequence_Length': 0,
                'Status': 'Timeout'
            })
    except Exception as e:
        # エラーの場合、全IDをエラーとして返す
        for protein_id in protein_ids_batch:
            results.append({
                'UniProt_ID': protein_id,
                'Gene_Name': '',
                'Protein_Name': '',
                'Sequence': '',
                'Sequence_Length': 0,
                'Status': f'Error: {str(e)[:50]}'
            })

    return results

def load_protein_ids(input_file, max_count=50000):
    """
    CSVファイルからタンパク質IDを読み込む
    """
    protein_ids = []
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i >= max_count:
                break
            protein_ids.append(row['UniProt_ID'])
    return protein_ids

def fetch_protein_details_batch(protein_ids, output_file="protein_details.csv",
                                 batch_size=100, max_workers=20):
    """
    バッチAPIを使用して複数のタンパク質情報を並列取得
    """
    print(f"全{len(protein_ids)}件のタンパク質情報を取得します")
    print(f"バッチサイズ: {batch_size}, 並列数: {max_workers}")
    print("=" * 60)

    # バッチに分割
    batches = []
    for i in range(0, len(protein_ids), batch_size):
        batches.append(protein_ids[i:i+batch_size])

    print(f"{len(batches)}個のバッチに分割しました")

    processed_count = 0
    success_count = 0
    all_results = []

    # CSVファイルを開く
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'UniProt_ID', 'Gene_Name', 'Protein_Name',
            'Sequence_Length', 'Sequence', 'Status'
        ])
        writer.writeheader()

        # 並列処理
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_batch = {
                executor.submit(get_proteins_batch, batch, get_session()): batch
                for batch in batches
            }

            for future in as_completed(future_to_batch):
                batch = future_to_batch[future]
                try:
                    batch_results = future.result()

                    for result in batch_results:
                        all_results.append(result)
                        processed_count += 1

                        if result['Status'] == 'Success':
                            success_count += 1
                        elif result['Status'] not in ['Success', 'Not Found']:
                            print(f"  エラー: {result['UniProt_ID']} - {result['Status']}")

                        # CSVに書き込み
                        writer.writerow(result)

                    f.flush()  # バッファをフラッシュ

                    # 進捗表示（バッチごと）
                    print(f"進捗: {processed_count}/{len(protein_ids)} 件処理完了 "
                          f"({processed_count/len(protein_ids)*100:.1f}%) | "
                          f"取得成功: {success_count} 件")

                except Exception as e:
                    print(f"  バッチエラー: {e}")

    print("=" * 60)
    print(f"処理完了: {processed_count}/{len(protein_ids)} 件")
    print(f"取得成功: {success_count} 件 ({success_count/len(protein_ids)*100:.1f}%)")
    print(f"結果を {output_file} に保存しました")

    return all_results

if __name__ == "__main__":
    print("タンパク質詳細情報の取得を開始（最適化版）...")
    print("=" * 60)

    # IDリストを読み込み（最初の5万件）
    input_file = "human_protein_ids_separated.csv"
    print(f"{input_file} から最初の50,000件のIDを読み込み中...")
    protein_ids = load_protein_ids(input_file, max_count=50000)
    print(f"{len(protein_ids)} 件のIDを読み込みました")

    # タンパク質情報を取得（バッチAPI使用）
    output_file = "protein_details_50k_optimized.csv"
    fetch_protein_details_batch(
        protein_ids,
        output_file=output_file,
        batch_size=100,  # 1バッチあたり100ID
        max_workers=20   # 並列数（バッチAPIなので少なめ）
    )

    print("=" * 60)
    print("全処理完了")
