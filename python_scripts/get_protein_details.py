#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UniProtからタンパク質の詳細情報（配列、遺伝子名など）を取得してCSV出力
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

def get_protein_info(protein_id, session):
    """
    UniProt APIから1つのタンパク質の詳細情報を取得
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"

    try:
        response = session.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        # アミノ酸配列を取得
        sequence = data.get('sequence', {}).get('value', '')

        # 遺伝子名を取得
        gene_name = ''
        if 'genes' in data and len(data['genes']) > 0:
            gene = data['genes'][0]
            if 'geneName' in gene:
                gene_name = gene['geneName'].get('value', '')

        # タンパク質名を取得
        protein_name = ''
        if 'proteinDescription' in data:
            rec_name = data['proteinDescription'].get('recommendedName', {})
            if 'fullName' in rec_name:
                protein_name = rec_name['fullName'].get('value', '')
            else:
                # submittedNameの場合
                sub_names = data['proteinDescription'].get('submittedName', [])
                if sub_names and len(sub_names) > 0:
                    protein_name = sub_names[0].get('fullName', {}).get('value', '')

        return {
            'UniProt_ID': protein_id,
            'Gene_Name': gene_name,
            'Protein_Name': protein_name,
            'Sequence': sequence,
            'Sequence_Length': len(sequence),
            'Status': 'Success'
        }

    except requests.exceptions.Timeout:
        return {
            'UniProt_ID': protein_id,
            'Gene_Name': '',
            'Protein_Name': '',
            'Sequence': '',
            'Sequence_Length': 0,
            'Status': 'Timeout'
        }
    except Exception as e:
        return {
            'UniProt_ID': protein_id,
            'Gene_Name': '',
            'Protein_Name': '',
            'Sequence': '',
            'Sequence_Length': 0,
            'Status': f'Error: {str(e)[:50]}'
        }

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

def fetch_protein_details(protein_ids, output_file="protein_details.csv", max_workers=50):
    """
    複数のタンパク質情報を並列取得してCSVに保存
    """
    print(f"全{len(protein_ids)}件のタンパク質情報を取得します")
    print(f"{max_workers}並列で処理を開始します")
    print("=" * 60)

    results = []
    processed_count = 0
    success_count = 0

    # CSVファイルを開く（追記モード）
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'UniProt_ID', 'Gene_Name', 'Protein_Name',
            'Sequence_Length', 'Sequence', 'Status'
        ])
        writer.writeheader()

        # 並列処理
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # セッションを各スレッドで使い回す
            future_to_id = {
                executor.submit(get_protein_info, pid, get_session()): pid
                for pid in protein_ids
            }

            for future in as_completed(future_to_id):
                protein_id = future_to_id[future]
                try:
                    result = future.result()
                    results.append(result)
                    processed_count += 1

                    if result['Status'] == 'Success':
                        success_count += 1
                    else:
                        print(f"  エラー: {protein_id} - {result['Status']}")

                    # CSVに書き込み
                    writer.writerow(result)
                    f.flush()  # バッファをフラッシュして保存

                    # 進捗表示（100件ごと）
                    if processed_count % 100 == 0:
                        print(f"進捗: {processed_count}/{len(protein_ids)} 件処理完了 "
                              f"({processed_count/len(protein_ids)*100:.1f}%) | "
                              f"取得成功: {success_count} 件")

                except Exception as e:
                    print(f"  エラー: {protein_id} - {e}")
                    processed_count += 1

    print("=" * 60)
    print(f"処理完了: {processed_count}/{len(protein_ids)} 件")
    print(f"取得成功: {success_count} 件 ({success_count/len(protein_ids)*100:.1f}%)")
    print(f"結果を {output_file} に保存しました")

    return results

if __name__ == "__main__":
    print("タンパク質詳細情報の取得を開始...")
    print("=" * 60)

    # IDリストを読み込み（最初の5万件）
    input_file = "human_protein_ids_separated.csv"
    print(f"{input_file} から最初の50,000件のIDを読み込み中...")
    protein_ids = load_protein_ids(input_file, max_count=50000)
    print(f"{len(protein_ids)} 件のIDを読み込みました")

    # タンパク質情報を取得
    output_file = "protein_details_50k.csv"
    fetch_protein_details(protein_ids, output_file=output_file, max_workers=50)

    print("=" * 60)
    print("全処理完了")
