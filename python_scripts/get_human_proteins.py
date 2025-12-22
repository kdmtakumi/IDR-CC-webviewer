#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UniProtからヒトの全タンパク質情報を取得してCSV出力
"""
import requests
import csv
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import sys

def get_human_protein_ids(max_count=None):
    """
    UniProt REST APIを使ってヒトの全タンパク質IDリストを取得
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
    }

    all_ids = []
    cursor = None
    page_size = 500
    page = 1

    if max_count:
        print(f"UniProt APIからヒトのタンパク質IDを取得中（最大{max_count}件）...")
    else:
        print("UniProt APIからヒトのタンパク質IDを取得中（全件）...")

    while True:
        if max_count and len(all_ids) >= max_count:
            break
        params = {
            'query': 'organism_id:9606',  # 9606 = Homo sapiens
            'format': 'list',
            'size': page_size
        }

        if cursor:
            params['cursor'] = cursor

        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()

            ids = response.text.strip().split('\n')
            ids = [id.strip() for id in ids if id.strip()]

            if not ids:
                break

            # max_countを超えないように調整
            if max_count:
                remaining = max_count - len(all_ids)
                if remaining < len(ids):
                    ids = ids[:remaining]

            all_ids.extend(ids)
            print(f"ページ {page}: {len(ids)} 件取得（累計: {len(all_ids)} 件）", flush=True)

            # 目標件数に達したら終了
            if max_count and len(all_ids) >= max_count:
                break

            # Linkヘッダーから次のカーソルを取得
            link_header = response.headers.get('Link', '')
            if 'cursor=' in link_header:
                import re
                cursor_match = re.search(r'cursor=([^&>]+)', link_header)
                if cursor_match:
                    cursor = cursor_match.group(1)
                    page += 1
                    time.sleep(0.5)
                else:
                    break
            else:
                break

        except Exception as e:
            print(f"エラー: {e}")
            break

    print(f"取得完了: {len(all_ids)} 件のタンパク質ID")
    return all_ids

def get_protein_info(protein_id, session):
    """
    UniProt REST APIを使って、タンパク質の情報を取得（セッション再利用版）
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"

    try:
        response = session.get(url, timeout=10)
        response.raise_for_status()

        data = response.json()

        # ID
        primary_accession = data.get('primaryAccession', '')
        uniprotkb_id = data.get('uniProtkbId', '')
        full_id = f"{primary_accession} · {uniprotkb_id}"

        # Protein name
        protein_desc = data.get('proteinDescription', {})
        recommended_name = protein_desc.get('recommendedName', {}).get('fullName', {}).get('value', '')

        # Amino acid sequence
        sequence = data.get('sequence', {})
        aa_sequence = sequence.get('value', '')

        return {
            'ID': full_id,
            'Protein_Name': recommended_name,
            'Amino_Acid_Sequence': aa_sequence
        }

    except Exception as e:
        print(f"  エラー: {protein_id} - {e}")
        return None

def get_all_human_proteins(output_file="human_proteins.csv", max_count=None, max_workers=50):
    """
    ヒトの全タンパク質情報を取得してCSV出力（並列処理版・セッション再利用）
    """
    if max_count:
        print(f"ヒトのタンパク質IDリストを取得中（最大{max_count}件）...")
    else:
        print("ヒトのタンパク質IDリストを取得中（全件）...")
    protein_ids = get_human_protein_ids(max_count=max_count)
    print(f"\n全{len(protein_ids)} 件のタンパク質IDを取得しました")
    print(f"{max_workers}並列で処理を開始します")
    print("=" * 60)

    results = []
    lock = threading.Lock()
    completed_count = 0

    # スレッドローカルなセッションを作成
    thread_local = threading.local()

    def get_session():
        if not hasattr(thread_local, "session"):
            thread_local.session = requests.Session()
            thread_local.session.headers.update({
                'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
            })
        return thread_local.session

    def fetch_and_store(protein_id):
        nonlocal completed_count
        session = get_session()
        info = get_protein_info(protein_id, session)

        with lock:
            nonlocal completed_count
            completed_count += 1

            if info:
                results.append(info)

            # 100件ごとに進捗を表示
            if completed_count % 100 == 0:
                print(f"進捗: {completed_count}/{len(protein_ids)} 件処理完了 ({completed_count/len(protein_ids)*100:.1f}%) | 取得: {len(results)} 件", flush=True)

        return info

    # 並列処理で取得
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fetch_and_store, pid): pid for pid in protein_ids}

        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                protein_id = futures[future]
                print(f"  エラー: {protein_id} - {e}")

    # CSV出力
    if results:
        fieldnames = ['ID', 'Protein_Name', 'Amino_Acid_Sequence']

        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print("\n" + "=" * 60)
        print(f"完了: {len(results)} 件のタンパク質を {output_file} に保存しました")
        print("=" * 60)
    else:
        print("\nタンパク質が見つかりませんでした")

# 実行
if __name__ == "__main__":
    get_all_human_proteins()  # max_count=Noneで全件取得
