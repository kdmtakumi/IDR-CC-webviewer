#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UniProtからヒトの全タンパク質IDを取得してCSV出力
"""
import requests
import csv
import time

def get_human_protein_ids():
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

    print("UniProt APIからヒトのタンパク質IDを取得中（全件）...")

    while True:
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

            all_ids.extend(ids)
            print(f"ページ {page}: {len(ids)} 件取得（累計: {len(all_ids)} 件）")

            # 次のページのカーソルを取得
            if 'Link' in response.headers:
                link_header = response.headers['Link']
                if 'cursor=' in link_header:
                    cursor_start = link_header.find('cursor=') + 7
                    cursor_end = link_header.find('&', cursor_start)
                    if cursor_end == -1:
                        cursor_end = link_header.find('>', cursor_start)
                    cursor = link_header[cursor_start:cursor_end]
                else:
                    break
            else:
                break

            page += 1
            time.sleep(0.1)  # API rate limit対策

        except Exception as e:
            print(f"エラー: {e}")
            break

    print(f"取得完了: {len(all_ids)} 件のタンパク質ID")
    return all_ids

def save_to_csv(protein_ids, output_file="human_protein_ids.csv"):
    """
    タンパク質IDをCSVファイルに保存
    """
    print(f"\n{output_file} に保存中...")

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['UniProt_ID'])  # ヘッダー

        for protein_id in protein_ids:
            writer.writerow([protein_id])

    print(f"保存完了: {len(protein_ids)} 件のIDを {output_file} に保存しました")

if __name__ == "__main__":
    print("ヒトのタンパク質IDリストを取得中...")
    print("=" * 60)

    # タンパク質IDを取得
    protein_ids = get_human_protein_ids()

    # CSVに保存
    save_to_csv(protein_ids)

    print("=" * 60)
    print("処理完了")
