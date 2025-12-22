#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UniProtからヒトの全タンパク質IDを取得してCSV出力
- organism_id:9606 (205,205件)
- taxonomy_id:9606のみ (差分89件)
を分けて取得
"""
import requests
import csv
import time

def get_protein_ids(query, label):
    """
    UniProt REST APIを使ってタンパク質IDリストを取得
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
    }

    all_ids = []
    cursor = None
    page_size = 500
    page = 1

    print(f"\n{label}を取得中...")
    print("-" * 60)

    while True:
        params = {
            'query': query,
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

    print(f"取得完了: {len(all_ids)} 件")
    return all_ids

def save_to_csv(ids_dict, output_file="human_protein_ids_separated.csv"):
    """
    タンパク質IDをCSVファイルに保存
    """
    print(f"\n{output_file} に保存中...")

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['UniProt_ID', 'Category'])  # ヘッダー

        for category, protein_ids in ids_dict.items():
            for protein_id in protein_ids:
                writer.writerow([protein_id, category])

    total = sum(len(ids) for ids in ids_dict.values())
    print(f"保存完了: {total} 件のIDを {output_file} に保存しました")
    for category, ids in ids_dict.items():
        print(f"  - {category}: {len(ids)} 件")

if __name__ == "__main__":
    print("ヒトのタンパク質IDリストを取得中...")
    print("=" * 60)

    # 1. organism_id:9606 のタンパク質IDを取得 (205,205件)
    organism_ids = get_protein_ids(
        query="organism_id:9606",
        label="organism_id:9606 (標準のヒトタンパク質)"
    )

    # 2. taxonomy_id:9606 のタンパク質IDを取得 (205,294件)
    taxonomy_ids = get_protein_ids(
        query="taxonomy_id:9606",
        label="taxonomy_id:9606 (全ヒトタンパク質)"
    )

    # 3. 差分を計算 (taxonomy_id:9606のみに存在するID)
    organism_set = set(organism_ids)
    taxonomy_set = set(taxonomy_ids)
    diff_ids = list(taxonomy_set - organism_set)

    print("\n" + "=" * 60)
    print("取得結果:")
    print(f"  organism_id:9606     : {len(organism_ids):>7} 件")
    print(f"  taxonomy_id:9606     : {len(taxonomy_ids):>7} 件")
    print(f"  差分 (taxonomy_idのみ): {len(diff_ids):>7} 件")
    print("=" * 60)

    # CSVに保存
    ids_dict = {
        'organism_id:9606': organism_ids,
        'taxonomy_id_only': diff_ids
    }
    save_to_csv(ids_dict)

    # 差分IDを別ファイルにも保存
    if diff_ids:
        print(f"\n差分IDを taxonomy_id_only.csv に保存中...")
        with open('taxonomy_id_only.csv', 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['UniProt_ID'])
            for protein_id in diff_ids:
                writer.writerow([protein_id])
        print(f"保存完了: {len(diff_ids)} 件")

    print("\n" + "=" * 60)
    print("処理完了")
