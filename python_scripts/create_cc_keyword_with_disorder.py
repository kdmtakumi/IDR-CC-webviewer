#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UniProtキーワード"Coiled coil"を持つヒトタンパク質のIDR情報を取得
"""
import requests
import csv
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

def get_cc_keyword_proteins():
    """
    UniProt APIからCoiled coilキーワードを持つヒトタンパク質IDを取得
    """
    print("UniProtからCoiled coilキーワードを持つタンパク質ID取得中...")

    url = "https://rest.uniprot.org/uniprotkb/search"
    query = "organism_id:9606 AND keyword:KW-0175"  # KW-0175 = Coiled coil

    all_ids = []
    next_url = None

    params = {
        'query': query,
        'format': 'list',
        'size': 500
    }

    session = requests.Session()

    while True:
        if next_url:
            response = session.get(next_url)
        else:
            response = session.get(url, params=params)

        response.raise_for_status()

        # IDを取得
        ids = response.text.strip().split('\n')
        all_ids.extend([id.strip() for id in ids if id.strip()])

        # 次のページがあるかチェック
        if 'Link' in response.headers:
            links = response.headers['Link'].split(',')
            next_url = None
            for link in links:
                if 'rel="next"' in link:
                    next_url = link[link.find('<')+1:link.find('>')]
                    break

            if not next_url:
                break
        else:
            break

        print(f"  取得中... {len(all_ids)} 件")
        time.sleep(0.5)

    print(f"完了: {len(all_ids)} 件のCoiled coilタンパク質ID取得")
    return all_ids

def merge_with_disorder_data(cc_ids, disorder_file="../human_protein_database_with_IDR-CCinformation.csv"):
    """
    Coiled coilキーワードタンパク質とIDR情報を統合
    """
    print(f"\nIDR情報読み込み中: {disorder_file}")

    cc_id_set = set(cc_ids)
    matched_proteins = []

    with open(disorder_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        for row in reader:
            if row['UniProt_ID'] in cc_id_set:
                matched_proteins.append(row)

    print(f"マッチ: {len(matched_proteins)} 件")

    return matched_proteins, fieldnames

def main():
    # Coiled coilキーワードを持つタンパク質IDを取得
    cc_ids = get_cc_keyword_proteins()

    # IDR/CC情報と統合
    matched_proteins, fieldnames = merge_with_disorder_data(cc_ids)

    # CSV出力
    output_file = "../coiled_coil_with_disorder_human_KW(CC).csv"

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(matched_proteins)

    # 統計
    has_idr = sum(1 for p in matched_proteins if p['Has_IDR'] == 'Yes')
    has_cc_domain = sum(1 for p in matched_proteins if int(p.get('Num_CC_Domains', 0) or 0) >= 1)
    has_both = sum(1 for p in matched_proteins
                   if p['Has_IDR'] == 'Yes' and int(p.get('Num_CC_Domains', 0) or 0) >= 1)

    print(f"\n{'='*70}")
    print(f"完了: {output_file} に保存")
    print(f"総件数: {len(matched_proteins):,}")
    print(f"  IDR保有: {has_idr:,} ({has_idr/len(matched_proteins)*100:.1f}%)")
    print(f"  CCドメイン保有: {has_cc_domain:,} ({has_cc_domain/len(matched_proteins)*100:.1f}%)")
    print(f"  両方保有: {has_both:,} ({has_both/len(matched_proteins)*100:.1f}%)")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
