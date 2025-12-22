#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UniProtからヒトのタンパク質総数を詳細確認
"""
import requests

def check_detailed_count():
    """
    UniProt REST APIでヒトのタンパク質総数を詳細確認
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
    }

    # より詳細なクエリで確認
    queries = {
        "全ヒトタンパク質 (organism_id)": "organism_id:9606",
        "全ヒトタンパク質 (taxonomy_id)": "taxonomy_id:9606",
        "reviewed + unreviewed": "(organism_id:9606 AND reviewed:true) OR (organism_id:9606 AND reviewed:false)",
        "アクティブなエントリのみ": "organism_id:9606 AND active:true",
        "全エントリ（削除済み含む？）": "organism_id:9606",
        "フラグメント含む": "organism_id:9606 AND fragment:*",
        "フラグメント除外": "organism_id:9606 NOT fragment:*",
    }

    print("UniProtのヒトタンパク質数を詳細確認中...")
    print("=" * 70)

    for label, query in queries.items():
        params = {
            'query': query,
            'format': 'json',
            'size': 1
        }

        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()

            if 'x-total-results' in response.headers:
                total = response.headers['x-total-results']
                print(f"{label:40}: {total:>10} 件")
            else:
                print(f"{label:40}: ヘッダー情報なし")

        except Exception as e:
            print(f"{label:40}: エラー - {e}")

    print("=" * 70)
    print("\nWebサイト表示: 205,294 件")
    print("API経由取得:   205,205 件")
    print("差分:              89 件")

if __name__ == "__main__":
    check_detailed_count()
