#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UniProtからヒトのタンパク質総数を確認
"""
import requests

def check_total_count():
    """
    UniProt REST APIでヒトのタンパク質総数を取得
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
    }

    # クエリごとに総数を確認
    queries = {
        "全ヒトタンパク質": "organism_id:9606",
        "Swiss-Prot (reviewed)": "organism_id:9606 AND reviewed:true",
        "TrEMBL (unreviewed)": "organism_id:9606 AND reviewed:false",
    }

    print("UniProtのヒトタンパク質数を確認中...")
    print("=" * 60)

    for label, query in queries.items():
        params = {
            'query': query,
            'format': 'json',
            'size': 1  # 1件だけ取得（総数を知りたいだけ）
        }

        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()

            # ヘッダーから総数を取得
            if 'x-total-results' in response.headers:
                total = response.headers['x-total-results']
                print(f"{label}: {total:>10} 件")
            else:
                data = response.json()
                print(f"{label}: ヘッダー情報なし")
                print(f"  レスポンス keys: {data.keys()}")

        except Exception as e:
            print(f"{label}: エラー - {e}")

    print("=" * 60)

if __name__ == "__main__":
    check_total_count()
