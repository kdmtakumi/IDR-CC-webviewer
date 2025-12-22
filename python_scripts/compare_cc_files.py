#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2つのCoiled coilファイルの共通・差分を分析
"""
import csv

def compare_files():
    """
    coiled_coil_with_disorder_human_KW(CC).csv と
    proteins_with_both_IDR_and_CC.csv の比較
    """
    # ファイル1: UniProtキーワード"Coiled coil"を持つタンパク質
    file1 = "../final/coiled_coil_with_disorder_human_KW(CC).csv"
    # ファイル2: IDRとCCを両方持つタンパク質
    file2 = "../final/proteins_with_both_IDR_and_CC.csv"

    print("ファイル読み込み中...")

    # ファイル1のID読み込み
    ids1 = set()
    with open(file1, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # 古い形式と新しい形式に対応
            uid = row.get('UniProt_ID') or row.get('ID')
            if uid:
                # "Q5T1B0 · AXDN1_HUMAN" -> "Q5T1B0" の形式に正規化
                uid = uid.split('·')[0].strip()
                ids1.add(uid)

    # ファイル2のID読み込み
    ids2 = set()
    with open(file2, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ids2.add(row['UniProt_ID'])

    # 集合演算
    common = ids1 & ids2  # 共通
    only_in_kw = ids1 - ids2  # KWのみ
    only_in_both = ids2 - ids1  # IDR+CCのみ

    print(f"\n{'='*70}")
    print(f"ファイル1: coiled_coil_with_disorder_human_KW(CC).csv")
    print(f"  UniProtキーワード'Coiled coil'を持つタンパク質")
    print(f"  件数: {len(ids1):,}")

    print(f"\nファイル2: proteins_with_both_IDR_and_CC.csv")
    print(f"  IDRとCCドメインを両方持つタンパク質")
    print(f"  件数: {len(ids2):,}")

    print(f"\n{'='*70}")
    print(f"共通するタンパク質: {len(common):,}")
    print(f"  (両ファイルに存在)")

    print(f"\nKWファイルのみに存在: {len(only_in_kw):,}")
    print(f"  (UniProtキーワードは持つが、CCドメインまたはIDRがない)")

    print(f"\nIDR+CCファイルのみに存在: {len(only_in_both):,}")
    print(f"  (CCドメインとIDRは持つが、UniProtキーワードはない)")

    # パーセンテージ
    print(f"\n{'='*70}")
    print(f"KWファイルの観点:")
    print(f"  共通率: {len(common)/len(ids1)*100:.1f}% ({len(common):,}/{len(ids1):,})")
    print(f"  KW専有率: {len(only_in_kw)/len(ids1)*100:.1f}% ({len(only_in_kw):,}/{len(ids1):,})")

    print(f"\nIDR+CCファイルの観点:")
    print(f"  共通率: {len(common)/len(ids2)*100:.1f}% ({len(common):,}/{len(ids2):,})")
    print(f"  IDR+CC専有率: {len(only_in_both)/len(ids2)*100:.1f}% ({len(only_in_both):,}/{len(ids2):,})")

    print(f"{'='*70}")

    # サンプル表示
    if only_in_kw:
        print(f"\nKWのみの例（最初の10件）:")
        sample = list(only_in_kw)[:10]
        for uid in sample:
            print(f"  {uid}")

    if only_in_both:
        print(f"\nIDR+CCのみの例（最初の10件）:")
        sample = list(only_in_both)[:10]
        for uid in sample:
            print(f"  {uid}")

if __name__ == "__main__":
    compare_files()
