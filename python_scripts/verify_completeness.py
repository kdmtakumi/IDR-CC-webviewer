#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
全データの被り・漏れ・総数を検証
"""
import csv

def verify_completeness():
    """
    4つのCSVファイルを検証
    """
    files = [
        ("protein_details_50k_complete.csv", 1, 50000),
        ("protein_details_50k_to_100k_complete.csv", 50001, 100000),
        ("protein_details_100k_to_150k_complete.csv", 100001, 150000),
        ("protein_details_150k_to_end_complete.csv", 150001, 205294),
    ]

    # 元データの全IDを読み込み
    print("元データ（human_protein_ids_separated.csv）読み込み中...")
    original_ids = []
    with open("human_protein_ids_separated.csv", 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            original_ids.append(row['UniProt_ID'])

    print(f"元データ総数: {len(original_ids)} 件")
    print(f"元データユニーク数: {len(set(original_ids))} 件")
    print("=" * 70)

    # 各ファイルのIDを収集
    all_retrieved_ids = []

    for filename, start, end in files:
        print(f"\n{filename} を検証中...")
        print(f"  期待範囲: {start}〜{end}件目")

        file_ids = []
        success_count = 0

        with open(filename, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                file_ids.append(row['UniProt_ID'])
                if row['Status'] == 'Success':
                    success_count += 1

        expected_count = end - start + 1
        print(f"  期待件数: {expected_count} 件")
        print(f"  実際の件数: {len(file_ids)} 件")
        print(f"  成功件数: {success_count} 件")
        print(f"  ユニーク数: {len(set(file_ids))} 件")

        # 件数チェック
        if len(file_ids) != expected_count:
            print(f"  ⚠️  件数不一致！ (差分: {len(file_ids) - expected_count})")
        else:
            print(f"  ✓ 件数OK")

        # 重複チェック
        if len(file_ids) != len(set(file_ids)):
            duplicates = len(file_ids) - len(set(file_ids))
            print(f"  ⚠️  ファイル内に重複あり: {duplicates} 件")
        else:
            print(f"  ✓ ファイル内重複なし")

        # 元データとの照合
        expected_ids = original_ids[start-1:end]
        if file_ids == expected_ids:
            print(f"  ✓ 元データと完全一致")
        else:
            print(f"  ⚠️  元データと不一致")
            # 差分を確認
            file_set = set(file_ids)
            expected_set = set(expected_ids)
            missing = expected_set - file_set
            extra = file_set - expected_set
            if missing:
                print(f"    漏れ: {len(missing)} 件")
            if extra:
                print(f"    余分: {len(extra)} 件")

        all_retrieved_ids.extend(file_ids)

    # 全体の検証
    print("\n" + "=" * 70)
    print("全体の検証:")
    print(f"  全ファイル合計: {len(all_retrieved_ids)} 件")
    print(f"  ユニーク数: {len(set(all_retrieved_ids))} 件")

    # 重複チェック
    if len(all_retrieved_ids) != len(set(all_retrieved_ids)):
        duplicates = len(all_retrieved_ids) - len(set(all_retrieved_ids))
        print(f"  ⚠️  ファイル間で重複あり: {duplicates} 件")

        # 重複IDを特定
        from collections import Counter
        id_counts = Counter(all_retrieved_ids)
        dup_ids = [id for id, count in id_counts.items() if count > 1]
        print(f"  重複ID例: {dup_ids[:5]}")
    else:
        print(f"  ✓ ファイル間重複なし")

    # 元データとの完全一致チェック
    if set(all_retrieved_ids) == set(original_ids):
        print(f"  ✓ 元データと完全一致（被り・漏れなし）")
    else:
        original_set = set(original_ids)
        retrieved_set = set(all_retrieved_ids)
        missing = original_set - retrieved_set
        extra = retrieved_set - original_set

        if missing:
            print(f"  ⚠️  漏れ: {len(missing)} 件")
            print(f"    漏れID例: {list(missing)[:5]}")
        if extra:
            print(f"  ⚠️  余分: {len(extra)} 件")
            print(f"    余分ID例: {list(extra)[:5]}")

    # 順序チェック
    if all_retrieved_ids == original_ids:
        print(f"  ✓ 元データと順序も完全一致")
    else:
        print(f"  ⚠️  順序が異なる（IDは一致していても順序が違う可能性）")

    print("=" * 70)

if __name__ == "__main__":
    verify_completeness()
