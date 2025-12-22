#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 11:23:09 2025

@author: tomofumi
"""
import requests
import csv
import time

def get_protein_details(protein_id):
    """
    UniProt REST APIを使って、タンパク質の詳細情報を取得
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
    }

    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()

        data = response.json()

        # 基本情報を取得
        primary_accession = data.get('primaryAccession', '')

        # Entry Name (Mnemonic ID)
        uniprotkb_id = data.get('uniProtkbId', '')

        # Protein name
        protein_desc = data.get('proteinDescription', {})
        recommended_name = protein_desc.get('recommendedName', {}).get('fullName', {}).get('value', '')

        # Gene name
        genes = data.get('genes', [])
        gene_name = genes[0].get('geneName', {}).get('value', '') if genes else ''

        # Status (reviewed/unreviewed)
        entry_type = data.get('entryType', '')

        # Organism
        organism = data.get('organism', {}).get('scientificName', '')

        # Sequence length (Amino acids)
        sequence = data.get('sequence', {})
        amino_acids = sequence.get('length', '')

        # Protein existence
        protein_existence = data.get('proteinExistence', '')

        # Annotation score
        annotation_score = data.get('annotationScore', '')

        # Features (ft_region) をチェック
        features = data.get('features', [])

        has_coiled_coil = False
        has_disordered = False
        disordered_positions = []

        for feature in features:
            feature_type = feature.get('type', '')
            description = str(feature.get('description', ''))

            if feature_type == 'Coiled coil' or 'coiled coil' in description.lower():
                has_coiled_coil = True

            if feature_type == 'Region' and 'disorder' in description.lower():
                has_disordered = True
                # Positionを取得
                location = feature.get('location', {})
                start = location.get('start', {}).get('value', '')
                end = location.get('end', {}).get('value', '')
                if start and end:
                    disordered_positions.append(f"{start}-{end}")

        # 両方を持つ場合のみ情報を返す
        if has_coiled_coil and has_disordered:
            # IDを "A0JNW5 · BLT3B_HUMAN" の形式にする
            full_id = f"{primary_accession} · {uniprotkb_id}"

            # Disordered positionsをカンマ区切りで結合
            disordered_str = ', '.join(disordered_positions)

            return {
                'ID': full_id,
                'Protein': recommended_name,
                'Gene': gene_name,
                'Organism': organism,
                'Amino acids': amino_acids,
                'Disordered positions': disordered_str
            }

        return None

    except Exception as e:
        print(f"  エラー: {protein_id} - {e}")
        return None

def get_protein_ids_from_keyword(start_from=2001, max_results=2000):
    """
    REST APIを使ってKW-0175のタンパク質IDリストを取得（cursorベースのページネーション使用）

    Args:
        start_from: 取得を開始する位置（デフォルト: 2001）
        max_results: 取得する件数（デフォルト: 2000）
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36'
    }

    all_ids = []
    cursor = None
    page_size = 500
    page = 1
    total_fetched = 0

    print(f"UniProt APIからKW-0175のタンパク質IDを取得中（{start_from}件目から{start_from + max_results - 1}件目まで）...")

    # start_from の位置までスキップしながらカーソルを進める
    while total_fetched < start_from + max_results - 1:
        params = {
            'query': 'keyword:KW-0175',
            'format': 'list',
            'size': page_size
        }

        # cursorパラメータを使ってページネーション
        if cursor:
            params['cursor'] = cursor

        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()

            # IDリストを取得
            ids = response.text.strip().split('\n')
            ids = [id.strip() for id in ids if id.strip()]

            if not ids:
                break

            total_fetched += len(ids)

            # start_from以降のデータのみを保存
            if total_fetched > start_from - 1:
                # 開始位置を計算
                if len(all_ids) == 0:
                    # 最初に範囲に入るページ
                    skip_count = start_from - 1 - (total_fetched - len(ids))
                    ids_to_add = ids[skip_count:]
                else:
                    ids_to_add = ids

                # max_resultsを超えないように調整
                remaining = max_results - len(all_ids)
                if remaining < len(ids_to_add):
                    ids_to_add = ids_to_add[:remaining]

                all_ids.extend(ids_to_add)
                print(f"ページ {page}: {len(ids_to_add)} 件取得（累計: {len(all_ids)} 件、全体位置: {total_fetched - len(ids) + len(ids_to_add) + start_from - 1}）")
            else:
                # スキップ中
                print(f"ページ {page}: スキップ中（{total_fetched}/{start_from - 1}）")

            # 目標件数に達したら終了
            if len(all_ids) >= max_results:
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

    print(f"取得完了: {len(all_ids)} 件のタンパク質ID（{start_from}件目から{start_from + len(all_ids) - 1}件目まで）")
    return all_ids

def scrape_coiled_coil_with_disorder(output_file="coiled_coil_with_disorder_2001-4000.csv"):
    """
    KW-0175のタンパク質からCoiled coilとDisorderedの両方を持つものを検索してCSV出力
    """
    print("KW-0175のタンパク質IDリストを取得中（2001-4000件目）...")
    protein_ids = get_protein_ids_from_keyword(start_from=2001, max_results=2000)
    print(f"\n全{len(protein_ids)} 件のタンパク質IDを取得しました")
    print("=" * 60)

    results = []

    for i, protein_id in enumerate(protein_ids, 1):
        # 実際の番号は2001から開始
        actual_number = i + 2000
        print(f"進捗: {actual_number}/{actual_number + len(protein_ids) - i} - {protein_id}")

        details = get_protein_details(protein_id)
        if details:
            results.append(details)
            print(f"  ✓ 該当: {details['ID']}")

        # 50件ごとに進捗ログを出力
        if i % 50 == 0:
            print("=" * 60)
            print(f"■ 進捗報告: {i}/{len(protein_ids)} 件処理完了 ({i/len(protein_ids)*100:.1f}%)")
            print(f"■ 該当件数: {len(results)} 件")
            print("=" * 60)

        time.sleep(0.2)

    # CSV出力
    if results:
        fieldnames = ['ID', 'Protein', 'Gene', 'Organism', 'Amino acids', 'Disordered positions']

        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print("\n" + "=" * 60)
        print(f"完了: {len(results)} 件のタンパク質を {output_file} に保存しました")
        print("=" * 60)
    else:
        print("\n該当するタンパク質が見つかりませんでした")

# 実行
scrape_coiled_coil_with_disorder()
