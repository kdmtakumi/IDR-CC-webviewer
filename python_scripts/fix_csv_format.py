#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CSVを標準フォーマットに修正（全フィールドを引用符で囲む）
"""
import csv

def fix_csv_format(input_file="human_protein_details_all.csv", output_file="human_protein_details_all_fixed.csv"):
    """
    CSVを標準フォーマットに変換
    - すべてのフィールドをダブルクォートで囲む
    - 改行文字を削除
    """
    print(f"修正中: {input_file}")

    with open(input_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', newline='', encoding='utf-8') as outfile:

        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames

        # QUOTE_ALL: すべてのフィールドを引用符で囲む
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
        writer.writeheader()

        count = 0
        for row in reader:
            # すべてのフィールドから改行文字を削除
            for key in row:
                if isinstance(row[key], str):
                    row[key] = row[key].replace('\n', '').replace('\r', '')

            writer.writerow(row)
            count += 1

    print(f"完了: {output_file}")
    print(f"総行数: {count}")

if __name__ == "__main__":
    fix_csv_format()
