#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CSVファイル内の改行文字を削除
"""
import csv

def fix_newlines(input_file="human_protein_details_all.csv", output_file="human_protein_details_all_fixed.csv"):
    """
    Sequenceフィールド内の改行文字を削除
    """
    print(f"修正中: {input_file}")

    fixed_count = 0
    total_count = 0

    with open(input_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', newline='', encoding='utf-8') as outfile:

        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            total_count += 1

            # Sequenceフィールドから改行文字を削除
            if '\n' in row['Sequence'] or '\r' in row['Sequence']:
                row['Sequence'] = row['Sequence'].replace('\n', '').replace('\r', '')
                fixed_count += 1

            # その他のフィールドからも改行を削除（念のため）
            for key in row:
                if isinstance(row[key], str):
                    row[key] = row[key].replace('\n', '').replace('\r', '')

            writer.writerow(row)

    print(f"完了: {output_file}")
    print(f"総行数: {total_count}")
    print(f"修正行数: {fixed_count}")

if __name__ == "__main__":
    fix_newlines()
