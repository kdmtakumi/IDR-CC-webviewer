#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 11:23:09 2025

@author: tomofumi
"""
import requests
import csv

def get_uniprot_ids_from_pdb(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}" #指定したpdb_idからuniprot_idを返す
    response = requests.get(url)
    if response.status_code != 200: #クリアしたら200、404→Not Found、500→サーバーエラー、403→アクセス拒否
        return []
    data = response.json()
    return list(data.get(pdb_id.lower(), {}).get("UniProt", {}).keys())

def get_uniprot_info(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json" #指定したuniprot_idからjson形式で詳細な情報を返す
    response = requests.get(url)
    if response.status_code != 200:
        return None

    data = response.json()
    organism = data.get("organism", {}).get("scientificName", "False")
    genes = data.get("genes", [])
    gene_name = genes[0].get("geneName", {}).get("value", "False") if genes else "False"

    protein_desc = data.get("proteinDescription", {})
    recommended_name = protein_desc.get("recommendedName", {}).get("fullName", {}).get("value", "False")

    localization_list = []
    for comment in data.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                loc_value = loc.get("location", {}).get("value")
                if loc_value:
                    localization_list.append(loc_value)

    return {
        "organism": organism,
        "gene_name": gene_name,
        "recommended_name": recommended_name,
        "localizations": list(set(localization_list)) if localization_list else ["False"]
    }

def process_csv(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ["pdb_id", "uniprot_ids", "organisms", "gene_names", "recommended_names", "localizations"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            pdb_id = row['pdb_id']
            uniprot_ids = get_uniprot_ids_from_pdb(pdb_id)

            #uniprot_idが存在しなかった場合に飛ばす
            if not uniprot_ids:
                writer.writerow({
                    "pdb_id": pdb_id,
                    "uniprot_ids": "False",
                    "organisms": "False",
                    "gene_names": "False",
                    "recommended_names": "False",
                    "localizations": "False"
                })
                continue

            #リストの初期化
            all_uids = []
            all_organisms = []
            all_genes = []
            all_names = []
            all_locs = []

            for uid in uniprot_ids:
                info = get_uniprot_info(uid)
                all_uids.append(uid)
                if info:
                    all_organisms.append(info["organism"])
                    all_genes.append(info["gene_name"])
                    all_names.append(info["recommended_name"])
                    all_locs.extend(info["localizations"])
                else:
                    all_organisms.append("False")
                    all_genes.append("False")
                    all_names.append("False")
                    all_locs.append("False")

            writer.writerow({
                "pdb_id": pdb_id,
                "uniprot_ids": "; ".join(all_uids),
                "organisms": "; ".join(set(all_organisms)),
                "gene_names": "; ".join(all_genes),
                "recommended_names": "; ".join(all_names),
                "localizations": "; ".join(set(all_locs))
            })

# 実行
process_csv("pdb_ids.csv", "output.csv")


