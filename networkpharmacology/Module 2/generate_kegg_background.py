# /// script
# requires-python = ">=3.10"
# dependencies = ["requests", "pandas"]
# ///

import requests
import pandas as pd
import os

print("Fetching mmu gene list from KEGG...")
r_list = requests.get("http://rest.kegg.jp/list/mmu")
mmu_to_symbol = {}
if r_list.status_code == 200:
    for line in r_list.text.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) >= 4:
            mmu_id = parts[0]
            # parts[3] is like "Aatk, Aatyk, mKIAA0641; apoptosis-associated tyrosine kinase"
            symbols_part = parts[3].split(';')[0]
            symbols = [s.strip() for s in symbols_part.split(',')]
            for sym in symbols:
                mmu_to_symbol[mmu_id] = sym
        elif len(parts) == 2:
            mmu_id = parts[0]
            symbols_part = parts[1].split(';')[0]
            symbols = [s.strip() for s in symbols_part.split(',')]
            for sym in symbols:
                mmu_to_symbol[mmu_id] = sym
else:
    print("Failed to fetch mmu list.")

# Reverse mapping: symbol -> mmu_id
symbol_to_mmu = {v: k for k, v in mmu_to_symbol.items()}

print("Fetching KO mapping from KEGG...")
r_link = requests.get("http://rest.kegg.jp/link/ko/mmu")
mmu_to_ko = {}
if r_link.status_code == 200:
    for line in r_link.text.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) == 2:
            mmu_id = parts[0]
            ko_id = parts[1].replace('ko:', '')
            mmu_to_ko[mmu_id] = ko_id
else:
    print("Failed to fetch KO links.")

# Load our DEGs
deg_file = "Simplified_Retina_DEGs_No_NA.txt"
if not os.path.exists(deg_file):
    print(f"Cannot find {deg_file}")
    exit(1)

df = pd.read_csv(deg_file, sep='\t')
gene_symbols = df['Gene_Name'].tolist()

# Map
output_file = "TBtools_Query2Knum.txt"
mapped_count = 0

with open(output_file, 'w') as out:
    for sym in gene_symbols:
        # Match case-insensitively just in case
        matched_mmu = None
        if sym in symbol_to_mmu:
            matched_mmu = symbol_to_mmu[sym]
        else:
            # try case-insensitive
            for s, m in symbol_to_mmu.items():
                if s.lower() == sym.lower():
                    matched_mmu = m
                    break
        
        if matched_mmu and matched_mmu in mmu_to_ko:
            ko_id = mmu_to_ko[matched_mmu]
            out.write(f"{sym}\t{ko_id}\n")
            mapped_count += 1

print(f"Background mapping created successfully! Mapped {mapped_count} out of {len(gene_symbols)} genes.")
