# /// script
# requires-python = ">=3.10"
# dependencies = ["requests", "pandas"]
# ///

import requests
import pandas as pd
import os

# 1. Download KEGG Backend (ko00001.keg)
print("Downloading KEGG Backend (ko00001.keg)...")
r = requests.get("https://rest.kegg.jp/get/br:ko00001")
if r.status_code == 200:
    with open("KEGG_Backend_ko00001.keg", "w") as f:
        f.write(r.text)
    print("Saved KEGG_Backend_ko00001.keg")
else:
    print("Failed to download ko00001")

# 2. Generate Selection Set
deg_file = "Simplified_Retina_DEGs_No_NA.txt"
if os.path.exists(deg_file):
    df = pd.read_csv(deg_file, sep='\t')
    genes = df['Gene_Name'].dropna().tolist()
    with open("TBtools_Selection_Set.txt", "w") as f:
        for g in genes:
            f.write(f"{g}\n")
    print("Saved TBtools_Selection_Set.txt")
else:
    print(f"Cannot find {deg_file}")
