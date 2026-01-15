#!/usr/bin/env python
# coding: utf-8

# In[2]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import urllib.request

# ------------------------------------------------------------------------------
# 1. Download structure
# ------------------------------------------------------------------------------


def download_pdb_file(pdb_id, output_dir="."):
    """Download PDB file"""
    print("\n=== Downloading PDB file ===")
    pdb_file = f"{protein_directory}/{pdb_id}.pdb"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        urllib.request.urlretrieve(url, pdb_file)
        print(f"Downloaded: {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"Download error: {e}")
        exit()


# ------------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    pdb_id = "4OHU"
    #pdb_id = os.getenv("PARAM_PDB_ID")
    #protein_directory = os.getenv("PARAM_PROTEIN_DIRECTORY")
    protein_directory = "molecular_docking/protein_files"
    download_pdb_file(pdb_id)

