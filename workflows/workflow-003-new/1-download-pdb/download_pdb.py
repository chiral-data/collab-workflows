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
    pdb_file = f"./{pdb_id}.pdb"
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
    # pdb_id = "4OHU"
    pdb_id = os.getenv("PARAM_PDB_ID")
    download_pdb_file(pdb_id)
