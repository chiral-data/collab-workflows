#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 3: Unpack File - Extract ligand ZIP file
Input: ligands.zip
Output: ligands.sdf (SDF files in extracted directory)
"""

import os
import zipfile

INPUT_FILE = "ligands.zip"
OUTPUT_DIR = "."

def main():
    """Main execution function"""
    print("=== Node 3: Extract ligand ZIP file ===")
    
    if not os.path.exists(INPUT_FILE):
        print(f"❌ Error: {INPUT_FILE} not found.")
        exit(1)
    
    print(f"Source: {INPUT_FILE}")
    print(f"Destination: {OUTPUT_DIR}")
    
    try:
        with zipfile.ZipFile(INPUT_FILE, "r") as z:
            z.extractall(OUTPUT_DIR)
            # Show list of extracted files
            file_list = z.namelist()
            sdf_files = [f for f in file_list if f.endswith(".sdf")]
            print(f"✅ Extraction complete: {len(file_list)} files")
            if sdf_files:
                print(f"   Number of SDF files: {len(sdf_files)}")
                print(f"   First 5 files: {sdf_files[:5]}")
    except Exception as e:
        print(f"❌ Error occurred during extraction: {e}")
        exit(1)


if __name__ == "__main__":
    main()

