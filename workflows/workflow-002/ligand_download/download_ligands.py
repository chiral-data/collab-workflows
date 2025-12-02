#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download File - Download ligand files
Input: None (implicit)
Output: ligands.zip
"""

import io
import os
import shutil
import requests

FILE_URL = (
    "https://zenodo.org/records/17374422/files/constructed_library.zip?download=1"
)
# Get script directory and set output directory relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "ligands.zip")

def main():
    """Main execution function"""
    print("=== Download ligand files ===")
    print(f"Download URL: {FILE_URL}")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    response = requests.get(FILE_URL, stream=True)
    response.raise_for_status()
    
    # Get file size
    total_size = int(response.headers.get("content-length", 0))
    block_size = 8192
    bytes_downloaded = 0
    update_interval_mb = 0.2  # Progress update interval (MB)
    
    print(f"Total size: {total_size / (1024 * 1024):.2f} MB")
    
    # Write to file
    with open(OUTPUT_FILE, "wb") as f:
        for chunk in response.iter_content(chunk_size=block_size):
            if chunk:
                f.write(chunk)
                bytes_downloaded += len(chunk)
                
                # Show progress
                if bytes_downloaded // (update_interval_mb * 1024 * 1024) > (
                    bytes_downloaded - len(chunk)
                ) // (update_interval_mb * 1024 * 1024):
                    progress_pct = (
                        (bytes_downloaded / total_size) * 100 if total_size > 0 else 0
                    )
                    print(
                        f"Progress: {bytes_downloaded / (1024 * 1024):.2f} MB / {total_size / (1024 * 1024):.2f} MB ({progress_pct:.1f}%)"
                    )
    
    print(f"✅ Download complete: {OUTPUT_FILE} ({bytes_downloaded / (1024 * 1024):.2f} MB)")
    
    # Move ligands.zip to ligand_unpack/input/ directory
    next_node_input_dir = os.path.join(SCRIPT_DIR, "..", "ligand_unpack", "input")
    os.makedirs(next_node_input_dir, exist_ok=True)
    target_file = os.path.join(next_node_input_dir, "ligands.zip")
    
    try:
        shutil.copy2(OUTPUT_FILE, target_file)
        print(f"✅ Copied ligands.zip to {target_file}")
    except Exception as e:
        print(f"⚠ Warning: Could not copy ligands.zip to {target_file}: {e}")


if __name__ == "__main__":
    main()

