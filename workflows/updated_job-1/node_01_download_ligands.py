#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 1: Download File - Download ligand files
Input: None (implicit)
Output: ligands.zip
"""

import io
import requests

FILE_URL = (
    "https://zenodo.org/records/17374422/files/constructed_library.zip?download=1"
)
OUTPUT_FILE = "ligands.zip"

def main():
    """Main execution function"""
    print("=== Node 1: Download ligand files ===")
    print(f"Download URL: {FILE_URL}")
    
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


if __name__ == "__main__":
    main()

