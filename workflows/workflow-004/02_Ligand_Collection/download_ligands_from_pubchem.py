#!/usr/bin/env python3
"""
Node 2: Ligand Collection
Downloads ligand structures from PubChem
"""

import os
import json
import shutil
import sys
import urllib.request
from pathlib import Path
from datetime import datetime

def download_pubchem_cid(cid, outdir: Path):
    """Download single ligand from PubChem"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/SDF?record_type=3d"
    outfile = outdir / f"{cid}.sdf"

    try:
        urllib.request.urlretrieve(url, outfile)
        print(f"[Ligand Collection] Downloaded CID {cid} → {outfile}")
        return {
            "cid": cid,
            "sdf_file": f"{cid}.sdf",
            "download_status": "success",
            "file_size_bytes": outfile.stat().st_size if outfile.exists() else 0
        }
    except Exception as e:
        print(f"[Ligand Collection] Failed to download CID {cid}: {e}")
        return {
            "cid": cid,
            "sdf_file": f"{cid}.sdf",
            "download_status": "failed",
            "error": str(e),
            "file_size_bytes": 0
        }


if __name__ == "__main__":
    cid_string = os.environ.get("PARAM_LIGAND_IDS", '["3672", "2662"]')
    
    try:
        # Handle case where user might pass single string instead of list
        if cid_string.startswith("["):
            cids = json.loads(cid_string)
        else:
            cids = [cid_string]
    except Exception as e:
        print(f"[Ligand Collection] Warning: Could not parse JSON, assuming single ID. Error: {e}")
        cids = [cid_string]

    outdir = Path("./")
    outdir.mkdir(parents=True, exist_ok=True)
    
    print(f"[Ligand Collection] Working directory: {outdir.resolve()}")
    print(f"[Ligand Collection] PARAM_LIGAND_IDS: {cid_string}")
    print(f"[Ligand Collection] Parsed CIDs: {cids}")
    
    # Download all ligands and collect metadata
    ligands_data = []
    for cid in cids:
        ligand_info = download_pubchem_cid(cid, outdir)
        ligands_data.append(ligand_info)
    
    # Create outputs directory
    outputs_dir = outdir / "outputs"
    outputs_dir.mkdir(exist_ok=True)
    
    success_count = sum(1 for lig in ligands_data if lig["download_status"] == "success")
    if success_count == 0:
        print("[Ligand Collection] CRITICAL ERROR: 0 ligands downloaded successfully.")
        sys.exit(1)
            
    # Copy SDF files to outputs folder for HTML visualization
    for ligand in ligands_data:
        if ligand["download_status"] == "success":
            src = outdir / ligand["sdf_file"]
            dst = outputs_dir / ligand["sdf_file"]
            if src.exists():
                shutil.copy(src, dst)
                print(f"[Ligand Collection] Copied {src.name} to outputs/")
    
    # Pass-through: Ensure receptor.pdb reaches Node 3
    for pdb in outdir.glob("*.pdb"):
        try:
            shutil.copy(pdb, outputs_dir / pdb.name)
            print(f"[Ligand Collection] ✓ Passed {pdb.name} forward to outputs/")
        except shutil.SameFileError:
            pass
        
    print(f"[Ligand Collection] ✓ Consolidated metadata and all files to {outputs_dir}")
    
    # Generate JSON output
    data = {
        "ligand_ids": cids,
        "ligands": ligands_data,
        "total_count": len(ligands_data),
        "successful_downloads": sum(1 for lig in ligands_data if lig["download_status"] == "success"),
        "failed_downloads": sum(1 for lig in ligands_data if lig["download_status"] == "failed"),
        "download_timestamp": datetime.now().isoformat()
    }
    
    json_file = outputs_dir / "ligands_metadata.json"
    with open(json_file, "w") as f:
        json.dump(data, f, indent=2)
    
    print(f"[Ligand Collection] Metadata saved to {json_file}")
    print("[Ligand Collection] ✓ Download complete")