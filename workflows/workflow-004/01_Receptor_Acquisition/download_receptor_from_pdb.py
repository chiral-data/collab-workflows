#!/usr/bin/env python3
"""
Node 1: Receptor Acquisition
Downloads PDB structure and generates metadata
"""

import os
import json
import shutil
from pathlib import Path
from datetime import datetime
from Bio.PDB import PDBList, PDBParser

def download_receptor(pdb_id: str, outdir: Path):
    """Download receptor from PDB and generate metadata"""
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.strip().lower()

    print(f"[Receptor Acquisition] Downloading receptor: {pdb_id}")
    pdbl = PDBList()

    try:
        # BioPython saves as pdbXXXX.ent or similar
        fn = pdbl.retrieve_pdb_file(pdb_id, pdir=str(outdir), file_format="pdb")
        saved = Path(fn)
        target = outdir / f"{pdb_id}.pdb"

        # Rename to clean name
        if saved.exists() and saved != target:
            shutil.move(str(saved), str(target))
            print(f"[Receptor Acquisition] Renamed {saved.name} to {target.name}")

        outputs_dir = outdir / "outputs"
        outputs_dir.mkdir(exist_ok=True)
        
        # Copy to outputs for Silva collection
        shutil.copy(str(target), str(outputs_dir / target.name))
        print(f"[Receptor Acquisition] Saved receptor to {target} and copied to {outputs_dir}")
        
        # Parse PDB to get metadata
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, str(target))
        
        num_atoms = sum(1 for _ in structure.get_atoms())
        num_residues = sum(1 for _ in structure.get_residues())
        
        # Generate JSON metadata
        data = {
            "pdb_id": pdb_id.upper(),
            "pdb_file": f"{pdb_id}.pdb",
            "download_timestamp": datetime.now().isoformat(),
            "file_size_bytes": target.stat().st_size,
            "num_atoms": num_atoms,
            "num_residues": num_residues,
            "num_chains": len(list(structure.get_chains()))
        }
        
        json_file = outputs_dir / "receptor_metadata.json"
        with open(json_file, "w") as f:
            json.dump(data, f, indent=2)
        
        print(f"[Receptor Acquisition] Metadata saved to {json_file}")
        print("[Receptor Acquisition] ✓ Download complete")
        return target
        
    except Exception as e:
        print(f"[Receptor Acquisition] ERROR: Failed to download {pdb_id}: {e}")
        raise


if __name__ == "__main__":
    pdb_id = os.environ.get("PARAM_PDB_ID", "5kir")
    output_dir = Path("./")
    download_receptor(pdb_id, output_dir)