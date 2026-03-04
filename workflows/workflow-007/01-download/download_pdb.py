# Part 1: Download PDB file
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import json
import os
import urllib.request

# =============================================================================
# CONFIGURATION FROM ENVIRONMENT VARIABLES
# =============================================================================

# PDB ID to analyze (from job parameter)
pdb_id = os.environ.get("PARAM_PDB_ID", "4TOS")

# Output directory (silva 0.4.0+)
output_dir = "outputs"

# =============================================================================

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Download the pdb file
print(f"Downloading PDB file: {pdb_id}...", flush=True)
pdb_filename = f"{pdb_id.upper()}.pdb"
output_path = os.path.join(output_dir, pdb_filename)
url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
urllib.request.urlretrieve(url, output_path)
print(f"Downloaded {output_path}", flush=True)

# Write config.json for downstream jobs
config_path = os.path.join(output_dir, "config.json")
with open(config_path, "w") as f:
    json.dump({"pdb_id": pdb_id.upper()}, f)
print(f"Config saved to {config_path}", flush=True)
