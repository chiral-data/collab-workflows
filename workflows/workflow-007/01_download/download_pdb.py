# Part 1: Download PDB file
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import os
import urllib.request

# =============================================================================
# CONFIGURATION FROM ENVIRONMENT VARIABLES
# =============================================================================

# PDB ID to analyze (from global workflow parameter)
pdb_id = os.environ.get("PARAM_PDB_ID", "4TOS")

# =============================================================================

# Download the pdb file
print(f"Downloading PDB file: {pdb_id}...", flush=True)
pdb_filename = f"{pdb_id.upper()}.pdb"
url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
urllib.request.urlretrieve(url, pdb_filename)
print(f"Downloaded {pdb_filename}", flush=True)
