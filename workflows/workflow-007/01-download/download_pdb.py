# Part 1: Download PDB files
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import json
import os
import urllib.request

# =============================================================================
# CONFIGURATION FROM ENVIRONMENT VARIABLES
# =============================================================================

# PDB IDs to analyze (comma-separated, from job parameter)
# 1FME — HIV-1 protease (classic drug target, clear binding pocket)
# 2RH1 — Beta-2 adrenergic receptor (GPCR with well-defined ligand pocket)
# 3HTB — CDK2 kinase (common benchmark for pocket detection)
# 4TOS — Tankyrase 1 (PARP family, Wnt signaling drug target with clear inhibitor binding pocket)

raw_ids = os.environ.get("PARAM_PDB_IDS", "1FME, 2RH1, 4TOS")
pdb_ids = [pid.strip().upper() for pid in raw_ids.split(",") if pid.strip()]

# Output directory (silva 0.4.0+)
output_dir = "outputs"

# =============================================================================

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

print(f"Downloading {len(pdb_ids)} PDB file(s): {', '.join(pdb_ids)}", flush=True)

for pdb_id in pdb_ids:
    pdb_filename = f"{pdb_id}.pdb"
    output_path = os.path.join(output_dir, pdb_filename)
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id}...", flush=True)
    urllib.request.urlretrieve(url, output_path)
    print(f"  Saved {output_path}", flush=True)

# Write config.json for downstream jobs
config_path = os.path.join(output_dir, "config.json")
with open(config_path, "w") as f:
    json.dump({"pdb_ids": pdb_ids}, f)
print(f"Config saved to {config_path} ({len(pdb_ids)} protein(s))", flush=True)
