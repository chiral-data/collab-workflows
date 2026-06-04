#!/usr/bin/env python3
"""
Node 01: Validate Inputs

- Validate ligand file format (.smiles/.smi or .sdf)
- Fetch receptor PDB and co-crystal ligand from RCSB
- Confirm structural resolution is below the cutoff
"""

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

# Try multiple SDF source URLs in order until one succeeds
LIGAND_SDF_URLS = [
    "https://www.rcsb.org/ccd/download?id={ligand_id}&type=ideal&format=SDF",
    "https://files.rcsb.org/ligands/{ligand_id}_ideal.sdf",
    "https://files.rcsb.org/ligands/{ligand_id}_model.sdf",
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand_id}/SDF",
]

# Non-ligand HETATM residues to skip when identifying the native ligand
EXCLUDE_RESIDUES = {
    "HOH", "WAT", "H2O",              # water
    "SO4", "PO4", "NO3",              # anions
    "CL", "NA", "MG", "ZN", "CA",    # ions
    "K", "FE", "MN", "CU", "CO",
    "GOL", "EDO", "PEG", "MPD",      # cryo/glycerol additives
    "IPA", "DMS", "ACT", "ACE",
    "MSE", "UNX", "UNL",
}


def fetch_url(url, label=""):
    """Fetch bytes from a URL; exit on failure."""
    try:
        print(f"  Fetching {label or url} ...", flush=True)
        with urllib.request.urlopen(url, timeout=30) as resp:
            return resp.read()
    except urllib.error.HTTPError as e:
        print(f"  HTTP {e.code} for {url}", flush=True)
        return None
    except Exception as e:
        print(f"  Error fetching {url}: {e}", flush=True)
        return None


# ── Receptor ──────────────────────────────────────────────────────────────────

def fetch_receptor(pdb_id):
    """Download PDB file; write receptor.pdb; return text."""
    url = RCSB_PDB_URL.format(pdb_id=pdb_id)
    data = fetch_url(url, f"PDB {pdb_id}")
    if data is None:
        print(f"ERROR: Could not download PDB {pdb_id} from RCSB", flush=True)
        sys.exit(1)
    with open("receptor.pdb", "wb") as f:
        f.write(data)
    print(f"  Saved receptor.pdb ({len(data):,} bytes)", flush=True)
    return data.decode("utf-8", errors="replace")


def get_resolution(pdb_id, pdb_text):
    """Return resolution in Angstroms, or None if not determinable."""
    # Try RCSB REST API first
    url = RCSB_ENTRY_URL.format(pdb_id=pdb_id)
    data = fetch_url(url, f"entry metadata {pdb_id}")
    if data:
        try:
            entry = json.loads(data)
            refine = entry.get("refine", [{}])
            if refine and "ls_d_res_high" in refine[0]:
                return float(refine[0]["ls_d_res_high"])
            info = entry.get("rcsb_entry_info", {})
            combo = info.get("resolution_combined")
            if combo:
                return float(combo[0])
        except (KeyError, IndexError, ValueError, TypeError):
            pass

    # Fall back to REMARK 2 in PDB text
    for line in pdb_text.splitlines():
        if line.startswith("REMARK   2 RESOLUTION."):
            m = re.search(r"RESOLUTION\.\s+([\d.]+)\s+ANGSTROMS", line)
            if m:
                return float(m.group(1))
    return None


# ── Native ligand ─────────────────────────────────────────────────────────────

def find_native_ligand(pdb_text):
    """Return the residue name of the primary HETATM ligand (most atoms)."""
    counts = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        res_name = line[17:20].strip()
        if res_name in EXCLUDE_RESIDUES:
            continue
        key = (res_name, line[21].strip(), line[22:26].strip())
        counts[key] = counts.get(key, 0) + 1

    if not counts:
        return None

    best = max(counts, key=lambda k: counts[k])
    res_name, chain, seq = best
    print(f"  Native ligand: {res_name} (chain {chain}, res {seq}, {counts[best]} atoms)", flush=True)
    return res_name


def fetch_native_ligand_sdf(ligand_id):
    """Download ligand SDF from multiple sources; write native_ligand.sdf."""
    for url_tmpl in LIGAND_SDF_URLS:
        url = url_tmpl.format(ligand_id=ligand_id)
        data = fetch_url(url, f"ligand {ligand_id} SDF")
        if data and len(data) > 50:  # guard against empty/error responses
            with open("native_ligand.sdf", "wb") as f:
                f.write(data)
            print(f"  Saved native_ligand.sdf ({len(data):,} bytes)", flush=True)
            return True

    print(f"  WARNING: Could not download SDF for ligand {ligand_id}", flush=True)
    return False


# ── Ligand file validation ────────────────────────────────────────────────────

_SMILES_CHARS = re.compile(r'^[A-Za-z0-9@+\-\[\]()\\/=#.:%*~^]+$')


def validate_smiles_file(path):
    """Validate .smiles/.smi file; write validated_ligands.smiles."""
    with open(path) as f:
        lines = f.readlines()

    valid, invalid = [], 0
    for i, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        smiles = parts[0]
        name = parts[1] if len(parts) > 1 else f"mol_{i}"
        if _SMILES_CHARS.match(smiles):
            valid.append(f"{smiles} {name}")
        else:
            print(f"  Line {i}: skipping malformed SMILES: {smiles[:50]}", flush=True)
            invalid += 1

    if not valid:
        print("ERROR: No valid SMILES found in ligand file", flush=True)
        sys.exit(1)

    with open("validated_ligands.smiles", "w") as f:
        f.write("\n".join(valid) + "\n")

    print(f"  {len(valid)} valid, {invalid} skipped", flush=True)
    return len(valid), invalid


def validate_sdf_file(path):
    """Basic SDF validation; copy to validated_ligands.sdf and write stub smiles."""
    import shutil
    with open(path) as f:
        content = f.read()

    mol_count = content.count("$$$$")
    if mol_count == 0:
        print("ERROR: SDF file contains no $$$$ terminators", flush=True)
        sys.exit(1)

    shutil.copy(path, "validated_ligands.sdf")
    with open("validated_ligands.smiles", "w") as f:
        f.write(f"# SDF input: {mol_count} molecules — see validated_ligands.sdf\n")

    print(f"  {mol_count} SDF records found", flush=True)
    return mol_count, 0


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Node 01: Validate Inputs")
    parser.add_argument("--pdb-id", default="1OKL")
    parser.add_argument("--ligand-input", default="inputs/ligands.smiles")
    parser.add_argument("--resolution-cutoff", type=float, default=3.0)
    args = parser.parse_args()

    pdb_id = args.pdb_id.strip().upper()
    ligand_path = args.ligand_input
    cutoff = args.resolution_cutoff

    print(f"PDB ID             : {pdb_id}", flush=True)
    print(f"Ligand input       : {ligand_path}", flush=True)
    print(f"Resolution cutoff  : {cutoff} A", flush=True)
    print(flush=True)

    # 1. Validate ligand file
    print("--- Step 1: Validate ligand file ---", flush=True)
    if not os.path.exists(ligand_path):
        print(f"ERROR: Ligand file not found: {ligand_path}", flush=True)
        sys.exit(1)

    ext = os.path.splitext(ligand_path)[1].lower()
    if ext in (".smiles", ".smi"):
        valid_count, invalid_count = validate_smiles_file(ligand_path)
        ligand_format = "smiles"
    elif ext == ".sdf":
        valid_count, invalid_count = validate_sdf_file(ligand_path)
        ligand_format = "sdf"
    else:
        print(f"ERROR: Unsupported format '{ext}'. Use .smiles, .smi, or .sdf", flush=True)
        sys.exit(1)

    # 2. Fetch receptor PDB
    print("\n--- Step 2: Fetch receptor from RCSB PDB ---", flush=True)
    pdb_text = fetch_receptor(pdb_id)

    # 3. Check resolution
    print("\n--- Step 3: Check structural resolution ---", flush=True)
    resolution = get_resolution(pdb_id, pdb_text)
    resolution_pass = True
    if resolution is None:
        print(f"  WARNING: Could not determine resolution for {pdb_id}; proceeding", flush=True)
    else:
        print(f"  Resolution: {resolution} A  (cutoff: {cutoff} A)", flush=True)
        if resolution >= cutoff:
            print(f"ERROR: Resolution {resolution} A exceeds cutoff {cutoff} A", flush=True)
            sys.exit(1)
        print(f"  Resolution check passed", flush=True)
        resolution_pass = True

    # 4. Identify and fetch native ligand
    print("\n--- Step 4: Fetch native co-crystal ligand ---", flush=True)
    native_ligand_id = find_native_ligand(pdb_text)
    ligand_fetched = False
    if native_ligand_id:
        ligand_fetched = fetch_native_ligand_sdf(native_ligand_id)
    else:
        print("  WARNING: No non-solvent HETATM residue found in PDB file", flush=True)
        open("native_ligand.sdf", "w").close()

    # 5. Write validation report
    report = {
        "pdb_id": pdb_id,
        "resolution_angstroms": resolution,
        "resolution_cutoff": cutoff,
        "resolution_pass": resolution_pass,
        "native_ligand_id": native_ligand_id,
        "native_ligand_sdf_fetched": ligand_fetched,
        "ligand_format": ligand_format,
        "valid_ligand_count": valid_count,
        "invalid_ligand_count": invalid_count,
    }
    with open("validation_report.json", "w") as f:
        json.dump(report, f, indent=2)
    print("\nWrote validation_report.json", flush=True)
    print(json.dumps(report, indent=2), flush=True)

    print("\nNode 01 completed", flush=True)


if __name__ == "__main__":
    main()
