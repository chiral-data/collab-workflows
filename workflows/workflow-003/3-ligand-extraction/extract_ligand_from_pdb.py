import os
import requests
from Bio.PDB import PDBIO, PDBParser, Select
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize


def extract_ligand_from_pdb(pdb_path, ligand_resname, out_pdb):
    """
    Extracts the 2TK ligand from a PDB file (e.g., 4OHU) and saves it as a PDB file,
    ensuring aromaticity is preserved.
    """
    
    # Load the PDB file (keep hydrogens and skip sanitization during the initial load)
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError("Could not load PDB file.")
    
    # Sanitize the molecule to ensure proper aromaticity and bond assignments
    Chem.SanitizeMol(mol)

    # Manually assign aromaticity to the molecule (this step is critical)
    rdmolops.AssignStereochemistry(mol, force=True, cleanIt=True)  # Assigns stereochemistry and aromaticity

    # Split the PDB into fragments (proteins, ligands, water, etc.)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

    ligand_frags = []

    # Iterate over fragments and extract the 2TK ligand
    for frag in frags:
        resnames = set()
        for atom in frag.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info:
                resnames.add(info.GetResidueName().strip())

        # Check if this fragment corresponds to the 2TK ligand (by residue name)
        if ligand_resname.upper() in resnames:
            ligand_frags.append(frag)

    if not ligand_frags:
        raise ValueError(f"No ligand with residue name {ligand_name} found.")

    # Write the extracted ligand to a PDB file
    writer = Chem.PDBWriter(out_pdb)
    for lig in ligand_frags:
        writer.write(lig)
    writer.close()

    print(f"Saved 2TK ligand to {out_pdb}")
    return ligand_frags

def fix_and_align():
    
    print("Fixing and correcting the pose of extracted ligand...")
    # Load and remove the hydrogens that we can't map anyways yet.
    ideal_mol = Chem.MolFromMolFile(f"{ligand_name}_ideal.sdf", removeHs=True)
    pose_mol = Chem.MolFromPDBFile(f"{ligand_name}_fromPDB.pdb", removeHs=True)
    
    # Disconnect any organometal
    rdMolStandardize.DisconnectOrganometallicsInPlace(pose_mol)

    # Remove disconnected fragments
    fragmenter = rdMolStandardize.FragmentRemover()
    pose_mol_f = fragmenter.remove(pose_mol)

    # Choose largest fragment
    chooser = rdMolStandardize.LargestFragmentChooser()
    pose_mol_lf = chooser.choose(pose_mol_f)

    # Assign bond orders from the template to the pose molecule
    corrected_pose = AllChem.AssignBondOrdersFromTemplate(ideal_mol, pose_mol_lf)

    # Add hydrogens
    corrected_pose_with_H = Chem.AddHs(corrected_pose, addCoords=True)

    # Save the corrected pose to an SDF file
    ligand_corrected_pose_file = f"{ligand_name}_corrected_pose.sdf"
    writer = Chem.SDWriter(ligand_corrected_pose_file)
    writer.write(corrected_pose_with_H)
    writer.close()
    print("Extracted ligand fixed and aligned!")

def download_ideal_ligand():

    ideal_ligand_filename = f"{ligand_name}_ideal.sdf"
    print(f"Downloading ligand {ligand_name}...")
    ligand_url = f"https://files.rcsb.org/ligands/download/{ideal_ligand_filename}"
    ligand_request = requests.get(ligand_url)
    ligand_request.raise_for_status() # Check for errors

    ideal_filepath = f"{ideal_ligand_filename}"

    with open(ideal_filepath, "w") as f:
        f.write(ligand_request.text)
    print(f"Saved ligand to {ideal_filepath}")

# Example usage: Extracting 2TK from 4OHU PDB while preserving aromaticity and saving as PDB
#extract_ligand_from_pdb(f"{pdb_id}.pdb", ligand_resname=ligand_name, out_pdb=f"{ligand_name}_fromPDB.pdb")

if __name__ == "__main__":
    # pdb_id = "4OHU"
    pdb_id = os.getenv("PARAM_PDB_ID")
    # ligand_name = "2TK"
    ligand_name = os.getenv("PARAM_LIGAND_NAME")
    extract_ligand_from_pdb(f"{pdb_id}.pdb", ligand_name, out_pdb=f"{ligand_name}_fromPDB.pdb")
    download_ideal_ligand()
    fix_and_align()
