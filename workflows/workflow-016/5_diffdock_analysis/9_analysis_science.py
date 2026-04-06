#!/usr/bin/env python3
"""
Node 5: Interface Analysis for Antibody-Antigen Complex

Analyzes the best-scoring docked pose to identify the antibody-antigen binding interface.
Computes atomic contacts, classifies interaction types, and generates detailed reports.
"""

import os
import sys
import json
import argparse
import warnings
from collections import defaultdict
import numpy as np

warnings.filterwarnings('ignore')

# Amino acid classification
HYDROPHOBIC_RESIDUES = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'}
POSITIVE_RESIDUES = {'LYS', 'ARG', 'HIS'}
NEGATIVE_RESIDUES = {'ASP', 'GLU'}
POLAR_RESIDUES = {'SER', 'THR', 'ASN', 'GLN', 'CYS', 'TYR'}


class PDBAtom:
    """Simple PDB atom representation."""
    def __init__(self, line):
        self.raw_line = line
        try:
            self.atom_num = int(line[6:11].strip())
            self.atom_name = line[12:16].strip()
            self.res_name = line[17:20].strip()
            self.chain_id = line[21] if len(line) > 21 else 'A'
            
            # Try to parse residue number (might be malformed)
            try:
                self.res_num = int(line[22:26].strip())
            except (ValueError, IndexError):
                self.res_num = None
            
            self.x = float(line[30:38].strip())
            self.y = float(line[38:46].strip())
            self.z = float(line[46:54].strip())
            self.coord = np.array([self.x, self.y, self.z])
            self.element = line[76:78].strip() if len(line) > 78 else 'C'
        except Exception as e:
            raise ValueError(f"Could not parse PDB line: {line.rstrip()}\nError: {e}")


def parse_pdb_manual(filename):
    """
    Manually parse PDB file to handle both well-formed and malformed PDBs.
    For well-formed PDBs: use residue numbers from file.
    For malformed PDBs (DiffDock-PP): each atom is treated as its own residue.
    Returns dict of {chain_id: {res_num: {atom_name: PDBAtom}}}
    """
    atoms_by_chain = defaultdict(lambda: defaultdict(dict))
    
    with open(filename) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            
            try:
                atom = PDBAtom(line)
            except ValueError:
                continue
            
            chain_id = atom.chain_id
            res_num = atom.res_num
            
            # If residue number is not valid/present, use atom number as residue identifier
            # This handles DiffDock-PP CA-only PDB where each atom is a residue
            if res_num is None or res_num == 0:
                res_num = atom.atom_num
            
            atoms_by_chain[chain_id][res_num][atom.atom_name] = atom
    
    return atoms_by_chain


def combine_structures(receptor_pdb, ligand_pdb, output_path):
    """Combine receptor and ligand PDB files into a single complex."""
    with open(output_path, 'w') as out:
        with open(receptor_pdb) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    out.write(line)
        
        out.write("TER\n")
        
        with open(ligand_pdb) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    out.write(line)
        
        out.write("END\n")
    
    return {}


def classify_contact(r_resname, l_resname, distance):
    """
    Classify interaction type based on amino acid properties and distance.
    Returns: (interaction_type, is_hbond)
    """
    is_hbond = False
    
    # Hydrogen bond check (distance < 3.5 Å + polar residues)
    if distance < 3.5:
        if r_resname in POLAR_RESIDUES or l_resname in POLAR_RESIDUES:
            is_hbond = True
    
    # Hydrophobic interaction
    if r_resname in HYDROPHOBIC_RESIDUES and l_resname in HYDROPHOBIC_RESIDUES:
        return 'hydrophobic', is_hbond
    
    # Electrostatic interaction
    r_positive = r_resname in POSITIVE_RESIDUES
    r_negative = r_resname in NEGATIVE_RESIDUES
    l_positive = l_resname in POSITIVE_RESIDUES
    l_negative = l_resname in NEGATIVE_RESIDUES
    
    if (r_positive and l_negative) or (r_negative and l_positive):
        return 'electrostatic', is_hbond
    
    # Hydrogen bond (already checked via polar atoms)
    if is_hbond:
        return 'hydrogen_bond', is_hbond
    
    # Polar or other interactions
    if r_resname in POLAR_RESIDUES or l_resname in POLAR_RESIDUES:
        return 'polar', is_hbond
    
    return 'other', False


def calculate_interface(receptor_pdb, ligand_pdb, distance_cutoff=5.0, hbond_cutoff=3.5):
    """
    Calculate interface contacts manually without strict PDB parsing.
    Handles structures with single atom per residue (e.g., DiffDock-PP CA-only output).
    """
    contacts = []
    interface_residues_receptor = defaultdict(set)
    interface_residues_ligand = defaultdict(set)
    
    try:
        rec_atoms = parse_pdb_manual(receptor_pdb)
        lig_atoms = parse_pdb_manual(ligand_pdb)
    except Exception as e:
        raise RuntimeError(f"Failed to parse PDB files: {e}")
    
    # Iterate through all residue pairs
    for rec_chain, rec_residues in rec_atoms.items():
        for rec_resnum, rec_atoms_dict in rec_residues.items():
            # Get residue name from any atom in residue
            if not rec_atoms_dict:
                continue
            r_resname = next(iter(rec_atoms_dict.values())).res_name
            
            for lig_chain, lig_residues in lig_atoms.items():
                for lig_resnum, lig_atoms_dict in lig_residues.items():
                    if not lig_atoms_dict:
                        continue
                    
                    l_resname = next(iter(lig_atoms_dict.values())).res_name
                    
                    # Calculate atom-atom distances
                    # If either structure only has CA atoms, use all available atoms
                    for r_atom_name, r_atom in rec_atoms_dict.items():
                        for l_atom_name, l_atom in lig_atoms_dict.items():
                            try:
                                distance_vec = r_atom.coord - l_atom.coord
                                distance = np.sqrt(np.sum(distance_vec ** 2))
                            except:
                                continue
                            
                            if distance < distance_cutoff:
                                # Classify the contact
                                interaction_type, is_hbond = classify_contact(
                                    r_resname, l_resname, distance
                                )
                                
                                contact = {
                                    'receptor_chain': rec_chain,
                                    'receptor_resnum': rec_resnum,
                                    'receptor_resname': r_resname,
                                    'receptor_atom': r_atom_name,
                                    'ligand_chain': lig_chain,
                                    'ligand_resnum': lig_resnum,
                                    'ligand_resname': l_resname,
                                    'ligand_atom': l_atom_name,
                                    'distance': float(distance),
                                    'interaction_type': interaction_type,
                                    'is_hydrogen_bond': is_hbond
                                }
                                
                                contacts.append(contact)
                                
                                # Track interface residues
                                interface_residues_receptor[rec_chain].add(
                                    (rec_resnum, r_resname)
                                )
                                interface_residues_ligand[lig_chain].add(
                                    (lig_resnum, l_resname)
                                )
    
    return contacts, interface_residues_receptor, interface_residues_ligand


def deduplicate_contacts(contacts):
    """
    Deduplicate contacts by residue-residue pairs, keeping the closest distance.
    """
    contact_map = {}
    
    for contact in contacts:
        key = (
            contact['receptor_chain'],
            contact['receptor_resnum'],
            contact['ligand_chain'],
            contact['ligand_resnum']
        )
        
        if key not in contact_map:
            contact_map[key] = contact
        else:
            if contact['distance'] < contact_map[key]['distance']:
                contact_map[key] = contact
    
    return list(contact_map.values())


def generate_interface_analysis_txt(contacts, interface_res_rec, interface_res_lig, output_path):
    """Generate human-readable interface analysis text file."""
    
    # Count interactions by type
    interaction_counts = defaultdict(int)
    hbond_count = 0
    
    for contact in contacts:
        interaction_counts[contact['interaction_type']] += 1
        if contact['is_hydrogen_bond']:
            hbond_count += 1
    
    # Count unique residues
    unique_rec_residues = sum(len(residues) for residues in interface_res_rec.values())
    unique_lig_residues = sum(len(residues) for residues in interface_res_lig.values())
    
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("ANTIBODY-ANTIGEN INTERFACE ANALYSIS\n")
        f.write("=" * 80 + "\n\n")
        
        # Summary statistics
        f.write("SUMMARY STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total atomic contacts (distance < 5.0 Å): {len(contacts)}\n")
        f.write(f"Unique receptor interface residues: {unique_rec_residues}\n")
        f.write(f"Unique ligand interface residues: {unique_lig_residues}\n")
        f.write(f"Hydrogen bonds (distance < 3.5 Å): {hbond_count}\n\n")
        
        # Interaction type breakdown
        f.write("INTERACTION TYPES BREAKDOWN\n")
        f.write("-" * 80 + "\n")
        total_interactions = len(contacts)
        for itype in sorted(interaction_counts.keys()):
            count = interaction_counts[itype]
            percentage = (count / total_interactions * 100) if total_interactions > 0 else 0
            f.write(f"{itype.upper():20s}: {count:4d} ({percentage:5.1f}%)\n")
        f.write("\n")
        
        # Receptor interface residues
        f.write("RECEPTOR (ANTIBODY) INTERFACE RESIDUES\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Chain':<6} {'ResNum':<8} {'ResName':<8}\n")
        for chain_id in sorted(interface_res_rec.keys()):
            residues = sorted(interface_res_rec[chain_id])
            for resnum, resname in residues:
                f.write(f"{chain_id:<6} {resnum:<8} {resname:<8}\n")
        f.write("\n")
        
        # Ligand interface residues
        f.write("LIGAND (ANTIGEN) INTERFACE RESIDUES\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Chain':<6} {'ResNum':<8} {'ResName':<8}\n")
        for chain_id in sorted(interface_res_lig.keys()):
            residues = sorted(interface_res_lig[chain_id])
            for resnum, resname in residues:
                f.write(f"{chain_id:<6} {resnum:<8} {resname:<8}\n")
        f.write("\n")
        
        # Top contacts
        f.write("TOP 20 CLOSEST CONTACTS\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Dist(Å)':<8} {'Receptor':<20} {'Ligand':<20} {'Type':<15}\n")
        f.write("-" * 80 + "\n")
        
        sorted_contacts = sorted(contacts, key=lambda x: x['distance'])[:20]
        for contact in sorted_contacts:
            rec_str = f"{contact['receptor_chain']}{contact['receptor_resnum']}{contact['receptor_resname']}"
            lig_str = f"{contact['ligand_chain']}{contact['ligand_resnum']}{contact['ligand_resname']}"
            f.write(
                f"{contact['distance']:<8.2f} {rec_str:<20} {lig_str:<20} "
                f"{contact['interaction_type']:<15}\n"
            )


def generate_contact_residues_json(interface_res_rec, interface_res_lig, output_path):
    """Generate machine-readable contact residues JSON file."""
    
    receptor_residues = []
    for chain_id in sorted(interface_res_rec.keys()):
        for resnum, resname in sorted(interface_res_rec[chain_id]):
            receptor_residues.append({
                'chain': chain_id,
                'residue_number': resnum,
                'residue_name': resname
            })
    
    ligand_residues = []
    for chain_id in sorted(interface_res_lig.keys()):
        for resnum, resname in sorted(interface_res_lig[chain_id]):
            ligand_residues.append({
                'chain': chain_id,
                'residue_number': resnum,
                'residue_name': resname
            })
    
    output = {
        'receptor_interface_residues': receptor_residues,
        'ligand_interface_residues': ligand_residues,
        'total_receptor_residues': len(receptor_residues),
        'total_ligand_residues': len(ligand_residues)
    }
    
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description='Analyze antibody-antigen binding interface'
    )
    parser.add_argument(
        '--receptor_pdb',
        required=True,
        help='Path to processed antibody (receptor) PDB file'
    )
    parser.add_argument(
        '--best_pose_pdb',
        required=True,
        help='Path to best docking pose (rank1.pdb) from Node 4'
    )
    parser.add_argument(
        '--confidence_json',
        default=None,
        help='Path to confidence scores JSON for reference'
    )
    parser.add_argument(
        '--output_dir',
        default='5_diffdock_analysis/outputs',
        help='Output directory for results'
    )
    parser.add_argument(
        '--distance_cutoff',
        type=float,
        default=5.0,
        help='Distance cutoff for contacts in Angstroms'
    )
    parser.add_argument(
        '--hbond_cutoff',
        type=float,
        default=3.5,
        help='Distance cutoff for hydrogen bonds in Angstroms'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Validate input files
    if not os.path.exists(args.receptor_pdb):
        raise FileNotFoundError(f"Receptor PDB not found: {args.receptor_pdb}")
    if not os.path.exists(args.best_pose_pdb):
        raise FileNotFoundError(f"Best pose PDB not found: {args.best_pose_pdb}")
    
    try:
        # 1. Combine structures
        final_complex_path = os.path.join(args.output_dir, 'final_complex.pdb')
        print(f"Combining antibody and best pose into: {final_complex_path}")
        combine_structures(args.receptor_pdb, args.best_pose_pdb, final_complex_path)
        
        # 2. Calculate interface
        print(f"Calculating interface contacts (cutoff: {args.distance_cutoff} Å)...")
        contacts, interface_res_rec, interface_res_lig = calculate_interface(
            args.receptor_pdb,
            args.best_pose_pdb,
            distance_cutoff=args.distance_cutoff,
            hbond_cutoff=args.hbond_cutoff
        )
        
        # 3. Deduplicate by residue pairs
        contacts = deduplicate_contacts(contacts)
        
        print(f"Found {len(contacts)} unique residue-pair contacts")
        
        # 4. Generate text report
        analysis_txt = os.path.join(args.output_dir, 'interface_analysis.txt')
        generate_interface_analysis_txt(contacts, interface_res_rec, interface_res_lig, analysis_txt)
        print(f"Generated: {analysis_txt}")
        
        # 5. Generate contact residues JSON
        contact_residues_json = os.path.join(args.output_dir, 'contact_residues.json')
        generate_contact_residues_json(interface_res_rec, interface_res_lig, contact_residues_json)
        print(f"Generated: {contact_residues_json}")
        
        # 6. Count interaction types
        interaction_counts = defaultdict(int)
        hbond_count = 0
        for contact in contacts:
            interaction_counts[contact['interaction_type']] += 1
            if contact['is_hydrogen_bond']:
                hbond_count += 1
        
        # 7. Generate data.json
        data = {
            'status': 'success',
            'node': 'interface_analysis',
            'total_contacts': len(contacts),
            'unique_receptor_residues': sum(len(r) for r in interface_res_rec.values()),
            'unique_ligand_residues': sum(len(r) for r in interface_res_lig.values()),
            'hydrogen_bonds': hbond_count,
            'interaction_types': dict(interaction_counts),
            'receptor_pdb': os.path.basename(args.receptor_pdb),
            'best_pose_pdb': os.path.basename(args.best_pose_pdb),
            'final_complex_pdb': os.path.basename(final_complex_path),
            'output_files': {
                'interface_analysis_txt': os.path.basename(analysis_txt),
                'contact_residues_json': os.path.basename(contact_residues_json),
                'final_complex_pdb': os.path.basename(final_complex_path)
            }
        }
        
        # Add confidence JSON reference if provided
        if args.confidence_json and os.path.exists(args.confidence_json):
            with open(args.confidence_json) as f:
                conf_data = json.load(f)
                data['confidence_reference'] = conf_data.get('scores', [])[0:1]
        
        data_json_path = os.path.join(args.output_dir, 'data.json')
        with open(data_json_path, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Generated: {data_json_path}")
        
        # Print summary
        print("\n" + "=" * 80)
        print("INTERFACE ANALYSIS COMPLETE")
        print("=" * 80)
        print(json.dumps(data, indent=2))
        
        return 0
        
    except Exception as e:
        error_data = {
            'status': 'error',
            'node': 'interface_analysis',
            'error': str(e)
        }
        data_json_path = os.path.join(args.output_dir, 'data.json')
        with open(data_json_path, 'w') as f:
            json.dump(error_data, f, indent=2)
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
