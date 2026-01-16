#!/usr/bin/env python3
"""
DiffDock-PP Protein-Only Extractor

This script extracts poses from DiffDock-PP pickle files while:
1. Only including ATOM records (standard protein atoms)
2. Excluding HETATM records (cofactors, ligands, waters)
3. Maintaining proper bond lengths through rigid body transformation

Usage:
    python extract_poses_protein_only.py predictions.pkl output_dir/
"""

import os
import pickle
import argparse
import numpy as np
import torch
from pathlib import Path


class ProteinOnlyDiffDockExtractor:
    """Extract only protein atoms from DiffDock-PP poses"""
    
    def __init__(self, predictions_file, output_dir, run_name=None):
        self.predictions_file = predictions_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.run_name = run_name or self._detect_run_name()
    
    def _detect_run_name(self):
        """Extract run name from predictions file"""
        filename = Path(self.predictions_file).stem
        return filename.replace('_predictions', '')
    
    def _filter_protein_atoms(self, pdb_file):
        """Filter PDB file to only include protein ATOM records"""
        protein_atoms = []
        ca_coords = []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                # Only include ATOM records (exclude HETATM like GDP)
                if line.startswith('ATOM'):
                    protein_atoms.append(line)
                    
                    # Extract CA coordinates for transformation
                    if line[12:16].strip() == 'CA':
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ca_coords.append([x, y, z])
        
        return protein_atoms, np.array(ca_coords)
    
    def _write_protein_only_pdb(self, protein_atoms, output_file, new_ca_coords=None, chain='A'):
        """Write protein-only PDB with optional coordinate transformation"""
        with open(output_file, 'w') as f:
            ca_index = 0
            atom_num = 1
            
            for line in protein_atoms:
                if line.startswith('ATOM'):
                    # Update atom number and chain
                    updated_line = (line[:6] + f"{atom_num:5d}" + line[11:21] + 
                                  f"{chain}" + line[22:])
                    
                    # If we have new CA coordinates, apply transformation
                    if new_ca_coords is not None and line[12:16].strip() == 'CA':
                        if ca_index < len(new_ca_coords):
                            x, y, z = new_ca_coords[ca_index]
                            updated_line = (updated_line[:30] + 
                                          f"{x:8.3f}{y:8.3f}{z:8.3f}" + 
                                          updated_line[54:])
                            ca_index += 1
                    
                    f.write(updated_line)
                    atom_num += 1
            
            f.write("END\n")
    
    def _calculate_rigid_body_transform(self, original_coords, new_coords):
        """Calculate rigid body transformation using Kabsch algorithm"""
        # Center coordinates
        centroid_orig = np.mean(original_coords, axis=0)
        centroid_new = np.mean(new_coords, axis=0)
        
        orig_centered = original_coords - centroid_orig
        new_centered = new_coords - centroid_new
        
        # Calculate rotation matrix using SVD
        H = orig_centered.T @ new_centered
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        
        # Ensure proper rotation (det(R) = 1)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # Calculate translation
        t = centroid_new - R @ centroid_orig
        
        return {'rotation': R, 'translation': t}
    
    def _apply_transformation_to_all_atoms(self, protein_atoms, original_ca, new_ca, output_file, chain):
        """Apply rigid body transformation to all protein atoms"""
        if len(original_ca) != len(new_ca):
            print(f"Warning: CA coordinate mismatch ({len(original_ca)} vs {len(new_ca)})")
            min_len = min(len(original_ca), len(new_ca))
            original_ca = original_ca[:min_len]
            new_ca = new_ca[:min_len]
        
        # Calculate transformation
        transform = self._calculate_rigid_body_transform(original_ca, new_ca)
        R = transform['rotation']
        t = transform['translation']
        
        with open(output_file, 'w') as f:
            atom_num = 1
            
            for line in protein_atoms:
                if line.startswith('ATOM'):
                    # Extract original coordinates
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    orig_coord = np.array([x, y, z])
                    
                    # Apply transformation
                    new_coord = R @ orig_coord + t
                    nx, ny, nz = new_coord
                    
                    # Write transformed atom
                    updated_line = (line[:6] + f"{atom_num:5d}" + line[11:21] + 
                                  f"{chain}" + line[22:30] + 
                                  f"{nx:8.3f}{ny:8.3f}{nz:8.3f}" + line[54:])
                    f.write(updated_line)
                    atom_num += 1
            
            f.write("END\n")
        
        print(f"Applied rigid body transformation to all atoms: {output_file}")
    
    def extract_poses(self, max_poses=5):
        """Extract protein-only poses with proper transformations"""
        print(f"Loading predictions from {self.predictions_file}")
        
        try:
            with open(self.predictions_file, 'rb') as f:
                predictions_data = pickle.load(f)
        except Exception as e:
            print(f"ERROR: Could not load pickle file: {e}")
            return []
        
        if isinstance(predictions_data, list) and len(predictions_data) > 0:
            pose_list = predictions_data[0] if isinstance(predictions_data[0], list) else predictions_data
            
            # Filter valid poses
            valid_poses = [(hetero_data, conf) for hetero_data, conf in pose_list 
                          if isinstance((hetero_data, conf), tuple) and conf != float('inf')]
            
            num_poses = min(max_poses, len(valid_poses))
            print(f"Extracting {num_poses} protein-only poses from {len(valid_poses)} valid poses")
            
            # Check for original PDB files
            receptor_pdb = f"{self.run_name}_r_b.pdb"
            ligand_pdb = f"{self.run_name}_l_b.pdb"
            
            if not (os.path.exists(receptor_pdb) and os.path.exists(ligand_pdb)):
                print(f"Warning: Original PDB files not found ({receptor_pdb}, {ligand_pdb})")
                return []
            
            # Filter original structures to protein-only
            receptor_atoms, receptor_ca = self._filter_protein_atoms(receptor_pdb)
            ligand_atoms, ligand_ca = self._filter_protein_atoms(ligand_pdb)
            
            print(f"Filtered to protein-only: {len(receptor_atoms)} receptor atoms, {len(ligand_atoms)} ligand atoms")
            
            extracted_files = []
            
            for i in range(num_poses):
                hetero_data, confidence_score = valid_poses[i]
                pose_num = i + 1
                
                print(f"\nExtracting pose {pose_num} (confidence: {confidence_score:.3f})")
                
                # Get new CA coordinates from DiffDock-PP
                new_receptor_ca = self._tensor_to_coords(hetero_data['receptor'].pos)
                new_ligand_ca = self._tensor_to_coords(hetero_data['ligand'].pos)
                
                # Apply transformation to all atoms (protein-only)
                receptor_file = self.output_dir / f"pose_{pose_num}_receptor_protein.pdb"
                ligand_file = self.output_dir / f"pose_{pose_num}_ligand_protein.pdb"
                complex_file = self.output_dir / f"pose_{pose_num}_complex_protein.pdb"
                
                self._apply_transformation_to_all_atoms(
                    receptor_atoms, receptor_ca, new_receptor_ca, receptor_file, 'A'
                )
                extracted_files.append(receptor_file)
                
                self._apply_transformation_to_all_atoms(
                    ligand_atoms, ligand_ca, new_ligand_ca, ligand_file, 'B'
                )
                extracted_files.append(ligand_file)
                
                # Create complex
                self._combine_protein_structures(receptor_file, ligand_file, complex_file, confidence_score)
                extracted_files.append(complex_file)
            
            return extracted_files
        
        print("ERROR: Could not parse predictions data")
        return []
    
    def _tensor_to_coords(self, tensor):
        """Convert tensor to numpy coordinates"""
        if hasattr(tensor, 'cpu'):
            return tensor.cpu().numpy()
        return tensor
    
    def _combine_protein_structures(self, receptor_file, ligand_file, output_file, confidence):
        """Combine receptor and ligand into protein-only complex"""
        with open(output_file, 'w') as out_f:
            out_f.write(f"REMARK   Confidence Score: {confidence:.3f}\n")
            out_f.write("REMARK   Protein-only structure (HETATM records excluded)\n")
            
            # Write receptor
            with open(receptor_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        out_f.write(line)
            
            out_f.write("TER\n")
            
            # Write ligand with updated atom numbers
            receptor_atom_count = self._count_atoms_in_file(receptor_file)
            
            with open(ligand_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        atom_num = int(line[6:11]) + receptor_atom_count
                        new_line = line[:6] + f"{atom_num:5d}" + line[11:]
                        out_f.write(new_line)
            
            out_f.write("END\n")
        
        print(f"Created protein-only complex: {output_file}")
    
    def _count_atoms_in_file(self, pdb_file):
        """Count ATOM records in PDB file"""
        count = 0
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    count += 1
        return count


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description="Extract protein-only poses from DiffDock-PP pickle files"
    )
    parser.add_argument('predictions_file', help='DiffDock-PP predictions pickle file')
    parser.add_argument('output_dir', help='Output directory for extracted poses')
    parser.add_argument('--max_poses', type=int, default=5, help='Maximum poses to extract')
    parser.add_argument('--run_name', help='Run name (auto-detected if not provided)')
    
    args = parser.parse_args()
    
    extractor = ProteinOnlyDiffDockExtractor(
        args.predictions_file, args.output_dir, args.run_name
    )
    
    extracted_files = extractor.extract_poses(args.max_poses)
    
    print(f"\n=== Protein-Only Extraction Complete ===")
    print(f"Extracted {len(extracted_files)} protein-only files:")
    for file in extracted_files:
        print(f"  {file}")
    print("==========================================")
    
    return len(extracted_files) > 0


if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)