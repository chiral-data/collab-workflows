#!/usr/bin/env python3
"""
Extract individual compounds from DrugCentral SDF library
Each compound becomes a separate SDF file for virtual screening
"""

import os
from pathlib import Path

def extract_compounds():
    """Extract individual compounds from multi-compound SDF file"""
    
    sdf_file = "drugcentral_structures.sdf"
    output_dir = "fda_drugs_sdf"
    
    if not os.path.exists(sdf_file):
        print(f"Error: SDF file {sdf_file} not found")
        return False
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    compound_count = 0
    current_compound = []
    compound_name = None
    compound_id = None
    
    with open(sdf_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            
            # Start of new compound (first line is compound name)
            if not current_compound and line and not line.startswith('  '):
                compound_name = line.strip()
                current_compound = [line]
            
            # Continue reading compound data
            elif current_compound:
                current_compound.append(line)
                
                # Extract compound ID from data field
                if line.startswith('>  <ID>'):
                    # Next line should contain the ID
                    continue
                elif len(current_compound) > 1 and current_compound[-2].startswith('>  <ID>') and line.strip().isdigit():
                    compound_id = line.strip()
                
                # End of compound record ($$$$)
                elif line == '$$$$':
                    if compound_name and compound_id:
                        # Clean compound name for filename
                        safe_name = "".join(c for c in compound_name if c.isalnum() or c in ('-', '_')).strip()
                        if not safe_name:
                            safe_name = f"compound_{compound_id}"
                        
                        # Write individual SDF file
                        output_file = f"{output_dir}/{safe_name}_{compound_id}.sdf"
                        with open(output_file, 'w') as out_f:
                            out_f.write('\n'.join(current_compound) + '\n')
                        
                        compound_count += 1
                        if compound_count % 100 == 0:
                            print(f"Extracted {compound_count} compounds...")
                    
                    # Reset for next compound
                    current_compound = []
                    compound_name = None
                    compound_id = None
    
    print(f"✅ Extracted {compound_count} individual SDF files to {output_dir}/")
    return True

def main():
    print("=== Extracting FDA Drug Library Compounds ===")
    print("Converting multi-compound SDF to individual files\n")
    
    if extract_compounds():
        print("\n✅ Extraction complete!")
        print("Individual SDF files ready for PDBQT conversion")
    else:
        print("\n❌ Extraction failed!")

if __name__ == "__main__":
    main()