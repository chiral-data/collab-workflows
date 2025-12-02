#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unpack File - Extract ligand ZIP file
Input: ligands.zip
Output: ligands.sdf (SDF files in extracted directory)
"""

import os
import re
import zipfile
import glob
from rdkit import Chem

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
INPUT_FILE = os.path.join(INPUT_DIR, "ligands.zip")
# Fallback to outputs for backward compatibility
if not os.path.exists(INPUT_FILE):
    INPUT_FILE = os.path.join(OUTPUT_DIR, "ligands.zip")

def main():
    """Main execution function"""
    print("=== Extract ligand ZIP file ===")
    
    if not os.path.exists(INPUT_FILE):
        print(f"❌ Error: {INPUT_FILE} not found.")
        exit(1)
    
    print(f"Source: {INPUT_FILE}")
    print(f"Destination: {OUTPUT_DIR}")
    
    try:
        with zipfile.ZipFile(INPUT_FILE, "r") as z:
            z.extractall(OUTPUT_DIR)
            # Show list of extracted files
            file_list = z.namelist()
            sdf_files = [f for f in file_list if f.endswith(".sdf")]
            print(f"✅ Extraction complete: {len(file_list)} files")
            if sdf_files:
                print(f"   Number of SDF files: {len(sdf_files)}")
                print(f"   First 5 files: {sdf_files[:5]}")
        
        # Validate and remove invalid SDF files
        print(f"\nValidating SDF files...")
        validate_and_clean_sdf_files(OUTPUT_DIR)
        
        # Rename SDF files based on ligand names
        print(f"\nRenaming SDF files based on ligand names...")
        rename_sdf_files(OUTPUT_DIR)
        
    except Exception as e:
        print(f"❌ Error occurred during extraction: {e}")
        exit(1)


def validate_sdf_file(sdf_file_path):
    """
    Validate SDF file by checking if it contains valid molecules.
    
    Args:
        sdf_file_path: Path to SDF file
        
    Returns:
        True if file is valid (contains at least one valid molecule), False otherwise
    """
    try:
        # Check if file is empty
        if os.path.getsize(sdf_file_path) == 0:
            return False
        
        # Try to read molecules from SDF file using RDKit
        supplier = Chem.SDMolSupplier(sdf_file_path)
        
        # Check if at least one valid molecule exists
        valid_molecule_found = False
        for mol in supplier:
            if mol is not None:
                valid_molecule_found = True
                break
        
        return valid_molecule_found
        
    except Exception as e:
        # If any error occurs during reading, consider the file invalid
        return False


def validate_and_clean_sdf_files(output_dir):
    """
    Validate all SDF files and remove invalid ones.
    
    Args:
        output_dir: Directory containing extracted SDF files
    """
    # Find all SDF files recursively
    sdf_pattern = os.path.join(output_dir, "**", "*.sdf")
    sdf_files = glob.glob(sdf_pattern, recursive=True)
    
    if not sdf_files:
        print("   No SDF files found to validate.")
        return
    
    valid_count = 0
    invalid_count = 0
    
    for sdf_file in sdf_files:
        if validate_sdf_file(sdf_file):
            valid_count += 1
        else:
            # Remove invalid file
            try:
                os.remove(sdf_file)
                invalid_count += 1
                print(f"   ❌ Removed invalid SDF file: {os.path.basename(sdf_file)}")
            except Exception as e:
                print(f"   ⚠ Error removing {os.path.basename(sdf_file)}: {e}")
    
    print(f"✅ Validation complete: {valid_count} valid file(s), {invalid_count} invalid file(s) removed")


def sanitize_filename(name):
    """
    Sanitize ligand name to be used as a filename.
    Removes or replaces characters that are not allowed in filenames.
    Preserves stereochemistry notation like (+/-), (±), etc.
    
    Args:
        name: Original ligand name
        
    Returns:
        Sanitized filename-safe string
    """
    # Preserve stereochemistry notation at the start (e.g., (+/-), (±), (+)-, (-)-)
    # Extract and temporarily store it
    stereochemistry = ""
    stereochemistry_pattern = r'^\([+\-/±]+\)\s*-?\s*'
    match = re.match(stereochemistry_pattern, name)
    if match:
        stereochemistry = match.group(0)
        # Remove the stereochemistry part from the name for processing
        name = name[len(stereochemistry):]
    
    # Remove or replace invalid characters for filenames
    # Replace spaces, slashes, and other special characters with underscores
    sanitized = re.sub(r'[<>:"/\\|?*]', '_', name)
    sanitized = re.sub(r'\s+', '_', sanitized)  # Replace spaces with underscores
    sanitized = re.sub(r'[()]', '', sanitized)  # Remove parentheses (except preserved stereochemistry)
    sanitized = re.sub(r'_+', '_', sanitized)  # Replace multiple underscores with single
    sanitized = sanitized.strip('_')  # Remove leading/trailing underscores
    
    # Remove leading + and - characters (but only if not part of preserved stereochemistry)
    # This handles cases like "+-", "-+", "+", "-" at the start
    sanitized = re.sub(r'^[+\-]+', '', sanitized)
    
    # Remove leading/trailing underscores again after removing +/-
    sanitized = sanitized.strip('_')
    
    # Sanitize the stereochemistry part if it exists
    if stereochemistry:
        # Preserve the stereochemistry notation format: "(-)-" -> "-_-", "(+)-" -> "+_-", etc.
        # Replace parentheses with underscores to make it filename-safe
        paren_content = re.search(r'\(([+\-/±]+)\)', stereochemistry)
        if paren_content:
            # Get the content inside parentheses
            content = paren_content.group(1)
            # Replace slashes with underscores for multi-character cases
            if '/' in content:
                content = content.replace('/', '_')
            # Format: "(-)-" -> "-_-", "(+)-" -> "+_-", "(+/-)-" -> "+_-_-"
            # Opening parenthesis becomes "-", closing parenthesis becomes "_", keep trailing hyphen
            stereochemistry_clean = f"-{content}_-"
        else:
            # Fallback: replace parentheses with underscores
            stereochemistry_clean = stereochemistry.replace('(', '-').replace(')', '_')
            stereochemistry_clean = re.sub(r'\s+', '_', stereochemistry_clean)
        
        stereochemistry_clean = re.sub(r'_+', '_', stereochemistry_clean)  # Multiple underscores
        stereochemistry_clean = stereochemistry_clean.strip('_')
        # Combine stereochemistry with the main name
        # Directly combine without additional separator (the trailing "-" in stereochemistry serves as separator)
        sanitized = f"{stereochemistry_clean}{sanitized}" if sanitized else stereochemistry_clean
    
    # Limit filename length (keep it reasonable)
    if len(sanitized) > 200:
        sanitized = sanitized[:200]
    
    # If empty after sanitization, use a default name
    if not sanitized:
        sanitized = "unnamed_ligand"
    
    return sanitized


def extract_ligand_name(sdf_file_path):
    """
    Extract ligand name from SDF file.
    Looks for '> <Name>' tag first, then '> <Synonyms>' tag if name is not found.
    Returns the first synonym if name is not available.
    
    Args:
        sdf_file_path: Path to SDF file
        
    Returns:
        Ligand name string, or None if not found
    """
    try:
        with open(sdf_file_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
            
            # Priority 1: Look for '> <Name>' tag
            for i, line in enumerate(lines):
                if line.strip() == '> <Name>':
                    # Get the next non-empty line
                    if i + 1 < len(lines):
                        name = lines[i + 1].strip()
                        if name:
                            return name
                    
                    # Sometimes there's an empty line after the tag
                    if i + 2 < len(lines):
                        name = lines[i + 2].strip()
                        if name:
                            return name
            
            # Priority 2: Look for '> <Synonyms>' tag if name not found
            for i, line in enumerate(lines):
                if line.strip() == '> <Synonyms>':
                    # Get the next non-empty line
                    if i + 1 < len(lines):
                        synonyms_line = lines[i + 1].strip()
                        if synonyms_line:
                            # Split by semicolon and get the first synonym
                            synonyms = [s.strip() for s in synonyms_line.split(';')]
                            if synonyms and synonyms[0]:
                                return synonyms[0]
                    
                    # Sometimes there's an empty line after the tag
                    if i + 2 < len(lines):
                        synonyms_line = lines[i + 2].strip()
                        if synonyms_line:
                            # Split by semicolon and get the first synonym
                            synonyms = [s.strip() for s in synonyms_line.split(';')]
                            if synonyms and synonyms[0]:
                                return synonyms[0]
            
            # Fallback: use first line if it looks like a name
            if lines:
                first_line = lines[0].strip()
                if first_line and not first_line.startswith(' ') and len(first_line) < 100:
                    return first_line
                    
    except Exception as e:
        print(f"   Warning: Could not read {sdf_file_path}: {e}")
    
    return None


def rename_sdf_files(output_dir):
    """
    Rename all SDF files in the output directory based on ligand names.
    
    Args:
        output_dir: Directory containing extracted SDF files
    """
    # Find all SDF files recursively
    sdf_pattern = os.path.join(output_dir, "**", "*.sdf")
    sdf_files = glob.glob(sdf_pattern, recursive=True)
    
    if not sdf_files:
        print("   No SDF files found to rename.")
        return
    
    renamed_count = 0
    failed_count = 0
    
    for sdf_file in sdf_files:
        try:
            # Extract ligand name
            ligand_name = extract_ligand_name(sdf_file)
            
            if not ligand_name:
                print(f"   ⚠ Could not extract name from: {os.path.basename(sdf_file)}")
                failed_count += 1
                continue
            
            # Sanitize the name for use as filename
            sanitized_name = sanitize_filename(ligand_name)
            
            # Get directory of the SDF file
            file_dir = os.path.dirname(sdf_file)
            new_file_path = os.path.join(file_dir, f"{sanitized_name}.sdf")
            
            # If the new filename is the same as the old one, skip
            if os.path.normpath(sdf_file) == os.path.normpath(new_file_path):
                continue
            
            # Handle duplicate names by adding a number suffix
            if os.path.exists(new_file_path):
                base_name = sanitized_name
                counter = 1
                while os.path.exists(new_file_path):
                    new_file_path = os.path.join(file_dir, f"{base_name}_{counter}.sdf")
                    counter += 1
            
            # Rename the file
            os.rename(sdf_file, new_file_path)
            renamed_count += 1
            
        except Exception as e:
            print(f"   ⚠ Error renaming {os.path.basename(sdf_file)}: {e}")
            failed_count += 1
    
    print(f"✅ Renaming complete: {renamed_count} files renamed")
    if failed_count > 0:
        print(f"   ⚠ {failed_count} files could not be renamed")


if __name__ == "__main__":
    main()

