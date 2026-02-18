#!/usr/bin/env python3
"""
Simple script to convert FASTA files to Boltz-2 format
"""

def convert_to_boltz_format(input_file, output_file, seq_id):
    """Convert regular FASTA to Boltz-2 format"""
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Get sequence without original header
    sequence = ''.join(line.strip() for line in lines[1:] if not line.startswith('>'))
    
    # Create Boltz-2 format header
    boltz_header = f">{seq_id}|protein"
    
    # Write Boltz-2 format file
    with open(output_file, 'w') as f:
        f.write(f"{boltz_header}\n{sequence}\n")
    
    print(f"Converted: {input_file} -> {output_file}")
    print(f"Sequence length: {len(sequence)} amino acids")

def extract_spike_rbd():
    """Extract RBD region from spike protein and create Boltz-2 file"""
    input_file = "1_mRNA/sequences/P0DTC2.fasta"
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Get full spike sequence
    full_sequence = ''.join(line.strip() for line in lines[1:])
    
    # Extract RBD (amino acids 331-524)
    rbd_sequence = full_sequence[330:524]
    
    # Write RBD in Boltz-2 format
    output_file = "1_mRNA/inputs/spike_rbd.fasta"
    with open(output_file, 'w') as f:
        f.write(f">spike_rbd|protein\n{rbd_sequence}\n")
    
    print(f"Extracted RBD: {len(rbd_sequence)} amino acids -> {output_file}")

def main():
    print("Converting sequences to Boltz-2 format...")
    
    # Create inputs directories
    import os
    os.makedirs("1_mRNA/inputs", exist_ok=True)
    os.makedirs("2_antibody/inputs", exist_ok=True)
    
    # Extract Spike RBD
    print("\n=== Spike RBD (Nobel 2023) ===")
    extract_spike_rbd()
    
    # Convert FMC63 protein
    print("\n=== FMC63 Antibody (Nobel 2018 connection) ===")
    convert_to_boltz_format(
        "2_antibody/sequences/FMC63-28Z.fasta",
        "2_antibody/inputs/fmc63.fasta", 
        "fmc63_antibody"
    )
    
    print("\n✅ All files converted to Boltz-2 format!")

if __name__ == "__main__":
    main()