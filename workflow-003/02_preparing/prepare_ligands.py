import os
import sys
import subprocess

def convert_ligands(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    ligands = [f for f in os.listdir(input_dir) if f.endswith('.pdb')]

    for lig in ligands:
        inp = os.path.join(input_dir, lig)
        out = os.path.join(output_dir, lig.replace('.pdb', '.pdbqt'))
        print(f"Converting {lig} → {os.path.basename(out)}")

        cmd = [
            "python2", "/opt/mgltools/bin/prepare_ligand4.py",
            "-l", inp,
            "-o", out
        ]
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    convert_ligands(sys.argv[1], sys.argv[2])