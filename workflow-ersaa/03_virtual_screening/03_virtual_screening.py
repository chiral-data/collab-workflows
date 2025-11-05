import os
import subprocess
import pandas as pd

def run_p2rank(protein_pdb, output_dir):
    """Run P2Rank to predict binding pockets."""
    os.makedirs(output_dir, exist_ok=True)
    subprocess.run(["prank", "predict", protein_pdb, "--output_dir", output_dir], check=True)
    print("P2Rank prediction completed.")

def choose_pocket(p2rank_dir):
    """Let user choose a pocket from P2Rank results."""
    csv_path = os.path.join(p2rank_dir, "predictions.csv")
    df = pd.read_csv(csv_path)
    print("\nDetected pockets:\n", df[["rank", "center_x", "center_y", "center_z", "score"]])
    pocket = input("Enter pocket rank to use (default = 1): ") or "1"
    row = df[df["rank"] == int(pocket)].iloc[0]
    return float(row["center_x"]), float(row["center_y"]), float(row["center_z"])

def run_vina(protein_pdbqt, ligand_pdbqt, center, box_size, output_dir):
    """Run AutoDock Vina for one ligand."""
    os.makedirs(output_dir, exist_ok=True)
    x, y, z = center
    out_file = os.path.join(output_dir, ligand_pdbqt.replace(".pdbqt", "_out.pdbqt"))
    log_file = out_file.replace(".pdbqt", ".log")

    cmd = [
        "vina",
        "--receptor", protein_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(x), "--center_y", str(y), "--center_z", str(z),
        "--size_x", str(box_size), "--size_y", str(box_size), "--size_z", str(box_size),
        "--out", out_file,
        "--log", log_file
    ]
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    protein_pdb = "./input/protein_prepared.pdb"
    protein_pdbqt = "./input/protein_prepared.pdbqt"
    ligands_dir = "./input/ligands_pdbqt"
    output_dir = "./results"
    p2rank_output = "./p2rank_output"

    run_p2rank(protein_pdb, p2rank_output)
    center = choose_pocket(p2rank_output)

    box_size = input("Enter box size (default = 80): ") or "80"
    box_size = int(box_size)

    for ligand in os.listdir(ligands_dir):
        if ligand.endswith(".pdbqt"):
            ligand_path = os.path.join(ligands_dir, ligand)
            print(f"Docking {ligand}...")
            run_vina(protein_pdbqt, ligand_path, center, box_size, output_dir)

    print("\nDocking finished. To analyze results, run:")
    print("  python rank_vina.py")