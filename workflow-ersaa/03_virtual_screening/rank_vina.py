import os
import pandas as pd

def extract_affinities(log_dir):
    """Extract and rank binding affinities from AutoDock Vina log files."""
    data = []
    for file in os.listdir(log_dir):
        if file.endswith(".log"):
            with open(os.path.join(log_dir, file)) as f:
                for line in f:
                    if line.strip().startswith("1 "):  # first mode
                        parts = line.split()
                        if len(parts) > 1:
                            data.append({
                                "Ligand": file.replace(".log", ""),
                                "Affinity (kcal/mol)": float(parts[1])
                            })
                        break
    if data:
        df = pd.DataFrame(data).sort_values("Affinity (kcal/mol)")
        df.to_excel(os.path.join(log_dir, "binding_affinities.xlsx"), index=False)
        print("Results saved to binding_affinities.xlsx")
        print(df.head())
    else:
        print("No affinity values found.")

if __name__ == "__main__":
    log_dir = "./results"
    extract_affinities(log_dir)