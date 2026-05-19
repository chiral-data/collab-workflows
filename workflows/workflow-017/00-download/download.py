import os
import urllib.request

os.makedirs("outputs", exist_ok=True)

url = "https://files.rcsb.org/download/1BRS.pdb"
print("Downloading 1BRS.pdb from RCSB ...", flush=True)
urllib.request.urlretrieve(url, "1BRS.pdb")

with open("1BRS.pdb") as f:
    lines = f.readlines()


def write_chain(lines, chain_id, output):
    selected = [l for l in lines
                if (l.startswith(("ATOM", "HETATM")) and len(l) > 21 and l[21] == chain_id)
                or l.startswith("END")]
    with open(output, "w") as f:
        f.writelines(selected)
    atom_count = sum(1 for l in selected if l.startswith("ATOM"))
    print(f"Saved {output} ({atom_count} ATOM records)", flush=True)


write_chain(lines, "A", os.path.join("outputs", "protein2_barnase.pdb"))
write_chain(lines, "D", os.path.join("outputs", "protein1_barstar.pdb"))
os.remove("1BRS.pdb")
