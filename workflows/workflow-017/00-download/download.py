import os
import tempfile
import urllib.request

url = "https://files.rcsb.org/download/1BRS.pdb"
print("Downloading 1BRS.pdb from RCSB ...", flush=True)
tmp_fd, tmp_path = tempfile.mkstemp(suffix=".pdb")
os.close(tmp_fd)
urllib.request.urlretrieve(url, tmp_path)

with open(tmp_path) as f:
    lines = f.readlines()
os.remove(tmp_path)


def write_chain(lines, chain_id, output):
    selected = [l for l in lines
                if (l.startswith(("ATOM", "HETATM")) and len(l) > 21 and l[21] == chain_id)
                or l.startswith("END")]
    with open(output, "w") as f:
        f.writelines(selected)
    atom_count = sum(1 for l in selected if l.startswith("ATOM"))
    print(f"Saved {output} ({atom_count} ATOM records)", flush=True)


write_chain(lines, "A", "protein2_barnase.pdb")
write_chain(lines, "D", "protein1_barstar.pdb")
