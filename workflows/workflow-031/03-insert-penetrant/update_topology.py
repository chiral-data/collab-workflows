"""
Append penetrant ITP include and molecule count to topol.top → topol_penetrant.top.
"""
import sys, pathlib, re

PENETRANT   = sys.argv[1]   # "O2" or "H2O"
N_PENETRANT = int(sys.argv[2])
MOL_NAME    = "O2" if PENETRANT == "O2" else "SOL"

top = pathlib.Path("topol.top").read_text()

# Insert #include "penetrant.itp" before [ system ]
itp_line = '#include "penetrant.itp"\n'
if itp_line not in top:
    top = re.sub(r'(\[ system \])', itp_line + r'\1', top)

# Append molecule entry under [ molecules ]
mol_entry = f"{MOL_NAME}             {N_PENETRANT}\n"
if mol_entry not in top:
    top = top.rstrip() + "\n" + mol_entry

pathlib.Path("topol_penetrant.top").write_text(top)
print(f"  Updated topology: +{N_PENETRANT} {MOL_NAME}")
