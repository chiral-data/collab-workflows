import glob
import os
import sys
import yaml

VALID_CHARS = {
    "protein": set("ACDEFGHIKLMNPQRSTVWYBXZJUO"),
    "rna":     set("ACGU"),
    "dna":     set("ACGT"),
}
SEQUENCE_TYPES = set(VALID_CHARS)
SUPPORTED_TYPES = SEQUENCE_TYPES | {"ligand"}

yaml_files = glob.glob("./inputs/*.yaml") + glob.glob("./inputs/*.yml")
if not yaml_files:
    print("Error: no YAML file found in ./inputs/", file=sys.stderr)
    sys.exit(1)

data = yaml.safe_load(open(yaml_files[0]))
sequences = data.get("sequences")
if not sequences:
    print("Error: 'sequences' key missing or empty", file=sys.stderr)
    sys.exit(1)

for i, s in enumerate(sequences):
    for field in ("id", "type"):
        if field not in s:
            print(f"Error: chain {i} missing required field '{field}'", file=sys.stderr)
            sys.exit(1)

    seq_type = s["type"].lower()
    if seq_type not in SUPPORTED_TYPES:
        print(f"Error: chain '{s['id']}' has unsupported type '{s['type']}'. Must be one of: {sorted(SUPPORTED_TYPES)}", file=sys.stderr)
        sys.exit(1)

    if seq_type == "ligand":
        if not s.get("ccd") and not s.get("smiles"):
            print(f"Error: ligand chain '{s['id']}' must have 'ccd' or 'smiles'", file=sys.stderr)
            sys.exit(1)
    else:
        if "sequence" not in s:
            print(f"Error: chain '{s['id']}' ({seq_type}) missing required field 'sequence'", file=sys.stderr)
            sys.exit(1)
        bad = set(s["sequence"].upper()) - VALID_CHARS[seq_type]
        if bad:
            print(f"Error: chain '{s['id']}' ({seq_type}) contains invalid characters: {bad}", file=sys.stderr)
            sys.exit(1)
        s["sequence"] = s["sequence"].upper()

    s["type"] = seq_type

print(f"Validated {len(sequences)} chain(s):", flush=True)
for s in sequences:
    print(f"  chain {s['id']} ({s['type']}): {len(s['sequence'])} residues", flush=True)

os.makedirs("./outputs", exist_ok=True)
with open("./outputs/validated.yaml", "w") as f:
    yaml.dump({"sequences": sequences}, f, default_flow_style=False)

print("Wrote ./outputs/validated.yaml", flush=True)
