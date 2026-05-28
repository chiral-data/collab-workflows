import json
import os
import torch
import yaml

from esm.models.esmfold2 import (
    DNAInput,
    ESMFold2InputBuilder,
    LigandInput,
    ProteinInput,
    RNAInput,
    StructurePredictionInput,
)
from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

print("Loading ESMFold2 model...", flush=True)
model = ESMFold2Model.from_pretrained(
    "biohub/ESMFold2",
    dtype=torch.bfloat16,
    device_map="auto",
).eval()
print("Model loaded.", flush=True)

data = yaml.safe_load(open("./inputs/validated.yaml"))
chain_inputs = []
for s in data["sequences"]:
    t = s["type"]
    if t == "protein":
        chain_inputs.append(ProteinInput(id=s["id"], sequence=s["sequence"]))
    elif t == "rna":
        chain_inputs.append(RNAInput(id=s["id"], sequence=s["sequence"]))
    elif t == "dna":
        chain_inputs.append(DNAInput(id=s["id"], sequence=s["sequence"]))
    elif t == "ligand":
        chain_inputs.append(LigandInput(id=s["id"], ccd=s.get("ccd"), smiles=s.get("smiles")))

spi = StructurePredictionInput(sequences=chain_inputs)

num_loops            = int(os.environ.get("PARAM_NUM_LOOPS", 3))
num_sampling_steps   = int(os.environ.get("PARAM_NUM_SAMPLING_STEPS", 50))
num_diffusion_samples = int(os.environ.get("PARAM_NUM_DIFFUSION_SAMPLES", 1))

print(f"Running prediction: loops={num_loops}, sampling_steps={num_sampling_steps}, samples={num_diffusion_samples}", flush=True)

result = ESMFold2InputBuilder().fold(
    model, spi,
    num_loops=num_loops,
    num_sampling_steps=num_sampling_steps,
    num_diffusion_samples=num_diffusion_samples,
    seed=0,
)

os.makedirs("./outputs", exist_ok=True)

with open("./outputs/structure.cif", "w") as f:
    f.write(result.complex.to_mmcif())

confidence = {
    "plddt_mean":        float(result.plddt.mean()),
    "ptm":               float(result.ptm),
    "iptm":              float(result.iptm),
    "plddt_per_residue": result.plddt.cpu().tolist(),
    "pae":               result.pae.cpu().tolist() if hasattr(result, "pae") and result.pae is not None else None,
}
with open("./outputs/confidence.json", "w") as f:
    json.dump(confidence, f, indent=2)

print(f"Done — pLDDT: {confidence['plddt_mean']:.3f}  pTM: {confidence['ptm']:.3f}  ipTM: {confidence['iptm']:.3f}", flush=True)
