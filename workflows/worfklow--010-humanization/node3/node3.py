"""
Node 3 – ABB3 Structure Prediction
====================================
Loads each .pt file from Node 1 (plain ABB3) or Node 2 (ABB3-LM),
runs the model forward pass, reconstructs atom37 coordinates,
and writes one PDB file per antibody pair.

Inputs
------
Node 1 outputs (plain ABB3):
    results/node1/<pair_id>.pt  →  contains "ab_input"

Node 2 outputs (ABB3-LM):
    results/node2/<pair_id>.pt  →  contains "ab_input" + "plm_embedding"

Outputs
-------
results/node3/<pair_id>.pdb
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict

import torch


# ---------------------------------------------------------------------------
# Inference helper
# ---------------------------------------------------------------------------

def run_inference(
    ab_input: Dict[str, Any],
    checkpoint_path: Path,
    out_pdb: Path,
    device: torch.device,
) -> None:
    """Run ABB3 inference for a single antibody and write a PDB file."""
    from abodybuilder3.lightning_module import LitABB3
    from abodybuilder3.utils import add_atom37_to_output, output_to_pdb

    # Load model from checkpoint
    module = LitABB3.load_from_checkpoint(
        checkpoint_path,
        weights_only=False,
        map_location=device,
    )
    model = module.model.to(device).eval()

    # Build batched input (add batch dimension where needed)
    ab_input_batch = {
        k: (v.unsqueeze(0).to(device) if k not in {"single", "pair"} else v.to(device))
        for k, v in ab_input.items()
    }

    # Forward pass
    with torch.no_grad():
        output = model(ab_input_batch, ab_input_batch["aatype"])

    # Reconstruct full atom37 representation
    output = add_atom37_to_output(output, ab_input["aatype"].to(device))

    # Convert to PDB string and write
    pdb_string = output_to_pdb(output, ab_input)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text(pdb_string)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Node 3: ABB3 / ABB3-LM structure prediction."
    )
    parser.add_argument(
        "--inputs", required=True,
        help="Directory of input .pt files (from Node 1 or Node 2)"
    )
    parser.add_argument(
        "--checkpoint", required=True,
        help="Path to ABB3 model checkpoint (.ckpt)"
    )
    parser.add_argument(
        "--outputs", required=True,
        help="Directory to write predicted PDB files"
    )
    parser.add_argument(
        "--device", default="auto", choices=["auto", "cpu", "cuda"],
        help="Compute device (default: auto)"
    )
    args = parser.parse_args()

    inputs_dir = Path(args.inputs)
    outputs_dir = Path(args.outputs)
    checkpoint = Path(args.checkpoint)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    print(f"[Node 3] Device:     {device}")
    print(f"[Node 3] Checkpoint: {checkpoint}")

    if not checkpoint.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint}")

    pt_files = sorted(inputs_dir.glob("*.pt"))
    if not pt_files:
        raise RuntimeError(f"No .pt files found in {inputs_dir}")

    print(f"[Node 3] Found {len(pt_files)} input file(s).")

    # Allow ml_collections ConfigDict to be deserialized safely
    try:
        import ml_collections
        torch.serialization.add_safe_globals(
            [ml_collections.config_dict.config_dict.ConfigDict]
        )
    except ImportError:
        pass  # not available; hope the checkpoint doesn't need it

    for pt_file in pt_files:
        obj: Dict[str, Any] = torch.load(pt_file, map_location="cpu")
        pair_id: str = obj["id"]
        ab_input: Dict[str, Any] = obj["ab_input"]

        # ABB3-LM: attach PLM embedding as the "single" residue feature
        if "plm_embedding" in obj:
            ab_input["single"] = obj["plm_embedding"]
            print(f"[Node 3] {pair_id}: using PLM embedding (ABB3-LM mode)")
        else:
            print(f"[Node 3] {pair_id}: no PLM embedding found (plain ABB3 mode)")

        out_pdb = outputs_dir / f"{pair_id}.pdb"
        print(f"[Node 3] Predicting structure → {out_pdb.name} …")

        run_inference(
            ab_input=ab_input,
            checkpoint_path=checkpoint,
            out_pdb=out_pdb,
            device=device,
        )

        print(f"[Node 3] Saved → {out_pdb}")

    print(f"[Node 3] Done. {len(pt_files)} PDB(s) written to {outputs_dir}")


if __name__ == "__main__":
    main()