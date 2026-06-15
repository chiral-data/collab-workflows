"""
Node 3 – ABB3 Structure Prediction
====================================
Loads each .pt file from Node 1 (plain ABB3) or Node 2 (ABB3-LM),
runs the model forward pass, reconstructs atom37 coordinates,
and writes one PDB file per antibody pair.

Inputs
------
Node 1 outputs (plain ABB3):
    inputs/<pair_id>.pt  →  contains "ab_input"

Node 2 outputs (ABB3-LM):
    inputs/<pair_id>.pt  →  contains "ab_input" + "plm_embedding"

Outputs
-------
outputs/<pair_id>.pdb
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
    model,
    out_pdb: Path,
    device: torch.device,
) -> None:
    """Run ABB3 inference for a single antibody and write a PDB file."""
    from abodybuilder3.utils import add_atom37_to_output, output_to_pdb

    # string_to_input already adds a batch dim to "single" (3D) and "pair" (4D);
    # only 1D scalars and the injected 2D PLM embedding need unsqueeze(0).
    ab_input_batch = {
        k: (v.unsqueeze(0).to(device) if v.dim() < 3 else v.to(device))
        if isinstance(v, torch.Tensor) else v
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
        help="Path to plain ABB3 checkpoint (.ckpt); auto-replaced by --checkpoint_lm when PLM embeddings are present"
    )
    parser.add_argument(
        "--checkpoint_lm", default=None,
        help="Path to ABB3-LM checkpoint (.ckpt); used when input .pt files contain 'plm_embedding'"
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
    outputs_dir.mkdir(parents=True, exist_ok=True)

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    pt_files = sorted(inputs_dir.glob("*.pt"))
    if not pt_files:
        raise RuntimeError(f"No .pt files found in {inputs_dir}")

    print(f"[Node 3] Found {len(pt_files)} input file(s).")

    # Detect mode from first input file
    first_obj: Dict[str, Any] = torch.load(pt_files[0], map_location="cpu")
    use_lm = "plm_embedding" in first_obj

    if use_lm:
        if args.checkpoint_lm is None:
            raise ValueError("Input contains PLM embeddings (ABB3-LM mode) but --checkpoint_lm was not provided.")
        checkpoint = Path(args.checkpoint_lm)
        print(f"[Node 3] Mode:       ABB3-LM (PLM embeddings detected)")
    else:
        checkpoint = Path(args.checkpoint)
        print(f"[Node 3] Mode:       plain ABB3")

    print(f"[Node 3] Device:     {device}")
    print(f"[Node 3] Checkpoint: {checkpoint}")

    if not checkpoint.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint}")

    # Allow ml_collections ConfigDict to be deserialized safely (PyTorch >= 2.4)
    try:
        import ml_collections
        if hasattr(torch.serialization, "add_safe_globals"):
            torch.serialization.add_safe_globals(
                [ml_collections.config_dict.config_dict.ConfigDict]
            )
    except ImportError:
        pass

    # Load model ONCE outside the loop
    from abodybuilder3.lightning_module import LitABB3

    print("[Node 3] Loading model from checkpoint …")
    module = LitABB3.load_from_checkpoint(
        checkpoint,
        weights_only=False,
        map_location=device,
    )
    model = module.model.to(device).eval()

    for pt_file in pt_files:
        obj: Dict[str, Any] = torch.load(pt_file, map_location="cpu")
        pair_id: str = obj["id"]
        ab_input: Dict[str, Any] = obj["ab_input"]

        if "plm_embedding" in obj:
            ab_input["single"] = obj["plm_embedding"]
            print(f"[Node 3] {pair_id}: ABB3-LM mode")
        else:
            print(f"[Node 3] {pair_id}: plain ABB3 mode")

        out_pdb = outputs_dir / f"{pair_id}.pdb"
        print(f"[Node 3] Predicting structure → {out_pdb.name} …")

        run_inference(
            ab_input=ab_input,
            model=model,
            out_pdb=out_pdb,
            device=device,
        )

        print(f"[Node 3] Saved → {out_pdb}")

    print(f"[Node 3] Done. {len(pt_files)} PDB(s) written to {outputs_dir}")


if __name__ == "__main__":
    main()
