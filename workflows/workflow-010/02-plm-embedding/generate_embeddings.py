"""
Node 2 – PLM Embedding Generation (ProtT5)
==========================================
Loads each .pt file produced by Node 1, generates a per-residue
ProtT5 protein language model embedding for the concatenated
heavy+light sequence, and saves an enriched .pt file that Node 3
can consume for the ABB3-LM variant.

This node is OPTIONAL.  If you want to run plain ABB3 (no language
model), skip directly from Node 1 → Node 3.

Outputs
-------
outputs/<pair_id>.pt
    Same dict as Node 1 output, extended with:
        {
            ...                         # everything from node1
            "plm_embedding": Tensor,    # shape (L, 1024)  float32
        }

Notes
-----
* ProtT5 is large (~3 GB).  A GPU is strongly recommended.
* If a pre-computed embedding .pt file lives in --precomputed_dir,
  it is loaded directly instead of recomputing.

Fixes vs original
-----------------
* ProtT5 wrapper does not expose .to() / .eval() directly — those
  calls are now made on the inner HuggingFace model instead.
* Embedding API uses get_embeddings([heavy], [light]) with separate
  lists matching the upstream ABodyBuilder3 API contract.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, Optional

import torch


# ---------------------------------------------------------------------------
# PLM helpers
# ---------------------------------------------------------------------------

def load_or_compute_embedding(
    pair_id: str,
    heavy: str,
    light: str,
    precomputed_dir: Optional[Path],
    model,          # ProtT5 instance (lazy-loaded)
    device: torch.device,
) -> torch.Tensor:
    """
    Return a PLM embedding tensor of shape (L, D).

    Looks for a pre-computed file at
        <precomputed_dir>/<pair_id>.pt
    before falling back to on-the-fly ProtT5 inference.
    """
    if precomputed_dir is not None:
        candidate = precomputed_dir / f"{pair_id}.pt"
        if candidate.exists():
            cached = torch.load(candidate, map_location="cpu")
            embedding = cached["plm_embedding"]
            print(f"[Node 2]   Loaded cached embedding from {candidate}  shape={tuple(embedding.shape)}")
            return embedding

    if model is None:
        raise RuntimeError(
            "ProtT5 model is not loaded and no pre-computed embedding was found. "
            "Either provide --precomputed_dir or ensure ProtT5 can be imported."
        )

    print(f"[Node 2]   Computing ProtT5 embedding for {pair_id} …")

    # Use the upstream API: separate heavy and light lists.
    # This matches ABodyBuilder3's get_embeddings([heavy_seqs], [light_seqs]) signature.
    with torch.no_grad():
        embedding = model.get_embeddings([heavy], [light])[0]   # (L, D)

    print(f"[Node 2]   Computed embedding shape={tuple(embedding.shape)}")
    return embedding.cpu()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Node 2: Generate ProtT5 PLM embeddings for ABB3-LM."
    )
    parser.add_argument(
        "--inputs", required=True,
        help="Directory of Node 1 .pt files"
    )
    parser.add_argument(
        "--outputs", required=True,
        help="Directory to write enriched .pt files"
    )
    parser.add_argument(
        "--precomputed_dir", default=None,
        help="(Optional) Directory containing pre-computed <pair_id>.pt embedding files"
    )
    parser.add_argument(
        "--device", default="auto", choices=["auto", "cpu", "cuda"],
        help="Compute device for ProtT5 inference (default: auto)"
    )
    parser.add_argument(
        "--use-plm", default="0",
        help="Set to '1' to compute PLM embeddings; '0' to passthrough inputs unchanged"
    )
    args = parser.parse_args()

    inputs_dir = Path(args.inputs)
    outputs_dir = Path(args.outputs)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    pt_files = sorted(inputs_dir.glob("*.pt"))
    if not pt_files:
        raise RuntimeError(f"No .pt files found in {inputs_dir}")

    use_plm = args.use_plm.strip() == "1"

    if not use_plm:
        print(f"[Node 2] USE_PLM=0 — passthrough mode (copying {len(pt_files)} file(s))")
        for pt_file in pt_files:
            out_path = outputs_dir / pt_file.name
            import shutil
            shutil.copy2(pt_file, out_path)
            print(f"[Node 2]   Copied → {out_path}")
        print(f"[Node 2] Done. {len(pt_files)} file(s) passed through to {outputs_dir}")
        return

    print(f"[Node 2] USE_PLM=1 — computing ProtT5 embeddings")

    precomputed_dir = Path(args.precomputed_dir) if args.precomputed_dir else None

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    print(f"[Node 2] Device: {device}")
    print(f"[Node 2] Found {len(pt_files)} input file(s).")

    plm_model = None

    def get_plm_model():
        nonlocal plm_model
        if plm_model is not None:
            return plm_model
        print("[Node 2] Loading ProtT5 model (this may take a while) …")
        try:
            from abodybuilder3.language.model import ProtT5
            _hub = Path("/root/.cache/huggingface/hub/models--Rostlab--prot_t5_xl_half_uniref50-enc")
            _rev = (_hub / "refs" / "main").read_text().strip()
            weights_path = str(_hub / "snapshots" / _rev)
            plm_model = ProtT5(weights_dir=weights_path)
            if hasattr(plm_model, "model"):
                plm_model.model.to(device).eval()
            elif hasattr(plm_model, "encoder"):
                plm_model.encoder.to(device).eval()
            else:
                try:
                    plm_model.to(device)
                    plm_model.eval()
                except AttributeError:
                    pass
        except ImportError as e:
            raise ImportError(
                "abodybuilder3 (with ProtT5) is not installed or not on the Python path."
            ) from e
        return plm_model

    for pt_file in pt_files:
        obj: Dict[str, Any] = torch.load(pt_file, map_location="cpu")
        pair_id: str = obj["id"]
        heavy: str = obj["heavy"]
        light: str = obj["light"]

        print(f"[Node 2] Processing: {pair_id}")

        need_model = precomputed_dir is None or not (precomputed_dir / f"{pair_id}.pt").exists()
        model_instance = get_plm_model() if need_model else None

        embedding = load_or_compute_embedding(
            pair_id=pair_id,
            heavy=heavy,
            light=light,
            precomputed_dir=precomputed_dir,
            model=model_instance,
            device=device,
        )

        out_obj = {**obj, "plm_embedding": embedding}
        out_path = outputs_dir / f"{pair_id}.pt"
        torch.save(out_obj, out_path)
        print(f"[Node 2]   Saved → {out_path}")

    print(f"[Node 2] Done. {len(pt_files)} file(s) written to {outputs_dir}")


if __name__ == "__main__":
    main()