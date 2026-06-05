#!/usr/bin/env python3
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def render_plddt_pymol(pdb_path, plddt_per_residue, out_path):
    import pymol
    from pymol import cmd as pm

    pm.reinitialize()
    pm.load(pdb_path, "structure")

    # Assign pLDDT values to B-factors
    stored_plddt = plddt_per_residue[:]
    pm.alter("structure and name CA", "b = stored_plddt.pop(0)", space={"stored_plddt": stored_plddt})
    pm.spectrum("b", "red_yellow_cyan_blue", "structure", minimum=0, maximum=100)

    pm.set("ray_opaque_background", 0)
    pm.set("antialias", 2)
    pm.orient("structure")
    pm.png(out_path, width=1200, height=900, ray=1)
    pm.delete("all")
    print(f"PyMOL render saved: {out_path}", flush=True)


def render_plddt_fallback(pdb_path, plddt_per_residue, out_path):
    """Matplotlib color-scale bar as fallback when PyMOL is unavailable."""
    print("WARNING: PyMOL not available — generating pLDDT color scale figure as fallback", flush=True)
    fig, ax = plt.subplots(figsize=(8, 2))
    gradient = np.linspace(0, 100, 256).reshape(1, -1)
    ax.imshow(gradient, aspect="auto", cmap="RdYlBu", extent=[0, 100, 0, 1])
    ax.set_yticks([])
    ax.set_xlabel("pLDDT score")
    ax.set_title(f"pLDDT color scale (PDB: {os.path.basename(pdb_path)})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Fallback pLDDT figure saved: {out_path}", flush=True)


def render_pae_heatmap(pae_matrix, max_pae, chain_lengths, out_path):
    pae = np.array(pae_matrix)
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(pae, cmap="RdYlGn_r", vmin=0, vmax=max_pae, origin="upper")
    plt.colorbar(im, ax=ax, label="PAE (Å)")
    ax.set_title("Predicted Aligned Error")
    ax.set_xlabel("Scored residue")
    ax.set_ylabel("Aligned residue")

    # Annotate chain boundaries for multimers
    if len(chain_lengths) > 1:
        cumulative = 0
        for length in chain_lengths[:-1]:
            cumulative += length
            ax.axvline(x=cumulative - 0.5, color="black", linewidth=1)
            ax.axhline(y=cumulative - 0.5, color="black", linewidth=1)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"PAE heatmap saved: {out_path}", flush=True)


def main():
    os.makedirs("./outputs", exist_ok=True)

    color_by = os.environ.get("PARAM_COLOR_BY", "pLDDT")
    renderer = os.environ.get("PARAM_RENDERER", "auto")

    # Find top-ranked PDB
    pdb_files = glob.glob("./inputs/*_unrelaxed_rank_001_*.pdb")
    if not pdb_files:
        print("ERROR: No rank_001 PDB file found in ./inputs/", flush=True)
        sys.exit(1)
    pdb_path = pdb_files[0]
    print(f"Using structure: {pdb_path}", flush=True)

    # Load confidence summary
    summary_path = "./inputs/confidence_summary.json"
    if not os.path.exists(summary_path):
        print(f"ERROR: confidence_summary.json not found at {summary_path}", flush=True)
        sys.exit(1)

    with open(summary_path) as f:
        summary = json.load(f)

    rank1 = summary["models"][0]
    plddt_per_residue = rank1["plddt_per_residue"]
    pae_matrix = rank1["pae"]
    max_pae = rank1["max_pae"]
    chain_lengths = [len(plddt_per_residue)]  # monomer default

    # For multimers, reconstruct chain lengths from per-chain pLDDT if available
    if summary.get("mode") == "multimer":
        # Best effort: distribute residues equally if chain info not stored
        # (chain_lengths are already embedded in validated_sequences.fasta via node 03)
        pass

    # Structure render
    struct_out = "./outputs/structure_plddt.png"
    if renderer != "py3dmol":
        try:
            import pymol  # noqa: F401
            render_plddt_pymol(pdb_path, plddt_per_residue, struct_out)
        except ImportError:
            render_plddt_fallback(pdb_path, plddt_per_residue, struct_out)
    else:
        render_plddt_fallback(pdb_path, plddt_per_residue, struct_out)

    # PAE heatmap — always matplotlib
    pae_out = "./outputs/pae_matrix.png"
    render_pae_heatmap(pae_matrix, max_pae, chain_lengths, pae_out)


if __name__ == "__main__":
    main()
