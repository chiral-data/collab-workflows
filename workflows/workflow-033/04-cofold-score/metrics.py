"""Structural metrics for binder-design scoring.

Kabsch RMSD + mmCIF CA-coord extraction (via gemmi) for self-consistency.
Ported/extended from NVIDIA BioNeMo Agent Toolkit (Apache-2.0 / CC-BY-4.0).
"""
from __future__ import annotations

import numpy as np


def kabsch_rmsd(p, q) -> float:
    """Minimal RMSD after optimal superposition of two (N, 3) coord arrays."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if p.shape != q.shape or p.ndim != 2 or p.shape[1] != 3:
        raise ValueError(f"coordinate shape mismatch: {p.shape} vs {q.shape}")
    if p.shape[0] == 0:
        raise ValueError("no coordinates provided")
    pc = p - p.mean(axis=0)
    qc = q - q.mean(axis=0)
    h = pc.T @ qc
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    p_rot = pc @ rot.T
    return float(np.sqrt(np.sum((p_rot - qc) ** 2) / p.shape[0]))


def ca_rmsd_from_pdb(pdb_a: str, pdb_b: str,
                     chain_a: str | None = None,
                     chain_b: str | None = None) -> float:
    """CA-RMSD between two PDB structures. Truncates to common residue count."""
    from pdb_utils import ca_coords
    a = ca_coords(pdb_a, chain_a)
    b = ca_coords(pdb_b, chain_b)
    n = min(len(a), len(b))
    if n == 0:
        raise ValueError("no CA atoms found for RMSD computation")
    return kabsch_rmsd(a[:n], b[:n])


def ca_coords_from_cif(cif_content: str, chain_id: str) -> list[tuple[float, float, float]]:
    """Extract CA Cartesian coordinates for a chain from a Boltz2 mmCIF string.

    Uses auth_asym_id first (matches the 'id' field we sent in the request),
    falls back to label_asym_id.
    """
    import gemmi
    doc = gemmi.cif.read_string(cif_content)
    block = doc.sole_block()

    # Try auth_asym_id (matches the chain ID we requested)
    for chain_col in ("auth_asym_id", "label_asym_id"):
        try:
            table = block.find("_atom_site.", [
                chain_col, "label_atom_id",
                "Cartn_x", "Cartn_y", "Cartn_z",
            ])
            coords = [
                (float(row[2]), float(row[3]), float(row[4]))
                for row in table
                if row[0] == chain_id and row[1] == "CA"
            ]
            if coords:
                return coords
        except (KeyError, RuntimeError):
            continue

    raise ValueError(f"no CA atoms for chain {chain_id!r} found in mmCIF")


def mean_plddt_from_cif(cif_content: str, chain_id: str) -> float:
    """Mean pLDDT for a chain from a Boltz2 mmCIF (B_iso_or_equiv = pLDDT).

    Only CA atoms are averaged to get a per-residue mean (avoids atom-count bias).
    """
    import gemmi
    doc = gemmi.cif.read_string(cif_content)
    block = doc.sole_block()

    for chain_col in ("auth_asym_id", "label_asym_id"):
        try:
            table = block.find("_atom_site.", [
                chain_col, "label_atom_id", "B_iso_or_equiv",
            ])
            vals = [
                float(row[2])
                for row in table
                if row[0] == chain_id and row[1] == "CA"
            ]
            if vals:
                return sum(vals) / len(vals)
        except (KeyError, RuntimeError):
            continue

    raise ValueError(f"no CA atoms for chain {chain_id!r} in mmCIF (pLDDT)")


def self_consistency_rmsd(
    cif_content: str,
    backbone_pdb: str,
    binder_cif_chain: str,
    binder_pdb_chain: str,
) -> float:
    """Kabsch CA-RMSD between predicted binder (mmCIF) and RFdiffusion backbone (PDB).

    Low RMSD indicates the designed sequence is predicted to fold back into
    the diffused backbone geometry (Bennett et al. 2023 filter).
    """
    from pdb_utils import ca_coords
    predicted = ca_coords_from_cif(cif_content, binder_cif_chain)
    backbone  = ca_coords(backbone_pdb, binder_pdb_chain)
    n = min(len(predicted), len(backbone))
    if n == 0:
        raise ValueError("no CA atoms to compare for self-consistency RMSD")
    return kabsch_rmsd(predicted[:n], backbone[:n])
