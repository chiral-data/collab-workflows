#!/usr/bin/env python3
"""Node 05 — Analyze Ensemble.

Turn the sampled ensemble into numbers: a per-frame table, ensemble-level statistics, and a
multi-model PDB of representative frames for the report's 3D viewer.

Every quantity here is **geometric**. Nothing in this node is a confidence estimate, and
nothing is a thermodynamic quantity: upstream is explicit that BioEmu2's sampling frequency
is not a calibrated Boltzmann probability, so a cluster population is a count of frames, not
a free energy, and this node deliberately computes no ΔG.

The trajectory is bound to the System's own topology (`system.topology`) rather than parsed
from topology.cif. The XTC holds coordinates only, and mdtraj's PDBx/mmCIF support is the
weak link in that chain; the System is authoritative on both the sampling and the analysis
side, so atom order matches by construction.

Per-frame: radius of gyration, C-alpha RMSD to the ensemble medoid, physicality status.
Ensemble: Rg spread, pairwise C-alpha RMSD matrix, per-residue RMSF, cluster count at a
given RMSD threshold, physicality pass rate. Propagated frames are analysed as a second
ensemble only when node 04 actually propagated.
"""
import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("HF_HUB_OFFLINE", "1")

COLUMNS = [
    "frame",
    "ensemble",
    "radius_of_gyration_nm",
    "rmsd_to_medoid_nm",
    "physicality_status",
    "physicality_reason",
]

# A pairwise RMSD matrix is O(n^2); above this many frames it is summarised, not stored.
MAX_MATRIX_FRAMES = 200


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def read_json(path, label, default=None):
    if not os.path.exists(path):
        return default if default is not None else {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        print(f"WARNING: could not read {label}: {exc}", flush=True)
        return default if default is not None else {}


def load_trajectory(xtc_path, topology, label):
    import mdtraj

    try:
        traj = mdtraj.load_xtc(xtc_path, top=topology)
    except Exception as exc:
        fail(f"Could not load the {label} ensemble from '{xtc_path}': {exc}")
    if traj.n_frames == 0:
        fail(f"The {label} ensemble contains no frames.")
    return traj


def ca_indices(traj):
    """C-alpha atom indices, falling back to all atoms for a non-protein system."""
    try:
        indices = traj.topology.select("name CA")
    except Exception:
        indices = []
    if len(indices) >= 3:
        return indices, True
    print(
        "NOTE: fewer than 3 C-alpha atoms found; RMSD and RMSF fall back to all atoms. "
        "Expected for a nucleic-acid or ligand-only system.",
        flush=True,
    )
    return list(range(traj.n_atoms)), False


def pairwise_rmsd(traj, indices):
    """Symmetric pairwise RMSD in nm, each row superposed on the corresponding frame."""
    import mdtraj
    import numpy as np

    n = traj.n_frames
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        matrix[i] = mdtraj.rmsd(traj, traj, frame=i, atom_indices=indices)
    # mdtraj.rmsd is already symmetric up to float error; average to make it exactly so.
    return (matrix + matrix.T) / 2.0


def greedy_clusters(matrix, threshold):
    """Count clusters by greedy medoid assignment: the simplest defensible cluster count.

    Repeatedly take the unassigned frame with the most unassigned neighbours within
    `threshold` as a cluster centre, and absorb its neighbours. Order-independent given the
    matrix, and needs no extra dependency.
    """
    import numpy as np

    n = matrix.shape[0]
    unassigned = set(range(n))
    clusters = []
    while unassigned:
        members = list(unassigned)
        neighbour_counts = [
            (int((matrix[i, members] <= threshold).sum()), i) for i in members
        ]
        _, centre = max(neighbour_counts)
        absorbed = [i for i in members if matrix[centre, i] <= threshold]
        clusters.append({"centre": int(centre), "size": len(absorbed),
                         "members": [int(i) for i in absorbed]})
        unassigned -= set(absorbed)
    clusters.sort(key=lambda c: c["size"], reverse=True)
    return clusters


def analyze_ensemble(traj, label, physicality, threshold):
    """Return (per_frame_rows, summary_dict, medoid_index, rmsd_matrix_or_None)."""
    import mdtraj
    import numpy as np

    indices, used_ca = ca_indices(traj)
    rg = mdtraj.compute_rg(traj)

    matrix = None
    if traj.n_frames <= MAX_MATRIX_FRAMES:
        matrix = pairwise_rmsd(traj, indices)
        medoid = int(np.argmin(matrix.sum(axis=1)))
        rmsd_to_medoid = matrix[medoid]
        clusters = greedy_clusters(matrix, threshold)
        upper = matrix[np.triu_indices(traj.n_frames, k=1)]
        pairwise_stats = {
            "mean_nm": round(float(upper.mean()), 4) if upper.size else 0.0,
            "max_nm": round(float(upper.max()), 4) if upper.size else 0.0,
            "min_nm": round(float(upper.min()), 4) if upper.size else 0.0,
        }
    else:
        # Too many frames for an n^2 matrix: use frame 0 as the reference instead.
        print(
            f"{traj.n_frames} frames exceeds {MAX_MATRIX_FRAMES}; skipping the pairwise "
            "RMSD matrix and clustering, and measuring RMSD against frame 0.",
            flush=True,
        )
        medoid = 0
        rmsd_to_medoid = mdtraj.rmsd(traj, traj, frame=0, atom_indices=indices)
        clusters = None
        pairwise_stats = None

    # Per-residue RMSF about the mean structure, superposed on the medoid first.
    superposed = traj[:]
    try:
        superposed.superpose(traj, frame=medoid, atom_indices=indices)
        rmsf = mdtraj.rmsf(superposed, None, atom_indices=indices)
    except Exception as exc:
        print(f"WARNING: could not compute RMSF: {exc}", flush=True)
        rmsf = None

    residue_rmsf = []
    if rmsf is not None:
        for position, atom_index in enumerate(indices):
            atom = traj.topology.atom(int(atom_index))
            residue_rmsf.append(
                {
                    "chain": atom.residue.chain.index,
                    "resSeq": int(atom.residue.resSeq),
                    "resName": atom.residue.name,
                    "rmsf_nm": round(float(rmsf[position]), 4),
                }
            )

    statuses = physicality.get("status", [])
    reasons = physicality.get("reason", [])

    rows = []
    for frame in range(traj.n_frames):
        rows.append(
            {
                "frame": frame,
                "ensemble": label,
                "radius_of_gyration_nm": round(float(rg[frame]), 4),
                "rmsd_to_medoid_nm": round(float(rmsd_to_medoid[frame]), 4),
                "physicality_status": statuses[frame] if frame < len(statuses) else "",
                "physicality_reason": reasons[frame] if frame < len(reasons) else "",
            }
        )

    n_passed = sum(1 for s in statuses if s == "passed")
    n_checked = sum(1 for s in statuses if s in ("passed", "rejected"))

    summary = {
        "label": label,
        "n_frames": int(traj.n_frames),
        "n_atoms": int(traj.n_atoms),
        "rmsd_atom_selection": "C-alpha" if used_ca else "all atoms",
        "radius_of_gyration_nm": {
            "mean": round(float(rg.mean()), 4),
            "std": round(float(rg.std()), 4),
            "min": round(float(rg.min()), 4),
            "max": round(float(rg.max()), 4),
        },
        "medoid_frame": medoid,
        "pairwise_rmsd": pairwise_stats,
        "rmsd_cluster_threshold_nm": threshold,
        "n_clusters": len(clusters) if clusters is not None else None,
        "clusters": clusters[:20] if clusters is not None else None,
        "physicality": {
            "n_passed": n_passed,
            "n_checked": n_checked,
            "pass_rate": round(n_passed / n_checked, 3) if n_checked else None,
        },
        "residue_rmsf": residue_rmsf,
    }
    return rows, summary, medoid, matrix


def parse_physicality(path):
    """physicality.json is a list of {frame, status, reason}; return parallel lists."""
    raw = read_json(path, "physicality report", default=[])
    entries = raw.get("frames", raw) if isinstance(raw, dict) else raw
    if not isinstance(entries, list):
        return {"status": [], "reason": []}
    ordered = sorted(
        (e for e in entries if isinstance(e, dict)),
        key=lambda e: e.get("frame", 0),
    )
    return {
        "status": [str(e.get("status", "")) for e in ordered],
        "reason": ["" if e.get("reason") is None else str(e.get("reason")) for e in ordered],
    }


def select_frames(summary, n_frames, top_n):
    """Representative frames: cluster centres by size first, then medoid-nearest."""
    clusters = summary.get("clusters")
    chosen = []
    if clusters:
        chosen = [c["centre"] for c in clusters][:top_n]
    if len(chosen) < top_n:
        for frame in range(n_frames):
            if frame not in chosen:
                chosen.append(frame)
            if len(chosen) >= top_n:
                break
    return chosen[: min(top_n, n_frames)]


def main():
    os.makedirs("./outputs", exist_ok=True)

    threshold = float(os.environ.get("PARAM_RMSD_CLUSTER_THRESHOLD_NM", "0.2"))
    top_n = int(os.environ.get("PARAM_TOP_N", "10"))

    if threshold <= 0:
        fail(f"rmsd_cluster_threshold_nm must be positive, got {threshold}")
    if top_n < 1:
        fail(f"top_n must be at least 1, got {top_n}")

    for name in ("samples.xtc", "system.yaml"):
        if not os.path.exists(f"./inputs/{name}"):
            fail(f"Required input not found: ./inputs/{name}")

    from bioemu2 import System

    try:
        system = System.from_yaml(Path("./inputs/system.yaml"))
    except Exception as exc:
        fail(f"Could not load system YAML: {exc}")
    topology = system.topology

    sampling_manifest = read_json("./inputs/sampling_manifest.json", "sampling manifest")
    propagation_manifest = read_json("./inputs/propagation_manifest.json", "propagation manifest")
    physicality = parse_physicality("./inputs/physicality.json")

    equilibrium = load_trajectory("./inputs/samples.xtc", topology, "sampled")
    print(f"Sampled ensemble: {equilibrium.n_frames} frame(s), {equilibrium.n_atoms} atoms",
          flush=True)

    rows, summary, medoid, matrix = analyze_ensemble(
        equilibrium, "sampled", physicality, threshold
    )

    summaries = {"sampled": summary}
    all_rows = list(rows)

    propagated_enabled = bool(propagation_manifest.get("enabled"))
    propagated_path = "./inputs/propagated_samples.xtc"
    if propagated_enabled and os.path.exists(propagated_path):
        propagated = load_trajectory(propagated_path, topology, "propagated")
        print(f"Propagated ensemble: {propagated.n_frames} frame(s)", flush=True)
        prop_rows, prop_summary, _, _ = analyze_ensemble(
            propagated,
            "propagated",
            {"status": [], "reason": []},  # node 04 runs no physicality check
            threshold,
        )
        prop_summary["lag_time_ns"] = propagation_manifest.get("lag_time_ns")
        summaries["propagated"] = prop_summary
        all_rows.extend(prop_rows)
    else:
        print(
            "Propagation was disabled — analysing the sampled ensemble only. Node 04's "
            "output is an identity copy and is not treated as time-lagged.",
            flush=True,
        )

    with open("./outputs/ensemble_metrics.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(all_rows)

    # Multi-model PDB of representative frames, for the report's 3D viewer.
    selected = select_frames(summary, equilibrium.n_frames, top_n)
    try:
        equilibrium[selected].save_pdb("./outputs/ensemble_top.pdb")
    except Exception as exc:
        fail(f"Could not write the multi-model PDB: {exc}")

    payload = {
        "ensembles": summaries,
        # The report's heatmap reads this directly; it is None above MAX_MATRIX_FRAMES, where
        # an n^2 matrix would bloat both this file and the HTML.
        "rmsd_matrix": (
            [[round(float(v), 4) for v in row] for row in matrix]
            if matrix is not None
            else None
        ),
        "viewer_frames": [int(f) for f in selected],
        "propagation_enabled": propagated_enabled,
        "mode": (sampling_manifest.get("settings") or {}).get("mode"),
        "sampler": (sampling_manifest.get("settings") or {}).get("sampler"),
        "seed": (sampling_manifest.get("settings") or {}).get("seed"),
        "bioemu2_version": sampling_manifest.get("bioemu2_version"),
        "msa_source": sampling_manifest.get("msa_source"),
        "msa_min_depth": sampling_manifest.get("msa_min_depth"),
        "metrics_are_geometric_only": True,
        "caveat": (
            "Sampling frequency is not a calibrated Boltzmann probability. Cluster sizes are "
            "frame counts, not populations or free energies. Sampled frames are independent "
            "draws, not a trajectory."
        ),
    }
    with open("./outputs/ensemble_summary.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    rg = summary["radius_of_gyration_nm"]
    print(
        f"Rg {rg['mean']} +/- {rg['std']} nm (range {rg['min']}-{rg['max']})",
        flush=True,
    )
    if summary["pairwise_rmsd"]:
        print(
            f"Pairwise C-alpha RMSD: mean {summary['pairwise_rmsd']['mean_nm']} nm, "
            f"max {summary['pairwise_rmsd']['max_nm']} nm",
            flush=True,
        )
        if summary["pairwise_rmsd"]["max_nm"] < 1e-3:
            print(
                "WARNING: the ensemble has collapsed — every frame is effectively identical. "
                "Check mode (equilibrium vs structure) and the MSA depth.",
                flush=True,
            )
    if summary["n_clusters"] is not None:
        print(
            f"{summary['n_clusters']} cluster(s) at {threshold} nm; medoid frame {medoid}",
            flush=True,
        )
    print(f"Wrote {len(selected)} representative frame(s) to ensemble_top.pdb", flush=True)


if __name__ == "__main__":
    main()
