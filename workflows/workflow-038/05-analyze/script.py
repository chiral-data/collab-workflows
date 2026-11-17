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
    """C-alpha atom indices, with the selection actually used reported back.

    `select("name CA")` is wrong here, and wrong in a way that does not look wrong: mdtraj
    matches **atom names**, and a calcium ion is an atom named CA in a residue named CA. Any
    system with structural calcium — calmodulin, a metalloprotease, haemoglobin's relatives —
    would have had its ions folded into the C-alpha set, shifting every RMSD, every RMSF and
    the medoid itself by an amount nothing in the output would reveal.

    Returns (indices, description, excluded_chains).
    """
    topology = traj.topology
    try:
        indices = list(topology.select("protein and name CA"))
    except Exception as exc:
        print(f"WARNING: C-alpha selection failed ({exc}); trying a direct scan.", flush=True)
        indices = []
    how = "protein and name CA"

    if len(indices) < 3:
        # mdtraj's `protein` keyword works from its own residue-name table, so it drops
        # modified residues that BioEmu2 models perfectly well. Fall back to C-alpha atoms in
        # multi-atom residues: that keeps a modified residue and still excludes a lone ion,
        # which is the whole point of not writing "name CA".
        indices = [
            atom.index
            for atom in topology.atoms
            if atom.name == "CA" and atom.residue.n_atoms > 1
        ]
        how = "name CA in residues with more than one atom"
        if len(indices) >= 3:
            print(
                f"NOTE: mdtraj's 'protein' selection found fewer than 3 C-alpha atoms; using "
                f"{len(indices)} C-alpha atom(s) from multi-atom residues instead. Expected "
                "for a chain of modified or non-standard residues.",
                flush=True,
            )

    if len(indices) < 3:
        print(
            "NOTE: fewer than 3 C-alpha atoms found; RMSD and RMSF fall back to all atoms. "
            "Expected for a nucleic-acid or ligand-only system.",
            flush=True,
        )
        return list(range(traj.n_atoms)), "all atoms", []

    selected = set(indices)
    excluded = sorted(
        {
            str(getattr(chain, "chain_id", None) or chain.index)
            for chain in topology.chains
            if not any(atom.index in selected for atom in chain.atoms)
        }
    )
    if excluded:
        print(
            f"Chains contributing no C-alpha atom, and so excluded from RMSD, RMSF and "
            f"clustering: {excluded}. They remain in the frames and in the radius of "
            "gyration, which is computed over every atom.",
            flush=True,
        )
    return indices, how, excluded


def pairwise_pass(traj, indices, keep_matrix):
    """One streaming pass over every frame pair.

    Returns (row_sums, stats, matrix_or_None). The row sums give the medoid, and the stats
    give the ensemble spread, both in O(n) memory. Only the heatmap needs the matrix itself,
    so above `MAX_MATRIX_FRAMES` it is not retained — but the medoid and the spread still
    are. They used to be dropped together, which silently turned "RMSD to the medoid" into
    "RMSD to whichever frame happened to be sampled first", under the same column name.
    """
    import mdtraj
    import numpy as np

    n = traj.n_frames
    matrix = np.zeros((n, n), dtype=float) if keep_matrix else None
    row_sums = np.zeros(n, dtype=float)
    total = 0.0
    largest = 0.0
    smallest = float("inf")

    for i in range(n):
        row = mdtraj.rmsd(traj, traj, frame=i, atom_indices=indices)
        row_sums[i] = float(row.sum())
        if matrix is not None:
            matrix[i] = row
        if i + 1 < n:
            # Upper triangle only, so each pair is counted once.
            upper = row[i + 1 :]
            total += float(upper.sum())
            largest = max(largest, float(upper.max()))
            smallest = min(smallest, float(upper.min()))

    n_pairs = n * (n - 1) // 2
    stats = {
        "mean_nm": round(total / n_pairs, 4) if n_pairs else 0.0,
        "max_nm": round(largest, 4) if n_pairs else 0.0,
        "min_nm": round(smallest, 4) if n_pairs else 0.0,
        "n_pairs": n_pairs,
    }
    if matrix is not None:
        # mdtraj.rmsd is already symmetric up to float error; average to make it exactly so.
        matrix = (matrix + matrix.T) / 2.0
    return row_sums, stats, matrix


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
    """Return (rows, summary, medoid_index, matrix_or_None, superposed_traj)."""
    import mdtraj
    import numpy as np

    indices, selection, excluded_chains = ca_indices(traj)
    rg = mdtraj.compute_rg(traj)

    keep_matrix = traj.n_frames <= MAX_MATRIX_FRAMES
    if not keep_matrix:
        print(
            f"{traj.n_frames} frames exceeds {MAX_MATRIX_FRAMES}; the pairwise RMSD matrix is "
            "not retained, so there is no heatmap and no clustering. The medoid and the "
            "ensemble spread are still measured over every pair.",
            flush=True,
        )
    row_sums, pairwise_stats, matrix = pairwise_pass(traj, indices, keep_matrix)

    medoid = int(np.argmin(row_sums))
    rmsd_to_medoid = (
        matrix[medoid] if matrix is not None else mdtraj.rmsd(traj, traj, frame=medoid, atom_indices=indices)
    )
    clusters = greedy_clusters(matrix, threshold) if matrix is not None else None

    # Per-residue RMSF about the mean structure, superposed on the medoid first. mdtraj's
    # rmsf(target, None) sets trajectory_prealigned=True, so this superposition is not
    # optional: without it the RMSF measures rigid-body drift rather than flexibility.
    superposed = traj[:]
    rmsf = None
    try:
        superposed.superpose(traj, frame=medoid, atom_indices=indices)
        rmsf = mdtraj.rmsf(superposed, None, atom_indices=indices)
    except Exception as exc:
        print(f"WARNING: could not compute RMSF: {exc}", flush=True)

    residue_rmsf = []
    if rmsf is not None:
        for position, atom_index in enumerate(indices):
            atom = traj.topology.atom(int(atom_index))
            chain = atom.residue.chain
            residue_rmsf.append(
                {
                    "chain": str(getattr(chain, "chain_id", None) or chain.index),
                    "chain_index": int(chain.index),
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
        "rmsd_atom_selection": selection,
        "excluded_chains": excluded_chains,
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
        # Truncated for the report only. select_frames() is given the full list, since a
        # top_n above this cut-off would otherwise fall back to consecutive frames without
        # saying so.
        "clusters": clusters[:20] if clusters is not None else None,
        "physicality": {
            "n_passed": n_passed,
            "n_checked": n_checked,
            "pass_rate": round(n_passed / n_checked, 3) if n_checked else None,
        },
        "residue_rmsf": residue_rmsf,
    }
    return rows, summary, medoid, matrix, superposed, clusters


def parse_physicality(path, expected_frames, manifest):
    """Read node 03's physicality.json into parallel lists.

    The file is a bare JSON list of `{"frame", "status", "reason"}`, written by
    `PhysicalityReport.save()` in bioemu2/physicality.py. Status is one of `passed`,
    `rejected` or `not_checked`; `not_checked` is what `no_filter` and a system with no
    eligible protein atoms both produce.

    Node 03 writes this file on every successful run and fails if `save_frames` did not, so
    a missing or short file means the plumbing dropped it, not that there was nothing to
    check. Reporting empty columns in that case would present "we did not check" and "we
    could not read the check" as the same thing.
    """
    if not os.path.exists(path):
        if manifest:
            fail(
                f"Node 03's physicality report is missing from {path}, but its manifest is "
                "present, so node 03 ran and wrote one. Check that node 03 declares "
                "physicality.json as an output and this node declares it as an input. "
                "Continuing would blank the physicality column for every frame and look "
                "like a run where no check was requested."
            )
        fail(f"Required input not found: {path}")

    entries = read_json(path, "physicality report", default=[])
    if not isinstance(entries, list):
        fail(
            f"{path} should be a JSON list of per-frame records, got "
            f"{type(entries).__name__}. This is the format bioemu2's PhysicalityReport.save() "
            "writes, so a different shape means a bioemu2 change rather than a bad input."
        )
    ordered = sorted(
        (e for e in entries if isinstance(e, dict)),
        key=lambda e: e.get("frame", 0),
    )
    if len(ordered) != expected_frames:
        print(
            f"WARNING: physicality report covers {len(ordered)} frame(s) but the ensemble has "
            f"{expected_frames}. Frames beyond the report are reported as unknown.",
            flush=True,
        )
    return {
        "status": [str(e.get("status", "")) for e in ordered],
        "reason": ["" if e.get("reason") is None else str(e.get("reason")) for e in ordered],
    }


def parse_propagated_physicality(expected_frames):
    """Read node 04's physicality report for the propagated endpoints, if it ran.

    Propagation can produce an unphysical endpoint exactly as sampling can, and upstream's
    own multi-step entry point checks every frame at the end of a rollout for that reason.
    When node 04 propagated it writes this file; when it passed the ensemble through there is
    nothing new to check, so the file carries node 03's verdicts forward instead.
    """
    path = "./inputs/propagated_physicality.json"
    if not os.path.exists(path):
        print(
            f"NOTE: {path} is absent, so the propagated frames carry no physicality verdict "
            "and the report will say so rather than show them as unchecked passes.",
            flush=True,
        )
        return {"status": [], "reason": []}
    entries = read_json(path, "propagated physicality report", default=[])
    if not isinstance(entries, list):
        print(f"WARNING: {path} is not a JSON list; ignoring it.", flush=True)
        return {"status": [], "reason": []}
    ordered = sorted(
        (e for e in entries if isinstance(e, dict)), key=lambda e: e.get("frame", 0)
    )
    if len(ordered) != expected_frames:
        print(
            f"WARNING: propagated physicality report covers {len(ordered)} frame(s) but the "
            f"propagated ensemble has {expected_frames}.",
            flush=True,
        )
    return {
        "status": [str(e.get("status", "")) for e in ordered],
        "reason": ["" if e.get("reason") is None else str(e.get("reason")) for e in ordered],
    }


def select_frames(clusters, n_frames, top_n):
    """Representative frames: cluster centres by size first, then medoid-nearest.

    Takes the full cluster list, not the truncated copy that goes into the JSON: reading the
    truncated one silently capped the viewer at 20 distinct conformations and padded the rest
    with consecutive frames.
    """
    chosen = []
    if clusters:
        chosen = [c["centre"] for c in clusters][:top_n]
    from_clusters = len(chosen)
    if len(chosen) < top_n:
        for frame in range(n_frames):
            if frame not in chosen:
                chosen.append(frame)
            if len(chosen) >= top_n:
                break
    chosen = chosen[: min(top_n, n_frames)]
    return chosen, min(from_clusters, len(chosen))


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

    equilibrium = load_trajectory("./inputs/samples.xtc", topology, "sampled")
    print(f"Sampled ensemble: {equilibrium.n_frames} frame(s), {equilibrium.n_atoms} atoms",
          flush=True)

    physicality = parse_physicality(
        "./inputs/physicality.json", equilibrium.n_frames, sampling_manifest
    )

    rows, summary, medoid, matrix, superposed, clusters = analyze_ensemble(
        equilibrium, "sampled", physicality, threshold
    )

    summaries = {"sampled": summary}
    all_rows = list(rows)

    propagated_path = "./inputs/propagated_samples.xtc"
    # Three cases that must not collapse into one message. Node 04 always writes a manifest,
    # so an absent one is a plumbing failure, not a disabled propagation; and a manifest that
    # says propagation ran, with no XTC beside it, is a lost declared output. Reporting
    # either as "propagation was disabled" is how a real failure becomes an ordinary-looking
    # single-ensemble report.
    if not propagation_manifest:
        fail(
            "Node 04's propagation_manifest.json is missing or unreadable. Node 04 writes it "
            "on every path, including the identity pass-through, so this is a staging failure "
            "rather than a disabled propagation. Whether these frames are time-lagged cannot "
            "be determined, and presenting them either way would be a guess."
        )
    propagated_enabled = bool(propagation_manifest.get("enabled"))
    if propagated_enabled and not os.path.exists(propagated_path):
        fail(
            f"Node 04 recorded enabled=true but {propagated_path} is not here. The propagated "
            "ensemble is a declared output of node 04, so it was lost in staging. Analysing "
            "the sampled ensemble alone would silently drop the propagation that was paid for."
        )

    if propagated_enabled:
        propagated = load_trajectory(propagated_path, topology, "propagated")
        print(f"Propagated ensemble: {propagated.n_frames} frame(s)", flush=True)
        prop_physicality = parse_propagated_physicality(propagated.n_frames)
        prop_rows, prop_summary, _, _, _, _ = analyze_ensemble(
            propagated,
            "propagated",
            prop_physicality,
            threshold,
        )
        prop_summary["lag_time_ns"] = propagation_manifest.get("lag_time_ns")
        summaries["propagated"] = prop_summary
        all_rows.extend(prop_rows)
    else:
        reason = propagation_manifest.get("reason") or "no reason recorded"
        print(
            f"Propagation was disabled ({reason}) — analysing the sampled ensemble only. "
            "Node 04's output is an identity copy and is not treated as time-lagged.",
            flush=True,
        )

    with open("./outputs/ensemble_metrics.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(all_rows)

    # Multi-model PDB of representative frames, for the report's 3D viewer. The frames are
    # taken from the medoid-superposed copy: the viewer cycles through models in place, and
    # un-superposed frames differ by a rigid-body rotation that reads as the whole molecule
    # tumbling, which looks like conformational change and is not.
    selected, n_from_clusters = select_frames(clusters, equilibrium.n_frames, top_n)
    try:
        superposed[selected].save_pdb("./outputs/ensemble_top.pdb")
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
        # How many of those frames are genuine cluster centres. The rest are consecutive
        # frames used as filler, and the report must not label them as cluster centres.
        "viewer_frames_from_clusters": int(n_from_clusters),
        "viewer_frames_superposed_on_frame": int(medoid),
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
    print(
        f"Wrote {len(selected)} representative frame(s) to ensemble_top.pdb "
        f"({n_from_clusters} cluster centre(s), {len(selected) - n_from_clusters} filler "
        f"frame(s)), each superposed on frame {medoid}",
        flush=True,
    )


if __name__ == "__main__":
    main()
