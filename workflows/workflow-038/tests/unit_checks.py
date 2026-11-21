#!/usr/bin/env python3
"""The three checks that need no input file, only the node code and mdtraj.

Run inside the image, because `pairwise_pass` is mdtraj all the way down:

    python -I tests/unit_checks.py <out.json>

Writes {check_name: {"status": ..., "detail": ...}} for check_fixtures.py to report, so a
failure here reads in the same table as everything else rather than as a traceback.
"""
import importlib.util
import json
import os
import sys
import tempfile

WF = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(node):
    """A node's script as a module. They import mdtraj and bioemu2 inside functions, so
    importing one costs nothing and runs nothing: everything is behind `def main()`."""
    path = os.path.join(WF, node, "script.py")
    spec = importlib.util.spec_from_file_location(node.replace("-", "_"), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def unit_a3m_path_helpers():
    """The unpaired set must exclude the paired files, and the paired map must be keyed on
    the chain id.

    `*.a3m` also matches `<chain>.paired.a3m`. bioemu2's `load_a3m_msa` requires the
    mapping's keys to equal the system's protein chain ids exactly, so a leaked
    `A.paired.a3m` raises -- and keying the paired map on `Path.stem` gives `A.paired`,
    which is not a chain id either. Both mistakes are one character from correct.
    """
    sample = load("03-sample")
    fetch = load("02-fetch-msa")
    with tempfile.TemporaryDirectory() as directory:
        for name in ("A.a3m", "B.a3m", "A.paired.a3m", "B.paired.a3m"):
            with open(os.path.join(directory, name), "w") as handle:
                handle.write(">query\nAAA\n")

        unpaired = [os.path.basename(p) for p in sample.unpaired_a3m_paths(directory)]
        if unpaired != ["A.a3m", "B.a3m"]:
            return "fail", f"node 03's unpaired set is {unpaired}, not ['A.a3m', 'B.a3m']"

        node02 = [os.path.basename(p) for p in fetch.unpaired_a3m_paths(directory)]
        if node02 != unpaired:
            return "fail", f"node 02 sees {node02} where node 03 sees {unpaired}"

        paired = sample.paired_a3m_paths(directory)
        if sorted(paired) != ["A", "B"]:
            return "fail", (
                f"paired map is keyed {sorted(paired)}, not ['A', 'B']: these keys go "
                "straight to load_a3m_msa, which matches them against chain ids"
            )
        for chain, path in paired.items():
            if os.path.basename(str(path)) != f"{chain}.paired.a3m":
                return "fail", f"chain {chain} maps to {path}"
    return "pass", "unpaired ['A.a3m', 'B.a3m'], paired keyed ['A', 'B'], nodes 02 and 03 agree"


def _toy_trajectory(n_frames, n_atoms, seed):
    import mdtraj
    import numpy as np

    topology = mdtraj.Topology()
    chain = topology.add_chain()
    for index in range(n_atoms):
        residue = topology.add_residue("ALA", chain, resSeq=index + 1)
        topology.add_atom("CA", mdtraj.core.element.carbon, residue)
    rng = np.random.default_rng(seed)
    base = rng.normal(scale=1.0, size=(n_atoms, 3))
    frames = np.stack([base + rng.normal(scale=0.3, size=base.shape) for _ in range(n_frames)])
    return mdtraj.Trajectory(frames.astype("float32"), topology)


def unit_pairwise_pass():
    """The streaming pass must agree with the dense matrix, matrix retained or not.

    Node 05 dropped the matrix above 200 frames and set `medoid = 0` while keeping the
    column name `rmsd_to_medoid_nm`: the ensemble's own reference frame became "whichever
    frame was sampled first", under a label that said otherwise. The streaming pass exists
    so the medoid and the spread survive at any frame count, which is only true if it
    gives the same answer as the matrix on a case where both are computable.
    """
    import mdtraj
    import numpy as np

    analyze = load("05-analyze")
    traj = _toy_trajectory(n_frames=9, n_atoms=12, seed=7)
    indices = np.arange(traj.n_atoms)

    dense = np.stack([mdtraj.rmsd(traj, traj, frame=i, atom_indices=indices) for i in range(9)])
    upper = dense[np.triu_indices(9, k=1)]
    reference_medoid = int(np.argmin(dense.sum(axis=1)))

    for keep in (True, False):
        sums, stats, matrix = analyze.pairwise_pass(traj, indices, keep)
        medoid = int(np.argmin(sums))
        if medoid != reference_medoid:
            return "fail", (
                f"keep_matrix={keep}: medoid {medoid} against {reference_medoid} from the "
                "dense matrix"
            )
        if (matrix is None) == keep:
            return "fail", f"keep_matrix={keep} returned matrix={type(matrix).__name__}"
        for field, expected in (
            ("max_nm", upper.max()),
            ("min_nm", upper.min()),
            ("mean_nm", upper.mean()),
        ):
            if abs(stats[field] - expected) > 1e-3:
                return "fail", (
                    f"keep_matrix={keep}: {field}={stats[field]} against {expected:.4f} "
                    "measured over the dense upper triangle"
                )
        if stats["n_pairs"] != 36:
            return "fail", f"n_pairs={stats['n_pairs']} for 9 frames, expected 36"
    return "pass", (
        f"9 frames, 36 pairs: medoid {reference_medoid}, max {upper.max():.4f} nm, "
        "identical with and without the matrix"
    )


def unit_select_frames():
    """The viewer's frames, and how many of them are filler.

    The viewer takes cluster centres first and pads with consecutive frames. Reporting the
    pad as though it were a distinct conformation overstates what the ensemble found, and
    reading the truncated cluster list capped the viewer at 20 regardless of `top_n`.
    """
    analyze = load("05-analyze")
    clusters = [{"centre": 4, "size": 10}, {"centre": 1, "size": 3}]

    chosen, from_clusters = analyze.select_frames(clusters, n_frames=20, top_n=5)
    if chosen[:2] != [4, 1] or from_clusters != 2 or len(chosen) != 5:
        return "fail", f"two clusters, top_n=5 gave {chosen} with {from_clusters} from clusters"
    if len(set(chosen)) != len(chosen):
        return "fail", f"filler repeated a cluster centre: {chosen}"

    chosen, from_clusters = analyze.select_frames(clusters, n_frames=20, top_n=1)
    if chosen != [4] or from_clusters != 1:
        return "fail", f"top_n=1 gave {chosen} with {from_clusters} from clusters"

    chosen, from_clusters = analyze.select_frames(clusters, n_frames=3, top_n=8)
    if len(chosen) != 3:
        return "fail", f"3 frames, top_n=8 gave {len(chosen)} frames: {chosen}"
    if from_clusters > len(chosen):
        return "fail", f"{from_clusters} reported from clusters out of {len(chosen)} chosen"

    chosen, from_clusters = analyze.select_frames([], n_frames=6, top_n=3)
    if chosen != [0, 1, 2] or from_clusters != 0:
        return "fail", f"no clusters gave {chosen} with {from_clusters} reported from clusters"
    return "pass", (
        "2 centres then filler (5 frames, 2 real), capped at top_n, capped at n_frames, "
        "and 0 claimed when there are no clusters"
    )


CHECKS = {
    "unit_a3m_path_helpers": unit_a3m_path_helpers,
    "unit_pairwise_pass": unit_pairwise_pass,
    "unit_select_frames": unit_select_frames,
}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "unit_results.json"
    results = {}
    for name, runner in CHECKS.items():
        try:
            status, detail = runner()
        except Exception as exc:
            status, detail = "fail", f"raised {type(exc).__name__}: {exc}"
        results[name] = {"status": status, "detail": detail}
        print(f"  {status.upper():<4} {name}  {detail}")
    with open(out, "w") as handle:
        json.dump(results, handle, indent=2)
    return 1 if any(r["status"] == "fail" for r in results.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
