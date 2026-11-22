#!/usr/bin/env python3
"""Build the inputs nodes 03-06 read, deterministically and without a GPU.

Node 03 needs a GPU and a model; nodes 04, 05 and 06 need only an XTC, a system and some
JSON. Those three are where the measurement bugs lived, so this script supplies their
inputs directly. That is the same argument workflow-037's `make_designs.py` makes, and it
buys the same thing: the downstream nodes are exercised exactly as a real run exercises
them, and *deterministically*, which a sampled ensemble never could.

Runs inside the image, because it needs bioemu2 to turn a `system.yaml` into the topology
the XTC is bound to -- the same `System.topology` that nodes 04 and 05 use, so the atom
order matches by construction rather than by hope.

Two ensemble shapes, and the difference between them is the point:

  spread     a base conformation plus per-frame Gaussian displacement, each frame then
             rotated by a random rigid rotation.
  identical  the *same* conformation in every frame, each still rotated randomly.

The rotation is deliberate in both. A frame that differs from another only by a rigid
rotation is not a different conformation, so:

  * RMSD and RMSF must be invariant to it. `identical` therefore has to trip node 05's
    collapse warning despite every frame having different coordinates -- an implementation
    that forgot to superpose would measure a large spread here and report a healthy
    ensemble.
  * `ensemble_top.pdb` must hold the superposed frames. Written un-superposed, the 3D
    viewer animates the whole molecule tumbling, which reads as conformational change.
    With the rotation present that difference is visible; without it the bug is invisible.

Usage:
    python3 -I make_ensembles.py <out-dir> <system.yaml> <n_frames> <shape> [physicality]

`physicality` is `passed` (default) or `not_checked`.
"""
import json
import os
import sys
import zlib


def rotation(rng):
    """A uniformly random rotation matrix, via QR of a Gaussian matrix."""
    import numpy as np

    matrix = rng.normal(size=(3, 3))
    q, r = np.linalg.qr(matrix)
    # Fix the signs so q is a rotation rather than a reflection; a reflection would mirror
    # the structure, which is a different molecule and not a different pose.
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def base_conformation(n_atoms):
    """A compact, deterministic starting shape.

    An ascending spiral with roughly bond-length spacing. It is not a protein, and it does
    not need to be: every quantity node 05 computes -- radius of gyration, pairwise RMSD,
    RMSF about the mean -- is a function of coordinates alone. What it does need to be is
    compact, finite, and the same on every machine.
    """
    import numpy as np

    index = np.arange(n_atoms, dtype=float)
    turns = index * 0.6
    radius = 0.9 + 0.35 * np.sqrt(index / max(n_atoms, 1))
    return np.stack(
        [
            radius * np.cos(turns),
            radius * np.sin(turns),
            index * 0.015 - 0.0075 * n_atoms,
        ],
        axis=1,
    )


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        return 2
    out_dir, system_path, n_frames, shape = sys.argv[1:5]
    physicality_status = sys.argv[5] if len(sys.argv) > 5 else "passed"
    # The mode node 03 would have recorded. Node 05 reads it rather than inferring anything
    # from the frame count, so a fixture that wants to exercise a mode-dependent message has
    # to declare the mode here.
    mode = sys.argv[6] if len(sys.argv) > 6 else "equilibrium"
    n_frames = int(n_frames)
    if shape not in ("spread", "identical"):
        print(f"ERROR: unknown shape '{shape}'")
        return 2
    if mode not in ("equilibrium", "structure"):
        print(f"ERROR: unknown mode '{mode}'")
        return 2

    import mdtraj
    import numpy as np
    from pathlib import Path

    from bioemu2 import System

    os.makedirs(out_dir, exist_ok=True)
    system = System.from_yaml(Path(system_path))
    topology = system.topology
    n_atoms = topology.n_atoms

    # Seeded per shape, size and mode, so a rerun of the same case is byte-identical and a
    # different case is genuinely different.
    #
    # Not `hash()`. Python salts the hash of a str per process unless PYTHONHASHSEED is set,
    # so this used to draw a different ensemble on every run while claiming to be
    # deterministic -- and that hid a flaky check: `node05_medoid_measured_over_all_pairs`
    # asserted that `rmsd_to_medoid_nm` is zero at the medoid, which is only true to within
    # mdtraj's float32 RMSD noise (measured 0.0013-0.0018 nm for a frame against itself), so
    # it passed or failed depending on which frame the dice made the medoid.
    key = f"{shape}-{n_frames}-{n_atoms}-{mode}".encode()
    rng = np.random.default_rng(zlib.crc32(key))

    base = base_conformation(n_atoms)
    frames = np.empty((n_frames, n_atoms, 3), dtype=np.float32)
    for index in range(n_frames):
        conformation = base if shape == "identical" else base + rng.normal(
            scale=0.05, size=base.shape
        )
        centred = conformation - conformation.mean(axis=0)
        frames[index] = centred @ rotation(rng).T

    trajectory = mdtraj.Trajectory(xyz=frames, topology=topology)
    trajectory.save_xtc(os.path.join(out_dir, "samples.xtc"))
    trajectory[0].save(os.path.join(out_dir, "topology.cif"))

    # Exactly the shape PhysicalityReport.save() writes: a bare list, one record per frame.
    records = [
        {
            "frame": index,
            "status": physicality_status,
            "reason": None,
        }
        for index in range(n_frames)
    ]
    with open(os.path.join(out_dir, "physicality.json"), "w") as handle:
        json.dump(records, handle, indent=2)

    checked = physicality_status in ("passed", "rejected")
    n_passed = sum(1 for r in records if r["status"] == "passed") if checked else None
    manifest = {
        "settings": {
            "mode": mode,
            "num_samples": n_frames,
            "batch_size": None,
            "seed": 42,
            "sampler": "edm",
            "num_steps": None,
            "physicality_filter": "high_quality" if checked else "no_filter",
            "precision": "bf16-mixed",
            "matmul_precision": "high",
            "compile_denoiser": False,
            "checkpoint": "bioemu2-preview",
            "hf_ckpt_revision": "fixture",
        },
        "bioemu2_version": "fixture",
        "n_frames": n_frames,
        "n_atoms": n_atoms,
        "physicality_checked": checked,
        "physicality_not_checked_because": (
            None if checked else "physicality_filter=no_filter, so no thresholds were applied"
        ),
        "n_passed": n_passed,
        "pass_rate": (round(n_passed / n_frames, 3) if n_passed is not None else None),
        "rejection_reasons": {},
        "msa_source": "fixture",
        "msa_min_depth": 64,
        "msa_loaded": True,
        "msa_how": "fixture: supplied by tests/make_ensembles.py",
        "msa_files": [],
        "msa_paired": False,
        "gpu": "none (fixture)",
        "elapsed_seconds": 0.0,
        "frames_are_independent_samples": True,
        "caveat": "Fixture ensemble. These coordinates are synthetic, not sampled.",
    }
    with open(os.path.join(out_dir, "sampling_manifest.json"), "w") as handle:
        json.dump(manifest, handle, indent=2)

    print(
        f"wrote {n_frames} frame(s) x {n_atoms} atoms, shape={shape}, "
        f"physicality={physicality_status} -> {out_dir}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
