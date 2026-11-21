#!/usr/bin/env python3
"""Synthesise what node 03 produces, with the answer known before node 04 runs.

Node 03 needs a GPU and a minute per design, which would put nodes 04 and 05 out of reach
of a suite anyone actually runs. But those nodes read only a scores JSON, a rank_001 PDB
and a manifest, so a hand-built set exercises them exactly as a real prediction does --
and with a true answer that is established by construction rather than by trusting this
file.

That last part is the whole point. A fixture whose expected value comes from running the
code under test asserts nothing. So:

  identity      the designed backbone, copied. True sc_rmsd_A = 0, BY CONSTRUCTION.
  rotated       the designed backbone under a rigid rotation and a large translation.
                True sc_rmsd_A = 0, BY CONSTRUCTION, because RMSD after superposition is
                invariant under rigid motion. This is the keystone case: a node that
                computes RMSD without superposing returns something enormous here and so
                rejects every good design while reporting plausible-looking failures --
                which is exactly how 029's chain-order defect presented (d5ba7487).
  displaced:<d> the designed backbone plus a displacement field of root-mean-square
                magnitude exactly d, constructed to be orthogonal to every rigid motion,
                then rigidly moved. True sc_rmsd_A = d.

Why `displaced` needs the orthogonality, and how it is guaranteed
----------------------------------------------------------------
If the displacement has any component along a translation or a rotation, the optimal
superposition absorbs part of it and the true RMSD is somewhere below d -- by an amount
that depends on the geometry, so it cannot be written down in a fixture.

With P = Q + D and both centred, the cross-covariance is

    H = P^T Q = Q^T Q + D^T Q

and the optimal rotation is the identity exactly when H is symmetric positive definite.
Q^T Q always is. The antisymmetric part of D^T Q vanishes iff sum_i (Q_i x D_i) = 0, which
is precisely the statement that D is orthogonal to every infinitesimal rotation. So the
construction projects both components out of an arbitrary starting field -- the mean, then
the least-squares rotation generator -- and only then scales what is left to
root-mean-square d.

That is an argument, so it is also checked. This file runs its OWN Kabsch, which node 04
never sees, and refuses to write anything unless the measured RMSD equals d to 1e-9 and
the recovered rotation is the identity to 1e-9. If a geometry ever makes H non-positive --
d large against the structure's own extent -- the check fails loudly here instead of
quietly handing node 04 an expected value that is wrong.

Chain order, pLDDT and PAE
--------------------------
The scores file carries NO chain boundaries -- measured, in tests/reference/: a flat
plddt list and a square pae, with nothing naming a chain. Node 04 therefore slices by
position, in the order node 02 recorded. So this file can write the binder block first or
last, with different values in each block, and a node that hardcodes 029's
plddt[:binder_length] passes one case and fails the other.

PAE blocks are set independently per direction, because the interface figure is the mean
over BOTH off-diagonal blocks. Setting target->binder to 6 and binder->target to 14 makes
the right answer exactly 10 and makes averaging one block alone give 6 or 14.

Run
---
  python3 tests/make_predictions.py <inventory> <prediction-inputs> <out-dir> [options]

    --shape S            identity | rotated | displaced:<d>        (default identity)
    --chain-order O      comma-separated, from {binder,target}     (default binder)
    --plddt-binder V     fill value for the binder block           (default 90.0)
    --plddt-target V     fill value for the target block           (default 85.0)
    --iptm V             iptm value to write
    --omit KEYS          comma-separated keys to leave OUT: iptm,pae,plddt,max_pae,ptm
    --pae-target-binder V   PAE in the target->binder block
    --pae-binder-target V   PAE in the binder->target block
    --pae-diagonal V     PAE inside each chain's own block         (default 1.0)
    --length-delta N     write N extra (or fewer) binder residues than the design has
    --plddt-delta N      write a plddt array N entries longer than the chain lengths
    --chains N           force the PDB to carry N chains
    --no-score-file D    write no scores file for design D
    --fail-design D      mark design D as failed in the manifest, with no output files
    --model-type T       recorded in filenames and the manifest    (default alphafold2_ptm)
"""
import json
import math
import os
import sys

ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    "X": "UNK", "U": "SEC", "O": "PYL",
}

# A rotation that mixes all three axes and a translation far enough out that an
# unsuperposed RMSD is hundreds of Angstroms rather than a few. 0.7 rad about an axis
# that is not any of the basis vectors: nothing here can be accidentally the identity.
ROTATION_AXIS = (0.4679, -0.6245, 0.6255)   # unit to 1e-6
ROTATION_ANGLE_RAD = 0.7
TRANSLATION_A = (137.0, -211.0, 59.0)

# How exactly the construction has to hold in floating point, before anything is written.
CONSTRUCTION_TOLERANCE = 1e-9

# And how exactly it can still hold after being written. A PDB ATOM record carries three
# decimal places, so every coordinate is rounded to 0.001 A on the way out: about 0.0005 A
# per axis, 0.00087 A in three dimensions, which moves the RMSD by roughly 1e-4. Measured
# on this suite's own 30-residue case: constructed 2.000000000000, as written
# 1.999883143312.
#
# That is not an inaccuracy to be hidden in a loose tolerance. It is the precision the
# FILE can carry, and node 04 reads the file, so truth.json records both numbers and the
# check below is on the artifact rather than on the array in memory (§8.3).
PDB_ROUNDING_TOLERANCE_A = 2e-3


def die(message):
    print(f"make_predictions: ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def rotation_matrix():
    import numpy as np

    axis = np.array(ROTATION_AXIS, dtype=float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c, s = math.cos(ROTATION_ANGLE_RAD), math.sin(ROTATION_ANGLE_RAD)
    return np.array([
        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
    ])


def kabsch(mobile, reference):
    """(rmsd, rotation) after optimal superposition. This file's own, not node 04's.

    Kept separate on purpose: a construction verified with the code it is meant to test
    verifies nothing (§8.4). If this and node 04 ever disagree, one of them is wrong and
    the suite says which value each produced.
    """
    import numpy as np

    P = np.asarray(mobile, dtype=float)
    Q = np.asarray(reference, dtype=float)
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    aligned = (R @ Pc.T).T
    rmsd = float(np.sqrt(((aligned - Qc) ** 2).sum() / len(P)))
    return rmsd, R


def rigid_free_displacement(coords, magnitude):
    """A displacement field of root-mean-square `magnitude`, orthogonal to rigid motion.

    Starts from an alternating field along a fixed axis -- deterministic, and nothing like
    a rigid motion to begin with -- then removes the translation component (its mean) and
    the rotation component (the least-squares generator omega minimising
    ||D - omega x Q||), and scales the remainder.
    """
    import numpy as np

    Q = np.asarray(coords, dtype=float)
    Q = Q - Q.mean(axis=0)
    n = len(Q)
    if n < 4:
        die(f"a displaced shape needs at least 4 residues, got {n}")

    signs = np.array([1.0 if i % 2 == 0 else -1.0 for i in range(n)])
    # Two axes, so a chain that happens to lie along one of them still gets a field with
    # a component across it.
    base = np.outer(signs, np.array([1.0, 0.0, 0.0]))
    base += np.outer(np.roll(signs, 1), np.array([0.0, 1.0, 0.0]))

    D = base - base.mean(axis=0)

    # omega x Q is linear in omega: build the 3n x 3 operator and solve least squares.
    operator = np.zeros((3 * n, 3))
    for i, q in enumerate(Q):
        # cross(omega, q) = -cross(q, omega) = -[q]_x omega
        qx = np.array([
            [0.0, -q[2], q[1]],
            [q[2], 0.0, -q[0]],
            [-q[1], q[0], 0.0],
        ])
        operator[3 * i:3 * i + 3, :] = -qx
    omega, *_ = np.linalg.lstsq(operator, D.reshape(-1), rcond=None)
    D = D - np.cross(np.broadcast_to(omega, Q.shape), Q)
    D = D - D.mean(axis=0)

    scale = math.sqrt(float((D ** 2).sum() / n))
    if scale < 1e-12:
        die("the displacement field collapsed to zero after removing rigid motion; "
            "this geometry cannot carry a known displacement")
    return D * (magnitude / scale)


def build_coords(designed_ca, shape):
    """Return (coords, true_rmsd). Verified here, with this file's own Kabsch."""
    import numpy as np

    Q = np.asarray([c for c in designed_ca if c is not None], dtype=float)
    if len(Q) < 1:
        die("the design has no CA coordinates")

    if shape == "identity":
        P, truth = Q.copy(), 0.0
    elif shape == "rotated":
        P = (rotation_matrix() @ Q.T).T + np.array(TRANSLATION_A)
        truth = 0.0
    elif shape.startswith("displaced:"):
        magnitude = float(shape.split(":", 1)[1])
        if magnitude <= 0:
            die(f"displaced magnitude must be positive, got {magnitude}")
        D = rigid_free_displacement(Q, magnitude)
        Qc = Q - Q.mean(axis=0)
        # Check identity is really optimal BEFORE the rigid move, where the claim lives.
        measured, R = kabsch(Qc + D, Qc)
        off_identity = float(np.abs(R - np.eye(3)).max())
        if off_identity > CONSTRUCTION_TOLERANCE:
            die(
                f"the displacement field is not orthogonal to rigid motion: the optimal "
                f"rotation differs from the identity by {off_identity:.3e}. The true RMSD "
                "would then be below the requested magnitude by an amount that depends on "
                "the geometry, so no fixture could state it."
            )
        if abs(measured - magnitude) > CONSTRUCTION_TOLERANCE:
            die(
                f"constructed RMSD is {measured!r}, not the requested {magnitude!r} "
                f"(difference {abs(measured - magnitude):.3e}). Refusing to write an "
                "expected value this file cannot stand behind."
            )
        P = (rotation_matrix() @ (Qc + D).T).T + np.array(TRANSLATION_A)
        truth = magnitude
    else:
        die(f"unknown shape '{shape}'")

    # The same check again on what is actually written, after the rigid move, because that
    # is the array node 04 will read.
    final, _ = kabsch(P, Q)
    if abs(final - truth) > CONSTRUCTION_TOLERANCE:
        die(
            f"the coordinates about to be written superpose onto the design at RMSD "
            f"{final!r}, not the claimed {truth!r}."
        )
    return [list(map(float, row)) for row in P], truth


def read_pdb_ca(path):
    """CA coordinates per chain, in file order. Reads back what was just written."""
    chains, order = {}, []
    with open(path) as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            chain = line[21]
            if chain not in chains:
                chains[chain] = []
                order.append(chain)
            chains[chain].append(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
    return [(cid, chains[cid]) for cid in order]


def write_pdb(path, blocks):
    """blocks: [(chain_id, sequence, coords), ...]. CA only, renumbered from 1 per chain.

    Renumbered from 1 per chain because that is what ColabFold does, and it is the trap:
    the design's own numbering may start anywhere, so nothing downstream may match on
    residue number (ef74646e).
    """
    serial = 1
    with open(path, "w") as handle:
        for chain_id, sequence, coords in blocks:
            for i, (letter, xyz) in enumerate(zip(sequence, coords), start=1):
                name = ONE_TO_THREE.get(letter, "UNK")
                handle.write(
                    "ATOM  %5d  CA  %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 50.00           C\n"
                    % (serial, name, chain_id, i, xyz[0], xyz[1], xyz[2])
                )
                serial += 1
            handle.write("TER\n")
        handle.write("END\n")


def parse_args(argv):
    options = {
        "shape": "identity",
        "chain_order": "binder",
        "plddt_binder": 90.0,
        "plddt_target": 85.0,
        "iptm": None,
        "omit": "",
        "pae_target_binder": None,
        "pae_binder_target": None,
        "pae_diagonal": 1.0,
        "length_delta": 0,
        "plddt_delta": 0,
        "chains": None,
        "no_score_file": None,
        "fail_design": None,
        "model_type": "alphafold2_ptm",
    }
    floats = {"plddt_binder", "plddt_target", "iptm", "pae_target_binder",
              "pae_binder_target", "pae_diagonal"}
    ints = {"length_delta", "plddt_delta", "chains"}
    i = 0
    while i < len(argv):
        token = argv[i]
        if not token.startswith("--"):
            die(f"unexpected argument '{token}'")
        key = token[2:].replace("-", "_")
        if key not in options:
            die(f"unknown option '{token}'")
        if i + 1 >= len(argv):
            die(f"'{token}' needs a value")
        raw = argv[i + 1]
        options[key] = float(raw) if key in floats else int(raw) if key in ints else raw
        i += 2
    return options


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    inventory_path, inputs_path, out_dir = sys.argv[1:4]
    options = parse_args(sys.argv[4:])
    omit = {part.strip() for part in options["omit"].split(",") if part.strip()}
    order = [part.strip() for part in options["chain_order"].split(",") if part.strip()]
    for name in order:
        if name not in ("binder", "target"):
            die(f"chain_order may only contain 'binder' and 'target', got '{name}'")
    if "binder" not in order:
        die("chain_order must contain 'binder'")

    with open(inventory_path) as handle:
        inventory = json.load(handle)
    with open(inputs_path) as handle:
        declared = json.load(handle)

    os.makedirs(out_dir, exist_ok=True)
    by_name = {e["design"]: e for e in inventory["designs"]}

    predictions, truth_rows = [], []
    for record in declared.get("inputs") or []:
        design = record["design"]
        entry = by_name.get(design)
        if entry is None:
            die(f"'{design}' is declared as an input but is not in the inventory")

        base = {
            "design": design,
            "mode": record["mode"],
            "record_id": record["record_id"],
            "model_type": options["model_type"],
            "msa_mode": "single_sequence",
            "num_models": 1,
            "num_recycle": 3,
            "host_url_passed": False,
            "argv": ["make_predictions.py", "--shape", options["shape"]],
            "elapsed_s": 0.0,
            "log": None,
        }

        if design == options["fail_design"]:
            predictions.append({
                **base,
                "chain_order": record["chain_order"],
                "chain_lengths": record["chain_lengths"],
                "expected_chains": record["expected_chains"],
                "expected_residues": record["expected_residues"],
                "status": "failed",
                "reason": (
                    "synthetic failure, so the suite can ask whether a design with no "
                    "prediction reaches the report as a failed row rather than vanishing "
                    "from it (§8.5)"
                ),
                "outputs": [],
            })
            truth_rows.append({"design": design, "status": "failed"})
            continue

        binder_seq = entry["binder_sequence"]
        binder_coords, truth = build_coords(entry["binder_ca"], options["shape"])

        if options["length_delta"]:
            delta = options["length_delta"]
            if delta > 0:
                # Repeat the last residue: the point is the LENGTH, and a mismatch has to
                # be refused rather than truncated and scored (§8.1).
                binder_seq = binder_seq + binder_seq[-1] * delta
                binder_coords = binder_coords + [binder_coords[-1]] * delta
            else:
                binder_seq = binder_seq[:delta]
                binder_coords = binder_coords[:delta]
            truth = None

        target_seq = "".join(t["sequence"] for t in entry["target_chains"])
        # Placed far from the binder: these coordinates are never scored against anything,
        # and overlapping them would make a reader think they were meant to be.
        target_coords = [[500.0 + 3.8 * i, 0.0, 0.0] for i in range(len(target_seq))]

        blocks, lengths, written_order = [], [], []
        letters = iter("ABCDEFGH")
        for name in order:
            if name == "binder":
                blocks.append((next(letters), binder_seq, binder_coords))
                lengths.append(len(binder_seq))
            else:
                if not target_seq:
                    die(f"'{design}' has no target chain, so chain_order cannot "
                        "include 'target'")
                blocks.append((next(letters), target_seq, target_coords))
                lengths.append(len(target_seq))
            written_order.append(name)

        if options["chains"] is not None and options["chains"] != len(blocks):
            if options["chains"] == 1 and len(blocks) == 2:
                # One chain carrying both sequences: exactly what a cplx_ input with its
                # ':' stripped comes back as, which is 026's failure.
                merged = "".join(seq for _, seq, _ in blocks)
                merged_coords = [c for _, _, coords in blocks for c in coords]
                blocks = [("A", merged, merged_coords)]
            else:
                die(f"--chains {options['chains']} is not constructible from "
                    f"{len(blocks)} block(s)")

        suffix = f"rank_001_{options['model_type']}_model_1_seed_000"
        pdb_name = f"{design}_unrelaxed_{suffix}.pdb"
        pdb_path = os.path.join(out_dir, pdb_name)
        write_pdb(pdb_path, blocks)

        # Measure the FILE, not the array that went into it. The rounding to three decimal
        # places is real, node 04 reads the rounded numbers, and a check that skipped this
        # step would be asserting the construction rather than the fixture (§8.3).
        truth_as_written = None
        if truth is not None and len(blocks) == len(written_order):
            index = written_order.index("binder")
            _, written_binder = read_pdb_ca(pdb_path)[index]
            reference = [c for c in entry["binder_ca"] if c is not None]
            if len(written_binder) == len(reference):
                truth_as_written, _ = kabsch(written_binder, reference)
                if abs(truth_as_written - truth) > PDB_ROUNDING_TOLERANCE_A:
                    die(
                        f"'{pdb_name}' as written superposes onto the design at RMSD "
                        f"{truth_as_written!r}, which is further than "
                        f"{PDB_ROUNDING_TOLERANCE_A} from the constructed {truth!r}. PDB "
                        "coordinates carry three decimal places and nothing else should "
                        "differ, so this is not rounding."
                    )

        plddt = []
        for name, length in zip(written_order, lengths):
            value = (options["plddt_binder"] if name == "binder"
                     else options["plddt_target"])
            plddt.extend([value] * length)
        if options["plddt_delta"]:
            delta = options["plddt_delta"]
            plddt = (plddt + [50.0] * delta) if delta > 0 else plddt[:delta]

        total = sum(lengths)
        pae = [[options["pae_diagonal"]] * total for _ in range(total)]
        if len(written_order) == 2:
            first_len = lengths[0]
            for row in range(total):
                for col in range(total):
                    first_row = row < first_len
                    first_col = col < first_len
                    if first_row == first_col:
                        continue
                    source = written_order[0] if first_row else written_order[1]
                    key = ("pae_target_binder" if source == "target"
                           else "pae_binder_target")
                    value = options[key]
                    if value is not None:
                        pae[row][col] = value

        scores = {}
        if "plddt" not in omit:
            scores["plddt"] = plddt
        if "pae" not in omit:
            scores["pae"] = pae
        if "max_pae" not in omit:
            scores["max_pae"] = max((max(row) for row in pae), default=0.0)
        if "ptm" not in omit:
            scores["ptm"] = 0.75
        # Absent unless asked for, because that is the MEASURED monomer case: the real
        # file carries ptm and no iptm key at all (tests/reference/).
        if options["iptm"] is not None and "iptm" not in omit:
            scores["iptm"] = options["iptm"]

        outputs = [pdb_name]
        scores_name = f"{design}_scores_{suffix}.json"
        if design != options["no_score_file"]:
            with open(os.path.join(out_dir, scores_name), "w") as handle:
                json.dump(scores, handle, indent=2)
            outputs.append(scores_name)

        failed = design == options["no_score_file"]
        predictions.append({
            **base,
            "chain_order": written_order,
            "chain_lengths": lengths,
            # What node 02 DECLARED, which is the number of entries in chain_order -- not
            # what the PDB happens to carry. With --chains 1 the two differ on purpose,
            # and the disagreement is the fixture: a complex that came back one chain
            # means the ':' was dropped, and node 04 has to notice without being told.
            "expected_chains": len(written_order),
            "expected_residues": total,
            "returned_chains": {cid: len(seq) for cid, seq, _ in blocks},
            "status": "failed" if failed else "ok",
            "reason": (
                "synthetic: colabfold_batch exited 0 and produced no rank_001 scores "
                "file. Exit codes see nothing (§8.3)." if failed else None
            ),
            "outputs": sorted(outputs),
            **({} if failed else {
                "scores_rank_001": scores_name,
                "pdb_rank_001": pdb_name,
            }),
        })
        truth_rows.append({
            "design": design,
            "status": "failed" if failed else "ok",
            "sc_rmsd_A": truth,
            # The same number measured off the written file. These differ by about 1e-4
            # because a PDB carries three decimal places; node 04 reads the file, so a
            # check that wants better than 1e-3 should compare against this one.
            "sc_rmsd_A_as_written": truth_as_written,
            "plddt_binder": options["plddt_binder"],
            "plddt_target": (options["plddt_target"]
                             if "target" in written_order else None),
            "iptm": options["iptm"] if "iptm" not in omit else None,
            "pae_interaction": (
                None if len(written_order) != 2
                or options["pae_target_binder"] is None
                or options["pae_binder_target"] is None
                or "pae" in omit
                else (options["pae_target_binder"] + options["pae_binder_target"]) / 2.0
            ),
            "chain_order": written_order,
            "chain_lengths": lengths,
        })

    manifest = {
        "n_declared": len(predictions),
        "n_ok": len([p for p in predictions if p["status"] == "ok"]),
        "n_failed": len([p for p in predictions if p["status"] == "failed"]),
        "dry_run": False,
        "phase": 2 if any("target" in p.get("chain_order", []) for p in predictions) else 1,
        "model_type": options["model_type"],
        "msa_mode": "single_sequence",
        "num_models": 1,
        "num_recycle": 3,
        "predict_timeout_s": 1800,
        "host_url_passed": False,
        "network_expected": "none",
        "weights": "synthetic: no predictor ran",
        "chain_order_source": "02-build-prediction-inputs",
        "designed_chain_sources": inventory.get("designed_chain_sources"),
        "securedna_screened": False,
        "synthetic": True,
        "generator": "tests/make_predictions.py",
        "predictions": predictions,
    }
    with open(os.path.join(out_dir, "prediction_manifest.json"), "w") as handle:
        json.dump(manifest, handle, indent=2)
    with open(os.path.join(out_dir, "truth.json"), "w") as handle:
        json.dump({
            "shape": options["shape"],
            "options": options,
            "note": (
                "sc_rmsd_A here is established by construction and verified by this "
                "file's own Kabsch to 1e-9, which node 04 never sees. null means the "
                "case is about a refusal rather than a value."
            ),
            "designs": truth_rows,
        }, handle, indent=2)

    print(
        f"make_predictions: shape={options['shape']} order={','.join(order)} "
        f"-> {len(predictions)} prediction(s) in {out_dir}; "
        f"true sc_rmsd_A={truth_rows[0].get('sc_rmsd_A') if truth_rows else 'n/a'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
