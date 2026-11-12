#!/usr/bin/env python3
"""Node 02 — Model Fv Structures.

One Fv structure per annotated clone: ABodyBuilder2 for paired VH/VL, NanoBodyBuilder2 for VHH.

Why the Python API and not the `ABodyBuilder2` CLI: the CLI builds one antibody per process,
which would reload four sets of model weights for every clone in the panel. The predictor is
constructed once here and reused across the whole panel, which is the difference between
seconds and minutes per clone.

Refinement is not a flag upstream. `Antibody.save()` always calls `refine()` and then
`add_errors_as_bfactors()`; `save_single_unrefined()` skips both. So `refine_structures=false`
takes the unrefined path and this node adds the B-factors itself, so that both paths produce
files whose B-factor column means the same thing: ImmuneBuilder's per-residue error estimate
(the spread of its four-model ensemble), not a crystallographic B-factor.

The predictor numbers sequences internally with ANARCI in `numbering_scheme`, the same engine
node 01 used, so PDB residue ids line up with node 01's annotation and node 03 can map
structure residues to CDRs by (chain, residue id, insertion code).

Supports --dry-run: resolve and print every setting, then exit 0 without loading weights.

Per-clone failures are recorded, not fatal. A clone that fails to model still has humanness
from node 04, which depends on node 01 rather than on this node for exactly that reason.
"""
import csv
import json
import os
import sys
import time

# Nothing here should reach the network: the weights are baked into the image at build time.
os.environ.setdefault("HF_HUB_OFFLINE", "1")


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def warn(message):
    print(f"WARNING: {message}", flush=True)


def package_version(name):
    try:
        from importlib.metadata import version

        return version(name)
    except Exception:
        return "unknown"


def read_json(path, label, required=True):
    if not os.path.exists(path):
        if required:
            fail(f"{label} not found: {path}. Node 01 must run first.")
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        fail(f"Could not read {label} ('{path}'): {exc}")


def per_residue_errors(prediction):
    """ImmuneBuilder's per-residue error estimate, in the same units save() writes.

    `error_estimates` is the squared deviation of the four-model ensemble about its mean, per
    model and residue. Upstream's own B-factor call is `error_estimates.mean(0).sqrt()`, so this
    mirrors it exactly rather than inventing a summary.
    """
    try:
        return prediction.error_estimates.mean(0).sqrt().cpu().numpy()
    except Exception as exc:
        warn(f"Could not extract per-residue error estimates: {exc}")
        return None


def save_prediction(prediction, path, refine, clone_id):
    """Write the top-ranked model, with the error estimate in the B-factor column either way."""
    from ImmuneBuilder.util import add_errors_as_bfactors

    errors = per_residue_errors(prediction)

    if refine:
        # save() refines and adds B-factors itself.
        prediction.save(path)
        return "refined"

    best = prediction.ranking[0]
    prediction.save_single_unrefined(path, index=best)
    if errors is not None:
        try:
            add_errors_as_bfactors(
                path,
                errors,
                header=[
                    "REMARK   1 UNREFINED MODEL — OpenMM refinement skipped "
                    "(refine_structures=false)"
                ],
            )
        except Exception as exc:
            warn(
                f"{clone_id}: wrote the unrefined model but could not add error B-factors: "
                f"{exc}. The B-factor column in this file is not an error estimate."
            )
    return "unrefined"


def main():
    dry_run = "--dry-run" in sys.argv

    os.makedirs("./outputs", exist_ok=True)
    os.makedirs("./outputs/models", exist_ok=True)

    receptor_type = os.environ.get("PARAM_RECEPTOR_TYPE", "auto")
    if receptor_type not in {"auto", "antibody", "nanobody"}:
        fail(f"receptor_type must be auto, antibody or nanobody, got '{receptor_type}'")
    refine = os.environ.get("PARAM_REFINE_STRUCTURES", "true").lower() == "true"

    annotations = read_json("./inputs/annotations.json", "annotations.json")
    clones = annotations.get("clones") or []
    if not clones:
        fail("annotations.json contains no annotated clones.")
    numbering_scheme = (annotations.get("settings") or {}).get("numbering_scheme", "imgt")

    settings = {
        "receptor_type": receptor_type,
        "refine_structures": refine,
        "numbering_scheme": numbering_scheme,
        "n_clones": len(clones),
        "immunebuilder_version": package_version("ImmuneBuilder"),
        "torch_threads": os.environ.get("OMP_NUM_THREADS", "default"),
    }
    print(f"ImmuneBuilder version: {settings['immunebuilder_version']}", flush=True)
    print(f"Resolved settings: {json.dumps(settings, indent=2)}", flush=True)
    if not refine:
        print(
            "NOTE: refine_structures=false. Models are written unrefined, which is several "
            "times faster and leaves minor geometry violations in place. Node 03's SASA-derived "
            "metrics are sensitive to exactly that geometry, so use refined models for numbers "
            "you intend to act on.",
            flush=True,
        )

    if dry_run:
        kinds = {}
        for record in clones:
            kind = record["receptor_type"] if receptor_type == "auto" else receptor_type
            kinds[kind] = kinds.get(kind, 0) + 1
        print(f"Would model: {json.dumps(kinds)}", flush=True)
        print("Dry run requested — not loading model weights.", flush=True)
        return

    from ImmuneBuilder import ABodyBuilder2, NanoBodyBuilder2

    needed = {
        (record["receptor_type"] if receptor_type == "auto" else receptor_type)
        for record in clones
    }
    predictors = {}
    for kind in sorted(needed):
        print(f"Loading {kind} weights", flush=True)
        try:
            predictors[kind] = (
                ABodyBuilder2(numbering_scheme=numbering_scheme)
                if kind == "antibody"
                else NanoBodyBuilder2(numbering_scheme=numbering_scheme)
            )
        except Exception as exc:
            fail(
                f"Could not load {kind} weights: {exc}. They are baked into the image at build "
                "time, so this is an image-build problem rather than a network one."
            )

    rows, failures, errors_by_clone = [], [], {}
    started_all = time.time()

    for record in clones:
        clone_id = record["clone_id"]
        kind = record["receptor_type"] if receptor_type == "auto" else receptor_type
        chains = record.get("chains") or {}
        heavy = (chains.get("VH") or {}).get("domain_sequence")
        light = (chains.get("VL") or {}).get("domain_sequence")

        if not heavy:
            failures.append({"clone_id": clone_id, "reason": "no heavy chain in annotations"})
            continue
        if kind == "antibody" and not light:
            failures.append(
                {
                    "clone_id": clone_id,
                    "reason": "receptor_type=antibody but the annotation has no light chain",
                }
            )
            continue

        sequence_dict = {"H": heavy} if kind == "nanobody" else {"H": heavy, "L": light}
        path = f"./outputs/models/{clone_id}.pdb"
        started = time.time()
        try:
            prediction = predictors[kind].predict(sequence_dict)
            mode = save_prediction(prediction, path, refine, clone_id)
        except AssertionError as exc:
            # ImmuneBuilder asserts on short sequences and unrecognised chain types.
            failures.append({"clone_id": clone_id, "reason": f"ImmuneBuilder rejected it: {exc}"})
            print(f"  {clone_id}: FAILED — {exc}", flush=True)
            continue
        except Exception as exc:
            failures.append({"clone_id": clone_id, "reason": f"{type(exc).__name__}: {exc}"})
            print(f"  {clone_id}: FAILED — {exc}", flush=True)
            continue
        elapsed = round(time.time() - started, 2)

        errors = per_residue_errors(prediction)
        mean_error = float(errors.mean()) if errors is not None else None
        max_error = float(errors.max()) if errors is not None else None

        # CDR-H3 is both the loop that matters most and the one modelled worst, so its error
        # is surfaced separately rather than hidden in the mean.
        cdrh3_error = None
        if errors is not None:
            heavy_rows = (chains.get("VH") or {}).get("residues") or []
            offsets = [
                index
                for index, row in enumerate(heavy_rows)
                if row.get("region") == "CDR3"
            ]
            if offsets and max(offsets) < len(errors):
                cdrh3_error = float(
                    sum(errors[index] for index in offsets) / len(offsets)
                )

        if errors is not None:
            errors_by_clone[clone_id] = [round(float(value), 3) for value in errors]

        rows.append(
            {
                "clone_id": clone_id,
                "receptor_type": kind,
                "model_file": f"models/{clone_id}.pdb",
                "refinement": mode,
                "n_residues_modelled": len(errors) if errors is not None else "",
                "mean_error_angstrom": round(mean_error, 3) if mean_error is not None else "",
                "max_error_angstrom": round(max_error, 3) if max_error is not None else "",
                "cdrh3_mean_error_angstrom": (
                    round(cdrh3_error, 3) if cdrh3_error is not None else ""
                ),
                "seconds": elapsed,
                "status": "ok",
                "reason": "",
            }
        )
        print(
            f"  {clone_id}: {kind}, {mode}, mean error "
            f"{round(mean_error, 2) if mean_error is not None else '?'} Å, {elapsed}s",
            flush=True,
        )

    for failure in failures:
        rows.append(
            {
                "clone_id": failure["clone_id"],
                "status": "failed",
                "reason": failure["reason"],
                "refinement": "",
                "model_file": "",
            }
        )

    fields = [
        "clone_id",
        "receptor_type",
        "model_file",
        "refinement",
        "n_residues_modelled",
        "mean_error_angstrom",
        "max_error_angstrom",
        "cdrh3_mean_error_angstrom",
        "seconds",
        "status",
        "reason",
    ]
    with open("./outputs/model_summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    modelled = [row for row in rows if row["status"] == "ok"]
    manifest = {
        "settings": settings,
        "n_modelled": len(modelled),
        "n_failed": len(failures),
        "failures": failures,
        "per_residue_errors": errors_by_clone,
        "elapsed_seconds": round(time.time() - started_all, 1),
        "bfactor_meaning": (
            "The B-factor column holds ImmuneBuilder's per-residue error estimate in Angstrom "
            "(the spread of its four-model ensemble), not a crystallographic B-factor. Higher "
            "means the model is less sure; CDR-H3 is usually highest."
        ),
        "caveat": (
            "Upstream reports mean CDR-H3 RMSD of 2.81 A on a 34-antibody benchmark. Every "
            "structure-derived metric downstream inherits that error, and CDR-H3 carries most "
            "of it."
        ),
    }
    with open("./outputs/modelling_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    if not modelled:
        fail(
            f"No clone could be modelled ({len(failures)} failure(s)). See model_summary.csv "
            "for the per-clone reasons."
        )

    print(
        f"Modelled {len(modelled)}/{len(clones)} clone(s) in "
        f"{manifest['elapsed_seconds']}s ({len(failures)} failure(s))",
        flush=True,
    )
    if failures:
        warn(
            f"{len(failures)} clone(s) have no structure. They keep their node 01 annotation "
            "and node 04 humanness, and node 03 will report their metrics as unavailable."
        )


if __name__ == "__main__":
    main()
