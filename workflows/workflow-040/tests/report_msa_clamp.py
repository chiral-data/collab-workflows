#!/usr/bin/env python3
"""Compare the two msa_mode runs over the SAME precomputed a3m.

The question: does --msa-mode single_sequence silently discard a handed-in alignment?

What makes this worth a script rather than an eyeball is the trap in it. Both runs write a
coverage.png and both re-save the parsed a3m, because colabfold does that from
unserialize_msa()'s output BEFORE max_seq is applied. So the artifacts that look like
"an MSA was used" are present either way, and §11.1's spike read exactly those artifacts
and concluded the alignment was consumed. The only honest evidence is the model's own
numbers.

Run:  python3 -I tests/report_msa_clamp.py <probe-out-dir>
"""
import glob
import json
import os
import re
import sys

MODES = ["single_sequence", "mmseqs2_uniref_env"]


def scores(root, mode):
    hits = sorted(glob.glob(os.path.join(root, mode, "out", "*scores_rank_001*.json")))
    if not hits:
        return None, None
    with open(hits[0]) as handle:
        return json.load(handle), os.path.basename(hits[0])


def artifacts(root, mode):
    out = os.path.join(root, mode, "out")
    if not os.path.isdir(out):
        return []
    return sorted(os.listdir(out))


def logged_max_seq(root, mode):
    """colabfold logs the model config, which is where max_seq shows up if anywhere."""
    path = os.path.join(root, f"{mode}.log")
    if not os.path.exists(path):
        return None
    with open(path, errors="replace") as handle:
        text = handle.read()
    hit = re.search(r"max_seq[^\d]*(\d+)", text)
    return int(hit.group(1)) if hit else None


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]

    width = 78
    print("=" * width)
    print("Does --msa-mode single_sequence discard a precomputed a3m?")
    print(f"  probe: {root}")
    print("=" * width)

    rows = {}
    for mode in MODES:
        data, name = scores(root, mode)
        files = artifacts(root, mode)
        rows[mode] = {
            "scores": data,
            "scores_file": name,
            "coverage_png": any(f.endswith("coverage.png") for f in files),
            "a3m_resaved": any(f.endswith(".a3m") for f in files),
            "n_files": len(files),
            "max_seq": logged_max_seq(root, mode),
        }

    print(f"\n{'mode':22} {'iptm':>7} {'ptm':>7} {'max_pae':>8} "
          f"{'plddt':>7} {'coverage':>9} {'a3m':>5}")
    for mode in MODES:
        row = rows[mode]
        data = row["scores"] or {}
        plddt = data.get("plddt")
        mean_plddt = (sum(plddt) / len(plddt)) if plddt else None

        def show(value, places=3):
            return "-" if value is None else f"{value:.{places}f}"

        print(f"{mode:22} {show(data.get('iptm')):>7} {show(data.get('ptm')):>7} "
              f"{show(data.get('max_pae'), 2):>8} {show(mean_plddt, 1):>7} "
              f"{str(row['coverage_png']):>9} {str(row['a3m_resaved']):>5}")

    missing = [m for m in MODES if rows[m]["scores"] is None]
    if missing:
        print(f"\nNO SCORES for {missing} -- the run failed; read {root}/<mode>.log")
        return 1

    single = rows["single_sequence"]["scores"]
    full = rows["mmseqs2_uniref_env"]["scores"]

    print("\nthe artifacts that look like evidence, and are not")
    for mode in MODES:
        row = rows[mode]
        print(f"  {mode:22} coverage.png={row['coverage_png']}  "
              f"a3m re-saved={row['a3m_resaved']}  {row['n_files']} output file(s)")
    print("  Both are written from unserialize_msa()'s output, BEFORE max_seq is applied")
    print("  (colabfold/batch.py:284-287). So they say the a3m was PARSED, never that it")
    print("  reached the model. §11.1's spike read these and concluded the MSA was used.")

    deltas = {}
    for key in ("iptm", "ptm", "max_pae"):
        a, b = single.get(key), full.get(key)
        if a is not None and b is not None:
            deltas[key] = b - a

    print("\nthe only honest evidence: the model's own numbers")
    for key, delta in deltas.items():
        print(f"  {key:10} single_sequence {single[key]:.4f} -> "
              f"mmseqs2_uniref_env {full[key]:.4f}   delta {delta:+.4f}")

    identical = all(abs(d) < 1e-9 for d in deltas.values())

    print("\n" + "=" * width)
    if identical:
        print("NO CLAMP OBSERVED: the two runs agree to the digit, so msa_mode does not")
        print("change what the model sees when an a3m is supplied. Node 03 may pass")
        print("single_sequence on the complex path, which is the fail-safe choice (if the")
        print("a3m were ever ignored, the fallback is no egress rather than sending the")
        print("designed binder to an MSA server).")
    else:
        print("THE CLAMP IS REAL. The same a3m produces different predictions under the")
        print("two msa_mode values, so --msa-mode single_sequence DOES change what the")
        print("model is fed -- colabfold/batch.py:137-142 sets max_seq = 1 regardless of")
        print("how many sequences the alignment holds.")
        print()
        print("Consequence for 040: node 03 must NOT pass single_sequence on the complex")
        print("path. It would report interface numbers measured at effectively")
        print("single_sequence while config.json, the re-saved a3m and coverage.png all")
        print("said an alignment was used -- §4's exact failure, and silent.")
        print()
        print("mmseqs2_uniref_env is safe to pass WITH an a3m: use_msa and use_env only")
        print("feed get_msa_and_templates, which is not called when a3m_lines is not None")
        print("(batch.py:258-270). This probe ran both arms under --network none and")
        print("neither touched the network.")
    print("=" * width)
    return 0


if __name__ == "__main__":
    sys.exit(main())
