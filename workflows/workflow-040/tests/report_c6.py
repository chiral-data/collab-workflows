#!/usr/bin/env python3
"""Judge C6: a real prediction killed at a short budget.

Five claims, and the second one is the reason this control exists at all:

  1. the node stopped near the budget rather than running to completion;
  2. the predictor's workers WERE alive mid-prediction, holding GPU memory
     -- without this, claim 3 is vacuous (§8.4);
  3. nothing of the predictor survived inside the container after the node returned;
  4. the GPU came back to its baseline;
  5. every design is recorded FAILED with a reason naming the timeout and the budget,
     not silently dropped (§8.5) and not left with a metric defaulted to a number (§8.7).

Claim 3 is answerable only because run_c6.sh starts the container first and lets it
outlive the node. Under `docker run`, the container's exit kills everything in it and
frees the device whatever kill_tree() did, so the check would pass for the wrong reason.

Run:  python3 -I tests/report_c6.py <c6-run-dir>
"""
import json
import os
import re
import sys

PREDICTOR = re.compile(r"colabfold_batch|colabfold\.batch|alphafold", re.I)
# The sleep that holds the container open, the listing command itself, and the shell
# docker exec starts are not the predictor. Matching them would make claim 3 fail for the
# wrong reason.
HARMLESS = re.compile(r"\bsleep infinity\b|\bps -eo\b|/bin/sh|/bin/bash|\[", re.I)
# A listing that could not be taken is not an empty listing. The first version of this
# control ran `docker exec ps` against an image with no ps, so the after-listing held
# "OCI runtime exec failed: ... executable file not found" -- no process lines, and the
# judge read that as "nothing survived" and PASSED. That is §8.4 inside the checking code:
# absence of evidence rendered as evidence of absence. `docker top` prints a PID header,
# so a real listing is now something the judge can require rather than assume.
UNREADABLE = re.compile(
    r"OCI runtime|executable file not found|no such container|is not running|"
    r"Cannot connect to the Docker daemon|permission denied", re.I)
HEADER = re.compile(r"^\s*PID\b", re.I | re.M)


def listing_usable(text):
    """(ok, why). A process listing must be shown to be a listing before it is read."""
    if not text or not text.strip():
        return False, "the listing is empty -- not 'no processes', but no output at all"
    hit = UNREADABLE.search(text)
    if hit:
        return False, (f"the listing could not be taken: {hit.group(0)!r}\n"
                       f"first line: {text.strip().splitlines()[0][:140]}")
    if not HEADER.search(text):
        return False, ("the listing carries no PID header, so it is not a process "
                       "table:\n" + text.strip()[:200])
    return True, "a real process table, with its PID header"


def read(root, name):
    path = os.path.join(root, name)
    if not os.path.exists(path):
        return None
    with open(path, errors="replace") as handle:
        return handle.read()


def gpu_mib(text):
    """Total MiB in use across the devices in nvidia-smi csv output.

    The query is name,memory.used,memory.total, so the FIRST MiB figure on each line is
    the used figure and the second is the capacity. Taking both would report every device
    as full.
    """
    if not text:
        return None
    total, found = 0, False
    for line in text.strip().splitlines():
        hit = re.search(r"(\d+)\s*MiB", line)
        if hit:
            total += int(hit.group(1))
            found = True
    return total if found else None


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]

    budget = int((read(root, "budget.txt") or "0").strip() or 0)
    elapsed = int((read(root, "elapsed.txt") or "0").strip() or 0)
    node_rc = (read(root, "node-rc.txt") or "?").strip()
    sampler = read(root, "sampler.txt") or ""
    ps_after = read(root, "ps-after.txt") or ""
    before = gpu_mib(read(root, "gpu-before.txt"))
    after = gpu_mib(read(root, "gpu-after.txt"))

    manifest = None
    for name in ("prediction_manifest.json", "predictions.json"):
        path = os.path.join(root, "n3", "outputs", name)
        if os.path.exists(path):
            with open(path) as handle:
                manifest = json.load(handle)
            break

    width = 78
    print("=" * width)
    print("workflow-040 C6: a real prediction killed at a short budget")
    print(f"  run:      {root}")
    print(f"  budget:   {budget}s   node 03 exit code {node_rc}, wall clock {elapsed}s")
    print(f"  gpu:      {before} MiB before -> {after} MiB after")
    print("=" * width)

    results = []

    def claim(name, ok, detail):
        results.append((name, ok, detail))
        print(f"\n  [{'PASS' if ok else 'FAIL'}] {name}")
        for line in detail.splitlines():
            print(f"         {line}")

    # 1 -------------------------------------------------------------------------------
    # A real prediction of these inputs takes 45-67s per design, so finishing inside the
    # budget would mean it was never actually killed.
    n_designs = len(manifest.get("predictions", [])) if manifest else 0
    ceiling = max(budget * max(n_designs, 1) * 3, budget * 3) + 120
    claim(
        "stopped at the budget rather than completing",
        0 < elapsed <= ceiling and elapsed >= budget,
        f"ran {elapsed}s against a {budget}s budget over {n_designs} design(s); "
        f"a completed prediction of these inputs takes 45-67s EACH, so an elapsed "
        f"time at or just above the budget per design is the kill landing.\n"
        f"ceiling allowed {ceiling}s",
    )

    # 2 -------------------------------------------------------------------------------
    live = [line.strip() for line in sampler.splitlines()
            if PREDICTOR.search(line) and not HARMLESS.search(line)]
    sampled_gpu = [int(m.group(1)) for m in re.finditer(r"gpu:\s*(\d+)\s*MiB", sampler)]
    peak = max(sampled_gpu) if sampled_gpu else None
    sampler_readable, sampler_why = listing_usable(sampler)
    gpu_filled = peak is not None and before is not None and peak > before + 100
    claim(
        "the predictor was really running (so claim 3 is not vacuous)",
        sampler_readable and bool(live) and gpu_filled,
        f"{len(live)} predictor process line(s) sampled inside the container while the "
        f"node ran; peak GPU {peak} MiB against a {before} MiB baseline.\n"
        f"If nothing was ever alive and the device never filled, 'nothing survived' "
        f"would be true of a node that did no work.\n"
        f"sampler listing: {sampler_why}\n"
        + (f"first:  {live[0][:150]}\n" if live else "")
        + (f"last:   {live[-1][:150]}" if live else ""),
    )

    # 3 -------------------------------------------------------------------------------
    readable, why = listing_usable(ps_after)
    survivors = [line.strip() for line in ps_after.splitlines()
                 if PREDICTOR.search(line) and not HARMLESS.search(line)]
    claim(
        "nothing of the predictor survived the node",
        readable and not survivors,
        (f"{why}, and it holds no predictor process after the node returned.\n"
         f"The container is deliberately still up: this is the claim a `docker run` "
         f"could not make,\nbecause a container's exit would have killed them "
         f"regardless of kill_tree() (§8.4)."
         if readable and not survivors else
         f"THE LISTING IS NOT USABLE, so this claim cannot be made either way: {why}"
         if not readable else
         "SURVIVORS -- kill_tree() did not reach these:\n"
         + "\n".join(line[:150] for line in survivors[:8])),
    )

    # 4 -------------------------------------------------------------------------------
    claim(
        "the GPU was freed",
        before is not None and after is not None and after <= before + 100,
        f"{before} MiB before, {after} MiB after. 037's defect (755d4af2) was workers "
        f"outliving the parent and holding device memory that had to be freed by hand.",
    )

    # 5 -------------------------------------------------------------------------------
    if manifest is None:
        claim("every design recorded failed, with the timeout named", False,
              "no prediction manifest in n3/outputs -- nothing to read")
    else:
        rows = manifest.get("predictions", [])
        bad = []
        for row in rows:
            design = row.get("design") or row.get("name") or "?"
            status = row.get("status")
            reason = (row.get("reason") or "")
            if status == "ok":
                bad.append(f"{design}: status ok -- it was supposed to be killed")
                continue
            if "timed out" not in reason and "timeout" not in reason.lower():
                bad.append(f"{design}: reason does not name the timeout: {reason[:90]!r}")
            elif str(budget) not in reason:
                bad.append(f"{design}: reason does not name the {budget}s budget: "
                           f"{reason[:90]!r}")
        claim(
            "every design recorded failed, with the timeout and budget named",
            bool(rows) and not bad,
            (f"{len(rows)} design(s), all recorded failed with the timeout and the "
             f"{budget}s budget in the reason. A design whose prediction failed must "
             f"appear as failed, not vanish (§8.5), and the reason has to be a reason "
             f"rather than a sentinel (§8.7).\n"
             + "\n".join(f"  {(r.get('design') or '?')}: {(r.get('reason') or '')[:95]}"
                         for r in rows[:6])
             if not bad else "\n".join(bad[:8])),
        )

    failed = [name for name, ok, _ in results if not ok]
    print("\n" + "=" * width)
    if failed:
        print(f"C6: NOT MET -- {len(failed)} claim(s) failed: {failed}")
    else:
        print("C6: MET -- a real colabfold_batch holding GPU memory was killed at the")
        print("budget, nothing of it survived inside a container that outlived the node,")
        print("the device came back to baseline, and every design is recorded failed with")
        print("the timeout and the budget in its reason.")
        print()
        print("This closes the Phase 1 control set except C5, and it is the one claim the")
        print("Tier 1 image run reports as a GAP rather than a pass, because a pid")
        print("recorded inside a container means nothing on the host.")
    print("=" * width)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
