#!/usr/bin/env python3
"""Audit every node's declared ports against the files its code actually touches.

Why this exists
---------------
Node 02 wrote `cplx_*.fasta` and `pairedmsa_*.a3m` and declared neither. `silva validate`
reported "6 job(s) ok" -- it checks that the declarations are well formed, not that they
match the code -- and the node passed standing alone, because standing alone everything it
writes is in its own ./outputs/. On a real run Silva copies the DECLARED outputs and
nothing else, so the complex path would have degraded to monomer-only in silence.

The Tier 1 suite caught that one because its harness reads the globs out of job.toml. This
catches the class directly, in a second, without staging anything: it reads each script for
the literal `./outputs/<name>` and `./inputs/<name>` paths it mentions and checks each one
against the node's own declaration.

What it cannot see
------------------
Paths built at run time -- `f"./outputs/{prefix}_{design}.fasta"` -- are reported as
patterns to check by eye rather than silently passed, because a tool that quietly ignored
what it could not parse would be the §8.3 case: an absence that generates no finding. Run
it with --verbose to see them. It is a cheap second opinion, not a replacement for the
suite, and it decides nothing on its own.

Usage:  python3 tests/audit_ports.py [--verbose]
"""
import fnmatch
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WF = os.path.dirname(HERE)
NODES = ["01-ingest-designs", "01b-align-target", "02-build-prediction-inputs",
         "03-predict", "04-score", "05-filter-and-report"]

LITERAL = re.compile(r"""["']\./(inputs|outputs)/([A-Za-z0-9_.*-]+)["']""")
BUILT = re.compile(r"""["']\./(inputs|outputs)/[^"']*\{""")
JOIN = re.compile(r"""os\.path\.join\(\s*["']\./(inputs|outputs)["']\s*,\s*([^)]+)\)""")

# A node may write a file it does not declare -- but it has to say why, here, in writing.
# Keyed by "<node>:<what the code calls it>", valued with the reason. An undeclared write
# that is not in this table is a finding, because the alternative is a table of exceptions
# that grows by accident until it covers the next cplx_*.fasta.
DELIBERATELY_UNDECLARED = {
    "03-predict:log_name": (
        "<design>_predict.log, the predictor's own output, retained for diagnosis (§8.5) "
        "and read by nothing downstream. Declaring it would copy it into node 04's "
        "inputs, and a dependency that pulls inputs it does not read is 7fba161c's defect "
        "(§8.10). The manifest records log_dir so the pointer resolves."
    ),
}


def declared(node, key):
    path = os.path.join(WF, node, ".chiral", "job.toml")
    text = open(path).read()
    match = re.search(r"^%s\s*=\s*\[(.*?)\]" % key, text, re.S | re.M)
    if not match:
        return None
    return re.findall(r'"([^"]+)"', match.group(1))


def covered(name, patterns):
    return any(fnmatch.fnmatch(name, p) for p in patterns)


def main(argv):
    verbose = "--verbose" in argv
    problems = []
    for node in NODES:
        source = open(os.path.join(WF, node, "script.py")).read()
        print(f"\n{node}")
        for key in ("inputs", "outputs"):
            ports = declared(node, key)
            if ports is None:
                problems.append(f"{node}: job.toml declares no {key} list at all")
                print(f"  {key:<8} NO LIST IN job.toml")
                continue
            print(f"  {key:<8} declared {ports}")
            seen = sorted({m.group(2) for m in LITERAL.finditer(source)
                           if m.group(1) == key})
            for name in seen:
                if covered(name, ports):
                    if verbose:
                        print(f"           ok      {name}")
                else:
                    problems.append(
                        f"{node}: the code touches ./{key}/{name} and job.toml does not "
                        f"declare it. Silva copies the declaration, not the directory."
                    )
                    print(f"           MISSING {name}")
        # Run-time-built names. The expression cannot be resolved here, so each one is
        # judged on whether the node declares a glob that could plausibly cover it, and
        # anything left over has to appear in DELIBERATELY_UNDECLARED with a reason.
        out_ports = declared(node, "outputs") or []
        for match in JOIN.finditer(source):
            if match.group(1) != "outputs":
                continue
            expression = " ".join(match.group(2).split()).strip().rstrip(",")
            key = f"{node}:{expression}"
            if key in DELIBERATELY_UNDECLARED:
                print(f"  undeclared write, BY DECISION: {expression}")
                if verbose:
                    print(f"           {DELIBERATELY_UNDECLARED[key]}")
                continue
            print(f"  built at run time: ./outputs/<{expression}> -- "
                  f"covered by one of {out_ports}?")

        built = [" ".join(m.group(0).split()) for m in BUILT.finditer(source)]
        if verbose and built:
            print(f"  f-string paths, check by eye: {len(built)}")
            for item in built[:12]:
                print(f"           {item[:96]}")
        elif built:
            print(f"  f-string paths, not checked here: {len(built)} (--verbose to list)")

    print("\n" + "=" * 78)
    if problems:
        print(f"{len(problems)} undeclared port(s):")
        for line in problems:
            print(f"  - {line}")
        print("=" * 78)
        return 1
    print("Every literal ./inputs/ and ./outputs/ path is covered by its node's own")
    print("declaration. Run-time-built names are listed above and are not checked here.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
