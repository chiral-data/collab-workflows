#!/usr/bin/env python3
"""01b-align-target -- Phase 2. Align the target chain once, and cache it.

WHY THIS IS A NODE AND NOT A FUNCTION

Two reasons, and the second is the one that matters.

Caching: §4's decision is that Phase 2 aligns the TARGET ONLY, once per run, and pairs the
unaligned designed binder against it. Batches of designs are the entire point of 040, so
an alignment per run rather than per design is the difference between one network call and
a hundred.

Egress: this is the ONLY node in the workflow permitted to touch the network, and it says
so in its own job.toml. Exactly one sequence goes out: the fused target chain, which is
public native protein from the PDB.

What keeps the designed binders out of that payload is a CHECK, not their absence. This
docstring used to claim the opposite -- that the designed sequences were not inputs here,
and that this was "not a promise about the code's behaviour, it is a fact about its
inputs". That was false, and it was the most confident sentence in the file.
design_inventory.json carries `binder_sequence` for every design and it is this node's
only input. The binder stayed out of the payload because `target_of()` returns
target-derived fields only and `from_server()` reads `target['sequence']`: true, correct,
and one edit away from not being true, in the one node that can open a socket.
assert_no_designed_sequence() now enforces it immediately before the POST, and
tests/fixtures.json shows it refusing. Workflows 029 had to retrofit a declared egress
boundary after the fact (92cf05b5, 8a6dcebf); 040 starts with one, and now also starts
with a reason to believe it.

WHAT §4 DECIDED, AND WHY A FALLBACK HERE WOULD BE THE WORST POSSIBLE BUG

029 measured a native complex at both settings, changing only msa_mode (2361277f):

                        single_sequence   mmseqs2_uniref_env
    iPTM                0.100             0.860
    interface PAE       27.6 A            7.7 A
    receptor pLDDT      38.4              83.7

The mechanism is NOT interface co-evolution -- a designed binder has no homologs, so there
is nothing to pair against. It is the TARGET being modelled at ~40 pLDDT instead of ~84,
and 8e4772a7's rule doing the rest: interface PAE can be no better than the worse-modelled
of the two chains.

So a Phase 2 run whose target alignment silently failed, and which fell back to
single_sequence, would report interface numbers that are a measurement of the fallback
rather than of the design -- and nothing in the output would say so. THIS NODE NEVER FALLS
BACK. It either produces an alignment with homologs in it, or it fails and the run stops.
An alignment containing only the query is single_sequence wearing a hat, and is refused by
name (§8.2: a fallback that gives the right answer is still a defect, and this one does not
even give the right answer).

INPUTS
    ./inputs/design_inventory.json      from 01-ingest-designs

OUTPUTS
    ./outputs/target_msa.a3m            single-chain a3m: query first, then homologs
    ./outputs/target_alignment.json     provenance: source, counts, digest, the fusion

PARAMETERS
    PARAM_SCORE_INTERFACE     false   this node is Phase 2; it no-ops when false
    PARAM_MSA_SOURCE          file    'file' or 'server'
    PARAM_TARGET_MSA_FILE     -       with source=file: the a3m under ./inputs/
    PARAM_MSA_HOST_URL        https://api.colabfold.com
    PARAM_MSA_TIMEOUT_S       1200
    PARAM_MIN_MSA_SEQUENCES   2       the query plus at least one homolog

'file' is the default on purpose. The server path is real egress and has to be asked for.
"""
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

AA = set("ACDEFGHIKLMNPQRSTVWYX")
# a3m lowercase letters are insertions relative to the query; '-' is a deletion.
A3M_EXTRA = set("abcdefghijklmnopqrstuvwxyz-")

DEFAULT_HOST = "https://api.colabfold.com"


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr, flush=True)
    raise SystemExit(1)


def bool_param(name, key, default):
    raw = os.environ.get(name)
    if raw is None:
        return default
    value = raw.strip().lower()
    if value in ("1", "true", "yes", "on"):
        return True
    if value in ("0", "false", "no", "off"):
        return False
    fail(f"{name}={raw!r} is not a boolean. Use true or false.")


def int_param(name, default):
    raw = os.environ.get(name)
    if raw is None or not raw.strip():
        return default
    try:
        value = int(raw.strip())
    except ValueError:
        fail(f"{name}={raw!r} is not an integer.")
    if value <= 0:
        fail(f"{name}={raw!r} must be positive.")
    return value


def read_a3m(path):
    """[(name, sequence)] from an a3m or fasta. Blank lines and '#' headers skipped.

    Hand-rolled rather than biopython, for the reason §8.4 gives: biopython 1.85 only
    WARNS on a FASTA comment line, so a file with a ';' banner in it parses to something
    plausible and different. And the '#<lens>\\t<cards>' first line of a PAIRED a3m is not
    a comment at all -- it is cardinality information, and a reader that silently dropped
    it would turn a complex alignment into a single-chain one without saying so.
    """
    with open(path, errors="replace") as handle:
        return parse_a3m(handle, path)


def parse_a3m(lines, where):
    """The parser itself, over any iterable of lines. `where` only names the source.

    Split out of read_a3m so the SERVER path never has to put what it fetched on disk in
    order to read it. It used to write the downloaded bytes into ./outputs/_server.a3m and
    delete the file afterwards -- a download parked, however briefly, inside this node's
    own declared output directory, and left there for good if the node died in between.
    Nothing needs it to be a file, so it is not one.
    """
    records, name, chunks = [], None, []
    header = None
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if header is None:
                header = line
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name, chunks = line[1:].strip(), []
            continue
        if name is None:
            fail(f"{where}: a sequence line appears before any '>' header: {line[:60]!r}")
        chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return header, records


def ungapped(sequence):
    """The residues only: a3m insertions are lowercase and deletions are '-'."""
    return "".join(c for c in sequence if c.isupper() and c != "-")


def target_of(inventory):
    """The one target this run is about, fused and disclosed.

    One target per run (§4), so this refuses a batch whose designs disagree about what
    they are binding rather than aligning the first one and carrying on -- which would
    produce interface numbers for the wrong receptor and say nothing.

    A multi-chain or internally gapped target is FUSED, not trimmed, and the fusion is
    recorded and warned about. 8e4772a7 trimmed a gapped receptor to the segment carrying
    the epitope and measured ~42 pLDDT against 81 for the whole chain, so the fused
    junction is a real defect kept knowingly (§8.12).
    """
    live = [row for row in inventory["designs"] if not row.get("skipped")]
    if not live:
        fail("every design was skipped, so there is no target to align.")

    seen = {}
    for row in live:
        chains = row.get("target_chains") or []
        if not chains:
            fail(
                f"design {row['design']!r} declares no target chain, so there is nothing "
                f"to align. designed_chain_source={row.get('designed_chain_source')!r}. "
                "A single-chain design cannot be scored on the complex path (§5) -- node "
                "02 refuses it, and this node refuses to invent a receptor for it."
            )
        key = tuple((c["chain"], c["sequence"]) for c in chains)
        seen.setdefault(key, []).append(row["design"])

    if len(seen) > 1:
        detail = "; ".join(
            f"{designs} -> {[c for c, _ in key]} "
            f"({sum(len(s) for _, s in key)} residues)"
            for key, designs in seen.items()
        )
        fail(
            "the designs do not share one target, and §4 aligns one target per run: "
            f"{detail}. Aligning the first and scoring the rest against it would report "
            "an interface for the wrong receptor. Split the batch by target."
        )

    key = next(iter(seen))
    chains = (live[0].get("target_chains") or [])
    sequence = "".join(c["sequence"] for c in chains)

    bad = sorted(set(sequence) - AA)
    if bad:
        fail(f"the target sequence carries characters outside the alphabet: {bad}")
    if not sequence:
        fail("the fused target sequence is empty.")

    internal_gaps = [
        {"chain": c["chain"], **gap} for c in chains for gap in (c.get("gaps") or [])
    ]
    fusion = []
    if len(chains) > 1:
        fusion.append(
            f"{len(chains)} target chains ({', '.join(c['chain'] for c in chains)}) are "
            f"concatenated into one query of {len(sequence)} residues"
        )
    for gap in internal_gaps:
        fusion.append(
            f"chain {gap['chain']} is missing {gap['missing']} residue(s) between "
            f"{gap['after']} and {gap['before']}, and the two sides are fused"
        )

    return {
        "chains": [c["chain"] for c in chains],
        "chain_lengths": [len(c["sequence"]) for c in chains],
        "sequence": sequence,
        "n_residues": len(sequence),
        "internal_gaps": internal_gaps,
        "fusion": fusion,
        "designs": sorted(d for designs in seen.values() for d in designs),
    }


def designed_sequences(inventory):
    """Every designed sequence in the inventory, with the design that owns it.

    These are the strings that must never leave this node. They are collected from the
    same inventory the target is read from, because the thing a guard has to be given is
    the real payload-adjacent data, not a restatement of it (§8.10).
    """
    out = []
    for row in inventory.get("designs") or []:
        sequence = (row.get("binder_sequence") or "").strip().upper()
        if sequence:
            out.append((row.get("design"), sequence))
    return out


def assert_no_designed_sequence(query, forbidden):
    """Refuse to send a payload containing a designed sequence.

    THE CLAIM THIS REPLACES WAS FALSE, and it was load-bearing. job.toml and
    tests/README.md both said a designed sequence "cannot leak" from this node because it
    "is not an input here" -- structural, not a promise. It is an input:
    design_inventory.json carries `binder_sequence` for every design, and it is this
    node's only input. What actually kept the binder out of the payload was that
    `target_of()` returns a dict of target-derived fields and `from_server()` reads
    `target['sequence']` from it. True, and one edit away from not being true, in the only
    node in this workflow permitted to open a socket.

    §8.4's shape exactly: a check -- here, a documented safety property -- that holds for
    a reason other than the one given. So the property is now enforced where it matters,
    immediately before the POST, against the sequences read from the actual inventory.

    This is a guard, so it can be shown refusing: tests/fixtures.json case
    `egress-refuses-a-designed-sequence` crafts an inventory whose fused target CONTAINS a
    binder sequence and requires this to fail by name, before any socket is opened.
    """
    leaked = sorted({design for design, sequence in forbidden if sequence in query})
    if leaked:
        fail(
            "REFUSING TO SEND: the outgoing MSA query contains the designed sequence of "
            f"{leaked}. This node is the only egress in the workflow and a designed "
            "sequence must never cross it -- a de novo binder is unpublished work and the "
            "MSA server is a third party. Nothing has been sent; no socket was opened.\n"
            "  This is checked rather than assumed because the sequences ARE present in "
            "this node's input: design_inventory.json carries binder_sequence for every "
            "design. The earlier claim that a leak was structurally impossible was wrong."
        )


def from_file(target, path):
    """An alignment supplied as a file. Validated against the target, not trusted."""
    if not os.path.exists(path):
        fail(
            f"PARAM_MSA_SOURCE=file but there is no alignment at {path}. "
            "This node does not fall back to single_sequence: see the module docstring."
        )
    header, records = read_a3m(path)
    if not records:
        fail(f"{path} holds no sequences.")

    query_name, query = records[0]
    if ungapped(query) != target["sequence"]:
        fail(
            f"{path}: the first record ({query_name!r}, {len(ungapped(query))} residues) "
            f"is not the target's sequence ({target['n_residues']} residues). An "
            "alignment of something else would model the wrong receptor and the interface "
            "numbers would be about that. Expected "
            f"{target['sequence'][:40]}... got {ungapped(query)[:40]}..."
        )
    if header is not None and "," in header:
        fail(
            f"{path}: the first line is {header!r}, which is a PAIRED a3m's "
            "chain-length header. This node wants the TARGET's own single-chain "
            "alignment; node 02 is what pairs it with a binder."
        )
    return records, {"source": "file", "path": os.path.basename(path)}


def from_server(target, host, timeout_s, forbidden):
    """The MSA server. The ONLY egress in this workflow, and only the target goes out.

    What is sent is one sequence: the fused target chain, which is public native protein
    from the PDB.

    The designed binders ARE in this node's input -- design_inventory.json carries
    binder_sequence for every design -- so "they cannot be sent" is enforced by
    assert_no_designed_sequence() below, not by their absence. See that function for the
    claim it replaced.
    """
    query = f">target\n{target['sequence']}\n"

    # Before the socket, not after. A guard that runs after the POST is a log entry.
    assert_no_designed_sequence(query, forbidden)

    print(
        f"  NETWORK: submitting {target['n_residues']} target residues to {host}. "
        "This is the only egress in the workflow. The payload was checked for every "
        "designed sequence in the inventory before this POST and carries none; the "
        "binders ARE in this node's input, so that is a check, not an absence.",
        flush=True,
    )

    def post(url, data):
        request = urllib.request.Request(
            url, data=data.encode(), method="POST",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.loads(response.read().decode())

    try:
        ticket = post(
            f"{host}/ticket/msa",
            urllib.parse.urlencode({"q": query, "mode": "env"}),
        )
    except (urllib.error.URLError, OSError, ValueError) as error:
        fail(
            f"could not reach the MSA server at {host}: {error}. Phase 2 does not "
            "proceed without a target alignment (§4), and this node will not substitute "
            "single_sequence for one. Supply an alignment with PARAM_MSA_SOURCE=file, or "
            "run Phase 1."
        )

    ticket_id = ticket.get("id")
    if not ticket_id:
        fail(f"the MSA server returned no ticket id: {ticket}")

    deadline = time.monotonic() + timeout_s
    state = ticket.get("status")
    while state in ("RUNNING", "PENDING", None):
        if time.monotonic() > deadline:
            fail(
                f"the MSA server did not finish within {timeout_s}s (last status "
                f"{state!r}). Not falling back to single_sequence."
            )
        time.sleep(5)
        try:
            with urllib.request.urlopen(
                f"{host}/ticket/{ticket_id}", timeout=60
            ) as response:
                state = json.loads(response.read().decode()).get("status")
        except (urllib.error.URLError, OSError, ValueError) as error:
            fail(f"lost contact with the MSA server while polling: {error}")
        print(f"  server status: {state}", flush=True)

    if state != "COMPLETE":
        fail(f"the MSA server finished in state {state!r} rather than COMPLETE.")

    # The result is a tar of a3m files. Written to scratch and read back, so the bytes
    # that were received are the bytes that were parsed.
    import io
    import tarfile

    try:
        with urllib.request.urlopen(
            f"{host}/result/download/{ticket_id}", timeout=timeout_s
        ) as response:
            blob = response.read()
    except (urllib.error.URLError, OSError) as error:
        fail(f"could not download the alignment: {error}")

    records = []
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:*") as archive:
        for member in archive.getmembers():
            if not member.name.endswith(".a3m"):
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            # Parsed in memory. The member's own name is never used as a path and nothing
            # is extracted to disk, so a tar entry called ../../anything is just a member
            # this loop skips or reads -- not a file it writes.
            text = extracted.read().decode("utf-8", errors="replace")
            _, records = parse_a3m(
                text.splitlines(), f"the a3m member {member.name!r} of the server's tar"
            )
            if records:
                break

    if not records:
        fail("the server's archive held no a3m with sequences in it.")
    if ungapped(records[0][1]) != target["sequence"]:
        fail(
            "the server returned an alignment whose query is not the sequence that was "
            "sent. Refusing it rather than modelling something else."
        )
    return records, {"source": "server", "host": host, "ticket": ticket_id}


def main():
    os.makedirs("./outputs", exist_ok=True)

    score_interface = bool_param("PARAM_SCORE_INTERFACE", "score_interface", False)
    inventory_path = "./inputs/design_inventory.json"

    if not score_interface:
        # A no-op that says so. An empty output directory and a zero exit code are
        # indistinguishable from a node that ran and found nothing (§8.3).
        note = {
            "score_interface": False,
            "aligned": False,
            "reason": (
                "score_interface is false, so this run is Phase 1 and predicts the "
                "binder alone with msa_mode=single_sequence, which §4 says is correct "
                "rather than a compromise for a de novo binder. No alignment is needed "
                "and no network was touched."
            ),
            "network_used": False,
        }
        with open("./outputs/target_alignment.json", "w") as handle:
            json.dump(note, handle, indent=2)
        print(
            "01b-align-target: score_interface is false, so there is nothing to align. "
            "Phase 1 predicts the binder alone (§4). Recorded in "
            "target_alignment.json; no network was touched.",
            flush=True,
        )
        return 0

    if not os.path.exists(inventory_path):
        fail(f"no {inventory_path}; 01-ingest-designs must run first.")
    with open(inventory_path) as handle:
        inventory = json.load(handle)

    target = target_of(inventory)
    print(
        f"01b-align-target: one target, {target['n_residues']} residues over chain(s) "
        f"{target['chains']}, shared by {len(target['designs'])} design(s).",
        flush=True,
    )
    for line in target["fusion"]:
        print(
            f"  WARNING: {line}. The target is folded WHOLE rather than trimmed to a "
            "segment: 8e4772a7 measured ~42 pLDDT trimmed against 81 whole, so the fused "
            "junction is a real defect kept knowingly (§8.12).",
            flush=True,
        )

    source = (os.environ.get("PARAM_MSA_SOURCE") or "file").strip().lower()
    minimum = int_param("PARAM_MIN_MSA_SEQUENCES", 2)

    if source == "file":
        name = (os.environ.get("PARAM_TARGET_MSA_FILE") or "target_msa.a3m").strip()
        records, provenance = from_file(target, os.path.join("./inputs", name))
        provenance["network_used"] = False
    elif source == "server":
        host = (os.environ.get("PARAM_MSA_HOST_URL") or DEFAULT_HOST).strip()
        records, provenance = from_server(
            target, host, int_param("PARAM_MSA_TIMEOUT_S", 1200),
            designed_sequences(inventory),
        )
        provenance["network_used"] = True
    else:
        fail(
            f"PARAM_MSA_SOURCE={source!r} is not recognised. Use 'file' or 'server'. "
            "Guessing would mean guessing whether to touch the network."
        )

    # The refusal §4 exists for. An alignment of one sequence IS single_sequence, and a
    # Phase 2 run on it would report interface numbers measuring the prediction mode
    # rather than the design -- 029's 0.100 iPTM against 0.860, and nothing in the output
    # would say which it was.
    if len(records) < minimum:
        fail(
            f"the alignment holds {len(records)} sequence(s) and at least {minimum} are "
            f"required (the query plus one homolog). An alignment containing only the "
            "query is single_sequence under another name, and §4's whole argument is "
            "that interface numbers measured that way are about the prediction mode and "
            "not about the design: 029 measured iPTM 0.100 against 0.860 and receptor "
            "pLDDT 38.4 against 83.7 on the same complex. Refusing rather than "
            "proceeding and reporting them."
        )

    query_name, query = records[0]
    with open("./outputs/target_msa.a3m", "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")

    digest = hashlib.sha256(
        open("./outputs/target_msa.a3m", "rb").read()
    ).hexdigest()

    record = {
        "score_interface": True,
        "aligned": True,
        "target_chains": target["chains"],
        "target_chain_lengths": target["chain_lengths"],
        "target_sequence": target["sequence"],
        "target_n_residues": target["n_residues"],
        "target_internal_gaps": target["internal_gaps"],
        "target_fusion": target["fusion"],
        "designs_sharing_target": target["designs"],
        "n_sequences": len(records),
        "n_homologs": len(records) - 1,
        "query_name": query_name,
        "min_required": minimum,
        "a3m_sha256": digest,
        "a3m": "target_msa.a3m",
        # Stated rather than implied: the binder is never aligned and never sent.
        "binder_aligned": False,
        "binder_sent_anywhere": False,
        "why_binder_unaligned": (
            "A de novo designed binder has no homologs, so there is no binder MSA to "
            "pair against the target's. §4: the deficit single_sequence causes is in the "
            "TARGET being modelled at ~40 pLDDT instead of ~84, and it is addressable "
            "without aligning the binder."
        ),
        **provenance,
    }
    with open("./outputs/target_alignment.json", "w") as handle:
        json.dump(record, handle, indent=2)

    print(
        f"01b-align-target: {len(records)} sequence(s) "
        f"({len(records) - 1} homolog(s)) for a {target['n_residues']}-residue target, "
        f"from {record['source']}. sha256={digest[:16]}. The binder is not aligned "
        f"(§4) and no designed sequence was in the payload, which is checked rather "
        f"than assumed -- see assert_no_designed_sequence().",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
