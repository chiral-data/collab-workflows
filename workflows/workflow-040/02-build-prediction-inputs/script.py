#!/usr/bin/env python3
"""Write the prediction inputs: one design per file, always.

Contract
--------
Reads   design_inventory.json         from 01-ingest-designs
        target_alignment.json         from 01b-align-target
        target_msa.a3m                from 01b-align-target, Phase 2 only
Writes  mono_<design>.fasta           the binder sequence alone, exactly one record
        cplx_<design>.fasta           TARGET:BINDER, one record -- Phase 2 only
        pairedmsa_<design>.a3m        what the predictor consumes -- Phase 2 only
        prediction_inputs.json        what was written, and the facts node 04 needs

The mono input is written in BOTH phases. Phase 2 ADDS interface scoring, it does not
replace self-consistency: an interface metric says nothing about whether the binder folds
as designed, which is the question this workflow exists to answer.

Two files per design on the complex path, and they are not redundant. cplx_<design>.fasta
is the declaration -- what was asked for, in the form the spec validates, with the ':'
rule enforceable on it. pairedmsa_<design>.a3m is what colabfold_batch actually reads,
in §11.1's measured format: the target carries its alignment and the binder is present
only as its own query. A bare TARGET:BINDER fasta would make the predictor align BOTH
chains, which is exactly what §4 decided against. Node 03 is handed both and cross-checks
them, so each representation is validated against the other.

Validation happens here because it is free here. A bad alphabet, a length outside bounds
or a missing chain costs nothing to catch on CPU and costs GPU minutes to catch inside
the predictor -- which is why workflow-037 made spec validation its own CPU node
(37c23bb5).

One design per file, never a batch
----------------------------------
Not a style preference. Workflow-026 takes a multi-record FASTA without complaining and
then produces a confident report about a design that does not exist: its 02_predict
copies all ColabFold output flat into one directory, and its 03_analyze_confidence globs
*_scores_rank_*.json across every design and takes job_id from whichever file sorts
first. Three designs at five models makes the reported "top 5" ranks 1,1,1,2,2 drawn from
three different designs, filed under one design's name (818954e9).

So: one record per file, one invocation per file, one output directory per design, and the
mapping back by filename prefix -- never by directory listing order.

Record ids must be filesystem-safe, because ColabFold names its outputs after them
(54bf21b6). The real sample confirms the shape:
sequences_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb.

The record id is the design name and NOTHING else -- no biosafety tag, no run id, no
length. 037 can afford a tag in its exported FASTA because nothing downstream of it
reads the file back by name; here the id becomes a filename prefix that node 03 writes
and node 04 parses, so anything appended to it has to be stripped by two more nodes that
would then each have their own idea of the format. The biosafety notice travels in
prediction_inputs.json and in node 05's export instead (7709809b).

prediction_inputs.json, per design
----------------------------------
    design              keyed on, everywhere
    path                the file written, basename only
    mode                "mono", or "cplx" for the Phase 2 complex record
    record_id           what ColabFold will name its outputs after
    chain_order         ["binder"] in Phase 1; ["target", "binder"] in Phase 2
    chain_lengths       parallel to chain_order
    expected_chains     how many chains the prediction must come back with
    expected_residues   sum(chain_lengths); node 04 asserts len(plddt) equals this

chain_order and chain_lengths exist for one reason: the scores file carries NO chain
boundaries. Measured first-hand, not inferred -- see tests/reference/: the real monomer
file carries plddt as a flat list[35] and pae as 35x35, the real complex file carries
plddt as a flat list[63] spanning both chains and pae as 63x63, and neither names a chain
anywhere. So node 04 must slice the per-residue array by position, and the position
depends on the order chosen HERE. This node decides it, so this node records it.

In Phase 1 the binder is the whole array and the slice is trivial. Record it anyway. When
Phase 2 adds TARGET:BINDER the binder becomes the TAIL -- measured: a 34+29 complex
written target-first came back with 34 residues in chain A and 29 in chain B -- and
workflow-029 slices plddt[:binder_length] because its binder is chain A
(04_fold_and_filter/script.py:207). A node 04 that inherits that assumption reports the
target's confidence as the binder's, which on a well-modelled target makes every design
look good.

Colon rules, both directions
----------------------------
A mono input must contain NO colon: a colon makes the predictor treat one design as two
chains. A Phase 2 complex input must contain EXACTLY ONE, with both halves non-empty --
zero silently predicts a monomer and reports its confidence as an interface, which is
026's bug, and two predicts three chains. 026 validates sequences against
ACDEFGHIKLMNPQRSTVWYX and so rejects the only syntax that encodes a complex (d3d707ed);
040 must accept it in the right file and refuse it in the wrong one.

The check is on the bytes written, not on the variable intended to be written. A guard
that reads the value it just assigned is asserting the assignment, not the file (§8.4).

Skips, with a warning rather than silently
------------------------------------------
A design whose sequence is empty or all-unknown is skipped rather than written as junk,
and unknown residues are counted: ColabFold folds X as a sequence gap, which depresses
confidence for reasons that have nothing to do with the design (54bf21b6). Node 01 marks
those designs; this node honours the mark rather than re-deciding it.

Re-running
----------
Must be safe or refused. Workflow-029's chain swap was not idempotent, so re-running over
a non-empty output directory corrupted backbones it had already fixed -- and it fired for
real, during a GPU smoke test, on a directory left by an earlier run (4d1a523e).

Rewriting this node's own outputs is safe, because they are derived from the inventory and
the target alignment and nothing else: same inputs, same bytes. What is NOT safe is a
mono_*.fasta or cplx_*.fasta left behind for a design that is not in this run. Node 03
globs both, so a stale file is predicted, and node 04 then scores a design the inventory
has never heard of. A stale cplx_ file is the worse of the two, because it comes back as
an interface number filed under a design name that does not exist. Both are refused, and
the refusal leaves every existing file exactly as it was.

Fails when
----------
  * a binder sequence carries a residue outside the accepted alphabet
  * a binder is too short for a superposition to mean anything
  * score_interface is true and 01b-align-target produced no alignment. NOT a fallback
    to single_sequence: §4 measured what that does to the numbers (iPTM 0.100 against
    0.860, receptor pLDDT 38.4 against 83.7), so a run that cannot align the target
    stops rather than reporting interface metrics that are about the prediction mode.
  * the target alignment holds only its own query, which is single_sequence under
    another name
  * the designs in one batch do not share a target (01b refuses this first)
  * the paired a3m's declared chain lengths disagree with the record, or a chain has no
    single-chain row -- the two ways §11.1's format fails, one silently and one loudly
  * the output directory holds a mono_*.fasta or cplx_*.fasta for a design not in this
    run
  * the file it just wrote does not pass the colon and record-count checks
"""
import json
import math
import os
import sys

# 026's alphabet, in 026's order, including X. X is accepted and counted, not rejected: a
# design carrying modified residues is a real design, and refusing it here would discard it
# where recording it costs nothing (54bf21b6). The order is kept because the string is
# quoted back in the refusal and recorded in prediction_inputs.json, and a reader comparing
# it against 026's source should see the same string rather than a sorted paraphrase.
ALPHABET_STRING = "ACDEFGHIKLMNPQRSTVWYX"
ALPHABET = set(ALPHABET_STRING)

# A Kabsch superposition needs three non-collinear points to be defined at all, and four
# before the fit stops being exactly determined by its own inputs. An RMSD over fewer than
# four alpha carbons is a number without a meaning, so it is refused rather than reported.
# This is a geometric floor, not a design-quality threshold -- 037's binder_length default
# is 50-120 and nothing here second-guesses it.
MIN_BINDER_RESIDUES = 4

# Not a limit, a warning line. The only timings this workflow has measured are 35 residues
# in 61 s and a 63-residue complex in 80 s on an RTX 3060 Ti (tests/reference/), which says
# nothing about where memory runs out. AF2 at this size is practically bounded by VRAM
# rather than by anything in this node, so the node says so once and proceeds; §11.5 is the
# open question, and predict_timeout_s is the thing that actually stops a runaway.
LONG_SEQUENCE_WARN_RESIDUES = 1000

MONO_PREFIX = "mono_"
FASTA_SUFFIX = ".fasta"

# Phase 2. TWO files per design, and they are not redundant.
#
#   cplx_<design>.fasta     TARGET:BINDER, one record. The DECLARATION: it is what was
#                           asked for, in the form the spec validates (a ':' is required
#                           here and forbidden in a mono_ file -- 026's bug from both
#                           directions), and it is human-readable.
#   pairedmsa_<design>.a3m  What the predictor actually consumes. §11.1's measured format:
#                           the target carries its alignment and the binder is present
#                           only as its own query.
#
# The a3m could not simply replace the fasta, because the a3m has no ':' in it -- its
# chain boundary is the '#<len1>,<len2>' header -- so the colon checks would have nothing
# to read and the spec's rule would become unenforceable. And the fasta could not replace
# the a3m, because a bare TARGET:BINDER fasta makes colabfold_batch compute an MSA for
# BOTH chains, which is the thing §4 decided against: the designed binder must never be
# aligned or sent anywhere.
#
# So node 03 is given both and cross-checks them: the a3m's declared lengths must match
# the fasta's two halves. Two representations of one fact, each validated against the
# other, is the opposite of a single unchecked one.
CPLX_PREFIX = "cplx_"
PAIRED_MSA_PREFIX = "pairedmsa_"
A3M_SUFFIX = ".a3m"

# unserialize_msa()'s chain ids. 101 is the first chain, 102 the second; they are not
# residue numbers and not the design's own chain letters.
A3M_CHAIN_IDS = ["101", "102"]


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def numeric_param(env_name, display_name, default, cast):
    """Read a numeric parameter, failing the same way every other bad parameter does."""
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    try:
        value = cast(raw)
    except ValueError:
        kind = "a whole number" if cast is int else "a number"
        fail(f"{display_name} must be {kind}, got '{raw}'.")
    if not math.isfinite(value):
        fail(
            f"{display_name} must be a finite number, got '{raw}', "
            f"which parses as {value}."
        )
    return value


def bool_param(env_name, display_name, default):
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    lowered = raw.strip().lower()
    if lowered in ("true", "1", "yes", "on"):
        return True
    if lowered in ("false", "0", "no", "off"):
        return False
    fail(f"{display_name} must be true or false, got '{raw}'.")


def load_inventory(path):
    if not os.path.exists(path):
        fail(
            f"No '{os.path.basename(path)}' in ./inputs/. This node reads node 01's "
            "inventory, which is where the designed chain was resolved; it does not parse "
            "structures itself, so there is no fallback and there should not be one -- a "
            "second parser is a second chance to pick a different chain (§8.1)."
        )
    try:
        with open(path) as handle:
            inventory = json.load(handle)
    except Exception as exc:
        fail(f"Could not read '{path}': {exc}")
    if not isinstance(inventory.get("designs"), list):
        fail(f"'{path}' has no 'designs' list; node 01 did not write it, or it is a stub.")
    return inventory


def check_rerun(out_dir, expected_names):
    """Refuse a stale prediction input before writing anything at all.

    Order matters: the refusal has to happen before the first write, or a run that is
    refused has still half-replaced the directory it refused to touch (4d1a523e).

    Every prefix this node writes is checked, not just mono_. A stale cplx_ file is the
    worse case of the two: node 03 globs both, so a complex input left over from a design
    that is not in this run would be predicted and scored as an interface, and node 05
    would report it under a design name the inventory has never heard of.
    """
    if not os.path.isdir(out_dir):
        return
    patterns = [
        (MONO_PREFIX, FASTA_SUFFIX),
        (CPLX_PREFIX, FASTA_SUFFIX),
        (PAIRED_MSA_PREFIX, A3M_SUFFIX),
    ]
    existing, expected = [], set()
    for prefix, suffix in patterns:
        existing += sorted(
            name for name in os.listdir(out_dir)
            if name.startswith(prefix) and name.endswith(suffix)
        )
        expected |= {f"{prefix}{name}{suffix}" for name in expected_names}
    stale = [name for name in existing if name not in expected]
    if stale:
        fail(
            f"./outputs/ already holds {len(stale)} prediction input(s) for design(s) "
            f"that are not in this run: {', '.join(stale)}. Node 03 globs "
            f"{MONO_PREFIX}*{FASTA_SUFFIX} and {CPLX_PREFIX}*{FASTA_SUFFIX}, so these "
            "would be predicted, and node 04 would then score designs the inventory has "
            "never heard of. Nothing has been written and nothing has been removed -- "
            "clear ./outputs/ or run in a fresh directory. (029's node 04 re-ran over a "
            "populated directory during a GPU smoke test and swapped already-swapped "
            "chains back: 4d1a523e.)"
        )
    if existing:
        print(
            f"NOTE: ./outputs/ already holds {len(existing)} input file(s), all of them "
            "for designs in this run. Rewriting them: they are derived from the inventory "
            "and nothing else, so the same inventory produces the same bytes.",
            flush=True,
        )


def write_mono(out_dir, entry):
    """Write one single-record FASTA and return what was recorded about it.

    Returns a dict, or None if the design was skipped upstream.
    """
    design = entry["design"]
    sequence = entry["binder_sequence"]
    path = os.path.join(out_dir, f"{MONO_PREFIX}{design}{FASTA_SUFFIX}")

    with open(path, "w") as handle:
        handle.write(f">{design}\n{sequence}\n")

    return {
        "design": design,
        "path": os.path.basename(path),
        "mode": "mono",
        "record_id": design,
        # Recorded even though Phase 1's slice is the whole array. The value is what node
        # 04 slices by, and Phase 2 changes it; a node 04 that reads it from here cannot
        # inherit 029's head-slice by accident.
        "chain_order": ["binder"],
        "chain_lengths": [len(sequence)],
        "expected_chains": 1,
        "expected_residues": len(sequence),
        "binder_chain_in_design": entry["binder_chain"],
        "designed_chain_source": entry["designed_chain_source"],
        "n_unknown": entry["n_unknown"],
    }


def match_columns(row):
    """An a3m row reduced to its match columns: uppercase residues and '-' deletions.

    a3m encodes an insertion relative to the query as a LOWERCASE letter, which occupies
    no alignment column. Dropping them makes every row exactly as long as the query, which
    is the shape §11.1 measured working -- the captured reference file's homolog rows are
    len(chain1) match columns followed by len(chain2) gaps, with no lowercase anywhere.
    Keeping the insertions would be more faithful and is not worth the risk: the format
    had to be read out of unserialize_msa()'s source as it was, after two guesses were
    rejected outright.
    """
    return "".join(c for c in row if c == "-" or c.isupper())


def write_cplx(out_dir, entry, target):
    """cplx_<design>.fasta -- TARGET:BINDER, one record.

    TARGET FIRST, deliberately, and this is the decision node 04 depends on. The scores
    JSON carries no chain boundaries, so node 04 slices the per-residue plddt array by
    position. Workflow-029 slices plddt[:binder_length] because ITS binder is chain A; in
    this order the binder is the TAIL, and a node 04 that copied 029's slice would average
    the TARGET's pLDDT and report it as the binder's. The order is a fact this node
    decides, so this node records it (§8.1) and node 04 reads it rather than assuming.
    """
    design = entry["design"]
    binder = entry["binder_sequence"]
    path = os.path.join(out_dir, f"{CPLX_PREFIX}{design}{FASTA_SUFFIX}")

    with open(path, "w") as handle:
        handle.write(f">{design}\n{target['sequence']}:{binder}\n")

    return {
        "design": design,
        "path": os.path.basename(path),
        "mode": "cplx",
        "record_id": design,
        "chain_order": ["target", "binder"],
        "chain_lengths": [target["n_residues"], len(binder)],
        "expected_chains": 2,
        "expected_residues": target["n_residues"] + len(binder),
        "binder_chain_in_design": entry["binder_chain"],
        "designed_chain_source": entry["designed_chain_source"],
        "n_unknown": entry["n_unknown"],
        "target_chains": target["chains"],
        "target_n_residues": target["n_residues"],
        "target_fusion": target["fusion"],
    }


def write_paired_msa(out_dir, entry, target):
    """pairedmsa_<design>.a3m -- the target aligned, the binder as its own query only.

    §11.1's measured format, which had to be read out of colabfold.batch's
    unserialize_msa() after two hand-rolled attempts were rejected with
    'MSA 0 must contain at least one sequence'. The rule both missed:

        any row carrying residues in two or more chain segments goes to the PAIRED block
        for every chain, and only single-chain rows populate unpaired_chunks[j] -- so
        EVERY chain needs at least one single-chain row, including its own query.

    A concatenated query alone leaves every unpaired block empty and generate_input_feature
    raises. Hence rows 2 and 3 below, which look redundant and are load-bearing.
    """
    design = entry["design"]
    binder = entry["binder_sequence"]
    t_len, b_len = target["n_residues"], len(binder)
    path = os.path.join(out_dir, f"{PAIRED_MSA_PREFIX}{design}{A3M_SUFFIX}")

    lines = [f"#{t_len},{b_len}\t1,1"]
    # The paired row: residues in both segments, so unserialize_msa sends it to the paired
    # block for both chains. It is also a3m_lines[2], which gets sliced by the declared
    # lengths to recover each chain's query -- so the declared lengths and this row have
    # to agree exactly or the chains come back wrong without an error.
    lines += [f">{A3M_CHAIN_IDS[0]}", f"{target['sequence']}{binder}"]
    # One single-chain row per chain, each seeding its own unpaired block.
    lines += [f">{A3M_CHAIN_IDS[0]}", f"{target['sequence']}{'-' * b_len}"]
    lines += [f">{A3M_CHAIN_IDS[1]}", f"{'-' * t_len}{binder}"]

    # The target's homologs, padded into chain 1's segment. The binder gets none: it is a
    # de novo design with no homologs, and §4's decision is that it is never aligned.
    n_homologs = 0
    for name, row in target["alignment"][1:]:
        columns = match_columns(row)
        if len(columns) != t_len:
            # Skipped rather than padded or truncated: a row that does not line up with
            # the query is not an alignment of it, and silently reshaping it would put
            # the wrong residue in a column and change what the target looks like.
            print(
                f"  NOTE: {design}: dropped alignment row {name!r} -- {len(columns)} "
                f"match column(s) against a {t_len}-residue target query. A row that "
                "does not line up with the query is not an alignment of it.",
                flush=True,
            )
            continue
        safe = name.replace(">", "_").replace("\t", " ").strip() or f"hom{n_homologs}"
        lines += [f">{safe}", f"{columns}{'-' * b_len}"]
        n_homologs += 1

    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")

    return os.path.basename(path), n_homologs


def verify_paired_msa(path, record):
    """Read the a3m back and check the three things that fail silently.

    Each of these produces a prediction rather than an error if it is wrong, which is the
    §8.3 case: a mis-declared length gives two chains of the wrong sizes, and node 04 then
    slices a plddt array by numbers that do not describe it.
    """
    t_len, b_len = record["chain_lengths"]
    with open(path) as handle:
        raw = [line.rstrip("\n") for line in handle if line.strip()]

    if not raw or not raw[0].startswith("#"):
        fail(f"'{os.path.basename(path)}' has no '#<len1>,<len2>' header.")
    declared = raw[0][1:].split("\t")[0]
    if declared != f"{t_len},{b_len}":
        fail(
            f"'{os.path.basename(path)}' declares chain lengths {declared!r} but the "
            f"record says {t_len},{b_len}. unserialize_msa slices the query row by the "
            "DECLARED lengths, so a mismatch returns two chains of the wrong sizes and "
            "no error at all."
        )

    rows = [(raw[i][1:], raw[i + 1]) for i in range(1, len(raw) - 1, 2)
            if raw[i].startswith(">")]
    if not rows:
        fail(f"'{os.path.basename(path)}' holds no sequence rows.")

    query = rows[0][1]
    if len(query) != t_len + b_len:
        fail(
            f"'{os.path.basename(path)}' query row is {len(query)} column(s), not "
            f"{t_len + b_len}. It is sliced by the declared lengths, so this silently "
            "mis-assigns residues between the chains."
        )

    # The §11.1 rule, checked rather than assumed: every chain needs a single-chain row.
    singles = {0: 0, 1: 0}
    for _, row in rows:
        if len(row) != t_len + b_len:
            continue
        head, tail = row[:t_len], row[t_len:]
        head_has = any(c != "-" for c in head)
        tail_has = any(c != "-" for c in tail)
        if head_has and not tail_has:
            singles[0] += 1
        elif tail_has and not head_has:
            singles[1] += 1
    missing = [i for i, n in singles.items() if n == 0]
    if missing:
        fail(
            f"'{os.path.basename(path)}' has no single-chain row for chain(s) "
            f"{[A3M_CHAIN_IDS[i] for i in missing]}. unserialize_msa only populates "
            "unpaired_chunks[j] from single-chain rows, so generate_input_feature raises "
            "'MSA 0 must contain at least one sequence' -- which is exactly how §11.1's "
            "two hand-rolled attempts failed."
        )

    # A colon here would be read as a chain break on top of the declared lengths.
    if any(":" in row for _, row in rows):
        fail(f"'{os.path.basename(path)}' contains a colon; its chain split is the header.")

    return {"n_rows": len(rows), "declared": declared,
            "single_chain_rows": {A3M_CHAIN_IDS[i]: n for i, n in singles.items()}}


def load_target(alignment_path, a3m_path):
    """The target and its alignment, from 01b. Refused rather than guessed at.

    01b writes target_alignment.json on every run, including a Phase 1 one where it says
    it was a no-op. So this can tell "no alignment was asked for" from "an alignment was
    asked for and failed", which an absent file could not (§8.3).
    """
    if not os.path.exists(alignment_path):
        fail(
            "score_interface = true but there is no target_alignment.json. "
            "01b-align-target must run first; it is a dependency of this node, so an "
            "absent file means that node did not run rather than that it found nothing."
        )
    with open(alignment_path) as handle:
        record = json.load(handle)

    if not record.get("aligned"):
        fail(
            "score_interface = true but 01b-align-target reports no alignment: "
            f"{record.get('reason') or 'no reason recorded'}. This node will not write a "
            "complex input without one. §4 measured what predicting a complex at "
            "single_sequence does -- iPTM 0.100 against 0.860, receptor pLDDT 38.4 "
            "against 83.7 -- so the numbers would be about the prediction mode rather "
            "than the design, and nothing downstream would say so."
        )
    if not os.path.exists(a3m_path):
        fail(
            f"target_alignment.json says it aligned the target, but {a3m_path} is not "
            "here. Refusing rather than proceeding with an alignment that is only a "
            "claim in a JSON file."
        )

    header, records = read_a3m_rows(a3m_path)
    if header is not None and "," in header:
        fail(f"{a3m_path} carries a paired a3m header; 01b writes a single-chain one.")
    if not records:
        fail(f"{a3m_path} holds no sequences.")

    sequence = record["target_sequence"]
    if match_columns(records[0][1]).replace("-", "") != sequence:
        fail(
            f"{a3m_path}'s query is not the target sequence target_alignment.json "
            "records. One of the two is describing a different protein, and pairing the "
            "binder against the wrong receptor would produce an interface number for it."
        )
    if len(records) < 2:
        fail(
            f"{a3m_path} holds {len(records)} sequence(s). An alignment of the query "
            "alone is single_sequence under another name; 01b refuses this too and "
            "refuses it first, but this node is reachable on its own."
        )

    return {
        "sequence": sequence,
        "n_residues": record["target_n_residues"],
        "chains": record.get("target_chains") or [],
        "fusion": record.get("target_fusion") or [],
        "alignment": records,
        "n_sequences": len(records),
        "a3m_sha256": record.get("a3m_sha256"),
        "msa_source": record.get("source"),
    }


def read_a3m_rows(path):
    """(header, [(name, sequence)]). Hand-rolled: see 01b for why not biopython."""
    header, name, chunks, out = None, None, [], []
    with open(path, errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                if header is None:
                    header = line
                continue
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(chunks)))
                name, chunks = line[1:].strip(), []
                continue
            chunks.append(line)
    if name is not None:
        out.append((name, "".join(chunks)))
    return header, out


def verify_written(path, expect_colons):
    """Check the FILE, not the string that was meant to be written.

    A guard that re-reads the variable it just assigned asserts the assignment. These
    checks exist because a colon miscount costs GPU minutes and produces a plausible
    number, so they read the bytes back off disk (§8.4).
    """
    with open(path) as handle:
        text = handle.read()

    records = [line for line in text.splitlines() if line.startswith(">")]
    if len(records) != 1:
        fail(
            f"'{os.path.basename(path)}' holds {len(records)} record(s), not 1. One design "
            "per file is the whole defence against 026's conflation bug, where a batch "
            "file made three designs' ranks report as one design's top 5 (818954e9)."
        )

    body = "\n".join(
        line for line in text.splitlines() if line and not line.startswith(">")
    )
    colons = body.count(":")
    if colons != expect_colons:
        if expect_colons == 0:
            fail(
                f"'{os.path.basename(path)}' contains {colons} colon(s) in its sequence. A "
                "colon makes the predictor treat one design as two chains, so this "
                "monomer input would come back as a complex and its confidence would be "
                "read as an interface -- which is 026's bug arriving from the other "
                "direction (d3d707ed)."
            )
        fail(
            f"'{os.path.basename(path)}' contains {colons} colon(s) in its sequence, not "
            f"{expect_colons}. Zero silently predicts a monomer and reports its confidence "
            "as an interface; two predicts three chains."
        )

    if ":" in records[0]:
        fail(
            f"The record id in '{os.path.basename(path)}' contains a colon. 026 reads a "
            "colon in a record id or description as a chain-break marker, which is why "
            "037's UNSCREENED tag carries none (7709809b)."
        )
    return len(records), colons


def main():
    score_interface = bool_param("PARAM_SCORE_INTERFACE", "score_interface", False)
    numeric_param("PARAM_NUM_MODELS", "num_models", 5, int)

    os.makedirs("./outputs", exist_ok=True)
    inventory = load_inventory("./inputs/design_inventory.json")
    designs = inventory["designs"]
    if not designs:
        fail("node 01's inventory lists no designs. There is nothing to predict.")

    eligible = [e for e in designs if not e.get("skipped")]
    skipped = [e for e in designs if e.get("skipped")]

    # Validate every design before writing any file, so a batch with one bad sequence does
    # not leave a half-written output directory behind.
    for entry in eligible:
        sequence = entry["binder_sequence"]
        bad = sorted({c for c in sequence if c not in ALPHABET})
        if bad:
            fail(
                f"Design '{entry['design']}' has residue code(s) {bad} outside the "
                f"accepted alphabet. 026 validates every sequence against "
                f"{ALPHABET_STRING} and refuses anything else, so this would be "
                "refused inside the predictor after the GPU had already been allocated "
                "(d3d707ed). Node 01 scores an unrecognised residue as X; a code that is "
                "neither a standard residue nor X means something between the two nodes "
                "disagrees about the alphabet."
            )
        if len(sequence) < MIN_BINDER_RESIDUES:
            fail(
                f"Design '{entry['design']}' has {len(sequence)} residue(s), fewer than "
                f"{MIN_BINDER_RESIDUES}. A Kabsch superposition needs three non-collinear "
                "points to be defined and four before the fit stops being exactly "
                "determined by its own inputs, so the self-consistency RMSD this workflow "
                "exists to report would be a number without a meaning."
            )
        if len(entry["binder_ca"]) != len(sequence):
            fail(
                f"Design '{entry['design']}' has {len(sequence)} residue(s) but "
                f"{len(entry['binder_ca'])} coordinate slot(s). Node 01 keeps those arrays "
                "index-aligned on purpose, because node 04's RMSD maps the two structures "
                "by position; a mismatch here means they were not."
            )

    # Loaded before anything is written, so a Phase 2 run with no usable alignment fails
    # with an empty output directory rather than a half-written one.
    target = None
    if score_interface:
        target = load_target(
            "./inputs/target_alignment.json", "./inputs/target_msa.a3m"
        )
        print(
            f"Phase 2: a {target['n_residues']}-residue target over chain(s) "
            f"{target['chains']}, {target['n_sequences']} alignment sequence(s) from "
            f"{target['msa_source']}. The binder is paired against it UNALIGNED and is "
            "not sent anywhere (§4).",
            flush=True,
        )
        for line in target["fusion"]:
            print(f"  WARNING: {line} (§8.12).", flush=True)

    check_rerun("./outputs", [e["design"] for e in eligible])

    records = []
    for entry in eligible:
        # The mono input is written in BOTH phases. Phase 2 adds interface scoring, it
        # does not replace self-consistency: an interface metric says nothing about
        # whether the binder folds as designed, which is the question this workflow was
        # built to answer.
        record = write_mono("./outputs", entry)
        n_records, n_colons = verify_written(
            os.path.join("./outputs", record["path"]), expect_colons=0
        )
        record["verified_records_in_file"] = n_records
        record["verified_colons_in_sequence"] = n_colons
        records.append(record)

        if not score_interface:
            continue

        cplx = write_cplx("./outputs", entry, target)
        n_records, n_colons = verify_written(
            os.path.join("./outputs", cplx["path"]), expect_colons=1
        )
        cplx["verified_records_in_file"] = n_records
        cplx["verified_colons_in_sequence"] = n_colons

        msa_name, n_homologs = write_paired_msa("./outputs", entry, target)
        cplx["paired_msa"] = msa_name
        cplx["n_target_homologs_in_msa"] = n_homologs
        cplx["paired_msa_verified"] = verify_paired_msa(
            os.path.join("./outputs", msa_name), cplx
        )
        cplx["binder_aligned"] = False
        cplx["target_a3m_sha256"] = target["a3m_sha256"]
        records.append(cplx)

    payload = {
        "n_inputs": len(records),
        "n_skipped": len(skipped),
        "phase": 2 if score_interface else 1,
        "mode": "mono+cplx" if score_interface else "mono",
        "score_interface": score_interface,
        # Recorded whether or not Phase 2 ran, so node 03 and node 05 can say which
        # question was asked rather than inferring it from which files exist.
        "target": None if not score_interface else {
            "chains": target["chains"],
            "n_residues": target["n_residues"],
            "n_alignment_sequences": target["n_sequences"],
            "msa_source": target["msa_source"],
            "a3m_sha256": target["a3m_sha256"],
            "fusion": target["fusion"],
            "binder_aligned": False,
        },
        "alphabet": ALPHABET_STRING,
        "min_binder_residues": MIN_BINDER_RESIDUES,
        # Carried forward from node 01 rather than recomputed, so one wording of the chain
        # provenance reaches the report.
        "designed_chain_sources": inventory.get("designed_chain_sources"),
        "n_target_numbering_gaps": inventory.get("n_target_numbering_gaps"),
        "securedna_screened": False,
        "biosafety_notice": (
            "These sequences are unscreened de novo designs. Passing validation makes "
            "synthesis more likely, not less."
        ),
        "inputs": records,
        "skipped": [
            {"design": e["design"], "reason": e.get("skip_reason")} for e in skipped
        ],
    }
    with open("./outputs/prediction_inputs.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    for record in records:
        extra = ""
        if record["mode"] == "cplx":
            verified = record["paired_msa_verified"]
            extra = (
                f"; {record['paired_msa']} declares {verified['declared']} over "
                f"{verified['n_rows']} row(s) with "
                f"{record['n_target_homologs_in_msa']} target homolog(s) and no binder "
                f"alignment, single-chain rows {verified['single_chain_rows']}"
            )
        print(
            f"  {record['path']}: 1 record, {record['expected_residues']} residue(s), "
            f"{record['verified_colons_in_sequence']} colon(s), "
            f"chain_order={record['chain_order']} lengths={record['chain_lengths']}"
            + (f", {record['n_unknown']} X" if record["n_unknown"] else "")
            + extra,
            flush=True,
        )
    for entry in skipped:
        print(f"  SKIPPED {entry['design']}: {entry.get('skip_reason')}", flush=True)

    long_ones = [r for r in records
                 if r["expected_residues"] > LONG_SEQUENCE_WARN_RESIDUES]
    if long_ones:
        print(
            f"WARNING: {len(long_ones)} design(s) exceed {LONG_SEQUENCE_WARN_RESIDUES} "
            "residues. This workflow has measured 35 residues in 61 s and a 63-residue "
            "complex in 80 s on an RTX 3060 Ti, and nothing beyond that, so it cannot say "
            "where this runs out of VRAM -- only that it is the dimension that decides "
            "(§11.5). predict_timeout_s is what stops a runaway, not this line.",
            flush=True,
        )
    if skipped:
        print(
            f"WARNING: {len(skipped)} design(s) were skipped by node 01 and have no "
            "prediction input. They stay in prediction_inputs.json under 'skipped' with "
            "the reason, so node 05 reports them as rows with a stated cause rather than "
            "omitting them (§8.5).",
            flush=True,
        )

    print(
        f"Wrote {len(records)} prediction input(s) and ./outputs/prediction_inputs.json",
        flush=True,
    )


if __name__ == "__main__":
    main()
